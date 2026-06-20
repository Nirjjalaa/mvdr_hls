#include "hls_math.h"
#include "hls_x_complex.h"
#include "ap_fixed.h"

#define NX  8
#define NY  8
#define N   64
#define PI  3.14159265358979323846f

typedef hls::x_complex<float> cplx;

// ================================================
// FUNCTION 1: STEERING VECTOR
// ================================================
void steering_vec(float az_deg, float el_deg, float d, cplx a[N])
{
#pragma HLS INLINE

    float az = az_deg * PI / 180.0f;
    float el = el_deg * PI / 180.0f;

    float u = hls::sin(az) * hls::cos(el);
    float v = hls::sin(el);

    float norm_val = 1.0f / hls::sqrt((float)N);

    int k = 0;

    loop_y:
    for (int y = 0; y < NY; y++) {
    #pragma HLS PIPELINE
        loop_x:
        for (int x = 0; x < NX; x++) {
            float phase = 2.0f * PI * d * (x*u + y*v);
            a[k].real( norm_val * hls::cos(phase) );
            a[k].imag( norm_val * hls::sin(phase) );
            k++;
        }
    }
}

// ================================================
// FUNCTION 2: GIVENS QR DECOMPOSITION
// ================================================
void givens_qr(cplx A[N][N], cplx Q[N][N], cplx R[N][N])
{
#pragma HLS INLINE

    // Initialize Q = identity, R = A
    init_i:
    for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE
        init_j:
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
            if (i == j) {
                Q[i][j].real(1.0f);
                Q[i][j].imag(0.0f);
            } else {
                Q[i][j].real(0.0f);
                Q[i][j].imag(0.0f);
            }
        }
    }

    // Givens rotations
    col_loop:
    for (int j = 0; j < N; j++) {
        row_loop:
        for (int i = N-1; i >= j+1; i--) {

            float a = R[i-1][j].real();
            float b = R[i][j].real();

            if (b < 1e-6f && b > -1e-6f) continue;

            float r = hls::sqrt(a*a + b*b);
            float c =  a / r;
            float s = -b / r;

            // Rotate rows of R
            rot_R:
            for (int k = j; k < N; k++) {
            #pragma HLS PIPELINE
                cplx tmp1, tmp2;
                tmp1.real( c * R[i-1][k].real() - s * R[i][k].real() );
                tmp1.imag( c * R[i-1][k].imag() - s * R[i][k].imag() );
                tmp2.real( s * R[i-1][k].real() + c * R[i][k].real() );
                tmp2.imag( s * R[i-1][k].imag() + c * R[i][k].imag() );
                R[i-1][k] = tmp1;
                R[i][k]   = tmp2;
            }

            // Rotate columns of Q
            rot_Q:
            for (int k = 0; k < N; k++) {
            #pragma HLS PIPELINE
                cplx tmp1, tmp2;
                tmp1.real( c * Q[k][i-1].real() - s * Q[k][i].real() );
                tmp1.imag( c * Q[k][i-1].imag() - s * Q[k][i].imag() );
                tmp2.real( s * Q[k][i-1].real() + c * Q[k][i].real() );
                tmp2.imag( s * Q[k][i-1].imag() + c * Q[k][i].imag() );
                Q[k][i-1] = tmp1;
                Q[k][i]   = tmp2;
            }
        }
    }
}

// ================================================
// FUNCTION 3: BACK SUBSTITUTION
// ================================================
void back_substitution(cplx R[N][N], cplx b[N], cplx x[N])
{
#pragma HLS INLINE

    bs_loop:
    for (int i = N-1; i >= 0; i--) {
        cplx sum;
        sum.real(0.0f);
        sum.imag(0.0f);

        inner_loop:
        for (int j = i+1; j < N; j++) {
        #pragma HLS PIPELINE
            sum.real( sum.real() + R[i][j].real()*x[j].real()
                                 - R[i][j].imag()*x[j].imag() );
            sum.imag( sum.imag() + R[i][j].real()*x[j].imag()
                                 + R[i][j].imag()*x[j].real() );
        }

        cplx diff;
        diff.real( b[i].real() - sum.real() );
        diff.imag( b[i].imag() - sum.imag() );

        float denom = R[i][i].real()*R[i][i].real()
                    + R[i][i].imag()*R[i][i].imag();

        x[i].real( (diff.real()*R[i][i].real()
                  + diff.imag()*R[i][i].imag()) / denom );
        x[i].imag( (diff.imag()*R[i][i].real()
                  - diff.real()*R[i][i].imag()) / denom );
    }
}

// ================================================
// TOP FUNCTION (entry point for Vitis HLS)
// ================================================
void mvdr_top(
    float  az_target,       // target azimuth
    float  el_target,       // target elevation
    float  az_i1,           // interferer 1 azimuth
    float  el_i1,           // interferer 1 elevation
    float  az_i2,           // interferer 2 azimuth
    float  el_i2,           // interferer 2 elevation
    cplx   R_in[N][N],      // input covariance matrix
    float  g_target[1],     // output: gain at target
    float  g_i1_out[1],     // output: gain at interferer 1
    float  g_i2_out[1]      // output: gain at interferer 2
)
{
// Interface pragmas — tells Vitis how to connect ports
#pragma HLS INTERFACE s_axilite port=az_target
#pragma HLS INTERFACE s_axilite port=el_target
#pragma HLS INTERFACE s_axilite port=az_i1
#pragma HLS INTERFACE s_axilite port=el_i1
#pragma HLS INTERFACE s_axilite port=az_i2
#pragma HLS INTERFACE s_axilite port=el_i2
#pragma HLS INTERFACE m_axi port=R_in     offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=g_target offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=g_i1_out offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=g_i2_out offset=slave bundle=gmem
#pragma HLS INTERFACE s_axilite port=return

    // Local arrays (on-chip BRAM)
    cplx a_t[N], a_i1[N], a_i2[N];
    cplx Q[N][N], Ru[N][N];
    cplx y[N], w[N];

#pragma HLS ARRAY_PARTITION variable=a_t   complete
#pragma HLS ARRAY_PARTITION variable=a_i1  complete
#pragma HLS ARRAY_PARTITION variable=a_i2  complete

    // Step 1: steering vectors
    steering_vec(az_target, el_target, 0.5f, a_t);
    steering_vec(az_i1,     el_i1,     0.5f, a_i1);
    steering_vec(az_i2,     el_i2,     0.5f, a_i2);

    // Step 2: Givens QRD
    givens_qr(R_in, Q, Ru);

    // Step 3: y = Q^H * a_target
    y_loop:
    for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE
        y[i].real(0.0f);
        y[i].imag(0.0f);
        for (int j = 0; j < N; j++) {
            y[i].real( y[i].real() + Q[j][i].real()*a_t[j].real()
                                   + Q[j][i].imag()*a_t[j].imag() );
            y[i].imag( y[i].imag() + Q[j][i].real()*a_t[j].imag()
                                   - Q[j][i].imag()*a_t[j].real() );
        }
    }

    // Step 4: back substitution
    back_substitution(Ru, y, w);

    // Step 5: normalize
    float norm_r = 0.0f, norm_i = 0.0f;
    for (int i = 0; i < N; i++) {
        norm_r += a_t[i].real()*w[i].real() + a_t[i].imag()*w[i].imag();
        norm_i += a_t[i].real()*w[i].imag() - a_t[i].imag()*w[i].real();
    }
    float norm_mag2 = norm_r*norm_r + norm_i*norm_i;
    for (int i = 0; i < N; i++) {
        float wr = w[i].real(), wi = w[i].imag();
        w[i].real( (wr*norm_r + wi*norm_i) / norm_mag2 );
        w[i].imag( (wi*norm_r - wr*norm_i) / norm_mag2 );
    }

    // Step 6: compute gains
    float gt_r=0, gt_i=0, gi1_r=0, gi1_i=0, gi2_r=0, gi2_i=0;
    for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE
        gt_r  += w[i].real()*a_t[i].real()  + w[i].imag()*a_t[i].imag();
        gt_i  += w[i].real()*a_t[i].imag()  - w[i].imag()*a_t[i].real();
        gi1_r += w[i].real()*a_i1[i].real() + w[i].imag()*a_i1[i].imag();
        gi1_i += w[i].real()*a_i1[i].imag() - w[i].imag()*a_i1[i].real();
        gi2_r += w[i].real()*a_i2[i].real() + w[i].imag()*a_i2[i].imag();
        gi2_i += w[i].real()*a_i2[i].imag() - w[i].imag()*a_i2[i].real();
    }

    g_target[0]  = hls::sqrt(gt_r*gt_r   + gt_i*gt_i);
    g_i1_out[0]  = hls::sqrt(gi1_r*gi1_r + gi1_i*gi1_i);
    g_i2_out[0]  = hls::sqrt(gi2_r*gi2_r + gi2_i*gi2_i);
}
