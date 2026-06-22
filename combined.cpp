#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

#define NX  8
#define NY  8
#define N   64
#define PI  3.14159265358979323846

typedef complex<double> cplx;

// ================================================
// FUNCTION 1: STEERING VECTOR
// ================================================
void steering_vec(double az_deg, double el_deg, double d, cplx a[N])
{
    double az = az_deg * PI / 180.0;
    double el = el_deg * PI / 180.0;

    double u = sin(az) * cos(el);
    double v = sin(el);

    double norm_val = 1.0 / sqrt((double)N);

    int k = 0;
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            double phase = 2.0 * PI * d * (x*u + y*v);
            a[k] = cplx(norm_val * cos(phase), norm_val * sin(phase));
            k++;
        }
    }
}

// ================================================
// FUNCTION 2: GIVENS QR DECOMPOSITION
// ================================================
void givens_qr(cplx A[N][N], cplx Q[N][N], cplx R[N][N])
{
    // Initialize Q = identity, R = A
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
            Q[i][j] = (i == j) ? cplx(1,0) : cplx(0,0);
        }

    // Givens rotations column by column, bottom up
    for (int j = 0; j < N; j++) {
        for (int i = N-1; i >= j+1; i--) {

            if (abs(R[i][j]) < 1e-12) continue;

            double a = R[i-1][j].real();
            double b = R[i][j].real();
            double r = sqrt(a*a + b*b);
            double c =  a / r;
            double s = -b / r;

            // Rotate rows of R
            for (int k = j; k < N; k++) {
                cplx tmp1 =  c * R[i-1][k] - s * R[i][k];
                cplx tmp2 =  s * R[i-1][k] + c * R[i][k];
                R[i-1][k] = tmp1;
                R[i][k]   = tmp2;
            }

            // Rotate columns of Q
            for (int k = 0; k < N; k++) {
                cplx tmp1 = c * Q[k][i-1] - s * Q[k][i];
                cplx tmp2 = s * Q[k][i-1] + c * Q[k][i];
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
    for (int i = N-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < N; j++) {
            x[i] = x[i] - R[i][j] * x[j];
        }
        x[i] = x[i] / R[i][i];
    }
}

// ================================================
// HELPER: MATRIX-VECTOR MULTIPLY → y = A*x
// ================================================
void mat_vec_mul(cplx A[N][N], cplx x[N], cplx y[N])
{
    for (int i = 0; i < N; i++) {
        y[i] = 0;
        for (int j = 0; j < N; j++)
            y[i] += A[i][j] * x[j];
    }
}

// ================================================
// HELPER: DOT PRODUCT → result = a^H * b
// ================================================
cplx dot_product(cplx a[N], cplx b[N])
{
    cplx result = 0;
    for (int i = 0; i < N; i++)
        result += conj(a[i]) * b[i];
    return result;
}

// ================================================
// MAIN: FULL MVDR PIPELINE
// ================================================
int main()
{
    // ---------- STEP 1: Steering vectors ----------
    cplx a_target[N], a_i1[N], a_i2[N];

    steering_vec( 40.0,   0.0,  0.5, a_target);
    steering_vec( 30.0,  10.0,  0.5, a_i1);
    steering_vec(-60.0, -25.0,  0.5, a_i2);

    cout << "Steering vectors computed." << endl;

    // ---------- STEP 2: Build covariance matrix R ----------
    // R = a_i1*a_i1^H * 25 + a_i2*a_i2^H * 25 + 0.5*I
    // (simplified: no snapshots, just signal model)
    // interf_power = 25 (amplitude 5 squared)
    // noise_power  = 0.5

    cplx R[N][N] = {};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = 25.0 * a_i1[i] * conj(a_i1[j])
                    + 25.0 * a_i2[i] * conj(a_i2[j]);
        }
        R[i][i] += cplx(0.5, 0);   // add noise on diagonal
    }

    // Diagonal loading
    double trace = 0;
    for (int i = 0; i < N; i++)
        trace += R[i][i].real();
    double delta = 0.01 * trace / N;
    for (int i = 0; i < N; i++)
        R[i][i] += cplx(delta, 0);

    cout << "Covariance matrix built." << endl;

    // ---------- STEP 3: Givens QRD ----------
    cplx Q[N][N] = {};
    cplx Ru[N][N] = {};

    givens_qr(R, Q, Ru);

    cout << "Givens QRD complete." << endl;

    // ---------- STEP 4: Solve for weights ----------
    // y = Q^H * a_target
    cplx y[N] = {};
    for (int i = 0; i < N; i++) {
        y[i] = 0;
        for (int j = 0; j < N; j++)
            y[i] += conj(Q[j][i]) * a_target[j];  // Q^H = conj transpose
    }

    // w = back_substitution(Ru, y)
    cplx w[N] = {};
    back_substitution(Ru, y, w);

    // ---------- STEP 5: Normalize weights ----------
    // w = w / (a_target^H * w)
    cplx norm_factor = dot_product(a_target, w);
    for (int i = 0; i < N; i++)
        w[i] = w[i] / norm_factor;

    cout << "Weights computed." << endl;

    // ---------- STEP 6: Gain check ----------
    cplx g_target = dot_product(w, a_target);
    cplx g_i1     = dot_product(w, a_i1);
    cplx g_i2     = dot_product(w, a_i2);

    cout << "\n===== MVDR GAINS =====" << endl;
    cout << "Target gain (should be ~1.0): " << abs(g_target) << endl;
    cout << "Interferer 1 gain:            " << abs(g_i1)     << endl;
    cout << "Interferer 2 gain:            " << abs(g_i2)     << endl;

    cout << "\nSuppression (dB):" << endl;
    cout << "I1 rejection: "
         << 20*log10(abs(g_i1)/abs(g_target)) << " dB" << endl;
    cout << "I2 rejection: "
         << 20*log10(abs(g_i2)/abs(g_target)) << " dB" << endl;

    return 0;
}
