#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

#define NX  8
#define NY  8
#define N   64
#define PI  3.14159265358979323846

// Use standard complex for testbench (not HLS types)
typedef complex<double> cplx;

// ================================================
// DECLARE TOP FUNCTION
// (Vitis links this to mvdr_top.cpp automatically)
// ================================================
#include "hls_x_complex.h"
typedef hls::x_complex<float> hls_cplx;

void mvdr_top(
    float az_target, float el_target,
    float az_i1,     float el_i1,
    float az_i2,     float el_i2,
    hls_cplx R_in[N][N],
    float g_target[1],
    float g_i1_out[1],
    float g_i2_out[1]
);

// ================================================
// HELPER: build steering vector (for testbench)
// ================================================
void sv(double az_deg, double el_deg, double d, hls_cplx a[N])
{
    double az = az_deg * PI / 180.0;
    double el = el_deg * PI / 180.0;
    double u  = sin(az) * cos(el);
    double v  = sin(el);
    double nv = 1.0 / sqrt((double)N);
    int k = 0;
    for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++) {
            double phase = 2.0 * PI * d * (x*u + y*v);
            a[k].real( (float)(nv * cos(phase)) );
            a[k].imag( (float)(nv * sin(phase)) );
            k++;
        }
}

// ================================================
// MAIN TESTBENCH
// ================================================
int main()
{
    cout << "========================================" << endl;
    cout << "   MVDR HLS TESTBENCH                  " << endl;
    cout << "========================================" << endl;

    // ---------- Build covariance matrix ----------
    hls_cplx R_in[N][N];
    hls_cplx a_t[N], a_i1[N], a_i2[N];

    // Steering vectors
    sv( 40.0,   0.0, 0.5, a_t);
    sv( 30.0,  10.0, 0.5, a_i1);
    sv(-60.0, -25.0, 0.5, a_i2);

    // R = 25*a_i1*a_i1^H + 25*a_i2*a_i2^H + 0.5*I
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float r = 25.0f * (a_i1[i].real()*a_i1[j].real()
                             + a_i1[i].imag()*a_i1[j].imag())
                    + 25.0f * (a_i2[i].real()*a_i2[j].real()
                             + a_i2[i].imag()*a_i2[j].imag());

            float im = 25.0f * (a_i1[i].imag()*a_i1[j].real()
                              - a_i1[i].real()*a_i1[j].imag())
                     + 25.0f * (a_i2[i].imag()*a_i2[j].real()
                              - a_i2[i].real()*a_i2[j].imag());

            R_in[i][j].real(r);
            R_in[i][j].imag(im);
        }
        R_in[i][i].real( R_in[i][i].real() + 0.5f );
    }

    // Diagonal loading
    float trace = 0;
    for (int i = 0; i < N; i++)
        trace += R_in[i][i].real();
    float delta = 0.01f * trace / N;
    for (int i = 0; i < N; i++)
        R_in[i][i].real( R_in[i][i].real() + delta );

    cout << "Covariance matrix ready." << endl;

    // ---------- Call top function ----------
    float g_target[1], g_i1_out[1], g_i2_out[1];

    cout << "Calling mvdr_top..." << endl;

    mvdr_top(
        40.0f,   0.0f,
        30.0f,  10.0f,
       -60.0f, -25.0f,
        R_in,
        g_target, g_i1_out, g_i2_out
    );

    // ---------- Print results ----------
    cout << "\n===== RESULTS =====" << endl;
    cout << "Target gain (should be ~1.0): " << g_target[0]  << endl;
    cout << "Interferer 1 gain:            " << g_i1_out[0]  << endl;
    cout << "Interferer 2 gain:            " << g_i2_out[0]  << endl;

    cout << "\nSuppression (dB):" << endl;
    cout << "I1 rejection: "
         << 20*log10(g_i1_out[0]/g_target[0]) << " dB" << endl;
    cout << "I2 rejection: "
         << 20*log10(g_i2_out[0]/g_target[0]) << " dB" << endl;

    // ---------- Pass/Fail check ----------
    cout << "\n===== PASS/FAIL =====" << endl;

    int pass = 1;

    if (fabs(g_target[0] - 1.0f) < 0.01f)
        cout << "Target gain:  PASS ✅" << endl;
    else {
        cout << "Target gain:  FAIL ❌ (got " << g_target[0] << ")" << endl;
        pass = 0;
    }

    if (g_i1_out[0] < 0.1f)
        cout << "I1 rejection: PASS ✅" << endl;
    else {
        cout << "I1 rejection: FAIL ❌ (got " << g_i1_out[0] << ")" << endl;
        pass = 0;
    }

    if (g_i2_out[0] < 0.1f)
        cout << "I2 rejection: PASS ✅" << endl;
    else {
        cout << "I2 rejection: FAIL ❌ (got " << g_i2_out[0] << ")" << endl;
        pass = 0;
    }

    cout << "\n==============================" << endl;
    if (pass)
        cout << "ALL TESTS PASSED ✅" << endl;
    else
        cout << "SOME TESTS FAILED ❌" << endl;
    cout << "==============================" << endl;

    return (pass ? 0 : 1);
}
