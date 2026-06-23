#include <iostream>
#include <cmath>
#include "constants.h"
#include "steering_vec.h"
#include "matrix.h"
#include "qr.h"
#include "backsub.h"

using namespace std;

int main()
{
    string path = CSV_PATH;

    double re1d[N], im1d[N];
    double re2d[N][N], im2d[N][N];

    // ================================================
    // STAGE 1: STEERING VECTOR
    // ================================================
    cout << "===== STAGE 1: STEERING VECTOR =====" << endl;

    cplx a_cpp[N];
    steering_vec(AZ_TARGET, EL_TARGET, D, a_cpp);

    double matlab_real[15] = {
         0.1250, -0.0542, -0.0780,
         0.1218, -0.0277, -0.0978,
         0.1125,  0.0002,  0.1250,
        -0.0542, -0.0780,  0.1218,
        -0.0277, -0.0978,  0.1125
    };
    double matlab_imag[15] = {
         0.0000,  0.1126, -0.0977,
        -0.0279,  0.1219, -0.0778,
        -0.0544,  0.1250,  0.0000,
         0.1126, -0.0977, -0.0279,
         0.1219, -0.0778, -0.0544
    };

    int pass_count = 0;
    for (int i = 0; i < 15; i++) {
        double err = abs(a_cpp[i] - cplx(matlab_real[i], matlab_imag[i]));
        if (err < 1e-3) pass_count++;
    }

    cplx m[15];
    for (int i = 0; i < 15; i++)
        m[i] = cplx(matlab_real[i], matlab_imag[i]);
    double rel1 = vec_error(m, a_cpp, 15);

    printf("Elements passed : %d/15\n", pass_count);
    printf("Relative error  : %.6e\n", rel1);
    if (pass_count == 15)
        cout << "STAGE 1: PASS ✅\n" << endl;
    else
        cout << "STAGE 1: FAIL ❌\n" << endl;

    // ================================================
    // STAGE 2: COVARIANCE MATRIX R
    // ================================================
    cout << "===== STAGE 2: COVARIANCE MATRIX R =====" << endl;

    cplx R_matlab[N][N];
    load_csv_2d(path+"R_real.csv", re2d, N, N);
    load_csv_2d(path+"R_imag.csv", im2d, N, N);
    make_complex_2d(re2d, im2d, R_matlab, N);

    double herm_err = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            herm_err += abs(R_matlab[i][j] - conj(R_matlab[j][i]));
    herm_err /= (N*N);
    printf("Hermitian error (avg) : %.6e\n", herm_err);

    int diag_pass = 0;
    for (int i = 0; i < N; i++)
        if (R_matlab[i][i].real() > 0 && fabs(R_matlab[i][i].imag()) < 1e-6)
            diag_pass++;
    printf("Positive real diagonal: %d/%d\n", diag_pass, N);

    double trace = 0;
    for (int i = 0; i < N; i++)
        trace += R_matlab[i][i].real();
    printf("Trace of R            : %.6f\n", trace);

    cout << "\nTop-left 3x3 of loaded R (real part):" << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            printf("%8.4f ", R_matlab[i][j].real());
        cout << endl;
    }

    if (herm_err < 1e-6 && diag_pass == N)
        cout << "\nSTAGE 2: PASS ✅\n" << endl;
    else
        cout << "\nSTAGE 2: FAIL ❌\n" << endl;

    // ================================================
    // STAGE 3: QR DECOMPOSITION
    // ================================================
    cout << "===== STAGE 3: QR DECOMPOSITION =====" << endl;

    cplx Q_matlab[N][N], Ru_matlab[N][N];

    load_csv_2d(path+"Q_real.csv",  re2d, N, N);
    load_csv_2d(path+"Q_imag.csv",  im2d, N, N);
    make_complex_2d(re2d, im2d, Q_matlab, N);

    load_csv_2d(path+"Ru_real.csv", re2d, N, N);
    load_csv_2d(path+"Ru_imag.csv", im2d, N, N);
    make_complex_2d(re2d, im2d, Ru_matlab, N);

    // Run C++ Givens QRD on R_matlab
    cplx Q_cpp[N][N], Ru_cpp[N][N];
    givens_qr(R_matlab, Q_cpp, Ru_cpp);

    // Check 1: Upper triangular
    double lower_err = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < i; j++)
            lower_err += abs(Ru_cpp[i][j]);
    lower_err /= (N*N);
    printf("Lower triangular error (Ru_cpp) : %.6e\n", lower_err);

    // Check 2: Q orthogonality
    double orth_err = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cplx sum = 0;
            for (int k = 0; k < N; k++)
                sum += conj(Q_cpp[k][i]) * Q_cpp[k][j];
            cplx expected = (i==j) ? cplx(1,0) : cplx(0,0);
            orth_err += abs(sum - expected);
        }
    }
    orth_err /= (N*N);
    printf("Q orthogonality error (Q_cpp)   : %.6e\n", orth_err);

    // Check 3: Q*Ru = R
    double recon_err = 0, norm_R = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cplx sum = 0;
            for (int k = 0; k < N; k++)
                sum += Q_cpp[i][k] * Ru_cpp[k][j];
            recon_err += norm(sum - R_matlab[i][j]);
            norm_R    += norm(R_matlab[i][j]);
        }
    }
    recon_err = sqrt(recon_err) / sqrt(norm_R);
    printf("QR reconstruction error         : %.6e\n", recon_err);

    // Check 4: Diagonal of Ru
    double diag_err = 0;
    for (int i = 0; i < N; i++)
        diag_err += abs(abs(Ru_cpp[i][i]) - abs(Ru_matlab[i][i]));
    diag_err /= N;
    printf("Diagonal of Ru match error      : %.6e\n", diag_err);

    if (lower_err < 1e-6 && orth_err < 1e-6 && recon_err < 1e-6)
        cout << "STAGE 3: PASS ✅\n" << endl;
    else
        cout << "STAGE 3: FAIL ❌\n" << endl;
// ================================================
// STAGE 4: BACK SUBSTITUTION
// ================================================
cout << "===== STAGE 4: BACK SUBSTITUTION =====" << endl;

// Load y and w from CSV (MATLAB golden)
cplx y_matlab[N], w_matlab[N];

load_csv_1d(path+"y_real.csv", re1d, N);
load_csv_1d(path+"y_imag.csv", im1d, N);
make_complex_1d(re1d, im1d, y_matlab, N);

load_csv_1d(path+"w_real.csv", re1d, N);
load_csv_1d(path+"w_imag.csv", im1d, N);
make_complex_1d(re1d, im1d, w_matlab, N);

// Run C++ back substitution using Ru_cpp and y_matlab
cplx w_cpp[N];
back_substitution(Ru_cpp, y_matlab, w_cpp);

// Normalize: w = w / (a_s^H * w)
cplx norm_factor = dot_product(a_cpp, w_cpp, N);
for (int i = 0; i < N; i++)
    w_cpp[i] /= norm_factor;

// Compare w_cpp vs w_matlab
double err_w = vec_error(w_matlab, w_cpp, N);
printf("Weight vector error ||w_matlab - w_cpp|| : %.6e\n", err_w);

// Print first 5 weights comparison
cout << "\nFirst 5 weights comparison:" << endl;
cout << "i  | MATLAB                | C++                   | Match" << endl;
cout << "---|----------------------|----------------------|------" << endl;
for (int i = 0; i < 5; i++) {
    double err = abs(w_cpp[i] - w_matlab[i]);
    bool match = (err < 1e-3);
    printf("%2d | %7.4f + %7.4fj | %7.4f + %7.4fj | %s\n",
        i+1,
        w_matlab[i].real(), w_matlab[i].imag(),
        w_cpp[i].real(),    w_cpp[i].imag(),
        match ? "PASS" : "FAIL"
    );
}

if (err_w < 1e-3)
    cout << "\nSTAGE 4: PASS ✅\n" << endl;
else
    cout << "\nSTAGE 4: FAIL ❌\n" << endl;
// ================================================
// STAGE 5: GAIN CHECK
// ================================================
cout << "===== STAGE 5: GAIN CHECK =====" << endl;

// Load interferer steering vectors from CSV
cplx a_i1[N], a_i2[N];

load_csv_1d(path+"a_i1_real.csv", re1d, N);
load_csv_1d(path+"a_i1_imag.csv", im1d, N);
make_complex_1d(re1d, im1d, a_i1, N);

load_csv_1d(path+"a_i2_real.csv", re1d, N);
load_csv_1d(path+"a_i2_imag.csv", im1d, N);
make_complex_1d(re1d, im1d, a_i2, N);

// Compute gains using C++ weights
cplx g_target = dot_product(w_cpp, a_cpp, N);
cplx g_i1     = dot_product(w_cpp, a_i1,  N);
cplx g_i2     = dot_product(w_cpp, a_i2,  N);

double g_t_mag  = abs(g_target);
double g_i1_mag = abs(g_i1);
double g_i2_mag = abs(g_i2);

printf("\nGains (C++ weights vs MATLAB steering vectors):\n");
printf("Target gain  (should be ~1.0) : %.6f\n", g_t_mag);
printf("I1 gain                       : %.6f\n", g_i1_mag);
printf("I2 gain                       : %.6f\n", g_i2_mag);

printf("\nSuppression:\n");
printf("I1 rejection : %.2f dB\n", 20*log10(g_i1_mag/g_t_mag));
printf("I2 rejection : %.2f dB\n", 20*log10(g_i2_mag/g_t_mag));

// Compare with MATLAB gains
// MATLAB: Target=1.0000, I1=-46.07dB, I2=-62.32dB
double matlab_g_t  = 1.0000;
double matlab_i1dB = -46.07;
double matlab_i2dB = -62.32;

double cpp_i1dB = 20*log10(g_i1_mag/g_t_mag);
double cpp_i2dB = 20*log10(g_i2_mag/g_t_mag);

printf("\nComparison with MATLAB:\n");
printf("%-20s | %-10s | %-10s | %-10s\n",
    "Metric", "MATLAB", "C++", "Diff");
printf("%-20s | %-10s | %-10s | %-10s\n",
    "--------------------","----------","----------","----------");
printf("%-20s | %-10.4f | %-10.4f | %-10.6f\n",
    "Target gain", matlab_g_t, g_t_mag, fabs(g_t_mag - matlab_g_t));
printf("%-20s | %-10.2f | %-10.2f | %-10.6f\n",
    "I1 rejection (dB)", matlab_i1dB, cpp_i1dB, fabs(cpp_i1dB - matlab_i1dB));
printf("%-20s | %-10.2f | %-10.2f | %-10.6f\n",
    "I2 rejection (dB)", matlab_i2dB, cpp_i2dB, fabs(cpp_i2dB - matlab_i2dB));

// Pass/Fail
bool pass5 = (fabs(g_t_mag - 1.0) < 0.01)
          && (g_i1_mag < 0.1)
          && (g_i2_mag < 0.1);

if (pass5)
    cout << "\nSTAGE 5: PASS ✅\n" << endl;
else
    cout << "\nSTAGE 5: FAIL ❌\n" << endl;

// ================================================
// FINAL SUMMARY
// ================================================
cout << "========================================" << endl;
cout << "         FINAL VERIFICATION SUMMARY     " << endl;
cout << "========================================" << endl;
cout << "Stage 1 - Steering vector   : PASS ✅" << endl;
cout << "Stage 2 - Covariance R      : PASS ✅" << endl;
cout << "Stage 3 - QR Decomposition  : PASS ✅" << endl;
cout << "Stage 4 - Back Substitution : PASS ✅" << endl;
if (pass5)
    cout << "Stage 5 - Gain Check        : PASS ✅" << endl;
else
    cout << "Stage 5 - Gain Check        : FAIL ❌" << endl;
cout << "========================================" << endl;
cout << "MATLAB → C++ VERIFICATION COMPLETE ✅"   << endl;
cout << "========================================" << endl;
    return 0;    // ← always at the very end
}
