#include <iostream>
#include <cmath>
#include "constants.h"
#include "steering_vec.h"
#include "matrix.h"
#include "qr.h"

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

    return 0;    // ← always at the very end
}
