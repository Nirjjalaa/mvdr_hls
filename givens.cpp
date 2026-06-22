#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

#define N 4

typedef complex<double> cplx;

void givens_qr(cplx A[N][N], cplx Q[N][N], cplx R[N][N])
{
    // Step 1: Initialize Q = identity, R = A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
            Q[i][j] = (i == j) ? cplx(1,0) : cplx(0,0);
        }
    }

    // Step 2: Givens rotations
    for (int j = 0; j < N; j++) {
        for (int i = N-1; i >= j+1; i--) {

            if (abs(R[i][j]) < 1e-12) continue;

            double a = R[i-1][j].real();
            double b = R[i][j].real();

            double r = sqrt(a*a + b*b);

            double c =  a / r;
            double s = -b / r;

            // Apply rotation to R: rows i-1 and i
            for (int k = j; k < N; k++) {
                cplx tmp1 =  c * R[i-1][k] - s * R[i][k];
                cplx tmp2 =  s * R[i-1][k] + c * R[i][k];
                R[i-1][k] = tmp1;
                R[i][k]   = tmp2;
            }

            // Apply rotation to Q: columns i-1 and i
            for (int k = 0; k < N; k++) {
                cplx tmp1 = c * Q[k][i-1] - s * Q[k][i];
                cplx tmp2 = s * Q[k][i-1] + c * Q[k][i];
                Q[k][i-1] = tmp1;
                Q[k][i]   = tmp2;
            }
        }
    }
}

void mat_mul(cplx A[N][N], cplx B[N][N], cplx C[N][N])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

int main()
{
    cplx A[N][N] = {
        {4, 2, 1, 3},
        {2, 5, 2, 1},
        {1, 2, 6, 2},
        {3, 1, 2, 4}
    };

    cplx Q[N][N] = {};
    cplx R[N][N] = {};

    givens_qr(A, Q, R);

    // Print R
    cout << "R matrix (should be upper triangular):" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double val = R[i][j].real();
            if (abs(val) < 1e-10) val = 0.0;  // clean near-zero values
            cout << val;
            if (j < N-1) cout << "\t";
        }
        cout << endl;
    }

    // Verify Q*R = A
    cplx QR[N][N] = {};
    mat_mul(Q, R, QR);

    cout << "\nVerification Q*R (should match original A):" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double val = QR[i][j].real();
            if (abs(val) < 1e-10) val = 0.0;
            cout << val;
            if (j < N-1) cout << "\t";
        }
        cout << endl;
    }

    cout << "\nOriginal A:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j].real();
            if (j < N-1) cout << "\t";
        }
        cout << endl;
    }

    return 0;
}
