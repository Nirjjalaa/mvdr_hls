#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

#define N 4

typedef complex<double> cplx;

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

int main()
{
    cplx R[N][N] = {};
    cplx b[N] = {};
    cplx x[N] = {};

    // Fill upper triangular R (4x4)
    R[0][0]=2; R[0][1]=1; R[0][2]=3; R[0][3]=1;
               R[1][1]=4; R[1][2]=2; R[1][3]=2;
                          R[2][2]=5; R[2][3]=3;
                                     R[3][3]=2;

    // Fill b
    b[0]=10; b[1]=8; b[2]=15; b[3]=4;

    back_substitution(R, b, x);

    cout << "Solution x:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[" << i << "] = "
             << x[i].real() << " + "
             << x[i].imag() << "j" << endl;
    }

    cout << "\nVerification R*x (should match b):" << endl;
    for (int i = 0; i < N; i++) {
        cplx sum = 0;
        for (int j = 0; j < N; j++) {
            sum += R[i][j] * x[j];
        }
        cout << "R*x[" << i << "] = " << sum.real()
             << " (b[" << i << "] = " << b[i].real() << ")" << endl;
    }

    return 0;
}
