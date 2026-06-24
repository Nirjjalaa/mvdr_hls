#include "qr.h"
#include <cmath>

void givens_qr(cplx A[][N], cplx Q[][N], cplx R[][N])
{
    // Initialize Q = identity, R = A
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
            Q[i][j] = (i == j) ? cplx(1,0) : cplx(0,0);
        }

    for (int j = 0; j < N; j++) {
        for (int i = N-1; i >= j+1; i--) {

            if (abs(R[i][j]) < 1e-12) continue;

            cplx a = R[i-1][j];
            cplx b = R[i][j];

            double r = sqrt(norm(a) + norm(b));  // norm(x) = |x|^2

            // Complex Givens coefficients (matching MATLAB exactly)
            double c;
            cplx alpha, s;

            if (abs(a) < 1e-12) {
                alpha = cplx(1, 0);
            } else {
                alpha = a / abs(a);   // phase of a
            }

            c = abs(a) / r;
            s = alpha * conj(b) / r;

            // 2x2 rotation matrix
            // G = [ c      s  ]
            //     [-conj(s)  c]
            cplx G00 = cplx(c, 0);
            cplx G01 = s;
            cplx G10 = -conj(s);
            cplx G11 = cplx(c, 0);

            // Apply to rows i-1 and i of R
            for (int k = j; k < N; k++) {
                cplx tmp1 = G00 * R[i-1][k] + G01 * R[i][k];
                cplx tmp2 = G10 * R[i-1][k] + G11 * R[i][k];
                R[i-1][k] = tmp1;
                R[i][k]   = tmp2;
            }

            // Apply G' to columns i-1 and i of Q
            // G' = [c        -conj(s)]
            //      [conj(s)   c      ]
            for (int k = 0; k < N; k++) {
                cplx tmp1 = Q[k][i-1] * cplx(c,0) + Q[k][i] * conj(s);
                cplx tmp2 = Q[k][i-1] * (-s)       + Q[k][i] * cplx(c,0);
                Q[k][i-1] = tmp1;
                Q[k][i]   = tmp2;
            }
        }
    }
}
