#include "backsub.h"

void back_substitution(cplx R[][N], cplx b[], cplx x[])
{
    for (int i = N-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < N; j++)
            x[i] -= R[i][j] * x[j];
        x[i] /= R[i][i];
    }
}
