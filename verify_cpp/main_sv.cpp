#include <iostream>
#include <cmath>
#include "constants.h"
#include "steering_vec.h"

using namespace std;

int main()
{
    cout << "===== STAGE 1: STEERING VECTOR =====" << endl;

    cplx a_cpp[N];
    steering_vec(AZ_TARGET, EL_TARGET, D, a_cpp);

    // MATLAB golden reference (first 15 elements)
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

    // Compare element by element
    cout << "\ni  | MATLAB               | C++                  | Error      | Match" << endl;
    cout << "---|---------------------|----------------------|------------|------" << endl;

    int pass_count = 0;
    for (int i = 0; i < 15; i++) {
        double err_r = fabs(a_cpp[i].real() - matlab_real[i]);
        double err_i = fabs(a_cpp[i].imag() - matlab_imag[i]);
        double err   = sqrt(err_r*err_r + err_i*err_i);
        bool match = (err < 1e-3);
        if (match) pass_count++;

        printf("%2d | %6.4f + %6.4fj | %6.4f + %6.4fj | %.4e | %s\n",
            i+1,
            matlab_real[i], matlab_imag[i],
            a_cpp[i].real(), a_cpp[i].imag(),
            err,
            match ? "PASS" : "FAIL"
        );
    }

    // Overall relative error
    double diff = 0, norm_m = 0;
    for (int i = 0; i < 15; i++) {
        cplx m(matlab_real[i], matlab_imag[i]);
        diff   += norm(a_cpp[i] - m);
        norm_m += norm(m);
    }
    double rel_err = sqrt(diff) / sqrt(norm_m);

    cout << "\n===== RESULT =====" << endl;
    printf("Relative error  = %.6e\n", rel_err);
    printf("Elements passed = %d / 15\n", pass_count);

    if (pass_count == 15 && rel_err < 1e-3)
        cout << "STAGE 1 STEERING VECTOR: PASS ✅" << endl;
    else
        cout << "STAGE 1 STEERING VECTOR: FAIL ❌" << endl;

    return 0;
}
