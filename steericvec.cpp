#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

#define NX 8
#define NY 8
#define N  64
#define PI 3.14159265358979323846

typedef complex<double> cplx;

void steering_vec(double az_deg, double el_deg, double d, cplx a[N])
{
    double az = az_deg * PI / 180.0;
    double el = el_deg * PI / 180.0;

    double u = sin(az) * cos(el);
    double v = sin(el);

    double norm = 1.0 / sqrt((double)N);

    int k = 0;

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {

            double phase = 2.0 * PI * d * (x*u + y*v);

            a[k] = cplx(norm * cos(phase), norm * sin(phase));
            k++;
        }
    }
}

int main()
{
    cplx a[N];

    // Same parameters as your MATLAB code
    steering_vec(40.0, 0.0, 0.5, a);

    // Print first 5 elements to compare with MATLAB
    cout << "Steering vector (first 5 elements):" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "a[" << i << "] = "
             << a[i].real() << " + "
             << a[i].imag() << "j" << endl;
    }

    // Print magnitude of all elements (should all be 1/sqrt(64) = 0.125)
    cout << "\nMagnitude check (all should be 0.125):" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "|a[" << i << "]| = " << abs(a[i]) << endl;
    }

    return 0;
}
