#include "steering_vec.h"
#include <cmath>

void steering_vec(double az_deg, double el_deg, double d, cplx a[N])
{
    double az = az_deg * PI / 180.0;
    double el = el_deg * PI / 180.0;

    double u = sin(az) * cos(el);
    double v = sin(el);

    double norm_val = 1.0 / sqrt((double)N);

    int k = 0;
    for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++) {
            double phase = 2.0 * PI * d * (x*u + y*v);
            a[k++] = cplx(norm_val * cos(phase),
                          norm_val * sin(phase));
        }
}
