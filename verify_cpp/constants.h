#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
using namespace std;

// ================================================
// ARRAY PARAMETERS
// ================================================
#define NX  8           // elements in x direction
#define NY  8           // elements in y direction
#define N   64          // total elements = NX * NY
#define K   5000        // number of snapshots

// ================================================
// PHYSICAL PARAMETERS
// ================================================
#define D   0.5         // element spacing (half wavelength)

// ================================================
// MATH
// ================================================
#define PI  3.14159265358979323846

// ================================================
// SIGNAL POWERS
// ================================================
#define TARGET_POWER  10.0    // amplitude of target
#define INTERF_POWER   5.0    // amplitude of interferers
#define NOISE_POWER    0.5    // noise variance

// ================================================
// ANGLES (degrees)
// ================================================
#define AZ_TARGET   40.0
#define EL_TARGET    0.0
#define AZ_I1       30.0
#define EL_I1       10.0
#define AZ_I2      -60.0
#define EL_I2      -25.0

// ================================================
// VERIFICATION THRESHOLDS
// ================================================
#define THRESH_TIGHT   1e-6    // for steering vec, weights
#define THRESH_LOOSE   1e-3    // for beam pattern

// ================================================
// COMPLEX TYPE
// ================================================
typedef complex<double> cplx;

// ================================================
// CSV PATH
// ================================================
#define CSV_PATH "C:/mvdr_hls/verify_data/"

#endif
