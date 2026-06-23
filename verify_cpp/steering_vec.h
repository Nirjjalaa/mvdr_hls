#ifndef STEERING_VECTOR_H
#define STEERING_VECTOR_H

#include "constants.h"

// Function declaration
// az_deg → azimuth in degrees
// el_deg → elevation in degrees
// d      → element spacing
// a[]    → output steering vector (64 elements)

void steering_vec(double az_deg, double el_deg, double d, cplx a[N]);

#endif
