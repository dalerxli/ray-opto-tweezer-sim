#include "optics.h"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES

// Defines for CUBA library
#include <cuba.h>

#define NDIM 2
#define NCOMP 3
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

#define STATEFILE NULL
#define SPIN NULL

#define KEY1 47
#define KEY2 1
#define KEY3 1

#define KEY 0

double intens(double r, double th);

Vector3d get_total_force(Sphere* s, Lens* l);
int integrand(const int *ndim, const double xx[],
    const int *ncomp, double ff[], void *userdata);

int main(int argc, char* argv[]);
