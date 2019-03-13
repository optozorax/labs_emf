#pragma once

#ifdef ALL_FLOAT
typedef float real;
typedef float sumreal;
#endif

#ifdef ALL_DOUBLE
typedef double real;
typedef double sumreal;
#endif

#ifdef ALL_FLOAT_WITH_DOUBLE
typedef float real;
typedef double sumreal;
#endif

#ifndef ALL_FLOAT
#ifndef ALL_DOUBLE
#ifndef ALL_FLOAT_WITH_DOUBLE
#error "Type isn't defined"
typedef double real;
typedef double sumreal;
#endif
#endif
#endif


bool isNear(double a, double b);
double random(void);
int intRandom(int min, int max);