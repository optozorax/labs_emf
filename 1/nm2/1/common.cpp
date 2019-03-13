#include <cmath>
#include "common.h"

//-----------------------------------------------------------------------------
bool isNear(double a, double b) {
	if (a != 0) {
		if (fabs(a - b)/a > 0.0001)
			return false;
	} else {
		if (fabs(b) > 0.0001)
			return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
double random(void) {
	return std::rand() / double(RAND_MAX);
}

//-----------------------------------------------------------------------------
int intRandom(int min, int max) {
	return min + random() * (max - min);
}