#include "functions.h"
#include "math.h"

double f1(double x) { return 2 * pow(x, 4) + 4 * pow(x, 3) + 3 * pow(x, 2) - 10 * pow(x, 1) - 15; }
double derivative_f1(double x) { return 8 * pow(x, 3) + 12 * pow(x, 2) + 6 * pow(x, 1) - 10; }