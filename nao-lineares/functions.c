#include "functions.h"
#include "math.h"

double f1(double x) { return 2 * pow(x, 4) + 4 * pow(x, 3) + 3 * pow(x, 2) - 10 * pow(x, 1) - 15; }
double derivative_f1(double x) { return 8 * pow(x, 3) + 12 * pow(x, 2) + 6 * pow(x, 1) - 10; }

double f2(double x) { return 1 * pow(x, 5) - 2 * pow(x, 4) - 9 * pow(x, 3) + 22 * pow(x, 2) + 4 * pow(x, 1) - 24; }
double derivative_f2(double x) { return 5 * pow(x, 4) - 8 * pow(x, 3) - 27 * pow(x, 2) + 44 * pow(x, 1) + 4; }

double f3(double x) { return 3 * pow(x, 4) + 2 * pow(x, 3) + 4 * pow(x, 2) - 25 * pow(x, 1) - 30; }
double derivative_f3(double x) { return 12 * pow(x, 3) + 6 * pow(x, 2) + 8 * pow(x, 1) - 25; }

double f4(double x) { return 2 * pow(x, 3) - 5 * pow(x, 2) - 1 * pow(x, 1) + 3; }
double derivative_f4(double x) { return 6 * pow(x, 2) - 10 * pow(x, 1) - 1; }