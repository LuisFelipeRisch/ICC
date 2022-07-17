#include <math.h>
#include "coefficients.h"

real_t p1(real_t x)
{
  return 0.0;
}

real_t q1(real_t x)
{
  return 0.0;
}

real_t r1(real_t x)
{
  return (6 * x) - (0.5 * pow(x, 2));
}

real_t p2(real_t x)
{
  return 0.0;
}

real_t q2(real_t x)
{
  return 1.0;
}

real_t r2(real_t x)
{
  return 0.0;
}

real_t p3(real_t x)
{
  return (x + 1);
}

real_t q3(real_t x)
{
  return 2.0;
}

real_t r3(real_t x)
{
  return (1 - pow(x, 2)) * (exp(-x));
}