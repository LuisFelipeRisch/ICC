#include <stdio.h>
#include <math.h>
#include "functions.h"
#include "utils.h"
#include "ZeroFuncao.h"

int main()
{
  double lowerLimit = 0.0;
  double upperLimit = 3.0;

  printf("Zero da f1 bisseccao: %.15f\n", bisseccao(lowerLimit, upperLimit, f1, EPS));
  printf("Zero da f1 newtonRaphson: %.15f\n", newtonRaphson(upperLimit, f1, derivative_f1, EPS));

  return 0;
}
