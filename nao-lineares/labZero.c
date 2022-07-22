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
  printf("Zero da f1 secante: %.15f\n", secante(upperLimit - 1.0, upperLimit, f1, EPS));

  printf("\n=========================================================\n");

  printf("Zero da f2 bisseccao: %.15f\n", bisseccao(lowerLimit, upperLimit, f2, EPS));
  printf("Zero da f2 newtonRaphson: %.15f\n", newtonRaphson(upperLimit, f2, derivative_f2, EPS));
  printf("Zero da f2 secante: %.15f\n", secante(upperLimit - 1.0, upperLimit, f2, EPS));

  printf("\n=========================================================\n");

  printf("Zero da f3 bisseccao: %.15f\n", bisseccao(lowerLimit, upperLimit, f3, EPS));
  printf("Zero da f3 newtonRaphson: %.15f\n", newtonRaphson(upperLimit, f3, derivative_f3, EPS));
  printf("Zero da f3 secante: %.15f\n", secante(upperLimit - 1.0, upperLimit, f3, EPS));

  printf("\n=========================================================\n");

  printf("Zero da f4 bisseccao: %.15f\n", bisseccao(lowerLimit, upperLimit, f4, EPS));
  printf("Zero da f4 newtonRaphson: %.15f\n", newtonRaphson(upperLimit, f4, derivative_f4, EPS));
  printf("Zero da f4 secante: %.15f\n", secante(upperLimit - 1.0, upperLimit, f4, EPS));

  return 0;
}
