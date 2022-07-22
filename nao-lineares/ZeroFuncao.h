#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include <float.h>

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 500
#define EPS 1.0e-8

typedef struct
{
  double *p;
  int grau;
} Polinomio;

// Métodos
double bisseccao(double xl, double xu, double f(double x), double epsilon);
double newtonRaphson(double x0, double f(double x), double derivative_f(double x), double epsilon);
double secante(double x0, double x1, double f(double x), double epsilon);

// Cálculo de Polinômios
void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx);
void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx);

#endif // __ZEROFUNCAO_H__
