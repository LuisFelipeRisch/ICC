#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "float.h"
#include "utils.h"
#include "ZeroFuncao.h"

double bisseccao(double xl, double xu, double f(double x), double epsilon)
{
	double xmOld, xmNew;
	uint counter = 1;

	xmNew = (xl + xu) / 2;

	if (f(xl) * f(xmNew) > 0)
		xu = xmNew;
	else if (f(xl) * f(xmNew) < 0)
		xl = xmNew;
	else
		return xmNew;

	do
	{
		xmOld = xmNew;
		xmNew = (xl + xu) / 2;

		if (f(xl) * f(xmNew) > 0)
			xu = xmNew;
		else if (f(xl) * f(xmNew) < 0)
			xl = xmNew;
		else
			return xmOld;

		counter++;
	} while (fabs((xmNew - xmOld) * 100 / xmNew) > epsilon && counter < MAXIT);

	return xmNew;
}

double newtonRaphson(double x0, double f(double x), double derivative_f(double x), double epsilon)
{
	double xi = x0;
	double xiOld;
	uint counter = 0;

	do
	{
		xiOld = xi;
		xi = xi - f(xi) / derivative_f(xi);
		counter++;

	} while (fabs((xi - xiOld) * 100 / xi) > epsilon && counter < MAXIT);

	return xi;
}

double secante(double x0, double x1, double f(double x), double epsilon)
{
	double xOld, xNew, aux;
	uint counter = 0;

	xOld = x0;
	xNew = x1;

	do
	{
		aux = xNew;
		xNew = xNew - (f(xNew) * (xNew - xOld)) / (f(xNew) - f(xOld));
		xOld = aux;
		counter++;

	} while (fabs((xNew - xOld) * 100 / xNew) > epsilon && counter < MAXIT);

	return xNew;
}

void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx)
{
}

void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx)
{
}
