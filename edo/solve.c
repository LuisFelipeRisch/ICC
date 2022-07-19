#include <stdio.h>
#include <stdlib.h>
#include "solve.h"
#include "utils.h"
#include "math.h"

triDiagonal_SL *allocTriDiagonal(uint n)
{
  triDiagonal_SL *triDiagonalSL = (triDiagonal_SL *)malloc(sizeof(triDiagonal_SL));

  triDiagonalSL->n = n;

  triDiagonalSL->mainDiagonal = allocArray(n);
  triDiagonalSL->upperDiagonal = allocArray(n);
  triDiagonalSL->lowerDiagonal = allocArray(n);
  triDiagonalSL->independetTerms = allocArray(n);

  return triDiagonalSL;
}

void freeTriDiagonal(triDiagonal_SL *triDiagonalSL)
{
  free(triDiagonalSL->mainDiagonal);
  free(triDiagonalSL->upperDiagonal);
  free(triDiagonalSL->lowerDiagonal);
  free(triDiagonalSL->independetTerms);

  free(triDiagonalSL);
}

real_t L2norm(real_t *residue, uint size)
{
  real_t sum = 0.0;
  for (uint i = 0; i < size; i++)
    sum += residue[i] * residue[i];
  return sqrt(sum);
}

void calculateResidue(triDiagonal_SL *triDiagonalSL, real_t *solution, real_t *residue)
{
  uint size = triDiagonalSL->n;

  residue[0] = triDiagonalSL->independetTerms[0] - (triDiagonalSL->mainDiagonal[0] * solution[0] + triDiagonalSL->upperDiagonal[0] * solution[1]);
  for (uint i = 1; i <= size - 2; i++)
    residue[i] = triDiagonalSL->independetTerms[i] - (triDiagonalSL->mainDiagonal[i] * solution[i] + triDiagonalSL->upperDiagonal[i] * solution[i + 1] + triDiagonalSL->lowerDiagonal[i - 1] * solution[i - 1]);
  residue[size - 1] = triDiagonalSL->independetTerms[size - 1] - (triDiagonalSL->mainDiagonal[size - 1] * solution[size - 1] + triDiagonalSL->lowerDiagonal[size - 2] * solution[size - 2]);
}

/* needed to be fixed */
void gaussElimination(triDiagonal_SL *triDiagonalSL)
{
  uint size = triDiagonalSL->n;

  for (uint i = 0; i < size - 1; i++)
  {
    real_t m = triDiagonalSL->lowerDiagonal[i] / triDiagonalSL->mainDiagonal[i];
    triDiagonalSL->lowerDiagonal[i] = 0.0;
    triDiagonalSL->mainDiagonal[i + 1] -= m * triDiagonalSL->upperDiagonal[i];
    triDiagonalSL->independetTerms[i + 1] -= m * triDiagonalSL->independetTerms[i];
  }
}

/* needed to be fixed */
void retroSubstitution(triDiagonal_SL *triDiagonalSL,
                       real_t *solution)
{
  uint size = triDiagonalSL->n;

  solution[size - 1] = triDiagonalSL->independetTerms[size - 1] / triDiagonalSL->mainDiagonal[size - 1];
  for (uint i = size - 2; i >= 0; i--)
    solution[i] = (triDiagonalSL->independetTerms[i] - triDiagonalSL->upperDiagonal[i] * solution[i + 1]) / triDiagonalSL->mainDiagonal[i];
}

void gaussSeidel(triDiagonal_SL *triDiagonalSL, real_t *solution, real_t tolerance)
{
  uint counter = 0;
  uint size = triDiagonalSL->n;
  real_t error = 1.0 + tolerance;
  real_t *c = triDiagonalSL->upperDiagonal;
  real_t *d = triDiagonalSL->mainDiagonal;
  real_t *a = triDiagonalSL->lowerDiagonal;
  real_t *b = triDiagonalSL->independetTerms;
  real_t *residue = allocArray(size);

  for (uint i = 0; i < size; i++)
    solution[i] = 0;

  while (error > tolerance && counter < MAX_ITERATION)
  {
    solution[0] = (b[0] - c[0] * solution[1]) / d[0];
    for (uint i = 1; i <= size - 2; i++)
      solution[i] = (b[i] - a[i - 1] * solution[i - 1] - c[i - 1] * solution[i + 1]) / d[i];
    solution[size - 1] = (b[size - 1] - a[size - 2] * solution[size - 2]) / d[size - 1];
    calculateResidue(triDiagonalSL, solution, residue);
    error = L2norm(residue, size);
    counter++;
  }
  printf("couter: %d\n", counter);
  free(residue);
}

void differentGaussSeidel(Edo_EQ *edoEquation, triDiagonal_SL *triDiagonalSL, real_t *solution, real_t tolerance)
{
  uint size;
  real_t d, di, ds, b;
  real_t lastDI, lastDS;
  real_t xi, h;
  real_t error = 1.0 + tolerance;
  uint counter = 0;

  size = edoEquation->n;
  h = (edoEquation->b - edoEquation->a) / (size + 1.0);

  for (uint i = 0; i < size; i++)
    solution[i] = 0;

  while (counter < MAX_ITERATION && error > tolerance)
  {
    real_t *residue = allocArray(size);
    for (uint i = 0; i < size; i++)
    {
      xi = edoEquation->a + (i + 1) * h;

      di = 1 - h * edoEquation->p(xi) / 2.0;
      d = -2 + h * h * edoEquation->q(xi);
      ds = 1 + h * edoEquation->p(xi) / 2.0;
      b = h * h * edoEquation->r(xi);

      if (i == 0)
        b -= edoEquation->ya * (1 - h * edoEquation->p(edoEquation->a + h) / 2.0);
      if (i == size - 1)
        b -= edoEquation->yb * (1 + h * edoEquation->p(edoEquation->b - h) / 2.0);

      // calculating solution
      if (i == 0)
        solution[0] = (b - ds * solution[1]) / d;
      else if (i == size - 1)
        solution[size - 1] = (b - lastDI * solution[size - 2]) / d;
      else
        solution[i] = (b - lastDI * solution[i - 1] - lastDS * solution[i + 1]) / d;

      // calculating residue - needed to be reviwed
      // if (i == 0)
      //   residue[0] = b - (d * solution[0] + ds * solution[1]);
      // else if (i == size - 1)
      //   residue[size - 1] = b - (d * solution[size - 1] + lastDI * solution[size - 2]);
      // else
      //   residue[i] = b - (d * solution[i] + lastDI * solution[i - 1] + ds * solution[i + 1]);

      lastDI = di;
      lastDS = ds;
    }
    calculateResidue(triDiagonalSL, solution, residue);
    error = L2norm(residue, size);
    free(residue);
    counter++;
  }
}

void printTriDiagonalMatrix(triDiagonal_SL *triDiagonalSL)
{
  uint size = triDiagonalSL->n;

  for (uint i = 0; i < size; i++)
  {
    for (uint j = 0; j < size; j++)
    {
      if (i == j)
        printf("%.15g ", triDiagonalSL->mainDiagonal[i]);
      else if (j == i + 1)
        printf("%.15g ", triDiagonalSL->upperDiagonal[i]);
      else if (j == i - 1)
        printf("%.15g ", triDiagonalSL->lowerDiagonal[i]);
      else
        printf("%d ", 0);
    }
    printf(" | %.15g", triDiagonalSL->independetTerms[i]);
    printf("\n");
  }
}