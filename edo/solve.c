#include <stdlib.h>
#include "solve.h"
#include "utils.h"

triDiagonal_SL *allocTriDiagonal(uint n)
{
  triDiagonal_SL *triDiagonalSL = (triDiagonal_SL *)malloc(sizeof(triDiagonal_SL));

  triDiagonalSL->n = n;

  triDiagonalSL->mainDiagonal = allocArray(n);
  triDiagonalSL->upperDiagonal = allocArray(n - 1);
  triDiagonalSL->lowerDiagonal = allocArray(n - 1);
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

void gaussElimination(triDiagonal_SL *triDiagonalSL)
{
  uint size = triDiagonalSL->n;

  for (uint i = 0; i < size; i++)
  {
    real_t m = triDiagonalSL->lowerDiagonal[i] / triDiagonalSL->mainDiagonal[i];
    triDiagonalSL->lowerDiagonal[i] = 0.0;
    triDiagonalSL->mainDiagonal[i + 1] -= m * triDiagonalSL->upperDiagonal[i];
    triDiagonalSL->independetTerms[i + 1] -= m * triDiagonalSL->independetTerms[i];
  }
}

void retroSubstitution(triDiagonal_SL *triDiagonalSL,
                       real_t *solution)
{
  uint size = triDiagonalSL->n;

  solution[size - 1] = triDiagonalSL->independetTerms[size - 1] / triDiagonalSL->mainDiagonal[size - 1];
  for (uint i = size - 2; i >= 0; i--)
    solution[i] = (triDiagonalSL->independetTerms[i] - triDiagonalSL->upperDiagonal[i] * solution[i + 1]) / triDiagonalSL->mainDiagonal[i];
}