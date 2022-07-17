#include "solve.h"

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