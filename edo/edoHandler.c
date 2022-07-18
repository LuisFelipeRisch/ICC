#include "utils.h"
#include "edoHandler.h"
#include "solve.h"
#include <stdlib.h>

Edo_EQ *allocEdoEquation(uint n)
{
  Edo_EQ *edoEquation = (Edo_EQ *)malloc(sizeof(Edo_EQ));

  edoEquation->n = n;

  return edoEquation;
}

void freeEdoEquation(Edo_EQ *edoEquation)
{
  free(edoEquation);
}

void buildTriDiagonalSL(Edo_EQ *edoEquation,
                        triDiagonal_SL *triDiagonalSL)
{
  uint numberOfPoints;
  real_t xi, h;

  numberOfPoints = edoEquation->n;
  h = (edoEquation->b - edoEquation->a) / (numberOfPoints + 1.0);

  for (uint i = 0; i < numberOfPoints; i++)
  {
    xi = edoEquation->a + (i + 1) * h;

    triDiagonalSL->lowerDiagonal[i] = 1 - h * edoEquation->p(xi) / 2.0;
    triDiagonalSL->mainDiagonal[i] = -2 + h * h * edoEquation->q(xi);
    triDiagonalSL->upperDiagonal[i] = 1 + h * edoEquation->p(xi) / 2.0;
    triDiagonalSL->independetTerms[i] = h * h * edoEquation->r(xi);
  }

  triDiagonalSL->independetTerms[0] -= edoEquation->ya * (1 - h * edoEquation->p(edoEquation->a + h) / 2.0);
  triDiagonalSL->independetTerms[numberOfPoints - 1] -= edoEquation->yb * (1 + h * edoEquation->p(edoEquation->b - h) / 2.0);
}

void setEdoFunction(Edo_EQ *edoFunction,
                    real_t (*p)(real_t),
                    real_t (*q)(real_t),
                    real_t (*r)(real_t))
{
  edoFunction->p = p;
  edoFunction->q = q;
  edoFunction->r = r;
}