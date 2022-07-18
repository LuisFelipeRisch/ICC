#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "edoHandler.h"
#include "solve.h"

int main(int argc, char const *argv[])
{
  Edo_EQ *edoEquation;
  triDiagonal_SL *triDiagonalSL;
  edoEquation = allocEdoEquation(5);
  triDiagonalSL = allocTriDiagonal(5);
  printf("hello word!");
  freeEdoEquation(edoEquation);
  freeTriDiagonal(triDiagonalSL);
  return 0;
}
