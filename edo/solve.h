#ifndef __SOLVE_H__
#define __SOLVE_H__

#include "utils.h"
#include "edoHandler.h"

triDiagonal_SL *allocTriDiagonal(uint n);
void freeTriDiagonal(triDiagonal_SL *triDiagonalSL);
void printTriDiagonalMatrix(triDiagonal_SL *triDiagonalSL);
void gaussElimination(triDiagonal_SL *triDiagonalSL);
void retroSubstitution(triDiagonal_SL *triDiagonalSL,
                       real_t *solution);
void gaussSeidel(triDiagonal_SL *triDiagonalSL, real_t *solution, real_t tolerance);

#endif