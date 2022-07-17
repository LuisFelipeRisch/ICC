#ifndef __SOLVE_H__
#define __SOLVE_H__

#include "utils.h"
#include "edoHandler.h"

typedef struct
{
  real_t *mainDiagonal, *lowerDiagonal, *upperDiagonal, *independetTerms;
  uint n;
} triDiagonal_SL;

#endif