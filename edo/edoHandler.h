#ifndef __EDO_HANDLER_H__
#define __EDO_HANDLER_H__

#include "utils.h"

typedef struct
{
  uint n;
  real_t a, b;
  real_t ya, yb;
  real_t (*p)(real_t), (*q)(real_t), (*r)(real_t);
} Edo_EQ;

Edo_EQ *allocEdoEquation(uint n);
void freeEdoEquation(Edo_EQ *edoEquation);

#endif