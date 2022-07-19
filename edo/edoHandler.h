#ifndef __EDO_HANDLER_H__
#define __EDO_HANDLER_H__

#include "solve.h"
#include "utils.h"

Edo_EQ *allocEdoEquation();
void freeEdoEquation(Edo_EQ *edoEquation);
void initEdoEquation(Edo_EQ *edoFunction,
                     uint n,
                     uint a,
                     uint b,
                     uint ya,
                     uint yb,
                     real_t (*p)(real_t),
                     real_t (*q)(real_t),
                     real_t (*r)(real_t));
void buildTriDiagonalSL(Edo_EQ *edoEquation,
                        triDiagonal_SL *triDiagonalSL);

#endif