#ifndef __UTILS_H__
#define __UTILS_H__

#define MAX_ITERATION 50
#define DEFAULT_TOLERANCE 1.0e-6

typedef double real_t;
typedef unsigned int uint;

typedef struct triDiagonal_SL
{
  real_t *mainDiagonal, *lowerDiagonal, *upperDiagonal, *independetTerms;
  uint n;
} triDiagonal_SL;

typedef struct
{
  uint n;
  real_t a, b;
  real_t ya, yb;
  real_t (*p)(real_t), (*q)(real_t), (*r)(real_t);
} Edo_EQ;

typedef struct
{
  uint n[3];
  real_t a, b;
  real_t ya, yb;
} Assignment_Objective;

real_t *allocArray(uint size);
void setAssignmentObjective(Assignment_Objective *array);
void printArray(real_t *array, uint size);

#endif
