#include <stdio.h>
#include <stdlib.h>
#include "solve.h"
#include "utils.h"
#include "edoHandler.h"
#include "coefficients.h"

#define ASSIGNMENT_SIZE 3
#define N_SIZE 3

int main(int argc, char const *argv[])
{
  Assignment_Objective assignmentObjective[ASSIGNMENT_SIZE];
  setAssignmentObjective(assignmentObjective);

  for (uint i = 0; i < ASSIGNMENT_SIZE; i++)
  {
    for (uint j = 0; j < N_SIZE; j++)
    {
      uint size = assignmentObjective[i].n[j];

      Edo_EQ *edoEquation = allocEdoEquation();
      triDiagonal_SL *triDiagonalSL = allocTriDiagonal(size);
      real_t *solution = allocArray(size);

      switch (i)
      {
      case 0:
        initEdoEquation(edoEquation,
                        size,
                        assignmentObjective[i].a,
                        assignmentObjective[i].b,
                        assignmentObjective[i].ya,
                        assignmentObjective[i].yb,
                        p1,
                        q1,
                        r1);
        break;
      case 1:
        initEdoEquation(edoEquation,
                        size,
                        assignmentObjective[i].a,
                        assignmentObjective[i].b,
                        assignmentObjective[i].ya,
                        assignmentObjective[i].yb,
                        p2,
                        q2,
                        r2);
        break;
      case 2:
        initEdoEquation(edoEquation,
                        size,
                        assignmentObjective[i].a,
                        assignmentObjective[i].b,
                        assignmentObjective[i].ya,
                        assignmentObjective[i].yb,
                        p3,
                        q3,
                        r3);
        break;
      default:
        break;
      }

      buildTriDiagonalSL(edoEquation, triDiagonalSL);
      gaussSeidel(triDiagonalSL, solution, DEFAULT_TOLERANCE);

      printTriDiagonalMatrix(triDiagonalSL);
      printf("SOLUTION: \n");
      printArray(solution, size);

      printf("\n\n ======================================== \n\n");

      freeTriDiagonal(triDiagonalSL);
      freeEdoEquation(edoEquation);
      free(solution);
    }
  }
  return 0;
}
