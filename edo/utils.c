#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "solve.h"

real_t *allocArray(uint size)
{
  return (real_t *)calloc(size, sizeof(real_t));
}

void setAssignmentObjective(Assignment_Objective *array)
{
  array[0].n[0] = 5;
  array[0].n[1] = 10;
  array[0].n[2] = 100;
  array[0].a = 0.0;
  array[0].b = 12.0;
  array[0].ya = 0.0;
  array[0].yb = 0.0;

  array[1].n[0] = 5;
  array[1].n[1] = 10;
  array[1].n[2] = 100;
  array[1].a = 0.0;
  array[1].b = 1.0;
  array[1].ya = 0.0;
  array[1].yb = 1.0;

  array[2].n[0] = 5;
  array[2].n[1] = 10;
  array[2].n[2] = 100;
  array[2].a = 0.0;
  array[2].b = 1.0;
  array[2].ya = -1.0;
  array[2].yb = 0.0;
}

void printArray(real_t *array, uint size)
{
  for (uint i = 0; i < size; i++)
    printf("%.15g ", array[i]);
}