#include <stdlib.h>
#include "utils.h"

real_t *allocArray(uint size)
{
  return (real_t *)calloc(size, sizeof(real_t));
}