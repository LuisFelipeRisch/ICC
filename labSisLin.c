#include <stdio.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(time(NULL));
  
  SistLinear_t *SL;
  real_t *x;
  double *tTotal;
  int sizes[] = {10}; 
  int sizesTam = 1; 

  for (int i = 0; i < sizesTam; i++)
  {
    // SL = lerSisLin(pontPont);
    SL = alocaSisLin(sizes[i], pontPont);
    iniSisLin(SL, diagDominante, COEF_MAX);
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    // gaussSeidel(SL, x, ERRO, tTotal);
    refinamento(SL, x, ERRO, tTotal);
    
    // printf("O sistema linear que obtivemos da eliminação de gauss foi: \n"); 
    prnSisLin(SL); 
    // printf("E foi essa resposta que obtivemos em %lf milisegundos\n", *tTotal);
    prnVetor(x, SL->n); 

    liberaSisLin(SL); 
    free(x); 
    free(tTotal);
  }
  

  // código do programa aqui
  
}

