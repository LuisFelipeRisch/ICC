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
  SistLinear_t *originalSL;

  real_t *x;
  double *tTotal;
  int sizes[] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000}; 
  int sizesTam = 9; 
  int qntGSiterations, qntREFiterations;
  int t_egp, t_gs, t_ref; 

  for (int i = 0; i < sizesTam; i++)
  {
    printf("n \t t_egp \t t_gs \t normaResiduo_gs \t t_ref \t it_ref \t normaResiduo_ref");
    SL = alocaSisLin(sizes[i], pontPont);
    originalSL = alocaSisLin(sizes[i], pontPont);
    iniSisLin(SL, diagDominante, COEF_MAX);
    copySL(originalSL, SL);
    
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    setZeroVet(x, SL->n); 

    eliminacaoGauss(SL, x, tTotal);
    t_egp = *tTotal; 
    
    copySL(SL, originalSL); 
    setZeroVet(x, SL->n); 

    qntGSiterations = gaussSeidel(SL, x, ERRO, tTotal); 
    t_gs = *tTotal; 

    qntREFiterations = refinamento(SL, x, ERRO, tTotal); 
    t_ref = *tTotal; 

    liberaSisLin(SL); 
    liberaSisLin(originalSL); 
    free(x); 
    free(tTotal);
    printf("\n\n===========================================================================================\n\n");
  }
  

  // código do programa aqui
  
}

