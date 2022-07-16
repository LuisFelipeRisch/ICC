#include <stdio.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(202201);
  
  SistLinear_t *SL;
  SistLinear_t *originalSL;

  real_t *x, *residueGS, *residueREF;
  double *tTotal;
  int sizes[] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000}; 
  int sizesTam = 9; 
  int qntGSiterations, qntREFiterations;
  double t_egp, t_gs, t_ref; 

  printf("==============================================\n");
  printf("TIPO SL: diagDominate\n");
  printf("TIPO DE ALOCAÇÃO: pontPont\n");
  printf("==============================================\n");

  for (int i = 0; i < sizesTam; i++)
  {
    printf("\n\n================================== N: %d ==================================\n\n", sizes[i]);
    SL = alocaSisLin(sizes[i], pontPont);
    originalSL = alocaSisLin(sizes[i], pontPont);
    iniSisLin(SL, diagDominante, COEF_MAX);
    copySL(originalSL, SL);
    
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueGS = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueREF = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    setZeroVet(x, SL->n); 

    eliminacaoGauss(SL, x, tTotal);
    t_egp = *tTotal; 
    
    copySL(SL, originalSL); 
    setZeroVet(x, SL->n); 

    qntGSiterations = gaussSeidel(SL, x, ERRO, tTotal); 
    t_gs = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueGS);

    qntREFiterations = refinamento(SL, x, ERRO, tTotal); 
    t_ref = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueREF);

    printf("t_egp: %10g; t_gs: %10g; it_gs: %d; normaResiduo_gs: %10g; t_ref: %10g; it_ref: %d; normaResiduo_ref: %10g", t_egp, t_gs, qntGSiterations, normaL2Residuo(residueGS, SL->n),t_ref, qntREFiterations, normaL2Residuo(residueREF, SL->n));
    liberaSisLin(SL); 
    liberaSisLin(originalSL); 
    free(x); 
    free(residueGS); 
    free(residueREF); 
    free(tTotal);
  }

  printf("\n==============================================\n");
  printf("TIPO SL: diagDominate\n");
  printf("TIPO DE ALOCAÇÃO: pontVet\n");
  printf("==============================================\n");

  for (int i = 0; i < sizesTam; i++)
  {
    printf("\n\n================================== N: %d ==================================\n\n", sizes[i]);
    SL = alocaSisLin(sizes[i], pontVet);
    originalSL = alocaSisLin(sizes[i], pontVet);
    iniSisLin(SL, diagDominante, COEF_MAX);
    copySL(originalSL, SL);
    
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueGS = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueREF = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    setZeroVet(x, SL->n); 

    eliminacaoGauss(SL, x, tTotal);
    t_egp = *tTotal; 
    
    copySL(SL, originalSL); 
    setZeroVet(x, SL->n); 

    qntGSiterations = gaussSeidel(SL, x, ERRO, tTotal); 
    t_gs = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueGS);

    qntREFiterations = refinamento(SL, x, ERRO, tTotal); 
    t_ref = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueREF);

    printf("t_egp: %10g; t_gs: %10g; it_gs: %d; normaResiduo_gs: %10g; t_ref: %10g; it_ref: %d; normaResiduo_ref: %10g", t_egp, t_gs, qntGSiterations, normaL2Residuo(residueGS, SL->n),t_ref, qntREFiterations, normaL2Residuo(residueREF, SL->n));
    liberaSisLin(SL); 
    liberaSisLin(originalSL); 
    free(x); 
    free(residueGS); 
    free(residueREF); 
    free(tTotal);
  }

  printf("\n==============================================\n");
  printf("TIPO SL: HILBERT\n");
  printf("TIPO DE ALOCAÇÃO: pontPont\n");
  printf("==============================================\n");

  for (int i = 0; i < sizesTam; i++)
  {
    printf("\n\n================================== N: %d ==================================\n\n", sizes[i]);
    SL = alocaSisLin(sizes[i], pontPont);
    originalSL = alocaSisLin(sizes[i], pontPont);
    iniSisLin(SL, hilbert, COEF_MAX);
    copySL(originalSL, SL);
    
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueGS = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueREF = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    setZeroVet(x, SL->n); 

    eliminacaoGauss(SL, x, tTotal);
    t_egp = *tTotal; 
    
    copySL(SL, originalSL); 
    setZeroVet(x, SL->n); 

    qntGSiterations = gaussSeidel(SL, x, ERRO, tTotal); 
    t_gs = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueGS);

    qntREFiterations = refinamento(SL, x, ERRO, tTotal); 
    t_ref = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueREF);

    printf("t_egp: %10g; t_gs: %10g; it_gs: %d; normaResiduo_gs: %10g; t_ref: %10g; it_ref: %d; normaResiduo_ref: %10g", t_egp, t_gs, qntGSiterations, normaL2Residuo(residueGS, SL->n),t_ref, qntREFiterations, normaL2Residuo(residueREF, SL->n));
    liberaSisLin(SL); 
    liberaSisLin(originalSL); 
    free(x); 
    free(residueGS); 
    free(residueREF); 
    free(tTotal);
  }

  printf("\n==============================================\n");
  printf("TIPO SL: HILBERT\n");
  printf("TIPO DE ALOCAÇÃO: pontVet\n");
  printf("==============================================\n");

  for (int i = 0; i < sizesTam; i++)
  {
    printf("\n\n================================== N: %d ==================================\n\n", sizes[i]);
    SL = alocaSisLin(sizes[i], pontVet);
    originalSL = alocaSisLin(sizes[i], pontVet);
    iniSisLin(SL, hilbert, COEF_MAX);
    copySL(originalSL, SL);
    
    x = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueGS = (real_t *) malloc(sizes[i] * sizeof(real_t));
    residueREF = (real_t *) malloc(sizes[i] * sizeof(real_t));
    tTotal = (double *) malloc(sizeof(double));
    setZeroVet(x, SL->n); 

    eliminacaoGauss(SL, x, tTotal);
    t_egp = *tTotal; 
    
    copySL(SL, originalSL); 
    setZeroVet(x, SL->n); 

    qntGSiterations = gaussSeidel(SL, x, ERRO, tTotal); 
    t_gs = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueGS);

    qntREFiterations = refinamento(SL, x, ERRO, tTotal); 
    t_ref = *tTotal; 
    copySL(SL, originalSL); 
    calculateResidue(SL, x, residueREF);

    printf("t_egp: %10g; t_gs: %10g; it_gs: %d; normaResiduo_gs: %10g; t_ref: %10g; it_ref: %d; normaResiduo_ref: %10g", t_egp, t_gs, qntGSiterations, normaL2Residuo(residueGS, SL->n),t_ref, qntREFiterations, normaL2Residuo(residueREF, SL->n));
    liberaSisLin(SL); 
    liberaSisLin(originalSL); 
    free(x); 
    free(residueGS); 
    free(residueREF); 
    free(tTotal);
  }
  

  // código do programa aqui
  
}

