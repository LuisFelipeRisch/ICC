#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
  *tTotal = timestamp(); 

  for (int i = 0; i < SL->n - 1; i++){
    int pivot = findPivot(SL, i); 
    if(pivot != i)
      exchangeLines(SL, i, pivot);

    for (int j = i + 1; j < SL->n; j++){
      real_t m = SL->A[j][i] / SL->A[i][i];
      SL->A[j][i] = 0.0;

      for(int k = i + 1; k < SL->n; k++)
        SL->A[j][k] -= SL->A[i][k] * m;
      SL->b[j] -= SL->b[i] * m;
    }
  } 

  retroSubs(SL, x);

  *tTotal = timestamp() - *tTotal; 
  return 0; 
}


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(real_t *residue, int tam)
{ 
  real_t sum = 0.0;

  for (int i = 0; i < tam; i++)
    sum += residue[i] * residue[i];
  
  return sqrt(sum); 
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  *tTotal = timestamp(); 
  int currIteration = 0; 
  real_t result;
  real_t *currIterationX = (real_t *) malloc(SL->n * sizeof(real_t)); 
  real_t *oldIterationX = (real_t *) malloc(SL->n * sizeof(real_t)); 
  real_t *diffIterationX = (real_t *) malloc(SL->n * sizeof(real_t)); 
  setZeroVet(currIterationX, SL->n);
  setZeroVet(oldIterationX, SL->n);
  setZeroVet(diffIterationX, SL->n);
  
  while (currIteration == 1 || (currIteration < MAXIT && findMaxInVet(diffIterationX, SL->n) > erro))
  {
    for (int i = 0; i < SL->n; i++)
    {
      result = SL->b[i];
      
      for (int j = 0; j < SL->n; j++)
        if(i != j)
          result += ( -1 * SL->A[i][j]) * currIterationX[j];

      result /= SL->A[i][i];
      oldIterationX[i] = currIterationX[i]; 
      currIterationX[i] = result;
    }
    absDiffBetweenVets(diffIterationX, currIterationX, oldIterationX, SL->n);
    currIteration++;
  }

  copyVet(x, currIterationX, SL->n);

  free(currIterationX);
  free(oldIterationX);
  free(diffIterationX);
  *tTotal = timestamp() - *tTotal; 
  return currIteration; 
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  *tTotal = timestamp(); 

  real_t *originalB = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t *residue = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t *w = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t l2Norm = 1;
  int counter = 0;

  copyVet(originalB, SL->b, SL->n);
  gaussSeidel(SL, x, erro, tTotal);

  while (counter < MAXIT && l2Norm > erro)
  {
    calculateResidue(SL, x, residue);
    prnVetor(residue, SL->n);
    l2Norm = normaL2Residuo(residue, SL->n);
    copyVet(SL->b, residue, SL->n);
    copyVet(w, x, SL->n);
    gaussSeidel(SL, w, erro, tTotal);
    sumArrays(x, x, w, SL->n); 
    copyVet(SL->b, originalB, SL->n);
    counter++;
  }

  free(originalB);
  free(residue);
  free(w);
  *tTotal = timestamp() - *tTotal; 
  return counter;
}

