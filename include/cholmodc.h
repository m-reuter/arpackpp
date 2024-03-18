/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE CHOLMODc.h.
   Interface to CHOLMOD routines.

   Author of this class:
      Martin Reuter
      Date 11/05/2012
      
   Arpack++ Author:
      Francisco Gomes
      
   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CHOLMODC_H
#define CHOLMODC_H

#include "arcomp.h"
#include "arerror.h"
#include "cholmod.h"

inline cholmod_sparse* CholmodCreateSparse_impl(int m, int n, int nnz,
  void* a, int* irow, int* pcol, char uplo, int itype, int xtype, int dtype)
{
  cholmod_sparse* A = (cholmod_sparse*)malloc(sizeof(cholmod_sparse));

  if (!A) {
    throw new ArpackError(ArpackError::INSUFICIENT_MEMORY,
        "CholmodCreateSparse_impl");
  }

  A->nrow = m;
  A->ncol = n;
  A->nzmax = nnz;
  A->p = pcol;
  A->i = irow;
  A->nz = NULL;
  A->x = a;
  A->z = NULL;
  A->stype = (uplo == 'L' ? -1 : (uplo == 'U' ? 1 : 0));
  A->itype = itype;
  A->xtype = xtype;
  A->dtype = dtype;
  A->sorted = 0;
  A->packed = 1;

  return A;
}

inline cholmod_dense* CholmodCreateDense_impl(int m, int n, void* a, int xtype, int dtype)
{
  cholmod_dense* A = (cholmod_dense*)malloc(sizeof(cholmod_dense));

  if (!A) {
      throw new ArpackError(ArpackError::INSUFICIENT_MEMORY,
          "CholmodCreateDense_impl");
  }

  A->nrow = m;
  A->ncol = n;
  A->nzmax = m * n;
  A->d = m;
  A->x = a;
  A->z = NULL;
  A->xtype = xtype;
  A->dtype = dtype;

  return A;
}

/* CholmodCreateSparse */

#if CHOLMOD_MAIN_VERSION >= 5

inline cholmod_sparse* CholmodCreateSparse(int m, int n, int nnz,
  float* a, int* irow, int* pcol, char uplo)
{
  return CholmodCreateSparse_impl(m, n, nnz, a, irow, pcol, uplo, CHOLMOD_INT, CHOLMOD_REAL, CHOLMOD_SINGLE);
}

inline cholmod_sparse* CholmodCreateSparse(int m, int n, int nnz,
  arcomplex<float>* a, int* irow, int* pcol, char uplo)
{
  return CholmodCreateSparse_impl(m, n, nnz, a, irow, pcol, uplo, CHOLMOD_INT, CHOLMOD_COMPLEX, CHOLMOD_SINGLE);
}

#endif

inline cholmod_sparse* CholmodCreateSparse(int m, int n, int nnz,
  double* a, int* irow, int* pcol, char uplo)
{
  return CholmodCreateSparse_impl(m, n, nnz, a, irow, pcol, uplo, CHOLMOD_INT, CHOLMOD_REAL, CHOLMOD_DOUBLE);
}

inline cholmod_sparse* CholmodCreateSparse(int m, int n, int nnz,
  arcomplex<double>* a, int* irow, int* pcol, char uplo)
{
  return CholmodCreateSparse_impl(m, n, nnz, a, irow, pcol, uplo, CHOLMOD_INT, CHOLMOD_COMPLEX, CHOLMOD_DOUBLE);
}

/* CholmodCreateDense */

#if CHOLMOD_MAIN_VERSION >= 5

inline cholmod_dense* CholmodCreateDense(int m, int n, float* a)
{
  return CholmodCreateDense_impl(m, n, a, CHOLMOD_REAL, CHOLMOD_SINGLE);
}

inline cholmod_dense* CholmodCreateDense(int m, int n, arcomplex<float>* a)
{
  return CholmodCreateDense_impl(m, n, a, CHOLMOD_COMPLEX, CHOLMOD_SINGLE);
}

#endif

inline cholmod_dense* CholmodCreateDense(int m, int n, double* a)
{
  return CholmodCreateDense_impl(m, n, a, CHOLMOD_REAL, CHOLMOD_DOUBLE);
}

inline cholmod_dense* CholmodCreateDense(int m, int n, arcomplex<double>* a)
{
  return CholmodCreateDense_impl(m, n, a, CHOLMOD_COMPLEX, CHOLMOD_DOUBLE);
}

/* CholmodGetDenseData */

template <typename ARTYPE> inline void CholmodGetDenseData(cholmod_dense* A, int n, ARTYPE* a)
{
  memcpy(a, A->x, n * sizeof(ARTYPE));
}

/* CholmodAdd */

/* cholmod_add does not support xtype CHOLMOD_COMPLEX, so those are not provided. */

inline cholmod_sparse* CholmodAdd(cholmod_sparse* A, float sigma, cholmod_sparse* B, cholmod_common* c)
{
  if (sigma == 0.f) return A;
  double alpha[2] = { 1.0, 0.0 };
  double beta[2] = { sigma, 0.0 };
  return cholmod_add(A, B, alpha, beta, 1, 0, c);
}

inline cholmod_sparse* CholmodAdd(cholmod_sparse* A, double sigma, cholmod_sparse* B, cholmod_common* c)
{
  if (sigma == 0.0) return A;
  double alpha[2] = { 1.0, 0.0 };
  double beta[2] = { sigma, 0.0 };
  return cholmod_add(A, B, alpha, beta, 1, 0, c);
}

#endif // CHOLMODC_H
