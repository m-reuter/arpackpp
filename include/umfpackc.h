/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE UMFPACKc.h.
   Interface to UMFPACK routines.

   Author of this class:
      Martin Reuter
      Date 2/28/2013
      
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

#ifndef UMFPACKC_H
#define UMFPACKC_H

#include "arcomp.h"
#include "arerror.h"
#include <umfpack.h>

/* umfpack_defaults */

template <typename T> inline void umfpack_defaults(double* Control)
{
    throw ArpackError(ArpackError::NOT_IMPLEMENTED, "umfpack_defaults");
}

template <> inline void umfpack_defaults<double>(double* Control)
{
    umfpack_di_defaults(Control);
}

template <> inline void umfpack_defaults<arcomplex<double>>(double* Control)
{
    umfpack_zi_defaults(Control);
}

/* umfpack_triplet_to_col */

inline int umfpack_triplet_to_col(int32_t n_row, int32_t n_col, int32_t nz,
    const int32_t Ti[], const int32_t Tj[], const double Tx[], int32_t Ap[], int32_t Ai[], double Ax[])
{
    return umfpack_di_triplet_to_col(n_row, n_col, nz, Ti, Tj, Tx, Ap, Ai, Ax, nullptr);
}

inline int umfpack_triplet_to_col(int32_t n_row, int32_t n_col, int32_t nz,
    const int32_t Ti[], const int32_t Tj[], const arcomplex<double> Tx[], int32_t Ap[], int32_t Ai[], arcomplex<double> Ax[])
{
    return umfpack_zi_triplet_to_col(n_row, n_col, nz, Ti, Tj, 
        (double*)(&Tx[0]), nullptr, Ap, Ai, (double*)(&Ax[0]), nullptr, nullptr);
}

/* umfpack_symbolic */

inline int umfpack_symbolic(int32_t n_row, int32_t n_col, int32_t Ap[], int32_t Ai[], double Ax[],
    void** Symbolic, const double* Control, double* Info)
{
    return umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info);
}

inline int umfpack_symbolic(int32_t n_row, int32_t n_col, int32_t Ap[], int32_t Ai[], arcomplex<double> Ax[],
    void** Symbolic, const double* Control, double* Info)
{
    return umfpack_zi_symbolic(n_row, n_col, Ap, Ai, (double*)(&Ax[0]), nullptr, Symbolic, Control, Info);
}

/* umfpack_numeric */

inline int umfpack_numeric(int32_t Ap[], int32_t Ai[], double Ax[],
    void* Symbolic, void** Numeric, const double* Control, double* Info)
{
    return umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
}

inline int umfpack_numeric(int32_t Ap[], int32_t Ai[], arcomplex<double> Ax[],
    void* Symbolic, void** Numeric, const double* Control, double* Info)
{
    return umfpack_zi_numeric(Ap, Ai, (double*)(&Ax[0]), nullptr, Symbolic, Numeric, Control, Info);
}

/* umfpack_solve */

inline int umfpack_solve(int sys, int32_t Ap[], int32_t Ai[], double Ax[],
    double* X, double* B, void* Numeric, const double* Control, double* Info)
{
    return umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
}

inline int umfpack_solve(int sys, int32_t Ap[], int32_t Ai[], arcomplex<double> Ax[],
    arcomplex<double>* X, arcomplex<double>* B, void* Numeric, const double* Control, double* Info)
{
    return umfpack_zi_solve(sys, Ap, Ai, (double*)(&Ax[0]), nullptr, (double*)(&X[0]), nullptr, (double*)(&B[0]), nullptr, Numeric, Control, Info);
}

/* umfpack_free_symbolic */

template <typename T> inline void umfpack_free_symbolic(void** Symbolic)
{
    throw ArpackError(ArpackError::NOT_IMPLEMENTED, "umfpack_free_symbolic");
}

template <> inline void umfpack_free_symbolic<double>(void** Symbolic)
{
    umfpack_di_free_symbolic(Symbolic);
}

template <> inline void umfpack_free_symbolic<arcomplex<double>>(void** Symbolic)
{
    umfpack_zi_free_symbolic(Symbolic);
}

/* umfpack_free_numeric */

template <typename T> inline void umfpack_free_numeric(void** Numeric)
{
    throw ArpackError(ArpackError::NOT_IMPLEMENTED, "umfpack_free_numeric");
}

template <> inline void umfpack_free_numeric<double>(void** Numeric)
{
    umfpack_di_free_numeric(Numeric);
}

template <> inline void umfpack_free_numeric<arcomplex<double>>(void** Numeric)
{
    umfpack_zi_free_numeric(Numeric);
}

#endif // UMFPACKC_H
