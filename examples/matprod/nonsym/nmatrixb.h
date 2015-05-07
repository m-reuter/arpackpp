/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixB.h
   Class template for the tridiagonal matrix derived from
   the standard central difference discretization of the 1-d 
   convection diffusion operator u" + rho*u' on the interval 
   [0, 1] with zero Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXB_H
#define NMATRIXB_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class NonSymMatrixB: public MatrixWithProduct<T> {

 private:

  T    rho;
  T    shift;
  T    *Ad, *Adl, *Adu, *Adu2;
  int  *ipiv;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorOP();

  void MultMv(T* v, T* w);

  void MultOPv(T* v, T* w);

  NonSymMatrixB(int nv, T rhov);

  NonSymMatrixB(int nv, T shiftv, T rhov);

  virtual ~NonSymMatrixB();

}; // NonSymMatrixB.


template<class T>
inline void NonSymMatrixB<T>::FactorDataDeallocate()
// Eliminates the data structure used on matrix factorization.

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void NonSymMatrixB<T>::FactorOP()
/*
  Factors (M-shift*I).
*/

{

  int  j, ierr;
  T    h, s, s1, s2, s3;

  const T one = 1.0;
  const T two = 2.0;

  if (decsize != ncols()) {
    decsize = ncols();
    FactorDataDeallocate();
    Ad   = new T[ncols()];
    Adl  = new T[ncols()];
    Adu  = new T[ncols()];
    Adu2 = new T[ncols()];
    ipiv = new int[ncols()];
  }

  h  = one/T(ncols()+1);
  s  = rho*h/two;
  s1 = -one - s;
  s2 = two - shift;
  s3 = -one + s;

  for (j=0; j<ncols()-1; j++) {
    Adl[j] = s1;
    Ad[j]  = s2;
    Adu[j] = s3;
  }
  Ad[ncols()-1]  = s2;

  gttrf(ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void NonSymMatrixB<T>::MultMv(T* v, T* w)
/*
  Performs the matrix-vector multiplication w <- A*v.
*/

{

  int  j;
  T    dd, dl, du, s, h;

  const T one = 1.0;
  const T two = 2.0;

  h  = one/T(ncols()+1);
  s  = rho*h/two;
  dd = two;
  dl = -one - s;
  du = -one + s;

  w[0] = dd*v[0] + du*v[1];
  for (j=1; j<ncols()-1; j++) {
    w[j] = dl*v[j-1] + dd*v[j] + du*v[j+1];
  }
  w[ncols()-1] = dl*v[ncols()-2] + dd*v[ncols()-1];

  return;

} //  MultMv.


template<class T>
void NonSymMatrixB<T>::MultOPv(T* v, T* w)
/*
  Performs the matrix-vector product w <- inv(M-shift*I)*v.
*/

{

  int  ierr;
  char *type = "N";

  copy(ncols(), v, 1, w, 1);
  gttrs(type, ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, w, ncols(), ierr);

} // MultOPv


template<class T>
inline NonSymMatrixB<T>::
NonSymMatrixB(int nval, T rhov): MatrixWithProduct<T>(nval)
/*
  Constructor
*/

{

  rho     = rhov;
  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = 0.0;

} // Constructor.


template<class T>
inline NonSymMatrixB<T>::
NonSymMatrixB(int nv, T shiftv, T rhov): MatrixWithProduct<T>(nv)
/*
  Constructor with shift.
*/

{

  rho     = rhov;
  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = shiftv;
  FactorOP();

} // Constructor with shift.


template<class T>
inline NonSymMatrixB<T>::~NonSymMatrixB()
/*
  Destructor.
*/
{

  FactorDataDeallocate();

} // Destructor.


#endif // NMATRIXB_H

