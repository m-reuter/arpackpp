/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SMatrixD.h
   Class template for the 1-dimensional mass matrix
   on the interval [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXD_H
#define SMATRIXD_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class ART>
class SymMatrixD: public MatrixWithProduct<ART> {

 private:

  ART  *Ad, *Adl, *Adu, *Adu2;
  int  *ipiv;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorM();

  void SolveM(ART* v);

  void MultMv(ART* v, ART* w);

  SymMatrixD(int nv);

  virtual ~SymMatrixD();

}; // SymMatrixD.


template<class ART>
inline void SymMatrixD<ART>::FactorDataDeallocate()
// Eliminates the data structure used on matrix factorization.

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class ART>
void SymMatrixD<ART>::FactorM()
// Factors M.

{

  int  i, ierr;
  ART  h, r1, r2;

  const ART one  = 1.0;
  const ART four = 4.0;
  const ART six  = 6.0;

  if (decsize != ncols()) {
    decsize = ncols();
    FactorDataDeallocate();
    Ad   = new ART[ncols()];
    Adl  = new ART[ncols()];
    Adu  = new ART[ncols()];
    Adu2 = new ART[ncols()];
    ipiv = new int[ncols()];
  }

  h  = one/ART(ncols()+1);
  r2 = h/six;
  r1 = r2*four;

  for (i=0; i<ncols(); i++) {
    Ad[i]  = r1;
    Adl[i] = r2;
  }

  copy(ncols(), Adl, 1, Adu, 1);
  gttrf(ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorM.


template<class ART>
inline void SymMatrixD<ART>::SolveM(ART* v)
// Solves M*w = v. v is overwritten with vector w.

{

  int  ierr;
  char *type = "N";

  gttrs(type, ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, v, ncols(), ierr);

} // SolveM.

template<class ART>
void SymMatrixD<ART>::MultMv(ART* v, ART* w)
//  Performs w <- M*v.

{

  int  j;
  ART  h;

  const ART one  = 1.0;
  const ART four = 4.0;
  const ART six  = 6.0;

  w[0] = four*v[0] + v[1];
  for (j=1; j<ncols()-1; j++) {
    w[j] = v[j-1] + four*v[j] + v[j+1];
  }
  w[ncols()-1] = v[ncols()-2] + four*v[ncols()-1];

  // Scaling the vector w by h.

  h = one / (ART(ncols()+1)*six);
  scal(ncols(), h, w, 1L);

  return;

} //  MultMv.


template<class ART>
inline SymMatrixD<ART>:: SymMatrixD(int nval): MatrixWithProduct<ART>(nval)
// Constructor.

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;

} // Constructor.


template<class ART>
inline SymMatrixD<ART>::~SymMatrixD()
// Destructor.

{

  FactorDataDeallocate();

} // Destructor.


#endif // SMATRIXD_H

