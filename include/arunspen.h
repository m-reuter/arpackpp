/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUNSPen.h.
   Arpack++ class ARumNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUNSPEN_H
#define ARUNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "umfpackc.h"
#include "arunsmat.h"


template<class ARTYPE, class ARFLOAT>
class ARumNonSymPencil
{

 protected:

  char                               part;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>* A;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>* B;
  ARumNonSymMatrix<ARTYPE, ARFLOAT>  AsB;
#ifdef ARCOMP_H
  ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> AsBc;
#endif

  virtual void Copy(const ARumNonSymPencil& other);

  void SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
                   int yind[], int ny, ARTYPE z[], int zind[], int& nz);

#ifdef ARCOMP_H
  void SparseSaxpy(arcomplex<ARFLOAT> a, ARFLOAT x[], int xind[], int nx,
                   ARFLOAT y[], int yind[], int ny,
                   arcomplex<ARFLOAT> z[], int zind[], int& nz);
#endif

  void SubtractAsB(ARTYPE sigma);

#ifdef ARCOMP_H
  void SubtractAsB(ARFLOAT sigmaR, ARFLOAT sigmaI);
#endif

 public:

#ifdef ARCOMP_H
  bool IsFactored() { return (AsB.IsFactored() || AsBc.IsFactored()); }
#else
  bool IsFactored() { return AsB.IsFactored(); }
#endif

  bool IsSymmetric() { return AsB.IsSymmetric(); }

  void FactorAsB(ARTYPE sigma);

#ifdef ARCOMP_H
  void FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp = 'R');
#endif

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

#ifdef ARCOMP_H
  void MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w);
#endif

  void MultInvAsBv(ARFLOAT* v, ARFLOAT* w);

  void DefineMatrices(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                      ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp);

  ARumNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARumNonSymPencil(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                   ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp);
  // Long constructor.

  ARumNonSymPencil(const ARumNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymPencil() { }
  // Destructor.

  ARumNonSymPencil& operator=(const ARumNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymPencil<ARTYPE, ARFLOAT>::
Copy(const ARumNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  part     = other.part;
  A        = other.A;
  B        = other.B;
  AsB      = other.AsB;
#ifdef ARCOMP_H
  AsBc     = other.AsBc;
#endif

} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
            int yind[], int ny, ARTYPE z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  // TODO

} // SparseSaxpy (ARTYPE).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SparseSaxpy(arcomplex<ARFLOAT> a, ARFLOAT x[], int xind[], int nx,
            ARFLOAT y[], int yind[], int ny, 
            arcomplex<ARFLOAT> z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  // TODO

} // SparseSaxpy (arcomplex<ARFLOAT>).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::SubtractAsB(ARTYPE sigma)
{

  // TODO

} // SubtractAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
SubtractAsB(ARFLOAT sigmaR, ARFLOAT sigmaI)
{

  // TODO

} // SubtractAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {

    // TODO

  }

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  AsB.FactorA();

} // FactorAsB (ARTYPE shift).


#ifdef ARCOMP_H
template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
FactorAsB(ARFLOAT sigmaR, ARFLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsBc.IsDefined()) {

    // TODO

  }

  // Subtracting sigma*B from A and storing the result on AsBc.

  SubtractAsB(sigmaR, sigmaI);

  // Decomposing AsBc.

  AsBc.FactorA();

} // FactorAsB (arcomplex<ARFLOAT> shift).
#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::
MultInvAsBv(arcomplex<ARFLOAT>* v, arcomplex<ARFLOAT>* w)
{

  AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

} // MultInvAsBv (arcomplex<ARFLOAT>).

#endif // ARCOMP_H.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymPencil<ARTYPE, ARFLOAT>::MultInvAsBv(ARFLOAT* v, ARFLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((ARTYPE*)v,(ARTYPE*)w);

  }
  else {                // shift is complex.

#ifdef ARCOMP_H

    int                i;
    arcomplex<ARFLOAT> *tv, *tw;

    tv = new arcomplex<ARFLOAT>[AsBc.ncols()];
    tw = new arcomplex<ARFLOAT>[AsBc.ncols()];

    for (i=0; i!=AsBc.ncols(); i++) tv[i] = arcomplex<ARFLOAT>(v[i], 0.0);

    AsBc.MultInvv(tv, tw);

    if (part=='I') {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = imag(tw[i]);
    }
    else {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = real(tw[i]);
    }

    delete[] tv;
    delete[] tw;

#endif // ARCOMP_H.

  }

} // MultInvAsBv (ARFLOAT).


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymPencil<ARTYPE, ARFLOAT>::
DefineMatrices(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
               ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARumNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymPencil<ARTYPE, ARFLOAT>::
ARumNonSymPencil(ARumNonSymMatrix<ARTYPE, ARFLOAT>& Ap, 
                 ARumNonSymMatrix<ARTYPE, ARFLOAT>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE, class ARFLOAT>
ARumNonSymPencil<ARTYPE, ARFLOAT>& ARumNonSymPencil<ARTYPE, ARFLOAT>::
operator=(const ARumNonSymPencil<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUNSPEN_H
