/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARCSPen.h.
   Arpack++ class ARchSymMPencil definition.
   (CHOLMOD wrapper)

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

#ifndef ARCSPEN_H
#define ARCSPEN_H

#include "blas1c.h"
#include "arcsmat.h"


template<class ARTYPE>
class ARchSymPencil
{

 protected:

  ARchSymMatrix<ARTYPE>* A;
  ARchSymMatrix<ARTYPE>* B;
  cholmod_factor *LAsB ; 
  bool    factoredAsB;
  cholmod_common c ;

  virtual void Copy(const ARchSymPencil& other);

 public:

  bool IsFactored() { return factoredAsB; }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp);

  ARchSymPencil() : factoredAsB(false), A(nullptr), B(nullptr), LAsB(nullptr) { cholmod_start(&c); }
  // Short constructor that does nothing.

  ARchSymPencil(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARchSymPencil(const ARchSymPencil& other) { cholmod_start (&c) ; Copy(other); }
  // Copy constructor.

  virtual ~ARchSymPencil() {  if (LAsB) cholmod_free_factor(&LAsB, &c);  cholmod_finish (&c) ;}
  // Destructor.

  ARchSymPencil& operator=(const ARchSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARchSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARchSymPencil<ARTYPE>::Copy(const ARchSymPencil<ARTYPE>& other)
{
  if (LAsB) cholmod_free_factor(&LAsB, &c);
  A        = other.A;
  B        = other.B;
  factoredAsB = other.factoredAsB;
  if (factoredAsB)
    LAsB = cholmod_copy_factor(other.LAsB, &c);

} // Copy.


template<class ARTYPE>
void ARchSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARchSymPencil::FactorAsB");
  }

  if (LAsB) cholmod_free_factor(&LAsB, &c);

  cholmod_sparse* AsB;

  AsB = CholmodAdd(A->A, -sigma, B->A, &c);
    
  LAsB = cholmod_analyze (AsB, &c);
  int info = cholmod_factorize (AsB, LAsB, &c);  

  factoredAsB = (info != 0);  
  if (c.status != CHOLMOD_OK) {

    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymPencil::FactorAsB");
    
    factoredAsB = false;
  }

  if (A->A != AsB) {
    cholmod_free_sparse(&AsB, &c);
  }

} // FactorAsB (ARTYPE shift).


template<class ARTYPE>
void ARchSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  ::copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.

template<class ARTYPE>
void ARchSymPencil<ARTYPE>::MultInvAsBv(ARTYPE* v, ARTYPE* w)
{
  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARchSymPencil::MultInvAsBv");
  }

  // Solving A.w = v (or AsI.w = v).
  
  //create b from v (data is not copied!!)
  cholmod_dense * b = CholmodCreateDense(A->n, 1, v);

  cholmod_dense *x = cholmod_solve (CHOLMOD_A, LAsB, b, &c);

  CholmodGetDenseData(x, A->n, w);

  free(b);
  cholmod_free_dense(&x, &c);

} // MultInvAsBv

template<class ARTYPE>
inline void ARchSymPencil<ARTYPE>::
DefineMatrices(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if (A->n != B->n) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARchSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARchSymPencil<ARTYPE>::
ARchSymPencil(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp)
{
  cholmod_start (&c);
  LAsB = nullptr; 
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARchSymPencil<ARTYPE>& ARchSymPencil<ARTYPE>::
operator=(const ARchSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSPEN_H
