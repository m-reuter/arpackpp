/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSPen.h.
   Arpack++ class ARumSymPencil definition.

   Modified to work with Umfpack v5.??
      Martin Reuter
      Date 02/28/2013

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

#ifndef ARUSPEN_H
#define ARUSPEN_H

#include "arch.h"
#include "arerror.h"
#include "arusmat.h"
#include "blas1c.h"


template<class ARTYPE>
class ARumSymPencil
{

 protected:

  ARumSymMatrix<ARTYPE>* A;
  ARumSymMatrix<ARTYPE>* B;
  ARumSymMatrix<ARTYPE>  AsB;

  void* Numeric;

  virtual void Copy(const ARumSymPencil& other);

  void Expand(ARumSymMatrix<ARTYPE>* A);

  void ClearMem();

 public:

  bool IsFactored() { return (Numeric != nullptr); }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp);

  ARumSymPencil(): A(nullptr), B(nullptr), Numeric(nullptr) { }
  // Short constructor that does nothing.

  ARumSymPencil(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARumSymPencil(const ARumSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymPencil() { }
  // Destructor.

  ARumSymPencil& operator=(const ARumSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::Expand(ARumSymMatrix<ARTYPE>* A)
{
    auto mat = A->A;

    if (mat->IsTriangular())
    {
        int n = mat->nrows();
        int ndiag = mat->DiagIndices();
        int nz = 2 * mat->nzeros();

        if (ndiag > 0)
        {
            // Don't count diagonal entries twice.
            nz -= ndiag;
        }

        int* ap = new int[n + 1];
        int* ai = new int[nz];
        ARTYPE* ax = new ARTYPE[nz];

        ARSparseMatrix<ARTYPE> full(n, n, ap, ai, ax, nz);

        int status = mat->Expand(full);

        if (status > 0)
        {
            throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                "ARumSymPencil::Expand");
        }

        // Memory is now owned by A.
        A->DefineMatrix(n, nz, ax, ai, ap, 'S', true, true);
    }

}


template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::ClearMem()
{

  //if (A) { delete A; A = nullptr; }
  //if (B) { delete B; B = nullptr; }

} // ClearMem.

template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::Copy(const ARumSymPencil<ARTYPE>& other)
{
  ClearMem();
  A        = other.A;
  B        = other.B;
  //AsB      = other.AsB;

} // Copy.

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

    // Quitting the function if A and B were not defined.

    if (!(A->IsDefined() && B->IsDefined())) {
        throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymPencil::FactorAsB");
    }

    // Defining matrix AsB.

    if (!AsB.IsDefined()) {

        Expand(A);
        Expand(B);

        int* count = new int[A->n];
        int* work = new int[A->m];

        int nnz = A->A->PrepareAdd(*B->A, count, work);

        delete[] count;
        delete[] work;

        int* ap = new int[A->n + 1];
        int* ai = new int[nnz];
        ARTYPE* ax = new ARTYPE[nnz];

        // Do not validate AsB since though the matrix is allocated, no
        // meaningful values are set.

        AsB.DefineMatrix(A->m, nnz, ax, ai, ap, 'S', false, true);
    }

    // Subtracting sigma*B from A and storing the result on AsB.

    A->A->Add(-sigma, *B->A, *AsB.A);

    // Decomposing AsB.

    void* Symbolic;

    auto ap = AsB.A->pcol();
    auto ai = AsB.A->irow();
    auto ax = AsB.A->values();

    if (umfpack_symbolic(A->n, A->n, ap, ai, ax, &Symbolic, AsB.control, AsB.info) != UMFPACK_OK) {
      throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB");
    }

    if (umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, AsB.control, AsB.info) != UMFPACK_OK) {
      throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB");
    }

    umfpack_free_symbolic<ARTYPE>(&Symbolic);

    AsB.factored = true;
}

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::MultInvAsBv(ARTYPE* v, ARTYPE* w)
{
  if (!Numeric) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARchSymPencil::MultInvAsBv");
  }

  auto ap = AsB.A->pcol();
  auto ai = AsB.A->irow();
  auto ax = AsB.A->values();

  // Solving A.w = v (or AsI.w = v).

  int status = umfpack_solve(UMFPACK_A, ap, ai, ax, w, v, Numeric, AsB.control, AsB.info);

  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::MultInvv");

} // MultInvAsBv

template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::
DefineMatrices(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if (A->n != B->n) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARumSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARumSymPencil<ARTYPE>::
ARumSymPencil(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp)
{
  //AsB.factored  = false;
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARumSymPencil<ARTYPE>& ARumSymPencil<ARTYPE>::
operator=(const ARumSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSPEN_H
