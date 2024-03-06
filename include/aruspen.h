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

#include "arusmat.h"
#include "blas1c.h"


template<class ARTYPE>
class ARumSymPencil
{

 protected:

  ARumSymMatrix<ARTYPE>* A;
  ARumSymMatrix<ARTYPE>* B;
  //ARumSymMatrix<ARTYPE> AsB;
  void*   Numeric;
  int*    Ap;
  int*    Ai;
  ARTYPE* Ax; 

  virtual void Copy(const ARumSymPencil& other);

  void ExpandAsB(ARTYPE sigma);

  void ClearMem();

 public:

  bool IsFactored() { return (Numeric != NULL); }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp);

  ARumSymPencil() { Numeric = NULL; Ap = NULL; Ai = NULL; Ax = NULL; }
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
inline void ARumSymPencil<ARTYPE>::ClearMem()
{

  if (Numeric) umfpack_free_numeric<ARTYPE>(&Numeric);
  if (Ai) delete [] Ai;
  Ai = NULL;
  if (Ap) delete [] Ap;
  Ap = NULL;
  if (Ax) delete [] Ax;
  Ax = NULL;

} // ClearMem.



template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::Copy(const ARumSymPencil<ARTYPE>& other)
{
  ClearMem();
  A        = other.A;
  B        = other.B;

} // Copy.


template<class ARTYPE>
void ARumSymPencil<ARTYPE>::ExpandAsB(ARTYPE sigma)
{

  ClearMem();
 
  int mynnz = 2*A->nnz+2*B->nnz;
  if (sigma == 0.0)
    mynnz = 2*A->nnz;
  
  // create triples (i,j,value)
  int * tripi = new int[mynnz];
  int * tripj = new int[mynnz];
  ARTYPE* tripx = new ARTYPE[mynnz];
  if (tripi == NULL || tripj == NULL || tripx ==NULL)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB out of memory (1)");
  
  int count = 0;
  int i,j;
  for (i=0; i < A->n; i++)
  {
    // create triplets from A
    for (j=A->pcol[i]; j<(A->pcol[i+1]); j++)
    {
      tripi[count] = i;
      tripj[count] = A->irow[j];
      tripx[count] = A->a[j];
      count++;
      if (i != A->irow[j]) // not on diag
      {
        tripj[count] = i;
        tripi[count] = A->irow[j];
        tripx[count] = A->a[j];
        count++;
      }
    }
  
   if (sigma != 0.0)
   {
    // create triplets from -sigma B
    for (j=B->pcol[i]; j<(B->pcol[i+1]); j++)
    {
      tripi[count] = i;
      tripj[count] = B->irow[j];
      tripx[count] = -sigma * B->a[j];
      count++;
      if (i != B->irow[j]) // not on diag
      {
        tripj[count] = i;
        tripi[count] = B->irow[j];
        tripx[count] = tripx[count-1];
        count++;
      }
    }
    }

  }

  // convert triples (A-sigma B) to Ax Ap Ai
  Ap = new int[A->n + 1];
  Ai = new int[count];
  Ax = new ARTYPE[count];
  if (!Ap || !Ai || !Ax )
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB out of memory (2)");
  
  int status = umfpack_triplet_to_col (A->n, A->n, count, tripi, tripj, tripx, Ap, Ai, Ax,  (int *)NULL) ;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB triplet to col");

  // cleanup
  delete [] tripi;
  delete [] tripj;
  delete [] tripx;

}

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumSymPencil::FactorAsB");
  }


  // Subtracting sigma*B from A and storing the result 
  ExpandAsB(sigma);

  // Decomposing AsB.
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  umfpack_defaults<ARTYPE>(Control) ;
  void *Symbolic ;
  int status = umfpack_symbolic (A->n, A->n, Ap, Ai, Ax, &Symbolic, Control, Info) ;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB symbolic");
  status =  umfpack_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ;
  if (status == 1)
  {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric (matrix singular)");
  }
  if (status < UMFPACK_OK)
  {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric");
  }
  umfpack_free_symbolic<ARTYPE>(&Symbolic) ;

} // FactorAsB (ARTYPE shift).


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

  // Solving A.w = v (or AsI.w = v).
  int status = umfpack_solve (UMFPACK_A, Ap, Ai, Ax, w, v, Numeric, NULL, NULL) ;
  if (status == 1)
  {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric (matrix singular)");
  }
  if (status < UMFPACK_OK)
  {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::MultInvAsBv");
 
  }

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
  Numeric = NULL;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;

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
