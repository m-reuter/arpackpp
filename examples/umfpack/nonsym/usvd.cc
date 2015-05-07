/*
   ARPACK++ v1.2/20/2000
   c++ interface to ARPACK code.

   MODULE USVD.cc.
   Example program that illustrates how to determine the condition
   number of a matrix using arpack++ to find its largest and smallest
   singular values.

   1) Problem description:

      In this example, Arpack++ is called to solve the symmetric problem:

                             (A'*A)*v = sigma*v

      where A is an m by n real matrix.
      This formulation is appropriate when m >= n.
      The roles of A and A' must be reversed in the case that m < n.

   2) Data structure used to represent the matrix:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      --------------------------------------------
      lnmatrxv.h       RectangularMatrix, a function that generates
                       matrix A in CSC format.
      arunsmat.h       The ARumNonSymMatrix class definition.
      arssym.h         The ARSymStdEig class definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arssym.h"
#include "lnmatrxv.h"
#include "arunsmat.h"
#include <cmath>


int main()
{

  // Defining variables;

  int     m;          // Number of rows in A.
  int     n;          // Number of columns in A.
  int     nnz;        // Number of nonzero elements in A.
  int*    irow;       // pointer to an array that stores the row
                      // indices of the nonzeros in A.
  int*    pcol;       // pointer to an array of pointers to the
                      // beginning of each column of A in valA.
  double* valA;       // pointer to an array that stores the
                      // nonzero elements of A.
  double  cond;       // Condition number of A.
  double* svalue = new double[6];

  // Creating a rectangular matrix with m = 200 and n = 100.

  n = 100;
  RectangularMatrix(n, m, nnz, valA, irow, pcol);

  // Using ARumNonSymMatrix to store matrix information and to
  // perform the product A'Ax (LU decomposition is not used).

  ARumNonSymMatrix<double, double> A(m, n, nnz, valA, irow, pcol);

  // Defining what we need: eigenvalues from both ends of the spectrum.

  ARSymStdEig<double, ARumNonSymMatrix<double, double> >
    dprob(n, 6L, &A, &ARumNonSymMatrix<double, double>::MultMtMv, "BE", 20L);

  // Finding eigenvalues.

  dprob.Eigenvalues(svalue);

  // Calculating singular values.

  svalue[0] = sqrt(svalue[0]);
  svalue[5] = sqrt(svalue[5]);

  // Obtaining the condition number.

  cond = svalue[5]/svalue[0];

  // Printing some information about the problem.

  std::cout << std::endl << "Testing ARPACK++ class ARSymStdEig" << std::endl;
  std::cout << "Obtaining singular values by solving (A'*A)*v = sigma*v" << std::endl;
  std::cout << "  greatest singular value: " << svalue[5] << std::endl;
  std::cout << "  smallest singular value: " << svalue[0] << std::endl;
  std::cout << "  condition number of A  : " << cond << std::endl;
  std::cout << "MATLAB solution:" << std::endl;
  std::cout << "  greatest singular value:  9.89757224207690 \n";
  std::cout << "  smallest singular value:  1.41683937261247 \n";
  std::cout << "  condition number of A  :  6.98566995906319 \n";

} // main.

