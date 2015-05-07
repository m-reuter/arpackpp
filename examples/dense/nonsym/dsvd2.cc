/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DSVD2.cc.
   Example program that illustrates how to determine the largest 
   singular values of an out-of-core dense matrix using arpack++.

   1) Problem description:

      In this example, Arpack++ is called to solve the symmetric problem:

                             (A'*A)*v = sigma*v

      where A is an m by n real matrix.
      This formulation is appropriate when m >= n.
      The roles of A and A' must be reversed in the case that m < n.

   2) Data structure used to represent the matrix:

      Because m >> n in this example, A is supposed to be stored 
      row-wise in a file. Since this program is intended to show 
      how to find the singular values of an out-of-core matrix, A 
      will be loaded by blocks. This means that only part of the 
      matrix will be stored in vector called A, and the contents
      of this vector will be changed several times while performing
      the matrix-vector product y <- Ax.

   3) Included header files:

      File             Contents
      -----------      --------------------------------------------
      ardnsmat.h       The ARdsNonSymMatrix class definition.
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
#include "ardnsmat.h"
#include <iostream>
#include <cmath>


int main()
{

  // Defining variables;

  int     i;
  double* svalue = new double[4];

  // Using ARdsNonSymMatrix to store matrix data and to
  // perform the product A'Ax (LU decomposition is not used).

  ARdsNonSymMatrix<double, double> A("matrix.dat", 250);

  // Defining what we need: eigenvalues with largest magnitude.

  ARSymStdEig<double, ARdsNonSymMatrix<double, double> >
    dprob(A.ncols(), 4L, &A, &ARdsNonSymMatrix<double, double>::MultMtMv);

  // Finding eigenvalues.

  dprob.Eigenvalues(svalue);

  // Calculating singular values.

  for (i=0; i<4; i++) {
    svalue[i] = sqrt(svalue[i]);
  }

  // Printing some information about the problem.

  std::cout << std::endl << "Testing ARPACK++ class ARSymStdEig" << std::endl;
  std::cout << "Obtaining singular values by solving (A'*A)*v = sigma*v" << std::endl;
  std::cout << std::endl << "greatest singular values: " << std::endl; 
  for (i=0; i<4; i++) {
    std::cout << "  sigma [" << i+1 << "]: " << svalue[i] << std::endl;
  }

} // main.
