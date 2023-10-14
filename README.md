# arpackpp (ARPACK++)

## Introduction

Arpackpp is a C++ interface to the ARPACK Fortran package, which
implements the implicit restarted Arnoldi method for iteratively solving
large-scale sparse eigenvalue problems.

Arpackpp is a collection of classes (C++ headers and examples) that
offers C++ programmers an interface to ARPACK. Furthermore, it interfaces
with LAPACK, SuperLU, Cholmod and UMFPACK to incorporate efficient
matrix solvers. Arpackpp preserves the full capability, performance,
accuracy and low memory requirements of the ARPACK Fortran package, but
takes advantage of the C++ object-oriented programming environment.

This GitHub project is designed to provide a common maintained version
of arpackpp. It is derived from the original package (ARPACK++ Version
1.2. by Gomes and Sorensen), which has not been actively maintained for
many years. Several updates have been included in this version (some of
them were previously hosted as patches at
http://reuter.mit.edu/software/arpackpatch/ ). This GitHub repository is
designed to collect fixes and updates (e.g. to more recent or future
releases of the involved libraries). Please consider contributing (see
todo list below).


## Features:

Features of original ARPACK++ package:

- Friendly interface that hides the complicated reverse communication
  interface of the Fortran Arpack package from the user.
- Easy interface using  matrices and vectors via the Standard Template
  Library (STL). 
- Provides an interface between ARPACK and solvers in SuperLU, LAPACK,
  UMFPACK and CHOLMOD to solve eigenvalue problems (specifically shift
  invert methods). 
- Use of templates for optimal performance.

Additional features of this GitHub arpackpp package:

- CMake support for building the examples
- Examples build on Mac OSX using CLang 
- Install scripts for getting and building the dependencies of examples 
- Support for [SuperLU](https://github.com/xiaoyeli/superlu) versions 5.0 and up
- Added initial support for CHOLMOD (for symmetric real problems) 
- Updated UMFPACK sym integration with [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)
- Fixed ARPACK++1.2 to run with g++ 4.4.6 and SuperLU 4.3 (patch see
  here: http://reuter.mit.edu/software/arpackpatch/ )


## TODO

- CMake: get rid of globbing and specify individual files, also add some testing
- UMFPACK complex examples do not build (need update like sym) 
- CHOLMOD complex examples not included (implement similar to real sym)
- Update documentation (install) to cover more scenarios (APT, Homebrew)


## Files

Files included in the main directory:

1. `README.md`

   This file.

2. `INSTALL.md`

   Compile and install notes.

3. `Makefile.inc` (historic)

   An include file used to compile arpackpp examples. You must change
   some directories and machine-dependent directives contained into
   this file prior to compiling the examples. See the description of
   the "makefiles" directory below.

4. `CmakeLists.txt`

   A Cmake file to compile arpackpp examples.

5. `install-*.sh`

   Shell scripts to download and install dependencies into a local
   ./external directory. Some dependencies can also be installed via
   a package-manager on your system. See [INSTALL.md](INSTALL.md)
   for details.

Arpackpp subdirectories:

1. `makefiles` (historic)

   This directory contains example Makefile.inc include files for
   some platforms. Choose one and copy it onto the
   arpackpp/Makefile.inc file.

2. `include`

   The directory that contains arpackpp library, i.e., all header
   files that define arpackpp class templates.

3. `examples`

   The directory where all arpackpp examples can be found. These
   examples are intended to illustrate how to use and compile
   arpackpp classes and are divided according to the type of problem
   being solved and also the kind of information that the user is
   supposed to supply. Look at the [examples/README](examples/README)
   file for further information.

   Note: additional header files are contained in examples/matrices
   and examples/matprod that are needed to build examples or your own
   code!

4. `doc`

   The directory that contains a the arpackpp user's manual and some
   instructions on how to install the libraries required by arpackpp.


## Dependencies

- BLAS/LAPACK, e.g. OpenBLAS https://github.com/xianyi/OpenBLAS

- ARPACK (arpack-ng) https://github.com/opencollab/arpack-ng

For efficient sparse matrix operations, any of these: 

- SuperLU https://github.com/xiaoyeli/superlu

- UMFPACK, CHOLMOD https://github.com/DrTimothyAldenDavis/SuiteSparse

Detailed description of dependencies:

1. ARPACK (Fortran):

   Arpackpp is a C++ interface to ARPACK Fortran code, so the original
   ARPACK library must be installed prior to using the C++ version. A
   maintained package (arpack new generation) can be obtained via the following
   GitHub repository (see also [install-arpack-ng.sh](install-arpack-ng.sh)):

   https://github.com/opencollab/arpack-ng

2. BLAS and LAPACK (Fortran versions):

   Some arpackpp examples require routines from BLAS and LAPACK, so these
   libraries need to be installed before compiling the examples.

   It is recommended that vendor-optimized versions of BLAS and LAPACK are
   installed using a package manager. To install from source, a good choice
   is OpenBLAS or FlexiBLAS. Follow their install instructions on

   https://github.com/xianyi/OpenBLAS
   
   or
   
   https://github.com/mpimd-csc/flexiblas

3. SuperLU:

   Some ARPACK++ classes call SuperLU library functions to solve
   eigenvalue problems that require complex or real (non)symmetric
   matrix decompositions. Thus, SuperLU must also be installed if you
   intend to use one of these classes. SuperLU is via the following
   GitHub repository webpage (see also [install-superlu.sh](install-superlu.sh)):

   https://github.com/xiaoyeli/superlu

4. UMFPACK, CHOLMOD:

   The UMFPACK package can also be used to solve eigenvalue problems
   that require real or complex (non)symmetric/non-Hermitian matrix
   decompositions.

   The CHOLMOD package is performing a Cholesky decomposition and some
   of the symmetric problems can now interface with it.
   
   Both UMFPACK and CHOLMOD are part of the SuiteSparse package which
   is available via the following GitHub repository (see also
   [install-suitesparse.sh](install-suitesparse.sh)):

   https://github.com/DrTimothyAldenDavis/SuiteSparse


## Documentation

   Arpackpp user's manual is available in the doc directory. It contains
   all information needed to declare and solve eigenvalue problems using
   arpack++ classes and functions. Arpackpp computational modes and data
   types are also described in the manual. Instructions on how to
   install the above mentioned libraries are given in the INSTALL.md
   file. Moreover, README files were include in many arpackpp
   directories to give additional information about arpackpp files and
   examples.


## Using arpackpp:

   As a collection of class templates, arpackpp need not to be compiled.
   Because templates are defined in header (.h) files, no object (.o) or
   library (.a) files have to be build, except those corresponding to
   other libraries required by arpackpp (see Dependencies above).
   Arpackpp header files are included in the "include" directory and can
   be moved to another directory if desired. An option in the form

   ```
   -I$(ARPACKPP_INC) \
   -I$(ARPACKPP_INC)/../examples/matrices \
   -I$(ARPACKPP_INC)/../examples/matprod
   ```

   should be added to the command line when compiling programs that use
   arpackpp. Here, ARPACKPP_INC is the name of the directory that
   contains all arpackpp header files. Note, depending on what type of
   problem you want so solve, you need to also include the example
   matrices and/or matprod directories (see examples). You can also use
   `cmake` (see below) with `make install` to install all headers to your
   system into a single directory.


## Compiling and running arpackpp examples:

   Arpackpp supports cmake for the compilation of the examples. To build
   all examples, including the ones that depend on SuperLU, do

   ```
   $ cmake -B build -D ENABLE_SUPERLU=ON
   $ cmake --build build
   ```

   For this to work all dependencies need to be installed (either on the
   system or in the external subdirectory). See [INSTALL.md](INSTALL.md)
   for details. Regular Makefiles (in-source build) are also still supported.


## Acknowledgements

ARPACK++ authors:

-  Francisco M. Gomes (chico AT ime.unicamp.br)
-  Danny Sorensen (LASTNAME AT caam.rice.edu)

arpackpp (2.0.0 and above) authors:

-  Martin Reuter (LASTNAME AT mit.edu)
