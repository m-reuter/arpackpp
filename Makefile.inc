# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# This file contains some definitions used to compile arpack++ examples
# with the g++ compiler under Sun Solaris.


# Defining the machine.

PLAT         = linux

# Defining the compiler.

CPP          = g++

# Defining ARPACK++ directories.
# ARPACKPP_INC is the directory that contains all arpack++ header files.
# SUPERLU_DIR and UMFPACK_DIR must be set to ARPACKPP_INC.

#ARPACKPP_DIR = $(HOME)/arpack++
ARPACKPP_DIR = ../../..
ARPACKPP_INC = $(ARPACKPP_DIR)/include
SUPERLU_DIR  = $(ARPACKPP_INC)
UMFPACK_DIR  = $(ARPACKPP_INC)

# Defining ARPACK, LAPACK, UMFPACK, SUPERLU, BLAS and FORTRAN libraries.
# See the arpack++ manual or the README file for directions on how to 
# obtain arpack, umfpack and SuperLU packages. 
# UMFPACK_LIB and SUPERLU_LIB must be declared only if umfpack and superlu 
# are going to be used. Some BLAS and LAPACK fortran routines are 
# distributed along with arpack fortran code, but the user should verify 
# if optimized versions of these libraries are available before installing 
# arpack. The fortran libraries described below are those required to link
# fortran and c++ code using gnu g++ and f77 compiler under linux.
# Other libraries should be defined if the user intends to compile
# arpack++ on another environment.

ARPACK_LIB   = -larpack
LAPACK_LIB   = -llapack
UMFPACK_LIB  = #-lumfpack
SUPERLU_LIB  = -lsuperlu
BLAS_LIB     = -lblas
FORTRAN_LIBS = -lgfortran

# Defining g++ flags and directories.

# CPP_WARNINGS = -fpermissive 
CPP_WARNINGS = -Wall -ansi -pedantic-errors
CPP_DEBUG    = -g
CPP_OPTIM    = -O
CPP_LIBS     = 
CPP_INC      = 

CPP_FLAGS    = $(CPP_DEBUG) -D$(PLAT) -I$(ARPACKPP_INC) -I$(CPP_INC) \
               $(CPP_WARNINGS)

# Putting all libraries together.

ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) \
               $(BLAS_LIB) $(LAPACK_LIB) $(FORTRAN_LIBS) 

# defining paths.

vpath %.h  $(ARPACK_INC)

