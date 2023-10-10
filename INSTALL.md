# arpackpp installation

The arpackpp library consists of header files and can be
installed without compiling. However, to compile the tests /
examples or a program that includes these headers, it is
necessary to install some libraries first.

## Install Headers Only

Installing the headers into the default include/arpackpp
directory can be done via `cmake`, for example:

```
$ cmake -B build -D ENABLE_TESTS=OFF
$ cmake --install build --prefix /path/to/install
```

Compiling the examples can be done without installing the
headers, but requires a few libraries ...

## System Libraries: GFORTRAN, BLAS, LAPACK, ARPACK

These libraries can be installed via a package manager,
for example when using APT:

```
$ sudo apt update
$ sudo apt install -y gfortran libopenblas-dev libarpack2-dev
```

In case you want to use SuperLU:

```
$ sudo apt install -y libsuperlu-dev
```

In case you want to use UMPACK/CHOLMOD:

```
$ sudo apt install -y libsuitesparse-dev
```

## Install Scripts:

All dependencies can be installed locally using one of the
provided `install-*.sh` shell scripts. The scripts require
`curl` or `wget` to download the source code archives from
the official Github repositories. By default, the scripts will
install static libraries into a subdirectory named _external_.

In case the dependencies should be installed globally (i.e.
/usr/local), the `-g` switch (or `--global-install`) can be
provided (this requires root access). Shared libraries can be
built using the `-s` switch (or `--shared-libs`).

To cleanup temporary build files and the CMake cache, the `-c`
switch (or `--cleanup`) can be used.

If compilers other than the system default should be used, specify
them before calling the install scripts:

```
$ export CC=/path/to/c-compiler
$ export FC=/path/to/fortran-compiler
```

### BLAS, LAPACK:

It is recommended that the BLAS provider is installed using
a package manager. A reasonable choice is OpenBLAS, which is
available for most Linux distributions, MacOS and Windows.

Most classes defined by arpackpp do not require BLAS and LAPACK
(classes for band and dense matrices are the only exception).
However, many examples included in the arpackpp "examples"
directory require some routines from these two packages.

In case a local install is required (for testing), OpenBLAS can
be obtained from GitHub: https://github.com/xianyi/OpenBLAS

The script

```
$ ./install-openblas.sh
```

will install OpenBLAS into the _external_ directory.

### ARPACK:

The actively maintained "new generation" package from GitHub
https://github.com/opencollab/arpack-ng
can be installed via

```
$ ./install-arpack-ng.sh
```
into the _external_ directory.

### SuperLU:

The SuperLU package from GitHub
https://github.com/xiaoyeli/superlu
can be installed via

```
$ ./install-superlu.sh
```
into the _external_ directory.

### UMFPACK and CHOLMOD:

These libraries are part of the SuiteSparse package
https://github.com/DrTimothyAldenDavis/SuiteSparse .
The script

```
$ ./install-suitesparse.sh
```
installs them into the _external_ directory.

## Remarks

Note that you should use the same BLAS for compiling the dependencies,
that you will use when compiling the arpackpp examples (or your own code).
If OpenBLAS is installed locally using the shell script, all other install
scripts will link against that local library by default.

Arpackpp uses the official
[FindBLAS](https://cmake.org/cmake/help/latest/module/FindBLAS.html) and
[FindLAPACK](https://cmake.org/cmake/help/latest/module/FindLAPACK.html)
modules by CMake to find these packages. The libraries found will be displayed
during the configuration. If you find that the wrong BLAS/LAPACK libraries were
chosen or you get an error that a library isn't found, CMake allows specifying
the libraries on the command line

```
$ cmake -B build -D ENABLE_SUPERLU=ON \
        -D BLAS_LIBRARIES=/path/to/your-blas-lib.a \
        -D LAPACK_LIBRARIES=/path/to/your-lapack-lib.a 
```

In case you have multiple BLAS providers installed, you can try enabling
a specific one using the `BLA_VENDOR` option:

```
$ cmake -B build -D ENABLE_SUPERLU=ON -D BLA_VENDOR=OpenBLAS
```

## Compile Examples (cmake):

Arpackpp supports `cmake` for the compilation of the examples. To build
some examples, including the ones that depend on SuperLU, do

```
$ cmake -B build -D ENABLE_SUPERLU=ON
$ cmake --build build --parallel
```

For this to work, all dependencies need to be installed (either on the
system or in the _external_ subdirectory). See above for details.

Compilation of CHOLMOD and UMFPACK examples can be switched-on via:

```
$ cmake -B build -D ENABLE_CHOLMOD=ON -D ENABLE_UMFPACK=ON
```

You can also use `ccmake` instead of `cmake` to see all variables and
manually overwrite specific paths to ensure the right libraries
are being used.

To compile a specific target only, for example "symsimp"
(which can be found in the examples/product/simple directory),
just write

```
$ cmake --build build --target symsimp
```

## Compile Examples (Makefiles in-source build):

Currently we still support standard Makefiles and in-source build:

Arpackpp example directories contain Makefiles that should be used
to compile the examples. For example, to compile example "symsimp",
you just need to write

```
$ make symsimp
```

File symsimp.cc will be compiled and linked to arpackpp libraries,
and an executable file named symsimp will be created.


## Compiler-dependent instructions

These compiler-dependent instructions were supplied originally, it
is unclear whether they are still relevant:

Some compiler-dependent functions and data types used by ARPACK++ were
grouped in the file include/arch.h. Thus, this file should be changed
to reflect the characteristics of your system. Because at the present
time the library was only compiled with the GNU g++ compiler and
tested in a SUN SparcStation, further work must be done in order to
allow the use of ARPACK++ in other environments.

Moreover, ARPACK++ also includes a file, include/arcomp,h, that contains
the definition of a class template called arcomplex, created to emulate
the g++ complex class when other compilers are being used. arcomplex is
the only complex type referenced by other ARPACK++ files, so you must
change the definition of this class in order to work with complex
numbers if g++ (or CC) is not being used.
