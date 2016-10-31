# Change Log

## upcomming

* Update install script to SuperLU 5.2.1
* Explicitly link with quadmath for Ubuntu build (fix #6)
* Add build targets: install and examples, thanks to [jschueller]
* Remove default build (header library does not need build)


## arpackpp - 2.3.0

* Add CMake support for building the examples
* Add install scripts for getting and building the dependencies of examples
* Update README.md, add INSTALL.md and CHANGES.md notes, clean up doc
* Add BSD 3-clause LICENSE thanks to original authors:
  - F.A.M. Gomes (UNICAMP)
  - D.C. Sorensen (Rice University)

-- Martin Reuter - Nov 22 2015


## arpackpp - 2.2.0

* Update interface for SuperLU 5.0 (thread-safe) thanks to [wo80]

-- Martin Reuter - Oct 14 2015


## arpackpp - 2.1.1

* Fix SuperLU 4.3 interface thanks to [wo80]

-- Martin Reuter - Oct 11 2015


## arpackpp - 2.1.0

* Add interface for SuiteSparse CHOLMOD (for real symmmetric only)
* Update UMFPACK interface to work with SuiteSparse UMFPACK (symetric only)

-- Martin Reuter - May 7 2015


## arpackpp - 2.0.0

* Move arpackpp to GitHub
* Support g++ 4.4.6
  - Use of this-> pointer for members of template classes
  - Includes corrected (e.g. math.h -> cmath)
  - Namespace std:: added
  - Replaced char * with const std::string& in most places
* Updated SuperLU interface to run with SuperLU 4.3
* From patch at http://reuter.mit.edu/software/arpackpatch/ 

-- Martin Reuter - May 7 2015


## arpackpp - 1.2.0

* ARPACK++ version 1.2. Feb 20, 2000
  - Authors: Francisco M Gomes and Danny Sorensen
  - http://www.ime.unicamp.br/%7echico/arpack++/arpack++.tar.gz
