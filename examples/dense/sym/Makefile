# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/dense/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining densemat directory.

DENSEMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: dsymreg dsymshf dsymgreg dsymgshf dsymgbkl dsymgcay

# compiling and linking each symmetric problem.

dsymreg:	dsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymreg dsymreg.o $(ALL_LIBS)

dsymshf:	dsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymshf dsymshf.o $(ALL_LIBS)

dsymgreg:	dsymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymgreg dsymgreg.o $(ALL_LIBS)

dsymgshf:	dsymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymgshf dsymgshf.o $(ALL_LIBS)

dsymgbkl:	dsymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymgbkl dsymgbkl.o $(ALL_LIBS)

dsymgcay:	dsymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsymgcay dsymgcay.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core dsymreg dsymshf dsymgreg dsymgshf dsymgbkl dsymgcay

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -c $<


