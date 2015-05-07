# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/umfpack/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining CSCMAT directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: usymreg usymshf usymgreg usymgshf usymgbkl usymgcay

# compiling and linking each symmetric problem.

usymreg:	usymreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymreg usymreg.o $(UMFPACK_LIB) $(ALL_LIBS)

usymshf:	usymshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymshf usymshf.o $(UMFPACK_LIB) $(ALL_LIBS)

usymgreg:	usymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymgreg usymgreg.o $(UMFPACK_LIB) $(ALL_LIBS)

usymgshf:	usymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymgshf usymgshf.o $(UMFPACK_LIB) $(ALL_LIBS)

usymgbkl:	usymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymgbkl usymgbkl.o $(UMFPACK_LIB) $(ALL_LIBS)

usymgcay:	usymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usymgcay usymgcay.o $(UMFPACK_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core usymreg usymshf usymgreg usymgshf usymgbkl usymgcay

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(UMFPACK_DIR) -c $<


