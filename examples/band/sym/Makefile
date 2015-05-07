# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/band/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining bandmat directory.

BANDMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: bsymreg bsymshf bsymgreg bsymgshf bsymgbkl bsymgcay

# compiling and linking each symmetric problem.

bsymreg:	bsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymreg bsymreg.o $(ALL_LIBS)

bsymshf:	bsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymshf bsymshf.o $(ALL_LIBS)

bsymgreg:	bsymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymgreg bsymgreg.o $(ALL_LIBS)

bsymgshf:	bsymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymgshf bsymgshf.o $(ALL_LIBS)

bsymgbkl:	bsymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymgbkl bsymgbkl.o $(ALL_LIBS)

bsymgcay:	bsymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsymgcay bsymgcay.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core bsymreg bsymshf bsymgreg bsymgshf bsymgbkl bsymgcay

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -c $<


