# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/band/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining bandmat directory.

BANDMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: bnsymreg bnsymshf bnsymgre bnsymgsh bnsymgsc bsvd

# compiling and linking each nonsymmetric problem.

bnsymreg:       bnsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bnsymreg bnsymreg.o $(ALL_LIBS)

bnsymshf:	bnsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bnsymshf bnsymshf.o $(ALL_LIBS)

bnsymgre:	bnsymgre.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bnsymgre bnsymgre.o $(ALL_LIBS)

bnsymgsh:	bnsymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bnsymgsh bnsymgsh.o $(ALL_LIBS)

bnsymgsc:	bnsymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bnsymgsc bnsymgsc.o $(ALL_LIBS)

bsvd:	bsvd.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bsvd bsvd.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core bnsymreg bnsymshf bnsymgre bnsymgsh bnsymgsc bsvd

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -c $<


