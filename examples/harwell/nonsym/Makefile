# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/harwell/nonsym directory makefile

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: hnsymstd hnsymgen

# compiling and linking each nonsymmetric problem.

hnsymstd:       hnsymstd.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hnsymstd hnsymstd.o $(SUPERLU_LIB) $(ALL_LIBS)

hnsymgen:	hnsymgen.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hnsymgen hnsymgen.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core hnsymstd hnsymgen

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<

