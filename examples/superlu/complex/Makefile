# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/superlu/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: lcompreg lcompshf lcompgre lcompgsh

# compiling and linking each complex problem.

lcompreg:	lcompreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lcompreg lcompreg.o $(SUPERLU_LIB) $(ALL_LIBS)

lcompshf:	lcompshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lcompshf lcompshf.o $(SUPERLU_LIB) $(ALL_LIBS)

lcompgre:	lcompgre.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lcompgre lcompgre.o $(SUPERLU_LIB) $(ALL_LIBS)

lcompgsh:	lcompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lcompgsh lcompgsh.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core lcompreg lcompshf lcompgre lcompgsh

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<

