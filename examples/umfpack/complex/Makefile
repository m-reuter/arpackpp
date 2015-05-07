# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/umfpack/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: ucompreg ucompshf ucompgre ucompgsh

# compiling and linking each complex problem.

ucompreg:	ucompreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o ucompreg ucompreg.o $(UMFPACK_LIB) $(ALL_LIBS)

ucompshf:	ucompshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o ucompshf ucompshf.o $(UMFPACK_LIB) $(ALL_LIBS)

ucompgre:	ucompgre.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o ucompgre ucompgre.o $(UMFPACK_LIB) $(ALL_LIBS)

ucompgsh:	ucompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o ucompgsh ucompgsh.o $(UMFPACK_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core ucompreg ucompshf ucompgre ucompgsh

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(UMFPACK_DIR) -c $<

