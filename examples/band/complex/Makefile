# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/band/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining bandmat directory.

BANDMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: bcompreg bcompshf bcompgre bcompgsh

# compiling and linking each complex problem.

bcompreg:       bcompreg.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bcompreg bcompreg.o $(ALL_LIBS)

bcompshf:	bcompshf.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bcompshf bcompshf.o $(ALL_LIBS)

bcompgre:	bcompgre.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bcompgre bcompgre.o $(ALL_LIBS)

bcompgsh:	bcompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -o bcompgsh bcompgsh.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core bcompreg bcompshf bcompgre bcompgsh

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(BANDMAT_DIR) -c $<


