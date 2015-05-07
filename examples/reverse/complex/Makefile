# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/reverse/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining objects.

EXMP_INC = $(ARPACKPP_DIR)/examples/matprod 
EXC_INC  = $(ARPACKPP_DIR)/examples/matprod/complex

# compiling and linking all examples.

all: rcompreg rcompshf rcompgre rcompgsh

# compiling and linking each complex problem.

rcompreg:       rcompreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o rcompreg rcompreg.o \
            $(ALL_LIBS)

rcompshf:       rcompshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o rcompshf rcompshf.o \
            $(ALL_LIBS)

rcompgre:       rcompgre.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o rcompgre rcompgre.o \
            $(ALL_LIBS)

rcompgsh:       rcompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o rcompgsh rcompgsh.o \
            $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *.o *~ core rcompreg rcompshf rcompgre rcompgsh

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -c $<

