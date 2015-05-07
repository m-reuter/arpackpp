# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining search directories.

EXMP_INC = $(ARPACKPP_DIR)/examples/matprod 
EXC_INC  = $(ARPACKPP_DIR)/examples/matprod/complex

# compiling and linking all examples.

all: compreg compshf compgreg compgshf

# compiling and linking each complex problem.

compreg:	compreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o compreg compreg.o $(ALL_LIBS)

compshf:	compshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o compshf compshf.o $(ALL_LIBS)

compgreg:       compgreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o compgreg compgreg.o $(ALL_LIBS)

compgshf:       compgshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -o compgshf compgshf.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core compgreg compgshf compreg compshf

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXC_INC) -c $<

