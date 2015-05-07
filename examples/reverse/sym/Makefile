# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/reverse/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining objects.

EXMP_INC = $(ARPACKPP_DIR)/examples/matprod 
EXS_INC  = $(ARPACKPP_DIR)/examples/matprod/sym 

# compiling and linking all examples.

all: rsymreg rsymshf rsymgreg rsymgshf rsymgbkl rsymgcay

# compiling and linking each symmetric problem.

rsymreg: 	rsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymreg rsymreg.o \
            $(ALL_LIBS)

rsymshf:        rsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymshf rsymshf.o \
            $(ALL_LIBS)

rsymgreg:       rsymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymgreg rsymgreg.o \
            $(ALL_LIBS)

rsymgshf:       rsymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymgshf rsymgshf.o \
            $(ALL_LIBS)

rsymgbkl:       rsymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymgbkl rsymgbkl.o \
            $(ALL_LIBS)

rsymgcay:       rsymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o rsymgcay rsymgcay.o \
            $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *.o *~ core rsymgbkl rsymgcay rsymgreg rsymgshf rsymreg rsymshf

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -c $<

