# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/simple directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# compiling and linking all examples.

all: symsimp

# compiling and linking a simple arpack++ problem.

symsimp: 	symsimp.o
	$(CPP) $(CPP_FLAGS) -o symsimp symsimp.o $(ALL_LIBS) 

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core symsimp

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -c $<

