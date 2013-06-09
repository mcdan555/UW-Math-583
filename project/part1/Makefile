
OBJECTS = functions.o quadrature.o test3.o
FFLAGS = 
NUM_PROCS ?= 4   # default if not specified on command line or env variable

.PHONY: test test2 clean 

test: test3.exe
	mpiexec -n $(NUM_PROCS) test3.exe

test.exe: $(OBJECTS)
	mpif90 $(FFLAGS) $(OBJECTS) -o test3.exe

%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

