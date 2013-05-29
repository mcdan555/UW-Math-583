# $UWHPSC/codes/fortran/newton/Makefile

OBJECTS = functions.o quadrature2.o test2.o
OBJECTS2 = functions.o quadrature3.o test3.o
MODULES = functions.mod quadrature2.mod 
MODULES2 = functions.mod quadrature3.mod 

FFLAGS = -fopenmp -g

.PHONY: test2 test3 clean 

test2: test2.exe
	./test2.exe

test2.exe: $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test2.exe
	
test3: test3.exe
	./test3.exe

test3.exe: $(MODULES2) $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test3.exe


%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod

