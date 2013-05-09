# $UWHPSC/codes/fortran/newton/Makefile

OBJECTS = quadrature.o test1.o
OBJECTS2 = quadrature.o test2.o
OBJECTS3 = quadrature_omp.o test2_omp.o
MODULES = quadrature.mod
MODULES2 = quadrature_omp.mod

FFLAGS = -g

.PHONY: test1 test2 test2_omp clean 

test1: test1.exe
	./test1.exe
	
test2: test2.exe
	./test2.exe
	
test2_omp: test2_omp.exe
	./test2_omp.exe

test1.exe: $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test1.exe

test2.exe: $(MODULES) $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test2.exe
	
test2_omp.exe: $(MODULES2) $(OBJECTS3)
	gfortran -fopenmp quadrature_omp.f90 test2_omp.f90 -o test2_omp.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod

