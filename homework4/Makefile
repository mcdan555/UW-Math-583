# $UWHPSC/codes/fortran/newton/Makefile

OBJECTS = quadrature.o test1.o
OBJECTS2 = quadrature.o test2.o
MODULES = quadrature.mod

FFLAGS = -g

.PHONY: test1 test2 clean 

test1: test1.exe
	./test1.exe
	
test2: test2.exe
	./test2.exe

test1.exe: $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test1.exe

test2.exe: $(MODULES) $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test2.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod

