! test2.f90 written by Daniel McInally for AMATH 583 HW5
! some code derived from test.f90 by 


program test2

	
	use omp_lib
	use quadrature2, only: simpson, error_table
	use functions, only: f, fevals, k
	
	implicit none
	integer :: j, nthreads
	integer :: nvals(12)
	real(kind=8) :: a,b, int_true

	fevals = 0
	k = 1.d3
	a = 0.d0
	b = 2.d0
	int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
	
	
	nthreads = 1 
	!$ nthreads = 4    ! for openmp
	!$ call omp_set_num_threads(nthreads)
	
	
	print 100, nthreads
100 format("This program uses ",i2," threads.")

	print 10, int_true
10	format("True integral: ", es22.14)
	print *, " "
	
	do j=1,12
		nvals(j) = 50 * 2**(j-1)
		enddo
		
	call error_table(f,a,b,nvals,int_true,simpson)
	
end program test2
		
	

	

