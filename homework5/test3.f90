! test2.f90 written by Daniel McInally for AMATH 583 HW5
! some code derived from test.f90 by 


program test3

	
	use omp_lib
	use quadrature3, only: trapezoid, error_table
	use functions, only: f, fevals, k
	
	implicit none
	integer :: j, nthreads
	integer :: nvals(12)
	real(kind=8) :: a,b, int_true
	
	real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate

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
	
	call system_clock(tclock1)
	call cpu_time(t1)  ! Time before
		
	call error_table(f,a,b,nvals,int_true,trapezoid)
	
	call cpu_time(t2)
	call system_clock(tclock2, clock_rate)
	
	elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
	print *, " "
	print 11, elapsed_time
11	format("Elapsed time = ",f12.8, " seconds")

    print 12, t2-t1
12 	format("CPU time = ",f12.8, " seconds")


    do j=0,nthreads-1
        print 101,  j, fevals(j)
101     format("fevals by thread ",i2,": ",i13)
        enddo

    print 102, sum(fevals)
102 format("Total number of fevals: ",i10)
	
end program test3
		
	

	

