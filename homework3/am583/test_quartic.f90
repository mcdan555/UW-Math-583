! $MYHPSC/homework3/am583/test_quartic.f90
! Created by Daniel McInally for AMATH 583, HW3

program test_quartic

	use newton, only: solve,tol
	use functions, only: f_quartic, fprime_quartic,epsilon

! 	Variable declaration
	
	implicit none
	integer :: iters,j,k
	logical :: debug
	real, dimension(3) :: testep,testtol
	real(kind=8) :: x,x0,xstar,fx,error

!	Initial test values
	
	debug = .false.
	testep=(/0.100D-03,0.100D-07,0.100D-11/)
	testtol=(/0.100D-04,0.100D-09,0.100D-13/)
	x0=0.4000000000000D+01

!   Initial output
	
	print *, 'Starting with initial guess ',x0
	print *
	print *, '    epsilon        tol    iters          x                 f(x)        x-xstar'
	
!	Double loop to calculate approximated values for different values of epsilon and tol
	
	do j=1,3
		do k=1,3
			
			epsilon=testep(j)
			xstar=1+epsilon**(1/4.0)
			tol=testtol(k)
			call solve(f_quartic, fprime_quartic, x0, x, iters, debug)
			error = x-xstar
			fx=f_quartic(x)
	
			print 11, epsilon, tol, iters, x, fx, error
11  		format(2d13.3, i4, d24.15, 2d13.3)

			enddo
		print *
		enddo
		
		
		
end program test_quartic
