! $MYHPSC/homework4/quarature.f90

! Daniel McInally AMATH 583, HW4


module quadrature

	! module parameters
	implicit none
	
contains

function linspace (a,b,n)

	! This function will generate the the linear space
	
	integer, intent(in) :: n
	real(kind=8) ::linspace(n),ls(n)
	real(kind=8) :: value,a,b
	integer :: j
	
	value = (b-a)/(n-1)
	do j=1,n
		ls(j)=a+(j-1)*value
		enddo
		
	linspace=ls
	
	
end function linspace
	
	
! ----------------------------------------------------------------------------


real(kind=8) function trapezoid (f,a,b,n)

	!Approximated the numerical intergral of the supplied function
	! Input:
	!	f:  The function
	!	a:  The lower limit of integration
	!	b:  The upper limit of integration
	!	n:	
	
	implicit none
	real(kind=8), external :: f
	integer, intent(in) :: n
	real(kind=8), intent(in) :: a,b
	
	!local variables
	real(kind=8) :: base
	integer :: j,iters
	real(kind=8), dimension(n) :: xj,fj
	
	
	base = (b-a) / (n-1)
	xj = linspace(a,b,n)
	
	do j=1,n
		fj(j)=f(xj(j))
		enddo
	
	trapezoid = base*sum(fj) - 0.5*base* (fj(1) +fj(n))
	
	
end function trapezoid

! -------------------------------------------------------------------------------------


subroutine error_table (f,a,b,nvals,int_true)

	implicit none
	real(kind=8), external :: f
	integer, dimension(:), intent(in) :: nvals
	real(kind=8), intent(in) :: int_true,a,b
	
	integer :: j,n
	real(kind=8) :: error, ratio, last_error, int_trap
	
	
	n=size(nvals)
	last_error = 0
	print *, "    n         trapezoid            error       ratio"
	do j=1,size(nvals)
		n=nvals(j)
		int_trap = trapezoid(f,a,b,n)
		error = abs(int_trap - int_true)
		ratio = last_error / error
		last_error = error
		print 11, n, int_trap, error, ratio
11		format(i8, es22.14, es13.3, es13.3)
		enddo
end subroutine error_table
	
	
end module quadrature
	

