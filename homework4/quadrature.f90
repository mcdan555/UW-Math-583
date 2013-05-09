! $MYHPSC/homework4/quarature.f90

! Daniel McInally AMATH 583, HW4


module quadrature

	! module parameters
	implicit none
	
contains

function linspace (a,b,n)

	! This function will generate the the linear space
	
	integer, intent(in) :: a,b,n
	real(kind=8) ::linspace(n),ls(n)
	real(kind=8) :: value
	integer :: j
	
	value = abs(b-a)/n
	
	do j=1,n
		ls(j)=value
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
	integer, intent(in) :: a,b,n
	
	!local variables
	real(kind=8) :: base
	integer :: j,iters
	real(kind=8), dimension(n) :: xj,fj
	
	
	base = (b-a) / (n-1)
	xj = linspace(a,b,n)
	
	do j=1,n
		fj(j)=f(xj(j))
		enddo
	
	trapezoid = base*sum(fj) - 0.5*base* (fj(a) +fj(b))
	
	
end function trapezoid
	
end module quadrature
	

