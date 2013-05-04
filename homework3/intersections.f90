! $MYHPSC/homework3/intesections.f90
! Created by Daniel McInally for AMATH 569 HW3


program intersections
use newton
use functions

!	This program will use newtons method to compute the intersections of
!	two fuctions whos difference and the derivitive of the difference is given 
!	in f and fp respectivly.

implicit none


integer, parameter :: n=4
integer :: k,iters
real(kind=8), dimension(n) :: x0,x
logical :: debug=.false.
!real(kind=8), external :: f, fp

x0 = (/-2.18026,-1.61399,-0.794267,1.44331/)


do k=1,n
	call solve(f,fp, x0(k), x(k), iters, debug)
	print *
	print *, "solve returns x=" ,x(k), "after", iters, "iterations"
	print *, "the value of f(x) is ", f(x(k))
	enddo

end program intersections


