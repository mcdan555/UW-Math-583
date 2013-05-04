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

x0 = (/-2.18026,-1.61399,-0.794267,1.44331/)


do k=1,n
	call solve(f,fp, x0(k), x(k), iters, debug)
	print *
	print *, "With initial guess x0= ",x0(k)
	print 11, x(k), iters
11	format('      solve returns x=   ',e22.15,'  after', i3, '  iterations')
	enddo

end program intersections


