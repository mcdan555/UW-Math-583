! quadrature2, modified by Daniel McInally for AMATH 583 for HW5


module quadrature3

    use omp_lib

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

real(kind=8) function trapezoid(f, a, b, n)

    ! Estimate the integral of f(x) from a to b using the
    ! Trapezoid Rule with n points.

    ! Input:
    !   f:  the function to integrate
    !   a:  left endpoint
    !   b:  right endpoint
    !   n:  number of points to use
    ! Returns:
    !   the estimate of the integral
     
    implicit none
    real(kind=8), intent(in) :: a,b
    real(kind=8), external :: f
    integer, intent(in) :: n

    ! Local variables:
    integer :: j
    real(kind=8) :: h, trap_sum, xj

    h = (b-a)/(n-1)
    trap_sum = 0.5d0*(f(a) + f(b))  ! endpoint contributions
    

    do j=2,n-1
        xj = a + (j-1)*h
        trap_sum = trap_sum + f(xj)
        enddo

    trapezoid = h * trap_sum

end function trapezoid

real(kind=8) function simpson(f, a, b, n)

	! Estimates the integral of f(x) from a to b using the
	! Simpson's rule
	
	implicit none
	real(kind=8), intent(in) :: a,b
	real(kind=8), external :: f
	integer, intent(in) :: n
	
	real(kind=8) :: h
	real(kind=8), dimension(n) :: xj,fj,xc,fc
	integer :: j
	
	h=(b-a)/(n-1)
	xj=linspace(a,b,n)
	
	!$omp parallel do private(xj) 
	do j=1,n
		fj(j)=f(xj(j))
		enddo
	
	xc = linspace(a+h/2,b-h/2,n-1)
	
	!$omp parallel do private(xc)
	do j=1,n
		fc(j)=f(xc(j))
		enddo
		
	simpson= (h/6.0) * (2.0*sum(fj) - (fj(1)+fj(n)) + 4.0*sum(fc))
	
end function simpson 
	
	

subroutine error_table(f,a,b,nvals,int_true,method)

    ! Compute and print out a table of errors when the quadrature
    ! rule specified by the input function method is applied for
    ! each value of n in the array nvals.

    implicit none
    real(kind=8), intent(in) :: a,b, int_true
    real(kind=8), external :: f, method
    integer, dimension(:), intent(in) :: nvals

    ! Local variables:
    integer :: j, n
    real(kind=8) :: ratio, last_error, error, int_approx

    print *, "      n         approximation        error       ratio"
    last_error = 0.d0   
    
    !$omp parallel do private(last_error,error,int_approx,n) schedule(dynamic)

    do j=size(nvals),1,-1
        n = nvals(j)
        int_approx = method(f,a,b,n)
        error = abs(int_approx - int_true)
        ratio = last_error / error
        last_error = error  ! for next n

        print 11, n, int_approx, error, ratio
 11     format(i8, es22.14, es13.3, es13.3)
        enddo

end subroutine error_table


end module quadrature3

