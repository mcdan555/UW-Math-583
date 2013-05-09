
! Created from Test1.f90 (Leveque 2013)
! Modified by Daniel McInally for AMATH 583, HW4


program test2

    use quadrature, only: trapezoid, error_table

    implicit none
    real(kind=8) :: a,b,int_true,k
    integer :: nvals(12), i

    a = 0.d0
    b = 2.d0
    k = 1000.d0
    
    int_true = (b-a) + (b**4 - a**4) / 4.d0 + (cos(k*a) - cos (k*b))/k

    print 10, int_true
 10 format("true integral: ", es22.14)
    print *, " "  ! blank line

    ! values of n to test:
    do i=1,12
        nvals(i) = 5 * 2**(i)
        enddo

    call error_table(f, a, b, nvals, int_true)

contains

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x
        
        f = 1.d0 + x**3 + sin(k*x)
        
    end function f

end program test2
