! $UWHPSC/codes/fortran/newton/functions.f90

! Modified by Daniel McInally for AMATH 583, HW 3

module functions

contains

real(kind=8) function f_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt = x**2 - 4.d0

end function f_sqrt


real(kind=8) function fprime_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt = 2.d0 * x

end function fprime_sqrt


real(kind=8) function f(x)
    implicit none
    real(kind=8), intent(in) :: x
    real, parameter :: Pi =3.1415927
	
    
    f = x*cos(Pi*x) - 1.0 + 0.6*x**2
    
end function f


real(kind=8) function fp(x)
    implicit none
    real(kind=8), intent(in) :: x
    real, parameter :: Pi =3.1415927
    
    fp = cos(Pi*x) - Pi*x*sin(Pi*x) + 1.2 * x
    
end function fp
    

end module functions
