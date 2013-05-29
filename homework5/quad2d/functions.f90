! Modified by Daniel McInally for AMATH 583 HW5


module functions

    use omp_lib
    implicit none
    integer :: fevals(0:7),gevals(0:7)
    real(kind=8) :: k
    save

contains


	real(kind=8) function g(x,y)
		implicit none
		real(kind=8), intent(in) :: x,y
		integer thread_num
		
		thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        gevals(thread_num) = gevals(thread_num) + 1
		
		g = sin(x+y)
		
	end function g
		

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        integer thread_num,ny,j
        real(kind=8) yj,h, trap_sum,c,d

        ! keep track of number of function evaluations by
        ! each thread:
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        fevals(thread_num) = fevals(thread_num) + 1
    
    c=1.0
    d=4.0
    ny=1000 
    h = (d-c)/(ny-1)
    trap_sum = 0.5d0*(g(x,c) + g(x,d))  ! endpoint contributions
    
    do j=2,ny-1
        yj = c + (j-1)*h
        trap_sum = trap_sum + g(yj,x)
        enddo

    f = h * trap_sum
        
    end function f

end module functions
