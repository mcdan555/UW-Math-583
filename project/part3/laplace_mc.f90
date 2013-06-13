! laplace_mc.f90   created by Daniel McInally for AMATH 583, Final Project Part 3



program laplace_mc
	
	use problem_description
	use mc_walk
	use random_util, only:init_random_seed
	
	
	
	implicit none
	
	
	! Local Variables
	integer :: seed, n_mc, max_steps, n_success, k, n_total
	real(kind=8) :: u_true,u_mc,i0,j0, error, u_mc_total, u_sum_old,u_sum_new,x0,y0
	
	

	! Set initial values
	
	x0 = 0.9
	y0 = 0.6
	
	i0 = nint((x0-ax)/dx)
	j0 = nint((y0-ay)/dy)
	
	x0 = ax + i0*dx
	y0 = ay + j0*dy
	
	u_true = utrue(x0,y0)
	print '("True Value is",es17.9)',u_true
	
	
	
	! Initialize random generator
	seed = 12345
    call init_random_seed(seed)
    
    
    ! max steps
    max_steps =100*max(nx, ny)
	!max_steps =10
    
    ! intial Monte-Carlo 
	n_mc = 10
	nwalks = 0
	open(unit=25, file='mc_laplace_error.txt', status='unknown')
	
	call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
	
	error = abs((u_mc-u_true)/u_true)
	
	print '(" ",i8," ",e22.15," ",e13.6)' &
			, n_success,u_mc,error
	write(25,'(i10,e23.15,e15.6)') n_success, u_mc, error
			
	! MC loop
	u_mc_total = u_mc
	n_total = n_success
	
	do k = 1,12
		u_sum_old = u_mc_total * n_total
		
		call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
		u_sum_new = u_mc * n_success
		n_total = n_total + n_success
		u_mc_total = (u_sum_old + u_sum_new) / n_total
	
		error = abs((u_mc_total - u_true)/u_true)
	
		print '(" ",i8," ",e22.15," ",e13.6)' &
				, n_total,u_mc_total,error
		write(25,'(i10,e23.15,e15.6)') n_total, u_mc_total, error
	
		n_mc=2*n_mc
		
		enddo
	
	print '("Final approximation to u(x0,y0):  ",es21.14)',u_mc_total
	print '("Total number of random walks:  ",i9)',nwalks
	

end program laplace_mc
