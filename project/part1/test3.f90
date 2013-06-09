!	Modified by Daniel McInally for AMATH 583 Final Project - Part1
!   Originally modified as test2.f90 for HW6


program test3

    use mpi

    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

    implicit none
    real(kind=8) :: a,b,int_true, int_approx,dxs,int_sub

    integer :: proc_num, num_procs, ierr, n, fevals_total,j,ns,jj,nerr,lastint,nextint
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(kind=8), dimension(:,2) :: ab_sub
    real(kind=8), allocatable, dimension(:) :: int_subs

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    ! All processes set these values so we don't have to broadcast:
    k = 1.d3   ! functions module variable 
    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
    n = 1000
    nerr = 0

    ! Each process keeps track of number of fevals:
    fevals_proc = 0

	! Initial input and output of true integral
    if (proc_num==0) then
		if (num_procs==1) then
			print *, " " ! Blank line
			print *, "*** Error, This program requires more than one process to operate"
			nerr = 1
		else
    		print *, " "  ! blank line
        	print '("Using ",i3," processes")', num_procs
        	print *, "How many subintervals? "
        	read *, ns
        	print '("true integral: ", es22.14)', int_true
        	print *, " "  ! blank line
        
        	call MPI_BCAST(ns, 1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        	endif
        endif
        	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for process 0 to print	
    	
    ! 	Halting all processes if error
    if (nerr==1) then
    	call MPI_FINALIZE(ierr)
    	stop
    	endif
    
    !  Process 0 (Master) code
   	if (proc_num==0) then 	
    	allocate(int_subs(ns))
    	allocate(ab_sub(ns,2))
    	dxs = (b-a) / ns
    	
    	
    	! Load the sub interval into matrix ab_sub
    	do j=1,ns
    		ab_sub(j,1) = a + (j-1)*dxs
    		ab_sub(j,2) = a + j*dxs
    		enddo	
    		
    	lastint = 0 ! Tracks the last value of the interval index
    		
    	! Send the first batch to worker processes
    	do j=1,min(num_procs-1, ns)
    		call MPI_SEND(ab_sub(j,:), 2, MPI_DOUBLE_PRECISION, j , j, &
    					MPI_COMM_WORLD, ierr)
    		lastint = lastint + 1		
    		enddo
    		
    	! Get results and send more work if needed
    	do j=1,ns
    		call MPI_RECV(int_sub, 1, MPI_DOUBLE_PRECISION, &
    					MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    		sender = status(MPI_SOURCE)
    		jj=status(MPI_TAG)
    		int_subs(jj)=int_sub
    		
    		! Send out more work
    		if (lastint<ns) then
    			nextint = lastint + 1
    			call MPI_SEND(ab_sub(nextint,:), 2, MPI_DOUBLE_PRECISION, sender , nextint, &
    					MPI_COMM_WORLD, ierr)
    			lastint = lastint + 1
    		 else
    		 	! If no more work send tag=0
    		 	call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender , 0, &
    					MPI_COMM_WORLD, ierr)
    		 endif
    		enddo		
        endif

    
    ! Note: The worker processes calculate portions of the integral
    if (proc_num /= 0) then
      
      	call MPI_RECV(ab_sub, 2, MPI_DOUBLE_PRECISION, &
      				0,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      	
      	j = status(MPI_TAG)			
      	
    	int_sub = trapezoid(f,ab_sub(1),ab_sub(2),n)
    	
    	call MPI_SEND(int_sub, 1, MPI_DOUBLE_PRECISION, 0 , j, &
    					MPI_COMM_WORLD, ierr)				
    	end if
    	
    	! print the number of function evaluations by each thread:
   	print '("fevals by Process ",i2,": ",i13)',  proc_num, fevals_proc

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all process to print
    	
    
    call MPI_REDUCE(fevals_proc,fevals_total,1, &
    		MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
    		MPI_COMM_WORLD,ierr)

    if (proc_num==0) then
    	int_approx=sum(int_subs)
        print '("Trapezoid approximation with ",i8," total points: ",es22.14)',&
            ns*n, int_approx
        print '("Total number of fevals: ",i10)', fevals_total
        endif

    call MPI_FINALIZE(ierr)

end program test3
