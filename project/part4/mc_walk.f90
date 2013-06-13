! module mc_walk  written by Daniel McInally for AMATH 583, Fianl Project, Part 3


module mc_walk

	implicit none
	integer :: nwalks


contains

subroutine random_walk(i0, j0, max_steps, ub, iabort)
	use problem_description
	
	
	implicit none
	
	real(kind=8), intent(in) :: i0,j0
	integer, intent(in) :: max_steps
	real(kind=8), intent(out) :: ub
	integer, intent(out) :: iabort
	
	!local variables
	real(kind=8) :: i,j,xb,yb
	real(kind=8) :: r(max_steps)
	integer :: istep
	
	
	! Starting point
	i = i0
	j = j0
	
	
	!generate random numbers
	call random_number(r)
	iabort=0
	
	
	
	! Taking the random steps
	do istep = 1,max_steps
		if (r(istep) < 0.25) then
			i = i-1
		elseif (r(istep) < 0.5) then
			i = i+1
		elseif (r(istep) < 0.75) then
			j = j-1
		else
			j = j+1
			endif
			
		
		! boundry check
		if (i*j*(nx+1-i)*(ny+1-j)==0) then
			xb = ax + i*dx
			yb = ay + j*dy
			
			ub = uboundary(xb, yb)
			
			exit
			endif
			
		if (istep==(max_steps-1)) then
			iabort=1
			exit
			endif
		
		enddo
	nwalks=nwalks+1
	
end subroutine random_walk

subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

	use mpi
	use problem_description
	
	implicit none
	real(kind=8), intent(in) :: i0,j0
	integer, intent(in) :: max_steps, n_mc
	real(kind=8), intent(out) :: u_mc
	integer, intent(out) :: n_success
	
	!local Variables
	real(kind=8) :: ub_sum,ub,i,j,a
	integer :: k,k2,iabort,ierr,proc_num,lastwalk,sender,numprocs
	integer, dimension(MPI_STATUS_SIZE) :: status
	logical ::done
	
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)
	
	ub_sum = 0.0
	n_success = 0
	
	if (proc_num==0) then
		! Single Processor
		if (numprocs==1) then
			do k = 1,n_mc
				call random_walk(i0, j0, max_steps, ub, iabort)
				if (iabort == 0) then
					ub_sum = ub_sum + ub
					n_success=n_success+1
					endif
				enddo
			u_mc = ub_sum / n_success
	  ! Multi Processor
	  	  else
			!send first round of processes
			lastwalk=0
			do k2 = 1,min(numprocs-1,n_mc)
				call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, k2, 1, & 
							MPI_COMM_WORLD, ierr)
				lastwalk=lastwalk+1
				enddo
				
			!receive results and send more work
			do k = 1,n_mc
				call MPI_RECV(ub,1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
						MPI_COMM_WORLD, status, ierr)
				iabort=status(MPI_TAG)
				sender=status(MPI_SOURCE)
				if (iabort == 0) then
					ub_sum = ub_sum + ub
					n_success=n_success+1
					endif
				!request more work
				if (lastwalk<n_mc) then
					call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, 1, &
									MPI_COMM_WORLD, ierr)
					lastwalk=lastwalk+1
				 else
				 	call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, 0, &
									MPI_COMM_WORLD, ierr)	
					endif
				enddo
			endif
				
				!tell extra procs work is not needed
			if (n_mc<(numprocs-1)) then
				do k =n_mc+1,numprocs-1
					call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, k, 0, &
									MPI_COMM_WORLD, ierr)	
					enddo
				endif	
			u_mc = ub_sum / n_success	
			endif
		
		! WORKER PROCESSES
		if (proc_num /= 0) then
			done=.false.
			do while (.not.(done))
				call MPI_RECV(a,0 , MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, &
						MPI_COMM_WORLD, status, ierr)
				k2=status(MPI_TAG)
				
				if (k2==0) then
					done=.true.
					exit
					endif
				call random_walk(i0, j0, max_steps, ub, iabort)
				call MPI_SEND(ub, 1, MPI_DOUBLE_PRECISION, 0, iabort, &
									MPI_COMM_WORLD, ierr)
				enddo
			endif			

end subroutine many_walks

end module mc_walk
