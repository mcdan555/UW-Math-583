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
	
end subroutine random_walk

subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

	use problem_description
	
	implicit none
	real(kind=8), intent(in) :: i0,j0
	integer, intent(in) :: max_steps, n_mc
	real(kind=8), intent(out) :: u_mc
	integer, intent(out) :: n_success
	
	!local Variables
	real(kind=8) :: ub_sum,ub,i,j
	integer :: k,iabort
	
	ub_sum = 0.0
	n_success = 0
	
	do k = 1,n_mc
		call random_walk(i0, j0, max_steps, ub, iabort)
		nwalks=nwalks+1
		if (iabort == 0) then
			ub_sum = ub_sum + ub
			n_success=n_success+1
			endif
		enddo
		
	u_mc = ub_sum / n_success

end subroutine many_walks

end module mc_walk
