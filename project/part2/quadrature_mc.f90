! $MYHPSC/project/part2/quadrature_mc.f90


module quadrature_mc


contains

	function quad_mc(g, a, b, ndim, npoints)
    	use random_util, only: init_random_seed
		
		implicit none
		! Parameters
		integer, intent(in) :: ndim, npoints
		real(kind=8), intent(in) :: a(ndim),b(ndim)
		real(kind=8), external :: g
		! Internal variables
		real(kind=8) :: quad_mc,gsum,x(ndim),k(ndim,npoints),V
		real(kind=8), allocatable :: r(:)
		integer :: n,seed,i,j
		
		
		!Seed the random number and load the random number array (number 0 to 1)
		n=ndim*npoints
		!seed=7
		!call init_random_seed(seed)
		allocate(r(n))
		call random_number(r)
		do i=1,npoints
			do j=1,ndim
				k(j,i)=a(j)*(1-r(j+(i-1)*ndim))+b(j)*(r(j+(i-1)*ndim))
				enddo
			enddo
		
		
		!Integrate
		gsum=0
		
		do i=1,npoints
			gsum = gsum + g(k(:,i),ndim)
			enddo
		
		! Determine Volume and compute final mc approximation	
		V=product(b-a)
		
		quad_mc=V/npoints*gsum
		
	end function quad_mc

end module quadrature_mc
