program laplsolv
	!-----------------------------------------------------------------------
	! Serial program for solving the heat conduction problem
	! on a square using the Jacobi method.
	! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
	! Modified by Berkant Savas (besav@math.liu.se) April 2006
	!-----------------------------------------------------------------------
	use omp_lib
	integer, parameter                  :: nmax=1000, maxiter=228
	integer										:: n, threads
	!double precision,parameter          :: tol=1.0E-3
	double precision          :: tol=1.0E-3
	double precision,dimension(0:nmax+1,0:nmax+1) :: T1,T2
	!double precision,allocatable 			:: T1(:,:),T2(:,:)
	double precision                    :: error,x
	integer                             :: i,j,k, solution
	character(len=20)                   :: str

	!Timing variables
	real(8)										:: start_t, end_t

	!CLI argument buffer
	character*100								::	buffer
	character(len=30)							:: filename

	call getarg(1, buffer)
	read(buffer,*) threads
	call getarg(2, buffer)
	read(buffer,*) n
	call getarg(3, buffer)
	read(buffer,*) filename

	if(n > nmax) then
		n = nmax
	end if
	if(threads > 8) then
		threads = 8
	end if
	call omp_set_num_threads(threads)
	!allocate( T1(0:n+1,0:n+1) )
	!allocate( T2(0:n+1,0:n+1) )

	! Set boundary conditions and initial values for the unknowns
	T1=0.0D0
	T1(0:n+1 , 0)     = 1.0D0
	T1(0:n+1 , n+1)   = 1.0D0
	T1(n+1   , 0:n+1) = 2.0D0
	T2=0.0D0
	T2(0:n+1 , 0)     = 1.0D0
	T2(0:n+1 , n+1)   = 1.0D0
	T2(n+1   , 0:n+1) = 2.0D0
	solution = 0

	start_t = OMP_get_wtime()

	! Solve the linear system of equations using the Jacobi method
	do k=1,maxiter
		solution = 0;
		error=0.0D0
!$omp parallel shared (T1,T2,error, solution,tol)
!!$omp parallel do shared(T1,T2) reduction(MAX:error)
!$omp  do schedule(guided, 1) reduction(MAX:error)
		do j=1,n
			T2(1:n,j)= ( T1(0:n-1,j) + T1(2:n+1,j) + T1(1:n,j+1) + T1(1:n,j-1) ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
		end do
!!$omp end parallel  do
!$omp end do
!$omp master
		if (error<tol) then
			!solution = 2
			solution = solution + 2
			!exit
		end if
		error=0.0D0
!$omp flush (error)
!$omp end master

!!$omp parallel do shared(T1,T2) reduction(MAX:error)
!!$omp do shared(T1,T2) reduction(MAX:error)
!$omp do schedule(guided, 1) reduction(MAX:error)
		do j=1,n
			T1(1:n,j)= ( T2(0:n-1,j) + T2(2:n+1,j) + T2(1:n,j+1) + T2(1:n,j-1) ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
		end do
!!$omp end parallel  do
!$omp end  do
!$omp master
		if (error<tol) then
			solution = solution + 1
			!solution = 1
			!exit
		end if
!$omp end master
!$omp end parallel
		if(solution > 0) then
			exit
		end if
	end do

	end_t = OMP_get_wtime()
	start_t = end_t - start_t
	k = k*2 - (min(2,solution) - 1)

	write(unit=*,fmt='(a,f7.3,a,i5.3)') 'Time:',start_t,'Number of Iterations:',k
	if (solution < 2) then
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T1(1,1)
	else
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T2(1,1)
	end if

	open(unit=7,action='write',file=filename,status='unknown')
	write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
	write(unit=7, fmt='(A,I2,A,I6,A,I6,A,F7.3)') 'threads=', threads, ' n =', n, ' iterations =', k, ' time-to-solve =', start_t
	if(solution == 1) then
		do i=0,n+1
			write (unit=7,fmt=str) T1(i,0:n+1)
		end do
	else
		do i=0,n+1
			write (unit=7,fmt=str) T2(i,0:n+1)
		end do
	end if
	close(unit=7)

end program laplsolv
