program laplsolv
	!-----------------------------------------------------------------------
	! Serial program for solving the heat conduction problem
	! on a square using the Jacobi method.
	! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
	! Modified by Berkant Savas (besav@math.liu.se) April 2006
	!-----------------------------------------------------------------------

	integer, parameter                  :: nmax=1000, maxiter=1000
	integer										:: n, threads
	double precision,parameter          :: tol=1.0E-3
	double precision,dimension(0:nmax+1,0:nmax+1) :: T1,T2
	!double precision,dimension(n)       :: tmp1,tmp2
	double precision                    :: error,x
	!real                                :: time1,time0
	!integer                             :: time1,time0
	integer                             :: i,j,k
	character(len=20)                   :: str

	!Timing variables
	!real(kind(0.0D0)) 						:: r_time !start_d, end_d
	double precision	 						:: r_time !start_d, end_d
	integer(4)									:: start_t, end_t, count_rate, count_max
	!CLI argument buffer
	character*100								::	buffer
	character(len=30)							:: filename

	call getarg(1, buffer)
	read(buffer,*) threads
	call getarg(2, buffer)
	read(buffer,*) n
	call getarg(3, buffer)
	read(buffer,*) filename

	if(n > 10000) then
		n = 10000
	end if
	if(threads > 8) then
		threads = 8
	end if


	call omp_set_num_threads(threads)

	! Set boundary conditions and initial values for the unknowns
	T1=0.0D0
	T1(0:n+1 , 0)     = 1.0D0
	T1(0:n+1 , n+1)   = 1.0D0
	T1(n+1   , 0:n+1) = 2.0D0
	T2=0.0D0
	T2(0:n+1 , 0)     = 1.0D0
	T2(0:n+1 , n+1)   = 1.0D0
	T2(n+1   , 0:n+1) = 2.0D0


	! Solve the linear system of equations using the Jacobi method
	!call cpu_time(time0)
	call system_clock(start_t)
!$omp parallel do shared(T1,T2,start_t)
	do k=1,maxiter

		error=0.0D0

!!$omp parallel do shared(T1,T2,start_t) reduction(MAX:error)
!$omp do  reduction(MAX:error)
		do j=1,n
			!tmp2=T(1:n,j)
			T2(1:n,j)= ( T1(0:n-1,j) + T1(2:n+1,j) + T1(1:n,j+1) + T1(1:n,j-1) ) / 4.0D0 !tmp1 ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
			!tmp1=tmp2
		end do
!!$omp end parallel  do
!$omp end do
		if (error<tol) then
			solution = 2
			exit
		end if

!$omp barrier
		error=0.0D0


!!$omp parallel do shared(T1,T2, start_t) reduction(MAX:error)
!$omp do reduction(MAX:error)
		do j=1,n
			T1(1:n,j)= ( T2(0:n-1,j) + T2(2:n+1,j) + T2(1:n,j+1) + T2(1:n,j-1) ) / 4.0D0 !tmp1 ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
		end do
!!$omp end parallel  do
!$omp end  do

		if (error<tol) then
			solution = 1
			exit
		end if

	end do


	!call cpu_time(time1)
	call system_clock(end_t, count_rate, count_max)
	if (end_t > start_t)  then
		count_max = 0
	end if
	r_time = ( dble(end_t) + dble(count_max) - dble(start_d)) / dble(count_rate)

	write(unit=*,fmt=*) 'Time:',r_time,'Number of Iterations:',k*2 + (2 - solution)
	if (solution < 2) then
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T1(1,1)
	else
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T2(1,1)
	end if

	open(unit=7,action='write',file=filename,status='unknown')
	write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
	write(unit=7, fmt=*) 'threads=', threads, ' n=', n, ' iterations=', k, ' time-to-solve=', r_time
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
