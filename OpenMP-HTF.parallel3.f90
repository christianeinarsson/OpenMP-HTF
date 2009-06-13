program laplsolv
	!-----------------------------------------------------------------------
	! Serial program for solving the heat conduction problem
	! on a square using the Jacobi method.
	! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
	! Modified by Berkant Savas (besav@math.liu.se) April 2006
	!-----------------------------------------------------------------------
	use omp_lib
	integer, parameter                  :: nmax=1000, maxiter=1000
	integer										:: n, threads
	double precision,parameter          :: tol=1.0E-3
	!double precision          :: tol=1.0E-3
	double precision,dimension(0:nmax+1,0:nmax+1) :: T1,T2
	double precision                    :: error,x
	integer                             :: i,j,k, solution
	character(len=20)                   :: str

	!Timing variables
	real(8)										:: start, end1, end2, end3

	!CLI argument buffer
	character*100								::	buffer
	character(len=30)							:: filename

	call getarg(1, buffer)
	read(buffer,*) threads
	call getarg(2, buffer)
	read(buffer,*) filename

	n = nmax
	if(threads > 8) then
		threads = 8
	end if
	call omp_set_num_threads(threads)

	! Set boundary conditions and initial values for the unknowns
	start = OMP_get_wtime()
!$omp parallel sections default(shared)
!$omp section
	T1(0:n , 1:n)		= 0.0D0
!$omp section
	T1(0:n , 0)			= 1.0D0
!$omp section
	T1(0:n , n+1)		= 1.0D0
!$omp section
	T1(n+1 , 0:n+1)	= 2.0D0
!$omp section
	T2(0:n , 1:n)		= 0.0D0
!$omp section
	T2(0:n , 0)     	= 1.0D0
!$omp section
	T2(0:n , n+1)   	= 1.0D0
!$omp section
	T2(n+1 , 0:n+1) 	= 2.0D0
!$omp end  parallel sections
	solution = 0

	end1 = OMP_get_wtime()

	! Solve the linear system of equations using the Jacobi method
!	!$omp parallel default(shared)
	do k=1,maxiter
		if (mod(k,2) == 0) then
			call jacobi(T1,T2,n,solution,tol)
			if( solution == 1) then
				solution = 2
				exit
			end if
		else
			call jacobi(T2,T1,n,solution,tol)
			if( solution == 1) then
				solution = 1
				exit
			end if
		end if
	end do
!	!$omp end parallel
	end2 = OMP_get_wtime()

	write(*,*) 'Writing solution...'
	open(unit=7,action='write',file=filename,status='unknown')
	write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
	write(unit=7, fmt='(A,I2,A,I6,A,I6,A,F7.3)') 'threads=', threads, ' n =', n, ' iterations =', k, ' time-to-solve =', end1-start
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

	end3 = OMP_get_wtime()

	write(unit=*,fmt='(a,f7.4,a)') 'Init time: ',end1-start, ' s'
	write(unit=*,fmt='(a,f7.4,a)') 'Solve time:',end2-end1 , ' s'
	write(unit=*,fmt='(a,f7.4a)') 'Write time:', end3-end2 , ' s'
	write(unit=*,fmt='(a,f7.4,a)') 'Total time:',end3-start, ' s'

	write(unit=*,fmt='(a,i5.3)') 'Number of Iterations:',k
	if (solution == 1) then
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T1(1,1)
	else
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T2(1,1)
	end if

end program laplsolv

subroutine jacobi (Told, Tnew, n, solution, tol)
	integer												:: n, solution
	double precision,dimension(0:n+1,0:n+1) 	:: Tnew,Told
	double precision									:: tol, error

!!$omp parallel default(shared)
!!$omp master
	error = 0.0D0
!!$omp end master

!!$omp   do  schedule(guided, 1) reduction(MAX:error)
!$omp  parallel do default(shared) schedule(guided, 1) reduction(MAX:error)
	do j=1,n
		Tnew(1:n,j)= ( Told(0:n-1,j) + Told(2:n+1,j) + Told(1:n,j+1) + Told(1:n,j-1) ) / 4.0D0
		error=max(error,maxval(abs(Tnew(1:n,j)-Told(1:n,j))))
	end do
!$omp end parallel do
!!$omp end do

!!$omp master
	!write(*,*) ' error: ', error
	if (error<tol) then
		solution = 1
	else
		solution = 0
	end if
!!$omp end master

!!$omp end parallel
end subroutine jacobi

