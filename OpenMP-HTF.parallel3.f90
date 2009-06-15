program laplsolv
	!-----------------------------------------------------------------------
	! Parallel program for solving the heat conduction problem
	! on a square using the Jacobi method.
	! Serial code written by Fredrik Berntsson (frber@math.liu.se) March 2003
	! Serial code modified by Berkant Savas (besav@math.liu.se) April 2006
	! Serial code parallelized by Christian Einarsson and Joel Purra June 2009
	!-----------------------------------------------------------------------
	use omp_lib
	integer, parameter								:: n=1000, maxiter=1000			!Problemsize and an restriction on itterations
	double precision,parameter						:: tol=1.0E-3						!Computationalerror threshold
	double precision,dimension(0:n+1,0:n+1)	:: T1,T2								!Matrices containing problem/solution
	integer												:: i,k, solution, threads		!Iteration counters, etc.
	character(len=30)									:: str, filename					!Strings storing outputformat and output filename
	real(8)												:: start, end1, end2, end3		!Timing variables
	character*100										::	clibuffer						!buffer for rading CLI-arguments

	!Getting CLI-parameters
	call getarg(1, clibuffer)
	read(clibuffer,*) threads
	call getarg(2, clibuffer)
	read(clibuffer,*) filename

	!Setting number of threads.
	if(threads > 8) then
		threads = 8
	end if
	call omp_set_num_threads(threads)

	!Output problem info
	write(*,*) '* * *  Solving for...  * * * '
	write(*,*) '  threads: ', threads
	write(*,*) '  n:       ', n

	!Record start time
	start = OMP_get_wtime()

	! Set boundary conditions and initial values for the unknowns.
	!  This part is divided in independent prallel sections, mostly becuse it's
	!  posible, only a smal speedup is gained. The probable reason is that the
	!  sections are to small.
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

	!Record then end time for initialization and start time for the solving.
	end1 = OMP_get_wtime()

	! Solve the linear system of equations using the Jacobi method.
	!	The loop repetedly calls the subroutine jacobi while alternating
	!	between the matrices. First jacobi reads from T2 and writes to T1
	!  in the next itteration its the otherway around.
	!  The loop is exited when one solution meets the given tolerance
	!	requirement. The number of the matrix containing the solution are
	!  stored in the solution variable.
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

	!Record end time fore equation solving and start time for writing
	end2 = OMP_get_wtime()

	!Write the solved problem to a file.
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

	!Record time for global end.
	end3 = OMP_get_wtime()

	!Write information to console.
	write(unit=*,fmt='(a,f12.4,a)') 'Init time: ',end1-start, ' s'
	write(unit=*,fmt='(a,f12.4,a)') 'Solve time:',end2-end1 , ' s'
	write(unit=*,fmt='(a,f12.4,a)') 'Write time:', end3-end2 , ' s'
	write(unit=*,fmt='(a,f12.4,a)') 'Total time:',end3-start, ' s'
	write(unit=*,fmt='(a,i5.3)') 'Number of Iterations:',k
	if (solution == 1) then
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T1(1,1)
	else
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T2(1,1)
	end if

end program laplsolv


! Subroutine for doing one itteration of the algorithm in parallel.

subroutine jacobi (Told, Tnew, n, solution, tol)
	integer												:: n, solution
	double precision,dimension(0:n+1,0:n+1) 	:: Tnew,Told
	double precision									:: tol, error, j

	error = 0.0D0
	! Share the work of the do loop over multiple threads using guded scheduling.
	! By using two matrices one fore reading and one for writing all the
	! individual calculations becomes independant, thus posible to do in parallel.
!$omp  parallel do default(shared) schedule(guided, 1) reduction(MAX:error)
	do j=1,n

		Tnew(1:n,j)= ( Told(0:n-1,j) + Told(2:n+1,j) + Told(1:n,j+1) + Told(1:n,j-1) ) / 4.0D0
		error=max(error,maxval(abs(Tnew(1:n,j)-Told(1:n,j))))
	end do
!$omp end parallel do

	! Check if the error is below the given tolerance level.
	if (error<tol) then
		solution = 1
	else
		solution = 0
	end if

end subroutine jacobi
