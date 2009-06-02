program laplsolv
	!-----------------------------------------------------------------------
	! Serial program for solving the heat conduction problem
	! on a square using the Jacobi method.
	! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
	! Modified by Berkant Savas (besav@math.liu.se) April 2006
	!-----------------------------------------------------------------------
	integer, parameter                  :: n=1000, maxiter=1000
	double precision,parameter          :: tol=1.0E-3
	double precision,dimension(0:n+1,0:n+1) :: T1,T2
	double precision,dimension(n)       :: tmp1,tmp2
	double precision                    :: error,x
	real                                :: time1,time0
	integer                             :: i,j,k
	character(len=20)                   :: str

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
	call cpu_time(time0)
	do k=1,maxiter

		!tmp1=T(1:n,0)
		error=0.0D0
		!
		!Två lopar första läser T1 och skriver T2 den andra tvärt om
		!Kolla feltollerans mellan looparan
		!Parallelisera loparna se till att det görs reduction max på error
		!T1,T2 shared
		!Barrier mellan looparna
		!inga tempvariabler behövs.
		!Räkna upp k emellan
		! hur styra output för att slippa koppiera resultatet? k modulo 2 vid läsning?

!$omp parallel do shared(T1,T2) reduction(MAX:error)
		do j=1,n
			!tmp2=T(1:n,j)
			T2(1:n,j)= ( T1(0:n-1,j) + T1(2:n+1,j) + T1(1:n,j+1) + T1(1:n,j-1) ) / 4.0D0 !tmp1 ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
			!tmp1=tmp2
		end do
!$omp end parallel  do
		if (error<tol) then
			solution = 2
			exit
		end if

!$omp barrier
		!Barrier bör nog användas här...

		!får inte inkrementera k, hur välja mellan T1 och T2?
		!k = k + 1

!$omp parallel do shared(T1,T2) reduction(MAX:error)
		do j=1,n
			!tmp2=T(1:n,j)
			T1(1:n,j)= ( T2(0:n-1,j) + T2(2:n+1,j) + T2(1:n,j+1) + T2(1:n,j-1) ) / 4.0D0 !tmp1 ) / 4.0D0
			error=max(error,maxval(abs(T2(1:n,j)-T1(1:n,j))))
			!tmp1=tmp2
		end do
!$omp end parallel  do

		if (error<tol) then
			solution = 1
			exit
		end if

	end do

	call cpu_time(time1)

	!solution = modulo(k,2)

	write(unit=*,fmt=*) 'Time:',time1-time0,'Number of Iterations:',k*2 + solution - 1
	if (solution < 2) then
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T1(1,1)
	else
		write(unit=*,fmt=*) 'Temperature of element Tx(1,1)  =',T2(1,1)
	end if
	! Uncomment the next part if you want to write the whole solution
	! to a file. Useful for plotting.

	open(unit=7,action='write',file='result.dat',status='unknown')
	write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
	if(solution < 2) then
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
