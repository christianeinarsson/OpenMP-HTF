
integer, parameter  		:: n = 1000
integer						:: j,i
real 							:: T(0,n+1),Tnew(0,n+1), error


!$omp parallel defaults(shared)
!$omp do
do j=1,n
	do i=1,n
		Tnew(i,j) = 0.0
		T(i,j) = 0.0
	end do
end do
!$omp end do
do while (error > tol)
!$omp do private(j) reduction(MAX: error)
	do j=1,n
		do i=1,n
			Tnew(i,j) = ( T(i-1,j) + T( , ) + T( , ) + T( , ) ) / 4
			error = max(error, abs(Tnew(i,j) - T(i,j) ))
		end do
	end do
!$omp end do
!$omp do
	do j=1,n
		do i=1,n
			T(i,j)=Tnew(i,j)
		end do
	end do
!$omp end do
end do
!$omp end parallel
