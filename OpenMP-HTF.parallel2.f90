program openmphtf
	use omp_lib

	integer, parameter  								:: n = 1000
	double precision,parameter			         :: tol=1.0E-3

	integer												:: j,i
	!real 												:: T(0,n+1),Tnew(0,n+1), error
	double precision, dimension(0:n+1, 0:n+1) :: Tnew, T

	real(8)													:: start1,start2,start3,end1,end2,end3, error


	call omp_set_num_threads(1)
	!call omp_set_dynamic(1)
	start1 = OMP_get_wtime()
		T=0.0D0
		T(0:n+1 , 0)     = 1.0D0
		T(0:n+1 , n+1)   = 1.0D0
		T(n+1   , 0:n+1) = 2.0D0
		Tnew=0.0D0
		Tnew(0:n+1 , 0)     = 1.0D0
		Tnew(0:n+1 , n+1)   = 1.0D0
		Tnew(n+1   , 0:n+1) = 2.0D0
		error = 0
	start2 = OMP_get_wtime()

!		!$omp do
!			do j=0,n+1
!				do i=0,n+1
!					if( i == n+1) then
!						Tnew(i,j) = 2.0D0
!						T(i,j) = 2.0D0
!					else if(j==0 .OR. j == n + 1) then
!						Tnew = 1.0D0
!					else
!						Tnew(i,j) = 0.0D0
!						T(i,j) = 0.0D0
!					end if
!				end do
!			end do
!		!$omp end do
	!$omp parallel default	(shared)
		!$omp do private(j) reduction(MAX: error)
				do j=1,n
					do i=1,n
						Tnew(i,j) = ( T(i-1,j) + T( i,j+1) + T( i , j-1 ) + T(i+2 , j) ) / 4
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
	!$omp end parallel
write(*,*) "KUK"
		do while (error > tol)
		error = 0.0D0
	!$omp parallel default	(shared)
		!$omp do private(j) reduction(MAX: error)
				do j=1,n
					do i=1,n
						Tnew(i,j) = ( T(i-1,j) + T( i,j+1) + T( i , j-1 ) + T(i+2 , j) ) / 4
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
			!write(*,*) "error: ", error
	!$omp end parallel
		end do

	end1 = OMP_get_wtime()

	start1 = start2 - start1
	end1 = end1 - start2

	write(*,*) "init time: ", start1
	write(*,*) "solve time: ", end1

	write(*,*) "T(1,1): ", T(1,1)

end program openmphtf

subroutine jacobi(T, Tnew, n, error)
	real :: error
	integer :: n
	double precision, dimension (0:n+1,0:n+1) :: T, Tnew
	!$omp parallel default	(shared)
		!$omp do private(j) reduction(MAX: error)
				do j=1,n
					do i=1,n
						Tnew(i,j) = ( T(i-1,j) + T( i,j+1) + T( i , j-1 ) + T(i+2 , j) ) / 4
						error = max(error, abs(Tnew(i,j) - T(i,j) ))
					end do
				end do
			!$omp end do
	!$omp end parallel
end subroutine jacobi




