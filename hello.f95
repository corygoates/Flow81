program helloWorld
    implicit none
    integer :: N
    real, dimension (:,:), allocatable :: matrix
    integer i,j

    N = 2000
    allocate(matrix(N,N))

    do i=1,N
        do j=1,N
            matrix(i,j) = i*j
            matrix(i,j) = matrix(i,j)-i+j
        end do
    end do

    do i=1,N
        do j=1,N
            print *,matrix(i,j)
        end do
    end do
    
end program helloWorld