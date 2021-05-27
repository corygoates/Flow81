module mesh_m

    implicit none

    type boundary_condition_t

        character :: type ! N for Neumann, D for Dirichlet
        real(kind=8) :: a, b, c ! Boundary condition formulated as a+b*x+c*x^2

    end type boundary_condition_t
    
    type mesh_t

        character(100) :: input_file
        integer :: n_rows, n_cols, N_cp, max_iter
        integer,allocatable,dimension(:) :: W_ind, E_ind, S_ind, N_ind
        real(kind=8) :: x_min, x_max, y_min, y_max
        real(kind=8) :: dx, dy, gamma, omega, conv
        real(kind=8),allocatable,dimension(:) :: x_cp, y_cp, AP, AW, AS, AE, AN, Sp, Su, phi
        type(boundary_condition_t) :: W_bc, S_bc, E_bc, N_bc

    end type mesh_t

    integer allocated

contains


subroutine mesh_allocate(t)
    type(mesh_t) :: t
    if(allocated.eq.1) call mesh_deallocate(t)

    allocate(t%W_ind(t%N_cp))
    allocate(t%E_ind(t%N_cp))
    allocate(t%S_ind(t%N_cp))
    allocate(t%N_ind(t%N_cp))
    allocate(t%x_cp(t%n_cols))
    allocate(t%y_cp(t%n_rows))
    allocate(t%AP(t%N_cp))
    allocate(t%AW(t%N_cp))
    allocate(t%AS(t%N_cp))
    allocate(t%AE(t%N_cp))
    allocate(t%AN(t%N_cp))
    allocate(t%Su(t%N_cp))
    allocate(t%Sp(t%N_cp))
    allocate(t%phi(t%N_cp))

    allocated = 1

end subroutine mesh_allocate


subroutine mesh_deallocate(t)
    type(mesh_t) :: t

    deallocate(t%W_ind)
    deallocate(t%E_ind)
    deallocate(t%S_ind)
    deallocate(t%N_ind)
    deallocate(t%x_cp)
    deallocate(t%y_cp)
    deallocate(t%AP)
    deallocate(t%AW)
    deallocate(t%AS)
    deallocate(t%AE)
    deallocate(t%AN)
    deallocate(t%Su)
    deallocate(t%Sp)
    deallocate(t%phi)

    allocated = 0
end subroutine mesh_deallocate


subroutine mesh_load_input(t, inp)
    type(mesh_t) :: t
    character(100) :: inp
    integer :: i,j

    t%input_file = inp

    ! Read input file
    open(1,file=t%input_file)
        read(1,*)
        read(1,*)
        read(1,*) t%x_min, t%x_max, t%y_min, t%y_max
        read(1,*)
        read(1,*) t%n_rows, t%n_cols
        read(1,*)
        read(1,*) t%gamma
        read(1,*)
        read(1,*)
        read(1,*) t%W_bc%type, t%W_bc%a, t%W_bc%b, t%W_bc%c
        read(1,*)
        read(1,*) t%S_bc%type, t%S_bc%a, t%S_bc%b, t%S_bc%c
        read(1,*)
        read(1,*) t%E_bc%type, t%E_bc%a, t%E_bc%b, t%E_bc%c
        read(1,*)
        read(1,*) t%N_bc%type, t%N_bc%a, t%N_bc%b, t%N_bc%c
        read(1,*)
        read(1,*)
        read(1,*) t%omega, t%conv, t%max_iter
    close(1)

    ! Parse some derived quantities
    t%dx = (t%x_max-t%x_min)/t%n_cols
    t%dy = (t%y_max-t%y_min)/t%n_rows
    t%N_cp = t%n_rows*t%n_cols

    ! Allocate
    call mesh_allocate(t)

    ! Initialize control point locations
    do i = 1, t%n_rows
        t%y_cp(i) = t%y_min+(i-0.5)*t%dx
    end do
    do j = 1, t%n_cols
        t%x_cp(j) = t%x_min+(j-0.5)*t%dx
    end do

    write(*,*) "Successfully initialized mesh."

end subroutine mesh_load_input


subroutine mesh_set_up_matrices(t)
    ! Sets up the matrices for actually solving the problem

    implicit none
    type(mesh_t) :: t
    integer :: i, j, ind
    real(kind=8) :: phi_bc, dphi_dn

    ! Loop through control points
    ind = 0
    do i = 1, t%n_rows
        do j = 1, t%n_cols

            ! Calculate indices for this point
            ind = ind + 1
            t%W_ind(ind) = ind-1
            t%E_ind(ind) = ind+1
            t%S_ind(ind) = ind-t%n_cols
            t%N_ind(ind) = ind+t%n_cols

            ! Initialize solution
            t%phi(ind) = 0.0

            ! Calculate central coefficients
            t%AW(ind) = t%gamma*t%dy/t%dx
            t%AE(ind) = t%gamma*t%dy/t%dx
            t%AS(ind) = t%gamma*t%dx/t%dy
            t%AN(ind) = t%gamma*t%dx/t%dy
            t%Su(ind) = 0
            t%Sp(ind) = 0

            ! Handle west boundary
            if (j .eq. 1) then
                if (t%W_bc%type .eq. 'D') then
                    phi_bc = t%W_bc%a+t%W_bc%b*t%y_cp(i)+t%W_bc%c*t%y_cp(i)**2.0
                    t%Su(ind) = t%Su(ind) + 2.0*t%gamma*t%dy/t%dx
                    t%Sp(ind) = t%Sp(ind) - 2.0*t%gamma*t%dy*phi_bc/t%dx
                else
                    dphi_dn = t%W_bc%a+t%W_bc%b*t%y_cp(i)+t%W_bc%c*t%y_cp(i)**2.0
                    t%Sp(ind) = t%Sp(ind) - t%dy*t%gamma*dphi_dn
                end if
                t%AW(ind) = 0.0
                t%W_ind(ind) = 1 ! Not needed in this case
            end if

            ! Handle east boundary
            if (j .eq. t%n_cols) then
                if (t%E_bc%type .eq. 'D') then
                    phi_bc = t%E_bc%a+t%E_bc%b*t%y_cp(i)+t%E_bc%c*t%y_cp(i)**2.0
                    t%Su(ind) = t%Su(ind) + 2.0*t%gamma*t%dy/t%dx
                    t%Sp(ind) = t%Sp(ind) - 2.0*t%gamma*t%dy*phi_bc/t%dx
                else
                    dphi_dn = t%E_bc%a+t%E_bc%b*t%y_cp(i)+t%E_bc%c*t%y_cp(i)**2.0
                    t%Sp(ind) = t%Sp(ind) + t%dy*t%gamma*dphi_dn
                end if
                t%AE(ind) = 0.0
                t%E_ind(ind) = 1 ! Not needed in this case
            end if

            ! Handle south boundary
            if (i .eq. 1) then
                if (t%S_bc%type .eq. 'D') then
                    phi_bc = t%S_bc%a+t%S_bc%b*t%x_cp(j)+t%S_bc%c*t%x_cp(j)**2.0
                    t%Su(ind) = t%Su(ind) + 2.0*t%gamma*t%dx/t%dy
                    t%Sp(ind) = t%Sp(ind) - 2.0*t%gamma*t%dx*phi_bc/t%dy
                else
                    dphi_dn = t%S_bc%a+t%S_bc%b*t%x_cp(j)+t%S_bc%c*t%x_cp(j)**2.0
                    t%Sp(ind) = t%Sp(ind) - t%dx*t%gamma*dphi_dn
                end if
                t%AS(ind) = 0.0
                t%S_ind(ind) = 1 ! Not needed in this case
            end if

            ! Handle north boundary
            if (i .eq. t%n_rows) then
                if (t%N_bc%type .eq. 'D') then
                    phi_bc = t%N_bc%a+t%N_bc%b*t%x_cp(j)+t%N_bc%c*t%x_cp(j)**2.0
                    t%Su(ind) = t%Su(ind) + 2.0*t%gamma*t%dx/t%dy
                    t%Sp(ind) = t%Sp(ind) - 2.0*t%gamma*t%dx*phi_bc/t%dy
                else
                    dphi_dn = t%N_bc%a+t%N_bc%b*t%x_cp(j)+t%N_bc%c*t%x_cp(j)**2.0
                    t%Sp(ind) = t%Sp(ind) + t%dx*t%gamma*dphi_dn
                end if
                t%AN(ind) = 0.0
                t%N_ind(ind) = 1 ! Not needed in this case
            end if

            ! Calculate P coefficient
            t%AP(ind) = t%AW(ind) + t%AE(ind) + t%AN(ind) + t%AS(ind) - t%Sp(ind)

        end do
    end do

end subroutine mesh_set_up_matrices

subroutine mesh_sor(t)

    type(mesh_t) :: t
    real(kind=8) :: max_corr, corr
    integer :: i, iteration

    write(*,*)
    write(*,*) "Running SOR..."

    ! Iterative SOR
    iteration = 0
    max_corr = 1.0+t%conv
    do while (max_corr>t%conv .and. iteration<t%max_iter)

        max_corr = 0.0
        
        ! Loop through control points
        iteration = iteration + 1
        do i = 1, t%N_cp

            ! Calculate correction
            corr = t%AW(i)*t%phi(t%W_ind(i)) + t%AW(i)*t%phi(t%E_ind(i))
            corr = corr + t%AS(i)*t%phi(t%S_ind(i)) + t%AN(i)*t%phi(t%N_ind(i))
            corr = corr + t%Su(i) - t%AP(i)*t%phi(i)
            corr = corr*t%omega/t%AP(i)

            ! Check magnitude
            if (abs(corr)>max_corr) then
                max_corr  = abs(corr)
            end if

            ! Apply correction
            t%phi(i) = t%phi(i) + corr

        end do

        ! Output progress
        write(*,*) iteration, max_corr

    end do
    write(*,*) "Complete! Final maximum correction: ", max_corr

end subroutine mesh_sor
    
end module mesh_m