module ns_case_m

    use boundary_condition_m

    implicit none

    type ns_case_t

        character(100) :: input_file, results_file
        real(kind=8) :: x_min, x_max, y_min, y_max
        real(kind=8) :: omegaU, omegaV, omegaP, mu, rho
        real(kind=8) :: dx, dy
        integer :: nx, ny, outer_iterations
        type(bc_t) :: Uw_bc, Vw_bc, Us_bc, Vs_bc, Ue_bc, Ve_bc, Un_bc, Vn_bc
        real(kind=8),allocatable,dimension(:,:) :: U, U_old, V, V_old, P ! Values
        real(kind=8),allocatable,dimension(:,:) :: AuP, AuW, AuS, AuE, AuN ! x-mom coefficients
        real(kind=8),allocatable,dimension(:,:) :: AvP, AvW, AvS, AvE, AvN ! y-mom coefficients
        real(kind=8),allocatable,dimension(:) :: x_cp, y_cp

    end type ns_case_t

    integer allocated

contains

subroutine ns_case_allocate(t)
    type(ns_case_t) :: t

    ! Check allocation
    if(allocated.eq.1) call ns_case_deallocate(t)

    ! Allocate arrays
    allocate(t%P(1:t%nx,1:t%ny))
    allocate(t%U(1:t%nx+1,0:t%ny+1))
    allocate(t%V(0:t%nx+1,1:t%ny+1))
    allocate(t%U_old(1:t%nx+1,0:t%ny+1))
    allocate(t%V_old(0:t%nx+1,1:t%ny+1))
    allocate(t%x_cp(t%nx))
    allocate(t%y_cp(t%ny))
    allocate(t%AuW(2:t%nx,1:t%ny))
    allocate(t%AuS(2:t%nx,1:t%ny))
    allocate(t%AuE(2:t%nx,1:t%ny))
    allocate(t%AuN(2:t%nx,1:t%ny))
    allocate(t%AuP(2:t%nx,1:t%ny))

    allocated = 1

end subroutine ns_case_allocate


subroutine ns_case_deallocate(t)
    type(ns_case_t) :: t

    ! Clean up memory
    deallocate(t%P)
    deallocate(t%U)
    deallocate(t%V)
    deallocate(t%U_old)
    deallocate(t%V_old)
    deallocate(t%x_cp)
    deallocate(t%y_cp)

    allocated = 0

end subroutine ns_case_deallocate


subroutine ns_case_load_input(t, input_file)
    type(ns_case_t) :: t
    character(100) :: input_file

    ! Read input file
    open(1,file=input_file)
        read(1,*) ! Header

        read(1,*) ! Geometry
        read(1,*)
        read(1,*) t%x_min, t%x_max, t%y_min, t%y_max
        read(1,*)
        read(1,*) t%nx, t%ny

        read(1,*) ! Flow properties
        read(1,*)
        read(1,*) t%mu, t%rho

        read(1,*) ! Boundary conditions
        read(1,*)
        read(1,*)
        read(1,*) t%Uw_bc%type, t%Uw_bc%a, t%Uw_bc%b, t%Uw_bc%c
        read(1,*)
        read(1,*) t%Vw_bc%type, t%Vw_bc%a, t%Vw_bc%b, t%Vw_bc%c
        read(1,*)
        read(1,*)
        read(1,*) t%Us_bc%type, t%Us_bc%a, t%Us_bc%b, t%Us_bc%c
        read(1,*)
        read(1,*) t%Vs_bc%type, t%Vs_bc%a, t%Vs_bc%b, t%Vs_bc%c
        read(1,*)
        read(1,*)
        read(1,*) t%Ue_bc%type, t%Ue_bc%a, t%Ue_bc%b, t%Ue_bc%c
        read(1,*)
        read(1,*) t%Ve_bc%type, t%Ve_bc%a, t%Ve_bc%b, t%Ve_bc%c
        read(1,*)
        read(1,*)
        read(1,*) t%Un_bc%type, t%Un_bc%a, t%Un_bc%b, t%Un_bc%c
        read(1,*)
        read(1,*) t%Vn_bc%type, t%Vn_bc%a, t%Vn_bc%b, t%Vn_bc%c

        read(1,*) ! Solver Options
        read(1,*)
        read(1,*) t%outer_iterations
        read(1,*)
        read(1,*) t%omegaU, t%omegaV, t%omegaP
        read(1,*)
        read(1,*) t%results_file

    close(1)

    call ns_case_allocate(t)
    call ns_case_initialize(t)

end subroutine ns_case_load_input


subroutine ns_case_initialize(t)
    type(ns_case_t) :: t

    ! Initialize grid
    t%dx = (t%x_max-t%x_min)/t%nx
    t%dy = (t%y_max-t%y_min)/t%ny

    ! Guess values for U, V, and P
    t%P = 0.0
    t%U = 0.0
    t%V = 0.0

    ! Initialize boundary conditions
    call ns_case_apply_bc(t)

    ! Initialize old values
    t%U_old = t%U
    t%V_old = t%V

end subroutine ns_case_initialize


subroutine ns_case_apply_bc(t)
    type(ns_case_t) :: t
    integer :: i

    ! Place scalar control points
    do i=1,t%nx
        t%x_cp(i) = (i+0.5)*t%dx
    end do

    do i=1,t%ny
        t%y_cp(i) = (i+0.5)*t%dy
    end do

    ! U boundary conditions
    do i=0,t%ny+1

        ! West
        if (t%Uw_bc%type .eq. 'D') then
            t%U(1,i) = bc_get_value(t%Uw_bc, t%y_cp(i))
        end if

        ! East
        if (t%Ue_bc%type .eq. 'D') then
            t%U(t%nx+1,i) = bc_get_value(t%Ue_bc, t%y_cp(i))
        end if

    end do

    do i=1,t%nx+1

        ! South
        if (t%Us_bc%type .eq. 'D') then
            t%U(i,0) = bc_get_value(t%Us_bc, t%x_cp(i))
        end if

        ! North
        if (t%Ue_bc%type .eq. 'D') then
            t%U(i,t%ny+1) = bc_get_value(t%Un_bc, t%x_cp(i))
        end if

    end do

    ! V
    t%V(:,0) = 0.0 ! South
    t%V(:,t%ny+1) = 0.0 ! North
    t%V(1,:) = 0.0 ! West
    t%V(t%nx+1,:) = 0.0 ! East

end subroutine ns_case_apply_bc


subroutine ns_case_run_simple(t)
    type(ns_case_t) :: t
    integer :: outer_i

    ! Outer loop of SIMPLE algorithm
    do outer_i = 1,t%outer_iterations

        ! Handle x-momentum
        call ns_case_x_mom(t)

        ! Handle y-momentum

        ! Solve pressure correction

        ! Update pressure

        ! Correct velocities based on pressure corrections

        ! Prepare for next iteration
        t%U_old = t%U
        t%V_old = t%V

    end do

end subroutine ns_case_run_simple

subroutine ns_case_x_mom(t)
    type(ns_case_t) :: t
    integer :: iter, i, j
    real(kind=8) :: mw, ms, me, mn

    ! Calculate x-momentum coefficients
    do i=2,t%nx
        do j=1,t%ny

            ! Mass flux terms
            mw = 0.5*t%rho*(t%U(i-1,j)+t%U(i,j))
            ms = 0.5*t%rho*(t%U(i,j-1)+t%U(i,j))
            me = 0.5*t%rho*(t%U(i,j)+t%U(i+1,j))
            me = 0.5*t%rho*(t%U(i,j)+t%U(i,j+1))

            ! Coefs
            t%AuW(i,j) = max(mw, 0.0)+t%mu*t%dy/t%dx
            t%AuS(i,j) = max(ms, 0.0)+t%mu*t%dx/t%dy
            t%AuE(i,j) = max(-me, 0.0)+t%mu*t%dy/t%dx
            t%AuN(i,j) = max(-mn, 0.0)+t%mu*t%dx/t%dy
            t%AuP(i,j) = t%AuW(i,j)+t%AuS(i,j)+t%AuE(i,j)+t%AuN(i,j)+ms-mw+mn-ms

        end do
    end do

    ! Apply x-momentum boundary conditions

    ! Iterate x-momentum
    do iter=1,10
        do i=2,t%nx
            do j=1,t%ny

                t%U(i,j) = (1.0-t%omegaU)*t%U_old(i,j)+&
                t%omegaU/t%AuP(i,j)*(t%AuW(i,j)*t%U(i-1,j)+t%AuS(i,j)*t%U(i,j-1)+t%AuE(i,j)*t%U(i+1,j)+t%AuN(i,j)*t%U(i,j+1)+&
                (t%P(i-1,j)-t%P(i,j))*t%dy)

            end do
        end do
    end do

end subroutine ns_case_x_mom

end module ns_case_m