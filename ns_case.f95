module ns_case_m

    use boundary_condition_m

    implicit none

    type ns_case_t

        character(200) :: input_file, results_file
        real(kind=8) :: x_min, x_max, y_min, y_max
        real(kind=8) :: omegaU, omegaV, omegaP, mu, rho
        real(kind=8) :: dx, dy
        integer :: nx, ny, outer_i, outer_iterations
        type(bc_t) :: Uw_bc, Vw_bc, Us_bc, Vs_bc, Ue_bc, Ve_bc, Un_bc, Vn_bc
        real(kind=8),allocatable,dimension(:,:) :: U, U_old, V, V_old, P ! Values
        real(kind=8),allocatable,dimension(:,:) :: P_prime ! Corrections
        real(kind=8),allocatable,dimension(:,:) :: u_out, v_out, P_out, v_mag_out ! Output values
        real(kind=8),allocatable,dimension(:,:) :: AuP, AuW, AuS, AuE, AuN ! x-mom coefficients
        real(kind=8),allocatable,dimension(:,:) :: AvP, AvW, AvS, AvE, AvN ! y-mom coefficients
        real(kind=8),allocatable,dimension(:,:) :: AP, AW, AS, AE, AN, S ! pressure coefficients
        real(kind=8),allocatable,dimension(:) :: x_u, x_v, y_u, y_v ! Velocity control point locations

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

    allocate(t%P_prime(1:t%nx,1:t%ny))

    allocate(t%u_out(0:t%nx,0:t%ny))
    allocate(t%v_out(0:t%nx,0:t%ny))
    allocate(t%P_out(0:t%nx,0:t%ny))
    allocate(t%v_mag_out(0:t%nx,0:t%ny))

    allocate(t%AuW(2:t%nx,1:t%ny))
    allocate(t%AuS(2:t%nx,1:t%ny))
    allocate(t%AuE(2:t%nx,1:t%ny))
    allocate(t%AuN(2:t%nx,1:t%ny))
    allocate(t%AuP(2:t%nx,1:t%ny))

    allocate(t%AvW(1:t%nx,2:t%ny))
    allocate(t%AvS(1:t%nx,2:t%ny))
    allocate(t%AvE(1:t%nx,2:t%ny))
    allocate(t%AvN(1:t%nx,2:t%ny))
    allocate(t%AvP(1:t%nx,2:t%ny))

    allocate(t%AW(1:t%nx,1:t%ny))
    allocate(t%AS(1:t%nx,1:t%ny))
    allocate(t%AE(1:t%nx,1:t%ny))
    allocate(t%AN(1:t%nx,1:t%ny))
    allocate(t%AP(1:t%nx,1:t%ny))
    allocate(t%S(1:t%nx,1:t%ny))

    allocate(t%x_u(1:t%nx+1))
    allocate(t%y_u(0:t%ny+1))
    allocate(t%x_v(0:t%nx+1))
    allocate(t%y_v(1:t%ny+1))

    allocated = 1

end subroutine ns_case_allocate


subroutine ns_case_deallocate(t)
    type(ns_case_t) :: t

    ! Clean up memory
    deallocate(t%P)
    deallocate(t%P_prime)
    deallocate(t%U)
    deallocate(t%V)
    deallocate(t%U_old)
    deallocate(t%V_old)

    deallocate(t%u_out)
    deallocate(t%v_out)
    deallocate(t%P_out)
    deallocate(t%v_mag_out)

    deallocate(t%AuW)
    deallocate(t%AuS)
    deallocate(t%AuE)
    deallocate(t%AuN)
    deallocate(t%AuP)

    deallocate(t%AvW)
    deallocate(t%AvS)
    deallocate(t%AvE)
    deallocate(t%AvN)
    deallocate(t%AvP)

    deallocate(t%AW)
    deallocate(t%AS)
    deallocate(t%AE)
    deallocate(t%AN)
    deallocate(t%AP)
    deallocate(t%S)

    deallocate(t%x_u)
    deallocate(t%y_u)
    deallocate(t%x_v)
    deallocate(t%y_v)

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
        read(1,"(a)") t%results_file

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
    t%P_prime = 0.0
    t%U = 0.0
    t%V = 0.0

    ! Set up control point locations
    call ns_case_setup_control_points(t)

    ! Initialize boundary conditions
    call ns_case_set_dirichlet_bc(t)

    ! Initialize old values
    t%U_old = t%U
    t%V_old = t%V

    ! Initialize coefficients
    t%AP = 1.0
    t%AuP = 1.0
    t%AvP = 1.0

end subroutine ns_case_initialize


subroutine ns_case_setup_control_points(t)
    type(ns_case_t) :: t
    integer :: i,j

    ! U
    ! Get y location
    do j=0,t%ny+1
        if (j.eq.0) then
            t%y_u(j) = t%y_min
        else if (j.eq.t%ny+1) then
            t%y_u(j) = t%y_max
        else
            t%y_u(j) = (j-0.5)*t%dy
        end if
    end do

    ! Get x location
    do i=1,t%nx+1
        t%x_u(i) = (i-1.0)*t%dx
    end do

    ! V
    ! Get y location
    do j=1,t%ny+1
        t%y_v(j) = (j-1.0)*t%dy
    end do

    ! Get x location
    do i=0,t%nx+1
        if (i .eq. 0) then
            t%x_v(i) = t%x_min
        else if (i .eq. t%nx+1) then
            t%x_v(i) = t%x_max
        else
            t%x_v(i) = (i-0.5)*t%dx
        end if
    end do

end subroutine ns_case_setup_control_points


subroutine ns_case_set_dirichlet_bc(t)
    type(ns_case_t) :: t
    integer :: i,j

    ! U boundary conditions
    ! West
    if (t%Uw_bc%type .eq. 'D') then
        do j=0,t%ny+1
            t%U(1,j) = bc_get_value(t%Uw_bc, t%y_u(j))
        end do
    end if

    ! East
    if (t%Ue_bc%type .eq. 'D') then
        do j=0,t%ny+1
            t%U(t%nx+1,j) = bc_get_value(t%Ue_bc, t%y_u(j))
        end do
    end if

    ! South
    if (t%Us_bc%type .eq. 'D') then
        do i=1,t%nx+1
            t%U(i,0) = bc_get_value(t%Us_bc, t%x_u(i))
        end do
    end if

    ! North
    if (t%Un_bc%type .eq. 'D') then
        do i=1,t%nx+1
            t%U(i,t%ny+1) = bc_get_value(t%Un_bc, t%x_u(i))
        end do
    end if

    ! V boundary conditions
    ! West
    if (t%Vw_bc%type .eq. 'D') then
        do j=1,t%ny+1
            t%V(0,j) = bc_get_value(t%Vw_bc, t%y_v(j))
        end do
    end if

    ! East
    if (t%Ve_bc%type .eq. 'D') then
        do j=1,t%ny+1
            t%V(t%nx+1,j) = bc_get_value(t%Ve_bc, t%y_v(j))
        end do
    end if


    ! South
    if (t%Vs_bc%type .eq. 'D') then
        do i=0,t%nx+1
            t%V(i,1) = bc_get_value(t%Vs_bc, t%x_v(i))
        end do
    end if

    ! North
    if (t%Vn_bc%type .eq. 'D') then
        do i=0,t%nx+1
            t%V(i,t%ny+1) = bc_get_value(t%Vn_bc, t%x_v(i))
        end do
    end if

end subroutine ns_case_set_dirichlet_bc


subroutine ns_case_set_neumann_u_bc(t)
    type(ns_case_t) :: t
    integer :: i,j

    ! West
    if (t%Uw_bc%type .eq. 'N') then
        do j=1,t%ny
            t%U(1,j) = t%U(2,j)-bc_get_value(t%Uw_bc, t%y_u(j))*t%dx
        end do
    end if

    ! South
    if (t%Us_bc%type .eq. 'N') then
        do i=2,t%nx
            t%U(i,0) = t%U(i,1)-0.5*bc_get_value(t%Us_bc, t%x_u(i))*t%dy
        end do
    end if

    ! East
    if (t%Ue_bc%type .eq. 'N') then
        do j=1,t%ny
            t%U(t%nx+1,j) = t%U(t%nx,j)+bc_get_value(t%Ue_bc, t%y_u(j))*t%dx
        end do
    end if

    ! North
    if (t%Un_bc%type .eq. 'N') then
        do i=2,t%nx
            t%U(i,t%ny+1) = t%U(i,t%ny+1)+0.5*bc_get_value(t%Un_bc, t%x_u(i))*t%dy
        end do
    end if

end subroutine ns_case_set_neumann_u_bc


subroutine ns_case_set_neumann_v_bc(t)
    type(ns_case_t) :: t
    integer :: i,j

    ! West
    if (t%Uw_bc%type .eq. 'N') then
        do j=2,t%ny
            t%V(0,j) = t%V(1,j)-0.5*bc_get_value(t%Uw_bc, t%y_v(j))*t%dx
        end do
    end if

    ! South
    if (t%Us_bc%type .eq. 'N') then
        do i=1,t%nx
            t%V(i,1) = t%V(i,2)-bc_get_value(t%Vs_bc, t%x_v(i))*t%dy
        end do
    end if

    ! East
    if (t%Ue_bc%type .eq. 'N') then
        do j=2,t%ny
            t%V(t%nx+1,j) = t%V(t%nx,j)+0.5*bc_get_value(t%Vw_bc, t%y_v(j))*t%dx
        end do 
    end if

    ! North
    if (t%Un_bc%type .eq. 'N') then
        do i=1,t%nx
            t%V(i,t%ny+1) = t%V(i,t%nx)+bc_get_value(t%Vn_bc, t%x_v(i))*t%dy
        end do
    end if

end subroutine ns_case_set_neumann_v_bc


subroutine ns_case_run_simple(t)
    type(ns_case_t) :: t
    integer :: outer_i

    ! Initialize display
    write(*,*)
    write(*,*) " Outer Iteration         Mass Imbalance"
    write(*,*) "-------------------------------------------------"

    ! Outer loop of SIMPLE algorithm
    do outer_i = 0,t%outer_iterations
        t%outer_i = outer_i

        ! Handle x-momentum
        call ns_case_x_mom(t)

        ! Handle y-momentum
        call ns_case_y_mom(t)

        ! Pressure correction
        call ns_case_pressure_corr(t)

        ! Correct velocities based on pressure corrections
        call ns_case_velocity_corr(t)

        ! Prepare for next iteration
        t%U_old = t%U
        t%V_old = t%V

        ! Check mass imbalance
        call ns_case_calc_mass_imbal(t)

        ! Update display
        if (modulo(outer_i, 50) .eq. 0) then
            write(*,*) t%outer_i, sqrt(sum(t%S**2))
        end if

    end do

    ! Write results
    call ns_case_write_results(t)

end subroutine ns_case_run_simple


subroutine ns_case_x_mom(t)
    type(ns_case_t) :: t
    integer :: iter, i, j
    real(kind=8) :: mw, ms, me, mn

    ! Calculate x-momentum coefficients
    do i=2,t%nx
        do j=1,t%ny

            ! Mass flux terms
            mw = 0.5*t%rho*t%dy*(t%U_old(i-1,j)+t%U_old(i,j))
            ms = 0.5*t%rho*t%dx*(t%V_old(i-1,j)+t%V_old(i,j))
            me = 0.5*t%rho*t%dy*(t%U_old(i,j)+t%U_old(i+1,j))
            mn = 0.5*t%rho*t%dx*(t%V_old(i-1,j+1)+t%V_old(i,j+1))

            ! Coefs
            t%AuW(i,j) = max(mw, 0.0)+t%mu*t%dy/t%dx
            t%AuE(i,j) = max(-me, 0.0)+t%mu*t%dy/t%dx

            ! Southern boundary
            if (j .eq. 1) then
                t%AuS(i,j) = max(ms, 0.0)+2.0*t%mu*t%dx/t%dy
                t%AuN(i,j) = max(-mn, 0.0)+t%mu*t%dx/t%dy

            ! Northern boundary
            else if (j .eq. t%ny) then
                t%AuS(i,j) = max(ms, 0.0)+t%mu*t%dx/t%dy
                t%AuN(i,j) = max(-mn, 0.0)+2.0*t%mu*t%dx/t%dy

            ! Center cells
            else
                t%AuS(i,j) = max(ms, 0.0)+t%mu*t%dx/t%dy
                t%AuN(i,j) = max(-mn, 0.0)+t%mu*t%dx/t%dy
            end if

            t%AuP(i,j) = t%AuW(i,j)+t%AuS(i,j)+t%AuE(i,j)+t%AuN(i,j)!+me-mw+mn-ms

        end do
    end do

    ! Iterate x-momentum
    t%AuP = t%AuP/t%omegaU
    do iter=1,10
        do i=2,t%nx
            do j=1,t%ny

                t%U(i,j) = (1.0-t%omegaU)*t%U_old(i,j) + &
                           1.0/t%AuP(i,j)*(&
                           t%AuW(i,j)*t%U(i-1,j) + &
                           t%AuS(i,j)*t%U(i,j-1) + &
                           t%AuE(i,j)*t%U(i+1,j) + &
                           t%AuN(i,j)*t%U(i,j+1) + &
                           (t%P(i-1,j)-t%P(i,j))*t%dy)

            end do
        end do

        ! Update Neumann boundary condition
        call ns_case_set_neumann_u_bc(t)

    end do

    ! Ensure mass conservation on boundaries
    call ns_case_ensure_mass_cons(t)

end subroutine ns_case_x_mom


subroutine ns_case_y_mom(t)
    type(ns_case_t) :: t
    integer :: iter, i, j
    real(kind=8) :: mw, ms, me, mn

    ! Calculate y-momentum coefficients
    do i=1,t%nx
        do j=2,t%ny

            ! Mass flux terms
            mw = 0.5*t%rho*t%dy*(t%U_old(i,j-1)+t%U_old(i,j))
            ms = 0.5*t%rho*t%dx*(t%V_old(i,j-1)+t%V_old(i,j))
            me = 0.5*t%rho*t%dy*(t%U_old(i+1,j-1)+t%U_old(i+1,j))
            mn = 0.5*t%rho*t%dx*(t%V_old(i,j)+t%V_old(i,j+1))

            ! Coefs
            t%AvS(i,j) = max(ms, 0.0)+t%mu*t%dx/t%dy
            t%AvN(i,j) = max(-mn, 0.0)+t%mu*t%dx/t%dy

            ! Western boundary
            if (i .eq. 1) then
                t%AvW(i,j) = max(mw, 0.0)+2.0*t%mu*t%dy/t%dx
                t%AvE(i,j) = max(-me, 0.0)+t%mu*t%dy/t%dx

            ! Eastern boundary
            else if (i .eq. t%nx) then
                t%AvW(i,j) = max(mw, 0.0)+t%mu*t%dy/t%dx
                t%AvE(i,j) = max(-me, 0.0)+2.0*t%mu*t%dy/t%dx

            ! Center cells
            else
                t%AvW(i,j) = max(mw, 0.0)+t%mu*t%dy/t%dx
                t%AvE(i,j) = max(-me, 0.0)+t%mu*t%dy/t%dx

            end if

            t%AvP(i,j) = t%AvW(i,j)+t%AvS(i,j)+t%AvE(i,j)+t%AvN(i,j)!+me-mw+mn-ms

        end do
    end do


    ! Iterate y-momentum
    t%AvP = t%AvP/t%omegaV
    do iter=1,10
        do i=1,t%nx
            do j=2,t%ny

                t%V(i,j) = (1.0-t%omegaV)*t%V_old(i,j) + &
                           1.0/t%AvP(i,j)*( &
                           t%AvW(i,j)*t%V(i-1,j) + &
                           t%AvS(i,j)*t%V(i,j-1) + &
                           t%AvE(i,j)*t%V(i+1,j) + &
                           t%AvN(i,j)*t%V(i,j+1) + &
                           (t%P(i,j-1)-t%P(i,j))*t%dx)

            end do
        end do

        ! Update Neumann boundary condition
        call ns_case_set_neumann_v_bc(t)

    end do

    ! Ensure mass conservation on boundaries
    call ns_case_ensure_mass_cons(t)

end subroutine ns_case_y_mom


subroutine ns_case_pressure_corr(t)
    type(ns_case_t) :: t
    integer :: iter, i, j

    ! Calculate coefficients
    do i=1,t%nx
        do j=1,t%ny
            t%AW(i,j) = t%rho*t%dy**2/t%AuP(i,j)
            t%AS(i,j) = t%rho*t%dx**2/t%AvP(i,j)
            t%AE(i,j) = t%rho*t%dy**2/t%AuP(i+1,j)
            t%AN(i,j) = t%rho*t%dx**2/t%AvP(i,j+1)
        end do
    end do

    ! Set boundary coefficients
    t%AW(1,:) = 0.0
    t%AS(:,1) = 0.0
    t%AE(t%nx,:) = 0.0
    t%AN(:,t%ny) = 0.0

    ! Calculate cell center coefficients
    t%AP = t%AW+t%AS+t%AE+t%AN
    t%AP(1,1) = 1.0e30 ! Reference pressure cell

    ! Calculate mass source
    call ns_case_calc_mass_imbal(t)

    ! Iterate
    t%P_prime = 0.0
    do iter=1,1000
        do i=1,t%nx
            do j=1,t%ny
                t%P_prime(i,j) = t%P_prime(i,j) + &
                                 1.7/t%AP(i,j)*( &
                                 t%AW(i,j)*t%P_prime(i-1,j) + &
                                 t%AS(i,j)*t%P_prime(i,j-1) + &
                                 t%AE(i,j)*t%P_prime(i+1,j) + &
                                 t%AN(i,j)*t%P_prime(i,j+1) - &
                                 t%S(i,j) - t%AP(i,j)*t%P_prime(i,j))
            end do
        end do
    end do

    t%P = t%P+t%omegaP*t%P_prime

end subroutine ns_case_pressure_corr


subroutine ns_case_velocity_corr(t)
    type(ns_case_t) :: t
    integer :: i, j

    ! Update x-velocity
    do i=2,t%nx
        do j=1,t%ny
            t%U(i,j) = t%U(i,j) + t%dy/t%AuP(i,j)*(t%P_prime(i-1,j)-t%P_prime(i,j))
        end do
    end do

    ! Update y-velocity
    do i=1,t%nx
        do j=2,t%ny
            t%V(i,j) = t%V(i,j) + t%dx/t%AvP(i,j)*(t%P_prime(i,j-1)-t%P_prime(i,j))
        end do
    end do

end subroutine ns_case_velocity_corr


subroutine ns_case_ensure_mass_cons(t)
    type(ns_case_t) :: t
    real(kind=8) :: total_mass_imbal, forced_mass_imbal, corr
    integer :: N_neu, N_dir

    ! Calculate total mass imbalance
    total_mass_imbal = sum(t%U(1,:)-t%U(t%nx+1,:))*t%dy+sum(t%V(:,1)-t%V(:,t%ny+1))*t%dx

    ! Calculate forced mass imbalance (due to user-prescribed Dirichlet boundary conditions)
    forced_mass_imbal = 0.0
    if (t%Uw_bc%type .eq. 'D') then
        forced_mass_imbal = forced_mass_imbal+sum(t%U(1,:))*t%dy
    end if
    if (t%Ue_bc%type .eq. 'D') then
        forced_mass_imbal = forced_mass_imbal-sum(t%U(t%nx+1,:))*t%dy
    end if
    if (t%Vs_bc%type .eq. 'D') then
        forced_mass_imbal = forced_mass_imbal+sum(t%V(:,1))*t%dx
    end if
    if (t%Vn_bc%type .eq. 'D') then
        forced_mass_imbal = forced_mass_imbal-sum(t%V(:,t%ny+1))*t%dx
    end if

    ! Determine number of Neumann boundaries which can influence mass conservation
    N_neu = 0
    if (t%Uw_bc%type .eq. 'N') then
        N_neu = N_neu + 1
    end if
    if (t%Ue_bc%type .eq. 'N') then
        N_neu = N_neu + 1
    end if
    if (t%Vs_bc%type .eq. 'N') then
        N_neu = N_neu + 1
    end if
    if (t%Vn_bc%type .eq. 'N') then
        N_neu = N_neu + 1
    end if

    ! Determine which type of correction is needed
    ! If there is forced mass imbalance (i.e. finite mass flux
    ! specified by user), then the non-user-specified boundaries
    ! can be corrected multiplicatively. Otherwise, an 
    ! additive correction needs to be made.
    if (abs(forced_mass_imbal)>1e-12 .and. abs(forced_mass_imbal-total_mass_imbal)>1e-12) then

        ! Mass imbalance can be fixed multiplicatively
        ! Calculate correction
        corr = -forced_mass_imbal/(total_mass_imbal-forced_mass_imbal)

        ! Apply correction
        if (t%Uw_bc%type .eq. 'N') then
            t%U(1,:) = t%U(1,:)*corr
        end if
        if (t%Ue_bc%type .eq. 'N') then
            t%U(t%nx+1,:) = t%U(t%nx+1,:)*corr
        end if
        if (t%Vs_bc%type .eq. 'N') then
            t%V(:,1) = t%V(:,1)*corr
        end if
        if (t%Vn_bc%type .eq. 'N') then
            t%V(:,t%ny+1) = t%V(:,t%ny+1)*corr
        end if

    else

        ! Mass imbalance needs to be fixed additively
        ! Calculate correction
        corr = (total_mass_imbal-forced_mass_imbal)/real(N_neu)

        ! Apply correction
        if (t%Uw_bc%type .eq. 'N') then
            t%U(1,:) = t%U(1,:)-corr
        end if
        if (t%Ue_bc%type .eq. 'N') then
            t%U(t%nx+1,:) = t%U(t%nx+1,:)+corr
        end if
        if (t%Vs_bc%type .eq. 'N') then
            t%V(:,1) = t%V(:,1)-corr
        end if
        if (t%Vn_bc%type .eq. 'N') then
            t%V(:,t%ny+1) = t%V(:,t%ny+1)+corr
        end if

    end if


end subroutine ns_case_ensure_mass_cons


subroutine ns_case_calc_mass_imbal(t)
    type(ns_case_t) :: t
    integer :: i, j

    ! Calculate mass sources
    do i=1,t%nx
        do j=1,t%ny
            t%S(i,j) = t%rho*((t%U(i+1,j)-t%U(i,j))*t%dy+(t%V(i,j+1)-t%V(i,j))*t%dx)
        end do
    end do

end subroutine ns_case_calc_mass_imbal


subroutine ns_case_write_results(t)
    type(ns_case_t) :: t
    integer :: i,j
    real(kind=8) :: x,y

    ! Interpolate values to scalar cell corners
    ! Center values
    do i=1,t%nx-1
        do j=1,t%ny-1
            t%P_out(i,j) = 0.25*(t%P(i,j)+t%P(i+1,j)+t%P(i,j+1)+t%P(i+1,j+1))
            t%u_out(i,j) = 0.5*(t%U(i+1,j)+t%U(i+1,j+1))
            t%v_out(i,j) = 0.5*(t%V(i,j+1)+t%V(i+1,j+1))
        end do
    end do

    ! West boundary
    t%P_out(0,1:t%ny-1) = 0.5*(t%P(1,1:t%ny-1)+t%P(1,2:t%ny))
    t%u_out(0,1:t%ny-1) = 0.5*(t%U(1,1:t%ny-1)+t%U(1,2:t%ny))
    t%v_out(0,1:t%ny-1) = t%V(0,2:t%ny)

    ! East boundary
    t%P_out(t%nx,1:t%ny-1) = 0.5*(t%P(t%nx,1:t%ny-1)+t%P(t%nx,2:t%ny))
    t%u_out(t%nx,1:t%ny-1) = 0.5*(t%U(t%nx+1,1:t%ny-1)+t%U(t%nx+1,2:t%ny))
    t%v_out(t%nx,1:t%ny-1) = t%V(t%nx+1,2:t%ny)

    ! South boundary
    t%P_out(1:t%nx-1,0) = 0.5*(t%P(1:t%nx-1,1)+t%P(2:t%nx,1))
    t%u_out(1:t%nx-1,0) = t%U(2:t%nx,0)
    t%v_out(1:t%nx-1,0) = 0.5*(t%V(1:t%nx-1,1)+t%V(2:t%nx,1))

    ! North boundary
    t%P_out(1:t%nx-1,t%ny) = 0.5*(t%P(1:t%nx-1,t%ny)+t%P(2:t%nx,t%ny))
    t%u_out(1:t%nx-1,t%ny) = t%U(2:t%nx,t%ny+1)
    t%v_out(1:t%nx-1,t%ny) = 0.5*(t%V(1:t%nx-1,t%ny+1)+t%V(2:t%nx,t%ny+1))

    ! Northwest corner
    t%P_out(0,t%ny) = t%P(1,t%ny)

    ! Southwest corner
    t%P_out(0,0) = t%P(1,1)

    ! Southeast corner
    t%P_out(t%nx,0) = t%P(t%nx,1)

    ! Northeast corner
    t%P_out(t%nx,t%ny) = t%P(t%nx,t%ny)

    ! Calculate magnitudes
    t%v_mag_out = sqrt(t%u_out**2+t%v_out**2)

    ! Write to file
    write(*,*)
    write(*,*) 'Writing results to ', t%results_file
    open(1, file=t%results_file)
    write(1,*) 'x,y,z,P,u,v,w,V'
    do i=0,t%nx
        do j=0,t%ny
            x = i*t%dx
            y = j*t%dy
            write(1,100) x,y,0.0,t%P_out(i,j),t%u_out(i,j),t%v_out(i,j),0.0,t%v_mag_out(i,j)
        end do
    end do

    close(1)

    100   format(f8.4,', ',f8.4,', ',f8.4,', ',f8.4,', ',f8.5,', ',f8.5,', ',f8.5,', ',f8.5)

end subroutine ns_case_write_results

end module ns_case_m