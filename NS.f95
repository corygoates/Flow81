program navier_stokes
    ! Solves the incompressible Navier-Stokes equations on a rectangular grid with even spacing and velocity boundary conditions

    use boundary_condition_m

    ! Declare variables
    character(100) :: input_file, results_file
    real(kind=8) :: x_min, x_max, y_min, y_max
    real(kind=8) :: omegaU, omegaV, omegaP, mu, rho
    real(kind=8) :: dx, dy
    integer :: nx, ny, outer_iterations, outer_i, inner_i
    type(bc_t) :: Uw_bc, Vw_bc, Us_bc, Vs_bc, Ue_bc, Ve_bc, Un_bc, Vn_bc

    ! Array variables
    real(kind=8),allocatable,dimension(:,:) :: U, U_old, V, V_old, P ! Values
    real(kind=8),allocatable,dimension(:,:) :: AuP, AuW, AuS, AuE, AuN ! x-mom coefficients
    real(kind=8),allocatable,dimension(:,:) :: AvP, AvW, AvS, AvE, AvN ! y-mom coefficients

    ! Allocate arrays
    allocate(P(nx,ny))
    allocate(U(1:nx+1,0:ny+1))
    allocate(V(0:nx+1,1:ny+1))
    allocate(U_old(1:nx+1,0:ny+1))
    allocate(V_old(0:nx+1,1:ny+1))

    ! Get input file
    call getarg(1, input_file)

    ! Read input file
    open(1,file=input_file)
        read(1,*) ! Header

        read(1,*) ! Geometry
        read(1,*)
        read(1,*) x_min, x_max, y_min, y_max
        read(1,*)
        read(1,*) nx, ny

        read(1,*) ! Flow properties
        read(1,*)
        read(1,*) mu, rho

        read(1,*) ! Boundary conditions
        read(1,*)
        read(1,*)
        read(1,*) Uw_bc%type, Uw_bc%a, Uw_bc%b, Uw_bc%c
        read(1,*)
        read(1,*) Vw_bc%type, Vw_bc%a, Vw_bc%b, Vw_bc%c
        read(1,*)
        read(1,*)
        read(1,*) Us_bc%type, Us_bc%a, Us_bc%b, Us_bc%c
        read(1,*)
        read(1,*) Vs_bc%type, Vs_bc%a, Vs_bc%b, Vs_bc%c
        read(1,*)
        read(1,*)
        read(1,*) Ue_bc%type, Ue_bc%a, Ue_bc%b, Ue_bc%c
        read(1,*)
        read(1,*) Ve_bc%type, Ve_bc%a, Ve_bc%b, Ve_bc%c
        read(1,*)
        read(1,*)
        read(1,*) Un_bc%type, Un_bc%a, Un_bc%b, Un_bc%c
        read(1,*)
        read(1,*) Vn_bc%type, Vn_bc%a, Vn_bc%b, Vn_bc%c

        read(1,*) ! Solver Options
        read(1,*)
        read(1,*) outer_iter
        read(1,*)
        read(1,*) omegaU, omegaV, omegaP
        read(1,*)
        read(1,*) results_file

    close(1)

    ! Initialize grid
    dx = (x_max-x_min)/nx
    dy = (y_max-y_min)/ny

    ! Guess values for U, V, and P
    P = 0.0
    U = 0.0
    V = 0.0
    U_old = 0.0
    V_old = 0.0

    ! Apply boundary conditions

    ! Outer loop of SIMPLE algorithm
    do outer_i = 1,outer_iterations

        ! Solve momentum for U and V

        ! Iterate x-momentum

        ! Iterate y-momentum

        ! Solve pressure correction

        ! Update pressure

        ! Correct velocities based on pressure corrections

        ! Prepare for next iteration
        U_old = U
        V_old = V

    end do

    ! Clean up memory
    deallocate(P)
    deallocate(U)
    deallocate(V)
    deallocate(U_old)
    deallocate(V_old)


end program navier_stokes