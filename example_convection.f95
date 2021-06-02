program example_convection
    implicit none

    ! Number of cells
    integer, parameter :: nx=50, ny=50

    ! Variables
    real(kind=8),dimension(nx,ny) :: AE, AW, AN, AS, Su, Sp, AP ! Upwinding coefficients
    real(kind=8),dimension(nx,ny) :: A_E, A_W, A_N, A_S, S_u, S_p, A_P ! Second-order coefficients
    real(kind=8),dimension(0:nx+1,0:ny+1) :: phi
    real(kind=8) :: x_min, y_min, x_max, y_max, u, v, dx, dy, B, rho, m_e, m_w, m_n, m_s, phi_wbc, phi_sbc, EQ3, EQ2, phi_new
    real(kind=8) :: max_diff, new_diff, x, y
    integer :: i, j, iterations, max_iter
    character(100) :: filename

    ! Welcome message
    write(*,*)
    write(*,*) "Example Convection Problem"
    write(*,*) "--------------------------"

    ! Solver options
    B = 0.9
    max_iter = 100000

    ! Domain
    x_min = 0.0
    x_max = 1.0
    y_min = 0.0
    y_max = 1.0

    ! Flow options
    u = 1.0
    v = 1.0
    rho = 1.0

    ! Boundary conditions
    phi_wbc = 1.0
    phi_sbc = 0.0

    ! Grid
    dx = (x_max-x_min)/nx
    dy = (y_max-y_min)/ny

    ! Initialize phi
    phi = 0.0
    phi(0,0:nx+1) = phi_sbc ! Southern boundary
    phi(0:ny+1,0) = phi_wbc ! Western boundary
    phi(0,0) = 0.5*(phi_sbc+phi_wbc) ! Southwest corner

    ! Mass flows
    m_e = rho*u*dy
    m_w = rho*u*dy
    m_n = rho*v*dx
    m_s = rho*v*dx

    ! Initialize all coefficients to interior values
    ! Upwind
    AE = 0.0
    AW = m_w
    AN = 0.0
    AS = m_s
    Su = 0.0
    Sp = 0.0

    ! Second-order
    A_E = -0.5*m_e
    A_W = 0.5*m_w
    A_N = -0.5*m_n
    A_S = 0.5*m_s
    S_u = 0.0
    S_p = 0.0

    ! West boundary; Dirichlet condition
    ! Upwind
    AE(2:ny-1,1) = 0.0
    AW(2:ny-1,1) = 0.0
    Sp(2:ny-1,1) = -m_e
    Su(2:ny-1,1) = m_w*phi_wbc

    ! Second-order
    A_E(2:ny-1,1) = -0.5*m_e
    A_W(2:ny-1,1) = 0.0
    S_p(2:ny-1,1) = -m_e
    S_u(2:ny-1,1) = m_w*phi_wbc

    ! East boundary; Neumann condition
    ! Upwind
    AE(2:ny-1,nx) = 0.0
    Aw(2:ny-1,nx) = m_w
    Sp(2:ny-1,nx) = -(m_e-m_w)

    ! Second-order
    A_E(2:ny-1,nx) = 0.0
    A_W(2:ny-1,nx) = 0.5*m_w
    S_p(2:ny-1,nx) = -(m_e-m_w)

    ! South boundary; Dirichlet condition
    ! Upwind
    AN(1,2:nx-1) = 0.0
    AS(1,2:nx-1) = 0.0
    Sp(1,2:nx-1) = -m_n
    Su(1,2:nx-1) = m_s*phi_sbc

    ! Second-order
    A_N(1,2:nx-1) = -0.5*m_n
    A_S(1,2:nx-1) = 0.0
    S_p(1,2:nx-1) = -m_n
    S_u(1,2:nx-1) = m_s*phi_sbc

    ! North boundary; Neumann condition
    ! Upwind
    AN(ny,2:nx-1) = 0.0
    AS(ny,2:nx-1) = m_s
    Sp(ny,2:nx-1) = -(m_n-m_s)

    ! Second-order
    A_N(ny,2:nx-1) = 0.0
    A_S(ny,2:nx-1) = 0.5*m_s
    S_p(ny,2:nx-1) = -(m_n-m_s)

    ! Northwest corner
    Su(ny,1) = m_w*phi_wbc
    Sp(ny,1) = -m_e-(m_n-m_s)
    AS(ny,1) = m_s
    AE(ny,1) = 0.0
    AW(ny,1) = 0.0
    AN(ny,1) = 0.0

    S_u(ny,1) = m_w*phi_wbc
    S_p(ny,1) = -m_e-(m_n-m_s)
    A_W(ny,1) = 0.0
    A_N(ny,1) = 0.0
    A_S(ny,1) = 0.5*m_s
    A_E(ny,1) = -0.5*m_e

    ! Southwest corner
    Su(1,1) = m_w*phi_wbc+m_s*phi_sbc
    Sp(1,1) = -m_e-m_n
    AE(1,1) = 0.0
    AW(1,1) = 0.0
    AN(1,1) = 0.0
    AS(1,1) = 0.0

    S_u(1,1) = m_w*phi_wbc+m_s*phi_sbc
    S_p(1,1) = -m_e-m_n
    A_W(1,1) = 0.0
    A_E(1,1) = -0.5*m_e
    A_S(1,1) = 0.0
    A_N(1,1) = -0.5*m_n

    ! Southeast corner
    Su(1,nx) = m_s*phi_sbc
    Sp(1,nx) = -m_n-(m_e-m_w)
    AW(1,nx) = m_w
    AN(1,nx) = 0.0
    AE(1,nx) = 0.0
    AS(1,nx) = 0.0

    S_u(1,nx) = m_s*phi_sbc
    S_p(1,nx) = -m_n-(m_e-m_w)
    A_E(1,nx) = 0.0
    A_W(1,nx) = 0.5*m_w
    A_S(1,nx) = 0.0
    A_N(1,nx) = -0.5*m_n

    ! Northeast corner
    Sp(ny,nx) = -(m_n-m_s)-(m_e-m_w)
    AN(ny,nx) = 0.0
    AE(ny,nx) = 0.0
    Aw(ny,nx) = m_w
    AS(ny,nx) = m_s

    S_p(ny,nx) = -(m_n-m_s)-(m_e-m_w)
    A_N(ny,nx) = 0.0
    A_S(ny,nx) = 0.5*m_s
    A_E(ny,nx) = 0.0
    A_W(ny,nx) = 0.5*m_w

    ! Calculate AP and A_P
    AP = AN+AE+AS+AW-SP+m_e-m_w+m_n-m_s
    A_P = A_N+A_E+A_S+A_W-S_P+m_e-m_w+m_n-m_s
    
    ! Iterate
    iterations = 0
    max_diff = 1.0
    write(*,*)
    write(*,*) "   Iteration      Max Correction"
    write(*,*) "--------------------------------"
    do while (.not. max_diff .eq. 0.0 .and. iterations<max_iter)

        ! Reset for iteration
        iterations = iterations+1
        max_diff = 0.0

        ! Loop through control points
        do i = 1,ny
            do j = 1,nx

                ! Upwind
                EQ3 = AE(i,j)*phi(i,j+1)+AW(i,j)*phi(i,j-1)+AN(i,j)*phi(i+1,j)+AS(i,j)*phi(i-1,j)+Su(i,j)

                ! Second order
                EQ2 = A_P(i,j)*phi(i,j)-A_E(i,j)*phi(i,j+1)-A_W(i,j)*phi(i,j-1)-A_N(i,j)*phi(i+1,j)-A_S(i,j)*phi(i-1,j)-S_u(i,j)

                ! Blend
                phi_new = EQ3/AP(i,j)-B/AP(i,j)*(EQ2-(AP(i,j)*phi(i,j)-EQ3))

                ! Check convergence
                new_diff = abs(phi_new-phi(i,j))
                if (new_diff > max_diff) then
                    max_diff = new_diff
                end if

                ! Assign
                phi(i,j) = phi_new

            end do
        end do
        write(*,*) iterations, max_diff
    end do

    ! Finalize north and east boundaries
    phi(ny+1,:) = phi(ny,:)
    phi(:,nx+1) = phi(:,nx)
    phi(ny+1,nx+1) = 0.5*(phi(ny+1,nx)+phi(ny,nx+1))

    ! Write results to file
    ! Open file
    filename = 'results.csv'
    open(1,file=filename)

        ! Header
        write(1,*) 'x,y,z,phi'

        ! Control points
        do i = 1,ny
            do j =1,nx
                x = x_min+0.5*dx*(2*j-1)
                y = y_min+0.5*dy*(2*i-1)
                write(1,*) x, ',', y, ',', 0.0, ',', phi(i,j)
            end do
        end do

        ! West and east boundaries
        do i = 1,ny
            y = y_min+0.5*dy*(2*i-1)
            write(1,*) x_min, ',', y, ',', 0.0, ',', phi(i,0) ! West
            write(1,*) x_max, ',', y, ',', 0.0, ',', phi(i,nx+1) ! East
        end do

        ! North and south boundaries
        do i = 1,nx
            x = x_min+0.5*dx*(2*i-1)
            write(1,*) x, ',', y_min, ',', 0.0, ',', phi(0,i) ! South
            write(1,*) x, ',', y_max, ',', 0.0, ',', phi(ny+1,i) ! North
        end do

        ! Corners
        write(1,*) x_min, ',', y_min, ',', 0.0, ',', phi(0,0) ! SW
        write(1,*) x_min, ',', y_max, ',', 0.0, ',', phi(ny+1,0) ! SE
        write(1,*) x_max, ',', y_min, ',', 0.0, ',', phi(0,nx+1) ! NW
        write(1,*) x_max, ',', y_max, ',', 0.0, ',', phi(ny+1,nx+1) ! NE

    close(1)

end program example_convection