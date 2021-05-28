program example_convection
    implicit none

    ! Number of cells
    integer, parameter :: nx=25, ny=25

    ! Variables
    real(kind=8),dimension(nx,ny) :: AE, AW, AN, AS, Su, Sp, AP ! Upwinding coefficients
    real(kind=8),dimension(nx,ny) :: A_E, A_W, A_N, A_S, S_u, S_p, A_P ! Second-order coefficients
    real(kind=8),dimension(0:nx+1,0:ny+1) :: phi
    real(kind=8) :: x_min, y_min, x_max, y_max, u, v, dx, dy, B, rho, m_e, m_w, m_n, m_s, phi_wbc, phi_sbc
    integer :: i

    ! Welcome message
    write(*,*)
    write(*,*) "Example Convection Problem"
    write(*,*) "--------------------------"

    ! Solver options
    B = 0.0

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

    ! Mass flows
    m_e = rho*u*dy
    m_w = rho*u*dy
    m_n = rho*v*dx
    m_s = rho*v*dx

    ! Initialize all points to interior values
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
    Su = 0.0
    Sp = 0.0

    ! West boundary; Dirichlet condition
    do i = 2,ny-1

        ! Upwind
        AE(i,1) = 0.0
        AW(i,1) = 0.0
        Sp(i,1) = -m_e
        Su(i,1) = m_w*phi_wbc

        ! Second-order
        A_E(i,1) = -0.5*m_e
        A_W(i,1) = 0.0
        S_p(i,1) = -m_e
        S_u(i,1) = m_w*phi_wbc

    end do

    ! East boundary; Neumann condition
    do i = 2,ny-1

        ! Upwind
        AE(i,nx) = 0.0
        Aw(i,nx) = m_w
        Sp(i,nx) = -(m_e-m_w)

        ! Second-order
        A_E(i,nx) = 0.0
        A_W(i,nx) = 0.5*m_w
        S_p(i,nx) = -(m_e-m_w)

    end do

    ! South boundary; Dirichlet condition
    do i = 2,ny-1

        ! Upwind
        AN(i,1) = 0.0
        AS(i,1) = 0.0
        Sp(i,1) = -m_n
        Su(i,1) = m_s*phi_sbc

        ! Second-order
        A_N(i,1) = -0.5*m_n
        A_S(i,1) = 0.0
        S_p(i,1) = -m_n
        S_u(i,1) = m_s*phi_sbc

    end do

    ! North boundary; Neumann condition
    do i = 2,ny-1

        ! Upwind
        AN(i,nx) = 0.0
        AS(i,nx) = m_s
        Sp(i,nx) = -(m_n-m_s)

        ! Second-order
        A_N(i,nx) = 0.0
        A_S(i,nx) = 0.5*m_s
        S_p(i,nx) = -(m_n-m_s)

    end do

    ! Northwest corner
    Su(1,ny) = m_w*phi_wbc
    Sp(1,ny) = -m_e-(m_n-m_s)
    AW(1,ny) = 0.0
    AN(1,ny) = 0.0

    S_u(1,ny) = m_w*phi_wbc
    S_p(1,ny) = -m_e-(m_n-m_s)
    A_W(1,ny) = 0.0
    A_N(1,ny) = 0.0

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
    AW(1,nx) = 0.0
    AN(1,nx) = 0.0

    S_u(1,nx) = m_s*phi_sbc
    S_p(1,nx) = -m_n-(m_e-m_w)
    A_W(1,nx) = 0.0
    A_N(1,nx) = 0.0

    ! Northeast corner
    
end program example_convection