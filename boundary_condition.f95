module boundary_condition_m

    implicit none

    type bc_t

        character :: type ! N for Neumann, D for Dirichlet
        real(kind=8) :: a, b, c ! Boundary condition formulated as a+b*x+c*x^2

    end type bc_t

contains

real(kind=8) function bc_get_value(t, x)
    type(bc_t), intent(in) :: t
    real(kind=8), intent(in) :: x

    bc_get_value = t%a+t%b*x+t%c*x**2

end function bc_get_value

end module boundary_condition_m