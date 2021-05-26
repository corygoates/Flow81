module mesh_m

    implicit none

    type boundary_condition_t

        character :: type ! N for Neumann, D for Dirichlet
        real(kind=8) :: a, b, c ! Boundary condition formulated as a+b*x+c*x^2

    end type boundary_condition_t
    
    type mesh_t

        character(100) :: input_file
        integer :: n_rows, n_cols
        real(kind=8) :: x_min, x_max, y_min, y_max
        real(kind=8) :: dx, dy, gamma
        real(kind=8),allocatable,dimension(:,:) :: phi, x_cp, y_cp
        type(boundary_condition_t) :: W_bc, S_bc, E_bc, N_bc

    end type mesh_t

contains

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
    close(1)

    ! Parse some derived quantities
    t%dx = (t%x_max-t%x_min)/t%n_cols
    t%dy = (t%y_max-t%y_min)/t%n_rows

    ! Allocate mesh
    allocate(t%phi(t%n_rows,t%n_cols))
    allocate(t%x_cp(t%n_rows,t%n_cols))
    allocate(t%y_cp(t%n_rows,t%n_cols))

    ! Initialize control point locations
    do i = 1, t%n_rows
        do j = 1, t%n_cols
            t%x_cp(i,j) = t%x_min+(j-0.5)*t%dx
            t%y_cp(i,j) = t%y_min+(i-0.5)*t%dx
        end do
        
    end do

    write(*,*) "Successfully initialized mesh."

end subroutine mesh_load_input

subroutine set_up_matrices(mesh)
    type(mesh_t) :: mesh

end subroutine set_up_matrices
    
end module mesh_m