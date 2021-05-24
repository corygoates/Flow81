program diffusion_solver
    ! Solves the diffusion problem on a rectangular grid with even spacing

    integer :: M, N
    real(kind=8) :: x_min, x_max, y_min, y_max
    real(kind=8) :: dx, dy, gamma
    character(len=32) :: input_file, misc
    character :: w_bc_type, s_bc_type, e_bc_type, n_bc_type
    real(kind=8) :: w_bc_val, s_bc_val, e_bc_val, n_bc_val

    ! Welcome message
    write (*,*)
    write (*,*) "Flow81 diffusion solver"
    write (*,*) "-----------------------"
    write (*,*)

    ! Read in input file
    call getarg(1, input_file)
    open(1,file=input_file)
        read(1,*)
        read(1,*)
        read(1,*) x_min, x_max, y_min, y_max
        read(1,*)
        read(1,*) M, N
        read(1,*)
        read(1,*) gamma
        read(1,*)
        read(1,*)
        read(1,*) w_bc_type, w_bc_val
        read(1,*)
        read(1,*) s_bc_type, s_bc_val
        read(1,*)
        read(1,*) e_bc_type, e_bc_val
        read(1,*)
        read(1,*) n_bc_type, n_bc_val
    close(1)

    ! Parse some derived quantities
    dx = (x_max-x_min)/N
    dy = (y_max-y_min)/M

    ! Display inputs
    write (*,*) "Problem boundaries"
    write (*,*) "    x_min: ", x_min
    write (*,*) "    x_max: ", x_max
    write (*,*) "    y_min: ", y_min
    write (*,*) "    y_max: ", y_max
    write (*,*)
    write (*,*) "Grid is ", M, " rows by ", N, " columns."
    write (*,*) "    dx: ", dx
    write (*,*) "    dy: ", dy

end program diffusion_solver