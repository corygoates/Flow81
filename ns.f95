program navier_stokes
    ! Solves the incompressible Navier-Stokes equations on a rectangular grid with even spacing and velocity boundary conditions

    use boundary_condition_m
    use ns_case_m

    character(100) :: input_file
    type(ns_case_t) :: case

    ! Load input
    call getarg(1, input_file)
    call ns_case_load_input(case, input_file)

    ! Run case
    call ns_case_run_simple(case)

    ! Output

end program navier_stokes