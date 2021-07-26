program navier_stokes
    ! Solves the incompressible Navier-Stokes equations on a rectangular grid with even spacing and velocity boundary conditions

    use boundary_condition_m
    use ns_case_m

    character(100) :: input_file
    type(ns_case_t) :: case

    ! Print welcome message
    write(*,*) ""
    write(*,*) "-------------------------------------------------"
    write(*,*) "|                                               |"
    write(*,*) "|  Welcome to the Flow81 Navier-Stokes Solver!  |"
    write(*,*) "|                                               |"
    write(*,*) "|      Created by Cory Goates, (c) 2021         |"
    write(*,*) "|                                               |"
    write(*,*) "-------------------------------------------------"
    write(*,*) ""

    ! Load input
    call getarg(1, input_file)
    write(*,*) "Got input file: ", input_file
    write(*,*) "Loading input..."
    call ns_case_load_input(case, input_file)

    ! Run case
    write(*,*) ""
    write(*,*) "Running case..."
    call ns_case_run_simple(case)

    ! Output

    ! Print quit message
    write(*,*) ""
    write(*,*) "-------------------------------------------------"
    write(*,*) "|                                               |"
    write(*,*) "|         Flow81 finished successfully          |"
    write(*,*) "|                                               |"
    write(*,*) "-------------------------------------------------"
    write(*,*) ""

end program navier_stokes