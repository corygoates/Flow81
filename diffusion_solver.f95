program diffusion_solver
    ! Solves the diffusion problem on a rectangular grid with even spacing

    use mesh_m
    use io_functions_m

    character(100) :: input_file, results_file
    type(mesh_t) :: mesh

    ! Welcome message
    write (*,*)
    write (*,*) "Flow81 diffusion solver"
    write (*,*) "-----------------------"
    write (*,*)

    ! Read in input file
    call getarg(1, input_file)
    call mesh_load_input(mesh, input_file)

    ! Display inputs
    call display_inputs(mesh)

    ! Set up problem
    call mesh_set_up_matrices(mesh)

    ! Solve
    call mesh_sor(mesh)
    write(*,*) mesh%phi(mesh%N_cp/2+1)
    write(*,*) abs(mesh%phi(mesh%N_cp/2+1)-0.25d0)/0.25d0

    ! Write results to csv
    call getarg(2, results_file)
    call mesh_write_results_to_csv(mesh, results_file)

    ! Clean up
    call mesh_deallocate(mesh)

end program diffusion_solver