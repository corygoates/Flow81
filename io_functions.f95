module io_functions_m
    use mesh_m
    implicit none
    
contains

subroutine display_inputs(mesh)
    type(mesh_t) :: mesh

    ! Boundary geometry
    write (*,*) "Problem boundaries"
    write (*,*) "    x_min: ", mesh%x_min
    write (*,*) "    x_max: ", mesh%x_max
    write (*,*) "    y_min: ", mesh%y_min
    write (*,*) "    y_max: ", mesh%y_max
    write (*,*)

    ! Grid params
    write (*,*) "Grid is ", mesh%n_rows, " rows by ", mesh%n_cols, " columns."
    write (*,*) "    dx: ", mesh%dx
    write (*,*) "    dy: ", mesh%dy
    write (*,*) "    Total control points: ", mesh%N_cp
    write (*,*)

    ! Boundary conditions
    write (*,*) "Boundary conditions"

    ! West
    write (*,*) "    West type: ", mesh%W_bc%type
    if (mesh%W_bc%type == "D") then
        write (*,*) "        phi=", mesh%W_bc%a, "+", mesh%W_bc%b, "x+", mesh%W_bc%c, "x^2"
    else
        write (*,*) "        dp/dn=", mesh%W_bc%a, "+", mesh%W_bc%b, "x+", mesh%W_bc%c, "x^2"
    end if
    write (*,*)

    ! South
    write (*,*) "    South type: ", mesh%S_bc%type
    if (mesh%S_bc%type == "D") then
        write (*,*) "        phi=", mesh%S_bc%a, "+", mesh%S_bc%b, "x+", mesh%S_bc%c, "x^2"
    else
        write (*,*) "        dp/dn=", mesh%S_bc%a, "+", mesh%S_bc%b, "x+", mesh%S_bc%c, "x^2"
    end if
    write (*,*)

    ! East
    write (*,*) "    East type: ", mesh%E_bc%type
    if (mesh%E_bc%type == "D") then
        write (*,*) "        phi=", mesh%E_bc%a, "+", mesh%E_bc%b, "x+", mesh%E_bc%c, "x^2"
    else
        write (*,*) "        dp/dn=", mesh%E_bc%a, "+", mesh%E_bc%b, "x+", mesh%E_bc%c, "x^2"
    end if
    write (*,*)

    ! North
    write (*,*) "    North type: ", mesh%N_bc%type
    if (mesh%N_bc%type == "D") then
        write (*,*) "        phi=", mesh%N_bc%a, "+", mesh%N_bc%b, "x+", mesh%N_bc%c, "x^2"
    else
        write (*,*) "        dp/dn=", mesh%N_bc%a, "+", mesh%N_bc%b, "x+", mesh%N_bc%c, "x^2"
    end if
    write (*,*)

end subroutine display_inputs
    
end module io_functions_m