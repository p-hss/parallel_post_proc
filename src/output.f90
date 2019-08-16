!----------------------------------------------------------------------------------------
! SUBROUTINE: outputting the simulation time to terminal
!----------------------------------------------------------------------------------------
subroutine final_out
    use global
    use mpi
    implicit none
    real*8 :: simulation_time
    real*8, dimension(num_procs) :: sort_time_global, comm_time_global

    call MPI_REDUCE(comm_total_time(:), comm_time_global(:), num_procs, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    call MPI_REDUCE(sort_total_time(:), sort_time_global(:), num_procs, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    if (proc_id .eq. root_process) then 
        call cpu_time(t_stop)
        simulation_time = t_stop - t_start

        !Final output 
        print*, "--------------------------------------------------------"

        write(*,'(A30, I3, A6, I3, A6)')"Simulation CPU time: ", int((simulation_time)/60), "[min]", &
            mod(int((simulation_time)),60), "[sec]"
        write(*,'(A30, I3, A6, I3, A6)')"Input reading CPU time: ", int((io_total_time)/60), "[min]", &
            mod(int((io_total_time)),60), "[sec]"
        write(*,'(A30, I3, A6, I3, A6)')"Input comm. CPU time: ", int(comm_time_global(1)/60), "[min]", &
            mod(int(comm_time_global(1)),60), "[sec]"
        write(*,'(A30, I3, A6, I3, A6)')"Input sorting CPU time: ", int(sort_time_global(1)/60), "[min]", &
            mod(int(sort_time_global(1)),60), "[sec]"
        print*, "========================================================"
    end if
end subroutine final_out

!----------------------------------------------------------------------------------------
! SUBROUTINE: writing paricle positions to file in ./data/lp_pos.dat 
!             is called after the main loop at the end of the main.f90 code.
!----------------------------------------------------------------------------------------
subroutine output_time_averages 
    use global
    implicit none
    integer :: t_pdf
    real*8, dimension(4) :: zeta_tave, zeta_tvar, xi_tave, xi_tvar !index is statistical moment
    real*8, dimension(5) :: gamma_a_tave, gamma_a_var_tave
    real*8, dimension(10) :: gamma_l_tave, gamma_l_var_tave
    real*8 :: theta_tave, theta_var_tave, phi_tave, phi_var_tave  
    real*8 :: time_average, time_variance
    real*8, dimension(N) :: eval_tave, eval_var_tave
    real*8, dimension(max_frame) :: buffer
    character(1) :: creturn

    if (proc_id .eq. root_process) then 
        zeta_tave(1) = time_average(zeta_mean(:))
        zeta_tvar(1) = time_variance(zeta_mean(:))
        zeta_tave(2) = time_average(zeta_var(:))
        zeta_tvar(2) = time_variance(zeta_var(:))
        zeta_tave(3) = time_average(zeta_skew(:))
        zeta_tvar(3) = time_variance(zeta_skew(:))
        zeta_tave(4) = time_average(zeta_kurt(:))
        zeta_tvar(4) = time_variance(zeta_kurt(:))

        xi_tave(1) = time_average(xi_mean(:))
        xi_tvar(1) = time_variance(xi_mean(:))
        xi_tave(2) = time_average(xi_var(:))
        xi_tvar(2) = time_variance(xi_var(:))
        xi_tave(3) = time_average(xi_skew(:))
        xi_tvar(3) = time_variance(xi_skew(:))
        xi_tave(4) = time_average(xi_kurt(:))
        xi_tvar(4) = time_variance(xi_kurt(:))

        buffer(:) =  eval_mean(1,:)
        eval_tave(1) = time_average(buffer(:))
        eval_var_tave(1) = time_variance(buffer(:))
        buffer(:) =  eval_mean(2,:)
        eval_tave(2) = time_average(buffer(:))
        eval_var_tave(2) = time_variance(buffer(:))
        buffer(:) =  eval_mean(3,:)
        eval_tave(3) = time_average(buffer(:))
        eval_var_tave(3) = time_variance(buffer(:))
        
        gamma_l_tave(1) = time_average(gamma_line(:,1))
        gamma_l_tave(2) = time_average(gamma_line(:,2))
        gamma_l_tave(3) = time_average(gamma_line(:,3))
        gamma_l_tave(4) = time_average(gamma_line(:,4))
        gamma_l_tave(5) = time_average(gamma_line(:,5))
        gamma_l_tave(6) = time_average(gamma_line(:,6)) !v b
        gamma_l_tave(7) = time_average(gamma_line(:,7)) !v l
        gamma_l_tave(8) = time_average(gamma_line(:,8)) !v T1
        gamma_l_tave(9) = time_average(gamma_line(:,9)) !T1 b
        gamma_l_var_tave(1) = time_variance(gamma_line(:,1))
        gamma_l_var_tave(2) = time_variance(gamma_line(:,2))
        gamma_l_var_tave(3) = time_variance(gamma_line(:,3))
        gamma_l_var_tave(4) = time_variance(gamma_line(:,4))
        gamma_l_var_tave(5) = time_variance(gamma_line(:,5))
        gamma_l_var_tave(6) = time_variance(gamma_line(:,6))
        gamma_l_var_tave(7) = time_variance(gamma_line(:,7))
        gamma_l_var_tave(8) = time_variance(gamma_line(:,8))
        gamma_l_var_tave(9) = time_variance(gamma_line(:,9))
        gamma_a_tave(1) = time_average(gamma_surf(:,1))
        gamma_a_tave(2) = time_average(gamma_surf(:,2))
        gamma_a_tave(3) = time_average(gamma_surf(:,3))
        gamma_a_tave(4) = time_average(gamma_surf(:,4))
        gamma_a_tave(5) = time_average(gamma_surf(:,5))
        gamma_a_var_tave(1) = time_variance(gamma_surf(:,1))
        gamma_a_var_tave(2) = time_variance(gamma_surf(:,2))
        gamma_a_var_tave(3) = time_variance(gamma_surf(:,3))
        gamma_a_var_tave(4) = time_variance(gamma_surf(:,4))
        gamma_a_var_tave(5) = time_variance(gamma_surf(:,5))

        buffer(:) = xi_mean(:) - 2._8*zeta_mean(:)  
        theta_tave = time_average(buffer(:))
        theta_var_tave = time_variance(buffer(:))

        buffer(:) = -(xi_mean(:) + zeta_mean(:))  
        phi_tave = time_average(buffer(:))
        phi_var_tave = time_variance(buffer(:))

        10 format(A3,A1,I8.0,A1,F10.3,A1,F10.3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
        open(unit=102, file=results_table, action='write', form='formatted')
            write(102,10) id, "&", max_lp, "&", (average_start_frame-start_frame)*dt, "&",&
                (max_frame-start_frame)*dt, "&", zeta_tave(1), "&", sqrt(zeta_tvar(1)), "&",&
                xi_tave(1), "&", sqrt(xi_tvar(1)), "\\"
        close(102)

        t_pdf=int(start_frame+(max_frame-start_frame)/2)

        11 format(A3,A1,F10.3,A1,ES10.3E2,3(A1,ES10.3E2,A1,ES10.3E2),A2)
        open(unit=102, file=zeta_pdf_table, action='write', form='formatted')
            write(102,11) id, "&",zeta_tave(1), "&", sqrt(zeta_tvar(1)), "&",&
                                  zeta_tave(2), "&", sqrt(zeta_tvar(2)), "&",&
                                  zeta_tave(3), "&", sqrt(zeta_tvar(3)), "&",&
                                  zeta_tave(4), "&", sqrt(zeta_tvar(4)), "\\"
        close(102)

        open(unit=102, file=xi_pdf_table, action='write', form='formatted')
            write(102,11) id, "&",xi_tave(1), "&", sqrt(xi_tvar(1)), "&",&
                                  xi_tave(2), "&", sqrt(xi_tvar(2)), "&",&
                                  xi_tave(3), "&", sqrt(xi_tvar(3)), "&",&
                                  xi_tave(4), "&", sqrt(xi_tvar(4)), "\\"
        close(102)

        20 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
        open(unit=102, file=eval_results_table, action='write', form='formatted')
            write(102,20) id, "&", eval_tave(1), "&", sqrt(eval_var_tave(1)),&
                "&", eval_tave(2), "&", sqrt(eval_var_tave(2)),&
                "&", eval_tave(3), "&", sqrt(eval_var_tave(3)),&
                "&", theta_tave,   "&", sqrt(theta_var_tave),&
                "&", phi_tave,     "&", sqrt(phi_var_tave), "\\"
        close(102)

        21 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
        open(unit=102, file=eval_512_results_table, action='write', form='formatted')
            write(102,21) id, "&", eval_tave(1), "&", sqrt(eval_var_tave(1)),&
                "&", eval_tave(2), "&", sqrt(eval_var_tave(2)),&
                "&", eval_tave(3), "&", sqrt(eval_var_tave(3)),"\\"
        close(102)
       

        if(mhd == 1)then
            30 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
            open(unit=102, file=line_angle_results_table, action='write', form='formatted')
                write(102,30) id, "&", gamma_l_tave(1), "&", sqrt(gamma_l_var_tave(1)),&  
                "&", gamma_l_tave(2), "&", sqrt(gamma_l_var_tave(2)),&
                "&", gamma_l_tave(3), "&", sqrt(gamma_l_var_tave(3)),&
                "&", gamma_l_tave(4), "&", sqrt(gamma_l_var_tave(4)),&
                "&", gamma_l_tave(5), "&", sqrt(gamma_l_var_tave(5)),"\\"
            close(102)

            31 format(A3,A1,ES8.3E1,A1,ES8.3E1,A1,F10.3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
            open(unit=102, file=line_angle_extended_results_table, action='write', form='formatted')
                write(102,31) id, "&", gamma_l_tave(6), "&", sqrt(gamma_l_var_tave(6)),&  
                                  "&", gamma_l_tave(7), "&", sqrt(gamma_l_var_tave(7)),&
                                  "&", gamma_l_tave(8), "&", sqrt(gamma_l_var_tave(8)),&
                                  "&", gamma_l_tave(9), "&", sqrt(gamma_l_var_tave(9)),"\\"
            close(102)
        end if

        if(mhd == 0)then
            40 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,A1,A1,A1,A2)
            open(unit=102, file=line_angle_results_table, action='write', form='formatted')
                write(102,40) id, "&", gamma_l_tave(1), "&", sqrt(gamma_l_var_tave(1)),&  
                "&", gamma_l_tave(2), "&", sqrt(gamma_l_var_tave(2)),&
                "&", gamma_l_tave(3), "&", sqrt(gamma_l_var_tave(3)),&
                "&", gamma_l_tave(4), "&", sqrt(gamma_l_var_tave(4)),&
                    "&", "-", "&", "-","\\"
            close(102)
        end if

        if(mhd == 1)then
            50 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A2)
            open(unit=102, file=surf_angle_results_table, action='write', form='formatted')
                write(102,50) id, "&", gamma_a_tave(1), "&", sqrt(gamma_a_var_tave(1)),&  
                "&", gamma_a_tave(2), "&", sqrt(gamma_a_var_tave(2)),&
                "&", gamma_a_tave(3), "&", sqrt(gamma_a_var_tave(3)),&
                "&", gamma_a_tave(4), "&", sqrt(gamma_a_var_tave(4)),&
                "&", gamma_a_tave(5), "&", sqrt(gamma_a_var_tave(5)),"\\"
            close(102)
        end if

        if(mhd == 0)then
            60 format(A3,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,F10.3,A1,ES8.3E1,A1,A1,A1,A1,A2)
            open(unit=102, file=surf_angle_results_table, action='write', form='formatted')
                write(102,60) id, "&", gamma_a_tave(1), "&", sqrt(gamma_a_var_tave(1)),&  
                "&", gamma_a_tave(2), "&", sqrt(gamma_a_var_tave(2)),&
                "&", gamma_a_tave(3), "&", sqrt(gamma_a_var_tave(3)),&
                "&", gamma_a_tave(4), "&", sqrt(gamma_a_var_tave(4)),&
                    "&", "-", "&", "-","\\"
            close(102)

        end if
    end if

end subroutine output_time_averages 

!--------------------------------------------------------------------------------------------
! FUNCTION: computing the averge of an input time series variable of length max_frame
!--------------------------------------------------------------------------------------------
function time_average(variable)
    use global
    implicit none
    integer :: i, n_steps
    real*8, dimension(max_frame) :: variable
    real*8 :: variable_sum, time_average 

    n_steps = max_frame - average_start_frame
    variable_sum = 0

    do i = 0, n_steps-1
       variable_sum = variable_sum + variable(average_start_frame + i) 
    end do

    time_average = variable_sum/n_steps
end function time_average

!----------------------------------------------------------------------------------------
! FUNCTION: computing the variance of an input time series variable  of length max_frame
!----------------------------------------------------------------------------------------
function time_variance(variable)
    use global
    implicit none
    integer :: i, n_steps
    real*8, dimension(max_frame) :: variable
    real*8 :: variable_sum, time_variance, variable_2_sum  

    n_steps = max_frame - average_start_frame
    variable_sum = 0
    variable_2_sum = 0

    do i = 0, n_steps-1
       variable_sum = variable_sum + variable(average_start_frame + i) 
       variable_2_sum = variable_2_sum + variable(average_start_frame + i)*variable(average_start_frame + i) 
    end do

    variable_sum = variable_sum/n_steps
    variable_2_sum = variable_2_sum/n_steps
    time_variance = variable_2_sum - variable_sum*variable_sum
end function time_variance
