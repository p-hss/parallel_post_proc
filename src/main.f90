!----------------------------------------------------------------------------------
! Program for processing the mhdt data: 
!   - computing the mean material line streching coefficients 
!----------------------------------------------------------------------------------
program main 
    use global
    use mpi

    implicit none

    integer :: i,lpi  
    integer  status(MPI_STATUS_SIZE)

    if (proc_id .eq. root_process) call cpu_time(t_start)
    !Let process 0 be the root process.
    root_process = 0
    
    !initiallize
    call MPI_INIT (ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, proc_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

    !reading input parameters
    call input_parameters

    iframe = start_frame !time frame number

    !allocating accoring to pariticle number
    call alloc
   
    !initialize material elemets
    call initialize                
    !reading the Position.***.*** files from MHDT  
    call input 

    do iframe = start_frame + 1 , max_frame!time frame number !max_time

        !reading the Position.***.*** files from MHDT  
        call input
        !integrate B:
        call b_int
        !evolve material lines in time:
        call line_evo
        !growth rates: zeta and xi
        call growth_rates(iframe)
        !strain rate statistics
        call strain_stats(iframe)
        !strain angles  histograms
        call angle_histograms(iframe)
        !fluid element deformaton
        call cauchy_green(iframe)
        
        do lpi = 1, proc_particles 
            !resetting variables 
            lp_vgr_local(lpi,1,:,:) = lp_vgr_local(lpi,2,:,:)
            B_local(lpi,1,:,:) = B_local(lpi,5,:,:)
            A_local(lpi,1,:) = A_local(lpi,2,:)
            A_length_local(lpi,1) = A_length_local(lpi,2)
            if(mhd == 1) mag_field_local(lpi,1,:) = mag_field_local(lpi,2,:)
            do i=1,3
                le_length_local(lpi,1,i) = le_length_local(lpi,2,i)
                le_local(lpi,1,:,i) = le_local(lpi,2,:,i)
            end do
        end do
    end do

    call output_time_averages 
    call final_out
    call dealloc
    call MPI_FINALIZE (ierr)
end program main 
