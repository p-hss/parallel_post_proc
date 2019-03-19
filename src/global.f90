module global
    implicit none
    !dimensions 
    integer, parameter                      :: n = 3 !number of spatial dimensions
    integer                                 :: iframe !loop indeces
    integer                                 :: max_i,max_j
    integer                                 :: max_lp 
    integer                                 :: nioproc !number of output files per io process
    character*3                             :: id 
    integer                                 :: mhd
    real*8, parameter                       :: pi = 3.14159265359
    !time
    integer                                 :: start_frame  
    integer                                 :: max_frame
    integer                                 :: corr_start, corr_end !dimension(max_lp,max_frame)
    integer                                 :: corr_start_time !kolmo time 
    integer                                 :: corr_end_time !kolmo time 
    integer                                 :: histo_time !kolmo time frame for histogram evaluation
    integer                                 :: histo_frame !kolmo time frame for histogram evaluation
    integer                                 :: average_start_time !kolmo time 
    integer                                 :: average_start_frame
    !MPI variables                                                     
    integer                                 :: num_procs, proc_id, proc_particles
    integer                                 :: root_process, ierr
    !time variables                                                     
    real*8                                  :: dt, dt_init, dt_fin 
    real*8                                  :: t_diss 
    real*8                                  :: t_start, t_stop
    real*8                                  :: io_t_start, io_t_stop, io_total_time 
    real*8, dimension(:), allocatable       :: sort_t_start, sort_t_stop, sort_total_time 
    real*8, dimension(:), allocatable       :: comm_t_start, comm_t_stop, comm_total_time 
    real*8                                  :: t_kolmo 
    !histograms
    real*8                                  :: bin_size  
    real*8                                  :: bin_size_w  !for the cauchy green histogram 
    real*8                                  :: bin_size_angle !angle histogram
    !line elemets
    integer*4                               :: lp_number 
    integer, parameter                      :: seed = 8264
    real*8                                  :: le3_sign 
    real*8, dimension(:,:,:,:), allocatable :: lp_vgr_local 
    real*8, dimension(:,:,:), allocatable   :: mag_field_local
    real*8, dimension(:,:,:,:), allocatable :: le_local 
    real*8, dimension(:,:,:), allocatable   :: le_initial_local
    real*8, dimension(:,:,:), allocatable   :: le_length_local 
    real*8, dimension(:,:,:,:), allocatable :: B_local
    real*8, dimension(:,:,:), allocatable   :: A_local
    real*8, dimension(:,:), allocatable     :: A_length_local 
    real*8, dimension(:), allocatable       :: zeta_mean,zeta_var,zeta_skew,zeta_kurt,&
                                               xi_mean, xi_var,xi_skew,xi_kurt 
    real*8, dimension(:,:), allocatable     :: eval_mean, eval_var ! strain rates dimension(N,max_frame)
    real*8, dimension(:,:), allocatable     :: gamma_line, gamma_line_var, gamma_surf, gamma_surf_var
    real*8, dimension(:), allocatable       :: theta_var, phi_var

    character*26                            :: readdir
    character*13                            :: outputdir
    character*140                           :: file_out, vel_file, line_evo_file, surf_evo_file,&
                                               angle_file, cg_file, strain_file,&
                                               line_histo_file, surf_histo_file,&
                                               cg_histo_file, mhd_angle_histo_file, angle_histo_file,&
                                               corr_file, mean_length_file,&
                                               mhd_line_evo_file
    character*100                            :: results_table, line_angle_results_table,&
                                               surf_angle_results_table, eval_results_table, pos_file,&
                                               zeta_pdf_table, xi_pdf_table

    integer, parameter                      :: verbose = 0
end module global

subroutine  input_parameters
    use global
    use mpi
    implicit none

    if (proc_id .eq. root_process) then 
        10 format(A40, A8)
        20 format(A40, F8.3)
        30 format(A40, I8)

        read(*, *) max_lp 
        read(*, *) nioproc 
        read(*, *) id 
        read(*, *) mhd 
        read(*, *) start_frame 
        read(*, *) max_frame
        read(*, *) corr_start_time
        read(*, *) corr_end_time
        read(*, *) average_start_time 
        read(*, *) histo_time
        read(*, *) t_kolmo
        read(*, *) bin_size 
        read(*, *) bin_size_w
        read(*, *) bin_size_angle

        proc_particles=max_lp/num_procs

        print*, "========================================================"
        print*, "          Line Element Deformation Simulation" 
        print*, "========================================================"
        write(*,10) "simulation ID: ", id 
        write(*,30) "MHD: ",  mhd 
        print*, "--------------------------------------------------------"
        write(*,20) "Kolmogorov time scale: ", t_kolmo
        write(*,30) "number of lagrange particles: ", max_lp 
        write(*,30) "number of PEs: ", num_procs
        write(*,30) "number of particles per PE: ", proc_particles
        write(*,30) "number of input files/time frame: ", nioproc
        print*, "--------------------------------------------------------"
        write(*,30) "MHDT start frame: ", start_frame 
        write(*,30) "MHDT end frame: ", max_frame
        write(*,30) "histogram time: ", histo_time 
        write(*,30) "start time for correlation measure: ", corr_start_time
        write(*,30) "end time for correlation measure: ", corr_end_time
        write(*,30) "start time for averages: ", average_start_time 
        print*, "--------------------------------------------------------"
        write(*,20) "zeta and xi histogram bin size: ", bin_size 
        write(*,20) "cauchy green histogram bin size: ", bin_size_w 
        write(*,20) "angle histograms bin size: ", bin_size_angle 
        print*, "========================================================"
    end if

    call MPI_BCAST (max_lp, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (nioproc, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (proc_particles, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (id, 3, MPI_CHAR, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (mhd, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (start_frame, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (max_frame, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (corr_start_time, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (corr_end_time , 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (average_start_frame, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (histo_time, 1, MPI_INT, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (t_kolmo, 1, MPI_DOUBLE_PRECISION, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (bin_size, 1, MPI_DOUBLE_PRECISION, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (bin_size_w, 1, MPI_DOUBLE_PRECISION, root_process, MPI_COMM_WORLD, ierr)
    call MPI_BCAST (bin_size_angle, 1, MPI_DOUBLE_PRECISION, root_process, MPI_COMM_WORLD, ierr)

end subroutine 

!----------------------------------------------------------------------------------------
! subroutine: initializing arrays after number of particles [lp_number] is known 
!----------------------------------------------------------------------------------------
subroutine initialize 
    use global
    use mpi

    implicit none
    integer :: i
    logical :: inquire, there 
    character*256 :: table_directory 

        !initializing output directories
        write(readdir, '("../results/tmp/output_",A,"/")')id
        !write(readdir, '("../linedeform_mhd/diagoutput/output_",A,"/")')id
        write(outputdir, '("data/sim_",A,"/")')id

    if (proc_id .eq. root_process) then 
        inquire(file=outputdir, exist=there)
        if (there .eqv. .false.)then
            call system('mkdir ' // outputdir)
        end if
    end if

    write(line_evo_file, '(A, "line_evo_",A,".dat")')outputdir, id
    write(mhd_line_evo_file, '(A, "mhd_line_evo_",A,".dat")')outputdir, id
    write(surf_evo_file, '(A, "surf_evo_",A,".dat")')outputdir, id
    write(mean_length_file, '(A,"mean_length_evo_",A,".dat")')outputdir, id
    write(angle_file, '(A,"gamma_",A,".dat")')outputdir, id
    write(line_histo_file, '(A, "line_histo_",A,"_tkolmo_",I2.2,".dat")')outputdir, id, histo_time
    write(surf_histo_file, '(A, "surf_histo_",A,"_tkolmo_",I2.2,".dat")')outputdir, id, histo_time
    write(cg_histo_file, '(A,"cg_histo_",A,"_tkolmo_",I2.2,".dat")')outputdir, id, histo_time
    write(cg_file, '(A,"cg_",A,".dat")')outputdir, id
    write(angle_histo_file, '(A, "angle_histo_",A,"_tkolmo_",I2.2,".dat")')outputdir, id, histo_time
    write(mhd_angle_histo_file, '(A, "mhd_angle_histo_",A,"_tkolmo_",I2.2,".dat")')outputdir, id, histo_time
    write(strain_file, '(A, "strain_",A,".dat")')outputdir, id
    write(corr_file, '(A, "covariance_",A,".dat")')outputdir, id

    if (proc_id .eq. root_process) then 
        table_directory = trim(outputdir) // "tables" 
        inquire(file=table_directory, exist=there)
        if (there .eqv. .false.)then
            call system('mkdir ' // table_directory)
        end if
    end if

    write(results_table, '(A, "tables/tab_results_",A,".txt")')outputdir, id 
    write(zeta_pdf_table, '(A, "tables/tab_zeta_results_pdf_",A,".txt")')outputdir, id 
    write(xi_pdf_table, '(A, "tables/tab_xi_results_pdf_",A,".txt")')outputdir, id 
    write(line_angle_results_table, '(A, "tables/tab_results_angle_line_",A,".txt")')outputdir, id 
    write(surf_angle_results_table, '(A, "tables/tab_results_angle_surf_",A,".txt")')outputdir, id 
    write(eval_results_table, '(A, "tables/tab_results_eval_",A,".txt")')outputdir, id 
    write(vel_file, '(A, "mean_vel_",A,".dat")')outputdir, id 

    !initial values
    le_local = 0
    le_initial_local = 0
    le_length_local = 0
    xi_var = 0
    xi_mean = 0
    zeta_var = 0
    zeta_mean = 0
    eval_mean = 0
    eval_var = 0
    gamma_line = 0
    gamma_line_var = 0
    gamma_surf = 0
    gamma_surf_var = 0
    theta_var = 0
    phi_var = 0
    io_total_time=0
    sort_total_time=0
    comm_total_time=0

    B_local(:,1,:,:) = 0

    do i=1,3!unit matrix
        B_local(:,1,i,i) = 1
    end do

    !call initialize_random 
    call initialize_ordered 
end subroutine initialize 

!----------------------------------------------------------------------------------------
! subroutine: initializes line elements #1 parrallel to x-axis 
!             initializes line elements #2 parrallel to y-axis 
!             initializes line elements #3 parrallel to z-axis 
!----------------------------------------------------------------------------------------
subroutine initialize_ordered 
    use global
    implicit none
    integer :: i

    do i=1, proc_particles
        le_initial_local(i,1,1) = 1._8 !x-coordidate
        le_initial_local(i,2,1) = 0._8 !y-coordidate
        le_initial_local(i,3,1) = 0._8 !z-coordidate

        le_initial_local(i,1,2) = 0._8 !x-coordidate
        le_initial_local(i,2,2) = 1._8 !y-coordidate
        le_initial_local(i,3,2) = 0._8 !z-coordidate

        le_initial_local(i,1,3) = 0._8 !x-coordidate
        le_initial_local(i,2,3) = 0._8 !y-coordidate
        le_initial_local(i,3,3) = 1._8 !z-coordidate
    end do
        
    do i=1,3
        le_length_local(:,1,i) = sqrt(le_initial_local(:,1,i)**2 + le_initial_local(:,2,i)**2 + le_initial_local(:,3,i)**2)
        le_local(:,1,:,i) = le_initial_local(:,:,i) 
    end do
         
    DO i=1,proc_particles
        A_local(i,1,1) = le_local(i,1,2,1)*le_local(i,1,3,2) - le_local(i,1,3,1)*le_local(i,1,2,2) 
        A_local(i,1,2) = le_local(i,1,3,1)*le_local(i,1,1,2) - le_local(i,1,1,1)*le_local(i,1,3,2) 
        A_local(i,1,3) = le_local(i,1,1,1)*le_local(i,1,2,2) - le_local(i,1,2,1)*le_local(i,1,1,2) 
    END DO

    A_length_local(:,1) = sqrt(A_local(:,1,1)**2 + A_local(:,1,2)**2 + A_local(:,1,3)**2)
end subroutine initialize_ordered

!----------------------------------------------------------------------------------------
! subroutine: initializes line elements with random orientation and length unity
!----------------------------------------------------------------------------------------
subroutine initialize_random 
    use global
    implicit none
    integer :: i
    real*8 :: rand

    !initializing 3 random oriented, orthogonal line elements per particle 
    call srand(int(seed,4)) 
    do i=1, proc_particles
        !line element #1
        le_initial_local(i,1,1) = rand() 
        le_initial_local(i,2,1) = rand() 
        le_initial_local(i,3,1) = rand() 

        !line element #2
        le_initial_local(i,1,2) = rand() 
        le_initial_local(i,2,2) = rand() 
        le_initial_local(i,3,2) = -(le_initial_local(i,1,1)*le_initial_local(i,1,2) +&
                            le_initial_local(i,2,1)*le_initial_local(i,2,2))/le_initial_local(i,3,1)
        
        !random orientation for the 3rd line element 
        le3_sign =  sign(1.0,rand()-0.5)

        !line element #3
        le_initial_local(i,1,3) = le3_sign*(le_initial_local(i,2,1)*le_initial_local(i,3,2)&
                                 - le_initial_local(i,3,1)*le_initial_local(i,2,2))
        le_initial_local(i,2,3) = le3_sign*(le_initial_local(i,3,1)*le_initial_local(i,1,2)&
                                 - le_initial_local(i,1,1)*le_initial_local(i,3,2))
        le_initial_local(i,3,3) = le3_sign*(le_initial_local(i,1,1)*le_initial_local(i,2,2)&
                                 - le_initial_local(i,2,1)*le_initial_local(i,1,2))
    end do
        
    do i=1,3
        !computing the norm
        le_length_local(:,1,i) = sqrt(le_initial_local(:,1,i)**2 + le_initial_local(:,2,i)**2 + le_initial_local(:,3,i)**2)

        !normalizing the components
        le_initial_local(:,1,i) = le_initial_local(:,1,i)/le_length_local(:,1,i)
        le_initial_local(:,2,i) = le_initial_local(:,2,i)/le_length_local(:,1,i)
        le_initial_local(:,3,i) = le_initial_local(:,3,i)/le_length_local(:,1,i)

        !computing the unity norm
        le_length_local(:,1,i) = sqrt(le_initial_local(:,1,i)**2 + le_initial_local(:,2,i)**2 + le_initial_local(:,3,i)**2)

        !setting the normalized components
        le_local(:,1,:,i) = le_initial_local(:,:,i) 
    end do
         
    DO i=1,proc_particles
        A_local(I,1,1) = le_local(i,1,2,1)*le_local(i,1,3,2) - le_local(i,1,3,1)*le_local(i,1,2,2) 
        A_local(I,1,2) = le_local(i,1,3,1)*le_local(i,1,1,2) - le_local(i,1,1,1)*le_local(i,1,3,2) 
        A_local(I,1,3) = le_local(i,1,1,1)*le_local(i,1,2,2) - le_local(i,1,2,1)*le_local(i,1,1,2) 
    END DO

    A_length_local(:,1) = sqrt(A_local(:,1,1)**2 + A_local(:,1,2)**2 + A_local(:,1,3)**2)
end subroutine initialize_random

!----------------------------------------------------------------------------------------
! subroutine: allocating the arrays
!----------------------------------------------------------------------------------------
subroutine alloc 
    use global
    implicit none

    allocate(lp_vgr_local(max_lp/num_procs,2,3,3),&
             sort_t_start(num_procs), sort_t_stop(num_procs), sort_total_time(num_procs),&
             comm_t_start(num_procs), comm_t_stop(num_procs), comm_total_time(num_procs),& 
             mag_field_local(max_lp/num_procs,2,3),&
             le_local(max_lp/num_procs,2,3,3), le_initial_local(max_lp/num_procs,3,3), le_length_local(max_lp/num_procs,2,3),&
             B_local(max_lp/num_procs,5,3,3),&
             A_local(max_lp/num_procs,2,3), A_length_local(max_lp/num_procs,2),&
             zeta_mean(max_frame), zeta_var(max_frame),&
             zeta_skew(max_frame), zeta_kurt(max_frame),&
             xi_mean(max_frame), xi_var(max_frame),&
             xi_skew(max_frame), xi_kurt(max_frame),&
             eval_mean(N, max_frame), eval_var(N, max_frame), gamma_line(max_frame,5), &
             gamma_line_var(max_frame,5),gamma_surf(max_frame,5), &
             gamma_surf_var(max_frame,5), theta_var(max_frame),&
             phi_var(max_frame))
end subroutine 
                                                                                                         
!----------------------------------------------------------------------------------------
! subroutine: deallocating the arrays
!----------------------------------------------------------------------------------------
subroutine dealloc 
    use global
    implicit none

    deallocate(lp_vgr_local, mag_field_local,&
              sort_t_start, sort_t_stop, sort_total_time,&
              comm_t_start, comm_t_stop, comm_total_time,& 
              le_local, le_initial_local, le_length_local, &
              A_local, A_length_local, &
              B_local, zeta_mean, zeta_var, zeta_skew, zeta_kurt,&
              xi_mean, xi_var, xi_skew, xi_kurt,&
              eval_mean, eval_var, gamma_line, gamma_line_var,gamma_surf, gamma_surf_var, &
              theta_var, phi_var)

end subroutine
