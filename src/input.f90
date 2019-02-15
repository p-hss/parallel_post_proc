!----------------------------------------------------------------------------------------
! subroutine: read data from binary MHDT files
!----------------------------------------------------------------------------------------
subroutine input 
    use mpi
    use global
    implicit none

    integer :: index_1, index_2, index_3
    integer :: io,i,j,lp,k
    integer :: lp_loc, fcount
    integer :: io_start, io_stop 
    real*8  :: time
    real*8, dimension(:,:,:), allocatable   :: lp_vgr_global_unsort, lp_vgr_local_unsort 
    real*8, dimension(:,:), allocatable     :: lp_pos_unsort, lp_vel_unsort, mag_field_global_unsort, mag_field_local_unsort
    integer*4, dimension(:), allocatable    :: lp_ID_list, lp_ID_list_local

    allocate(lp_ID_list(max_lp), lp_ID_list_local(max_lp/num_procs), lp_pos_unsort(max_lp,3),&
            lp_vel_unsort(max_lp,3), lp_vgr_global_unsort(max_lp,3,3),&
            lp_vgr_local_unsort(max_lp/num_procs,3,3), mag_field_global_unsort(max_lp,3),&
            mag_field_local_unsort(max_lp/num_procs,3))

    if (proc_id .eq. root_process) call cpu_time(io_t_start)

    if (proc_id .eq. root_process) then 

        io_start=1
        io_stop=0
        index_1=0
        index_2=iframe

        do index_3=0, nioproc-1
            
            write(file_out, '(A,"Position.",I3.3,".",I4.4,"_",I4.4)')readdir, index_1, index_2, index_3

            open(unit=201, file=file_out, form="unformatted", status='old', action='read')
                
                !counting number of records in input file 
                do
                    read(201,iostat=io) 
                    if(io==-1) exit
                    fcount= fcount + 1
                end do

                fcount=0

                rewind(201)

                read(201,iostat=io) time
                read(201,iostat=io) lp_number
                
                !counting number of particles per input file
                io_stop=io_stop + lp_number

                if(iframe==start_frame)then
                    lp=1
                    dt_init=time
                end if

                if(iframe==start_frame+1)then
                    dt_fin=time
                    dt=(dt_fin - dt_init)/t_kolmo
                end if

                !determine frames for line element analysis
                histo_frame=int(histo_time/dt+start_frame)
                corr_start=int(corr_start_time/dt+start_frame)
                corr_end=int(corr_end_time/dt+start_frame)
                average_start_frame=int(average_start_time/dt+start_frame)

                read(201,iostat=io) lp_ID_list(io_start:io_stop)
                !reading positions
                read(201,iostat=io) lp_pos_unsort(io_start:io_stop,1)
                read(201,iostat=io) lp_pos_unsort(io_start:io_stop,2)
                read(201,iostat=io) lp_pos_unsort(io_start:io_stop,3)
                !reading velocities 
                read(201,iostat=io) lp_vel_unsort(io_start:io_stop,1)
                read(201,iostat=io) lp_vel_unsort(io_start:io_stop,2)
                read(201,iostat=io) lp_vel_unsort(io_start:io_stop,3)
                !reading velocity gradient: 
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,1,1)!dx vx
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,1,2)!dy vx
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,1,3)!dz vx
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,2,1)!dx vy
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,2,2)!dy vy
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,2,3)!dz vy
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,3,1)!dx vz
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,3,2)!dy vz
                read(201,iostat=io) lp_vgr_global_unsort(io_start:io_stop,3,3)!dz vz

                if(mhd == 1)then
                    read(201,iostat=io) mag_field_global_unsort(io_start:io_stop,1)!m x
                    read(201,iostat=io) mag_field_global_unsort(io_start:io_stop,2)!m y
                    read(201,iostat=io) mag_field_global_unsort(io_start:io_stop,3)!m z
                end if

            close(unit=201)
            
            if(io > 0)then
                print*, "reading: something went wrong" 
            else if(io < 0)then
                print*, "reading: end of file reached" 
            end if

            io_start=io_stop + 1
        end do

    end if

    if (proc_id .eq. root_process)then
        call cpu_time(io_t_stop)
        io_total_time = io_total_time + io_t_stop - io_t_start
    end if

    if (proc_id .eq. root_process) call cpu_time(comm_t_start)

    call MPI_SCATTER(lp_ID_list(:), max_lp/num_procs, MPI_INTEGER,&
                     lp_ID_list_local(:), max_lp/num_procs, MPI_INTEGER,&
                     root_process, MPI_COMM_WORLD, ierr)

    do j=1,3
        do i=1,3
            call MPI_SCATTER(lp_vgr_global_unsort(:,i,j), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                             lp_vgr_local_unsort(:,i,j), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                             root_process, MPI_COMM_WORLD, ierr)
        end do
    end do

    do i=1,3
        call MPI_SCATTER(mag_field_global_unsort(:,i), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                         mag_field_local_unsort(:,i), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                         root_process, MPI_COMM_WORLD, ierr)
    end do

    if (proc_id .eq. root_process)then
        call cpu_time(comm_t_stop)
        comm_total_time = comm_total_time + comm_t_stop - comm_t_start
    end if

    if (proc_id .eq. root_process) call cpu_time(sort_t_start)

    do lp=1,proc_particles

        !map unsorted index lp_loc to sorted block index lp 
        lp_loc=minloc(abs(lp_ID_list_local - lp), 1)

        lp_vgr_local(lp,2,1,1)=lp_vgr_local_unsort(lp_loc,1,1)
        lp_vgr_local(lp,2,1,2)=lp_vgr_local_unsort(lp_loc,1,2)
        lp_vgr_local(lp,2,1,3)=lp_vgr_local_unsort(lp_loc,1,3)
        lp_vgr_local(lp,2,2,1)=lp_vgr_local_unsort(lp_loc,2,1)
        lp_vgr_local(lp,2,2,2)=lp_vgr_local_unsort(lp_loc,2,2)
        lp_vgr_local(lp,2,2,3)=lp_vgr_local_unsort(lp_loc,2,3)
        lp_vgr_local(lp,2,3,1)=lp_vgr_local_unsort(lp_loc,3,1)
        lp_vgr_local(lp,2,3,2)=lp_vgr_local_unsort(lp_loc,3,2)    
        lp_vgr_local(lp,2,3,3)=lp_vgr_local_unsort(lp_loc,3,3)    

        lp_vgr_local(lp,2,:,:)=lp_vgr_local(lp,2,:,:)*t_kolmo

        if(iframe==start_frame) then
            lp_vgr_local(:,1,:,:) = lp_vgr_local(:,2,:,:)
        end if

        if(mhd == 1)then
            mag_field_local(lp,2,1)=mag_field_local_unsort(lp_loc,1)
            mag_field_local(lp,2,2)=mag_field_local_unsort(lp_loc,2)
            mag_field_local(lp,2,3)=mag_field_local_unsort(lp_loc,3)
        end if
    end do

    if (proc_id .eq. root_process)then
         call cpu_time(sort_t_stop)
        sort_total_time = sort_total_time + sort_t_stop - sort_t_start
    end if
        
    if(iframe==start_frame+1)then
        call MPI_BCAST(histo_frame, 1, MPI_INTEGER,&
                    root_process, MPI_COMM_WORLD, ierr)
    end if

    if(iframe==start_frame+1)then
        call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION,&
                    root_process, MPI_COMM_WORLD, ierr)
    end if

    deallocate(lp_ID_list, lp_ID_list_local,  lp_pos_unsort, lp_vel_unsort, lp_vgr_global_unsort,&
                lp_vgr_local_unsort, mag_field_global_unsort, mag_field_local_unsort)

end subroutine input
