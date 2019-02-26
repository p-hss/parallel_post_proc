!----------------------------------------------------------------------------------------
! subroutine: intergrating B in time using RK4. Resulting in B(t+dt) 
!----------------------------------------------------------------------------------------
subroutine b_int()
    use global
    implicit none
     
    integer :: i,j!loop indeces
    integer :: lpi

    do lpi = 1, proc_particles!loop over all particles per processor

        !Runge-Kutta time integration scheme
        b_local(lpi,2,:,:) = b_local(lpi,1,:,:) &
                        + 0.5_8*dt*matmul(lp_vgr_local(lpi,1,:,:), b_local(lpi,1,:,:))   

        b_local(lpi,3,:,:) = b_local(lpi,1,:,:) &
                        + 0.25_8*dt*matmul(lp_vgr_local(lpi,1,:,:) + lp_vgr_local(lpi,2,:,:),b_local(lpi,2,:,:))   

        b_local(lpi,4,:,:) = b_local(lpi,1,:,:) &
                        + 0.5_8*dt*matmul(lp_vgr_local(lpi,1,:,:)  + lp_vgr_local(lpi,2,:,:),b_local(lpi,3,:,:))   

        b_local(lpi,5,:,:) = b_local(lpi,1,:,:) &
                        + 1.0_8/6.0_8*dt*(matmul(lp_vgr_local(lpi,1,:,:), b_local(lpi,1,:,:))&
                        + matmul(lp_vgr_local(lpi,1,:,:) + lp_vgr_local(lpi,2,:,:), b_local(lpi,2,:,:) + b_local(lpi,3,:,:))&
                        + matmul(lp_vgr_local(lpi,2,:,:), b_local(lpi,4,:,:)))  

    end do

end subroutine b_int

!----------------------------------------------------------------------------------------
! subroutine: evolving the material line element length using the updated matrix B
!----------------------------------------------------------------------------------------
subroutine line_evo()
    use global
    use mpi
    implicit none
    integer :: i,lpi, frame
    real*8  :: vector_length 

    do lpi = 1, proc_particles!loop over all particles per processor
        
        !updating line elements 
        do i=1,3
            !evolving the line element in time
            le_local(lpi,2,:,i) = matmul(b_local(lpi,5,:,:), le_initial_local(lpi,:,i)) !Carefull b-index =  5?

            !computing the lenght of a line element
            le_length_local(lpi,2,i) = vector_length(le_local(lpi,2,:,i))
        end do 
        
        !a =  l1 x l2
        a_local(lpi,2,1) =  le_local(lpi,2,2,1)*le_local(lpi,2,3,2) - le_local(lpi,2,3,1)*le_local(lpi,2,2,2) 
        a_local(lpi,2,2) =  le_local(lpi,2,3,1)*le_local(lpi,2,1,2) - le_local(lpi,2,1,1)*le_local(lpi,2,3,2) 
        a_local(lpi,2,3) =  le_local(lpi,2,1,1)*le_local(lpi,2,2,2) - le_local(lpi,2,2,1)*le_local(lpi,2,1,2) 

        a_length_local(lpi,2) = sqrt(a_local(lpi,2,1)**2 + a_local(lpi,2,2)**2 + a_local(lpi,2,3)**2)

    end do

end subroutine line_evo

!----------------------------------------------------------------------------------------
! subroutine: copmuting growth rates: Zeta for lines 
!                                   : Xi for surfaces
!----------------------------------------------------------------------------------------
subroutine growth_rates(frame)
    use global
    use mpi
    implicit none

    integer                     :: lpi, i, frame, source, destination, N_local, N_global
    real*8                      :: vector_length, theta_mean, phi_mean
    real*8, dimension(n,n)      :: mat_buffer !storing the VG for eingenv* computation
    real*8, dimension(n,n)      :: strain_rate, rotation_rate !tensors
    real*8, dimension(n)        :: line_element_dir !storing the VG for eingenv* computation
    real*8, dimension(proc_particles)   :: zeta, xi !storing the VG for eingenv* computation
    !real*8                      :: mean_a, mean_b, var_a, var_b
    character(1)                :: creturn
    real*8                      :: zeta_sum_local, zeta_isum
    real*8                      :: delta, M1, M2, M3, M4
    real*8                      :: mean_global, var_global, skew_global, kurt_global
    real*8                      :: sum_mean, sum_var, sum_skew, sum_kurt
    real*8                      :: sum_mean_global, sum_var_global, sum_skew_global, sum_kurt_global
    character*10                :: mode
    integer, dimension(MPI_STATUS_SIZE) :: status
    real*8, dimension(:), allocatable       :: zeta_global

    if (proc_id .eq. root_process) then 
        allocate(zeta_global(max_lp))
        zeta_global=0
    end if

    !===================compute the line stretching rate zeta===================
    do lpi = 1, proc_particles!loop over all particles per processor

        !============method #1===================
        !mat_buffer(:,:) = lp_vgr_local(lpi,2,:,:) 
        !call mat_decomp(mat_buffer, strain_rate, rotation_rate)
        !line_element_dir(:) = matmul(strain_rate(:,:), le_local(lpi,2,:,1)/le_length_local(lpi,2,1))
        !zeta(lpi) = dot_product(line_element_dir(:), le_local(lpi,2,:,1)/le_length_local(lpi,2,1))
        !============method #2===================
        zeta(lpi) = (log(le_length_local(lpi,2,1)) - log(le_length_local(lpi,1,1)))/dt 
        !call mean(lpi, zeta(lpi), sum_mean)
    end do

    !===================Gather all zetas for running stats===================
    call MPI_GATHER(zeta(:), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                    zeta_global(:), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                    root_process, MPI_COMM_WORLD, ierr)

    if (proc_id .eq. root_process) then 
        delta=0; M1=0; M2=0; M3=0; M4=0

        do lpi = 1, max_lp
            call running_stats(lpi, zeta_global(lpi), delta, M1, M2, M3, M4)
        end do

        zeta_mean(frame)= M1
        zeta_var(frame) = M2/max_lp
        zeta_skew(frame)= sqrt(real(max_lp))*M3/M2**3/2     
        zeta_kurt(frame)= max_lp*M4/M2**2 
    end if

    !===================Compute the rest of the moments in parallel ===================
    !call MPI_REDUCE(sum_mean, sum_mean_global, 1, MPI_DOUBLE_PRECISION,&
    !                MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    !mean_global=sum_mean_global/(real(max_lp))
    !zeta_mean(frame)=mean_global

    !call MPI_BCAST(mean_global, 1, MPI_DOUBLE_PRECISION,&
    !            root_process, MPI_COMM_WORLD, ierr)

    !do lpi = 1, proc_particles
    !    call var( lpi, zeta(lpi), sum_var,  mean_global)
    !    !call skew(lpi, zeta(lpi), sum_skew, mean_global)
    !    call kurt(lpi, zeta(lpi), sum_kurt, mean_global)
    !end do

    !call MPI_REDUCE(sum_var, sum_var_global, 1, MPI_DOUBLE_PRECISION,&
    !                MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    !!call MPI_REDUCE(sum_skew, sum_skew_global, 1, MPI_DOUBLE_PRECISION,&
    !!                MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    !call MPI_REDUCE(sum_kurt, sum_kurt_global, 1, MPI_DOUBLE_PRECISION,&
    !                MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    !var_global=sum_var_global/(real(max_lp))
    !zeta_var(frame)=var_global
    !!zeta_skew(frame)=sum_skew_global/(real(max_lp)*var_global**(3./2.))
    !zeta_kurt(frame)=sum_kurt_global/(real(max_lp)*var_global**(4./2.))

    if(frame == histo_frame)then
        mode="zeta"
        call histogram1D(zeta, mode)
    end if

    !===================compute the surface stretching rate xi===================
    do lpi = 1, proc_particles!loop over all particles per processor
        xi(lpi) = (log(a_length_local(lpi,2)) - log(a_length_local(lpi,1)))/dt 
        call mean(lpi, xi(lpi), sum_mean)
     end do

    call MPI_REDUCE(sum_mean, sum_mean_global, 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    mean_global=sum_mean_global/(real(max_lp))
    xi_mean(frame)=mean_global

    call MPI_BCAST(mean_global, 1, MPI_DOUBLE_PRECISION,&
                root_process, MPI_COMM_WORLD, ierr)

    do lpi = 1, proc_particles
        call var( lpi, xi(lpi), sum_var, mean_global)
        call skew(lpi, xi(lpi), sum_skew, mean_global)
        call kurt(lpi, xi(lpi), sum_kurt, mean_global)
    end do

    call MPI_REDUCE(sum_var, sum_var_global, 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sum_skew, sum_skew_global, 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sum_kurt, sum_kurt_global, 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    var_global=sum_var_global/(real(max_lp-1))
    xi_var(frame)=var_global
    xi_skew(frame)=sum_skew_global/(real(max_lp)*var_global**(3./2.))
    xi_kurt(frame)=sum_kurt_global/(real(max_lp)*var_global**(4./2.))

    if(frame == histo_frame)then
        mode="xi"
        call histogram1D(xi, mode)
    end if

    if (proc_id .eq. root_process) then 
        theta_mean = -2._8*zeta_mean(frame) + xi_mean(frame)
        phi_mean = -(zeta_mean(frame) + xi_mean(frame))
        !================================Output to file==========================================
        if(frame == start_frame +1)then 
            open(unit=102, file=line_evo_file, action='write', form='formatted')
            write(102,*) 0, 0, 0, 0
        else
            open(unit=102, file=line_evo_file, action='write', form='formatted', position='append')
        end if 
            write(102,*) (frame-start_frame)*dt, zeta_mean(frame), zeta_var(frame),&
                        theta_mean, zeta_skew(frame), zeta_kurt(frame)
            close(102)

        if(frame == start_frame +1)then 
            open(unit=103, file=surf_evo_file, status='unknown')
            write(103,*) 0, 0, 0, 0
        else
            open(unit=103, file=surf_evo_file, action='write', form='formatted', position='append')
        end if 
            write(103,*) (frame-start_frame)*dt, xi_mean(frame), xi_var(frame),&
                        phi_mean, xi_skew(frame), xi_kurt(frame)
            close(103)
        !================================Output to Terminal==========================================
        creturn = achar(13)
        30 format(A1,I3.3,A5,4(A15, ES12.3))
        write(*,30,advance='no') creturn,& 
            int((frame-start_frame)/real((max_frame-start_frame))*100)," % ",&
           "mean[zeta]:", zeta_mean(frame), "var[zeta]:", zeta_var(frame),&
           "skew[zeta]:", zeta_skew(frame), "kurt[zeta]:", zeta_kurt(frame)

        if(iframe==max_frame) print*,""
    end if

    if (proc_id .eq. root_process) then 
        deallocate(zeta_global)
    end if
end subroutine growth_rates

subroutine mean(i, x, sum_mean)
    use mpi, 
    use global
    implicit none
    integer, intent(in)     :: i
    real*8, intent(inout)   :: sum_mean
    real*8, intent(in)      :: x

    if(i==1) sum_mean=0

    sum_mean=sum_mean + x 

end subroutine mean

subroutine var(i, x, sum_var, mean)
    use mpi,
    use global
    implicit none
    integer, intent(in)     :: i
    real*8, intent(inout)   :: sum_var
    real*8, intent(in)      :: x, mean
    real*8                  :: tmp
    real*8, save            :: isum

    if(i==1) sum_var=0

    sum_var=sum_var + (x - mean)**2

end subroutine var

subroutine skew(i, x, sum_skew, mean)
    use mpi
    use global
    implicit none

    integer, intent(in)     :: i
    real*8, intent(inout)   :: sum_skew
    real*8, intent(in)      :: x, mean
    real*8                  :: tmp
    real*8, save            :: isum

    if(i==1) sum_skew=0

    sum_skew=sum_skew + (x - mean)**3

end subroutine skew

subroutine kurt(i, x, sum_kurt, mean)
    use mpi
    use global
    implicit none

    integer, intent(in)     :: i
    real*8, intent(inout)   :: sum_kurt
    real*8, intent(in)      :: x, mean
    real*8                  :: tmp
    real*8, save            :: isum

    if(i==1) sum_kurt=0

    sum_kurt=sum_kurt + (x - mean)**4

end subroutine kurt

!----------------------------------------------------------------------------------------
! subroutine: computes statistical moments I-IV in one pass 
!----------------------------------------------------------------------------------------
subroutine running_stats(i, xi, delta, M1, M2, M3, M4)
    implicit none

    integer, intent(in)     :: i 
    real*8                  :: xi
    real*8                  :: n, delta, M1, M2, M3, M4

    n=real(i,8)

    delta = xi - M1

    M1 = M1*(n-1)/n + xi/n

    M4 = M4 + delta**4*(n-1)*(n**2-3._8*n+3)/n**3 + 6._8*delta**2*M2/n**2 - 4._8*delta*M3/n

    M3 = M3 + delta**3*((n-1)*(n-2)/n**2)-(3._8*delta*M2)/n

    M2 = M2 + delta**2*(n-1)/n

end subroutine running_stats


!function parallel_mean(M1_a, N_a, M1_b, N_b) result(M1)
!    implicit none
!    real*8                  :: M1_a, M1_b, M1
!    integer                 :: N_a, N_b
!
!    M1 = (M1_a*N_a + M1_b*N_b)/(real(N_a + N_b))
!
!end function parallel_mean
!
!function parallel_variance(M1_a, N_a, M2_a, M1_b, N_b, M2_b) result(M2)
!    implicit none
!    real*8                  :: delta, M1_a, M1_b, M2_a, M2_b, M_a, M_b, M2
!    integer                 :: N_b, N_a
!
!    delta = M1_b - M1_a
!
!    M_a = M2_a*(N_a - 1)
!
!    M_b = M2_b*(N_b - 1)
!
!    M2 = M_a + M_b + delta**2*N_a*N_b/(N_a + N_b)
!
!    M2 = M2/(N_a + N_b - 1)
!end function parallel_variance


!----------------------------------------------------------------------------------------
! subroutine: copmuting growth rates: Zeta for lines 
!                                   : Xi for surfaces
!----------------------------------------------------------------------------------------
!subroutine growth_rates(frame)
!    use global
!    use mpi
!    implicit none
!
!    integer                     :: lpi, i, frame, source, destination, N_a, N_b
!    real*8                      :: vector_length, theta_mean, phi_mean
!    real*8, dimension(n,n)      :: mat_buffer !storing the VG for eingenv* computation
!    real*8, dimension(n,n)      :: strain_rate, rotation_rate !tensors
!    real*8, dimension(n)        :: line_element_dir !storing the VG for eingenv* computation
!    real*8, dimension(:), allocatable    :: lp_vgr_global
!    real*8                      :: delta, M1, M2, M3, M4
!    real*8                      :: mean_a, mean_b, var_a, var_b
!    real*8                      :: zeta, xi
!    integer*4                   :: seed_1
!    character(1)                :: creturn
!    real*8                      :: zeta_mean_global, parallel_variance, parallel_mean
!    integer, dimension(MPI_STATUS_SIZE) :: status
!
!    allocate(lp_vgr_global(max_lp,3,3))
!
!    delta=0; M1=0; M2=0; M3=0; M4=0
!    mean_a=0; mean_b=0; var_a=0; var_b=0
!
!
!    do j=1,3
!        do i=1,3
!            call MPI_GATHER(lp_vgr_local_unsort(:,2,i,j), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
!                             lp_vgr_global(:,i,j), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
!                             root_process, MPI_COMM_WORLD, ierr)
!        end do
!    end do
!
!    if (proc_id == root_process) then 
!        !compute the line stretching rate zeta
!        do lpi = 1, max_lp!loop over all particles per processor
!
!            !method #1
!            mat_buffer(:,:) = lp_vgr_global(lpi,2,:,:) 
!            call mat_decomp(mat_buffer, strain_rate, rotation_rate)
!
!            line_element_dir(:) = matmul(strain_rate(:,:), le_local(lpi,2,:,1)/le_length_local(lpi,2,1))
!
!            zeta = dot_product(line_element_dir(:), le_local(lpi,2,:,1)/le_length_local(lpi,2,1))
!
!            !method #2
!            !zeta = (log(le_length_local(lpi,2,1)) - log(le_length_local(lpi,1,1)))/dt 
!
!            call running_stats(lpi, zeta, delta, M1, M2, M3, M4)
!
!        end do
!    end if
!
!!    if (proc_id == root_process) then 
!!        !print*,""
!!        mean_a = M1
!!        var_a = M2/proc_particles
!!    end if
!!
!!    !merging statistical moments
!!    do source = 1, num_procs-1
!!
!!        destination = 0
!!
!!        if (proc_id == source) then 
!!
!!            mean_b=M1
!!            var_b=M2/proc_particles
!!
!!            call MPI_SEND( mean_b, 1, MPI_DOUBLE_PRECISION, destination, 1, MPI_COMM_WORLD, ierr)
!!            call MPI_SEND(  var_b, 1, MPI_DOUBLE_PRECISION, destination, 1, MPI_COMM_WORLD, ierr)
!!
!!        end if
!!
!!        if (proc_id == root_process) then 
!!
!!            call MPI_RECV( mean_b, 1, MPI_DOUBLE_PRECISION, source, 1, MPI_COMM_WORLD, status, ierr)
!!            call MPI_RECV(  var_b, 1, MPI_DOUBLE_PRECISION, source, 1, MPI_COMM_WORLD, status, ierr)
!!
!!            N_a = source*max_lp/num_procs
!!            N_b = max_lp/num_procs
!!
!!            mean_a = (mean_a*N_a + mean_b*N_b)/(N_a + N_b)
!!            var_a = parallel_variance(mean_a, N_a, var_a, mean_b, N_b, var_b)
!!        end if
!!
!!    end do
!
!    zeta_mean(frame)= M1
!    zeta_var(frame) = M2/proc_particles
!    zeta_skew(frame)= sqrt(real(proc_particles))*M3/M2**3/2     
!    zeta_kurt(frame)= proc_particles*M4/M2**2 
!
!    delta=0; M1=0; M2=0; M3=0; M4=0
!
!    !compute the surface stretching rate xi
!    do lpi = 1, proc_particles!loop over all particles per processor
!        
!        !xi(lpi) = (log(a_length(lpi,2)) - log(a_length(lpi,1)))/dt 
!        xi = (log(a_length_local(lpi,2)) - log(a_length_local(lpi,1)))/dt 
!
!        !call running_stats(lpi, xi(lpi) , delta, M1, M2, M3, M4)
!        call running_stats(lpi, xi , delta, M1, M2, M3, M4)
!
!    end do
!
!    xi_mean(frame)= M1
!    xi_var(frame) = M2/proc_particles
!    xi_skew(frame)= sqrt(real(proc_particles))*M3/M2**3/2     
!    xi_kurt(frame)= proc_particles*M4/M2**2 
!
!    theta_mean = -2._8*zeta_mean(frame) + xi_mean(frame)
!    phi_mean = -(zeta_mean(frame) + xi_mean(frame))
!
!    if (proc_id .eq. root_process) then 
!        !================================Output to file==========================================
!        if(frame == start_frame +1)then 
!            open(unit=102, file=line_evo_file, action='write', form='formatted')
!            write(102,*) 0, 0, 0, 0
!        else
!            open(unit=102, file=line_evo_file, action='write', form='formatted', position='append')
!        end if 
!
!            write(102,*) (frame-start_frame)*dt, zeta_mean(frame), zeta_var(frame),&
!                        theta_mean, zeta_skew(frame), zeta_kurt(frame)
!            close(102)
!
!        if(frame == start_frame +1)then 
!            open(unit=103, file=surf_evo_file, status='unknown')
!            write(103,*) 0, 0, 0, 0
!        else
!            open(unit=103, file=surf_evo_file, action='write', form='formatted', position='append')
!        end if 
!
!            write(103,*) (frame-start_frame)*dt, xi_mean(frame), xi_var(frame),&
!                        phi_mean, xi_skew(frame), xi_kurt(frame)
!
!            close(103)
!
!        !================================Output to Terminal==========================================
!        creturn = achar(13)
!        30 format(A1,I3.3,A5,4(A15, ES12.3))
!
!        write(*,30,advance='no') creturn,& 
!            int((frame-start_frame)/real((max_frame-start_frame))*100)," % ",&
!            "mean[zeta]:", zeta_mean(frame), "var[zeta]:", zeta_var(frame),&
!            "skew[zeta]:", zeta_skew(frame), "kurt[zeta]:", zeta_kurt(frame)
!
!        if(iframe==max_frame) print*,""
!    end if
!
!    deallocate(lp_vgr_global)
!
!end subroutine growth_rates
