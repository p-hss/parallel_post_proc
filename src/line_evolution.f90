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
    character*10                :: mode
    integer, dimension(MPI_STATUS_SIZE) :: status
    real*8, dimension(:), allocatable       :: zeta_global, xi_global

    if (proc_id .eq. root_process) then 
        allocate(zeta_global(max_lp), xi_global(max_lp))
        zeta_global=0; xi_global=0
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

    if(frame == histo_frame)then
        mode="zeta"
        call histogram1D(zeta, mode)
    end if

    !===============compute the surface stretching rate xi===================
    do lpi = 1, proc_particles!loop over all particles per processor
        xi(lpi) = (log(a_length_local(lpi,2)) - log(a_length_local(lpi,1)))/dt 
     end do

    !===================Gather all xis for running stats=====================
    call MPI_GATHER(xi(:), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                    xi_global(:), max_lp/num_procs, MPI_DOUBLE_PRECISION,&
                    root_process, MPI_COMM_WORLD, ierr)

    if (proc_id .eq. root_process) then 
        delta=0; M1=0; M2=0; M3=0; M4=0

        do lpi = 1, max_lp
            call running_stats(lpi, xi_global(lpi), delta, M1, M2, M3, M4)
        end do

        xi_mean(frame)= M1
        xi_var(frame) = M2/max_lp
        xi_skew(frame)= sqrt(real(max_lp))*M3/M2**3/2     
        xi_kurt(frame)= max_lp*M4/M2**2 
    end if

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
        deallocate(zeta_global, xi_global)
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

    M1 = M1 + delta/n

    M4 = M4 + delta**4*(n-1)*(n**2-3._8*n+3)/n**3 + 6._8*delta**2*M2/n**2 - 4._8*delta*M3/n

    M3 = M3 + delta**3*((n-1)*(n-2)/n**2)-(3._8*delta*M2)/n

    M2 = M2 + delta**2*(n-1)/n

end subroutine running_stats
