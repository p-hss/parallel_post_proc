
!----------------------------------------------------------------------------------------
! subroutine: the cauchy green tensor and its eigenvalues for fluid element deformation 
!             analysis.
!----------------------------------------------------------------------------------------
subroutine cauchy_green(frame)
    use mpi
    use global
    implicit none
    integer :: lpi, frame, i,j
    real*8, dimension(n,n)  :: b_t
    real*8, dimension(n,n)  :: cg_tensor, cg_evec
    real*8, dimension(n)    :: cg_eval, cg_eval_mean, cg_eval_sorted
    real*8, dimension(n)    :: sum_mean, sum_mean_global, mean_global
    real*8, dimension(proc_particles)   :: eval1, eval2 !storing the VG for eingenv* computation
    character*10            :: mode

    do lpi = 1, proc_particles!loop over all particles per processor

        !compute transpose
        do i = 1, n 
            do j = 1, n 
                b_t(i,j) = b_local(lpi,5,j,i)
            end do
        end do
        
        !compute cauchy green tensor
        cg_tensor(:,:) = matmul(b_local(lpi,5,:,:), b_t(:,:))

        call eigen_val_vec(cg_tensor(:,:), cg_eval(:), cg_evec(:,:))

        !sorting
        cg_eval_sorted=0
        cg_eval_sorted(1)=maxval(cg_eval)
        cg_eval_sorted(3)=minval(cg_eval)

        if(cg_eval(1) .lt. cg_eval_sorted(1) .and.&
           cg_eval(1) .gt. cg_eval_sorted(3))then
                    cg_eval_sorted(2) = cg_eval(1)
        else if(cg_eval(2) .lt. cg_eval_sorted(1) .and.&
                cg_eval(2) .gt. cg_eval_sorted(3))then
                    cg_eval_sorted(2) = cg_eval(2)
        else if(cg_eval(3) .lt. cg_eval_sorted(1) .and.&
                cg_eval(3) .gt. cg_eval_sorted(3))then
                    cg_eval_sorted(2) = cg_eval(3)
        end if

        eval1(lpi)=log(cg_eval_sorted(1))
        eval2(lpi)=log(cg_eval_sorted(2))
    
        call mean(lpi, eval1(lpi), sum_mean(1))
        call mean(lpi, eval2(lpi), sum_mean(2))
        call mean(lpi, -(eval1(lpi) + eval2(lpi)), sum_mean(3))
        
    end do

    if(frame == histo_frame)then
        mode="cauchy"
        call histogram2D(eval1, eval2, mode)
    end if

    call MPI_REDUCE(sum_mean(1), sum_mean_global(1), 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sum_mean(2), sum_mean_global(2), 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sum_mean(3), sum_mean_global(3), 1, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    mean_global(1)=sum_mean_global(1)/(real(max_lp))
    mean_global(2)=sum_mean_global(2)/(real(max_lp))
    mean_global(3)=sum_mean_global(3)/(real(max_lp))

    if (proc_id .eq. root_process) then 
        if(frame == start_frame+1)then 
            open(unit=104, file=cg_file, action='write', form='formatted')
        else
            open(unit=104, file=cg_file, action='write', form='formatted', position='append')
        end if 
                write(104,*) (frame-start_frame)*dt, mean_global(1),&
                                mean_global(2), mean_global(3)
            if(verbose == 8)then
                write(*,*) (frame-start_frame)*dt, mean_global(1),&
                                mean_global(2), mean_global(3)
            end if
        close(104)
    end if

end subroutine cauchy_green 

!----------------------------------------------------------------------------------------
! subroutine: comuptes histograms   
!----------------------------------------------------------------------------------------
subroutine histogram1D(variable, mode)
    use global
    use mpi
    implicit none
    integer :: lpi, frame, i,j
    integer :: max_val, count_excess, count_excess_global
    integer :: bin, number_bin
    integer, dimension(2) :: max_val_w 
    real*8, dimension(proc_particles), intent(in) :: variable
    real*8, dimension(:), allocatable :: histo, histo_global 
    character*10 :: mode

    max_val = 10
    number_bin = int(2*max_val/bin_size) 

    allocate(histo(number_bin), histo_global(number_bin))

    histo = 0
    histo_global = 0
    count_excess = 0

    do lpi=1, proc_particles
        if(variable(lpi) < max_val)then
            bin = ceiling((variable(lpi) + max_val)/bin_size)
            histo(bin) =  histo(bin) + 1
        else
            count_excess = count_excess + 1
        end if
    end do

    call MPI_REDUCE(histo(:), histo_global(:), number_bin, MPI_DOUBLE_PRECISION,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    call MPI_REDUCE(count_excess, count_excess_global, 1, MPI_INTEGER,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    histo_global(:) = histo_global(:)/((max_lp-real(count_excess_global))*bin_size)

    if (proc_id .eq. root_process) then 
        if(Trim(mode) == "zeta")then
            open(unit=102, file=line_histo_file, status='unknown')
        elseif(Trim(mode) == "xi")then
            open(unit=102, file=surf_histo_file, status='unknown')
        end if
        do i=1, number_bin
            write(102,*) i*bin_size - max_val, histo_global(i)
        end do
        close(102)
    end if

    deallocate(histo, histo_global)


end subroutine histogram1D

!----------------------------------------------------------------------------------------
! subroutine: comuptes histograms   
!----------------------------------------------------------------------------------------
subroutine histogram2D(variable1, variable2, mode)
    use global
    use mpi
    implicit none
    integer :: lpi,i,j
    integer :: count_excess, count_excess_global
    integer :: bin, number_bin
    integer, dimension(2) :: bin_w, number_bin_w 
    integer, dimension(2) :: max_val_w 
    real*8, dimension(proc_particles) :: variable1, variable2
    real*8, dimension(:,:), allocatable :: histo, histo_global
    character*10 :: mode

    max_val_w(1) = 20
    max_val_w(2) = 20
    count_excess = 0

    number_bin_w(1) = int(2*max_val_w(1)/bin_size_w) 
    number_bin_w(2) = int(2*max_val_w(2)/bin_size_w) 

    allocate(histo(number_bin_w(1),number_bin_w(2)), histo_global(number_bin_w(1),number_bin_w(2)))

    histo=0; histo_global=0; bin_w=0

    do lpi=1, proc_particles
        if(variable1(lpi) < max_val_w(1) .and. &
           variable1(lpi) > -max_val_w(1) .and. &
           variable2(lpi) < max_val_w(2) .and. &
           variable2(lpi) > -max_val_w(2) )then
            bin_w(1) = ceiling((variable1(lpi) + max_val_w(1))/bin_size_w)
            bin_w(2) = ceiling((variable2(lpi) + max_val_w(2))/bin_size_w)
            histo(bin_w(1), bin_w(2)) = histo(bin_w(1), bin_w(2))  + 1
        else
            count_excess = count_excess + 1
        end if
    end do

    call MPI_REDUCE(count_excess, count_excess_global, 1, MPI_INTEGER,&
                    MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

    do i=1,number_bin_w(2)
        call MPI_REDUCE(histo(:,i), histo_global(:,i), number_bin_w(1), MPI_DOUBLE_PRECISION,&
                        MPI_SUM, root_process, MPI_COMM_WORLD, ierr)
    end do

    histo_global(:,:)=histo_global(:,:)/((max_lp-real(count_excess_global))*bin_size_w)

    if (proc_id .eq. root_process) then 
        if(Trim(mode) == "cauchy")then
            open(unit=102, file=cg_histo_file, status='unknown')
        end if
        do j=1, number_bin_w(2)
            do i=1, number_bin_w(1)
                if(histo_global(i,j) > 0.00000000001)then
                    write(102,*) i*bin_size_w - max_val_w(1), j*bin_size_w - max_val_w(2),  histo_global(i,j)
                end if
            end do
        end do
        close(102)
    end if
    deallocate(histo, histo_global)

end subroutine histogram2D

!----------------------------------------------------------------------------------------
! subroutine:   computes the maximum positive and negative strain rates
!               and determines the angle with respect to the line elements 
!----------------------------------------------------------------------------------------
subroutine strain_stats(frame)
    use global
    implicit none
    integer                     :: i, lpi, frame
    real*8                      :: max_pos_eval, max_neg_eval, max_mid_eval 
    real*8                      :: scalar_prod, vector_length, angle_histo
    real*8                      :: max_mid 
    real*8, dimension(3)        :: pos_strain_dir, mid_strain_dir, neg_strain_dir 
    real*8                      :: pos_strain_dir_length, mid_strain_dir_length, neg_strain_dir_length
    real*8, dimension(5,max_lp) :: angle 
    real*8, dimension(5)        :: gamma_line_l2 
    real*8, dimension(5)        :: gamma_surf_l2 
    real*8, dimension(n)        :: vorticity, vorticity_normed
    real*8, dimension(n)        :: eval
    real*8, dimension(n,n)      :: mat_buffer !storing the VG for eingenv* computation
    real*8, dimension(n,n)      :: strain_rate, rotation_rate !tensors
    real*8, dimension(n,n)      :: evec
    real*8, dimension(n)        :: eval_l2_sum


    do lpi = 1 , proc_particles!loop over all particle
        mat_buffer(:,:) = lp_vgr_local(lpi,2,:,:) 

        !computing the symmetric part of the VG 
        call mat_decomp(mat_buffer, strain_rate, rotation_rate)
        !computing the eigen vectors and values
        call eigen_val_vec(strain_rate, eval, evec)
        !computes the maximum pos and neg strain rate and the directions
        call max_sort(eval, evec, pos_strain_dir, mid_strain_dir, neg_strain_dir, max_pos_eval, max_mid_eval, max_neg_eval)
        
        !eigenvalue statistics 
        eval_mean(1, frame) = eval_mean(1, frame)*(lpi-1)/real(lpi) + max_pos_eval/real(lpi)
        eval_mean(2, frame) = eval_mean(2, frame)*(lpi-1)/real(lpi) + max_mid_eval/real(lpi)
        eval_mean(3, frame) = eval_mean(3, frame)*(lpi-1)/real(lpi) + max_neg_eval/real(lpi)

        eval_l2_sum(1) = eval_l2_sum(1)*(lpi-1)/real(lpi) + max_pos_eval*max_pos_eval/real(lpi)
        eval_l2_sum(2) = eval_l2_sum(2)*(lpi-1)/real(lpi) + max_mid_eval*max_mid_eval/real(lpi)
        eval_l2_sum(3) = eval_l2_sum(3)*(lpi-1)/real(lpi) + max_neg_eval*max_neg_eval/real(lpi)

        eval_var(1, frame) = eval_var(1, frame)*(lpi-1)/real(lpi) + max_pos_eval/real(lpi)
        eval_var(2, frame) = eval_var(2, frame)*(lpi-1)/real(lpi) + max_mid_eval/real(lpi)
        eval_var(3, frame) = eval_var(3, frame)*(lpi-1)/real(lpi) + max_neg_eval/real(lpi)
        
        !unity normalization for max positive strain rate
        pos_strain_dir(:) = pos_strain_dir(:)/vector_length(pos_strain_dir(:))
        !unity normalization for intermediate strain rate
        mid_strain_dir(:) = mid_strain_dir(:)/vector_length(mid_strain_dir(:))
        !unity normalization for max negativ strain rate
        neg_strain_dir(:) = neg_strain_dir(:)/vector_length(neg_strain_dir(:))

        vorticity(1) = - 0.5_8*(rotation_rate(2,3) - rotation_rate(3,2)) 
        vorticity(2) = - 0.5_8*(rotation_rate(3,1) - rotation_rate(1,3)) 
        vorticity(3) = - 0.5_8*(rotation_rate(1,2) - rotation_rate(2,1)) 
        
        !++++++++++++++++++++++++++++++line stats+++++++++++++++++++++++++++++++

        !mean angle between max pos strain rate and line element 
        angle(1,lpi) =  acos(abs(dot_product(le_local(lpi,2,:,1), pos_strain_dir(:)))&
                          /(le_length_local(lpi,2,1)))
        !mean angle between intermediate strain rate and line element 
        angle(2,lpi) =  acos(abs(dot_product(le_local(lpi,2,:,1), mid_strain_dir(:)))&
                          /(le_length_local(lpi,2,1)))
        !angle between max strain rate and line element 
        angle(3,lpi) =  acos(abs(dot_product(le_local(lpi,2,:,1), neg_strain_dir(:)))&
                          /(le_length_local(lpi,2,1)))
        !angle between vorticity and line element 
        angle(4,lpi) =  acos(abs(dot_product(le_local(lpi,2,:,1), vorticity(:)))&
                          /(le_length_local(lpi,2,1)*vector_length(vorticity(:))))
        if(mhd == 1)then
            !angle between magnetic field and line element 
            angle(5,lpi) =  acos(abs(dot_product(le_local(lpi,2,:,1), mag_field_local(lpi,2,:)))&
                              /(le_length_local(lpi,2,1)* vector_length(mag_field_local(lpi,2,:))))
        else
            angle(5,lpi) = 0
        end if

        !mean angle
        gamma_line(frame,:) = gamma_line(frame,:)*(lpi-1)/real(lpi) + angle(:,lpi)/real(lpi)
        gamma_line_l2(:) = gamma_line_l2(:)*(lpi-1)/real(lpi) + angle(:,lpi)*angle(:,lpi)/real(lpi)

        !++++++++++++++++++++++++++++++surf stats+++++++++++++++++++++++++++++++

        !mean angle between max pos strain rate and line element 
        angle(1,lpi) =  acos(abs(dot_product(a_local(lpi,2,:), pos_strain_dir(:)))&
                          /(a_length_local(lpi,2)))
        !mean angle between intermediate strain rate and line element 
        angle(2,lpi) =  acos(abs(dot_product(a_local(lpi,2,:), mid_strain_dir(:)))&
                          /(a_length_local(lpi,2)))
        !angle between max strain rate and line element 
        angle(3,lpi) =  acos(abs(dot_product(a_local(lpi,2,:), neg_strain_dir(:)))&
                          /(a_length_local(lpi,2)))
        !angle between vorticity and line element 
        angle(4,lpi) =  acos(abs(dot_product(a_local(lpi,2,:), vorticity(:)))&
                          /(a_length_local(lpi,2)*vector_length(vorticity(:))))
        if(mhd == 1)then
            !angle between magnetic field and line element 
            angle(5,lpi) =  acos(abs(dot_product(a_local(lpi,2,:), mag_field_local(lpi,2,:)))&
                              /(a_length_local(lpi,2)* vector_length(mag_field_local(lpi,2,:))))
        else
            angle(5,lpi) = 0
        end if

        !mean angle
        gamma_surf(frame,:) = gamma_surf(frame,:)*(lpi-1)/real(lpi) + angle(:,lpi)/real(lpi)
        gamma_surf_l2(:) = gamma_surf_l2(:)*(lpi-1)/real(lpi) + angle(:,lpi)*angle(:,lpi)/real(lpi)

    end do

    eval_var(:,frame) = eval_l2_sum(:) - eval_mean(:,frame)*eval_mean(:,frame)   
    gamma_line_var(frame,:) = gamma_line_l2(:) - gamma_line(frame,:)*gamma_line(frame,:) 
    gamma_surf_var(frame,:) = gamma_surf_l2(:) - gamma_surf(frame,:)*gamma_surf(frame,:) 

    if (proc_id .eq. root_process) then 
        if(frame == start_frame+1)then 
            open(unit=102, file=angle_file, action='write', form='formatted')
        else
            open(unit=102, file=angle_file, action='write', form='formatted', position='append')
        end if 

            write(102,*) (frame-start_frame)*dt,&
                gamma_line(frame,:), gamma_line_var(frame,:),&
                gamma_surf(frame,:), gamma_surf_var(frame,:) 
            close(102)

        if(frame == start_frame+1)then 
            open(unit=104, file=strain_file, action='write', form='formatted')
        else
            open(unit=104, file=strain_file, action='write', form='formatted', position='append')
        end if 
                write(104,*) (frame-start_frame)*dt, eval_mean(1, frame),&
                                eval_mean(2, frame), eval_mean(3, frame)
        close(104)
    end if

end subroutine strain_stats

!----------------------------------------------------------------------------------------
! subroutine:   computes the maximum positive and negative strain rates
!               and determines the angle with respect to the line elements 
!----------------------------------------------------------------------------------------
subroutine angle_histograms(frame)
    use global
    use mpi
    implicit none
    integer                     :: i, lpi, frame
    real*8                      :: max_pos_eval, max_neg_eval, max_mid_eval 
    real*8                      :: scalar_prod, vector_length, angle
    real*8                      :: max_mid
    real*8, dimension(3)        :: pos_strain_dir, mid_strain_dir, neg_strain_dir 
    real*8                      :: pos_strain_dir_length, mid_strain_dir_length, neg_strain_dir_length
    real*8, dimension(5)        :: gamma_line_l2 
    real*8, dimension(5)        :: gamma_surf_l2 
    real*8, dimension(n)        :: eval
    real*8, dimension(n,n)      :: mat_buffer !storing the VG for eingenv* computation
    real*8, dimension(n,n)      :: strain_rate, rotation_rate !tensors
    real*8, dimension(n,n)      :: evec
    real*8, dimension(n)        :: eval_l2_sum
    real*8, dimension(n)        :: vorticity, vorticity_normed
    real*8, parameter           :: max_val = 1
    integer                     :: bin
    integer                     :: number_bin
    real*8, dimension(:,:), allocatable :: angle_histo, angle_histo_global

        !1 - vorticity_angle_histo
        !2 - positive_strain_angle_histo
        !3 - middle_strain_angle_histo
        !4 - negative_strain_angle_histo
        !5 - pos_strain_vorticity_angle_histo
        !6 - mid_strain_vorticity_angle_histo
        !7 - neg_strain_vorticity_angle_histo
        !8 - magnetic_angle_histo
        !9 - magnetic_pos_strain_angle_histo
        !10 - magnetic_mid_strain_angle_histo
        !11 - magnetic_neg_strain_angle_histo
        !12 - magnetic_vorticity_angle_histo

    if(frame == histo_frame)then
        number_bin = int(ceiling(2*max_val/bin_size_angle))

        allocate(angle_histo(number_bin,12), angle_histo_global(number_bin,12))
        angle_histo=0
        angle_histo_global=0

        eval_l2_sum = 0

        do lpi = 1 , proc_particles!loop over all particle
            mat_buffer(:,:) = lp_vgr_local(lpi,2,:,:) 

            !computing the symmetric part of the VG 
            call mat_decomp(mat_buffer, strain_rate, rotation_rate)
            !computing the eigen vectors and values
            call eigen_val_vec(strain_rate, eval, evec)
            !computes the maximum pos and neg strain rate and the directions
            call max_sort(eval, evec, pos_strain_dir, mid_strain_dir, neg_strain_dir, max_pos_eval, max_mid_eval, max_neg_eval)
            
            !unity normalization for max positive strain rate
            pos_strain_dir(:) = pos_strain_dir(:)/vector_length(pos_strain_dir(:))
            !unity normalization for intermediate strain rate
            mid_strain_dir(:) = mid_strain_dir(:)/vector_length(mid_strain_dir(:))
            !unity normalization for max negativ strain rate
            neg_strain_dir(:) = neg_strain_dir(:)/vector_length(neg_strain_dir(:))

            vorticity(1) = - 0.5_8*(rotation_rate(2,3) - rotation_rate(3,2)) 
            vorticity(2) = - 0.5_8*(rotation_rate(3,1) - rotation_rate(1,3)) 
            vorticity(3) = - 0.5_8*(rotation_rate(1,2) - rotation_rate(2,1)) 
            
            !++++++++++++++++++++++++++++++vorticity+++++++++++++++++++++++++++++++
            angle = abs( dot_product(vorticity(:), le_local(lpi,2,:,1))&
                /(vector_length(vorticity(:))*vector_length(le_local(lpi,2,:,1))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,1) = angle_histo(bin,1) + 1
            
            !+++++++++++++++++++++++max positive strain rate angle binning++++++++++++++++++++++++
            angle = abs( dot_product(pos_strain_dir(:), le_local(lpi,2,:,1))&
                /(vector_length(pos_strain_dir(:))*vector_length(le_local(lpi,2,:,1))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,2) = angle_histo(bin,2) + 1
            
            !+++++++++++++++++++++++max middle strain rate angle binning++++++++++++++++++++++++
            angle = abs( dot_product(mid_strain_dir(:), le_local(lpi,2,:,1))&
                /(vector_length(mid_strain_dir(:))*vector_length(le_local(lpi,2,:,1))))
            
            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,3) = angle_histo(bin,3) + 1

            !++++++++++++++++++++++max negative strain rate angle binning+++++++++++++++++++++++++
            angle = abs( dot_product(neg_strain_dir(:), le_local(lpi,2,:,1))&
                /(vector_length(pos_strain_dir(:))*vector_length(le_local(lpi,2,:,1))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,4) = angle_histo(bin,4) + 1

            !++++++++++++++++++++++ (max) positive strain voritcity angle binning+++++++++++++++++++++++++
            angle = abs( dot_product(pos_strain_dir(:), vorticity(:))&
                /(vector_length(pos_strain_dir(:))*vector_length(vorticity(:))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,5) = angle_histo(bin,5) + 1

            !++++++++++++++++++++++ middle strain voritcity angle binning+++++++++++++++++++++++++
            angle = abs( dot_product(mid_strain_dir(:), vorticity(:))&
                /(vector_length(mid_strain_dir(:))*vector_length(vorticity(:))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,6) = angle_histo(bin,6) + 1

            !++++++++++++++++++++++ negative strain voritcity angle binning+++++++++++++++++++++++++
            angle = abs( dot_product(neg_strain_dir(:), vorticity(:))&
                /(vector_length(neg_strain_dir(:))*vector_length(vorticity(:))))

            bin = ceiling((angle + max_val)/bin_size_angle)
            angle_histo(bin,7) = angle_histo(bin,7) + 1

            if(mhd == 1)then
                !++++++++++++++++++++++magnetic field line element angle binning+++++++++++++++++++++++++
                angle = abs( dot_product(mag_field_local(lpi,2,:), le_local(lpi,2,:,1))&
                    /(vector_length(mag_field_local(lpi,2,:))*vector_length(le_local(lpi,2,:,1))))

                bin = ceiling((angle + max_val)/bin_size_angle)
                angle_histo(bin,8) = angle_histo(bin,8) + 1

                !++++++++++++++++++++++magnetic pos strain field angle binning+++++++++++++++++++++++++
                angle = abs( dot_product(mag_field_local(lpi,2,:), pos_strain_dir(:))&
                    /(vector_length(mag_field_local(lpi,2,:))*vector_length(pos_strain_dir(:))))

                bin = ceiling((angle + max_val)/bin_size_angle)
                angle_histo(bin,9) = angle_histo(bin,9) + 1

                !++++++++++++++++++++++magnetic mid strain field angle binning+++++++++++++++++++++++++
                angle = abs( dot_product(mag_field_local(lpi,2,:), mid_strain_dir(:))&
                    /(vector_length(mag_field_local(lpi,2,:))*vector_length(mid_strain_dir(:))))

                bin = ceiling((angle + max_val)/bin_size_angle)
                angle_histo(bin,10) = angle_histo(bin,10) + 1

                !++++++++++++++++++++++magnetic neg strain field angle binning+++++++++++++++++++++++++
                angle = abs( dot_product(mag_field_local(lpi,2,:), neg_strain_dir(:))&
                    /(vector_length(mag_field_local(lpi,2,:))*vector_length(neg_strain_dir(:))))

                bin = ceiling((angle + max_val)/bin_size_angle)
                angle_histo(bin,11) = angle_histo(bin,11) + 1

                !++++++++++++++++++++++magnetic voritcity angle binning+++++++++++++++++++++++++
                angle = abs( dot_product(mag_field_local(lpi,2,:), vorticity(:))&
                    /(vector_length(mag_field_local(lpi,2,:))*vector_length(vorticity(:))))

                bin = ceiling((angle + max_val)/bin_size_angle)
                angle_histo(bin,12) = angle_histo(bin,12) + 1

            end if
        end do

        !----------------1D angle histogram----------------------------------------------
        do i=1, 7

            call MPI_REDUCE(angle_histo(:,i), angle_histo_global(:,i), number_bin, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

            if (proc_id .eq. root_process) then 
                angle_histo_global(:,i) = angle_histo_global(:,i)/(max_lp*bin_size_angle)
            end if

        end do


        if(mhd == 1)then

            do i=8, 12
                call MPI_REDUCE(angle_histo(:,i), angle_histo_global(:,i), number_bin, MPI_DOUBLE_PRECISION,&
                                MPI_SUM, root_process, MPI_COMM_WORLD, ierr)

                if (proc_id .eq. root_process) then 
                    angle_histo_global(:,i) = angle_histo_global(:,i)/(max_lp*bin_size_angle)
                end if
            end do

        end if


        if (proc_id .eq. root_process) then 
            open(unit=102, file=angle_histo_file, status='unknown')
                do i=(number_bin +2)/2, number_bin

                    write(102,*) i*bin_size_angle - (max_val + bin_size_angle/2.),&
                                angle_histo_global(i,2),&
                                angle_histo_global(i,3),&
                                angle_histo_global(i,4),&
                                angle_histo_global(i,1),&
                                angle_histo_global(i,5),&
                                angle_histo_global(i,6),&
                                angle_histo_global(i,7)
                end do
            close(102)

            if(mhd == 1)then
                open(unit=102, file=mhd_angle_histo_file, status='unknown')
                    do i=(number_bin +2)/2, number_bin
                        write(102,*) i*bin_size_angle - (max_val + bin_size_angle/2.),&
                                    angle_histo_global(i,8),&
                                    angle_histo_global(i,9),&
                                    angle_histo_global(i,10),&
                                    angle_histo_global(i,11),&
                                    angle_histo_global(i,12),&
                                    angle_histo_global(i,5),&
                                    angle_histo_global(i,6),&
                                    angle_histo_global(i,7)

                    end do
                close(102)
            end if
        end if

        deallocate(angle_histo, angle_histo_global)

    end if

end subroutine angle_histograms

function vector_length(vector)
    use global
    implicit none
    real*8 :: vector_length
    real*8, dimension(n) :: vector 

    vector_length = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)
end function

!----------------------------------------------------------------------------------------
! subroutine: sorts the eigenvectors and values in positive, middle and negative 
!----------------------------------------------------------------------------------------
subroutine max_sort(eval_in, evec_in, pos_in, mid_in, neg_in, max_pos, max_mid, max_neg)
    use global
    implicit none

    real*8, dimension(3) :: pos_in !max positive direction
    real*8, dimension(3) :: neg_in !max negative direction
    real*8, dimension(3) :: mid_in !mid strain direction
    real*8, dimension(3) :: eval_in 
    real*8, dimension(3,3) :: evec_in 

    real*8 :: max_pos 
    real*8 :: max_neg 
    real*8 :: max_mid 

    ! 2 negative eval_in, eval_in 1 and 2
    if(sign(1._8,eval_in(1)) == sign(1._8,eval_in(2)) .and. sign(1._8,eval_in(2)) .lt. 0)then
        max_pos = eval_in(3)
        pos_in(:) = evec_in(:,3)
        if(abs(eval_in(1)) > abs(eval_in(2)))then
           max_neg = eval_in(1) 
           neg_in(:) = evec_in(:,1)
           max_mid = eval_in(2) 
           mid_in(:) = evec_in(:,2)
        else
           max_neg = eval_in(2) 
           neg_in(:) = evec_in(:,2)
           max_mid = eval_in(1) 
           mid_in(:) = evec_in(:,1)
        end if
    end if
    ! 2 positive eval_in, eval_in 1 and 2
    if(sign(1._8,eval_in(1)) == sign(1._8,eval_in(2)) .and. sign(1._8,eval_in(2)) .gt. 0)then
        max_neg = eval_in(3)
        neg_in(:) = evec_in(:,3)
        if(abs(eval_in(1)) > abs(eval_in(2)))then
           max_pos = eval_in(1) 
           pos_in(:) = evec_in(:,1)
           max_mid = eval_in(2) 
           mid_in(:) = evec_in(:,2)
        else
           max_pos = eval_in(2) 
           pos_in(:) = evec_in(:,2)
           max_mid = eval_in(1) 
           mid_in(:) = evec_in(:,1)
        end if
    end if

    ! 2 negative eval_in, eval_in 2 and 3
    if(sign(1._8,eval_in(2)) == sign(1._8,eval_in(3)) .and. sign(1._8,eval_in(2)) .lt. 0)then
        max_pos = eval_in(1)
        pos_in(:) = evec_in(:,1)
        if(abs(eval_in(2)) > abs(eval_in(3)))then
           max_neg = eval_in(2) 
           neg_in(:) = evec_in(:,2)
           max_mid = eval_in(3) 
           mid_in(:) = evec_in(:,3)
        else
           max_neg = eval_in(3) 
           neg_in(:) = evec_in(:,3)
           max_mid = eval_in(2) 
           mid_in(:) = evec_in(:,2)
        end if
    end if
    ! 2 positive eval_in, eval_in 2 and 3
    if(sign(1._8,eval_in(2)) == sign(1._8,eval_in(3)) .and. sign(1._8,eval_in(2)) .gt. 0)then
        max_neg = eval_in(1)
        neg_in(:) = evec_in(:,1)
        if(abs(eval_in(2)) > abs(eval_in(3)))then
           max_pos = eval_in(2) 
           pos_in(:) = evec_in(:,2)
           max_mid = eval_in(3) 
           mid_in(:) = evec_in(:,3)
        else
           max_pos = eval_in(3) 
           pos_in(:) = evec_in(:,3)
           max_mid = eval_in(2) 
           mid_in(:) = evec_in(:,2)
        end if
    end if

    ! 2 negative eval_in, eval_in 1 and 3
    if(sign(1._8,eval_in(1)) == sign(1._8,eval_in(3)) .and. sign(1._8,eval_in(1)) .lt. 0)then
        max_pos = eval_in(2)
        pos_in(:) = evec_in(:,2)
        if(abs(eval_in(1)) > abs(eval_in(3)))then
           max_neg = eval_in(3) 
           neg_in(:) = evec_in(:,3)
           max_mid = eval_in(1) 
           mid_in(:) = evec_in(:,1)
        else
           max_neg = eval_in(1) 
           neg_in(:) = evec_in(:,1)
           max_mid = eval_in(3) 
           mid_in(:) = evec_in(:,3)
        end if
    end if
    ! 2 positive eval_in, eval_in 1 and 3
    if(sign(1._8,eval_in(1)) == sign(1._8,eval_in(3)) .and. sign(1._8,eval_in(1)) .gt. 0)then
        max_neg = eval_in(2)
        neg_in(:) = evec_in(:,2)
        if(abs(eval_in(1)) > abs(eval_in(3)))then
           max_pos = eval_in(1) 
           pos_in(:) = evec_in(:,1)
           max_mid = eval_in(3) 
           mid_in(:) = evec_in(:,3)
        else
           max_pos = eval_in(3) 
           pos_in(:) = evec_in(:,3)
           max_mid = eval_in(1) 
           mid_in(:) = evec_in(:,1)
        end if
    end if

    if(verbose == 5)then 
        print*, "max pos", max_pos, "max mid", max_mid, "max neg", max_neg
    end if
end subroutine max_sort

!----------------------------------------------------------------------------------------
! subroutine:   computes eigenvalues wr of a real double matrix, for real eigenvalues it
!               also needs to be symmertric. 
!----------------------------------------------------------------------------------------
subroutine eigen_val_vec(M, wr, vr)
    use global
    implicit none
    !--------------------lapack--------------------------- 
    integer, parameter :: lda = n !ld is leading dimension
    integer, parameter :: ldvr = n
    integer, parameter :: ldvl = n
    real*8, dimension(ldvl,n) :: vl 
    real*8, dimension(ldvr,n) :: vr !right side eigenvector 
    real*8, dimension(n) :: wr!real eigenvalues
    real*8, dimension(n) :: wi 
    integer :: lwork = 8*n
    real*8, dimension(n*4) :: work 
    character*1 :: jobvl = 'n' !compute left eigenvectors : n = no, V = yes
    character*1 :: jobvr = 'V' ! ""     right "" 
    integer :: info ! if = 0 then successful, if = - i, ith paramter had illegal value

    !--------------------eval--------------------------- 
    real*8, dimension(n,n) :: M ! input matrix 

    vl = 0 
    vr = 0 
    wr = 0
    wi = 0
    
    call dgeev(jobvl, jobvr, n, M, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

    if(verbose == 61)then 
        print*, ""
        print*, "----------------EValues--------------------"
        print*,"re eigenvalues: ",  wr(:)
        !print*,"im eigenvalues: ", wi(:)
        print*, ""
        print*, "-----------------EVectors--------------------"
        print*, vr(:,1)
        print*, vr(:,2)
        print*, vr(:,3)
    end if
End subroutine eigen_val_vec

subroutine mat_decomp(M, M_sym, M_asym)
    use global
    implicit none

    integer :: i, j 
    real*8, dimension(n,n) :: M ! input matrix 
    real*8, dimension(n,n) :: M_t ! transpose of M.
    real*8, dimension(n,n) :: M_sym ! symmetric part of M.
    real*8, dimension(n,n) :: M_asym ! anti symmetric part of M

    do i = 1, n 
        do j = 1, n 
            M_t(i,j) = M(j,i)
        end do
    end do

    M_sym = 0.5_8 *( M + M_t)
    M_asym = 0.5_8 *( M - M_t)

    if(verbose == 6)then 
        print*, ""
        print*, "-----------------Matrix decomposition--------------------"
        print*, "M"
        print*,"|", M(1,1:3) ,"|"
        print*,"|", M(2,1:3) ,"|"
        print*,"|", M(3,1:3) ,"|"

        print*, ""
        print*, "M_t"
        print*,"|", M_t(1,1:3) ,"|"
        print*,"|", M_t(2,1:3) ,"|"
        print*,"|", M_t(3,1:3) ,"|"

        print*, ""
        print*, "M_sym"
        print*,"|", M_sym(1,1:3) ,"|"
        print*,"|", M_sym(2,1:3) ,"|"
        print*,"|", M_sym(3,1:3) ,"|"

        M = M_sym + M_asym

        print*, ""
        print*, "M = M_sym + M_asym"
        print*,"|", M(1,1:3) ,"|"
        print*,"|", M(2,1:3) ,"|"
        print*,"|", M(3,1:3) ,"|"
    end if
end subroutine mat_decomp 

!----------------------------------------------------------------------------------------
! subroutine: measures the autocovariance of zeta and xi
!----------------------------------------------------------------------------------------
!subroutine correlation(corr_start_frame, corr_end_frame)
!    use global
!    implicit none
!    integer :: lpi,corr_end_frame, corr_start_frame, i,j
!    real*8, dimension(corr_end_frame - corr_start_frame)  :: zeta_cov, xi_cov 
!    real*8  ::  norm_sum 
!
!    zeta_cov(:) = 0
!    xi_cov(:) = 0
!
!    if(verbose == 10)then
!        print*, corr_start_frame, corr_end_frame
!    end if
!    
!    do i = 1,corr_end_frame - corr_start_frame !loop over all particle
!        do lpi = 1 , max_lp!loop over all particle
!            zeta_cov(i) = &
!            zeta_cov(i) + (zeta(lpi, corr_start_frame) - zeta_mean(corr_start_frame))*&
!                          (zeta(lpi, corr_start_frame+i-1) - zeta_mean(corr_start_frame+i-1))
!            xi_cov(i) = &
!            xi_cov(i) + (xi(lpi, corr_start_frame) - xi_mean(corr_start_frame))*&
!                        (xi(lpi, corr_start_frame+i-1) - xi_mean(corr_start_frame+i-1))
!        end do
!
!        zeta_cov(i) = zeta_cov(i)/max_lp 
!        xi_cov(i) = xi_cov(i)/max_lp 
!
!        if(i == 1)then 
!            open(unit=104, file=corr_file, action='write', form='formatted')
!        else
!            open(unit=104, file=corr_file, action='write', form='formatted', position='append')
!        end if 
!                write(104,*) (i-1)*dt, zeta_cov(i), xi_cov(i)
!            if(verbose == 10)then
!                write(*,*) (i-1)*dt, zeta_cov(i), xi_cov(i)
!            end if
!        close(104)
! 
!    end do
!end subroutine correlation
