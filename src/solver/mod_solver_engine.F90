module mod_solver_engine
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! Copyright 2021, Liheng Zheng
!
! This file is part of UBER.
!
!    UBER is free software: you can redistribute it and/or modify it under the
!    terms of the MIT License as published by Massachusetts Institute of
!    Technology. UBER is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
!    FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
!
!    You should have received a copy of the MIT License along with UBER. If not,
!    see <https://opensource.org/licenses/MIT>.
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, dummy_ParamFname => ParamFname, &
          & dummy_initialize_params => initialize_params
      use mod_domain, only: nBoundaries
      use mod_solver_batch
      use statistics
      implicit none
      private
      public solver_engine
! This module contains the engine of the solver, which solves Eq. (1) of mod_
! equation_typedef.F90 at a given spacetime, and returns the solution and estimated
! error. Statistics of the solving process are instantaneously sent to the standard
! output unit, and may also be exported in much greater detail into a file if the
! verbose mode is chosen. The verbose mode is enabled by specifying ISTAT='V',
! however it is not recommended for large-scale computations.

!     minimum number of batches. this number overrides MMin if MMin < MLeast.
      integer, parameter :: MLeast = 4
      
contains

      subroutine solver_engine(xt_start, t_stop, solution, error, istat)
      type(spacetime), intent(in) :: xt_start
      real(8), intent(in) :: t_stop
      real(8), intent(out) :: solution
      real(8), intent(out) :: error
      character(1), intent(in), optional :: istat
      type(stat_moments) :: moments, moments_bat
      type(batch_statistics) :: batch_stat
      real(8) :: running_meanF, running_varF
      real(8), allocatable :: w_meanF(:), w_medianU(:), w_medianV(:)
      real(8), allocatable :: w_quantileF(:,:)
      real(8), allocatable :: w_mean_local_time(:)
      real(8), allocatable :: w_mean_t_step(:), w_mean_s_step(:)
      real(8), allocatable :: boundary_bins(:)
      real(8) :: mean_local_time, mean_t_step, mean_s_step
      real(8) :: max_spatial_step
      real(8) :: sum_boundary_bins
      integer :: i, j, ntotal, fission_number
      logical :: print_msg, terminate, verbose

      if( present(istat) .and. (istat=='V' .or. istat=='v') ) then
         verbose = .true.
      else
         verbose = .false.
      end if

100   format(A,4(F10.4,A))
200   format(1X,A7,1X,A14,3X,A13)
300   format(1X,I7,1X,ES14.6,2X,F7.1)
310   format(1X,I7,1X,ES14.6,2X,A7)
      write(6,*)''
      write(6,'(A)')  ' The SDE solver is calculating functional expectation at:'
      write(6,fmt=100)' x1 =',xt_start%x1,', x2 =',xt_start%x2, &
                    & ', x3 =',xt_start%x3,', T =',xt_start%t,'.'
      write(6,*)''
      write(6,fmt=200)'Batch #','E[F(X)]','Rel. err. (%)'
      flush(6)

      allocate(w_meanF(MMax), w_medianU(MMax), w_medianV(MMax))
      allocate(w_quantileF(MMax, nGenerations))
      allocate(w_mean_local_time(MMax), w_mean_t_step(MMax), w_mean_s_step(MMax))
      allocate(boundary_bins(0:nBoundaries))

      w_meanF = 0.d0
      w_medianU = 0.d0
      w_medianV = 0.d0
      w_quantileF = 0.d0
      w_mean_local_time = 0.d0
      w_mean_t_step = 0.d0
      w_mean_s_step = 0.d0
      fission_number = 0
      boundary_bins = 0.d0
      max_spatial_step = 0.d0

      moments%meanF = 0.d0
      moments%devF = -1.d0
      moments%medianU = 0.d0
      moments%medianV = 0.d0
      moments%quantileF = -1.d0

      do j=1,MMax
         
         solution = NoData
         error = NoData

!        simulate one batch of ito processes (this routine allocates memory to
!        batch_stat)
         call solver_batch(xt_start, t_stop, moments, moments_bat, batch_stat)

         w_meanF(j) = moments_bat%meanF
         w_medianU(j) = moments_bat%medianU
         w_medianV(j) = moments_bat%medianV
         w_quantileF(j,:) = moments_bat%quantileF(1:nGenerations)
         w_mean_local_time(j) = batch_stat%mean_local_time
         w_mean_t_step(j) = batch_stat%mean_temporal_step
         w_mean_s_step(j) = batch_stat%mean_spatial_step
         fission_number = fission_number + batch_stat%fission_number
         boundary_bins = boundary_bins + batch_stat%boundary_bins
         max_spatial_step = max(max_spatial_step, batch_stat%max_spatial_step)

!        release the memory claimed in batch_stat
         deallocate( batch_stat%boundary_bins )

!        MLeast batches are forced (N.B., j must be at least 2 to call the routines
!        cfdint and avevar)
         if( j<MLeast ) cycle

!        evaluate the solution and its confidence interval (error)
         call cfdint(w_meanF(1:j), j, Ma, solution, error)

!        print message on screen
         print_msg = (j<100 .and. mod(j,10)==0) .or. &
                   & (j>=100 .and. mod(j,100)==0) .or. j==MMax
         if( print_msg ) then
            if( solution>0.d0 ) then
               write(6,fmt=300) j, solution, 100.d0*error/solution
            else
               write(6,fmt=310) j, solution, '-------'
            end if
            flush(6)
         end if

!        estimate the running statistical moments for use in the next iteration of
!        SOLVER_BATCH
         call avevar(w_meanF(1:j), j, running_meanF, running_varF)
         moments%meanF = running_meanF
         moments%devF = sqrt(running_varF)
         moments%medianU = selection((j+1)/2, j, w_medianU(1:j))
         moments%medianV = selection((j+1)/2, j, w_medianV(1:j))
         do i=1,nGenerations
            moments%quantileF(i) = selection((j+1)/2, j, w_quantileF(1:j,i))
         end do

!        when j>=MMin, determine the terminate condition
         if( j>=MMin ) then
            if( j<=MRed ) then
               terminate = (error<=0.5d0*Beta*solution)
            else
               terminate = (error<=Beta*solution)
            end if

            if( terminate ) then
               if( .not.print_msg ) then
                  if( solution>0.d0 ) then
                     write(6,fmt=300) j, solution, 100.d0*error/solution
                  else
                     write(6,fmt=310) j, solution, '-------'
                  end if
               end if
               write(6,'(A)')' Designated accuracy reached.'
               flush(6)
               exit
            end if
         end if
      end do

!     calculate the statistics of this solution
      if( j==(MMax+1) ) then
         write(6,'(A)')' Maximum number of batches completed.'
         j = MMax
      end if
      ntotal = j*NBAT

      sum_boundary_bins = sum(boundary_bins)
      call stat_sheet(j, verbose, mean_local_time, mean_t_step, mean_s_step)

      write(6,*)''
      write(6,'(A,I8)')   ' Total number of Ito processes: ',ntotal
      write(6,'(A,I8,A,F5.1,A)')' Fissioned Ito processes: ',fission_number,' ( ', &
                              & 100.d0*dble(fission_number)/dble(ntotal),'% )'
      write(6,'(A,F10.7)')' Mean local time: ',mean_local_time
      write(6,'(A,F10.7)')' Mean temporal step size: ',mean_t_step
      write(6,'(A,F10.7)')' Mean spatial step size: ',mean_s_step
      write(6,'(A,F10.7)')' Maximum spatial step size: ',max_spatial_step
      write(6,*)''
      write(6,'(A)')      ' Stopping locations:'
      write(6,'(A16,A23)') 'Location','Number ( Percent)'

400   format(A,I1,A,F10.1,A,F5.1,A)
      do i=0,nBoundaries
         write(6,fmt=400) ' boundary id [',i,']: ',boundary_bins(i), &
                        & ' ( ',100.d0*boundary_bins(i)/sum_boundary_bins,'% )'
      end do

      write(6,'(A)')'-------------------------------'
      flush(6)

      deallocate(w_meanF, w_medianU, w_medianV, w_quantileF, w_mean_local_time, &
               & w_mean_t_step, w_mean_s_step, boundary_bins)

      contains

         subroutine stat_sheet(n, fwrite, ave_local_time, ave_t_step, ave_s_step)
         integer, intent(in) :: n
         logical, intent(in) :: fwrite
         real(8), intent(out) :: ave_local_time, ave_t_step, ave_s_step
         real(8) :: aveF
         real(8) :: varF, var_local_time, var_t_step, var_s_step
         real(8) :: devF, dev_local_time, dev_t_step, dev_s_step
         real(8) :: percent_in_domain
         character(STRLEN) :: fname
         character(8) :: date
         integer :: fid, i

         call avevar(w_meanF(1:n), n, aveF, varF)
         devF = sqrt(varF)

         call avevar(w_mean_local_time(1:n), n, ave_local_time, var_local_time)
         dev_local_time = sqrt(var_local_time)

         call avevar(w_mean_t_step(1:n), n, ave_t_step, var_t_step)
         dev_t_step = sqrt(var_t_step)

         call avevar(w_mean_s_step(1:n), n, ave_s_step, var_s_step)
         dev_s_step = sqrt(var_s_step)

         if( fwrite ) then
!           export detailed solution statistics into a file
            percent_in_domain = 100.d0*boundary_bins(0)/sum_boundary_bins
            call solver_engine_name_str(fname)

            open(newunit=fid,file=fname,form='formatted',status='replace')
            write(fid,'(A)')' UBER solution statistics sheet'
            call date_and_time(date)
            write(fid,'(A)')' Produced on '//date(1:4)//'-'//date(5:6)//'-'// &
                          & date(7:8)
            write(fid,'(A)')' Data table format: (I4,2X,ES14.6,3F10.6)'
            write(fid,*)''
            write(fid,'(A,F10.4)')' x1 =',xt_start%x1
            write(fid,'(A,F10.4)')' x2 =',xt_start%x2
            write(fid,'(A,F10.4)')' x3 =',xt_start%x3
            write(fid,'(A,F10.4)')'  T =',xt_start%t
            write(fid,*)''
            write(fid,'(A,I8)')' Total number of Ito processes =',ntotal
            write(fid,'(A,F6.1,A)')' Percent of Ito processes stayed in domain =',&
                                 & percent_in_domain,'%'
            write(fid,'(A,F10.6)')' Maximum spatial step size =',max_spatial_step
            write(fid,*)''
            write(fid,'(6X,A14,3A10)') 'Solution','LocalTime','TimeStep','SpaceStep'
            write(fid,'(A4,2X,ES14.6,3F10.6)') 'Mean', aveF, ave_local_time, &
                                             & ave_t_step, ave_s_step
            write(fid,'(A5,1X,ES14.6,3F10.6)') 'StDev', devF, dev_local_time, &
                                             & dev_t_step, dev_s_step
            write(fid,fmt='(A6,2X)',advance='no') 'Batch#'
            do i=1,47
               write(fid,fmt='(A1)',advance='no') '-'
            end do
            write(fid,*)''

            do i=1,j
               write(fid,'(I4,2X,ES14.6,3F10.6)') i, w_meanF(i), &
                   & w_mean_local_time(i), w_mean_t_step(i), w_mean_s_step(i)
            end do

            do i=1,55
               write(fid,fmt='(A1)',advance='no') '-'
            end do
            write(fid,*)''

            close(fid)
         end if
         end subroutine stat_sheet

      end subroutine solver_engine



      subroutine solver_engine_name_str(name_str)
      character(STRLEN), intent(out) :: name_str
      character(12) :: date, time

      call date_and_time(date, time)
      name_str = 'STAT_'//date(1:8)//'_'//time(1:6)//'.txt'
      name_str = trim(name_str)

      end subroutine solver_engine_name_str
!
!
!
!      subroutine solver_engine_name_str(xt, name_str)
!      type(spacetime), intent(in) :: xt
!      character(STRLEN), intent(out) :: name_str
!      integer :: intx1, intx2, intx3, intt
!      integer :: decix1, decix2, decix3, decit
!
!      intx1 = int(xt%x1)
!      decix1 = nint(1.d3*(xt%x1 - aint(xt%x1)))
!      intx2 = int(xt%x2)
!      decix2 = nint(1.d3*(xt%x2 - aint(xt%x2)))
!      intx3 = int(xt%x3)
!      decix3 = nint(1.d3*(xt%x3 - aint(xt%x3)))
!      intt = int(xt%t)
!      decit = nint(1.d3*(xt%t - aint(xt%t)))
!
!500   format('X',I2.2,'p',I3.3,'Y',I2.2,'p',I3.3,'Z',I2.2,'p',I3.3,'T',I2.2,'p',I3.3)
!      write(unit=name_str,fmt=500) intx1, decix1, intx2, decix2, intx3, decix3, &
!                                 & intt, decit
!      name_str = trim(name_str)
!
!      end subroutine solver_engine_name_str

end module mod_solver_engine

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) removed the source distribution sorting function intended but unfinished in the
!    previous version
! 
! 2) removed the returned initial value from subroutine solver_engine
!
! 3) changed input argument ISTAT from logical to optional character(1)
!
! 4) changed MLeast from 5 to 4
!
! Liheng Zheng, 01/09/2020
!
! 1) added statistics of medianU, medianV and fission number per the changes in mod_
!    solver_batch, mod_solver_ito_process and mod_solver_tree
!
! Liheng Zheng, 02/27/2020
!
! 1) added quantiles of F per the changes in mod_solver_batch and mod_solver_ito_
!    process
!
! Liheng Zheng, 03/04/2020
!
! 1) changed naming convention of the diagnostic statistics file. it now uses date
!    and time as a part of its name.
!
! Liheng Zheng, 06/21/2024
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




