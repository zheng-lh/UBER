module mod_solver_batch
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) made relevant adaptions to the new ito process module in this revision
!
! 2) removed the source distribution sorting function intended but unfinished in the
!    previous version
!
! 3) moved declaration of NBAT to mod_typedef_params.F90
!
! 4) improved the fission strategy. in the previous version, only one fission
!    generation is allowed, which is born when the partial functional value of the
!    running ito process is larger than the functional expectation by CriticalRatio
!    times of deviation. in this revision, multiple generations are allowed. the
!    fission criterion for each generation is when the partial functional value of
!    the ito process segment in this generation is (GENERATION + 1)*CriticalRatio
!    times of deviation larger than the functional expectation. GENERATION is 0 for
!    the root generation.
!
! Liheng Zheng, 01/08/2020
!
! 1) modified the implementations of ito process and tree in accordance to their
!    02/25/2020 edition. especially, the fission strategy no longer uses the partial
!    functional value, but uses the projected functional value of the entire path,
!    in the fission criterion.
!
! 2) added a field fission_number to type batch_statistics
!
! Liheng Zheng, 02/27/2020
!
! 1) changed the fission strategy. now fission occurs when the projected functional
!    value of an ito process is greater than some statistical quantile of previous
!    functional values. the quantile is specified by CriticalRatio, whose meaning
!    is also changed (see mod_typedef_params.F90).
!
! Liheng Zheng, 03/04/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, dummy_ParamFname => ParamFname, &
          & dummy_initialize_params => initialize_params
      use mod_domain, dummy_domain_initial_condition => domain_initial_condition, &
          & dummy_initialize_domain => initialize_domain, &
          & dummy_finalize_domain => finalize_domain
      use mod_solver_solution_grid, only: &
          & initial_condition => solver_initial_condition
      use mod_solver_tree
      use mod_solver_ito_process
      use statistics
      implicit none
      private
      public stat_moments, batch_statistics, solver_batch
! This module simulates one batch of ito processes from a common starting spacetime,
! and calclates their functional expectation as well as other statistics.

      type stat_moments
!        mean and standard deviation of functional F
         real(8) :: meanF
         real(8) :: devF
!        median path integral U
         real(8) :: medianU
!        median path integral V
         real(8) :: medianV
!        quantiles of F
         real(8) :: quantileF(MaxGenerations)
      end type stat_moments

      type batch_statistics
         integer :: fission_number
         real(8) :: mean_local_time
         real(8) :: mean_temporal_step
         real(8) :: mean_spatial_step
         real(8) :: max_spatial_step
         real(8), allocatable :: boundary_bins(:)
      end type batch_statistics       

contains

      subroutine solver_batch(xt_start, t_stop, moments_in, moments_out, batch_stat)
!     common starting spacetime for this batch
      type(spacetime), intent(in) :: xt_start
!     common ending time for this batch
      real(8), intent(in) :: t_stop
!     input statistical moments from previous calculations
      type(stat_moments), intent(in) :: moments_in
!     output statistical moments from this batch
      type(stat_moments), intent(out) :: moments_out
      type(batch_statistics), intent(out) :: batch_stat
!     a parameter used to control fission
      real(8), parameter :: Frac = 1.618d-1
      type(boundary), pointer :: bptr => null()
      type(boundary_condition) :: bcondition
      type(tree) :: my_tree
      type(ito_proc_params) :: my_ito_params
      type(ito_process) :: my_ito_process
      type(spacetime) :: my_xt_start, my_xt_end
      integer :: my_bid, fission_number
      integer :: i, j, k, fistat
      real(8) :: nbat_f, weight
      real(8) :: my_Ut, my_Vt
      real(8) :: dF, dU, dV
      real(8) :: sumF, sum_pathlen, sum_nsteps, max_spatial_step
      real(8) :: total_tau, total_kt
      real(8), dimension(NBAT) :: w_dF, w_dU, w_dV
      real(8), allocatable :: boundary_bins(:)
      logical :: success, terminal, fissile
      !$OMP threadprivate(bptr)

      allocate(batch_stat%boundary_bins(0:nBoundaries))
      allocate(boundary_bins(0:nBoundaries))

      sumF = 0.d0
      sum_pathlen = 0.d0
      sum_nsteps = 0.d0
      total_tau = 0.d0
      total_kt = 0.d0
      fission_number = 0
      boundary_bins = 0.d0
      max_spatial_step = 0.d0
      w_dF = 0.d0
      w_dU = 0.d0
      w_dV = 0.d0

      !$OMP parallel do schedule(static) default(private) copyin(bptr) &
      !$OMP shared(xt_start, t_stop, moments_in, w_dF, w_dU, w_dV) &
      !$OMP shared(dtMin, dtMax, dS, Neighborhood, epsBnd, CriticalRatio) &
      !$OMP shared(nChildren, nGenerations, nBoundaries, boundaries) &
      !$OMP reduction(+: sumF, sum_pathlen, sum_nsteps, total_tau, total_kt) &
      !$OMP reduction(+: fission_number, boundary_bins) &
      !$OMP reduction(MAX: max_spatial_step)
      do i=1,NBAT

!        create a tree structure of stochastic paths and copy the common starting
!        spacetime to the root node of the tree, i.e., "plant" the tree at xt_start
         my_tree = tree()
         call my_tree%set_node_value(u=1.d0, v=0.d0, xt=xt_start)

!        set the parameters for the ito process
         my_ito_params%t_total = xt_start%t - t_stop
         my_ito_params%medianU = moments_in%medianU
         my_ito_params%medianV = moments_in%medianV
         if( my_ito_params%t_total<dtMin ) then
            print *,' mod_solver_batch procedure solver_batch:'
            print *,'  error: xt_start%t (=',xt_start%t,') <= t_stop (=',t_stop,')'
            stop
         end if

!        let the tree grow!
         j = 0
         do while( .true. ) ! loop over generations
            do while( .true. ) ! loop over paths in one generation

!              for the 0th generation, my_xt_start is a copy of xt_start that has
!              been stored in the root node, and pastU=1 and pastV=0; for higher
!              generations, my_xt_start is the my_xt_end from its parent path stored
!              in the current node, and pastU and pastV are path integrals up to
!              the parent.
               call my_tree%get_node_value(u=my_ito_params%pastU, &
                                         & v=my_ito_params%pastV, xt=my_xt_start)

!              determine whether the ito process is fissile. the false cases are:
!                 1) nChildren==1    - fission is disabled
!                 2) j==nGenerations - currently in the last generation
!                 3) devF<Frac*meanF - the distribution is already concentrated
!                    enough, or there has not been enough statistics (devF<0)
               if( nChildren==1 .or. j==nGenerations .or. &
                  & moments_in%devF<Frac*moments_in%meanF ) then
                  fissile = .false.
                  my_ito_params%criticalF = -1.d0
               else
                  fissile = .true.
                  my_ito_params%criticalF = moments_in%quantileF(j+1)
               end if

!              simulate one path segment of the ito process at the current node
               my_ito_process = ito_process()
               call my_ito_process%walk(my_xt_start, t_stop, fissile, my_ito_params)
               my_xt_end = my_ito_process%xt_end
               my_Ut = my_ito_process%Ut
               my_Vt = my_ito_process%Vt
               my_bid = my_ito_process%bid

!              node is terminal if the ito process ends at t=t_stop or stops on a
!              Dirichlet type boundary -- use my_xt_end%x1 to store the initial or
!              boundary condition in these cases. N.B., an ito process segment must
!              have ended in either of these cases if not fissile.
               if( my_bid/=0 .or. (my_xt_end%t - t_stop)<dtMin ) then
                  terminal = .true.
                  if( my_bid==0 ) then
                     my_xt_end%x1 = initial_condition(my_xt_end)
                  else
                     bptr => boundaries(my_bid)
                     if( bptr%b_type==1 ) then
                        bcondition = bptr%condition(my_xt_end)
                        my_xt_end%x1 = bcondition%g
                     else
                        print *,' mod_solver_batch procedure solver_batch:'
                        print *,'  Ito process should not have stopped on'
                        print *,'  non-Dirichlet-type boundaries.'
                        print *,'  bid =',my_bid
                        print *,'  bptr%b_type =',bptr%b_type
                        print *,'  Ito process fission tree:'
                        call my_tree%dump()
                        stop
                     end if
                     nullify(bptr)
                  end if
                  my_xt_end%x2 = NoData
                  my_xt_end%x3 = NoData
                  my_xt_end%t = NoData
               else
                  terminal = .false.
               end if

!              aggregate the path integrals (c.f. mod_solver_tree.F90), and store
!              them along with the ending spacetime and terminal status to this node
               my_Ut = my_ito_params%pastU * my_Ut
               my_Vt = my_ito_params%pastV + my_ito_params%pastU * my_Vt
               call my_tree%set_node_value(my_Ut, my_Vt, my_xt_end, terminal)

!              update statistics
               weight = 1.d0/dble(nChildren**j)
               total_tau = total_tau + my_ito_process%tau*weight
               total_kt = total_kt + my_ito_process%kt*weight
               sum_nsteps = sum_nsteps + dble(my_ito_process%nsteps)
               sum_pathlen = sum_pathlen + my_ito_process%pathlen
               max_spatial_step = max(max_spatial_step, my_ito_process%maxstep)
               if( terminal ) boundary_bins(my_bid) = boundary_bins(my_bid) + weight

!              move to next tree branch in this generation
               success = my_tree%move_to_next_path()
               if( .not.success ) exit
            end do

!           fission all fissile ito processes in this generation
            fistat = my_tree%fission()
            if( fistat==0 ) then
!              tree had not grown, which means it is already in the last generation
               exit
            else
!              tree grew, and the internal cursor in the tree had been moved to the
!              first branch in the new generation.
               j = j + 1
               if( j>nGenerations ) then
                  print *,' mod_solver_batch procedure solver_batch:'
                  print *,'  invalid generation (=',j,') > nGenerations (=', &
                          & nGenerations,')'
                  print *,'  Ito process fission tree:'
                  call my_tree%dump()
                  stop
               end if
!              count once for each tree
               if( j==1 ) fission_number = fission_number + 1
            end if
         end do

!        trace the tree for the functional value F, path integrals U and V of the
!        ito process from xt_start. the traces are equivalent to g(xt_end)*Ut + Vt,
!        Ut and Vt, respectively, as if the ito process had not fissioned, where g
!        is the initial or Dirichlet-type-boundary condition, and Ut and Vt are path
!        integrals as defined in mod_solver_ito_process.F90.
         call my_tree%trace(dF, dU, dV)
         sumF = sumF + dF

!        store dF, dU, dV into working arrays for later statistical analysis
         w_dF(i) = dF
         w_dU(i) = dU
         w_dV(i) = dV

!        destroy the tree to release memory
         call my_tree%destroy()

      end do
      !$OMP end parallel do

!     calcualte the statistics of functional value and path integrals
      nbat_f = dble(NBAT)
      moments_out%meanF = sumF/nbat_f
      moments_out%devF = -1.d0
      moments_out%medianU = selection((NBAT+1)/2, NBAT, w_dU)
      moments_out%medianV = selection((NBAT+1)/2, NBAT, w_dV)
      moments_out%quantileF = -1.d0
      do j=1,nGenerations
!        for the j-th quantile of dF's, select the (1 - CriticalRatio/2**(j-1)) *
!        nbat_f smallest value from w_dF
         k = nint( (1.d0 - CriticalRatio/dble(2**(j-1)))*nbat_f )
         moments_out%quantileF(j) = selection(k, NBAT, w_dF)
      end do

!     calculate the auxiliary statistics of ito processes
      batch_stat%mean_local_time = total_kt/nbat_f
      batch_stat%mean_temporal_step = total_tau/sum_nsteps
      batch_stat%mean_spatial_step = sum_pathlen/sum_nsteps
      batch_stat%fission_number = fission_number
      batch_stat%boundary_bins = boundary_bins
      batch_stat%max_spatial_step = max_spatial_step

      deallocate(boundary_bins)

      end subroutine solver_batch

end module mod_solver_batch






