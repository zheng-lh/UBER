module mod_solver
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 01/10/2020
!
! 1) made solution_prev (from mod_solver_solution_grid) public for its access from
!    submodule mod_equation:user_input. solution_prev will be needed in solving
!    non-diffusive integro-differential kinetic equations.
!
! Liheng Zheng, 10/24/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, only: space, spacetime, NoData
      use mod_domain, only: domain_query, domain_initial_condition
      use mod_solver_random_number, only: init_dcmt, free_dcmt
      use mod_solver_solution_grid, &
          & dummy_solver_initial_condition => solver_initial_condition
      use mod_solver_io
      use mod_solver_engine
      implicit none
      save
      private
      public initialize_solver, finalize_solver, solver, solver_engine, &
           & solution_prev
! This module contains the solver of Eq. (1) given in mod_equation_typedef.F90. It
! is used as:
!    ...
!    call initialize_solver(FNAME_IN, FNAME_OUT)
!    call solver(ncp=NCP, icp=ICP, istat=ISTAT)
!    call finalize_solver()
!    ...
! where FNAME_IN and FNAME_OUT are the character strings of the solution input file
! name and the output file name. The spacetime positions contained in FNAME_IN will
! be solved, and the solutions and errors will be writtin in FNAME_OUT.
!
! Alternatively, this module also provides the lower level solving procedure, which
! solves Eq. (1) at just one spacetime, whose usage is:
!    ...
!    call initialize_solver()
!    call solver_engine(xt_start, 0.d0, solution, error, istat=ISTAT)
!    call finalize_solver()
!    ...
! The optional input character ISTAT = 'V' turns on the verbose mode of solution
! statistics output, and ISTAT with any other character or being absent turns the
! verbose mode off.
!
! The two optional input integers NCP and ICP in procedure solver are designed to
! deal with periodic boundary conditions (at most one pair). The periodic boundary
! condition is not really considered a boundary condition internally; rather, it is
! dealt with by extending the domain to include multiple periods, so that the ito
! processes would not move out of the domain (except stopping on other Dirichlet
! type boundaries) within the given time duration. An example 2D domain with one
! periodic coordinate (section 2, 0 ~ 2*pi) is illustrated below:
!
!     +- - - - -+----------+- - - - -+
!     |         |          |         |
!               |          |
!     |    1    |     2    |    3    |
!               |          |
!     |         |          |         |
!     +- - - - -+----------+- - - - -+
!   -2*pi       0         2*pi      4*pi
!
! This domain is extended one extra period toward both ends (sections 1 and 3), and
! the initial and boundary conditions are copied to the extended sections
! accordingly. The enlarged domain (sections 1 to 3) along with the repeated initial
! and boundary conditions are what the user needs to provide explicitly, although
! the actual way of extension depends on the problem and time length of solution,
! and may be different from that illustrated here. By the way, the periodic
! boundaries of the enlarged domain, though ineffective, are recommended to take
! Dirichlet type boundary conditions with values of the initial condition for
! simplicity.
!
! In the 1st I/O mode, the requested solution positions are restricted within
! section 2. In this case, the input NCP and ICP are irrelevant and ignored. In the
! 2nd I/O mode, a solution grid is used. This solution grid must cover the entire
! enlarged domain. However, only a portion of the solution grid lies within the
! section to be solved, and NCP and ICP are used to tell the solver which portion it
! is. NCP indicates how many copies of the original domain (inclusive) are contained
! in the enlarged domain; and ICP indicates which copy is the original one. In the
! exemplifed case, NCP = 3 and ICP = 2 respectively. If only NCP is provided, ICP
! assumes the default value 1. After solving the ICP-th portion, the solution is
! copied to the other portions to populate the entire solution grid.
! 
! There are several implications of the aforementioned method. (1) The order how
! solution grid nodes are counted matters. Linear storage of the multi-dimensional
! grid (as in memory) must first fully cover section 1, then section 2, etc..
! Therefore in the example the 2D grid must be counted from the upper-left corner
! and in column major order. (2) The total number of grid nodes must be divisible
! by NCP. And the periodic boundaries must lie through the middle of grid cells
! rather than the grid nodes. 
!
! In the absence of NCP and ICP, they both take the default value 1.

      type(spacetime), allocatable :: xt_list(:)
      integer, allocatable :: mask_list(:)
      real(8), allocatable :: sol_list(:), err_list(:)

      interface initialize_solver
         module procedure initialize_solver, initialize_solver_engine
      end interface

contains

      subroutine initialize_solver(fname_in, fname_out)
      character(*), intent(in) :: fname_in, fname_out
      integer :: i, j, k, p
      integer :: marker, bid
      real(8), allocatable :: nodata_1d(:), nodata_2d(:,:), nodata_3d(:,:,:)

      print '(/A)',' Initializing the SDE solver...'

!     check the operating system, which must be unix
      call solver_io_os_check()

!     read the solution input file and open the output file. allocate the data lists
      call solver_io_read_input_file(fname_in)
      call solver_io_open_output_file(fname_out)

      allocate( xt_list(np), mask_list(np), sol_list(np), err_list(np) )
      xt_list = spacetime(0.d0, 0.d0, 0.d0, 0.d0)
      mask_list = 0
      sol_list = NoData
      err_list = NoData

!     assign values to xt_list
      select case( io_mode )
      case( 1 )
         do i=1,np
            xt_list(i)%x1 = x1_list(i)
            xt_list(i)%x2 = x2_list(i)
            xt_list(i)%x3 = x3_list(i)
            xt_list(i)%t = t_list(i)
         end do
      case( 2 )
!        xt_list%t deliberately remains 0
         do k=1,io_nx3
            do j=1,io_nx2
               do i=1,io_nx1
                  p = (k-1)*io_nx2*io_nx1 + (j-1)*io_nx1 + i
                  xt_list(p)%x1 = x1_list(i)
                  xt_list(p)%x2 = x2_list(j)
                  xt_list(p)%x3 = x3_list(k)
               end do
            end do
         end do
      case default
         print *,' mod_solver procedure initialize_solver:'
         print *,'  unrecognized I/O mode =',io_mode
         call solver_io_close_output_file()
         stop
      end select

!     assign mask_list, and initialize sol_list with initial condition and err_list
!     with 0. export the initial condition as a solution at t = 0 to the output file
      do i=1,np
         call domain_query(xt_list(i), marker, bid)
         if( marker>=0 ) then
            mask_list(i) = 1
            sol_list(i) = domain_initial_condition( xt_list(i)%space )
            err_list(i) = 0.d0
         end if
      end do

      call solver_io_export(sol_list, err_list)

!     if io_mode = 2, initialize solution_prev with the initial values in sol_list.
!     solution_prev is a positive semi-definite scalar field, and is used to hold
!     the solution at a previous time stamp.
      if( io_mode==2 ) then
         select case( io_npd )
         case( 1 )
            allocate( nodata_1d(io_nx1) )
            nodata_1d = NoData
            solution_prev = grid_solution(io_nx1, x1_list, nodata_1d)
            deallocate( nodata_1d )
         case( 2 )
            allocate( nodata_2d(io_nx1, io_nx2) )
            nodata_2d = NoData
            solution_prev = grid_solution(io_nx1, io_nx2, x1_list, x2_list, &
                                        & nodata_2d)
            deallocate( nodata_2d )
         case( 3 )
            allocate( nodata_3d(io_nx1, io_nx2, io_nx3) )
            nodata_3d = NoData
            solution_prev = grid_solution(io_nx1, io_nx2, io_nx3, x1_list, &
                                        & x2_list, x3_list, nodata_3d)
            deallocate( nodata_3d )
         case default
            print *,' mod_solver procedure initialize_solver:'
            print *,'  unrecognized dimension number io_npd =',io_npd
            call solver_io_close_output_file()
            stop
         end select

         call solution_prev%update(0.d0, sol_list)
      end if

!     initialize the DCMT pseudo-random number generator
      call init_dcmt()

      end subroutine initialize_solver



      subroutine initialize_solver_engine()
      print '(/A)',' Initializing the SDE solver...'
!     only initialize the DCMT pseudo-random number generator for solver_engine
      call init_dcmt()
      end subroutine initialize_solver_engine



      subroutine solver(ncp, icp, istat)
      integer, intent(in), optional :: ncp, icp
      character(1), intent(in), optional :: istat
      character(1) :: xstat
      integer :: i, j, k, marker, bid
      integer :: mcp, jcp, i_sol, n_sol
      real(8) :: t_stop

      if( present(istat) ) then
         xstat = istat
      else
         xstat = 'C'
      end if

      select case( io_mode )
      case( 1 )

         t_stop = 0.d0
         do i=1,np
            if( mask_list(i)==1 ) then
               call solver_engine(xt_list(i), t_stop, sol_list(i), err_list(i), &
                                & xstat)
            else
               sol_list(i) = NoData
               err_list(i) = NoData
            end if
         end do

         call solver_io_export(sol_list, err_list)

      case( 2 )

         if( present(ncp) ) then
            if( ncp<1 ) then
               print *,' mod_solver_solve procedure solver:'
               print *,'  input ncp (=',ncp,') must be no less than 1'
               stop
            elseif( mod(np, ncp)/=0 ) then
               print *,' mod_solver_solve procedure solver:'
               print *,'  solution grid node number np (=',np, &
                       & ') is indivisible by input ncp (=',ncp,')'
               stop
            end if

            mcp = ncp

            if( present(icp) ) then
               if( icp<=ncp ) then
                  jcp = icp
               else
                  print *,' mod_solver_solve procedure solver:'
                  print *,'  input icp (=',icp,') must not exceed ncp (=',ncp,')'
                  stop
               end if
            else
               jcp = 1
            end if
         else
            if( present(icp) ) then
               print *,' mod_solver_solve procedure solver:'
               print *,'  input argument ncp absent while icp (=',icp,') provided'
               stop
            end if
            mcp = 1
            jcp = 1
         end if

!        determine which portion of sol_list is to be solved by solver_engine
         n_sol = np/mcp
         i_sol = n_sol*(jcp - 1) + 1

         t_stop = 0.d0
         do k=1,io_nt
!           update xt_list%t and mask_list (because Dirichlet type boundary geometry
!           may be time-variable)
            xt_list(:)%t = t_list(k)
            do i=1,np
               call domain_query(xt_list(i), marker, bid)
               if( marker>=0 ) then
                  mask_list(i) = 1
               else
                  mask_list(i) = 0
               end if
            end do
            
!           solve for sol_list in the jcp-th copy of the mcp copies
            do i=i_sol,i_sol+n_sol-1
               if( mask_list(i)==1 ) then
                  call solver_engine(xt_list(i), t_stop, sol_list(i), err_list(i), &
                                   & xstat)
               else
                  sol_list(i) = NoData
                  err_list(i) = NoData
               end if
            end do

!           copy the solved portion to the rest of sol_list
            do j=1,mcp
               if( j==jcp ) cycle
               sol_list(n_sol*(j-1)+1 : n_sol*j) = sol_list(i_sol : i_sol+n_sol-1)
               err_list(n_sol*(j-1)+1 : n_sol*j) = err_list(i_sol : i_sol+n_sol-1)
            end do

!           update stopping time
            t_stop = t_list(k)

!           export solutions and errors and update solution_prev
            call solver_io_export(sol_list, err_list)
            if( k<io_nt ) call solution_prev%update(t_stop, sol_list)
         end do

      case default
         print *,' mod_solver procedure solver:'
         print *,'  unrecognized I/O mode =',io_mode
         stop
      end select

      end subroutine solver



      subroutine finalize_solver()
      logical :: fopened

!     close output files
      inquire(unit=Funit_out, opened=fopened)
      if( fopened ) then
         call solver_io_close_output_file()
         if( io_mode==2 ) call solution_prev%clean()
      end if

!     release memory claimed by the pseudo-random number generator
      call free_dcmt()

      if( allocated(xt_list) ) deallocate( xt_list )
      if( allocated(mask_list) ) deallocate( mask_list )
      if( allocated(sol_list) ) deallocate( sol_list )
      if( allocated(err_list) ) deallocate( err_list )

      end subroutine finalize_solver

end module mod_solver



