module mod_solver_solution_grid
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
      use mod_typedef_params, only: space, spacetime, dtMin
      use mod_domain, only: domain_initial_condition
      use mod_grid, only: grid_scalar_field
      implicit none
      save
      private
      public grid_solution, solution_prev, solver_initial_condition
! this module defines the grid_solution data type, and holds the solution_prev data
! structure in this type used to store the gridded solution from a previous time
! stamp. solution_prev is then used to provide the initial condition for the ito
! processes.

      type, extends(grid_scalar_field) :: grid_solution
         real(8), private :: t
      contains
         procedure :: get_t => grid_solution_get_t
         procedure :: update => grid_solution_update
         final :: grid_solution_finalization
      end type grid_solution

      type(grid_solution) :: solution_prev

      interface grid_solution
         module procedure grid_solution_initialization_1d, &
                        & grid_solution_initialization_2d, &
                        & grid_solution_initialization_3d
      end interface grid_solution

contains

      function grid_solution_initialization_1d(nx1, xx1, darr)
      integer, intent(in) :: nx1
      real(8), intent(in) :: xx1(:)
      real(8), intent(in) :: darr(:)
      type(grid_solution) :: grid_solution_initialization_1d

      grid_solution_initialization_1d%grid_scalar_field = grid_scalar_field(nx1, &
                                                        & xx1, darr, 'P')
      grid_solution_initialization_1d%t = 0.d0

      end function grid_solution_initialization_1d



      function grid_solution_initialization_2d(nx1, nx2, xx1, xx2, darr)
      integer, intent(in) :: nx1, nx2
      real(8), intent(in) :: xx1(:), xx2(:)
      real(8), intent(in) :: darr(:,:)
      type(grid_solution) :: grid_solution_initialization_2d

      grid_solution_initialization_2d%grid_scalar_field = grid_scalar_field(nx1, &
                                                        & nx2, xx1, xx2, darr, 'P')
      grid_solution_initialization_2d%t = 0.d0

      end function grid_solution_initialization_2d



      function grid_solution_initialization_3d(nx1, nx2, nx3, xx1, xx2, xx3, darr)
      integer, intent(in) :: nx1, nx2, nx3
      real(8), intent(in) :: xx1(:), xx2(:), xx3(:)
      real(8), intent(in) :: darr(:,:,:)
      type(grid_solution) :: grid_solution_initialization_3d

      grid_solution_initialization_3d%grid_scalar_field = grid_scalar_field(nx1, &
                                              & nx2, nx3, xx1, xx2, xx3, darr, 'P')
      grid_solution_initialization_3d%t = 0.d0

      end function grid_solution_initialization_3d



      subroutine grid_solution_get_t(this, t)
      class(grid_solution) :: this
      real(8), intent(out) :: t
      t = this%t
      end subroutine grid_solution_get_t



      subroutine grid_solution_update(this, t, sol_list)
      class(grid_solution) :: this
      real(8), intent(in) :: t
      real(8), intent(in) :: sol_list(:)
      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3 )

         if( size(sol_list)==nx1*nx2*nx3 ) then
            this%t = t
            this%aa = reshape(sol_list, (/nx1, nx2, nx3/))
            call this%set_mask()
         else
            print *,' mod_solver_grid_solution_typedef procedure'// &
                  & ' grid_solution_update:'
            print *,'  size mismatch in input list:'
            print *,'  nx1*nx2*nx3 =',nx1*nx2*nx3
            print *,'  size(sol_list) =',size(sol_list)
            stop
         end if

      end associate
      end subroutine grid_solution_update



      subroutine grid_solution_finalization(this)
      type(grid_solution) :: this
      call this%clean()
      end subroutine grid_solution_finalization



      function solver_initial_condition(xt) result(init)
!     this function provides initial condition for the ito processes. when the
!     stopping time is 0, initial condition is obtained from that provided by the
!     user in the domain submodule; when the stopping time is greater than 0,
!     solution from the last time stamp is regarded as the new initial condition.
      type(spacetime), intent(in) :: xt
      real(8) :: init
      real(8) :: t_stop

!     use dtMin as the epsilon in time
      if( xt%t<dtMin ) then
!        at t = 0, call domain_initial_condition directly
         init = domain_initial_condition(xt%space)
      else
!        at t = t_stop > 0, use the previous solution as the initial condition for
!        the next solution
         call solution_prev%get_t(t_stop)
         if( abs(xt%t - t_stop)<dtMin ) then
            call solution_prev%get(xt%space, init)
!           the following line is solely for the neutron-3D test purpose, it should
!           remain commented out at all other times
!@            call solution_prev%get(space(sqrt(xt%x1**2 + xt%x2**2 + xt%x3**2), &
!@                                 & 0.d0, 0.d0), init)
         else
            print *,' mod_solver_solution_grid procedure solver_initial_condition:'
            print '(2(A,F10.4))','  input xt%t =',xt%t,' is different from'//&
            & ' t_stop =',t_stop
            stop
         end if
      end if

      end function solver_initial_condition

end module mod_solver_solution_grid

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 01/02/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




