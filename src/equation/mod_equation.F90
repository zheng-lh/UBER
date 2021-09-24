module mod_equation
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
      use mod_typedef_params, only: space, spacetime, STRLEN, BUFFER, NoData
      use mod_equation_typedef
      use linear_algebra
      implicit none
      private
      public coefficients, components, equation_cmpnt, initialize_equation, &
           & finalize_equation
! This module is the interface between the equation object and the SDE solver. It
! assembles the user-provided Christoffel symbols and the equation coefficients into
! the constructing components of Ito processes defined in mod_equation_typdef.F90.

!     module procedures to be specified in submodule mod_equation:user_input
      interface
         module function equation_djlnG(x) result(djlnG)
!        the Christoffel symbols dj(lnG) as a function of space
         class(space), intent(in) :: x
         real(8) :: djlnG(3)
         end function

         module function equation_coeff(xt) result(coeff)
!        the equation coefficients (see Eq. (2) in mod_equation_typedef.F90) as a
!        function of spacetime
         class(spacetime), intent(in) :: xt
         type(coefficients) :: coeff
         end function

         module subroutine equation_set_data()
!        user-specified routine to set coefficient data arrays if gridded data from
!        external data files are used. a phony procedure must be in place even if
!        no gridded data are actually used. to use the gridded data, first declare
!        in the header of submodule mod_equation:user_input, for example (2D):
!           use mod_grid
!           ...
!           type(grid_tensor_field) :: TT
!           ...
!        then in this subroutine:  
!           ...
!           character(STRLEN) :: fname11, fname12, fname22
!           integer :: nx1, nx2
!           real(8), allocatable :: xx1(:), xx2(:)
!           real(8), allocatable, dimension(:,:) :: data_t11, data_t12, data_t22
!           ...
!           call grid_io_read_data_file(fname11, nx1, nx2, xx1, xx2, data_t11)
!           call grid_io_read_data_file(fname12, nx1, nx2, xx1, xx2, data_t12)
!           call grid_io_read_data_file(fname22, nx1, nx2, xx1, xx2, data_t12)
!           TT = grid_tensor_field(nx1, nx2, xx1, xx2, data_t11, data_t12, data_t22)
!           deallocate( xx1, xx2, data_t11, data_t12, data_t22 )
!           ...
         end subroutine

         module subroutine equation_clean_data()
!        user-specified routine to clean up data set in equation_set_data. for
!        example:
!           call TT%clean()
         end subroutine
      end interface

contains

      function equation_cmpnt(xt) result(cmpnt)
!     see comments in mod_equation_typedef.F90 for explanation of the calculations
!     in this routine
      class(spacetime), intent(in) :: xt
      type(components) :: cmpnt
      type(coefficients) :: coeff
      real(8) :: djlnG(3), divD(3), D(3,3), divh

      djlnG = equation_djlnG(xt%space)
      coeff = equation_coeff(xt)

      call symmetric_matrix(3, coeff%Dij, D)
      divD = coeff%djDij + matmul(D, djlnG)
      divh = coeff%dihi + dot_product(coeff%hi, djlnG)

      cmpnt%A = 2.d0*D
      cmpnt%B = divD - coeff%hi
      cmpnt%h = coeff%hi
      cmpnt%c = coeff%S - divh
      cmpnt%v = coeff%v

      end function equation_cmpnt



      subroutine initialize_equation()
      print *,''
      print '(A)',' Initializing Boltzmann equation...'
      call equation_set_data()
      print *,''
      end subroutine initialize_equation



      subroutine finalize_equation()
      call equation_clean_data()
      print '(A)',' Equation object cleaned.'
      print *,''
      end subroutine finalize_equation

end module mod_equation

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) rewritten based on the orginal mod_equation and mod_equation_inc
!
! 2) user-specified functions will be contained in a submodule, interfaces of these
!    functions are provided here
!
! Liheng Zheng, 11/18/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *





