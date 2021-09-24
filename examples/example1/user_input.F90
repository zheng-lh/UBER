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
!
! ========================= UBER user input source file ============================
!
! This file provides the user-specified functions and subroutines that define:
!
!     1) the computational domain, including the domain geometry, the initial
!        condition and the boundary conditions.
!
!     2) the Boltzmann equation, including the equation coefficients and the
!        Christoffel symbols (the gradient of the logrithm of the Jacobian
!        determinant) of the coordinate system.
!
! ------------------------------ Computational Domain ------------------------------
!
submodule(mod_domain) user_input
! This submodule has access to:
!
! Parameters:
!     integer, parameter :: STRLEN = 64
!     integer, parameter :: BUFFER = 128
!     real(8), parameter :: NoData = -1.d38
!     real(8), protected :: epsBnd
!
! Derived types:
!     type space
!     type spacetime
!     type boundary_condition
!     type boundary
!
! Variables:
!     integer :: nBoundaries
!     type(boundary), allocatable, target :: boundaries(:)
!
! Procedures:
!     module function domain_initial_condition(x) result(init)
!     class(space), intent(in) :: x
!     real(8) :: init
!     end function
!
!     module subroutine domain_set_boundaries()
!     ! for example:
!        ...
!        nBoundaries = 1
!        allocate( boundaries(nBoundaries) )
!        boundaries(1)%b_type = 2
!        boundaries(1)%equation => b1_equation
!        boundaries(1)%condition => b1_condition
!        ...
!     end subroutine
!
!     module subroutine domain_set_data()
!     ! for example, location of a 2D surface zz = f(xx1, xx2).
!     ! declare the following in the header of this submodule:
!        ...
!        use mod_grid
!        type(grid_scalar_field) :: zz
!        ...
!     ! then in this subroutine:
!        ...
!        character(STRLEN) :: fname
!        integer :: nx1, nx2
!        real(8), allocatable :: xx1(:), xx2(:)
!        real(8), allocatable :: data_zz(:,:)
!        ...
!        call grid_io_read_data_file(fname, nx1, nx2, xx1, xx2, data_zz)
!        zz = grid_scalar_field(nx1, nx2, xx1, xx2, data_zz)
!        deallocate( xx1, xx2, data_zz )
!        ...
!     end subroutine
!
!     module subroutine domain_clean_data()
!     ! continuing with the example in DOMAIN_SET_DATA()
!        ...
!        call zz%clean()
!        ...
!     end subroutine
!
! The user must provide all of the above procedures. Even if a procedure is not
! needed (e.g., DOMAIN_SET_DATA), a phony procedure that does nothing must be in
! place. The user may also define other procedures than those listed above in this
! submodule.
! 
! Uncomment the following lines to use the modules
!      use mod_grid
!      use functions
      implicit none

contains

      module function domain_initial_condition(x) result(init)
      class(space), intent(in) :: x
      real(8) :: init
      real(8), parameter :: Sigma = 1.d-1

      init = exp(-0.5d0*(x%x1/Sigma)**2)

      end function domain_initial_condition



      module subroutine domain_set_boundaries()

      nBoundaries = 2
      print '(A,I3)',' Number of boundary pieces =', nBoundaries
      allocate( boundaries(nBoundaries) )

      boundaries(1)%b_type = 2
      boundaries(1)%equation => b1_equation
      boundaries(1)%condition => b1_condition

      boundaries(2)%b_type = 2
      boundaries(2)%equation => b2_equation
      boundaries(2)%condition => b2_condition

      contains

         function b1_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b1_equation

         b1_equation = xt%x1 - epsBnd

         end function b1_equation


         function b1_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b1_condition

         if( abs(b1_equation(xt))<=Neighborhood ) then
            b1_condition%g = 0.d0
            b1_condition%n = (/1.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b1_condition


         function b2_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b2_equation
         real(8), parameter :: x1Max = 1.d0

         b2_equation = x1Max - xt%x1

         end function b2_equation


         function b2_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b2_condition

         if( abs(b2_equation(xt))<=Neighborhood ) then
            b2_condition%g = 0.5d0
            b2_condition%n = (/-1.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b2_condition


         subroutine domain_set_boundaries_error(xt)
         class(spacetime), intent(in) :: xt
         print *,' mod_domain:user_input procedure domain_set_boundaries:'
         print *,'  passed spacetime appears not near the boundary'
         print '(A8,4(F7.3,A1))','  xt = (',xt%x1,',',xt%x2,',',xt%x3,',',xt%t,')'
         stop
         end subroutine domain_set_boundaries_error

      end subroutine domain_set_boundaries



      module subroutine domain_set_data()
!     do nothing
      end subroutine domain_set_data



      module subroutine domain_clean_data()
!     do nothing
      end subroutine domain_clean_data

end submodule user_input

! ------------------------------- Boltzmann Equation -------------------------------

submodule(mod_equation) user_input
! This submodule has access to:
!
! Parameters:
!     integer, parameter :: STRLEN = 64
!     integer, parameter :: BUFFER = 128
!     real(8), parameter :: NoData = -1.d38
!
! Derived types:
!     type space
!     type spacetime
!     type coefficients
!     type components
!
! Procedures:
!     module function equation_djlnG(x) result(djlnG)
!     ! the Christoffel symbols dj(lnG) as a function of space
!     class(space), intent(in) :: x
!     real(8) :: djlnG(3)
!     end function
!
!     module function equation_coeff(xt) result(coeff)
!     ! the equation coefficients as a function of spacetime
!     class(spacetime), intent(in) :: xt
!     type(coefficients) :: coeff
!     end function
!
!     module subroutine equation_set_data()
!     ! for example, a rank-2 tensor coefficient field.
!     ! declare the following in the header of this submodule:
!        ...
!        use mod_grid
!        type(grid_tensor_field) :: TT
!        ...
!     ! then in this subroutine:
!        ...
!        character(STRLEN) :: fname11, fname12, fname22
!        integer :: nx1, nx2
!        real(8), allocatable :: xx1(:), xx2(:)
!        real(8), allocatable :: data_t11(:,:), data_t12(:,:), data_t22(:,:)
!        ...
!        call grid_io_read_data_file(fname11, nx1, nx2, xx1, xx2, data_t11)
!        call grid_io_read_data_file(fname12, nx1, nx2, xx1, xx2, data_t12)
!        call grid_io_read_data_file(fname22, nx1, nx2, xx1, xx2, data_t22)
!        TT = grid_tensor_field(nx1, nx2, xx1, xx2, data_t11, data_t12, data_t22)
!        deallocate( xx1, xx2, data_t11, data_t12, data_t22 )
!        ...
!     end subroutine
!
!     module subroutine domain_clean_data()
!     ! continuing with the example in EQUATION_SET_DATA()
!        ...
!        call TT%clean()
!        ...
!     end subroutine
!
! The user must provide all of the above procedures. Even if a procedure is not
! needed (e.g., EQUATION_SET_DATA), a phony procedure that does nothing must be in
! place. The user may also define other procedures than those listed above in this
! submodule.
!
! Uncomment the following lines to use the modules
!      use mod_grid
!      use functions 
      implicit none

contains

      module function equation_djlnG(x) result(djlnG)
      class(space), intent(in) :: x
      real(8) :: djlnG(3)

      djlnG = (/2.d0/x%x1, 0.d0, 0.d0/)

      end function equation_djlnG



      module function equation_coeff(xt) result(coeff)
      class(spacetime), intent(in) :: xt
      type(coefficients) :: coeff

      coeff%Dij = (/1.d-1, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
      coeff%djDij = (/0.d0, 0.d0, 0.d0/)
      coeff%hi = (/0.d0, 0.d0, 0.d0/)
      coeff%dihi = 0.d0
      coeff%S = 2.5d0
      coeff%v = 1.d-6/(1.d0 + xt%x1)

      end function equation_coeff



      module subroutine equation_set_data()
!     do nothing
      end subroutine equation_set_data



      module subroutine equation_clean_data()
!     do nothing
      end subroutine equation_clean_data

end submodule user_input
!
! ==================================================================================


