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
      use mod_grid
      use functions
      implicit none

      type(grid_scalar_field) :: bounce_loss_cone

contains

      module subroutine domain_set_data()
      character(STRLEN) :: fname = '../data/ex3_bounce_loss_cone.dat'
      integer :: nx1
      real(8), allocatable :: xx1(:)
      real(8), allocatable :: darr(:)

      call grid_io_read_data_file(fname, nx1, xx1, darr)
      bounce_loss_cone = grid_scalar_field(nx1, xx1, darr)
      deallocate( xx1, darr )

      end subroutine domain_set_data



      module function domain_initial_condition(x) result(init)
      class(space), intent(in) :: x
      real(8) :: init

      init = 0.d0

      end function domain_initial_condition



      module subroutine domain_set_boundaries()
!
!                 [2]
!        +- - - - -+---------+ pi/2     [i]: the i-th boundary piece
!        |         |         |
!   x1             |         |           x1: equatorial pitch angle a0
!   ^ [3]|         |         | [4]
!   |              |         |           x2: geomagnetic longitude phi
!   |    |         |         |
!   |    +- - - - -+---------+ pi/6
!   | -2*pi       [1]      2*pi
!   +--------> x2

      nBoundaries = 4
      print '(A,I3)',' Number of boundary pieces =', nBoundaries
      allocate( boundaries(nBoundaries) )

      boundaries(1)%b_type = 1
      boundaries(1)%equation => b1_equation
      boundaries(1)%condition => b1_condition

      boundaries(2)%b_type = 2
      boundaries(2)%equation => b2_equation
      boundaries(2)%condition => b2_condition

      boundaries(3)%b_type = 1
      boundaries(3)%equation => b3_equation
      boundaries(3)%condition => b3_condition

      boundaries(4)%b_type = 1
      boundaries(4)%equation => b4_equation
      boundaries(4)%condition => b4_condition

      contains

         function b1_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b1_equation
         type(space) :: phi
         real(8) :: blc_at_phi

         if( xt%x2>=0.d0 ) then
            phi%x1 = xt%x2
         else
            phi%x1 = xt%x2 + 2.d0*Pi
         end if

         phi%x2 = NoData
         phi%x3 = NoData

         call bounce_loss_cone%get(phi, blc_at_phi)

         b1_equation = xt%x1 - blc_at_phi

         end function b1_equation


         function b1_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b1_condition

         if( abs(b1_equation(xt))<=Neighborhood ) then
            b1_condition%g = 0.d0
            b1_condition%n = (/0.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b1_condition



         function b2_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b2_equation

!        N.B., xt%x1 = Pi/2 are singular points of the equation
         b2_equation = 0.5d0*Pi - xt%x1 - epsBnd

         end function b2_equation



         function b2_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b2_condition

         if( abs(b2_equation(xt))<=Neighborhood ) then
            b2_condition%g = 0.d0
            b2_condition%n = (/-1.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b2_condition



         function b3_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b3_equation

         b3_equation = xt%x2 + 2.d0*Pi

         end function b3_equation



         function b3_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b3_condition

         if( abs(b3_equation(xt))<=Neighborhood ) then
            b3_condition%g = 0.d0
            b3_condition%n = (/0.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b3_condition


         function b4_equation(xt)
         class(spacetime), intent(in) :: xt
         real(8) :: b4_equation

         b4_equation = 2.d0*Pi - xt%x2

         end function b4_equation


         function b4_condition(xt)
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: b4_condition

         if( abs(b4_equation(xt))<=Neighborhood ) then
            b4_condition%g = 0.d0
            b4_condition%n = (/0.d0, 0.d0, 0.d0/)
         else
            call domain_set_boundaries_error(xt)
         end if

         end function b4_condition



         subroutine domain_set_boundaries_error(xt)
         class(spacetime), intent(in) :: xt
         print *,' mod_domain:user_input procedure domain_set_boundaries:'
         print *,'  passed spacetime appears not near the boundary'
         print '(A8,4(F7.3,A1))','  xt = (',xt%x1,',',xt%x2,',',xt%x3,',',xt%t,')'
         stop
         end subroutine domain_set_boundaries_error

      end subroutine domain_set_boundaries



      module subroutine domain_clean_data()
      call bounce_loss_cone%clean()
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
      use mod_grid
      use functions 
      implicit none

      type(grid_scalar_field) :: bounce_loss_cone, dipole_lshell

contains

      module subroutine equation_set_data()
      character(STRLEN) :: fname1 = '../data/ex3_bounce_loss_cone.dat'
      character(STRLEN) :: fname2 = '../data/ex3_dipole_lshell.dat'
      integer :: nx1
      real(8), allocatable :: xx1(:)
      real(8), allocatable :: darr(:)

      call grid_io_read_data_file(fname1, nx1, xx1, darr)
      bounce_loss_cone = grid_scalar_field(nx1, xx1, darr)
      deallocate( xx1, darr )

      call grid_io_read_data_file(fname2, nx1, xx1, darr)
      dipole_lshell = grid_scalar_field(nx1, xx1, darr)
      deallocate( xx1, darr )

      end subroutine equation_set_data



      module function equation_djlnG(x) result(djlnG)
!
!        G(a0) = T(sin(a0))*sin(2*a0)
!        ln[G(a0)] = ln[T(sin(a0))] + ln[sin(a0)] + ln[cos(a0)] + ln(2)
!        dln(G)/da0 = (1/T)*(dT/dy)*cos(a0) + cot(a0) - tan(a0)
!
!     where
!
!        y = sin(a0)
!        T(y) = 1.380173 - 0.639693*y**(3/4)
!        dT/dy = -0.479770*y**(-1/4)
!
      class(space), intent(in) :: x
      real(8) :: djlnG(3)
      real(8) :: s, c, d1lnG

      s = sin(x%x1)
      c = cos(x%x1)
      d1lnG = -0.479770d0*s**(-0.25d0)*c/SchulzT(s) + c/s - s/c
      djlnG = (/d1lnG, 0.d0, 0.d0/)

      end function equation_djlnG



      module function equation_coeff(xt) result(coeff)
      class(spacetime), intent(in) :: xt
      type(coefficients) :: coeff
      real(8), parameter :: E = 0.304d0 ! MeV
      real(8), parameter :: L = 1.25d0  ! McIlwain's L
      type(space) :: phi
      real(8) :: blc_at_phi, lshell_at_phi
      real(8) :: s, c
!
!       Daa = 1.d-5*exp(92.55*(cos(a0)**4 - cos(aL)**4)) + 1.d-9 (s^-1)
!       dDaa/da0 = -3.702d-3*sin(a0)*cos(a0)**3*exp(92.55*(cos(a0)**4 - cos(aL)**4))
!                 (s^-1)
!
      real(8) :: Daa, daDaa, expf
!
!       Se/p**2 = 1.70d-12*(1.530 - E')**2/(Ls**2.7*sin(a0)) (c^3 cm^-3 MeV^-3 s^-1)
!
!     where Ls is dipole L-shell and E' = E/E0, therefore for 0.304 MeV electrons
!
!       Se/p**2 = 1.48d-12/(Ls**2.7*sin(a0)) (c^3 cm^-3 MeV^-3 s^-1)
!
      real(8) :: crand

      if( xt%x2>=0.d0 ) then
         phi%x1 = xt%x2
      else
         phi%x1 = xt%x2 + 2.d0*Pi
      end if

      phi%x2 = NoData
      phi%x3 = NoData

      call bounce_loss_cone%get(phi, blc_at_phi)
      call dipole_lshell%get(phi, lshell_at_phi)

      s = sin(xt%x1)
      c = cos(xt%x1)

      expf = exp(9.255d1*(c**4 - cos(blc_at_phi)**4))
      Daa = 1.d-5*expf + 1.d-9
      daDaa = -3.702d-3*s*c**3*expf

      crand = 1.48d-12/(s*lshell_at_phi**2.7d0)

      coeff%Dij = (/Daa, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
      coeff%djDij = (/daDaa, 0.d0, 0.d0/)
      coeff%hi = (/0.d0, Omega3(E, s, L), 0.d0/)
      coeff%dihi = 0.d0
      coeff%S = 0.d0
      coeff%v = crand

      end function equation_coeff



      module subroutine equation_clean_data()
      call bounce_loss_cone%clean()
      call dipole_lshell%clean()
      end subroutine equation_clean_data

end submodule user_input
!
! ==================================================================================


