module mod_domain_typedef
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
      use mod_typedef_params, only: spacetime
      implicit none
      public

!     for inhomogeneous boundary condition of the first type (Dirichlet)
!
!        u(t,x) = g1(t,x)                                                        (1)
!
!     g contains the value of g1, and n(3) are zeros; for homogeneous boundary
!     condition of the second type (Neumann) or the third type (Robin)
!
!        F(t,x,u) - g2(t,x)*u(t,x) = 0                                           (2)
!
!     where
!
!        F(t,x,u) = gamma(t,x).grad u(t,x) - eta(t,x)*u(t,x)                     (3)
!
!     is the outward flux on the boundary, with gamma(x) = n.A/2 (unnormalized
!     conormal vector, n is the unit inward normal of the boundary, and A is twice
!     the diffusion coefficient tensor), eta = n.h (h is the advection coefficient
!     vector), g stores the value of g2, and n(3) stores the vector n. physical
!     meaning of the boundary condition (2) is more easily seen when written as
!     F(t,x,u) = g2(t,x)*u(t,x), which says that the outward flux on the boundary is
!     proportional to the distribution function u, and the proportionality
!     coefficient is g2.
      type boundary_condition
         real(8) :: g
         real(8) :: n(3)
      end type boundary_condition

      type boundary
!        type 1: the first (Dirichlet) type boundary
!        type 2: the second (Neumann) or the third (Robin) type boundary
         integer :: b_type
!        equation > 0, if inside boundary
!                 = 0, if on boundary
!                 < 0, if outside boundary
!        for Dirichlet type boundaries, the equation function may be time-dependent,
!        meaning that the boundary geometry can be time-variable. for the second or
!        the third type boundaries, the boundary geometry must be static.
         procedure(f_x), pointer, nopass :: equation => null()
!        In case of the second or the third type boundary condition, the procedure
!        condition should be able to return the normal vector even when the passed
!        xt is not on the boundary but in its the neighborhood. this implies some
!        projection mechanism upon knowledge of the boundary shape.
         procedure(v_x), pointer, nopass :: condition => null()
      end type boundary

      abstract interface
         function f_x(xt)
         import :: spacetime
         class(spacetime), intent(in) :: xt
         real(8) :: f_x
         end function f_x

         function v_x(xt)
         import :: spacetime, boundary_condition
         class(spacetime), intent(in) :: xt
         type(boundary_condition) :: v_x
         end function v_x
      end interface

end module mod_domain_typedef

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) moved declarations of integer nBoundaries and type(boundary) boundaries(:) to
!    the mod_domain module.
!
! 2) moved declaration and evaluation of real(8) epsBnd to the mod_typedef_params
!    module.
!
! Liheng Zheng, 01/19/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




