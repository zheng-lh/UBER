module mod_domain
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) moved declarations of integer nBoundaries and type(boundary) boundaries(:) from
!    the mod_domain_typedef module to this module
!
! 2) changed module mod_domain_ibc to a user-specified submodule, and added the
!    appropriate interfaces here
!
! Liheng Zheng, 11/15/2019
!
! 1) added interfaces to user-specified procedures domain_set_data and domain_clean_
!    data
!
! Liheng Zheng, 01/19/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, only: space, spacetime, STRLEN, BUFFER, NoData, &
          & epsBnd, Neighborhood
      use mod_domain_typedef
      implicit none
      save
      private 
      public boundary_condition, boundary, nBoundaries, boundaries, &
           & domain_initial_condition, domain_query, initialize_domain, &
           & finalize_domain
! This module is the interface between the domain objects and the SDE solver, so
! that the solver needs only 'use mod_domain' rather than looking into the concrete
! mod_domain_* modules.

!     number of boundary pieces
      integer :: nBoundaries

!     boundary objects
      type(boundary), allocatable, target :: boundaries(:)

!     module procedures to be specified in submodule mod_domain:user_input
      interface
         module function domain_initial_condition(x) result(init)
!        user-specified function to provide the initial condition to the problem        
         class(space), intent(in) :: x
         real(8) :: init
         end function

         module subroutine domain_set_boundaries()
!        user-specified routine to set value to nBoundaries, allocate the boundary
!        array, and associate boundary equation and boundary condition to each
!        boundary piece.
!        for example:
!           nBoundaries = 1
!           allocate( boundaries(nBoundaries) )
!           boundaries(1)%b_type = 2
!           boundaries(1)%equation => b1_equation
!           boundaries(1)%condition => b1_condition
         end subroutine

         module subroutine domain_set_data()
!        user-specified routine to set up gridded data. use of mod_grid is required.
!        a phony procedure must be in place even if no gridded data is actually
!        needed.
         end subroutine

         module subroutine domain_clean_data()
!        user-specified routine to clean up data set in domain_set_data
         end subroutine
      end interface

contains

      subroutine domain_query(xt, marker, bid)
!     this subroutine tells whether a spacetime point xt is in the interior of the
!     domain, or is in the inside neighborhood of boundary, or is out of the domain.
!     in the latter two cases, the corresponding boundary number is returned as bid.
!     this subroutine may only be called after the domain boundaries having been
!     initialized.
      type(spacetime), intent(in) :: xt
!     marker = 1, if xt is in the interior of the domain, in this case bid = 0
!              0, if xt is in the inside neighborhood of the boundary piece bid
!             -1, if xt is out of the boundary piece bid
      integer, intent(out) :: marker, bid
      real(8) :: equation_val
      integer :: i
      type(boundary), pointer :: bptr => null()
      !$OMP threadprivate(bptr)

      marker = 1
      bid = 0

      do i=1,nBoundaries
         bptr => boundaries(i)
         equation_val = bptr%equation(xt)
!        it is important to use epsBnd rather than 0.d0 here, because in
!        mod_solver_ito_proccess, either projection or intersection procedure
!        finds the root on a boundary within accuracy epsBnd. with epsBnd as the
!        criterion, the root will be regarded on the domain.
         if( equation_val<-epsBnd ) then
            marker = -1
            bid = i
            exit
         elseif( equation_val<Neighborhood ) then
            marker = 0
!           the neighborhood of a Neumann or Robin type boundary overrides that of a
!           Dirichlet type boundary. this is in order to deal with the situation
!           when, close to a corner of the domain, xt is in the neighborhood of both
!           a Dirichlet type boundary and a Neumann or Robin type boundary. the
!           half-space reflection algorithm in mod_solver_ito_process requires
!           information of the Neumann or Robin type boundary whenever the ito
!           process is in vicinity thereof.
            if( bid==0 ) then
               bid = i
            else
               if( bptr%b_type==2 ) bid = i
            end if
         end if
      end do

      nullify(bptr)
      end subroutine domain_query



      subroutine initialize_domain()
      print *,''
      print '(A)',' Initializing computational domain...'
      call domain_set_data()
      call domain_set_boundaries()
      print *,''
      end subroutine initialize_domain



      subroutine finalize_domain()
      if( allocated(boundaries) ) deallocate( boundaries )
      call domain_clean_data()
      print '(A)',' Domain object cleaned.'
      print *,''
      end subroutine finalize_domain

end module mod_domain




