module mod_grid_tensor_field
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
      use mod_typedef_params, only: space, NoData
      use mod_grid_typedef
      use interpolation
      implicit none
      save
      private
      public grid_tensor_field
! this module defines the data type for a grid of a rank-2 symmetric tensor field
! a 2d grid example:
!     ...
!     type(grid_tensor_field) :: my_grid
!     integer :: nx1, nx2
!     real(8) :: xx1(nx1), xx2(nx2)
!     real(8) :: tarr11(nx1,nx2), tarr12(nx1,nx2), tarr22(nx1,nx2)
!     class(space) :: x
! T = [T11, T12, T22, T13, T23, T33], djTij = [djT1j, djT2j, djT3j]
!     real(8) :: T(6), djTij(3)
!     ...
! initializing the grid:
!     my_grid = grid_tensor_field(nx1, nx2, xx1, xx2, tarr11, tarr12, tarr22) 
! getting field value and divergence at position x:
!     call my_grid%get(x, T, djTij)
! forcibly releasing memory after use
!     call my_grid%clean()

      type, extends(grid) :: grid_tensor_field
         private
         real(8), dimension(:,:,:), allocatable :: tt11, tt12, tt22, &
                                                 & tt13, tt23, tt33
         real(8), dimension(:,:,:), allocatable :: djtt1j, djtt2j, djtt3j
         procedure(grid_tensor_get), pointer, public :: get => null()
      contains
         private
         procedure, public :: set_pointers => grid_tensor_set_pointers
         procedure, public :: data_alloc => grid_tensor_data_alloc
         procedure, public :: data_dealloc => grid_tensor_data_dealloc
         procedure, public :: set_mask => grid_tensor_set_mask
         procedure :: grid_tensor_set_data_1d
         procedure :: grid_tensor_set_data_2d
         procedure :: grid_tensor_set_data_3d
         procedure :: grid_tensor_interp_1d
         procedure :: grid_tensor_interp_2d
         procedure :: grid_tensor_interp_3d
         procedure, public :: clean => grid_tensor_clean
         generic :: set_data => grid_tensor_set_data_1d, grid_tensor_set_data_2d, &
                              & grid_tensor_set_data_3d
         final :: grid_tensor_finalization
      end type grid_tensor_field

      abstract interface
         subroutine grid_tensor_get(this, x, T, djTij)
         import :: grid_tensor_field, space
         class(grid_tensor_field) :: this
         class(space), intent(in) :: x
         real(8), intent(out) :: T(6), djTij(3)
         end subroutine
      end interface

      interface grid_tensor_field
         module procedure grid_tensor_initialization_1d, &
                        & grid_tensor_initialization_2d, &
                        & grid_tensor_initialization_3d
      end interface grid_tensor_field

contains

      function grid_tensor_initialization_1d(nx1, xx1, tarr)
      integer, intent(in) :: nx1
      real(8), intent(in) :: xx1(:)
      real(8), intent(in) :: tarr(:)
      type(grid_tensor_field) :: grid_tensor_initialization_1d

!     allocate the grid and the mask
      call grid_tensor_initialization_1d%grid_alloc(1, nx1)

!     associate all procedure pointers
      call grid_tensor_initialization_1d%set_pointers()

!     allocate data
      call grid_tensor_initialization_1d%data_alloc()

!     set data
      call grid_tensor_initialization_1d%set_data(xx1, tarr)

!     set mask
      call grid_tensor_initialization_1d%set_mask()

      end function grid_tensor_initialization_1d



      function grid_tensor_initialization_2d(nx1, nx2, xx1, xx2, tarr11, tarr12, &
                                           & tarr22)
      integer, intent(in) :: nx1, nx2
      real(8), intent(in) :: xx1(:), xx2(:)
      real(8), intent(in) :: tarr11(:,:), tarr12(:,:), tarr22(:,:)
      type(grid_tensor_field) :: grid_tensor_initialization_2d

!     allocate the grid and the mask
      call grid_tensor_initialization_2d%grid_alloc(2, nx1, nx2)

!     associate all procedure pointers
      call grid_tensor_initialization_2d%set_pointers()

!     allocate data
      call grid_tensor_initialization_2d%data_alloc()

!     set data
      call grid_tensor_initialization_2d%set_data(xx1, xx2, tarr11, tarr12, tarr22)

!     set mask
      call grid_tensor_initialization_2d%set_mask()

      end function grid_tensor_initialization_2d



      function grid_tensor_initialization_3d(nx1, nx2, nx3, xx1, xx2, xx3, &
                                           & tarr11, tarr12, tarr22, tarr13, &
                                           & tarr23, tarr33)
      integer, intent(in) :: nx1, nx2, nx3
      real(8), intent(in) :: xx1(:), xx2(:), xx3(:)
      real(8), intent(in) :: tarr11(:,:,:), tarr12(:,:,:), tarr22(:,:,:)
      real(8), intent(in) :: tarr13(:,:,:), tarr23(:,:,:), tarr33(:,:,:)
      type(grid_tensor_field) :: grid_tensor_initialization_3d

!     allocate the grid and the mask
      call grid_tensor_initialization_3d%grid_alloc(3, nx1, nx2, nx3)

!     associate all procedure pointers
      call grid_tensor_initialization_3d%set_pointers()

!     allocate data
      call grid_tensor_initialization_3d%data_alloc()

!     set data
      call grid_tensor_initialization_3d%set_data(xx1, xx2, xx3, tarr11, tarr12, &
                                                & tarr22, tarr13, tarr23, tarr33)

!     set mask
      call grid_tensor_initialization_3d%set_mask()

      end function grid_tensor_initialization_3d



      subroutine grid_tensor_set_pointers(this)
      class(grid_tensor_field) :: this

      select case( this%ndim )
      case( 1 )
         this%get => grid_tensor_interp_1d
      case( 2 )
         this%get => grid_tensor_interp_2d
      case( 3 )
         this%get => grid_tensor_interp_3d
      case default
         print *,' mod_grid_tensor_field procedure grid_tensor_set_pointers:'
         print *,'  wrong grid_tensor_field%ndim =',this%ndim
         stop
      end select

      end subroutine grid_tensor_set_pointers



      subroutine grid_tensor_data_alloc(this)
      class(grid_tensor_field) :: this
      associate( ndim=>this%ndim, nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3 )
         select case( ndim )
         case( 1 )
            allocate( this%tt11(nx1,nx2,nx3), this%djtt1j(nx1,nx2,nx3) )
            this%tt11 = NoData
            this%djtt1j = 0.d0
         case( 2 )
            allocate( this%tt11(nx1,nx2,nx3), this%tt12(nx1,nx2,nx3), &
                    & this%tt22(nx1,nx2,nx3), this%djtt1j(nx1,nx2,nx3), &
                    & this%djtt2j(nx1,nx2,nx3) )
            this%tt11 = NoData
            this%tt12 = NoData
            this%tt22 = NoData
            this%djtt1j = 0.d0
            this%djtt2j = 0.d0
         case( 3 )
            allocate( this%tt11(nx1,nx2,nx3), this%tt12(nx1,nx2,nx3), &
                    & this%tt22(nx1,nx2,nx3), this%tt13(nx1,nx2,nx3), &
                    & this%tt23(nx1,nx2,nx3), this%tt33(nx1,nx2,nx3), &
                    & this%djtt1j(nx1,nx2,nx3), this%djtt2j(nx1,nx2,nx3), &
                    & this%djtt3j(nx1,nx2,nx3) )
            this%tt11 = NoData
            this%tt12 = NoData
            this%tt22 = NoData
            this%tt13 = NoData
            this%tt23 = NoData
            this%tt33 = NoData
            this%djtt1j = 0.d0
            this%djtt2j = 0.d0
            this%djtt3j = 0.d0
         case default
            print *,' mod_grid_tensor_field procedure grid_tensor_data_alloc:'
            print *,'  wrong grid_tensor_field%ndim =',ndim
            stop
         end select

      end associate
      end subroutine grid_tensor_data_alloc



      subroutine grid_tensor_data_dealloc(this)
      class(grid_tensor_field):: this
      if( allocated(this%tt11) ) deallocate( this%tt11 )
      if( allocated(this%tt12) ) deallocate( this%tt12 )
      if( allocated(this%tt22) ) deallocate( this%tt22 )
      if( allocated(this%tt13) ) deallocate( this%tt13 )
      if( allocated(this%tt23) ) deallocate( this%tt23 )
      if( allocated(this%tt33) ) deallocate( this%tt33 )
      if( allocated(this%djtt1j) ) deallocate( this%djtt1j )
      if( allocated(this%djtt2j) ) deallocate( this%djtt2j )
      if( allocated(this%djtt3j) ) deallocate( this%djtt3j )
      end subroutine grid_tensor_data_dealloc



      subroutine grid_tensor_set_data_1d(this, yy1, tarr)
      class(grid_tensor_field) :: this
      real(8), intent(in) :: yy1(:)
      real(8), intent(in) :: tarr(:)
      associate( nx1=>this%nx1, xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & tt11=>this%tt11, djtt1j=>this%djtt1j )
         if( size(yy1)==nx1 .and. size(tarr)==nx1 ) then
            xx1 = yy1
            tt11(:,1,1) = tarr
         else
            print *,' mod_grid_tensor_field procedure grid_tensor_set_data_1d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(tarr) =',size(tarr)
            stop
         end if

         xx2 = 0.d0
         xx3 = 0.d0

         call this%grid_deriv(tt11, 1, djtt1j, fss002)

      end associate
      end subroutine grid_tensor_set_data_1d



      subroutine grid_tensor_set_data_2d(this, yy1, yy2, tarr11, tarr12, tarr22)
      class(grid_tensor_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:)
      real(8), intent(in) :: tarr11(:,:), tarr12(:,:), tarr22(:,:)
      real(8), dimension(:,:,:), allocatable :: d1tt11, d2tt12, d1tt21, d2tt22
      integer, dimension(2) :: t11shp, t12shp, t22shp
      associate( nx1=>this%nx1, nx2=>this%nx2, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & tt11=>this%tt11, tt12=>this%tt12, tt22=>this%tt22, &
               & djtt1j=>this%djtt1j, djtt2j=>this%djtt2j )

         t11shp = shape(tarr11)
         t12shp = shape(tarr12)
         t22shp = shape(tarr22)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. &
           & t11shp(1)==nx1 .and. t11shp(2)==nx2 .and. &
           & t12shp(1)==nx1 .and. t12shp(2)==nx2 .and. &
           & t22shp(1)==nx1 .and. t22shp(2)==nx2 ) then
            xx1 = yy1
            xx2 = yy2
            tt11(:,:,1) = tarr11
            tt12(:,:,1) = tarr12
            tt22(:,:,1) = tarr22
         else
            print *,' mod_grid_tensor_field procedure grid_tensor_set_data_2d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  shape(tarr11) = (',t11shp(1),t11shp(2),')'
            print *,'  shape(tarr12) = (',t12shp(1),t12shp(2),')'
            print *,'  shape(tarr22) = (',t22shp(1),t22shp(2),')'
            stop
         end if

         xx3 = 0.d0

         allocate( d1tt11(nx1, nx2, 1), d2tt12(nx1, nx2, 1), &
                   d1tt21(nx1, nx2, 1), d2tt22(nx1, nx2, 1) )
         call this%grid_deriv(tt11, 1, d1tt11, fss002)
         call this%grid_deriv(tt12, 2, d2tt12, fss002)
         call this%grid_deriv(tt12, 1, d1tt21, fss002)
         call this%grid_deriv(tt22, 2, d2tt22, fss002)
         djtt1j = d1tt11 + d2tt12
         djtt2j = d1tt21 + d2tt22
         deallocate( d1tt11, d2tt12, d1tt21, d2tt22 )

      end associate
      end subroutine grid_tensor_set_data_2d



      subroutine grid_tensor_set_data_3d(this, yy1, yy2, yy3, tarr11, tarr12, &
                                       & tarr22, tarr13, tarr23, tarr33)
      class(grid_tensor_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:), yy3(:)
      real(8), intent(in) :: tarr11(:,:,:), tarr12(:,:,:), tarr22(:,:,:)
      real(8), intent(in) :: tarr13(:,:,:), tarr23(:,:,:), tarr33(:,:,:)
      real(8), dimension(:,:,:), allocatable :: d1tt11, d2tt12, d3tt13, &
                                              & d1tt21, d2tt22, d3tt23, &
                                              & d1tt31, d2tt32, d3tt33
      integer, dimension(3) :: t11shp, t12shp, t22shp, t13shp, t23shp, t33shp
      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & tt11=>this%tt11, tt12=>this%tt12, tt22=>this%tt22, &
               & tt13=>this%tt13, tt23=>this%tt23, tt33=>this%tt33, &
               & djtt1j=>this%djtt1j, djtt2j=>this%djtt2j, djtt3j=>this%djtt3j )

         t11shp = shape(tarr11)
         t12shp = shape(tarr12)
         t22shp = shape(tarr22)
         t13shp = shape(tarr13)
         t23shp = shape(tarr23)
         t33shp = shape(tarr33)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. size(yy3)==nx3 .and. &
           & t11shp(1)==nx1 .and. t11shp(2)==nx2 .and. t11shp(3)==nx3 .and. &
           & t12shp(1)==nx1 .and. t12shp(2)==nx2 .and. t12shp(3)==nx3 .and. &
           & t22shp(1)==nx1 .and. t22shp(2)==nx2 .and. t22shp(3)==nx3 .and. &
           & t13shp(1)==nx1 .and. t13shp(2)==nx2 .and. t13shp(3)==nx3 .and. &
           & t23shp(1)==nx1 .and. t23shp(2)==nx2 .and. t23shp(3)==nx3 .and. &
           & t33shp(1)==nx1 .and. t33shp(2)==nx2 .and. t33shp(3)==nx3 ) then
            xx1 = yy1
            xx2 = yy2
            xx3 = yy3
            tt11 = tarr11
            tt12 = tarr12
            tt22 = tarr22
            tt13 = tarr13
            tt23 = tarr23
            tt33 = tarr33
         else
            print *,' mod_grid_tensor_field procedure grid_tensor_set_data_3d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2,' nx3 =',nx3
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  size(yy3) =',size(yy3)
            print *,'  shape(tarr11) = (',t11shp(1),t11shp(2),t11shp(3),')'
            print *,'  shape(tarr12) = (',t12shp(1),t12shp(2),t12shp(3),')'
            print *,'  shape(tarr22) = (',t22shp(1),t22shp(2),t22shp(3),')'
            print *,'  shape(tarr13) = (',t13shp(1),t13shp(2),t13shp(3),')'
            print *,'  shape(tarr23) = (',t23shp(1),t23shp(2),t23shp(3),')'
            print *,'  shape(tarr33) = (',t33shp(1),t33shp(2),t33shp(3),')'
            stop
         end if

         allocate( d1tt11(nx1,nx2,nx3), d2tt12(nx1,nx2,nx3), d3tt13(nx1,nx2,nx3), &
                   d1tt21(nx1,nx2,nx3), d2tt22(nx1,nx2,nx3), d3tt23(nx1,nx2,nx3), &
                   d1tt31(nx1,nx2,nx3), d2tt32(nx1,nx2,nx3), d3tt33(nx1,nx2,nx3) )
         call this%grid_deriv(tt11, 1, d1tt11, fss002)
         call this%grid_deriv(tt12, 2, d2tt12, fss002)
         call this%grid_deriv(tt13, 3, d3tt13, fss002)
         call this%grid_deriv(tt12, 1, d1tt21, fss002)
         call this%grid_deriv(tt22, 2, d2tt22, fss002)
         call this%grid_deriv(tt23, 3, d3tt23, fss002)
         call this%grid_deriv(tt13, 1, d1tt31, fss002)
         call this%grid_deriv(tt23, 2, d2tt32, fss002)
         call this%grid_deriv(tt33, 3, d3tt33, fss002)
         djtt1j = d1tt11 + d2tt12 + d3tt13
         djtt2j = d1tt21 + d2tt22 + d3tt23
         djtt3j = d1tt31 + d2tt32 + d3tt33
         deallocate( d1tt11, d2tt12, d3tt13, d1tt21, d2tt22, d3tt23, &
                   & d1tt31, d2tt32, d3tt33 )

      end associate
      end subroutine grid_tensor_set_data_3d



      subroutine grid_tensor_set_mask(this)
      class(grid_tensor_field) :: this
      associate( ndim=>this%ndim, mask=>this%mask, tt11=>this%tt11, &
               & tt12=>this%tt12, tt22=>this%tt22, tt13=>this%tt13, &
               & tt23=>this%tt23, tt33=>this%tt33 )
         mask = 0
         select case( ndim )
         case( 1 )
            where( tt11>NoData ) mask = 1
         case( 2 )
            where( tt11>NoData .and. tt12>NoData .and. tt22>NoData ) mask = 1
         case( 3 )
            where( tt11>NoData .and. tt12>NoData .and. tt22>NoData .and. &
                 & tt13>NoData .and. tt23>NoData .and. tt33>NoData ) mask = 1
         case default
            print *,' mod_grid_tensor_field procedure grid_tensor_set_mask:'
            print *,'  unrecognized ndim =',ndim
            stop
         end select
      end associate
      end subroutine grid_tensor_set_mask



      subroutine grid_tensor_interp_1d(this, x, p, djpij)
      class(grid_tensor_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: p(6)
      real(8), intent(out) :: djpij(3)
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(2), p11, djp1j
      associate( x1=>x%x1, j1=>indices%j1, xx1=>this%xx1, tt11=>this%tt11, &
               & djtt1j=>this%djtt1j )

         call this%search(x, indices, marker)

         if( marker==1 ) then
!           linear interpolation
            y(1) = tt11(j1,1,1)
            y(2) = tt11(j1+1,1,1)
            call linint(y, xx1(j1), xx1(j1+1), x1, p11)

            y(1) = djtt1j(j1,1,1)
            y(2) = djtt1j(j1+1,1,1)
            call linint(y, xx1(j1), xx1(j1+1), x1, djp1j)
         elseif( marker==0 .or. marker==-2 ) then
!           nearest value extrapolation
            p11 = this%avarray(indices, tt11)
            djp1j = this%avarray(indices, djtt1j)
            if( p11==NoData .or. djp1j==NoData ) then
               call grid_tensor_interp_error(x, indices, marker)
            end if
         else
!           the marker==-1 case is an error in 1d because the domain is not allowed
!           to be disconnected. but this case could happen in higher dimensions, see
!           the comments in mod_grid_scalar_field.F90.
            call grid_tensor_interp_error(x, indices, marker)
         end if
         
         if( p11<=0.d0 ) p11 = 0.d0
         p = (/p11, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
         djpij = (/djp1j, 0.d0, 0.d0/)

      end associate
      end subroutine grid_tensor_interp_1d



      subroutine grid_tensor_interp_2d(this, x, p, djpij)
      class(grid_tensor_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: p(6)
      real(8), intent(out) :: djpij(3)
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(4), p11, p12, p22, djp1j, djp2j
      associate( x1=>x%x1, x2=>x%x2, j1=>indices%j1, j2=>indices%j2, &
               & xx1=>this%xx1, xx2=>this%xx2, tt11=>this%tt11, tt12=>this%tt12, &
               & tt22=>this%tt22, djtt1j=>this%djtt1j, djtt2j=>this%djtt2j )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec(y, j1, j2, tt11(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, p11)

            call packvec(y, j1, j2, tt12(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, p12)

            call packvec(y, j1, j2, tt22(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, p22)

            call packvec(y, j1, j2, djtt1j(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, djp1j)

            call packvec(y, j1, j2, djtt2j(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, djp2j)
         elseif( marker==0 .or. marker==-2 ) then
            p11 = this%avarray(indices, tt11)
            p12 = this%avarray(indices, tt12)
            p22 = this%avarray(indices, tt22)
            djp1j = this%avarray(indices, djtt1j)
            djp2j = this%avarray(indices, djtt2j)
            if( p11==NoData .or. p12==NoData .or. p22==NoData .or. &
              & djp1j==NoData .or. djp2j==NoData ) then
               call grid_tensor_interp_error(x, indices, marker)
            end if
         elseif( marker==-1 ) then
            p11 = 0.d0
            p12 = 0.d0
            p22 = 0.d0
            djp1j = 0.d0
            djp2j = 0.d0
         else
            call grid_tensor_interp_error(x, indices, marker)
         end if
         
         if( p11<=0.d0 ) p11 = 0.d0; p12 = 0.d0
         if( p22<=0.d0 ) p12 = 0.d0; p22 = 0.d0
         p = (/p11, p12, p22, 0.d0, 0.d0, 0.d0/)
         djpij = (/djp1j, djp2j, 0.d0/)

      end associate
      end subroutine grid_tensor_interp_2d



      subroutine grid_tensor_interp_3d(this, x, p, djpij)
      class(grid_tensor_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: p(6)
      real(8), intent(out) :: djpij(3)
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(8), p11, p12, p22, p13, p23, p33, djp1j, djp2j, djp3j
      associate( x1=>x%x1, x2=>x%x2, x3=>x%x3, &
               & j1=>indices%j1, j2=>indices%j2, j3=>indices%j3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & tt11=>this%tt11, tt12=>this%tt12, tt22=>this%tt22, &
               & tt13=>this%tt13, tt23=>this%tt23, tt33=>this%tt33, &
               & djtt1j=>this%djtt1j, djtt2j=>this%djtt2j, djtt3j=>this%djtt3j )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec3d(y, j1, j2, j3, tt11)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p11)

            call packvec3d(y, j1, j2, j3, tt12)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p12)

            call packvec3d(y, j1, j2, j3, tt22)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p22)

            call packvec3d(y, j1, j2, j3, tt13)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p13)

            call packvec3d(y, j1, j2, j3, tt23)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p23)

            call packvec3d(y, j1, j2, j3, tt33)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, p33)

            call packvec3d(y, j1, j2, j3, djtt1j)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, djp1j)

            call packvec3d(y, j1, j2, j3, djtt2j)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, djp2j)

            call packvec3d(y, j1, j2, j3, djtt3j)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, djp3j)
         elseif( marker==0 .or. marker==-2 ) then
            p11 = this%avarray(indices, tt11)
            p12 = this%avarray(indices, tt12)
            p22 = this%avarray(indices, tt22)
            p13 = this%avarray(indices, tt13)
            p23 = this%avarray(indices, tt23)
            p33 = this%avarray(indices, tt33)
            djp1j = this%avarray(indices, djtt1j)
            djp2j = this%avarray(indices, djtt2j)
            djp3j = this%avarray(indices, djtt3j)
            if( p11==NoData .or. p12==NoData .or. p22==NoData .or. &
              & p13==NoData .or. p23==NoData .or. p33==NoData .or. &
              & djp1j==NoData .or. djp2j==NoData .or. djp3j==NoData ) then
               call grid_tensor_interp_error(x, indices, marker)
            end if
         elseif( marker==-1 ) then
            p11 = 0.d0
            p12 = 0.d0
            p22 = 0.d0
            p13 = 0.d0
            p23 = 0.d0
            p33 = 0.d0
            djp1j = 0.d0
            djp2j = 0.d0
            djp3j = 0.d0
         else
            call grid_tensor_interp_error(x, indices, marker)
         end if
         
         if( p11<=0.d0 ) p11 = 0.d0; p12 = 0.d0; p13 = 0.d0
         if( p22<=0.d0 ) p12 = 0.d0; p22 = 0.d0; p23 = 0.d0
         if( p33<=0.d0 ) p13 = 0.d0; p23 = 0.d0; p33 = 0.d0
         p = (/p11, p12, p22, p13, p23, p33/)
         djpij = (/djp1j, djp2j, djp3j/)

      end associate
      end subroutine grid_tensor_interp_3d



      subroutine grid_tensor_interp_error(x, indices, marker)
      class(space), intent(in) :: x
      type(grid_indices), intent(in) :: indices
      integer, intent(in) :: marker
      print *,' mod_grid_tensor_field procedure grid_tensor_interp_error:'
      print *,'  the interpolating points appears out of domain'
      print *,'  x = (',x%x1,',',x%x2,',',x%x3,')'
      print *,'  indices = (',indices%j1,',',indices%j2,',',indices%j3,')'
      print *,'  marker =',marker
      stop
      end subroutine grid_tensor_interp_error



      subroutine grid_tensor_clean(this)
      class(grid_tensor_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_tensor_clean



      subroutine grid_tensor_finalization(this)
      type(grid_tensor_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_tensor_finalization

end module mod_grid_tensor_field

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 11/18/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


