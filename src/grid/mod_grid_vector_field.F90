module mod_grid_vector_field
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
      public grid_vector_field
! this module defines the data type for a grid of a vector field
! a 2d grid example:
!     ...
!     type(grid_vector_field) :: my_grid
!     integer :: nx1, nx2
!     real(8) :: xx1(nx1), xx2(nx2)
!     real(8) :: varr1(nx1,nx2), varr2(nx1,nx2)
!     class(space) :: x
!     real(8) :: v(3), divi  ! v = [v1, v2, 0], divi = d1v1 + d2v2
!     ...
! initializing the grid:
!     my_grid = grid_vector_field(nx1, nx2, xx1, xx2, varr1, varr2) 
! getting field value and divergence at position x:
!     call my_grid%get(x, v, divi)
! forcibly releasing memory after use
!     call my_grid%clean()

      type, extends(grid) :: grid_vector_field
         private
         real(8), dimension(:,:,:), allocatable :: vv1, vv2, vv3, divi
         procedure(grid_vector_get), pointer, public :: get => null()
      contains
         private
         procedure, public :: set_pointers => grid_vector_set_pointers
         procedure, public :: data_alloc => grid_vector_data_alloc
         procedure, public :: data_dealloc => grid_vector_data_dealloc
         procedure, public :: set_mask => grid_vector_set_mask
         procedure :: grid_vector_set_data_1d
         procedure :: grid_vector_set_data_2d
         procedure :: grid_vector_set_data_3d
         procedure :: grid_vector_interp_1d
         procedure :: grid_vector_interp_2d
         procedure :: grid_vector_interp_3d
         procedure, public :: clean => grid_vector_clean
         generic :: set_data => grid_vector_set_data_1d, grid_vector_set_data_2d, &
                              & grid_vector_set_data_3d
         final :: grid_vector_finalization
      end type grid_vector_field

      abstract interface
         subroutine grid_vector_get(this, x, v, divi)
         import :: grid_vector_field, space
         class(grid_vector_field) :: this
         class(space), intent(in) :: x
         real(8), intent(out) :: v(3), divi
         end subroutine
      end interface

      interface grid_vector_field
         module procedure grid_vector_initialization_1d, &
                        & grid_vector_initialization_2d, &
                        & grid_vector_initialization_3d
      end interface grid_vector_field

contains

      function grid_vector_initialization_1d(nx1, xx1, varr)
      integer, intent(in) :: nx1
      real(8), intent(in) :: xx1(:)
      real(8), intent(in) :: varr(:)
      type(grid_vector_field) :: grid_vector_initialization_1d

!     allocate the grid and the mask
      call grid_vector_initialization_1d%grid_alloc(1, nx1)

!     associate all procedure pointers
      call grid_vector_initialization_1d%set_pointers()

!     allocate data
      call grid_vector_initialization_1d%data_alloc()

!     set data
      call grid_vector_initialization_1d%set_data(xx1, varr)

!     set mask
      call grid_vector_initialization_1d%set_mask()

      end function grid_vector_initialization_1d



      function grid_vector_initialization_2d(nx1, nx2, xx1, xx2, varr1, varr2)
      integer, intent(in) :: nx1, nx2
      real(8), intent(in) :: xx1(:), xx2(:)
      real(8), intent(in) :: varr1(:,:), varr2(:,:)
      type(grid_vector_field) :: grid_vector_initialization_2d

!     allocate the grid and the mask
      call grid_vector_initialization_2d%grid_alloc(2, nx1, nx2)

!     associate all procedure pointers
      call grid_vector_initialization_2d%set_pointers()

!     allocate data
      call grid_vector_initialization_2d%data_alloc()

!     set data
      call grid_vector_initialization_2d%set_data(xx1, xx2, varr1, varr2)

!     set mask
      call grid_vector_initialization_2d%set_mask()

      end function grid_vector_initialization_2d



      function grid_vector_initialization_3d(nx1, nx2, nx3, xx1, xx2, xx3, &
                                           & varr1, varr2, varr3)
      integer, intent(in) :: nx1, nx2, nx3
      real(8), intent(in) :: xx1(:), xx2(:), xx3(:)
      real(8), intent(in) :: varr1(:,:,:), varr2(:,:,:), varr3(:,:,:)
      type(grid_vector_field) :: grid_vector_initialization_3d

!     allocate the grid and the mask
      call grid_vector_initialization_3d%grid_alloc(3, nx1, nx2, nx3)

!     associate all procedure pointers
      call grid_vector_initialization_3d%set_pointers()

!     allocate data
      call grid_vector_initialization_3d%data_alloc()

!     set data
      call grid_vector_initialization_3d%set_data(xx1, xx2, xx3, &
                                                & varr1, varr2, varr3)

!     set mask
      call grid_vector_initialization_3d%set_mask()

      end function grid_vector_initialization_3d



      subroutine grid_vector_set_pointers(this)
      class(grid_vector_field) :: this

      select case( this%ndim )
      case( 1 )
         this%get => grid_vector_interp_1d
      case( 2 )
         this%get => grid_vector_interp_2d
      case( 3 )
         this%get => grid_vector_interp_3d
      case default
         print *,' mod_grid_vector_field procedure grid_vector_set_pointers:'
         print *,'  wrong grid_vector_field%ndim =',this%ndim
         stop
      end select

      end subroutine grid_vector_set_pointers



      subroutine grid_vector_data_alloc(this)
      class(grid_vector_field) :: this
      associate( ndim=>this%ndim, nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3 )
         select case( ndim )
         case( 1 )
            allocate( this%vv1(nx1, nx2, nx3) )
            this%vv1 = NoData
         case( 2 )
            allocate( this%vv1(nx1, nx2, nx3), this%vv2(nx1, nx2, nx3) )
            this%vv1 = NoData
            this%vv2 = NoData
         case( 3 )
            allocate( this%vv1(nx1, nx2, nx3), this%vv2(nx1, nx2, nx3), &
                    & this%vv3(nx1, nx2, nx3) )
            this%vv1 = NoData
            this%vv2 = NoData
            this%vv3 = NoData
         case default
            print *,' mod_grid_vector_field procedure grid_vector_data_alloc:'
            print *,'  wrong grid_vector_field%ndim =',ndim
            stop
         end select

         allocate( this%divi(nx1, nx2, nx3) )
         this%divi = 0.d0

      end associate
      end subroutine grid_vector_data_alloc



      subroutine grid_vector_data_dealloc(this)
      class(grid_vector_field):: this
      if( allocated(this%vv1) ) deallocate( this%vv1 )
      if( allocated(this%vv2) ) deallocate( this%vv2 )
      if( allocated(this%vv3) ) deallocate( this%vv3 )
      if( allocated(this%divi) ) deallocate( this%divi )
      end subroutine grid_vector_data_dealloc



      subroutine grid_vector_set_data_1d(this, yy1, varr)
      class(grid_vector_field) :: this
      real(8), intent(in) :: yy1(:)
      real(8), intent(in) :: varr(:)
      associate( nx1=>this%nx1, xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & vv1=>this%vv1, divi=>this%divi )
         if( size(yy1)==nx1 .and. size(varr)==nx1 ) then
            xx1 = yy1
            vv1(:,1,1) = varr
         else
            print *,' mod_grid_vector_field procedure grid_vector_set_data_1d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(varr) =',size(varr)
            stop
         end if

         xx2 = 0.d0
         xx3 = 0.d0

         call this%grid_deriv(vv1, 1, divi, fss002)

      end associate
      end subroutine grid_vector_set_data_1d



      subroutine grid_vector_set_data_2d(this, yy1, yy2, varr1, varr2)
      class(grid_vector_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:)
      real(8), intent(in) :: varr1(:,:), varr2(:,:)
      real(8), dimension(:,:,:), allocatable :: d1vv1, d2vv2
      integer, dimension(2) :: v1shp, v2shp
      associate( nx1=>this%nx1, nx2=>this%nx2, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & vv1=>this%vv1, vv2=>this%vv2, divi=>this%divi )

         v1shp = shape(varr1)
         v2shp = shape(varr2)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. &
           & v1shp(1)==nx1 .and. v1shp(2)==nx2 .and. &
           & v2shp(1)==nx1 .and. v2shp(2)==nx2 ) then
            xx1 = yy1
            xx2 = yy2
            vv1(:,:,1) = varr1
            vv2(:,:,1) = varr2
         else
            print *,' mod_grid_vector_field procedure grid_vector_set_data_2d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  shape(varr1) = (',v1shp(1),v1shp(2),')'
            print *,'  shape(varr2) = (',v2shp(1),v2shp(2),')'
            stop
         end if

         xx3 = 0.d0

         allocate( d1vv1(nx1, nx2, 1), d2vv2(nx1, nx2, 1) )
         call this%grid_deriv(vv1, 1, d1vv1, fss002)
         call this%grid_deriv(vv2, 2, d2vv2, fss002)
         divi = d1vv1 + d2vv2
         deallocate( d1vv1, d2vv2 )

      end associate
      end subroutine grid_vector_set_data_2d



      subroutine grid_vector_set_data_3d(this, yy1, yy2, yy3, varr1, varr2, varr3)
      class(grid_vector_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:), yy3(:)
      real(8), intent(in) :: varr1(:,:,:), varr2(:,:,:), varr3(:,:,:)
      real(8), dimension(:,:,:), allocatable :: d1vv1, d2vv2, d3vv3
      integer, dimension(3) :: v1shp, v2shp, v3shp
      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & vv1=>this%vv1, vv2=>this%vv2, vv3=>this%vv3, divi=>this%divi )

         v1shp = shape(varr1)
         v2shp = shape(varr2)
         v3shp = shape(varr3)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. size(yy3)==nx3 .and. &
           & v1shp(1)==nx1 .and. v1shp(2)==nx2 .and. v1shp(3)==nx3 .and. &
           & v2shp(1)==nx1 .and. v2shp(2)==nx2 .and. v2shp(3)==nx3 .and. &
           & v3shp(1)==nx1 .and. v3shp(2)==nx2 .and. v3shp(3)==nx3 ) then
            xx1 = yy1
            xx2 = yy2
            xx3 = yy3
            vv1 = varr1
            vv2 = varr2
            vv3 = varr3
         else
            print *,' mod_grid_vector_field procedure grid_vector_set_data_3d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2,' nx3 =',nx3
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  size(yy3) =',size(yy3)
            print *,'  shpae(varr1) = (',v1shp(1),v1shp(2),v1shp(3),')'
            print *,'  shpae(varr2) = (',v2shp(1),v2shp(2),v2shp(3),')'
            print *,'  shpae(varr3) = (',v3shp(1),v3shp(2),v3shp(3),')'
            stop
         end if

         allocate( d1vv1(nx1,nx2,nx3), d2vv2(nx1,nx2,nx3), d3vv3(nx1,nx2,nx3) )
         call this%grid_deriv(vv1, 1, d1vv1, fss002)
         call this%grid_deriv(vv2, 2, d2vv2, fss002)
         call this%grid_deriv(vv3, 3, d3vv3, fss002)
         divi = d1vv1 + d2vv2 + d3vv3
         deallocate( d1vv1, d2vv2, d3vv3 )

      end associate
      end subroutine grid_vector_set_data_3d



      subroutine grid_vector_set_mask(this)
      class(grid_vector_field) :: this
      associate( ndim=>this%ndim, mask=>this%mask, &
               & vv1=>this%vv1, vv2=>this%vv2, vv3=>this%vv3 )
         mask = 0
         select case( ndim )
         case( 1 )
            where( vv1>NoData ) mask = 1
         case( 2 )
            where( vv1>NoData .and. vv2>NoData ) mask = 1
         case( 3 )
            where( vv1>NoData .and. vv2>NoData .and. vv3>NoData ) mask = 1
         case default
            print *,' mod_grid_vector_field procedure grid_vector_set_mask:'
            print *,'  unrecognized ndim =',ndim
            stop
         end select
      end associate
      end subroutine grid_vector_set_mask



      subroutine grid_vector_interp_1d(this, x, w, diwi)
      class(grid_vector_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: w(3)
      real(8), intent(out) :: diwi
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(2), w1
      associate( x1=>x%x1, j1=>indices%j1, xx1=>this%xx1, vv1=>this%vv1, &
               & divi=>this%divi )

         call this%search(x, indices, marker)

         if( marker==1 ) then
!           linear interpolation
            y(1) = vv1(j1,1,1)
            y(2) = vv1(j1+1,1,1)
            call linint(y, xx1(j1), xx1(j1+1), x1, w1)

            y(1) = divi(j1,1,1)
            y(2) = divi(j1+1,1,1)
            call linint(y, xx1(j1), xx1(j1+1), x1, diwi)
         elseif( marker==0 .or. marker==-2 ) then
!           nearest value extrapolation
            w1 = this%avarray(indices, vv1)
            diwi = this%avarray(indices, divi)
            if( w1==NoData .or. diwi==NoData ) then
               call grid_vector_interp_error(x, indices, marker)
            end if
         else
!           the marker==-1 case is an error in 1d because the domain is not allowed
!           to be disconnected. but this case could happen in higher dimensions, see
!           the comments in mod_grid_scalar_field.F90.
            call grid_vector_interp_error(x, indices, marker)
         end if
         
         w = (/w1, 0.d0, 0.d0/)

      end associate
      end subroutine grid_vector_interp_1d



      subroutine grid_vector_interp_2d(this, x, w, diwi)
      class(grid_vector_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: w(3)
      real(8), intent(out) :: diwi
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(4), w1, w2
      associate( x1=>x%x1, x2=>x%x2, j1=>indices%j1, j2=>indices%j2, &
               & xx1=>this%xx1, xx2=>this%xx2, &
               & vv1=>this%vv1, vv2=>this%vv2, divi=>this%divi )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec(y, j1, j2, vv1(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, w1)

            call packvec(y, j1, j2, vv2(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, w2)

            call packvec(y, j1, j2, divi(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, diwi)
         elseif( marker==0 .or. marker==-2 ) then
            w1 = this%avarray(indices, vv1)
            w2 = this%avarray(indices, vv2)
            diwi = this%avarray(indices, divi)
            if( w1==NoData .or. w2==NoData .or. diwi==NoData ) then
               call grid_vector_interp_error(x, indices, marker)
            end if
         elseif( marker==-1 ) then
            w1 = 0.d0
            w2 = 0.d0
            diwi = 0.d0
         else
            call grid_vector_interp_error(x, indices, marker)
         end if
         
         w = (/w1, w2, 0.d0/)

      end associate
      end subroutine grid_vector_interp_2d



      subroutine grid_vector_interp_3d(this, x, w, diwi)
      class(grid_vector_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: w(3)
      real(8), intent(out) :: diwi
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(8), w1, w2, w3
      associate( x1=>x%x1, x2=>x%x2, x3=>x%x3, &
               & j1=>indices%j1, j2=>indices%j2, j3=>indices%j3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & vv1=>this%vv1, vv2=>this%vv2, vv3=>this%vv3, divi=>this%divi )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec3d(y, j1, j2, j3, vv1)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, w1)

            call packvec3d(y, j1, j2, j3, vv2)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, w2)

            call packvec3d(y, j1, j2, j3, vv3)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, w3)

            call packvec3d(y, j1, j2, j3, divi)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                          xx3(j3), xx3(j3+1), x1, x2, x3, diwi)
         elseif( marker==0 .or. marker==-2 ) then
            w1 = this%avarray(indices, vv1)
            w2 = this%avarray(indices, vv2)
            w3 = this%avarray(indices, vv3)
            diwi = this%avarray(indices, divi)
            if( w1==NoData .or. w2==NoData .or. w3==NoData .or. diwi==NoData ) then
               call grid_vector_interp_error(x, indices, marker)
            end if
         elseif( marker==-1 ) then
            w1 = 0.d0
            w2 = 0.d0
            w3 = 0.d0
            diwi = 0.d0
         else
            call grid_vector_interp_error(x, indices, marker)
         end if
         
         w = (/w1, w2, w3/)

      end associate
      end subroutine grid_vector_interp_3d



      subroutine grid_vector_interp_error(x, indices, marker)
      class(space), intent(in) :: x
      type(grid_indices), intent(in) :: indices
      integer, intent(in) :: marker
      print *,' mod_grid_vector_field procedure grid_vector_interp_error:'
      print *,'  the interpolating points appears out of domain'
      print *,'  x = (',x%x1,',',x%x2,',',x%x3,')'
      print *,'  indices = (',indices%j1,',',indices%j2,',',indices%j3,')'
      print *,'  marker =',marker
      stop
      end subroutine grid_vector_interp_error



      subroutine grid_vector_clean(this)
      class(grid_vector_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_vector_clean



      subroutine grid_vector_finalization(this)
      type(grid_vector_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_vector_finalization

end module mod_grid_vector_field

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 11/17/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


