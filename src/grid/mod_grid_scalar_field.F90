module mod_grid_scalar_field
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 11/16/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, only: space, NoData
      use mod_grid_typedef
      use interpolation
      implicit none
      save
      private
      public grid_scalar_field
! this module defines the data type for a grid of a scalar field
! a 2d grid example:
!     ...
!     type(grid_scalar_field) :: my_grid
!     integer :: nx1, nx2
!     real(8) :: xx1(nx1), xx2(nx2)
!     real(8) :: darr(nx1,nx2)
!     class(space) :: x
!     real(8) :: val
!     ...
! initializing the grid for a positive semi-definite field:
!     my_grid = grid_scalar_field(nx1, nx2, xx1, xx2, darr, 'P') 
! initializing the grid for an indefinite field:
!     my_grid = grid_scalar_field(nx1, nx2, xx1, xx2, darr) 
! or put any character other than 'P', e.g.,
!     my_grid = grid_scalar_field(nx1, nx2, xx1, xx2, darr, 'U') 
! getting field value at position x:
!     call my_grid%get(x, val)
! forcibly releasing memory after use
!     call my_grid%clean()

      type, extends(grid) :: grid_scalar_field
         real(8), dimension(:,:,:), allocatable :: aa
         procedure(grid_scalar_get), pointer :: get => null()
      contains
         private
         procedure, public :: set_pointers => grid_scalar_set_pointers
         procedure, public :: set_pointers_p => grid_scalar_set_pointers_p
         procedure, public :: data_alloc => grid_scalar_data_alloc
         procedure, public :: data_dealloc => grid_scalar_data_dealloc
         procedure, public :: set_mask => grid_scalar_set_mask
         procedure :: grid_scalar_set_data_1d
         procedure :: grid_scalar_set_data_2d
         procedure :: grid_scalar_set_data_3d
         procedure :: grid_scalar_interp_1d
         procedure :: grid_scalar_interp_2d
         procedure :: grid_scalar_interp_3d
         procedure :: grid_scalar_interp_1d_p
         procedure :: grid_scalar_interp_2d_p
         procedure :: grid_scalar_interp_3d_p
         procedure, public :: clean => grid_scalar_clean
         generic :: set_data => grid_scalar_set_data_1d, grid_scalar_set_data_2d, &
                              & grid_scalar_set_data_3d
         final :: grid_scalar_finalization
      end type grid_scalar_field

      abstract interface
         subroutine grid_scalar_get(this, x, val)
         import :: grid_scalar_field, space
         class(grid_scalar_field) :: this
         class(space), intent(in) :: x
         real(8), intent(out) :: val
         end subroutine
      end interface

      interface grid_scalar_field
         module procedure grid_scalar_initialization_1d, &
                        & grid_scalar_initialization_2d, &
                        & grid_scalar_initialization_3d
      end interface grid_scalar_field

contains

      function grid_scalar_initialization_1d(nx1, xx1, darr, psd)
      integer, intent(in) :: nx1
      real(8), intent(in) :: xx1(:)
      real(8), intent(in) :: darr(:)
!     specify psd=='P' if the scalar field is meant to be positive semi-definite.
!     this affects how the scalar field is interpolated.
      character(1), intent(in), optional :: psd
      type(grid_scalar_field) :: grid_scalar_initialization_1d

!     allocate the grid and the mask
      call grid_scalar_initialization_1d%grid_alloc(1, nx1)

!     associate all procedure pointers
      if( present(psd) .and. (psd=='P' .or. psd=='p') ) then
         call grid_scalar_initialization_1d%set_pointers_p()
      else
         call grid_scalar_initialization_1d%set_pointers()
      end if

!     allocate data
      call grid_scalar_initialization_1d%data_alloc()

!     set data
      call grid_scalar_initialization_1d%set_data(xx1, darr)

!     set mask
      call grid_scalar_initialization_1d%set_mask()

      end function grid_scalar_initialization_1d



      function grid_scalar_initialization_2d(nx1, nx2, xx1, xx2, darr, psd)
      integer, intent(in) :: nx1, nx2
      real(8), intent(in) :: xx1(:), xx2(:)
      real(8), intent(in) :: darr(:,:)
      character(1), intent(in), optional :: psd
      type(grid_scalar_field) :: grid_scalar_initialization_2d

!     allocate the grid and the mask
      call grid_scalar_initialization_2d%grid_alloc(2, nx1, nx2)

!     associate all procedure pointers
      if( present(psd) .and. (psd=='P' .or. psd=='p') ) then
         call grid_scalar_initialization_2d%set_pointers_p()
      else
         call grid_scalar_initialization_2d%set_pointers()
      end if

!     allocate data
      call grid_scalar_initialization_2d%data_alloc()

!     set data
      call grid_scalar_initialization_2d%set_data(xx1, xx2, darr)

!     set mask
      call grid_scalar_initialization_2d%set_mask()

      end function grid_scalar_initialization_2d



      function grid_scalar_initialization_3d(nx1, nx2, nx3, xx1, xx2, xx3, darr, &
                                           & psd)
      integer, intent(in) :: nx1, nx2, nx3
      real(8), intent(in) :: xx1(:), xx2(:), xx3(:)
      real(8), intent(in) :: darr(:,:,:)
      character(1), intent(in), optional :: psd
      type(grid_scalar_field) :: grid_scalar_initialization_3d

!     allocate the grid and the mask
      call grid_scalar_initialization_3d%grid_alloc(3, nx1, nx2, nx3)

!     associate all procedure pointers
      if( present(psd) .and. (psd=='P' .or. psd=='p') ) then
         call grid_scalar_initialization_3d%set_pointers_p()
      else
         call grid_scalar_initialization_3d%set_pointers()
      end if

!     allocate data
      call grid_scalar_initialization_3d%data_alloc()

!     set data
      call grid_scalar_initialization_3d%set_data(xx1, xx2, xx3, darr)

!     set mask
      call grid_scalar_initialization_3d%set_mask()

      end function grid_scalar_initialization_3d



      subroutine grid_scalar_set_pointers(this)
      class(grid_scalar_field) :: this

      select case( this%ndim )
      case( 1 )
         this%get => grid_scalar_interp_1d
      case( 2 )
         this%get => grid_scalar_interp_2d
      case( 3 )
         this%get => grid_scalar_interp_3d
      case default
         print *,' mod_grid_scalar_field procedure grid_scalar_set_pointers:'
         print *,'  wrong grid_scalar_field%ndim =',this%ndim
         stop
      end select

      end subroutine grid_scalar_set_pointers



      subroutine grid_scalar_set_pointers_p(this)
      class(grid_scalar_field) :: this

      select case( this%ndim )
      case( 1 )
         this%get => grid_scalar_interp_1d_p
      case( 2 )
         this%get => grid_scalar_interp_2d_p
      case( 3 )
         this%get => grid_scalar_interp_3d_p
      case default
         print *,' mod_grid_scalar_field procedure grid_scalar_set_pointers_p:'
         print *,'  wrong grid_scalar_field%ndim =',this%ndim
         stop
      end select

      end subroutine grid_scalar_set_pointers_p



      subroutine grid_scalar_data_alloc(this)
      class(grid_scalar_field) :: this
      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3 )
         allocate( this%aa(nx1, nx2, nx3) )
         this%aa = NoData
      end associate
      end subroutine grid_scalar_data_alloc



      subroutine grid_scalar_data_dealloc(this)
      class(grid_scalar_field) :: this
      if( allocated(this%aa) ) deallocate( this%aa )
      end subroutine grid_scalar_data_dealloc



      subroutine grid_scalar_set_data_1d(this, yy1, darr)
      class(grid_scalar_field) :: this
      real(8), intent(in) :: yy1(:)
      real(8), intent(in) :: darr(:)
      associate( nx1=>this%nx1, xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & aa=>this%aa )
         if( size(yy1)==nx1 .and. size(darr)==nx1 ) then
            xx1 = yy1
            aa(:,1,1) = darr
         else
            print *,' mod_grid_scalar_field procedure grid_scalar_set_data_1d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(darr) =',size(darr)
            stop
         end if

         xx2 = 0.d0
         xx3 = 0.d0
      end associate
      end subroutine grid_scalar_set_data_1d



      subroutine grid_scalar_set_data_2d(this, yy1, yy2, darr)
      class(grid_scalar_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:)
      real(8), intent(in) :: darr(:,:)
      integer, dimension(2) :: dshp
      associate( nx1=>this%nx1, nx2=>this%nx2, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, aa=>this%aa )

         dshp = shape(darr)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. &
           & dshp(1)==nx1 .and. dshp(2)==nx2 ) then
            xx1 = yy1
            xx2 = yy2
            aa(:,:,1) = darr
         else 
            print *,' mod_grid_scalar_field procedure grid_scalar_set_data_2d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  size(darr,1) =',dshp(1)
            print *,'  size(darr,2) =',dshp(2)
            stop
         end if

         xx3 = 0.d0

      end associate
      end subroutine grid_scalar_set_data_2d



      subroutine grid_scalar_set_data_3d(this, yy1, yy2, yy3, darr)
      class(grid_scalar_field) :: this
      real(8), intent(in) :: yy1(:), yy2(:), yy3(:)
      real(8), intent(in) :: darr(:,:,:)
      integer, dimension(3) :: dshp
      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, aa=>this%aa )

         dshp = shape(darr)

         if( size(yy1)==nx1 .and. size(yy2)==nx2 .and. size(yy3)==nx3 .and. &
           & dshp(1)==nx1 .and. dshp(2)==nx2 .and. dshp(3)==nx3 ) then
            xx1 = yy1
            xx2 = yy2
            xx3 = yy3
            aa = darr
         else 
            print *,' mod_grid_scalar_field procedure grid_scalar_set_data_3d:'
            print *,'  size mismatch in input arrays:'
            print *,'  nx1 =',nx1,' nx2 =',nx2,' nx3 =',nx3
            print *,'  size(yy1) =',size(yy1)
            print *,'  size(yy2) =',size(yy2)
            print *,'  size(yy3) =',size(yy3)
            print *,'  size(darr,1) =',dshp(1)
            print *,'  size(darr,2) =',dshp(2)
            print *,'  size(darr,3) =',dshp(3)
            stop
         end if

      end associate
      end subroutine grid_scalar_set_data_3d



      subroutine grid_scalar_set_mask(this)
      class(grid_scalar_field) :: this
      associate( mask=>this%mask, aa=>this%aa )
         mask = 0
         where( aa>NoData ) mask = 1
      end associate
      end subroutine grid_scalar_set_mask



      subroutine count_np_points(this, j1, np, j0, n)
!     this routine finds the n (in [1, np]) consecutive valued grid nodes either to
!     the left or to the right of the given index j1 (inclusive if the j1 node
!     itself is valued). the direction is determined by the end of the valued grid
!     on which j1 is. the leftmost index of these n points is returned as j0. this
!     routine is only intended to work for 1d grid, and it presumes that one and
!     only one of the {j1, j1+1} nodes is valued. 
      class(grid_scalar_field) :: this
      integer, intent(in) :: j1, np
      integer, intent(out) :: j0, n

      if( np<=0 ) then
         print *,' mod_grid_scalar_field procedure count_np_points:'
         print *,'  input integer np (=',np,') must be positive'
         stop
      end if

      associate( nx1=>this%nx1, mask=>this%mask )

         if( j1<nx1 ) then
            if( mask(j1+1,1,1)==1 ) then
!              find np points to the right of j1
               j0 = j1 + 1
               call count_to_the_right(j0, np, n)
            elseif( mask(j1,1,1)==1 ) then
!              find np points to the left of j1 and update j0
               j0 = j1
               call count_to_the_left(j0, np, n)
            else
               call count_np_points_error()
            end if
         elseif( j1==nx1 ) then
            if( mask(j1,1,1)==1 ) then
!              find np points to the left of j1 and update j0
               j0 = j1
               call count_to_the_left(j0, np, n)
            else
               call count_np_points_error()
            end if
         else
            call count_np_points_error()
         end if

      end associate
      contains

         subroutine count_to_the_right(j0, np, n)
!        count np points with mask = 1 starting from j0 to its right, the counting
!        result is n.
         integer :: j0, np, n
         integer :: i
         associate( nx1=>this%nx1, mask=>this%mask )
            n = 0
            do i=j0,j0+np-1
               if( i<=nx1 .and. mask(i,1,1)==1 ) then
                  n = n + 1
               else
                  exit
               end if
            end do
         end associate
         end subroutine count_to_the_right


         subroutine count_to_the_left(j0, np, n)
!        count np points with mask = 1 starting from j0 to its left, the counting
!        result is n. j0 is then updated to the index of the leftmost point.
         integer :: j0, np, n
         integer :: i
         associate( mask=>this%mask )
            n = 0
            do i=j0,j0-np+1,-1
               if( i>=1 .and. mask(i,1,1)==1 ) then
                  n = n + 1
               else
                  exit
               end if
            end do
            j0 = j0 - n + 1
         end associate
         end subroutine count_to_the_left


         subroutine count_np_points_error()
         print *,' mod_grid_scalar_field procedure count_np_points:'
         print *,'  wrong input index j1 =',j1
         call this%print_mask()
         stop
         end subroutine count_np_points_error

      end subroutine count_np_points



      subroutine grid_scalar_interp_1d(this, x, val)
!     use linear interpolation and nearest value extrapolation for a regular scalar
!     field (most probably used as an equation coefficient)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(2)
      associate( x1=>x%x1, j1=>indices%j1, xx1=>this%xx1, aa=>this%aa )
         
         call this%search(x, indices, marker)

         if( marker==1 ) then
            y(1) = aa(j1,1,1)
            y(2) = aa(j1+1,1,1)
            call linint(y, xx1(j1), xx1(j1+1), x1, val)
         elseif( marker==0 .or. marker==-2 ) then
!           nearest value extrapolation
            val = this%avarray(indices, aa)
            if( val==NoData ) call grid_scalar_interp_error(x, indices, marker)
         else
!           the marker==-1 case is an error in 1d because the domain is not allowed
!           to be disconnected. but this case could happen in higher dimensions, see
!           the comments below.
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_1d



      subroutine grid_scalar_interp_1d_p(this, x, val)
!     use linear interpolation and linear extrapolation in log-scale for 1d positive
!     semi-definite scalar field (most probably used as a solution)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker, j0, n
      real(8) :: y(2)
      associate( x1=>x%x1, j1=>indices%j1, xx1=>this%xx1, aa=>this%aa )
         
         call this%search(x, indices, marker)

         if( marker==1 ) then
            y(1) = aa(j1,1,1)
            y(2) = aa(j1+1,1,1)
            call loglinint(y, xx1(j1), xx1(j1+1), x1, val)
         elseif( marker==0 .or. marker==-2 ) then
!           linear extrapolation
            call count_np_points(this, j1, 2, j0, n)
            if( n==2 ) then
               y(1) = aa(j0,1,1)
               y(2) = aa(j0+1,1,1)
               call loglinint(y, xx1(j0), xx1(j0+1), x1, val)
            elseif( n==1 ) then
               val = aa(j0,1,1)
            else
               call grid_scalar_interp_error(x, indices, marker)
            end if
         else
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_1d_p



      subroutine grid_scalar_interp_2d(this, x, val)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(4)

      associate( x1=>x%x1, x2=>x%x2, j1=>indices%j1, j2=>indices%j2, &
               & nx1=>this%nx1, nx2=>this%nx2, xx1=>this%xx1, xx2=>this%xx2, &
               & aa=>this%aa )
         
         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec(y, j1, j2, aa(:,:,1))
            call blnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, val)
         elseif( marker==0 .or. marker==-2 ) then
!           nearest value extrapolation
            val = this%avarray(indices, aa)
!           the following case occurs when x is out of the grid bounds and its
!           nearest grid nodes all have no value, hence there is no way to do
!           extrapolation. this has two possible causes: 1) some grid nodes within
!           the domain are not evaluated; 2) x is out of the domain. either of these
!           causes is an error.
            if( val==NoData ) call grid_scalar_interp_error(x, indices, marker)
         elseif( marker==-1 ) then
!           given that x remains in domain, this case could happen if the local
!           domain is of a finger-shape and the width of the finger is narrower than
!           the size of a grid cell. if the field is used to provide equation
!           coefficients, a zero coefficient would prevent the stochastic process
!           within this pathological domain region from further developing.
            val = 0.d0
         else
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_2d



      subroutine grid_scalar_interp_2d_p(this, x, val)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(4)

      associate( x1=>x%x1, x2=>x%x2, j1=>indices%j1, j2=>indices%j2, &
               & nx1=>this%nx1, nx2=>this%nx2, xx1=>this%xx1, xx2=>this%xx2, &
               & aa=>this%aa )
         
         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec(y, j1, j2, aa(:,:,1))
            call logblnint(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), x1, x2, val)
         elseif( marker==0 .or. marker==-2 ) then
            val = this%avarray(indices, aa)
            if( val==NoData ) then
               call grid_scalar_interp_error(x, indices, marker)
            elseif( val<0.d0 ) then
               val = 0.d0
            end if
         elseif( marker==-1 ) then
            val = 0.d0
         else
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_2d_p



      subroutine grid_scalar_interp_3d(this, x, val)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(8)

      associate( x1=>x%x1, x2=>x%x2, x3=>x%x3, &
               & j1=>indices%j1, j2=>indices%j2, j3=>indices%j3, &
               & nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & aa=>this%aa )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec3d(y, j1, j2, j3, aa)
            call blnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                        & xx3(j3), xx3(j3+1), x1, x2, x3, val)
         elseif( marker==0 .or. marker==-2 ) then
            val = this%avarray(indices, aa)
            if( val==NoData ) call grid_scalar_interp_error(x, indices, marker)
         elseif( marker==-1 ) then
            val = 0.d0
         else
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_3d



      subroutine grid_scalar_interp_3d_p(this, x, val)
      class(grid_scalar_field) :: this
      class(space), intent(in) :: x
      real(8), intent(out) :: val
      type(grid_indices) :: indices
      integer :: marker
      real(8) :: y(8)

      associate( x1=>x%x1, x2=>x%x2, x3=>x%x3, &
               & j1=>indices%j1, j2=>indices%j2, j3=>indices%j3, &
               & nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & aa=>this%aa )

         call this%search(x, indices, marker)

         if( marker==1 ) then
            call packvec3d(y, j1, j2, j3, aa)
            call logblnint3d(y, xx1(j1), xx1(j1+1), xx2(j2), xx2(j2+1), &
                           & xx3(j3), xx3(j3+1), x1, x2, x3, val)
         elseif( marker==0 .or. marker==-2 ) then
            val = this%avarray(indices, aa)
            if( val==NoData ) then
               call grid_scalar_interp_error(x, indices, marker)
            elseif( val<0.d0 ) then
               val = 0.d0
            end if
         elseif( marker==-1 ) then
            val = 0.d0
         else
            call grid_scalar_interp_error(x, indices, marker)
         end if

      end associate
      end subroutine grid_scalar_interp_3d_p



      subroutine grid_scalar_interp_error(x, indices, marker)
      class(space), intent(in) :: x
      type(grid_indices), intent(in) :: indices
      integer, intent(in) :: marker
      print *,' mod_grid_scalar_field procedure grid_scalar_interp_error:'
      print *,'  the interpolating points appears out of domain'
      print *,'  x = (',x%x1,',',x%x2,',',x%x3,')'
      print *,'  indices = (',indices%j1,',',indices%j2,',',indices%j3,')'
      print *,'  marker =',marker
      stop
      end subroutine grid_scalar_interp_error      



      subroutine grid_scalar_clean(this)
      class(grid_scalar_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_scalar_clean



      subroutine grid_scalar_finalization(this)
      type(grid_scalar_field) :: this

      call this%grid_dealloc()
      call this%data_dealloc()
      if( associated(this%get) ) nullify(this%get)

      end subroutine grid_scalar_finalization

!
!
!      subroutine interp_1d_point_index(this, x, j0, n)
!!     this subroutine finds the n(<=Nq) grid nodes, starting from index j0, most
!!     closely surrounding the point x in the 1d grid. this routine is used to
!!     select the data points for polynomial interpolations in the following
!!     grid_scalar_interp_1d* procedures.
!      class(grid_scalar_field) :: this
!      class(space), intent(in) :: x
!      integer, intent(out) :: j0, n
!!     nominal number of points used in polynomial interpolation (1 ~ 10)
!!     linear interpolation and extrapolation
!      integer, parameter :: Nq = 2
!      type(grid_indices) :: indices
!!     j0 is the starting index of the n available points for interpolation
!      integer :: marker
!
!      associate( x1=>x%x1, j1=>indices%j1, nx1=>this%nx1, xx1=>this%xx1, &
!               & mask=>this%mask )
!         
!!        find the points used by interpolation
!         call this%search(x, indices, marker)
!
!         if( marker==1 ) then
!!           find nq points around j1
!            call find_np_points(j1, Nq, j0, n)
!         elseif( marker==0 .or. marker==-2 ) then
!            if( j1<nx1 ) then
!               if( mask(j1+1,1,1)==1 ) then
!!                 find nq points to the right of j1
!                  j0 = j1 + 1
!                  call count_to_the_right(j0, Nq, n)
!               elseif( mask(j1,1,1)==1 ) then
!!                 find nq points to the left of j1
!                  j0 = j1
!                  call count_to_the_left(j0, Nq, n)
!               else
!!                 the 1D mask domain is not allowed to be disconnected, and the
!!                 interpolation point x is not allowed to be out of the boundaries
!                  call grid_scalar_interp_error(x, indices, marker)
!               end if
!            else!( j1==nx1 )
!               if( mask(j1,1,1)==1 ) then
!                  j0 = j1
!                  call count_to_the_left(j0, Nq, n)
!               else
!                  call grid_scalar_interp_error(x, indices, marker)
!               end if
!            end if
!         else!( marker==-1 )
!            call grid_scalar_interp_error(x, indices, marker)
!         end if
!
!      end associate
!      contains
!
!         subroutine find_np_points(j, np, j0, n)
!!        given an index j and a positive integer np, this subroutine tries to find
!!        the nearest consecutive np points around j that are with mask = 1. the
!!        the starting index of these points is returned as j0, and the actually
!!        available number of points n (<=np). this subroutine assumes that the
!!        domain of mask = 1 points is simply connected, and that at least the j and
!!        j+1 points are with mask = 1. a few example results, assuming that all np
!!        points are available, are shown as below:
!!           (1) np = 1: j0 = j, available points {j}
!!           (2) np = 2: j0 = j, available points {j, j+1}
!!           (3) np = 3: j0 = j-1, available points {j-1, j, j+1}
!!           (4) np = 4: j0 = j-1, available points {j-1, j, j+1, j+2}
!         integer, intent(in) :: j, np
!         integer, intent(out) :: j0, n
!         integer :: jl, jr
!         integer :: i, i_nearest
!
!         associate( nx1=>this%nx1, mask=>this%mask )
!!           left- and right-bound of the np points around j
!            jl = j - (np-1)/2
!            jr = j + np/2
!
!!           search through these np points for the nearest one to j (if any) that
!!           is out of the grid bounds or with mask = 0
!            i_nearest = jl - 1
!            do i=jl,jr
!               if( i<1 .or. i>nx1 ) then
!                  if( abs(i-j)<=abs(i_nearest-j) ) i_nearest = i
!               elseif( mask(i,1,1)==0 ) then
!                  if( abs(i-j)<=abs(i_nearest-j) ) i_nearest = i
!               end if
!            end do
!
!            if( i_nearest==jl-1 ) then
!!              if there is no such a point, all these np points are available
!               j0 = jl
!               n = np
!            else!( i_nearest>jl-1 )
!               if( i_nearest<j ) then
!!                 if such a point is to the left of j, count np points to its right
!                  j0 = i_nearest + 1
!                  call count_to_the_right(j0, np, n)
!               else!( i_nearest>j+1 )
!!                 if such a point is to the right of j, count np points to its left
!                  j0 = i_nearest - 1
!                  call count_to_the_left(j0, np, n)
!               end if
!            end if
!         end associate
!         end subroutine find_np_points
!
!
!         subroutine count_to_the_right(j0, np, n)
!!        count np points with mask = 1 starting from j0 to its right, the counting
!!        result is n.
!         integer :: j0, np, n
!         integer :: i
!         associate( nx1=>this%nx1, mask=>this%mask )
!            n = 0
!            do i=j0,j0+np-1
!               if( i<=nx1 .and. mask(i,1,1)==1 ) then
!                  n = n + 1
!               else
!                  exit
!               end if
!            end do
!         end associate
!         end subroutine count_to_the_right
!
!
!         subroutine count_to_the_left(j0, np, n)
!!        count np points with mask = 1 starting from j0 to its left, the counting
!!        result is n. j0 is then updated to the index of the leftmost point.
!         integer :: j0, np, n
!         integer :: i
!         associate( mask=>this%mask )
!            n = 0
!            do i=j0,j0-np+1,-1
!               if( i>=1 .and. mask(i,1,1)==1 ) then
!                  n = n + 1
!               else
!                  exit
!               end if
!            end do
!            j0 = j0 - n + 1
!         end associate
!         end subroutine count_to_the_left
!
!      end subroutine interp_1d_point_index
!
end module mod_grid_scalar_field


