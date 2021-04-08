module mod_grid_typedef
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) added type-bound deferred procedures grid_set_pointers, grid_data_alloc and
!    grid_data_dealloc
!
! 2) removed the type-bound procedure grid_continuation. philosophy behind this move
!    is that any grid data preparation is now regarded as a pre-processing of input
!    to this model, rather than a part of this model.
!
! 3) removed the procedure pointers search and avarray, which are now procedures
!
! 4) added procedure print_grid and edited procedure print_mask
!
! Liheng Zheng, 11/16/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, only: space, STRLEN, NoData
      use interpolation
      implicit none
      public
! this module defines a prototype of a Cartesian grid object and its basic
! operations, including searching for grid indices, averaging over grid cell node
! values, and calculating the 1st-order grid derivatives.

      type :: grid_indices
         integer :: j1
         integer :: j2
         integer :: j3
      end type grid_indices

      type, abstract :: grid
!        grid dimension number: 1, 2 or 3
         integer :: ndim
         integer :: nx1, nx2, nx3
         real(8), dimension(:), allocatable :: xx1, xx2, xx3
!        mask value 1 indicates valid data is available on the grid node, 0
!        otherwise
!        for a 1d grid, dimension of mask is (nx1,   1,   1)
!        for a 2d grid, dimension of mask is (nx1, nx2,   1)
!        for a 3d grid, dimension of mask is (nx1, nx2, nx3)
         integer, dimension(:,:,:), allocatable :: mask
      contains
         procedure :: grid_alloc
         procedure :: grid_dealloc
         procedure(grid_op), deferred :: set_pointers
         procedure(grid_op), deferred :: data_alloc
         procedure(grid_op), deferred :: data_dealloc
         procedure(grid_op), deferred :: set_mask
         procedure :: search => grid_search
         procedure :: avarray => grid_avarray
         procedure :: print_grid => grid_print_grid
         procedure :: print_mask => grid_print_mask
         procedure :: grid_deriv
      end type grid

      abstract interface
         subroutine grid_op(this)
         import :: grid
         class(grid) :: this
         end subroutine grid_op
      end interface

contains

      subroutine grid_alloc(this, ndim, nx1, arg_nx2, arg_nx3)
      class(grid) :: this
      integer, intent(in) :: ndim
      integer, intent(in) :: nx1
      integer, intent(in), optional :: arg_nx2, arg_nx3
      integer :: nx2, nx3

      if( present(arg_nx2) ) then
         nx2 = arg_nx2
      else
         nx2 = 1
      end if

      if( present(arg_nx3) ) then
         nx3 = arg_nx3
      else
         nx3 = 1
      end if

      if( nx1<1 .or. nx2<1 .or. nx3<1 ) then
         print *,' mod_grid_typedef procedure grid_alloc:'
         print *,'  invalid input grid sizes: (nx1, nx2, nx3) = (', &
               & nx1,',',nx2,',',nx3,')'
         stop
      end if

      this%ndim = ndim
      this%nx1 = nx1
      this%nx2 = nx2
      this%nx3 = nx3

      allocate( this%xx1(nx1), this%xx2(nx2), this%xx3(nx3) )
      allocate( this%mask(nx1,nx2,nx3) )

      end subroutine grid_alloc



      subroutine grid_dealloc(this)
      class(grid) :: this
      if( allocated(this%xx1) ) deallocate( this%xx1 )
      if( allocated(this%xx2) ) deallocate( this%xx2 )
      if( allocated(this%xx3) ) deallocate( this%xx3 )
      if( allocated(this%mask) ) deallocate( this%mask )
      end subroutine grid_dealloc



      subroutine grid_search(this, x, indices, marker)
!     This routine searches in 1d, 2d or 3d grid for the indices (j1, j2, j3)
!     of the lowest-index corner of the grid cell that contains the given space
!     coordinates x = (x1, x2, x3). An integer flag MARKER is also returned to
!     tell the positional property of x.
!     MARKER =  1, if all corners of the grid cell containing x have mask = 1
!               0, if at least one but not all corners have mask = 1
!              -1, if no corner has mask = 1
!              -2, if x is out of grid bounds
      class(grid) :: this
      class(space), intent(in) :: x
      type(grid_indices), intent(out) :: indices
      integer, intent(out) :: marker
      integer :: i, j, k, ni, nj, nk, mask_sum, mask_sum_max

      indices = grid_indices(1, 1, 1)

      associate( nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, &
               & xx1=>this%xx1, xx2=>this%xx2, xx3=>this%xx3, &
               & mask=>this%mask, ndim=>this%ndim )

         select case( ndim )
         case( 1 )
            call locate(xx1, nx1, x%x1, indices%j1)
            if( indices%j1==0 .or. indices%j1==nx1 ) then
               marker = -2
               return
            end if
            ni = 1
            nj = 0
            nk = 0
         case( 2 )
            call locate(xx1, nx1, x%x1, indices%j1)
            call locate(xx2, nx2, x%x2, indices%j2)
            if( indices%j1==0 .or. indices%j1==nx1 .or. &
              & indices%j2==0 .or. indices%j2==nx2      ) then
               marker = -2
               return
            end if
            ni = 1
            nj = 1
            nk = 0
         case( 3 )
            call locate(xx1, nx1, x%x1, indices%j1)
            call locate(xx2, nx2, x%x2, indices%j2)
            call locate(xx3, nx3, x%x3, indices%j3)
            if( indices%j1==0 .or. indices%j1==nx1 .or. &
              & indices%j2==0 .or. indices%j2==nx2 .or. &
              & indices%j3==0 .or. indices%j3==nx3      ) then
               marker = -2
               return
            end if
            ni = 1
            nj = 1
            nk = 1
         case default
            print *,' mod_grid_typedef procedure grid_search:'
            print *,'  unrecognized ndim =',ndim
            stop
         end select

         mask_sum = 0
         do k=0,nk
            do j=0,nj
               do i=0,ni
                  mask_sum = mask_sum + &
                  & mask(indices%j1+i, indices%j2+j, indices%j3+k)
               end do
            end do
         end do
      end associate

      mask_sum_max = 2**(ni + nj + nk)
      if( mask_sum==mask_sum_max ) then
         marker = 1
      elseif( mask_sum<mask_sum_max .and. mask_sum>0 ) then
         marker = 0
      elseif( mask_sum==0 ) then
         marker = -1
      else
         print *,' mod_grid_typedef procedure grid_search:'
         print *,'  wrong mask_sum occurred'
         print '(A,3(F10.4,A1))','  x = (',x%x1,',',x%x2,',',x%x3,')'
         print '(A,3(I5,A1))','  indices = (',indices%j1,',', &
                                            & indices%j2,',',indices%j3,')'
         print *,' mask_sum =',mask_sum
         stop
      end if

      end subroutine grid_search



      function grid_avarray(this, indices, array) result(ave)
!     This function returns the average of the array values (with mask = 1) from
!     the grid cell with INDICES indicating its lowest-index corner. This function
!     serve as a simple extrapolation means in multiple dimensions. N.B., this
!     function allows for extrapolation beyond the grid bounds.
      class(grid) :: this
      type(grid_indices), intent(in) :: indices
      real(8), intent(in) :: array(:,:,:)
      real(8) :: ave
      integer :: i, j, k, ni, nj, nk, cnt
      type(grid_indices) :: my_indices
      associate( j1=>my_indices%j1, j2=>my_indices%j2, j3=>my_indices%j3, &
               & nx1=>this%nx1, nx2=>this%nx2, nx3=>this%nx3, mask=>this%mask )

      my_indices = indices

      if( j1==0 .or. j1==nx1 ) then
         ni = 0
         if( j1==0 ) j1 = 1
      else
         ni = 1
      end if

      if( j2==0 .or. j2==nx2 ) then
         nj = 0
         if( j2==0 ) j2 = 1
      else
         nj = 1
      end if

      if( j3==0 .or. j3==nx3 ) then
         nk = 0
         if( j3==0 ) j3 = 1
      else
         nk = 1
      end if

      cnt = 0
      ave = 0.d0

      do k=0,nk
         do j=0,nj
            do i=0,ni
               if( mask(j1+i, j2+j, j3+k)==1 ) then
                  cnt = cnt + 1
                  ave = ave + array(j1+i, j2+j, j3+k)
               end if
            end do
         end do
      end do

      if( cnt>0 ) then
         ave = ave/dble(cnt)
      else
         ave = NoData
      end if

      end associate
      end function grid_avarray



      subroutine grid_print_grid(this, fname)
      class(grid) :: this
      character(*), optional :: fname
      integer :: io_unit, i

      if( present(fname) ) then
         open(newunit=io_unit, file=trim(fname), form='formatted', status='replace')
      else
         io_unit = 6
      end if

      write(io_unit,'(/A)') ' ----------- Grid print out  -----------'
      write(io_unit,'(A,I2)') '  ndim =',this%ndim
      write(io_unit,'(3(A7,I6))') '  nx1 =',this%nx1,', nx2 =',this%nx2, &
                                & ', nx3 =',this%nx3
      write(io_unit,'(/A)') '  xx1 ='
      write(io_unit,'(10F8.2)') (this%xx1(i), i=1,this%nx1)
      write(io_unit,'(/A)') '  xx2 ='
      write(io_unit,'(10F8.2)') (this%xx2(i), i=1,this%nx2)
      write(io_unit,'(/A)') '  xx3 ='
      write(io_unit,'(10F8.2)') (this%xx3(i), i=1,this%nx3)
      write(io_unit,'(/A)') ' ---------------------------------------'

      if( present(fname) ) close(io_unit)

      end subroutine grid_print_grid



      subroutine grid_print_mask(this, fname)
      class(grid) :: this
      character(*), optional :: fname
      character(STRLEN) :: fmtstr
      integer :: io_unit, i, j, k

      if( present(fname) ) then
         open(newunit=io_unit, file=trim(fname), form='formatted', status='replace')
      else
         io_unit = 6
      end if

      write(unit=fmtstr,fmt='(I3,A)') this%nx2,'I3'
      write(io_unit,'(/A)') ' ----------- Mask print out  -----------'

      do k=1,this%nx3
         write(io_unit,'(''k = '',I3)') k
         write(io_unit,'(''i\j'','//trim(fmtstr)//')') (j, j=1,this%nx2)
         do i=1,this%nx1
            write(io_unit,'(I3,'//trim(fmtstr)//')') i, (this%mask(i,j,k), &
                                                       & j=1,this%nx2)
         end do
         write(io_unit,*)
      end do

      write(unit=io_unit,fmt='(A)') ' ---------------------------------------'
      
      if( present(fname) ) close(io_unit)

      end subroutine grid_print_mask



      subroutine grid_deriv(this, array, gdim, deriv, deriv_sub)
!     The grid_deriv routine calculates the 1st-order derivatives of array
!     along the given dimension gdim over all mask = 1 grid nodes, with user
!     specified numerical differentiation subroutine deriv_sub.
      class(grid) :: this
      integer, intent(in) :: gdim
      real(8), intent(in) :: array(:,:,:)
      real(8), intent(inout) :: deriv(:,:,:)
      integer, pointer :: p, q, r
      integer, target :: i, j, k
      integer :: ni, nj, nk
      integer :: cnt0, cnt1
      real(8), dimension(:), allocatable :: xa, ya, ta, xx
!     interface to subroutines to calculate derivatives. available options are:
!     interpolation: monotone_tangent, fss002
      interface
         subroutine deriv_sub(xx, n, yy, tt)
         integer :: n
         real(8) :: xx(n), yy(n), tt(n)
         end subroutine deriv_sub
      end interface

      deriv = 0.d0

      select case( gdim )
      case( 1 )
         ni = this%nx1
         nj = this%nx2
         nk = this%nx3
         p => i
         q => j
         r => k
         allocate( xx(ni) )
         xx = this%xx1
      case( 2 )
         ni = this%nx2
         nj = this%nx3
         nk = this%nx1
         q => i
         r => j
         p => k
         allocate( xx(ni) )
         xx = this%xx2
      case( 3 )
         ni = this%nx3
         nj = this%nx1
         nk = this%nx2
         r => i
         p => j
         q => k
         allocate( xx(ni) )
         xx = this%xx3
      case default
         print *,' mod_grid_typedef procedure grid_deriv:'
         print *,'  wrong gdim =',gdim
         stop
      end select

      do k=1,nk
         do j=1,nj

            cnt0 = -1
            do i=1,ni
               if( this%mask(p,q,r)==1 ) then
                  cnt0 = i - 1
                  exit
               end if
            end do
            if( cnt0==-1 ) cycle

            cnt1 = 0
            do i=cnt0+1,ni
               if( this%mask(p,q,r)==1 ) then
                  cnt1 = cnt1 + 1
               else
                  exit
               end if
            end do

            if( (cnt0 + cnt1)>ni ) then
               print *,' mod_grid_typedef procedure grid_deriv:'
               print *,'  cnt0 + cnt1 =',cnt0+cnt1,'> ni =',ni
               stop
            end if

            if( cnt1>1 ) then
               allocate( xa(cnt1), ya(cnt1), ta(cnt1) )
               xa = 0.d0
               ya = 0.d0
               ta = 0.d0
               xa(:) = xx((cnt0+1):(cnt0+cnt1))

               do i=cnt0+1,cnt0+cnt1
                  ya(i-cnt0) = array(p,q,r)
               end do

               if( cnt1>=3 ) then
                  call deriv_sub(xa, cnt1, ya, ta)
               else
                  ta(:) = (ya(2) - ya(1))/(xa(2) - xa(1))
               end if

               do i=cnt0+1,cnt0+cnt1
                  deriv(p,q,r) = ta(i-cnt0)
               end do

               deallocate(xa, ya, ta)
            end if

         end do
      end do

      nullify(p, q, r)
      deallocate( xx )

      end subroutine grid_deriv
!
!
!
!      subroutine grid_contn(this, array, gdim, direction, np)
!!     this routine continuates multi-dimensional data array along the specified
!!     dimension gdim toward either the first (begin) or the last (end) grid index,
!!     with polynomial extrapolation degree np - 1 (if there are at least np data
!!     points available; otherwise the actually available number of points are used).
!      class(grid) :: this
!      integer, intent(in) :: gdim
!      real(8), intent(inout) :: array(:,:,:)
!!     nubmer of points used in polynomial extrapolation (np <= 10)
!      integer, intent(in) :: np
!!     direction = 'B' for continuates toward the first index
!!               = 'E' for continuates toward the last index
!      character(1), intent(in) :: direction
!      integer, pointer :: p, q, r
!      integer, target :: i, j, k
!      integer :: ni, nj, nk
!      integer :: n, inc, ind
!      real(8) :: x, y, dy
!      real(8), dimension(:), allocatable :: xa, ya, xx
!      procedure(func1), pointer :: move
!
!      abstract interface
!         function func1(a)
!         integer :: a
!         integer :: func1
!         end function func1
!      end interface
!
!      select case( gdim )
!      case( 1 )
!         ni = this%nx1
!         nj = this%nx2
!         nk = this%nx3
!         p => i
!         q => j
!         r => k
!         allocate( xx(ni) )
!         xx = this%xx1
!      case( 2 )
!         ni = this%nx2
!         nj = this%nx3
!         nk = this%nx1
!         q => i
!         r => j
!         p => k
!         allocate( xx(ni) )
!         xx = this%xx2
!      case( 3 )
!         ni = this%nx3
!         nj = this%nx1
!         nk = this%nx2
!         r => i
!         p => j
!         q => k
!         allocate( xx(ni) )
!         xx = this%xx3
!      case default
!         print *,' mod_grid_typedef procedure grid_contn:'
!         print *,'  wrong gdim =',gdim
!         stop
!      end select
!
!      if( direction=='B' ) then
!         ind = 1
!         move => plus
!      elseif( direction=='E' ) then
!         ind = ni
!         move => minus
!      else
!         print *,' mod_grid_typedef procedure grid_contn:'
!         print *,'  unrecognized character direction ='//direction
!         stop
!      end if
!
!      do k=1,nk
!         do j=1,nj
!            i = ind
!!           cycle if the target point is out of domain
!            if( this%mask(p,q,r)==0 ) cycle
!!           count the actually available points for extrapolation
!            n = 0
!            do inc=1,np
!               i = ind + move(inc)
!               if( this%mask(p,q,r)==1 ) then
!                  n = n + 1
!               else
!                  exit
!               end if
!            end do
!
!            allocate(xa(n), ya(n))
!            xa = 0.d0
!            ya = 0.d0
!
!            do inc=1,n
!               i = ind + move(inc)
!               xa(inc) = xx(i)
!               ya(inc) = array(p,q,r)
!            end do
!
!            x = xx(ind)
!            call polint(xa, ya, n, x, y, dy)
!            array(p,q,r) = y
!
!            deallocate(xa, ya)
!         end do
!      end do
!
!      nullify(p, q, r, move)
!      deallocate( xx )
!      contains
!
!         function plus(a)
!         integer :: a
!         integer :: plus
!         plus = a
!         end function plus
!
!
!         function minus(a)
!         integer :: a
!         integer :: minus
!         minus = -a
!         end function minus
!
!      end subroutine grid_contn

end module mod_grid_typedef




