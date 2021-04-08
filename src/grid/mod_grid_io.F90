module mod_grid_io
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 01/13/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      implicit none
      private
      public grid_io_read_data_file
! This module contains the data input procedures for the grids. To use these
! procedures, specific data file formats are defined as follows. The data files
! shall be saved as binary files in stream I/O format, as opposed to the usual
! record-based I/O format of FORTRAN. In this way, it allows for inter-
! operability of these files by other data-processing software such as MATLAB.
!
! Example contents of an input data file containing a 2D coefficient (such as a
! component of a rank-2 tensor) are:
!
!    ndim, nx1, nx2, xx1(nx1), xx2(nx2), darr(nx1,nx2)
!      ^    ^           ^                    ^
!      |    |           |                    |
!      |    +-------+   +---------+          +--+
!      |            |             |             |
!    number of   number of   list of grid   coefficient data
!    dimensions  grid nodes  nodes in the   stored in column major
!    (2)         in the 1st  1st dimension
!                dimension
!
! The input data files are read by the generic procedure grid_io_read_data_file.
! N.B., the dummy output arguments XX? and DARR of this procedure are allocatable
! arrays, which must be deallocated explicitly outside of this module to release
! memory.

      interface grid_io_read_data_file
         module procedure grid_io_read_data_file_1d, &
                        & grid_io_read_data_file_2d, &
                        & grid_io_read_data_file_3d
      end interface

contains

      subroutine grid_io_read_data_file_1d(fname, nx1, xx1, darr)
      character(*), intent(in) :: fname
      integer, intent(out) :: nx1
      real(8), intent(out), allocatable :: xx1(:)
      real(8), intent(out), allocatable :: darr(:)
      integer :: funit, ndim

      open(newunit=funit, file=trim(fname), access='stream', status='old', &
         & form='unformatted', err=100)
      print '(/A)','  Reading data file '//trim(fname)//' ...'
      flush(6)
      goto 101
100   call grid_io_read_data_file_err(fname)
101   continue

      read(funit) ndim, nx1

      if( ndim/=1 ) then
         print *,' mod_grid_io procedure grid_io_read_data_file_1d:'
         print *,'  wrong input file ndim =',ndim
         close(funit)
         stop
      end if

      if( nx1>0 ) then
         allocate( xx1(nx1), darr(nx1) )
      else
         print *,' mod_grid_io procedure grid_io_read_data_file_1d:'
         print *,'  wrong input file nx1 =',nx1
         close(funit)
         stop
      end if

      read(funit) xx1, darr
      close(funit)

      end subroutine grid_io_read_data_file_1d



      subroutine grid_io_read_data_file_2d(fname, nx1, nx2, xx1, xx2, darr)
      character(*), intent(in) :: fname
      integer, intent(out) :: nx1, nx2
      real(8), intent(out), allocatable :: xx1(:), xx2(:)
      real(8), intent(out), allocatable :: darr(:,:)
      integer :: funit, ndim

      open(newunit=funit, file=trim(fname), access='stream', status='old', &
         & form='unformatted', err=100)
      print '(/A)','  Reading data file '//trim(fname)//' ...'
      flush(6)
      goto 101
100   call grid_io_read_data_file_err(fname)
101   continue

      read(funit) ndim, nx1, nx2

      if( ndim/=2 ) then
         print *,' mod_grid_io procedure grid_io_read_data_file_2d:'
         print *,'  wrong input file ndim =',ndim
         close(funit)
         stop
      end if
      
      if( nx1>0 .and. nx2>0 ) then
         allocate( xx1(nx1), xx2(nx2), darr(nx1,nx2) )
      else
         print *,' mod_grid_io procedure grid_io_read_data_file_2d:'
         print *,'  wrong input file nx1 =',nx1,'or nx2 =',nx2
         close(funit)
         stop
      end if

      read(funit) xx1, xx2, darr
      close(funit)

      end subroutine grid_io_read_data_file_2d



      subroutine grid_io_read_data_file_3d(fname, nx1, nx2, nx3, xx1, xx2, xx3, &
                                           & darr)
      character(*), intent(in) :: fname
      integer, intent(out) :: nx1, nx2, nx3
      real(8), intent(out), allocatable :: xx1(:), xx2(:), xx3(:)
      real(8), intent(out), allocatable :: darr(:,:,:)
      integer :: funit, ndim

      open(newunit=funit, file=trim(fname), access='stream', status='old', &
         & form='unformatted', err=100)
      print '(/A)','  Reading data file '//trim(fname)//' ...'
      flush(6)
      goto 101
100   call grid_io_read_data_file_err(fname)
101   continue

      read(funit) ndim, nx1, nx2, nx3

      if( ndim/=3 ) then
         print *,' mod_grid_io procedure grid_io_read_data_file_2d:'
         print *,'  wrong input file ndim =',ndim
         close(funit)
         stop
      end if

      if( nx1>0 .and. nx2>0 .and. nx3>0 ) then
         allocate( xx1(nx1), xx2(nx2), xx3(nx3), darr(nx1,nx2,nx3) )
      else
         print *,' mod_grid_io procedure grid_io_read_data_file_3d:'
         print *,'  wrong input file nx1 =',nx1,' or nx2 =',nx2,' or nx3 =',nx3
         close(funit)
         stop
      end if

      read(funit) xx1, xx2, xx3, darr
      close(funit)

      end subroutine grid_io_read_data_file_3d



      subroutine grid_io_read_data_file_err(fname)
      character(*), intent(in) :: fname
      print *,' mod_grid_io procedure grid_io_read_data_file:'
      print '(A)','  failed to open file '//trim(fname)
      stop
      end subroutine grid_io_read_data_file_err

end module mod_grid_io




