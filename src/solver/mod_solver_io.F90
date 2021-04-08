module mod_solver_io
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 01/03/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_typedef_params, only: STRLEN, BUFFER
      implicit none
      save
      public
      private Funit_tmp, fname_tmp
! This module contains the input/output procedures of the solver -- it is the link
! between the solver and the storage media. The solver reads in a solution input
! file to know the locations and times of requested solutions. The solution input
! comes with two modes. In the 1st mode, a list of spacetimes is provided to the
! solver -- their times may be different. And the solver seeks for solutions at this
! list of spacetimes. Contents of the input file would be:
!
!    mode, np, xt_list%x1, xt_list%x2, xt_list%x3, xt_list%t
!      ^    ^
!      |    |
!      1   number of entries in xt_list
!
! After the solving efforts, an output file is generated. The output file has a
! header that is a transcript of the contents of the input file. The solutions and
! errors are then appended to the header, as:
!
!    mode, np, xt_list%x1, xt_list%x2, xt_list%x3, xt_list%t, sol(np,2), err(np,2)
!    
! The contents of xt_list%x1 through err(np,2) can also be viewed as storing the
! following matrix in column major
!
!     /   x1,1     x2,1     x3,1     t,1     init,1     sol,1     0     err,1   \
!     |   x1,2     x2,2     x3,2     t,2     init,2     sol,2     0     err,2   |
!     |   x1,3     x2,3     x3,3     t,3     init,3     sol,3     0     err,3   |
!     |     .        .        .       .         .         .       .       .     |
!     |     .        .        .       .         .         .       .       .     |
!     \   x1,np    x2,np    x3,np    t,np    init,np    sol,np    0     err,np  /
!
! Note that there are two columns in the solution array sol(np,2), in which the
! first column stores the initial conditions, and the second column stores the
! solutions. In the error array err(np,2), the first column are zeros because the
! initial values are exact.
!
! In the 2nd input mode, a Cartesian grid of spatial positions and a list of time
! stamps are provided to the solver. Solutions are sought over the spatial grid at
! every time stamp in the t_list. As in the 1st mode, the output file has the
! solutions and the errors appended to a transcript of the input file (but with
! modification, see below), whose format is (a 2D example):
!
!    mode, ndim, nx1, nx2, nt, xx1, xx2, tt, sol(nx1,nx2,nt), err(nx1,nx2,nt)
!      ^     ^              ^             ^             ^
!      |     |              |             |             |
!      2  2 in this  output nt is   tt=(0,t_list)  sol(nx1,nx2,1) are 
!         example    1 + input nt                  initial conditions
!
! In this case, length of xt_list is np = nx1*nx2.
!
! The solution input file is read by the procedure solver_io_read_input_file, and
! the output file is processed by the procedures solver_io_open_output_file, solver_
! io_export, and solver_io_close_output_file.

      integer, parameter :: Funit_out = 99
      integer, parameter :: Funit_tmp = 98
      character(STRLEN) :: fname_out, fname_tmp

      integer, protected :: np
      integer, protected :: io_mode, io_npd, io_nx1, io_nx2, io_nx3, io_nt
      real(8), protected, allocatable :: x1_list(:), x2_list(:), x3_list(:), &
                                       & t_list(:)

contains

      subroutine solver_io_os_check()
#ifdef _WIN32
      print *,' mod_solver_io procedure solver_io_os_check:'
      print *,'  Operating system is Windows.'
      print *,'  This program requires execution on Unix/Linux.'
      stop
#endif
#ifdef __APPLE__
      print *,' mod_solver_io procedure solver_io_os_check:'
      print *,'  Operating system is Mac OSX.'
      print *,'  This program requires execution on Unix/Linux.'
      stop
#endif
      write(6,'(A)',advance='no') ' Operating system is '
      call execute_command_line('uname -ms');
      end subroutine solver_io_os_check



      subroutine solver_io_read_input_file(fname)
      character(*), intent(in) :: fname
      integer :: funit

      open(newunit=funit, file=trim(fname), access='stream', status='old', &
         & form='unformatted', err=100)
      print '(A)',' Reading solution input from file '//trim(fname)
      goto 101
100   print *,' mod_solver_io procedure solver_io_read_input_file:'
      print '(A)','  failed to open file '//trim(fname)
      stop
101   continue

      read(funit) io_mode
      print '(A,I2)','  I/O mode =',io_mode
      flush(6)

      read(funit) io_npd
      if( io_npd<=0 ) then
         print *,' mod_solver_io procedure solver_io_read_input_file:'
         print *,'  input np/ndim (=',io_npd,') must be > 0'
         close(funit)
         stop
      end if

      select case( io_mode )
      case( 1 )
         np = io_npd
         io_nx1 = np
         io_nx2 = np
         io_nx3 = np
         io_nt = np
         allocate( x1_list(np), x2_list(np), x3_list(np), t_list(np) )
         read(funit) x1_list, x2_list, x3_list, t_list
      case( 2 )
         select case( io_npd )
         case( 1 )
            read(funit) io_nx1
            io_nx2 = 1
            io_nx3 = 1
         case( 2 )
            read(funit) io_nx1, io_nx2
            io_nx3 = 1
         case( 3 )
            read(funit) io_nx1, io_nx2, io_nx3
         case default
            print *,' mod_solver_io procedure solver_io_read_input_file:'
            print *,'  input ndim (=',io_npd,') must be 1, 2 or 3'
            close(funit)
            stop
         end select

         if( io_nx1>0 .and. io_nx2>0 .and. io_nx3>0 ) then
            np = io_nx1 * io_nx2 * io_nx3
            allocate( x1_list(io_nx1), x2_list(io_nx2), x3_list(io_nx3) )
         else
            print *,' mod_solver_io procedure solver_io_read_input_file:'
            print *,'  wrong input file nx1 =',io_nx1,' or nx2 =',io_nx2, &
                  & ' or nx3 =',io_nx3
            close(funit)
            stop
         end if

         read(funit) io_nt
         if( io_nt>0 ) then
            allocate( t_list(io_nt) )
         else
            print *,' mod_solver_io procedure solver_io_read_input_file:'
            print *,'  wrong input file nt =',io_nt
            close(funit)
            stop
         end if

         select case( io_npd )
         case( 1 )
            read(funit) x1_list
            x2_list = 0.d0
            x3_list = 0.d0
         case( 2 )
            read(funit) x1_list, x2_list
            x3_list = 0.d0
         case( 3 )
            read(funit) x1_list, x2_list, x3_list
         end select

         read(funit) t_list
      case default
         print *,' mod_solver_io procedure solver_io_read_input_file:'
         print *,'  unrecognized I/O mode =',io_mode
         close(funit)
         stop
      end select

      close(funit)
      end subroutine solver_io_read_input_file



      subroutine solver_io_open_output_file(fname)
      character(*), intent(in) :: fname
      real(8), allocatable :: t_list_out(:)

!     open the solution output file 
      fname_out = trim(fname)
      open(unit=Funit_out, file=trim(fname_out), access='stream', &
         & status='replace', form='unformatted', err=100)
      print '(A)',' Solution output will be written to file '//trim(fname_out)
      flush(6)
      goto 101
100   print *,' mod_solver_io procedure solver_io_open_output_file:'
      print '(A)','  failed to open file '//trim(fname_out)
      stop
101   continue

!     transcribe the contents of the solution input file
      write(Funit_out) io_mode, io_npd

      select case( io_mode )
      case( 1 )
         write(Funit_out) x1_list, x2_list, x3_list, t_list
      case( 2 )
         allocate( t_list_out(0:io_nt) )
         t_list_out(0) = 0.d0
         t_list_out(1:io_nt) = t_list
         select case( io_npd )
         case( 1 )
            write(Funit_out) io_nx1, io_nt+1, x1_list, t_list_out
         case( 2 )
            write(Funit_out) io_nx1, io_nx2, io_nt+1, x1_list, x2_list, t_list_out
         case( 3 )
            write(Funit_out) io_nx1, io_nx2, io_nx3, io_nt+1, x1_list, x2_list, &
                       & x3_list, t_list_out
         case default
            print *,' mod_solver_io procedure solver_io_open_output_file:'
            print *,'  wront io_npd =',io_npd
            close(Funit_out)
            stop
         end select
         deallocate( t_list_out )
      case default
         print *,' mod_solver_io procedure soler_io_open_output_file:'
         print *,'  unrecognized I/O mode =',io_mode
         close(Funit_out)
         stop
      end select

      flush(Funit_out)

!     open a temporary output file to hold the errors. contents of this file will
!     be appended to the solution output file after the solving efforts.
      fname_tmp = trim(fname)//'_tmp.dat'
      open(unit=Funit_tmp, file=trim(fname_tmp), access='stream', &
         & status='replace', form='unformatted', err=200)
      goto 201
200   print *,' mod_solver_io procedure solver_io_open_output_file:'
      print '(A)','  failed to open file '//trim(fname_tmp)
      close(Funit_out)
      stop
201   continue

      end subroutine solver_io_open_output_file
      


      subroutine solver_io_export(sol_list, err_list)
      real(8), intent(in) :: sol_list(:), err_list(:)
      logical :: fopened

      inquire(unit=Funit_out, opened=fopened)
      if( .not.fopened ) then
         print *,' mod_solver_io procedure solver_io_export:'
         print *,'  output file '//trim(fname_out)//' is not opened'
         stop
      end if
      inquire(unit=Funit_tmp, opened=fopened)
      if( .not.fopened ) then
         print *,' mod_solver_io procedure solver_io_export:'
         print *,'  temporary file '//trim(fname_tmp)//' is not opened'
         stop
      end if

      if( size(sol_list)==np .and. size(err_list)==np ) then
         write(Funit_out) sol_list(:)
         flush(Funit_out)
         write(Funit_tmp) err_list(:)
         flush(Funit_tmp)
      else
         print *,' mod_solver_io procedure solver_io_export:'
         print *,'  wrong size of input list:'
         print *,'  np =',np
         print *,'  size(sol_list) =',size(sol_list)
         print *,'  size(err_list) =',size(err_list)
         close(Funit_out)
         close(Funit_tmp)
         stop
      end if

      end subroutine solver_io_export



      subroutine solver_io_close_output_file()
      character(BUFFER) :: cmd
      integer :: stat

!     release memory claimed by solver_io_open_output_file
      if( allocated(x1_list) ) deallocate( x1_list )
      if( allocated(x2_list) ) deallocate( x2_list )
      if( allocated(x3_list) ) deallocate( x3_list )
      if( allocated(t_list) ) deallocate( t_list )

!     close the output file and the temporary file
      close(Funit_out, err=100)
      goto 101
100   print *,' mod_solver_io procedure solver_io_close_output_file:'
      print *,'  failed to close file '//trim(fname_out)
101   continue

      close(Funit_tmp, err=200)
      goto 201
200   print *,' mod_solver_io procedure solver_io_close_output_file:'
      print *,'  failed to close file '//trim(fname_tmp)
201   continue

!     append the contents of fname_tmp to fname_out
      cmd = 'cat '//trim(fname_tmp)//' >> '//trim(fname_out)
      call execute_command_line(cmd, cmdstat=stat)

      if( stat==0 ) then
!        if successful, delete the temporary file
         cmd = 'rm -f '//trim(fname_tmp)
         call execute_command_line(cmd)
         print '(/A)',' Solution output file written.'
      else
         print *,' mod_solver_io procedure solver_io_close_output_file:'
         print *,'  unable to append '//trim(fname_tmp)//' to '//trim(fname_out)
         print *,'  keeping temporary file '//trim(fname_tmp)
      end if

      end subroutine solver_io_close_output_file

end module mod_solver_io




