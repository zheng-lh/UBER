module uber
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
      use mod_typedef_params, only: space, spacetime, STRLEN, initialize_params
      use mod_domain, only: initialize_domain, finalize_domain
      use mod_equation, only: initialize_equation, finalize_equation
      use mod_solver
      implicit none
      public
! This module provides the interface between the Boltzmann equation solver and the
! program calling it. Usage of the solver is:
!    ...
!    use uber
!    ...
!    call initialize_uber(FNAME_IN, FNAME_OUT)
!    call solver(ncp=NCP, icp=ICP, istat=ISTAT)
!    call finalize_uber()
!    ...
! for the higher level solver routine SOLVER, or:
!    ...
!    use uber
!    ...
!    call initialize_uber()
!    call solver_engine(xt_start, 0.d0, solution, error, istat=ISTAT)
!    call finalize_uber()
!    ...
! for the lower level solver SOLVER_ENGINE. Descriptions of the optional input
! arguments NCP, ICP, and ISTAT are given in mod_solver.F90.

contains

      subroutine initialize_uber(fname_in, fname_out)
      character(*), intent(in), optional :: fname_in, fname_out
      call initialize_params()
      call initialize_domain()
      call initialize_equation()
      if( present(fname_in) .and. present(fname_out)) then
         call initialize_solver(fname_in, fname_out)
      elseif( .not.present(fname_in) .and. .not.present(fname_out) ) then
         call initialize_solver()
      else
         print *,' UBER procedure initialize_uber:'
         print *,'  missing one of two input arguments (file names)'
         stop
      end if
      end subroutine initialize_uber



      subroutine finalize_uber()
      call finalize_solver()
      call finalize_equation()
      call finalize_domain()
      end subroutine finalize_uber

end module uber


