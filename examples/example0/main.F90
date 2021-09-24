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

      program main
      use uber
      implicit none
      real(8) :: sec0, sec1
      integer, dimension(8) :: digits0, digits1
      character(12) :: date, time, zone
      type(spacetime) :: xt_start
      real(8) :: solution, error

!     meter the starting time
      call date_and_time(date, time, zone, digits0)

      call initialize_uber()

      xt_start = spacetime(0.44d0, 0.d0, 0.d0, 5.d-2)
      call solver_engine(xt_start, 0.d0, solution, error, 'V')
      
      xt_start = spacetime(0.d0, 0.64d0, 0.d0, 5.d-2)
      call solver_engine(xt_start, 0.d0, solution, error, 'V')
      
      xt_start = spacetime(0.d0, 0.d0, 0.84d0, 5.d-2)
      call solver_engine(xt_start, 0.d0, solution, error, 'V')
      
      call finalize_uber()

!     meter the ending time
      call date_and_time(date, time, zone, digits1)
      sec0 = dble(digits0(5))*3.6d3 + dble(digits0(6))*6.d1 + &
           & dble(digits0(7)) + dble(digits0(8))*1.d-3
      sec1 = dble(digits1(5))*3.6d3 + dble(digits1(6))*6.d1 + &
           & dble(digits1(7)) + dble(digits1(8))*1.d-3
      if( sec1<=sec0 ) sec1 = sec1 + 24.d0*3600.d0

100   format(A,I2,':',I2,"'",I2,'.',I3,'"')
      write(6,fmt=100)      ' Time started: ', digits0(5:8)
      write(6,fmt=100)      ' Time ended  : ', digits1(5:8)
      write(6,'(A,F9.3,A)') ' Time used   : ', sec1 - sec0,' sec'
      flush(6)

      end program main

