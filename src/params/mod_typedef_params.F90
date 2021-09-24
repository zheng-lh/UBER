module mod_typedef_params
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
      implicit none
      save
      public

!     a point in space
      type space
         real(8) :: x1
         real(8) :: x2
         real(8) :: x3
      end type space

!     a point in spacetime
      type, extends(space) :: spacetime
         real(8) :: t
      end type spacetime

!     default length for character string and buffer
      integer, parameter :: STRLEN = 64
      integer, parameter :: BUFFER = 128

!     number of stochastic processes per batch. NBAT has to be declared as parameter
!     because it is used as the bound of iterations in OpenMP parallel DO. the
!     OpenMP clause SCHEDULE needs to know this bound a priori. note that, NBAT has
!     to be sufficiently large to reduce the effect of parallelization overhead.
      integer, parameter :: NBAT = 2048

!     maximum fission generations allowed by the program
      integer, parameter :: MaxGenerations = 4

!     default value of a data holder for which no data is assigned (this number is
!     close to the negative large bound of Fortran single precision) 
      real(8), parameter :: NoData = -1.d38

!     Neighborhood width = NeighborFactor * dS (see below)
      real(8), parameter, private :: NeighborFactor = 4.d0

!     boundary thickness = BoundaryFactor * Neighborhood (see below)
      real(8), parameter, private :: BoundaryFactor = 1.d-3

!     lower and upper limits of temporal step size
      real(8), protected :: dtMin, dtMax
!     designated root-mean-square spatial step size: dS = sqrt(tr(A)*dt)
      real(8), protected :: dS

!     neighborhood width
      real(8), protected :: Neighborhood

!     boundary thickness. in mod_solver_ito_process.F90, the projection and
!     intersection routines find their roots within accuracy of epsBnd, so that the
!     appeared thickness of a boundary piece would be 2*epsBnd.
      real(8), protected :: epsBnd

!     relative error tolerance
      real(8), protected :: Beta
!     confidence level (Ma = 1, 80%; 2, 90%; 3, 95%)
      integer, protected :: Ma
!     maximum and minimum batch numbers
      integer, protected :: MMax, MMin
!     when batch number is less than MRed, reduce error tolerance by half to prevent
!     a premature iteration termination
      integer, protected :: MRed

!     CriticalRatio is a real number between [0, 1) that specifies the upper
!     quantile of stochastic processes for which the first fission will occur. for
!     example, CriticalRatio = 0.25 specifies that, in terms of functional value,
!     the top 25% of stochastic processes will undergo the first fission. the j-th
!     generation fissions in the upper CriticalRatio/2**j quantile.
      real(8), protected :: CriticalRatio
!     number of children in the fission descendents. setting nChildren = 1 disables
!     fission
      integer, protected :: nChildren
!     maximum number of fission generations chosen by user (root generation = 0)
      integer, protected :: nGenerations

      character(STRLEN), protected :: ParamFname = 'UBER_params.in'

contains

      subroutine initialize_params()
      integer :: fid = 11

      open(fid,file=trim(ParamFname),form='formatted',status='old',err=100)
      goto 101
100   print *,' mod_typedef_params procedure initialize_params:'
      print *,'  failed to open input file '//trim(ParamFname)
      stop
101   continue

      read(fid,*) dtMin
      read(fid,*) dtMax
      read(fid,*) dS
      read(fid,*) Beta
      read(fid,*) Ma
      read(fid,*) MMax
      read(fid,*) MMin
      read(fid,*) MRed
      read(fid,*) CriticalRatio
      read(fid,*) nChildren
      read(fid,*) nGenerations

      close(fid)

      if( dtMin>dtMax .or. dtMin<0.d0 ) then
         print *,' mod_typedef_paras procedure initialize_params:'
         print *,'  invalid dtMin (=',dtMin,') or dtMax (=',dtMax,')'
         stop
      end if
      if( dS<=0.d0 ) then
         print *,' mod_typedef_paras procedure initialize_params:'
         print *,'  invalid dS value ',dS
         stop
      end if
      if( Beta<0.d0 ) then
         print *,' mod_typedef_paras procedure initialize_params:'
         print *,'  invalid Beta value ',Beta
         stop
      end if

      Neighborhood = NeighborFactor*dS
      epsBnd = BoundaryFactor*Neighborhood

      if( Ma/=1 .and. Ma/=2 .and. Ma/=3 ) then
         print *,' mod_typedef_paras procedure initialize_params:'
         print *,'  invalid Ma value ',Ma
         stop
      end if
      if( MMax<MMin .or. MMax<=0 ) then
         print *,' mod_typedef_paras procedure initialize_params:'
         print *,'  invalid MMax (=',MMax,') or MMin (=',MMin,')'
         stop
      end if
      if( CriticalRatio<0.d0 .or. CriticalRatio>=1.d0 ) then
         print *,' mod_typedef_params procedure initialize_params:'
         print *,'  invalid CriticalRatio value ',CriticalRatio
         stop
      end if

!     to save memory, set nGenerations to 0 if fission is disabled
      if( nChildren==1 ) nGenerations = 0
      nGenerations = min(nGenerations, MaxGenerations)

      write(6,'(/A)')      ' Parameters of this simulation:'
      write(6,'(A)')       ' * * * * * * * * * * * * * * * *'
      write(6,'(A,ES9.2)') ' Min temporal step size =',dtMin
      write(6,'(A,ES9.2)') ' Max temporal step size =',dtMax
      write(6,'(A,ES9.2)') ' RMS spatial step size =',dS
      write(6,'(A,ES9.2)') ' Neighborhood width =',Neighborhood
      write(6,'(A,ES9.2)') ' Boundary thickness =',epsBnd
      write(6,'(A,F5.1,A)')' Designated relative error =',100.d0*Beta,'%'
      write(6,'(A,I3)')    ' Confidence level selection =',Ma
      write(6,'(A,I5)')    ' Min number of stochastic batches =',MMin
      write(6,'(A,I5)')    ' Max number of stochastic batches =',MMax
      write(6,'(A,I5)')    ' Number of batches with reduced error tolerance =',MRed
      write(6,'(A,F6.3)')  ' Critical ratio =',CriticalRatio
      write(6,'(A,I3)')    ' Number of children at each fission =',nChildren
      write(6,'(A,I3)')    ' Max number of fission generations =',nGenerations
      write(6,'(A,I5)')    ' Number of Ito processes per batch =',NBAT
      write(6,'(A,ES9.2)') ' Default data points will be filled with ',NoData
      write(6,'(A/)')      ' * * * * * * * * * * * * * * * *'
      flush(6)

      end subroutine initialize_params

end module mod_typedef_params

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) redefined type coefficients and type components to incorporate the added term v
!    in the Boltzmann equation (1).
!
! 2) added a new parameter NoData
!
! 3) renamed FissionRatio to CriticalRatio to better reflect its meaning.
!
! 4) some revision of the comments.
!
! 5) added parameters STRLEN and BUFFER
!
! Liheng Zheng, 01/13/2020
! 
! 1) moved epsBnd from mod_domain_typedef to this module
!
! 2) moved declarations of type coefficients and type components to mod_equation_
!    typedef module 
!
! 3) moved declaration of NBAT from mod_solver_batch module to here
!
! Liheng Zheng, 01/19/2020
!
! 1) added an integer parameter MaxGenerations
!
! 2) changed the meaning of CriticalRatio (explained in comments)
!
! 3) added validity checks in INITIALIZE_PARAMS
!
! Liheng Zheng, 03/02/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



