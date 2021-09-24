module functions
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
! This modules provides the mathematical and physical constants, and the basic
! functions related to computations of the Earth's radiation belts. The physical
! constants are provided in three unit systems: *_mks for the meter-kilogram-second-
! ampere system (systeme international), *_cgs for the centimeter-gram-second system
! (Gaussian), and *_nat for the natural units. The natural unit system uses
! elementary charge as the unit for electric charge, speed of light in vacuum as the
! unit for speed, Earth mean radius as the unit for length, Gauss as the unit for
! magnetic intensity, and MeV as the unit for energy; other units are derived on
! these basis units. 
!
! The functions in this module could be roughly categorized in four groups:
!
!    1) relating to the magnetic dipole field, including dipole_loss_cone, and the
!       functions appearing in reference [1] SchulzY, SchulzInvy2, SchulzT, SchulzD
!       and SchulzQ;
!    2) relating to special relativity, including p_from_E, E_from_p, beta_sr, and
!       gamma_sr;
!    3) relating to particle periodic motions in dipole field, including the
!       frequencies Omega2, Omega3, and their corresponding periods Tbounce and
!       Tdrift;
!    4) relating to conversions between Universal Time and Julian days, including
!       jd2ut, ut2jd, and time_diff.
!
! References:
!   [1] M. Schulz (1991), The magnetosphere, in Geomagnetism, vol. 4, pp. 87-293, 
!       Academic Press.

      real(8), parameter :: Pi = 3.1415926535897932384626d0
      real(8), parameter :: Po180 = Pi/180.d0
!     Boltzmann constant
      real(8), parameter :: kB_mks = 1.380649d-23 ! J/K
      real(8), parameter :: kB_cgs = kB_mks*1.d7  ! erg/K
      real(8), parameter :: kB_nat = 8.617333d-11 ! MeV/K
!     reduced Planck constant
      real(8), parameter :: hbar_mks = 1.054572d-34  ! J*s
      real(8), parameter :: hbar_cgs = hbar_mks*1.d7 ! erg*s
      real(8), parameter :: hbar_nat = 6.582120d-22  ! MeV*s
!     elementary charge
      real(8), parameter :: e_mks = 1.602177d-19 ! C
      real(8), parameter :: e_cgs = 4.803204d-10 ! statcoul
      real(8), parameter :: e_nat = 1.d0
!     speed of light in vaccum
      real(8), parameter :: c_mks = 2.997925d8  ! m/s
      real(8), parameter :: c_cgs = c_mks*1.d2  ! cm/s
      real(8), parameter :: c_nat = 1.d0
!     electron mass
      real(8), parameter :: me_mks = 9.109384d-31 ! kg
      real(8), parameter :: me_cgs = me_mks*1.d3  ! g
      real(8), parameter :: me_nat = 0.5109990d0  ! MeV/c^2
!     proton mass
      real(8), parameter :: mp_mks = 1.672622d-27 ! kg
      real(8), parameter :: mp_cgs = mp_mks*1.d3  ! g
      real(8), parameter :: mp_nat = 9.382721d2   ! MeV/c^2
!     electron rest energy
      real(8), parameter :: E0e_mks = me_mks*c_mks**2 ! J
      real(8), parameter :: E0e_cgs = me_cgs*c_cgs**2 ! erg
      real(8), parameter :: E0e_nat = me_nat*c_nat**2 ! MeV
!     proton rest energy
      real(8), parameter :: E0p_mks = mp_mks*c_mks**2 ! J
      real(8), parameter :: E0p_cgs = mp_cgs*c_cgs**2 ! erg
      real(8), parameter :: E0p_nat = mp_nat*c_nat**2 ! MeV
!     energy associated with 1 eV
      real(8), parameter :: eV_mks = e_mks      ! J
      real(8), parameter :: eV_cgs = e_mks*1.d7 ! erg
      real(8), parameter :: eV_nat = 1.d-6      ! MeV
!     Earth mean radius
      real(8), parameter :: RE_mks = 6.371009d6  ! m
      real(8), parameter :: RE_cgs = RE_mks*1.d2 ! cm
      real(8), parameter :: RE_nat = 1.d0
!     Earth magnetic dipole moment
      real(8), parameter :: muE_mks = 0.311653d-4*RE_mks**3 ! T*m^3
      real(8), parameter :: muE_cgs = 0.311653d0*RE_cgs**3  ! G*cm^3
      real(8), parameter :: muE_nat = 0.311653d0*RE_nat**3  ! G*RE^3

contains

      function dipole_loss_cone(L)
!     dipole loss cone angle in radian
      real(8), intent(in) :: L
      real(8) :: dipole_loss_cone
      real(8) :: y

      if( L>1.d0 ) then
         y = 1.d0/sqrt(sqrt(4.d0*L**6 - 3.d0*L**5))
         dipole_loss_cone = asin(y)
      else
         print *,' function dipole_loss_cone:'
         print *,'  input L (=',L,') is out of domain (1, +inf).'
         stop
      end if

      end function dipole_loss_cone



      function SchulzY(y)
!     Eq. (161) in [1]. input argument y = sin(a0), where a0 is equatorial pitch
!     angle.
      real(8), intent(in) :: y
      real(8) :: SchulzY

      if( y>=0.d0 .and. y<=1.d0 ) then
         SchulzY = 2.760346d0 + 2.357194d0*y - 5.117540d0*y**(3.d0/4.d0)
      else
         print *,' function SchulzY:'
         print *,'  input argument y (=',y,') is out of domain [0, 1].'
         stop
      end if

      end function SchulzY



      function SchulzInvy2(x)
!     Eq. (132) in [1]. input argument x = sqrt(L*RE/muE)*K, where K is the second
!     adiabatic invariant defined by Eq. (130) in [1].
      real(8), intent(in) :: x
      real(8) :: SchulzInvy2

      if( x>=0.d0 ) then
         SchulzInvy2 = 1.0d0 + 1.35048d0*x - 0.030425d0*x**(4.d0/3.d0) + &
                     & 0.10066d0*x**(5.d0/3.d0) + (x/2.760346d0)**2
      else
         print *,' function SchulzInvy2:'
         print *,'  input argument x (=',x,') is out of domain [0, +inf).'
         stop
      end if

      end function SchulzInvy2



      function SchulzT(y)
!     Eq. (156) in [1]. input argument y = sin(a0).
      real(8), intent(in) :: y
      real(8) :: SchulzT

      if( y>=0.d0 .and. y<=1.d0 ) then
         SchulzT = 1.380173d0 - 0.639693d0*y**(3.d0/4.d0)
      else
         print *,' function SchulzT:'
         print *,'  input argument y (=',y,') is out of domain [0, 1].'
         stop
      end if

      end function SchulzT



      function SchulzD(y)
!     Eq. (170) in [1]. input argument y = sin(a0).
      real(8), intent(in) :: y
      real(8) :: SchulzD

      if( y>=0.d0 .and. y<=1.d0 ) then
         SchulzD = (5.520692d0 - 2.357194d0*y + 1.279385d0*y**(3.d0/4.d0))/12.d0
      else
         print *,' function SchulzD:'
         print *,'  input argument y (=',y,') is out of domain [0, 1].'
         stop
      end if

      end function SchulzD



      function SchulzQ(y)
!     Eq. (205) in [1]. input argument y = sin(a0).
      real(8), intent(in) :: y
      real(8) :: SchulzQ

      if( y>=0.d0 .and. y<=1.d0 ) then
         SchulzQ = -27.12667d0 - 45.39913*y**4 + 5.88256*y**8
      else
         print *,' function SchulzQ:'
         print *,'  input argument y (=',y,') is out of domain [0, 1].'
         stop
      end if

      end function SchulzQ



      function p_from_E(E, spec, units)
!     this routine calculates momentum from kinetic energy for particle species
!        spec = 'e' - electron (default)
!               'p' - proton
!     in the unit system 
!        units = 'M' - mks units
!                'C' - cgs units
!                'N' - natural units (default)
!     according to the formula
!
!        p = m0*c * sqrt(E'*(E' + 2)) = (1/c) * sqrt(E*(E + 2*E0)),             (1)
!
!     where E' = E/E0.      
      real(8), intent(in) :: E
      character(1), intent(in), optional :: spec, units
      real(8) :: p_from_E
      character(1) :: ispec, iunits
      real(8) :: c, E0

      if( present(spec) ) then
         ispec = spec
         if( ispec/='e' .and. ispec/='p' ) then
            print *,' function p_from_E:'
            print *,'  invalid input SPEC = ',ispec
            stop
         end if
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
         if( iunits/='M' .and. iunits/='C' .and. iunits/='N' ) then
            print *,' function p_from_E:'
            print *,'  invalid input UNITS = ',iunits
            stop
         end if
      else
         iunits = 'N'
      end if

      select case( iunits )
      case( 'N' )
         c = c_nat
         select case( ispec )
         case( 'e' )
            E0 = E0e_nat
         case( 'p' )
            E0 = E0p_nat
         end select
      case( 'M' )
         c = c_mks
         select case( ispec )
         case( 'e' )
            E0 = E0e_mks
         case( 'p' )
            E0 = E0p_mks
         end select
      case( 'C' )
         c = c_cgs
         select case( ispec )
         case( 'e' )
            E0 = E0e_cgs
         case( 'p' )
            E0 = E0p_cgs
         end select
      end select

      if( E>=0.d0 ) then
         p_from_E = sqrt(E*(E + 2.d0*E0))/c
      else
         print *,' function p_from_E:'
         print *,'  input argument E (=',E,') is out of domain [0, +inf).'
         stop
      end if

      end function p_from_E



      function E_from_p(p, spec, units)
!     this routine calculates kinetic energy from momentum for particle species SPEC
!     in unit system UNITS (c.f. function P_FROM_E), according to the formula
!
!        E = E0*( sqrt(1 + (p/(m0*c))**2) - 1 ) = sqrt(E0**2 + (p*c)**2) - E0.  (2)
!
      real(8), intent(in) :: p
      character(1), intent(in), optional :: spec, units
      real(8) :: E_from_p
      character(1) :: ispec, iunits
      real(8) :: c, E0

      if( present(spec) ) then
         ispec = spec
         if( ispec/='e' .and. ispec/='p' ) then
            print *,' function E_from_p:'
            print *,'  invalid input SPEC = ',ispec
            stop
         end if
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
         if( iunits/='M' .and. iunits/='C' .and. iunits/='N' ) then
            print *,' function E_from_p:'
            print *,'  invalid input UNITS = ',iunits
            stop
         end if
      else
         iunits = 'N'
      end if

      select case( iunits )
      case( 'N' )
         c = c_nat
         select case( ispec )
         case( 'e' )
            E0 = E0e_nat
         case( 'p' )
            E0 = E0p_nat
         end select
      case( 'M' )
         c = c_mks
         select case( ispec )
         case( 'e' )
            E0 = E0e_mks
         case( 'p' )
            E0 = E0p_mks
         end select
      case( 'C' )
         c = c_cgs
         select case( ispec )
         case( 'e' )
            E0 = E0e_cgs
         case( 'p' )
            E0 = E0p_cgs
         end select
      end select

      if( p>=0.d0 ) then
         E_from_p = sqrt(E0**2 + (p*c)**2) - E0
      else
         print *,' function E_from_p:'
         print *,'  input argument p (=',p,') is out of domain [0, +inf).'
         stop
      end if

      end function E_from_p



      function beta_sr(E, spec, units)
!     this routine calculates beta in special relativity for particle species SPEC
!     in unit system UNITS (c.f. function P_FROM_E) according to the formula
!
!        beta = v/c = sqrt(E'*(E' + 2))/(E' + 1),                               (3)
!
!     where E' = E/E0.
      real(8), intent(in) :: E
      character(1), intent(in), optional :: spec, units
      real(8) :: beta_sr
      character(1) :: ispec, iunits
      real(8) :: E0, E_prime

      if( present(spec) ) then
         ispec = spec
         if( ispec/='e' .and. ispec/='p' ) then
            print *,' function beta_sr:'
            print *,'  invalid input SPEC = ',ispec
            stop
         end if
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
         if( iunits/='M' .and. iunits/='C' .and. iunits/='N' ) then
            print *,' function beta_sr:'
            print *,'  invalid input UNITS = ',iunits
            stop
         end if
      else
         iunits = 'N'
      end if

      select case( iunits )
      case( 'N' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_nat
         case( 'p' )
            E0 = E0p_nat
         end select
      case( 'M' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_mks
         case( 'p' )
            E0 = E0p_mks
         end select
      case( 'C' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_cgs
         case( 'p' )
            E0 = E0p_cgs
         end select
      end select

      if( E>=0.d0 ) then
         E_prime = E/E0
         beta_sr = sqrt(E_prime*(E_prime + 2.d0))/(E_prime + 1.d0)
      else
         print *,' function beta_sr:'
         print *,'  input argument E (=',E,') is out of domain [0, +inf).'
         stop
      end if

      end function beta_sr



      function gamma_sr(E, spec, units)
!     this routine calculates gamma in special relativity for particle species SPEC
!     in unit system UNITS (c.f. function P_FROM_E) according to the formula
!
!        gamma = 1/sqrt(1 - beta**2) = 1 + E',                                  (4)
!
!     where E' = E/E0.
      real(8), intent(in) :: E
      character(1), intent(in), optional :: spec, units
      real(8) :: gamma_sr
      character(1) :: ispec, iunits
      real(8) :: E0, E_prime

      if( present(spec) ) then
         ispec = spec
         if( ispec/='e' .and. ispec/='p' ) then
            print *,' function gamma_sr:'
            print *,'  invalid input SPEC = ',ispec
            stop
         end if
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
         if( iunits/='M' .and. iunits/='C' .and. iunits/='N' ) then
            print *,' function gamma_sr:'
            print *,'  invalid input UNITS = ',iunits
            stop
         end if
      else
         iunits = 'N'
      end if

      select case( iunits )
      case( 'N' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_nat
         case( 'p' )
            E0 = E0p_nat
         end select
      case( 'M' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_mks
         case( 'p' )
            E0 = E0p_mks
         end select
      case( 'C' )
         select case( ispec )
         case( 'e' )
            E0 = E0e_cgs
         case( 'p' )
            E0 = E0p_cgs
         end select
      end select

      if( E>=0.d0 ) then
         E_prime = E/E0
         gamma_sr = 1.d0 + E_prime
      else
         print *,' function gamma_sr:'
         print *,'  input argument E (=',E,') is out of domain [0, +inf).'
         stop
      end if

      end function gamma_sr



      function Tbounce(E, y, L, spec, units)
!     particle bounce period (s) in Earth's dipole field calculated according to
!     Eq. (154) in [1]
!
!        Tb = (4*m*L*RE/p)*T(y) = (4*L*RE/(beta*c))*T(y).                       (5)
!
!     note that 4*L*RE*T(y) = oint{ ds/cos(a) }, where a is local pitch angle. input
!     arguments SPEC and UNITS are as explained in function P_FROM_E.
      real(8), intent(in) :: E, y, L
      character(1), intent(in), optional :: spec, units
      real(8) :: Tbounce
!     it doesn't matter whether RE and c are evaluated in MKS or CGS units here, as
!     long as they are consistent.
      real(8), parameter :: Const = 4.d0*RE_cgs/c_cgs
      character(1) :: ispec, iunits

      if( present(spec) ) then
         ispec = spec
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
      else
         iunits = 'N'
      end if

      if( L<=0.d0 ) then
         print *,' function Tbounce:'
         print *,'  input argument L (=',L,') is out of domain (0, +inf).'
         stop
      end if

      Tbounce = Const * L * SchulzT(y)/beta_sr(E, ispec, iunits)

      end function Tbounce



      function Omega2(E, y, L, spec, units)
!     particle bounce frequency (rad/s) in Earth's dipole field
      real(8), intent(in) :: E, y, L
      character(1), intent(in), optional :: spec, units
      real(8) :: Omega2
      character(1) :: ispec, iunits

      if( present(spec) ) then
         ispec = spec
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
      else
         iunits = 'N'
      end if

      Omega2 = 2.d0*Pi/Tbounce(E, y, L, ispec, iunits)

      end function Omega2



      function Omega3(E, y, L, spec, units)
!     particle drift frequency (rad/s) in Earth's dipole field calculated according
!     to Eq. (171) in [1]
!
!        Omega3 = (3*c*L*RE/(e*muE)) * (p**2/m) * D(y)/T(y)
!               = (3*c*L*RE/(e*muE)) * E0*beta**2*gamma * D(y)/T(y).            (6)
!
!     input arguments SPEC and UNITS are as explained in function P_FROM_E. N.B.
!     Eq. (6) is derived in CGS units; it would involve rationalization constants if
!     converted into SI units.
      real(8), intent(in) :: E, y, L
      character(1), intent(in), optional :: spec, units
      real(8) :: Omega3
      real(8), parameter :: Const = 3.d0*c_cgs*RE_cgs/(e_cgs*muE_cgs)
      real(8) :: E0_cgs
      character(1) :: ispec, iunits

      if( present(spec) ) then
         ispec = spec
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
      else
         iunits = 'N'
      end if

      if( L<=0.d0 ) then
         print *,' function Omega3:'
         print *,'  input argument L (=',L,') is out of domain (0, +inf).'
         stop
      end if

      select case( ispec )
      case( 'e' )
         E0_cgs = E0e_cgs
      case( 'p' )
         E0_cgs = E0p_cgs
      case default
         print *,' function Omega3:'
         print *,'  invalid input SPEC = ',ispec
         stop
      end select

      Omega3 = Const * L * E0_cgs * &
             & (beta_sr(E, ispec, iunits))**2 * gamma_sr(E, ispec, iunits) * &
             & SchulzD(y)/SchulzT(y)

      end function Omega3



      function Tdrift(E, y, L, spec, units)
!     particle drift period (s) in Earth's dipole field
      real(8), intent(in) :: E, y, L
      character(1), intent(in), optional :: spec, units
      real(8) :: Tdrift
      character(1) :: ispec, iunits

      if( present(spec) ) then
         ispec = spec
      else
         ispec = 'e'
      end if

      if( present(units) ) then
         iunits = units
      else
         iunits = 'N'
      end if

      Tdrift = 2.d0*Pi/Omega3(E, y, L, ispec, iunits)

      end function Tdrift
    
    
    
      subroutine jd2ut(jd,ny,nm,nd,ut)
!     This routine converts Julian Day number to Universal Time. Julian Day is the 
!     number of days since noon of 1 Jan 4713 B.C. This routine is adapted from
!     Mike Henderson's C function Lgm_JD_to_Date from LANL.
      implicit none
!     Julidan Day number
      real(8) :: jd
!     date and time in UT
      integer, intent(out) :: ny, nm, nd
      real(8), intent(out) :: ut
      integer :: I, A, B, C, D, E, G
      real(8) :: F
      
      jd = jd + 0.5d0
      
      I = int(jd)
      F = jd - dble(I)
      
      if( I>2299160 ) then
         A = int( (dble(I) - 1867216.25d0)/36524.25d0 )
         B = I + 1 + A - int( dble(A)/4.d0 )
      else
         B = I
      end if
      
      C = B + 1524
      D = int( (dble(C) - 122.1d0)/365.25d0 )
      E = int( 365.25d0*dble(D) )
      G = int( (dble(C) - dble(E))/30.6001d0 )
      
      nd = C - E - int( 30.6001d0*dble(G) )
      ut = 24.d0*F
      
      if( G<=13 ) then
         nm = G - 1
      else
         nm = G - 13
      end if
      
      if( nm>2 ) then
         ny = D - 4716
      else
         ny = D - 4715
      end if
      
      return
      end subroutine jd2ut
      
      
      
      subroutine ut2jd(ny,nm,nd,hh,mm,ss,jd)
!     This routine converts Universal Time to Julian Day number. This routine is
!     modified from function julday in "Numerical Recipes in FORTRAN", W. H. Press
!     et al., Cambridge University Press, 1992. Positive year signifies A.D.; 
!     negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      implicit none
!     integer number of year (ny), month (nm), day (nd), hour (hh) and minute (mm)
      integer, intent(in) :: ny, nm, nd, hh, mm
!     second (ss) is expressed in double precision, so that this routine has a
!     precision to subsecond.
      real(8), intent(in) :: ss
      real(8), intent(out) :: jd
!     Gregorian Calendar adopted 15 Oct. 1582
      integer, parameter :: IGREG = 15 + 31*(10+12*1582)
      integer :: ja, jm, jy, julday

!     calculate the integeral day number      
      jy = ny
      if( jy==0 ) then
         print *,'Error: there is no year zero'
         stop ' subroutine ut2jd stopped'
      end if
      if( jy<0 ) jy = jy + 1
      if( nm>2 ) then
         jm = nm + 1
      else
         jy = jy - 1
         jm = nm + 13
      end if
      julday = int(365.25*jy) + int(30.6001*jm) + nd + 1720995
!     test whether to change to Gregorian Calendar
      if( (nd+31*(nm+12*ny))>IGREG ) then
         ja = int(0.01*jy)
         julday = julday + 2 - ja + int(0.25*ja)
      end if
      
!     add to julday the decimals from hour, minute and second    
      jd = dble(julday) + (dble(hh)/24.d0 - 0.5d0) + dble(mm)/1440.d0 + ss/86400.d0
      
      return
      end subroutine ut2jd
    
    
    
      function time_diff(t0, t1)
!     This function returns the time difference between t0 and t1 in unit of day.
!     The format of input times is ISO 8601: YYYY-MM-DDThh:mm:ss, with accuracy to 1
!     second.
      integer, parameter :: Word_length = 19
      character(Word_length), intent(in) :: t0, t1
      real(8) :: time_diff
      integer :: year0, month0, day0, hour0, minute0, isecond0
      integer :: year1, month1, day1, hour1, minute1, isecond1
      real(8) :: second0, second1
!     Julian Day numbers of t0 and t1
      real(8) :: jd0, jd1
        
      time_diff = 0.d0
      if( LEN(t0)/=Word_length .or. LEN(t1)/=Word_length ) then
         print *,' function time_diff:'
         print *,'  invalid input time string format.'
         stop
      end if
        
      read(t0(1:4),'(I4)') year0
      read(t0(6:7),'(I2)') month0
      read(t0(9:10),'(I2)') day0
      read(t0(12:13),'(I2)') hour0
      read(t0(15:16),'(I2)') minute0
      read(t0(18:19),'(I2)') isecond0
      second0 = dble(isecond0)
        
      read(t1(1:4),'(I4)') year1
      read(t1(6:7),'(I2)') month1
      read(t1(9:10),'(I2)') day1
      read(t1(12:13),'(I2)') hour1
      read(t1(15:16),'(I2)') minute1
      read(t1(18:19),'(I2)') isecond1
      second1 = dble(isecond1)
        
      call ut2jd(year0,month0,day0,hour0,minute0,second0,jd0)
      call ut2jd(year1,month1,day1,hour1,minute1,second1,jd1)
      time_diff = jd1 - jd0
        
      end function time_diff

end module functions



!#define _TEST_
#ifdef _TEST_
! test for function Omega3 under different unit systems
      program test
      use functions
      implicit none
      real(8) :: L = 3.5d0
      real(8) :: a0(9), E, y
      integer :: i
      data a0 /10.d0, 20.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0/

100   format(F5.1,2X,ES14.6)

      E = E0e_mks
      print *,'MKS system'
      print '(A,ES14.6,A)','E =',E,' (J)'
      print '(A)','a0 (deg) Omega3 (rad/s)'
      do i=1,9
         y = sin(a0(i)*Po180)
         write(6,fmt=100) a0(i), Omega3(E, y, L, units='M')
      end do
      print *,''

      E = E0e_cgs
      print *,'CGS system'
      print '(A,ES14.6,A)','E =',E,' (erg)'
      print '(A)','a0 (deg) Omega3 (rad/s)'
      do i=1,9
         y = sin(a0(i)*Po180)
         write(6,fmt=100) a0(i), Omega3(E, y, L, units='C')
      end do
      print *,''

      E = E0e_nat
      print *,'Natural system'
      print '(A,ES14.6,A)','E =',E,' (MeV)'
      print '(A)','a0 (deg) Omega3 (rad/s)'
      do i=1,9
         y = sin(a0(i)*Po180)
         write(6,fmt=100) a0(i), Omega3(E, y, L, units='N')
      end do

      end program test
#endif

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! rewritten by Liheng Zheng on 03/25/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



