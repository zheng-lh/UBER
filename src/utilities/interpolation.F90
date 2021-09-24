module interpolation
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
      use fornberg
      implicit none
      public
      save
! This module contains the routines to perform interpolation and root finding.
! All routines in this module are public.
      
!     interface to the procedure pointers passed to the bisection root-finding
!     routines. the use of procedure pointers gives them more flexibility than taking
!     external procedures as argument. when parallilzing with OpenMP, the procedure
!     pointers passed to these root-finding routines must be declared threadprivate.
      abstract interface
!        a real function with a real argument
         real(8) function func_real(x)
            implicit none
            real(8), intent(in) :: x
         end function func_real
!        a real function with an integer argument
         real(8) function func_int(j)
            implicit none
            integer, intent(in) :: j
         end function func_int
      end interface

      interface hermint
         module procedure hermint1d
         module procedure hermint2d
         module procedure hermint3d
      end interface hermint
            
contains

      subroutine locate(xx,n,x,j)
!     A bisection index locator copied from "Numerical Recipes in FORTRAN",
!     W. H. Press et al., Cambridge University Press, 1992.    
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: j
      real(8), intent(in) :: x, xx(n)
!     Given an array xx(1:n), and given a value x, returns a value j such that x is
!     between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or
!     decreasing. j=0 or j=n is returned to indicate that x is out of range. If
!     x = xx(n), j=n is returned; if x = xx(1), j=1 is returned.
      integer :: jl, jm, ju
      
!     initialize lower and upper limits
      jl = 0
      ju = n + 1
      
!     if we are not yet done, compute a midpoint, and replace either the lower or
!     the upper limit, as appropriate. repeat until the test condition 10 is
!     satisfied.
10    if( (ju-jl)>1 ) then
         jm = (ju+jl)/2
!        the following condition has been changed from NR's original code. in the
!        original, if x happens to equal to some xx element, this loop gives the correct
!        location only when xx is decreasing. the modified condition avoids this problem.
         if(((xx(n)>xx(1)).and.(x>=xx(jm))).or.((xx(n)<xx(1)).and.(x<=xx(jm)))) then
            jl = jm
         else
            ju = jm
         end if
         goto 10
      end if
      
!     then set the output
      j = jl
      
      return
      end subroutine locate
      
      

      subroutine polint(xa,ya,n,x,y,dy)
!     Polynomial interpolation/extrapolation code adopted from "Numerical Recipes
!     in FORTRAN", W. H. Press et al., Cambridge University Press, 1992.
      implicit none
      integer :: n, NMAX
      real(8) :: dy, x, y, xa(n), ya(n)
      parameter (NMAX=10)
!     Given arrays xa and ya, each of length n, and given a value x, this routine returns
!     a value y, and an error estimate dy. if P(x) is the polynomial of degree N-1 such
!     that P(xa_i) = ya_i, i=1,...,n, then the returned value y=P(x).
      integer :: i, m, ns
      real(8) :: den, dif, dift, ho, hp, w, c(NMAX), d(NMAX)
      ns = 1
      dif = abs(x-xa(1))
      do i=1,n
         dift = abs(x-xa(i))
         if( dift<dif ) then
            ns = i
            dif = dift
         end if
         c(i) = ya(i)
         d(i) = ya(i)
      end do
      y = ya(ns)
      ns = ns - 1
      do m=1,n-1
         do i=1,n-m
            ho = xa(i) - x
            hp = xa(i+m) - x
            w = c(i+1) - d(i)
            den = ho - hp
            if( den==0. ) then
!              this error can occur only if two input xa's are (to within roundoff)
!              identical
               print *,' subroutine polint failure:'
               print *,' xa(n) =',xa(:)
               print *,' ya(n) =',ya(:)
               stop
            end if
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         end do
         if( 2*ns < n-m ) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns - 1
         end if
         y = y + dy
      end do
      return
      end subroutine polint
      
!     * * * * * * * * * * * * * * * * * * * * * * *
!
!     Bi-linear and Bi-cubic Interpolation
!
!     * * * * * * * * * * * * * * * * * * * * * * *
      
      subroutine linint(y,x1l,x1u,x1,ansy)
!     Linear interpolation
      implicit none
      real(8), intent(in) :: x1, x1l, x1u
      real(8), dimension(2), intent(in) :: y
      real(8), intent(out) :: ansy
      real(8) :: t

      if( x1u==x1l ) then
         print *,' Bad input in linint'
         stop ' subroutine linint stopped'
      end if

      t = (x1-x1l)/(x1u-x1l)
      ansy = (1.d0 - t)*y(1) + t*y(2)

      return
      end subroutine linint


      
      subroutine blnint(y,x1l,x1u,x2l,x2u,x1,x2,ansy)
!     Bilinear interpolation routine.
      implicit none
      real(8), intent(in) :: x1, x1l, x1u, x2, x2l, x2u
      real(8), dimension(4), intent(in) :: y
      real(8), intent(out) :: ansy
!     Bilinear interpolation within a grid square. Input quantities are y
!     (as described in blnint3d); x1l and x1u, the lower and upper coordinates of the
!     grid square in the 1-direction; x2l and x2u likewise for the 2-direction; and x1,
!     x2, the coordinates of the desired point for the interpolation. The interpolated 
!     function value is returned as ansy.
      real(8) :: t, u
       
      if( x1u==x1l .or. x2u==x2l ) then
         print *,' Bad input in blnint'
         stop ' subroutine blnint stopped'
      end if
      
      t = (x1-x1l)/(x1u-x1l)
      u = (x2-x2l)/(x2u-x2l)
      
      ansy = (1.0D0-t)*(1.0D0-u)*y(1) + t*(1.0D0-u)*y(2) + t*u*y(3) + (1.0D0-t)*u*y(4)
      
      return
      end subroutine blnint
      
      
      
      subroutine blnint3d(y,x1l,x1u,x2l,x2u,x3l,x3u,x1,x2,x3,ansy)
!     3D bilinear interpolation routine.
      implicit none
      real(8), intent(in) :: x1, x1l, x1u, x2, x2l, x2u, x3, x3l, x3u
!     from the lowest 1, 2 and 3 coordinates, y counts counterclockwisely in the x3l
!     plane, and then does the same way in the x3u plane. thus y(1) starts from
!     (x1l,x2l,x3l) and y(8) ends at (x1l,x2u,x3u)
      real(8), intent(in) :: y(8)
      real(8), intent(out) :: ansy
!     Bilinear interpolation within a grid cube. Input quantities are the grid node
!     values y(8); xil and xiu, the lower and upper coordinates of the grid cube in 
!     the i-direction; and (x1,x2,x3), the coordinates of the desired point for the
!     interpolation. The interpolated function value is returned as ansy.
      real(8) :: t, u, v
       
      if( x1u==x1l .or. x2u==x2l .or. x3u==x3l ) then
         print *,' Bad input in blnint3d'
         stop ' subroutine blnint3d stopped'
      end if
      
      t = (x1-x1l)/(x1u-x1l)
      u = (x2-x2l)/(x2u-x2l)
      v = (x3-x3l)/(x3u-x3l)
      
      ansy = (1.D0-t)*(1.D0-u)*(1.D0-v)*y(1) + &
                   t *(1.D0-u)*(1.D0-v)*y(2) + &
                   t *      u *(1.D0-v)*y(3) + &
             (1.D0-t)*      u *(1.D0-v)*y(4) + &
             (1.D0-t)*(1.D0-u)*      v *y(5) + &
                   t *(1.D0-u)*      v *y(6) + &
                   t *      u *      v *y(7) + &
             (1.D0-t)*      u *      v *y(8)
      
      return
      end subroutine blnint3d



      subroutine loglinint(y,x1l,x1u,x1,ansy)
!     Log linear interpolation in 1D.
      implicit none
      real(8), intent(in) :: x1l, x1u, x1
      real(8), dimension(2), intent(in) :: y
      real(8), intent(out) :: ansy
      real(8), dimension(2) :: logy

      if( ANY(y<=0.d0) ) then
         if( ANY(y>=0.d0) ) then
            call linint(y,x1l,x1u,x1,ansy)
         else
            logy = log(-y)
            call linint(logy,x1l,x1u,x1,ansy)
            ansy = -exp(ansy)
         end if
      else
         logy = log(y)
         call linint(logy,x1l,x1u,x1,ansy)
         ansy = exp(ansy)
      end if

      return
      end subroutine loglinint
      
      
      
      subroutine logblnint(y,x1l,x1u,x2l,x2u,x1,x2,ansy)
!     Log bilinear interpolation in 2D. This routine first takes natural log of the
!     grid cube values y, then bilinearly interpolates log(y), and returns exponential
!     of the interpolated value. If any element in y is less than or equal to zero,
!     a bilinear interpolation is instead performed. Interpolation in log scale is
!     more appropriate if data are positive and vary by orders of magnitude.
      implicit none
      real(8), intent(in) :: x1, x1l, x1u, x2, x2l, x2u
      real(8), dimension(4), intent(in) :: y
      real(8), intent(out) :: ansy
      real(8), dimension(4) :: logy
      
      if( ANY(y<=0.d0) ) then
         if( ANY(y>=0.d0) ) then
!           the entries of y change sign or vanish
            call blnint(y,x1l,x1u,x2l,x2u,x1,x2,ansy)
         else
!           all entries of y are less than zero
            logy = log(-y)
            call blnint(logy,x1l,x1u,x2l,x2u,x1,x2,ansy)
            ansy = -exp(ansy)
         end if
      else
         logy = log(y)
         call blnint(logy,x1l,x1u,x2l,x2u,x1,x2,ansy)
         ansy = exp(ansy)
      end if
      
      return
      end subroutine logblnint
      
      
      
      subroutine logblnint3d(y,x1l,x1u,x2l,x2u,x3l,x3u,x1,x2,x3,ansy)
!     Log bilinear interpolation in 3D.
      implicit none
      real(8), intent(in) :: x1, x1l, x1u, x2, x2l, x2u, x3, x3l, x3u
      real(8), intent(in) :: y(8)
      real(8), intent(out) :: ansy
      real(8) :: logy(8)
      
      if( ANY(y<=0.d0) ) then
         if( ANY(y>=0.d0) ) then
!           the entries of y change sign or vanish
            call blnint3d(y,x1l,x1u,x2l,x2u,x3l,x3u,x1,x2,x3,ansy)
         else
!           all entries of y are less than zero
            logy = log(-y)
            call blnint3d(logy,x1l,x1u,x2l,x2u,x3l,x3u,x1,x2,x3,ansy)
            ansy = -exp(ansy)
         end if
      else
         logy = log(y)
         call blnint3d(logy,x1l,x1u,x2l,x2u,x3l,x3u,x1,x2,x3,ansy)
         ansy = exp(ansy)
      end if
      
      return
      end subroutine logblnint3d


      
      subroutine bcucof(y,y1,y2,y12,d1,d2,c)
!     Procedure to calculate the coefficient matrix used in subroutine bcuint. Copied 
!     from "Numerical Recipes in FORTRAN", W. H. Press et al., Cambridge University 
!     Press, 1992.
      implicit none
      real(8), intent(in) :: d1, d2
      real(8), dimension(4), intent(in) :: y, y1, y2, y12
      real(8), dimension(4,4), intent(out) :: c
!     Given arrays y, y1, y2 and y12, each of length 4, containing the function,
!     gradients, and cross derivative at the four grid points of a rectangular grid
!     cell (numbered counterclockwise from the lower left), and given d1 and d2, 
!     the length of the grid cell in the 1- and 2- directions, this routine returns the
!     table c(1:4,1:4) that is used by routine bcuint for bicubic interpolation.
      integer :: i, j, k, l
      real(8) :: d1d2, xx, cl(16), wt(16,16), x(16)
      save wt
      data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4, &
               10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4, &
               4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2, &
               10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2, &
               0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2, &
               10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2, &
               5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1, &
               10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
               
      d1d2 = d1*d2
      
!     pack a temporary vector x
      do i=1,4
         x(i) = y(i)
         x(i+4) = y1(i)*d1
         x(i+8) = y2(i)*d2
         x(i+12) = y12(i)*d1d2
      end do
      
!     matrix multiply by the stored table
      do i=1,16
         xx = 0.0
         do k=1,16
            xx = xx + wt(i,k)*x(k)
         end do
         cl(i) = xx
      end do
      
!     unpack the result into the output table
      l = 0
      do i=1,4
         do j=1,4
            l = l + 1
            c(i,j) = cl(l)
         end do
      end do
      
      return
      end subroutine bcucof
          


      subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
!     Bicubic interpolation routine copied from "Numerical Recipes in FORTRAN",
!     W. H. Press et al., Cambridge University Press, 1992.
      implicit none
      real(8), intent(in) :: x1, x1l, x1u, x2, x2l, x2u
      real(8), dimension(4), intent(in) :: y, y1, y12, y2
      real(8), intent(out) :: ansy, ansy1, ansy2
!     Bicubic interpolation within a grid square. Input quantities are y, y1, y2, y12
!     (as described in bcucof); x1l and x1u, the lower and upper coordinates of the
!     grid square in the 1-direction; x2l and x2u likewise for the 2-direction; and x1,
!     x2, the coordinates of the desired point for the interpolation. The interpolated 
!     function value is returned as ansy, and the interpolated gradient values as ansy1
!     and ansy2. This routine calls bcucof.
      integer :: i
      real(8) :: t, u, c(4,4)
      
!     get the c's
      call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
      
      if( x1u==x1l .or. x2u==x2l ) then
         print *,' Bad input in bcuint'
         stop
      end if
      
      t = (x1-x1l)/(x1u-x1l)
      u = (x2-x2l)/(x2u-x2l)
      
      ansy = 0.0
      ansy1 = 0.0
      ansy2 = 0.0
      
      do i=4,1,-1
         ansy = t*ansy + ( (c(i,4)*u + c(i,3))*u + c(i,2) )*u + c(i,1)
         ansy2 = t*ansy2 + ( 3.0D0*c(i,4)*u + 2.0D0*c(i,3) )*u + c(i,2)
         ansy1 = u*ansy1 + ( 3.0D0*c(4,i)*t + 2.0D0*c(3,i) )*t + c(2,i)
      end do
      
      ansy1 = ansy1/(x1u-x1l)
      ansy2 = ansy2/(x2u-x2l)
      
      return
      end subroutine bcuint
      
      

      subroutine packvec(y,i,j,array)
!     Auxiliary routine used by interp.
      real(8), intent(out) :: y(4)
      integer, intent(in) :: i, j
      real(8), intent(in) :: array(:,:)
      
      y = 0.d0    
      y(1) = array(i,j)
      y(2) = array(i+1,j)
      y(3) = array(i+1,j+1)
      y(4) = array(i,j+1)
      
      return
      end subroutine packvec
      
      
      
      subroutine packvec3d(y,i,j,k,array)
!     axuiliary routine used by interp (3D)
      real(8), intent(out) :: y(8)
      integer, intent(in) :: i, j, k
      real(8), intent(in) :: array(:,:,:)
      
      y = 0.D0
      y(1) = array(i,j,k)
      y(2) = array(i+1,j,k)
      y(3) = array(i+1,j+1,k)
      y(4) = array(i,j+1,k)
      y(5) = array(i,j,k+1)
      y(6) = array(i+1,j,k+1)
      y(7) = array(i+1,j+1,k+1)
      y(8) = array(i,j+1,k+1)
      
      return
      end subroutine packvec3d

!     * * * * * * * * * * * * * * * * * * *
!
!     Monotonic Cubic Hermite Polynomial Interpolation
!
!     * * * * * * * * * * * * * * * * * * *

      subroutine monotone_tangent(xx,n,yy,tt)
!     This subroutine calculates the tangent (derivative) values tt(1:n) at each
!     (xx(i),yy(i)) for i=1,...,n used in monotonic cubic Hermite interpolation.
!     This method is based on
!        Fritsch and Carlson (1980), "Monotone piecewise cubic interpolation",
!        SIAM Journal on Numerical Analysis (SIAM) 17 (2): 238-246,
!        doi:10.1137/0717021
!     the wikipedia page
!        <https://en.wikipedia.org/wiki/Monotone_cubic_interpolation>
!     and Lunjin Chen's fortran subroutine hermit_g2.
      implicit none
      integer :: n
!     xx must be monotonically increasing and distinct
      real(8) :: xx(n), yy(n), tt(n)
      integer :: i
      real(8) :: eps, dx, r2, alpha, beta, tau
      real(8) :: delta(0:n-1)

      if( n<3 ) then
         print *,'Error: input arrays must have at least length 3.'
         stop ' Subroutine monotone_tangent stopped.'
      end if
      eps = 1.d-9*(maxval(yy) - minval(yy))

!     initialize tt
      delta(0) = 0.d0
      do i=1,n-1
         dx = xx(i+1) - xx(i)
         if( dx<eps ) then
            print *,' Error: input array xx must be monotonically increasing and distinct.'
            print *,' xx =',xx
            stop ' Subroutine monotone_tangent stopped.'
         end if
         delta(i) = (yy(i+1) - yy(i))/dx
         tt(i) = 0.5d0*(delta(i-1) + delta(i))
      end do

      tt(1) = delta(1)
      tt(n) = delta(n-1)

!     modify tt to meet the necessary conditions for monotonicity
      do i=1,n-1
!        step 1
         if( delta(i-1)*delta(i)<0.d0 .or. tt(i)*delta(i)<0.d0 .or. &
           & tt(i)*delta(i-1)<0.d0 ) then
            tt(i) = 0.d0
         end if
         if( abs(delta(i))<eps ) then
            tt(i) = 0.d0
            tt(i+1) = 0.d0
!        step 2
         else
            alpha = tt(i)/delta(i)
            beta = tt(i+1)/delta(i)
            r2 = alpha**2 + beta**2
            if( r2>9.d0 ) then
               tau = 3.d0/sqrt(r2)
               tt(i) = tau*alpha*delta(i)
               tt(i+1) = tau*beta*delta(i)
            end if
         end if
      end do

      return
      end subroutine monotone_tangent


!     Hermite polynomials and their first derivatives used in cubic Hermite interpolation
!     In each of these subroutines, x (in [0,1]) is the independent variable, y =
!     H(x) gives the value of the cubic Hermite polynomial, and dydx = dy/dx.
      subroutine hmp0(x,y,dydx)
      implicit none
      real(8), intent(in) :: x
      real(8), intent(out) :: y, dydx

      y = (1.d0 + 2.d0*x) * (1.d0 - x)**2
      dydx = 6.d0 * x * (x - 1.d0)

      return
      end subroutine hmp0



      subroutine hmp1(x,y,dydx)
      implicit none
      real(8), intent(in) :: x
      real(8), intent(out) :: y, dydx

      y = x * (1.d0 - x)**2
      dydx = (x - 1.d0) * (3.d0*x - 1.d0)

      return
      end subroutine hmp1

!     The following cubic Hermite interpolations are based on:
!     [1] Phillips, G. M., Explicit forms for certain Hermite approximations,
!         BIT 13 (1973), 177-180
!     [2] <https://en.wikipedia.org/wiki/Cubic_Hermite_spline>
!     [3] Lunjin Chen's fortran code in c40glf.f

      subroutine hermint1d(xx,yy,tt,n,x,y,dydx)
!     One dimensional cubic Hermite interpolation
      implicit none
!     Given data arrays (xx(i),yy(i)) and the tangent at each data point tt(i)
!     (i=1,...,n), this subroutine returns the interpolated value y and derivative
!     dydx at the given position x, so that the interpolating function satisfies
!     y(xx(i)) = yy(i) and dydx(xx(i)) = tt(i).
      integer, intent(in) :: n
      real(8), intent(in) :: xx(n), yy(n), tt(n)
      real(8), intent(in) :: x
      real(8), intent(out) :: y, dydx
      integer :: i
      real(8) :: t, h
      real(8) :: y0, y1, dydx0, dydx1, q0, q1, dqdx0, dqdx1, p0, p1, dpdx0, dpdx1

      if( x<xx(1) .or. x>xx(n) ) then
         print *,'Error: the interpolation point is out of the grid.'
         print *,'x =',x,', xx(1) =',xx(1),', xx(n) =',xx(n)
         stop ' Subroutine hermint1d stopped.'
      end if

      call locate(xx, n, x, i)
!     if x = xx(n), no interpolation needed
      if( i==n ) then
         y = yy(i)
         dydx = tt(i)
         return
      end if

      h = xx(i+1) - xx(i)
      t = (x - xx(i))/h
      y0 = yy(i)
      y1 = yy(i+1)
      dydx0 = tt(i)*h
      dydx1 = tt(i+1)*h

      call hmp0(t, q0, dqdx0)
      call hmp1(t, q1, dqdx1)
      call hmp0(1.d0-t, p0, dpdx0)
      call hmp1(1.d0-t, p1, dpdx1)

      y = y0*q0 + y1*p0 + dydx0*q1 - dydx1*p1
      dydx = y0*dqdx0 - y1*dpdx0 + dydx0*dqdx1 + dydx1*dpdx1
      dydx = dydx/h

      return
      end subroutine hermint1d



      subroutine hermint2d(xx,m,yy,n,zz,dzdx,dzdy,x,y,f,dfdx,dfdy)
!     Two dimensional cubic Hermite interpolation
      implicit none
!     Given coordinate arrays xx(m), yy(n), data array zz(m,n) on these
!     coordinates and its x-, y-derivatives dzdx(m,n) and dzdy(m,n), this subroutine
!     returns the interpolated value f and the two partial derivatives dfdx, dfdy at
!     the given position (x,y), so that the interpolating function satisfies that
!     f(xx(i),yy(j)) = zz(i,j), and dfdx(xx(i),yy(j)) = dzdx(i,j) and
!     dfdy(xx(i),yy(j)) = dzdy(i,j).
      integer, intent(in) :: m, n
      real(8), intent(in) :: xx(m), yy(n), zz(m,n), dzdx(m,n), dzdy(m,n)
      real(8), intent(in) :: x, y
      real(8), intent(out) :: f, dfdx, dfdy
      integer :: i, j
      real(8) :: tx, ty, hx, hy
      real(8) :: f00, f01, f10, f11
      real(8) :: dfdx00, dfdx01, dfdx10, dfdx11
      real(8) :: dfdy00, dfdy01, dfdy10, dfdy11
      real(8) :: qx0, qx1, dqdx0, dqdx1
      real(8) :: qy0, qy1, dqdy0, dqdy1
      real(8) :: px0, px1, dpdx0, dpdx1
      real(8) :: py0, py1, dpdy0, dpdy1

      if( x<xx(1) .or. x>xx(m) .or. y<yy(1) .or. y>yy(n) ) then
         print *,'Error: the interpolation point is out of the grid.'
         print *,'x =',x,', y =',y
         print *,'xx(1) =',xx(1),', xx(m) =',xx(m)
         print *,'yy(1) =',yy(1),', yy(n) =',yy(n)
         stop ' Subroutine hermint2d stopped.'
      end if

      call locate(xx, m, x, i)
      call locate(yy, n, y, j)
      if( i==m .and. j==n ) then
!        no need to interpolate in this case
         f = zz(i,j)
         dfdx = dzdx(i,j)
         dfdy = dzdy(i,j)
         return
      else
         if( i==m ) i = m - 1
         if( j==n ) j = n - 1
      end if

!     in a grid cell
!     y ^
!       |
!     - |  01   11
!    hy |     .
!     - |  00   10
!       +------------> x
!          |  hx |
      hx = xx(i+1) - xx(i)
      tx = (x - xx(i))/hx
      hy = yy(j+1) - yy(j)
      ty = (y - yy(j))/hy

      f00 = zz(i,j)
      f01 = zz(i,j+1)
      f10 = zz(i+1,j)
      f11 = zz(i+1,j+1)

      dfdx00 = dzdx(i,j)*hx
      dfdx01 = dzdx(i,j+1)*hx
      dfdx10 = dzdx(i+1,j)*hx
      dfdx11 = dzdx(i+1,j+1)*hx

      dfdy00 = dzdy(i,j)*hy
      dfdy01 = dzdy(i,j+1)*hy
      dfdy10 = dzdy(i+1,j)*hy
      dfdy11 = dzdy(i+1,j+1)*hy

      call hmp0(tx, qx0, dqdx0)
      call hmp0(ty, qy0, dqdy0)
      call hmp1(tx, qx1, dqdx1)
      call hmp1(ty, qy1, dqdy1)
      call hmp0(1.d0-tx, px0, dpdx0)
      call hmp0(1.d0-ty, py0, dpdy0)
      call hmp1(1.d0-tx, px1, dpdx1)
      call hmp1(1.d0-ty, py1, dpdy1)

      f = f00*qx0*qy0 + f01*qx0*py0 + f10*px0*qy0 + f11*px0*py0 &
      & + dfdx00*qx1*qy0 + dfdx01*qx1*py0 - dfdx10*px1*qy0 - dfdx11*px1*py0 &
      & + dfdy00*qx0*qy1 - dfdy01*qx0*py1 + dfdy10*px0*qy1 - dfdy11*px0*py1

      dfdx = f00*dqdx0*qy0 + f01*dqdx0*py0 - f10*dpdx0*qy0 - f11*dpdx0*py0 &
         & + dfdx00*dqdx1*qy0 + dfdx01*dqdx1*py0 &
         & + dfdx10*dpdx1*qy0 + dfdx11*dpdx1*py0 &
         & + dfdy00*dqdx0*qy1 - dfdy01*dqdx0*py1 &
         & - dfdy10*dpdx0*qy1 + dfdy11*dpdx0*py1

      dfdy = f00*qx0*dqdy0 - f01*qx0*dpdy0 + f10*px0*dqdy0 - f11*px0*dpdy0 &
         & + dfdx00*qx1*dqdy0 - dfdx01*qx1*dpdy0 &
         & - dfdx10*px1*dqdy0 + dfdx11*px1*dpdy0 &
         & + dfdy00*qx0*dqdy1 + dfdy01*qx0*dpdy1 &
         & + dfdy10*px0*dqdy1 + dfdy11*px0*dpdy1

      dfdx = dfdx/hx
      dfdy = dfdy/hy

      return
      end subroutine hermint2d



      subroutine hermint3d(xx,m,yy,n,zz,l,uu,dudx,dudy,dudz,x,y,z,f,dfdx,dfdy,dfdz)
!     Three dimensional cubic Hermite interpolation
      implicit none
      integer :: m, n, l
      real(8) :: xx(m), yy(n), zz(l)
      real(8) :: uu(m,n,l), dudx(m,n,l), dudy(m,n,l), dudz(m,n,l)
      real(8) :: x, y, z
      real(8) :: f, dfdx, dfdy, dfdz

      print *,' module interpolation procedure hermint3d is under construction'
      stop

      end subroutine hermint3d

!     * * * * * * * * * * * * * * * * * * * * * *
!
!     Root Finding
!
!     * * * * * * * * * * * * * * * * * * * * * *

      subroutine zbrac(func,x1,x2,xl,xu,succes)
!     A root bracketing routine adopted from "Numerical Recipes in FORTRAN",
!     W. H. Press et al., Cambridge University Press, 1992.    
      implicit none
      integer :: NTRY
!     xl, xu: lower and upper limit of x
      procedure(func_real), pointer :: func
      real(8) :: x1, x2, xl, xu, FACTOR
      parameter (FACTOR=1.6d0, NTRY=50)
!     Given a function func and an initial guessed range x1 to x2, the routine expands
!     the range geometrically until a root is bracketed by the returned values x1 and
!     x2 (in which case succes returns as .true.) or until the range becomes
!     unacceptably large (in which case succes returns as .false.).
      integer :: j
      real(8) :: f1, f2
      logical :: succes
      if( xl==xu ) stop ' xu has to be greater than xl in zbrac'
      if( (x1<xl.and.x2<xl) .or. (x1>xu.and.x2>xu) ) goto 100
      if( min(x1,x2)<xl ) x1 = xl
      if( max(x1,x2)>xu ) x2 = xu
      if( min(x1,x2)>=xu .or. max(x1,x2)<=xl ) goto 100
      if( x1==x2 ) stop ' you have to guess an initial range in zbrac'
      f1 = func(x1)
      f2 = func(x2)
      succes = .true.
      do j=1,NTRY
         if( f1*f2<0.D0 ) return
         if( abs(f1)<abs(f2) ) then
            if( x1==xl ) then
               succes = .false.
               return
            end if
            x1 = x1 + FACTOR*(x1 - x2)
            if( x1<xl ) x1 = xl
            f1 = func(x1)
         else
            if( x2==xu ) then
               succes = .false.
               return
            end if
            x2 = x2 + FACTOR*(x2 - x1)
            if( x2>xu ) x2 = xu
            f2 = func(x2)
         end if
      end do
100   succes = .false.
      return
      end subroutine zbrac
      
      
      
      real(8) function rtbis(func,x1,x2,xacc)
!     A bisection root finding routine adopted from "Numerical Recipes in FORTRAN",
!     W. H. Press et al., Cambridge University Press, 1992.    
      implicit none
!     maximum allowed number of bisections
      integer, parameter :: JMAX = 40
      procedure(func_real), pointer :: func
      real(8), intent(in) :: x1, x2, xacc
!     Using bisection, find the root of a function func known to lie between x1 and
!     x2. The root, returned as rtbis, will be refined until its accuracy is +/-xacc.
      integer :: j
      real(8) :: dx, f, fmid, xmid
      fmid = func(x2)
      f = func(x1)
      if( f*fmid>=0.D0 ) then
         print *,'Error in rtbis: root must be bracketed in rtbis'
         print *,'  x1 =',x1,',  x2 =',x2
         stop ' function rtbis stopped'
      end if
!     orient the search so that f>0 lies at x+dx
      if( f<0.D0 ) then
         rtbis = x1
         dx = x2 - x1
      else
         rtbis = x2
         dx = x1 - x2
      end if
!     bisection loop
      do j=1,JMAX
         dx = dx*0.5D0
         xmid = rtbis + dx
         fmid = func(xmid)
         if( fmid<=0.D0 ) rtbis = xmid
         if( abs(dx)<xacc .or. fmid==0.D0 ) return
      end do
      print '(A,I4,A)',' rtbis failed to find root in ',JMAX,' bisections'
      print *,' dx =',dx,',  xacc =',xacc
      print *,' root was asigned as the last bisection value'
      return
      end function rtbis
      
      
      
      integer function indbis(func,j1,j2)
!     This routine is an analogy of the function rtbis, except that rtbis finds the
!     root x of equation func(x) = 0, whereas this routine finds the index j that the
!     root of func(x(j)) = 0 lies within the indices j and (j + 1). There must be one
!     and only one root between [j1,j2].
      implicit none
      procedure(func_int), pointer :: func
      integer, intent(in) :: j1, j2
      integer :: jmid, i1, i2
      real(8) :: f1, fmid
      
      indbis = 0
      if( j2<=j1 .or. func(j1)*func(j2)>=0.d0 ) then
         print *,'indbis error:'
         print '(A,2(I4,A))',' one and only one root must be within indices [',j1,',',j2,']'
         stop ' function indbis stopped'
      end if
      
      i1 = j1
      i2 = j2
      do while( i2-i1>1 )
         jmid = i1 + (i2 - i1)/2
         f1 = func(i1)
         fmid = func(jmid)
         if( f1*fmid>=0.d0 ) then
            i1 = jmid
         else
            i2 = jmid
         end if
      end do      
      indbis = i1
      
      return
      end function indbis
      
      
            
!*    this is the subroutine version of the function rtbis. it is only used in debug.
!      subroutine subrtbis(func,x1,x2,xacc,ans,err)
!      implicit none
!*    maximum allowed number of bisections
!      integer, parameter :: JMAX = 40
!      procedure(func_real), pointer :: func
!      real(8), intent(in) :: x1, x2, xacc
!      real(8), intent(out) :: ans
!      integer, intent(out) :: err
!*    Using bisection, find the root of a function func known to lie between x1 and
!*    x2. The root, returned as rtbis, will be refined until its accuracy is +/-xacc.
!      integer :: j
!      real(8) :: dx, f, fmid, xmid
!      fmid = func(x2)
!      f = func(x1)
!      if( f*fmid>=0.D0 ) then
!         err = 1
!         print *,'Error in subrtbis: root must be bracketed within'
!         print *,'  x1 =',x1,',  x2 =',x2
!         return
!      end if
!      err = 0
!*    orient the search so that f>0 lies at x+dx
!      if( f<0.D0 ) then
!         ans = x1
!         dx = x2 - x1
!      else
!         ans = x2
!         dx = x1 - x2
!      end if
!*    bisection loop
!      do j=1,JMAX
!         dx = dx*0.5D0
!         xmid = ans + dx
!         fmid = func(xmid)
!         if( fmid<=0.D0 ) ans = xmid
!         if( abs(dx)<xacc .or. fmid==0.D0 ) return
!      end do
!      print '(A,I4,A)',' subrtbis failed to find root in ',JMAX,' bisections'
!      print *,' dx =',dx,',  xacc =',xacc
!      print *,' root was asigned as the last bisection value'
!      return
!      end subroutine subrtbis

!     * * * * * * * * * * * * * * * * * * * * *
!
!     Curve Fitting
!
!     * * * * * * * * * * * * * * * * * * * * *
      
      subroutine fit(x,y,ndata,a,b,siga,sigb,chi2,r)
!     Linear least-square regression routine adopted and modified from "Numerical
!     Recipes in FORTRAN", W. H. Press et al., Cambridge University Press, 1992.
      implicit none
!     Given a set of data points x(1:ndata), y(1:ndata) with assumed unit standard
!     deviations, this routine fits them to a straight line y = a + b*x by minimizing
!     the chi-square chi2. Returned are a, b and their respective probable uncertanties
!     siga and sigb, chi2, and the coefficient of correlation r.
      integer, intent(in) :: ndata
      real(8), intent(in) :: x(ndata), y(ndata)
      real(8), intent(out) :: a, b, siga, sigb, r, chi2
      integer :: i
      real(8) :: sigdat, ss, st2, sx, sxoss, sy, t
      
      if( ndata<=1 ) then
         print *,'Too few data to perform least-square regression'
         print *,'ndata =',ndata
         print *,'subroutine fit was not executed'
         return
      end if
      
      sx = 0.d0
      sy = 0.d0
      st2 = 0.d0
      b = 0.d0
      chi2 = 0.d0
      
      do i=1,ndata
         sx = sx + x(i)
         sy = sy + y(i)
      end do
      
      ss = DBLE(ndata)
      sxoss = sx/ss
      
      do i=1,ndata
         t = x(i) - sxoss
         st2 = st2 + t*t
         b = b + t*y(i)
      end do
      
      b = b/st2
      a = (sy - sx*b)/ss
      siga = sqrt((1.d0 + sx*sx/(ss*st2))/ss)
      sigb = sqrt(1.d0/st2)
      
      r = -sxoss/(st2*siga*sigb)
      
      do i=1,ndata
         chi2 = chi2 + (y(i) - a - b*x(i))**2
      end do
      
!     evaluate typical data sigma using chi2, and adjust the standard deviations
      sigdat = sqrt(chi2/DBLE(ndata-2))
      siga = siga*sigdat
      sigb = sigb*sigdat
      
      return
      end subroutine fit

end module interpolation
      

!#define _TEST_
#ifdef _TEST_
      program test
!     This is the monotonic cubic Hermite interpolation benchmark used in Fritsch and
!     Carlson (1980).
      use interpolation
      implicit none
      integer, parameter :: N = 11
      integer, parameter :: M = 1001
      integer :: i
      real(8) :: xx(N), yy(N), tt(N)
      real(8) :: x, y, t
      data xx /0.d0, 2.d0, 3.d0, 5.d0, 6.d0, 8.d0, 9.d0, 11.d0, 12.d0, 14.d0, 15.d0/
      data yy /10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.5d0, 15.d0, 50.d0, &
             & 60.d0, 85.d0/

      call monotone_tangent(xx,N,yy,tt)
      open(10,file='interpin.txt',status='replace',form='formatted')
      do i=1,N
         write(10,'(3(F16.9))') xx(i), yy(i), tt(i)
      end do
      close(10)

      open(11,file='interpout.txt',status='replace',form='formatted')
      do i=1,M
         x = dble(i-1)*15.d0/dble(M-1)
         call hermint1d(xx, yy, tt, N, x, y, t)
         write(11,'(3(F16.9))') x, y, t
      end do
      close(11)

      end program test

#endif
