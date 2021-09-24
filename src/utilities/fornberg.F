      module fornberg
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c Copyright 2021, Liheng Zheng
c
c This file is part of UBER.
c
c    UBER is free software: you can redistribute it and/or modify it under the
c    terms of the MIT License as published by Massachusetts Institute of
c    Technology. UBER is distributed in the hope that it will be useful, but
c    WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
c    FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
c
c    You should have received a copy of the MIT License along with UBER. If not,
c    see <https://opensource.org/licenses/MIT>.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c...  This module, named after Bengt Fornberg, the person who deviced these algorithms,
c...  contains the routines that calculate the first and second order derivatives of
c...  a function on a one-dimensional unevenly spaced grid.
c
c...  The routines were downloaded from Arizona State University website
c...  http://cococubed.asu.edu/code_pages/fdcoef.shtml
c
c...  They are based on Bengt Fornberg's articles "Generation of finite difference
c...  formulas on arbitrary spaced grids", Math. Comp., 51(184):699-706, 1988, and
c...  "Calculation of weights in finite difference formulas", SIAM Rev., 40(3):685b
c
c...  All possible rights belong to the original author.   
      contains 


      subroutine fss002(grid,ngrid,u,du)
      implicit none
      save

c..this routine computes second order accurate first derivatives 
c..on an arbitrary grid

c..input:
c..ngrid       = number of points in the grid, 
c..grid(ngrid) = array of independent values 
c..u(ngrid)    = function values at the grid points

c..output:
c..du(ngrid)   = first derivative values at the grid points

c..declare the pass
      integer          ngrid
      double precision grid(ngrid),u(ngrid),du(ngrid)


c..local variables
      integer          m,n
      parameter        (m=2, n=3)
      integer          i
      double precision coef(n)


c..first point; one sided
      call fdcoef(m,n,grid(1),grid(1),coef)
      du(1) = coef(1)*u(1) + coef(2)*u(2) + coef(3)*u(3)


c..middle of the grid; central differences
      do i=2,ngrid-1
       call fdcoef(m,n,grid(i),grid(i-1),coef)
       du(i) = coef(1)*u(i-1) + coef(2)*u(i) + coef(3)*u(i+1)
      enddo

c..last point; one sided
      call fdcoef(m,n,grid(ngrid),grid(ngrid-2),coef)
      du(ngrid) = coef(1)*u(ngrid-2)+coef(2)*u(ngrid-1)+coef(3)*u(ngrid)
      return
      end subroutine




      subroutine fss004(grid,ngrid,u,du)
      implicit none
      save

c..this routine computes fourth order accurate second derivatives 
c..on an arbitrary grid

c..input:
c..ngrid       = number of points in the grid, 
c..grid(ngrid) = array of independent values 
c..u(ngrid)    = function values at the grid points

c..output:
c..du(ngrid)   = second derivative values at the grid points

c..declare the pass
      integer          ngrid
      double precision grid(ngrid),u(ngrid),du(ngrid)


c..local variables
      integer          m,n
      parameter        (m=3, n=5)
      integer          i,j
      double precision coef(n)


c..first point; one sided
      call fdcoef(m,n,grid(1),grid(1),coef)
      du(1) = 0.0d0
      do j=1,5
       du(1) = du(1) + coef(j)*u(j)
      enddo


c..second point; one sided
      call fdcoef(m,n,grid(2),grid(1),coef)
      du(2) = 0.0d0
      do j=1,5
       du(2) = du(2) + coef(j)*u(j)
      enddo

c..middle of the grid
      do i=3,ngrid-2
       call fdcoef(m,n,grid(i),grid(i-2),coef)
       du(i) = 0.0d0
       do j=1,5
        du(i) = du(i) + coef(j)*u(i + j - 3)
       enddo
      enddo

c..second to last point; one sided
      call fdcoef(m,n,grid(ngrid-1),grid(ngrid-4),coef)
      du(ngrid-1) = 0.0d0
      do j=1,5
       du(ngrid-1) = du(ngrid-1) + coef(j)*u(ngrid + j - 5)
      enddo

c..last point; one sided
      call fdcoef(m,n,grid(ngrid),grid(ngrid-4),coef)
      du(ngrid) = 0.0d0
      do j=1,5
       du(ngrid) = du(ngrid) + coef(j)*u(ngrid + j - 5)
      enddo

      return
      end subroutine




      subroutine fdcoef(mord,nord,x0,grid,coef)
      implicit none
      save

c..this routine implements simple recursions for calculating the weights
c..of finite difference formulas for any order of derivative and any order
c..of accuracy on one-dimensional grids with arbitrary spacing.

c..from bengt fornberg's article
c..generation of finite difference formulas on arbitrary spaced grids. 
c..math. comp., 51(184):699-706, 1988. 


c..input:
c..mord       = the order of the derivative 
c..nord       = order of accuracy n
c..x0         = point at which to evaluate the coefficients
c..grid(nord) = array containing the grid

c..output:
c..coef(nord) = coefficients of the finite difference formula.


c..declare the pass
      integer          mord,nord
      double precision x0,grid(nord),coef(nord)


c..local variables
      integer          nu,nn,mm,nmmin,mmax,nmax
      parameter        (mmax=8, nmax=10)
      double precision weight(mmax,nmax,nmax),c1,c2,c3,c4,pfac


c..zero the weight array
      do nu=1,nord
       do nn=1,nord
        do mm=1,mord
         weight(mm,nn,nu) = 0.0d0
        enddo
       enddo
      enddo

      weight(1,1,1) = 1.0d0
      c1            = 1.0d0
      nmmin         = min(nord,mord)
      do nn = 2,nord
       c2 = 1.0d0
       do nu=1,nn-1
        c3 = grid(nn) - grid(nu)
        c2 = c2 * c3
        c4 = 1.0d0/c3
        pfac = grid(nn) - x0
        weight(1,nn,nu) = c4 * ( pfac * weight(1,nn-1,nu) )
        do mm=2,nmmin
         weight(mm,nn,nu) = c4 * ( pfac * weight(mm,nn-1,nu)
     1                      - dfloat(mm-1)*weight(mm-1,nn-1,nu) )
        enddo
       enddo
       pfac = (grid(nn-1) - x0)
       weight(1,nn,nn) = c1/c2 * (-pfac*weight(1,nn-1,nn-1))
       c4 = c1/c2
       do mm=2,nmmin
        weight(mm,nn,nn) = c4 * (dfloat(mm-1)*weight(mm-1,nn-1,nn-1) - 
     1                            pfac*weight(mm,nn-1,nn-1))
       enddo
       c1 = c2
      enddo

c..load the coefficients
      do nu = 1,nord
       coef(nu) = weight(mord,nord,nu)
      enddo
      return
      end subroutine


      end module
      
      
