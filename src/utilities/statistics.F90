module statistics
!     This module contains all the procedures to perform statistical analysis on the
!     SDE solutions. This module should not depend on any other module.
      implicit none
      public
      private ttable            
      
contains

      subroutine avevar(dat,n,ave,var)
!     Routine to compute mean and variance of dat. Adopted from "Numerical Recipes
!     in FORTRAN", W. H. Press et al., Cambridge University Press, 1992.
      implicit none
      integer :: n
      real(8) :: ave, var, dat(n)
!     Given array dat(1:n), returns its mean as ave and variance as var.
      integer :: j
      real(8) :: s, ep
      ave = 0.D0
      do j=1,n
         ave = ave + dat(j)
      end do
      ave = ave/DBLE(n)
      var = 0.D0
      ep = 0.D0
      do j=1,n
         s = dat(j) - ave
         ep = ep + s
         var = var + s*s
      end do
!     Corrected two-pass formula to reduce round-off error
      var = (var - ep**2/DBLE(n))/DBLE(n-1)
      return
      end subroutine avevar
      
      
      
      real(8) function cfdlev(a)
!     Returns '(1 - alpha)' according to 'a' from the table in the comments of
!     function TTABLE.
      integer, intent(in) :: a
      select case(a)
      case(1)
         cfdlev = 1.d0 - 0.2d0
      case(2)
         cfdlev = 1.d0 - 0.1d0
      case(3)
         cfdlev = 1.d0 - 0.05d0
      case default
         print *,' wrong argument a given to cfdlev: a =',a
         cfdlev = 0.d0
      end select
      return
      end function cfdlev
      
      
      
      real(8) function ttable(a,n)
!     Student's t-distribution table with significance probability alpha (indicated
!     by a) and degree of freedom n. The table value is returned as ttable. alpha
!     could only have values 0.2, 0.1 and 0.05. The table values are taken from
!     Table 30 in "Manual of Mathematics", High Education Press, Beijing, 1979, and
!     Wikipedia.
!         a| alpha
!       ---+-------
!         1|  0.2
!         2|  0.1
!         3|  0.05
!     table degree of freedom: 1, 2, 3, ..., 29, 30, 40, 50, 60, 80, 100, 120, +inf 
      implicit none
      integer, intent(in) :: a, n
      real(8) :: tt(37,3)
      save tt
      data tt /3.078D0, 1.886D0, 1.638D0, 1.533D0, 1.476D0, 1.440D0, 1.415D0, &
               1.397D0, 1.383D0, 1.372D0, 1.363D0, 1.356D0, 1.350D0, 1.345D0, &
               1.341D0, 1.337D0, 1.333D0, 1.330D0, 1.328D0, 1.325D0, 1.323D0, &
               1.321D0, 1.319D0, 1.318D0, 1.316D0, 1.315D0, 1.314D0, 1.313D0, &
               1.311D0, 1.310D0, 1.303D0, 1.299D0, 1.296D0, 1.292D0, 1.290D0, &
               1.289D0, 1.282D0, &
               6.314D0, 2.920D0, 2.353D0, 2.132D0, 2.015D0, 1.943D0, 1.895D0, &
               1.860D0, 1.833D0, 1.812D0, 1.796D0, 1.782D0, 1.771D0, 1.761D0, &
               1.753D0, 1.746D0, 1.740D0, 1.734D0, 1.729D0, 1.725D0, 1.721D0, &
               1.717D0, 1.714D0, 1.711D0, 1.708D0, 1.706D0, 1.703D0, 1.701D0, &
               1.699D0, 1.697D0, 1.684D0, 1.676D0, 1.671D0, 1.664D0, 1.660D0, &
               1.658D0, 1.645D0, &
              12.706D0, 4.303D0, 3.182D0, 2.776D0, 2.571D0, 2.447D0, 2.365D0, &
               2.306D0, 2.262D0, 2.228D0, 2.201D0, 2.179D0, 2.160D0, 2.145D0, &
               2.131D0, 2.120D0, 2.110D0, 2.101D0, 2.093D0, 2.086D0, 2.080D0, &
               2.074D0, 2.069D0, 2.064D0, 2.060D0, 2.056D0, 2.052D0, 2.048D0, &
               2.045D0, 2.042D0, 2.021D0, 2.009D0, 2.000D0, 1.990D0, 1.984D0, &
               1.980D0, 1.960D0/
            
      if( (a/=1 .and. a/=2 .and. a/=3) .or. n<1 ) then
         print *,' wrong argument given to ttable:'
         print *,' a =',a,' n =',n
         stop ' ttable stopped'
      end if
      
      select case(n)
      case(1:30)
         ttable = tt(n,a)
      case(31:34)
         ttable = tt(30,a)
      case(35:44)
         ttable = tt(31,a)
      case(45:54)
         ttable = tt(32,a)
      case(55:70)
         ttable = tt(33,a)
      case(71:90)
         ttable = tt(34,a)
      case(91:110)
         ttable = tt(35,a)
      case(111:130)
         ttable = tt(36,a)
      case(131:)
         ttable = tt(37,a)
      end select
      
      return
      end function ttable
      
      
      
      subroutine cfdint(dat,n,a,ave,del)
!     Given array dat(1:n), this routine returns its mean and confidence interval
!     with significance probability indicated by a (see function ttabble), so that
!     with confidence (1 - alpha(a)), the expectation of dat lies within
!     [ave - del, ave + del].
      implicit none
      integer, intent(in) :: n, a
      real(8), intent(in) :: dat(n)
      real(8), intent(out) :: ave, del
      real(8) :: var
      call avevar(dat,n,ave,var)
      del = ttable(a,n-1) * sqrt(var/DBLE(n))
      return
      end subroutine cfdint



      function selection(k, n, dat)
!     modified from function SELECT in "Numerical Recipes in FORTRAN", W. H. Press
!     et al., Cambridge University Press, 1992. this function can be used to find
!     median very efficiently.
      implicit none
      integer, intent(in) :: k, n
      real(8), intent(in) :: dat(n)
      real(8) :: selection, arr(n)
!     returns the k-th smallest value in the array arr(1:n). array arr(n) will be
!     rearranged to have this value in location arr(k), with all smaller elements to
!     arr(1:k-1) (in arbitrary order) and all larger elements in arr(k+1:n) (also in
!     arbitrary order).
      integer :: i, ir, j, l, mid
      real(8) :: a, temp

      arr = dat
      l = 1
      ir = n
!     active partition contains 1 or 2 elements      
1     if( (ir-l)<=1 ) then
!        active partition contains 2 elements         
         if( (ir-l)==1 ) then
            if( arr(ir)<arr(l) ) then
               temp = arr(l)
               arr(l) = arr(ir)
               arr(ir) = temp
            end if
         end if
         selection = arr(k)
         return
      else
!        choose median of left, center, and right elements as partition element a.
!        also rearrange so that arr(l+1)<=arr(l), arr(ir)>=arr(l).
         mid = (l + ir)/2
         temp = arr(mid)
         arr(mid) = arr(l+1)
         arr(l+1) = temp
         if( arr(l+1)>arr(ir) ) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
         end if
         if( arr(l)>arr(ir) ) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         end if
         if( arr(l+1)>arr(l) ) then
            temp = arr(l+1)
            arr(l+1) = arr(l)
            arr(l) = temp
         end if
!        initialize pointers for partitioning
         i = l + 1
         j = ir
!        partitioning element
         a = arr(l)
!        beginning of inner loop
3        continue
!        scan up to find element > a
         i = i + 1
         if( arr(i)<a ) goto 3
4        continue
!        scan down to find element < a
         j = j - 1
         if( arr(j)>a ) goto 4
!        pointer crossed. exit with partitioning complete.
         if( j<i ) goto 5
!        exchange elements
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
!        end of inner loop
         goto 3
!        insert partitioning element
5        arr(l) = arr(j)
         arr(j) = a
!        keep active the partition that contains the k-th element
         if( j>=k ) ir = j - 1
         if( j<=k ) l = i
      end if
      goto 1
      end function selection
      
end module statistics



!#define _TEST2_
#ifdef _TEST_
!     This is the example given in "Textbook for Advanced Mathematics", Vol. III,
!     Zhong Li, Jianying Zhou, Peking University Press, 1999, pp.272-273.
      program test
      use statistics
      implicit none
      integer, parameter :: N = 10
      real(8) :: dat(N)
      real(8) :: ave, del, var
      data dat /3350.d0, 3348.d0, 3351.d0, 3349.d0, 3351.d0, 3352.d0, 3350.d0, &
                3349.d0, 3347.d0, 3352.d0/
!     calculate the 95% confidence interval
      call cfdint(dat,N,3,ave,del)
      print *,' Elevation of the moutain is between (95% confidence):'
      print *,' [', ave - del,', ', ave + del, '] m'
      stop
      end program test
#endif
#ifdef _TEST2_
      program test2
      use statistics
      implicit none
      integer, parameter :: N = 9
      real(8) :: dat(N)
      real(8) :: median
      data dat /-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0/
      median = selection((N+1)/2, N, dat)
      print *,' median =',median
      end program test2
#endif






      
