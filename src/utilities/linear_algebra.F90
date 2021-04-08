module linear_algebra
!     This module contains all subroutines to perform linear algebraic calculations
!     for the REM. This module uses the LAPACK library (version 3.5.0), obtained from
!     "http://www.netlib.org/lapack/". This module shall not depend on any other module.

      implicit none
      public
      private :: minor, assembly, lower_triangle
 
contains
      
      subroutine print_matrix(lda,A,n)
!     prints the LDA x N matrix A on monitor
      implicit none
!     n has to be 1 ~ 9
      integer, intent(in) :: lda, n
      real(8), dimension(lda,n), intent(in) :: A
      integer :: i
      character(18) :: fmtstr
      
      write(unit=fmtstr,fmt='(A3,I1,A14)') '(A,',n,'(1X,ES10.3),A)'      
      do i=1,lda
         write(6,fmt=fmtstr) ' |',A(i,:),' |'
      end do
      print *,''
      
      return
      end subroutine print_matrix
      
      
      
      subroutine symmetric_matrix(na,V,A)
!     this function returns the na x na symmetric matrix A composed of the entires as
!     listed in V. Length of V should be na*(na+1)/2. The entries of V are A(1,1),
!     A(1,2), A(2,2), A(1,3), A(2,3), ...
      implicit none
      integer, intent(in) :: na
      real(8), intent(in) :: V(na*(na+1)/2)
      real(8), intent(out) :: A(na,na)
      integer :: i, j, p
      
      p = 0
      do j=1,na
         do i=1,j
            p = p + 1
            A(i,j) = V(p)
            if( i/=j ) A(j,i) = V(p)
         end do
      end do
      
      return
      end subroutine symmetric_matrix


 
      real(8) function trace(na,A)
!     this routine returns trace of an na x na matrix A.
      implicit none
      integer, intent(in) :: na
      real(8), intent(in) :: A(na,na)
      integer :: i
      trace = 0.d0
      do i=1,na
         trace = trace + A(i,i)
      end do
      return
      end function trace
      
      
      
      real(8) function det(na,mat)
!     this routine calculates determinant of an na x na matrix mat.
      implicit none
      integer, intent(in) :: na
      real(8), intent(in) :: mat(na,na)
      real(8) :: a(na,na)
      real(8) :: m, temp
      integer :: i, j, k, l
      logical :: detexists = .true.

!     make a copy of mat so that mat is not altered     
      a = mat

      l = 1
!     convert a to upper triangular form (Gaussian elimination)
      do k=1,na-1
         if( a(k,k)==0.d0 ) then
            detexists = .false.
            do i=k+1,na
               if( a(i,k)/=0.d0 ) then
                  do j=1,na
                     temp = a(i,j)
                     a(i,j) = a(k,j)
                     a(k,j) = temp
                  end do
                  detexists = .true.
                  l = -l
                  exit
               end if
            end do
            if( .not.detexists ) then
               det = 0.d0
               return
            end if
         end if
         do j=k+1,na
            m = a(j,k)/a(k,k)
            do i=k+1,na
               a(j,i) = a(j,i) - m*a(k,i)
            end do
         end do
      end do
      
!     calculate determinant by finding product of diagonal elements
      det = dble(l)
      do i=1,na
         det = det*a(i,i)
      end do
      
      return
      end function det
      
      
      
      subroutine minor(na,A,np,P,D)
!     Given a matrix A(1:na,1:na), this routine returns its principal minor P(1:np,1:np)
!     obtained from deleting the i-th row and i-th column of A where A(i,i) <= 0, and 
!     an integer vector D that shows the positions of A's non-positive diagonals. After
!     one call of this routine, matrix P must be explicitly deallocated to allow the
!     next call.
      implicit none
      integer, intent(in) :: na
      integer, intent(out) :: np
      real(8), intent(in) :: A(na,na)
!     D(i) = 0 where A(i,i) <= 0; D(i) = 1 where A(i,i) > 0.
      integer, intent(out) :: D(na)
      real(8), allocatable :: P(:,:), T(:)
      integer :: i
      logical :: mask(na,na)
      
      D = 1
      np = na
      do i=1,na
         if( A(i,i)<=0.d0 ) goto 100
      end do
!     if no diagonal entry of A is non-positive, simply copy A to P
      allocate(P(np,np))
      P = A
      return
      
100   continue
      mask = .true.
      do i=1,na
         if( A(i,i)<=0.d0 ) then
            mask(i,:) = .false.
            mask(:,i) = .false.
            D(i) = 0
            np = np - 1
         end if
      end do

!     PACK and RESHAPE can take zero-sized arrays      
      allocate(T(np*np),P(np,np))
      T = PACK(A,mask)
      P = RESHAPE(T,(/np,np/))
      deallocate(T)
      
      return
      end subroutine minor
      
      
      
      subroutine assembly(na,A,np,P,D)
!     This is the partial inverse routine of MINOR. It assembles the na x na matrix
!     A from the np x np (np <= na) principal minor P, by filling the rows of P into A
!     at which D(i) = 1, and rows of zeros into A at which D(i) = 0. When filling,
!     the rows of P are extended to length na by padding (na - np) zeros. See the
!     document "Factorization of Real Symmetric Matrices in the REM" for more detailed
!     explanations.
      implicit none
      integer, intent(in) :: na, np
      integer, intent(in) :: D(na)
      real(8), intent(in) :: P(np,np)
      real(8), intent(out) :: A(na,na)
      integer :: ia, ip
      
      if( SUM(D)/=np ) then
         print *,'Error: inconsistent P and D provided to subroutine assembly'
         print *,'np =',np
         print *,'P ='
         call print_matrix(np,P,np)
         print *,'D =',D
         stop ' subroutine assembly stopped'
      end if
      
      if( np==na ) then
         A = P
      elseif( np<na ) then
         A = 0.d0
         ip = 0
         do ia=1,na
            if( D(ia)==1 ) then
               ip = ip + 1
               A(ia,1:np) = P(ip,1:np)
            end if
         end do
      else
         print *,'Error: np > na in subroutine assembly'
         print *,'np =',np,', na =',na
         stop ' subroutine assembly stopped'
      end if
      
      return
      end subroutine assembly
      
      
      
      subroutine lower_triangle(na,A,L)
!     subroutine to extract the lower triangle L from na x na matrix A
      implicit none
      integer, intent(in) :: na
      real(8), dimension(na,na), intent(in) :: A
      real(8), dimension(na,na), intent(out) :: L
      integer :: i
      L = 0.d0
      do i=1,na
         L(i,1:i) = A(i,1:i)
      end do
      return
      end subroutine lower_triangle
      
      
      
      subroutine factorization(na,A,L,err)
!     This routine factorizes the na x na real symmetric matrix A as
!        A = L * L**T.
!     err is an integer flag that
!        err = 0 if A is positive semi-definite or slightly indefinite;
!        err = 1 if A is significantly indefinite.
!     This condition is determined by comparing the magnitude of A's smallest
!     negative eigenvalue (if it ever has) with its biggest eigenvalue.
!
!     The factorization is done in three steps:
!     1) Extract the principal minor P of A by deleting the rows and columns of A's
!        non-positive diagonals. P = A if all A's diagonals are positive.
!     2) Try Cholesky factorization on P. If successful, construct the matrix L from
!        the result of Cholesky factorization of P; if failed, proceed to step 3.
!     3) Spectral-factorize P since it is not positive definite as indicated by 
!        failure of step 2. Then construct the matrix L from the result of spectral
!        factorization.
!
!     For detailed explanations on how this routine works, see the document
!     "Factorization of Real Symmetric Matrices in the REM".
      implicit none
      integer, intent(in) :: na
      integer, intent(out) :: err
      real(8), intent(in) :: A(na,na)
      real(8), intent(out) :: L(na,na)
      character(1), parameter :: UPLO = 'L'
!     the optimal length of WORK for 3x3 matrices as returned from a (separate) test
!     program
      integer, parameter :: LWORK = 102
      integer :: INFO, i, np
      integer :: D(na)
      real(8) :: WORK(LWORK)
!     T is a temparory copy of P. S is the factorization result of T. W is the list
!     of eigenvalues of P, sorted from the smallest to the greatest. V is the matrix
!     square root of the matrix diag(W).
      real(8), allocatable :: W(:), P(:,:), T(:,:), S(:,:), V(:,:)

      L = 0.d0
      err = 0
!     extract the principal minor P
      call minor(na,A,np,P,D)      
      if( np==0 ) then
         deallocate(P)
         return
      end if

!     copy P to T since the LAPACK factorization routines destroy the input array
      allocate( T(np,np), S(np,np) )
      T = P
      
!     try Cholesky factorization on T
      call DPOTRF(UPLO, np, T, np, INFO)
      
      if( INFO==0 ) then
      
!        extract the Cholesky factorization result S from T
         call lower_triangle(np,T,S)
         
      elseif( INFO>0 ) then
!        Cholesky factorization failed, P is not positive definite

!        copy P to T again because T has been destroyed by DPOTRF
         T = P
!        spectral-factorize T
         allocate( W(np), V(np,np) )
         call DSYEV('V', UPLO, np, T, np, W, WORK, LWORK, INFO)
         
         if( INFO==0 ) then
!           check whether the matrix is too far away from positive semi-definite, give
!           warning if it is the case
            if( W(1)<0.d0 ) then
               if( abs(W(1))>=5.d-2*W(np) ) then
                  print *,'Warning: abs(W(1)) >= 0.05 * W(np)'
                  print *,'matrix is significantly indefinite'
                  print *,'W =',W
                  print *,'P ='
                  call print_matrix(np,P,np)
                  print *,'A ='
                  call print_matrix(na,A,na)
                  err = 1
               end if
            end if
!           calculate square root of the eigenvalue diagonal matrix
            V = 0.d0
            do i=1,np
               if( W(i)>0.d0 ) V(i,i) = sqrt(W(i))
            end do
!           this is the factorization result of T
            S = MATMUL(T,V)
         else
            print *,'Error: spectral factorization of T failed'
            print *,'INFO =',INFO
            print *,''
            print *,'P ='
            call print_matrix(np,P,np)
            print *,'A ='
            call print_matrix(na,A,na)
            stop ' subroutine factorization stopped'
         end if
         
         deallocate(W,V)
         
      else
         
         print *,'Error: wrong input argument list to DPOTRF'
         print *,'INFO =',INFO
         stop ' subroutine factorization stopped'
         
      end if

!     assemble L from S
      call assembly(na,L,np,S,D)
      deallocate(T,S,P)
      
      return
      end subroutine factorization
      
      
      
      subroutine cnan(na,A,pass)
!     This routine checks for NaN's in the na x na matrix A. pass = 1 is returned if
!     all entries of A are not NaN; pass = 0 otherwise.
      implicit none
      integer, intent(in) :: na
      real(8), intent(in) :: A(na,na)
      integer, intent(out) :: pass
      integer :: i, j
      
      do i=1,na
         do j=1,na
            if( ISNAN(A(i,j)) ) then
               pass = 0
               return
            end if
         end do
      end do
      
      pass = 1
      return
      end subroutine cnan
      
      
      
      subroutine cspd(na,A,pass)
!     This routine checks the positive semi-definiteness of the na x na real symmetric
!     matrix A. pass = 1 is returned if A is positive semi-definite; pass = 0
!     otherwise.
      implicit none
      integer, intent(in) :: na
      real(8), intent(in) :: A(na,na)
      integer, intent(out) :: pass
      integer :: i, np, INFO
      integer :: D(na)
      real(8), allocatable :: P(:,:)

!     diagonal entries can not be negative
      do i=1,na
         if( A(i,i)<0.d0 ) then
            pass = 0
            return
         end if
      end do

!     extract the principal minor P
      call minor(na,A,np,P,D)      
      if( np==0 ) then
         pass = 1
         deallocate(P)
         return
      end if

!     try Cholesky decomposition on P
      call DPOTRF('L', np, P, np, INFO)
      if( INFO==0 ) then
         pass = 1
      elseif( INFO>0 ) then
         if( abs(det(np,P))<=EPSILON(1.d0)*trace(np,P) ) then
            pass = 1
         else
            pass = 0
         end if
      else
         print *,'Error: wrong input argument list to DPOTRF'
         print *,'INFO =',INFO
         pass = 0
      end if
      
      deallocate(P)
      return
      end subroutine cspd
            
end module linear_algebra
      
      
      
