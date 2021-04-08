module mod_equation_typedef
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 01/19/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      implicit none
      public

!     Consider a Boltzmann equation up to three dimensions in the form (e.g.,
!     Eq. (180) in Schulz [1991]):
!
!        dt(u) + (1/G)*di(G*hi*u) = (1/G)*di[G*Dij*dj(u)] + S*u + v,             (1)
!
!     where dt = d/dt, di = d/dxi, and summation convention is assumed over repeated
!     indices. Eq. (1) is transformed into the Kolmogorov backward equation form as
!
!        dt(u) = Dij*dij(u) + [dj(Dij) + Dij*dj(lnG) - hi]*di(u)
!                + [S - di(hi) - hi*di(lnG)]*u + v,                              (2)
!
!     where dij = d2/dxidxj. The coefficients used to define Eq. (2) are expressed
!     by the following type. The term dj(lnG) related to the Jacobian determinant G
!     will be expressed as a user-provided function in the user-editable source file
!     user_input.F90.
      type coefficients
!        Dij = [D11, D12, D22, D13, D23, D33]
         real(8) :: Dij(6)
!        djDij = [dj D1j, dj D2j, dj D3j]
         real(8) :: djDij(3)
!        hi = [h1, h2, h3]
         real(8) :: hi(3)
         real(8) :: dihi
         real(8) :: S
         real(8) :: v
      end type coefficients

!     Comparing Eq. (2) with the standard form
!
!        dt(u) = (1/2)*Aij*dij(u) + Bi*di(u) + c*u + v,                          (3)
!
!     we have got the coefficients of Eq. (3) defined by the following type, which
!     are useful components in constructing the corresponding Ito stochastic
!     processes and the functional F in the Feynman-Kac formula u = E[F].
      type components
!        symmetric positive semi-definite matrix A = 2*D
         real(8) :: A(3,3)
!        B = div.D - h
         real(8) :: B(3)
!        h is used in expressing lambda in the third type boundary condition
         real(8) :: h(3)
!        c = S - div.h
         real(8) :: c
         real(8) :: v
      end type components

end module mod_equation_typedef




