module mod_solver_ito_process
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
      use mod_typedef_params, dummy_ParamFname => ParamFname, &
          & dummy_initialize_params => initialize_params
      use mod_domain, dummy_domain_initial_condition => domain_initial_condition, &
          & dummy_initialize_domain => initialize_domain, &
          & dummy_finalize_domain => finalize_domain
      use mod_equation, dummy_initialize_equation => initialize_equation, &
          & dummy_finalize_equation => finalize_equation
      use mod_solver_random_number, only: gausrand, exprand
      use mod_solver_solution_grid, only: &
          & initial_condition => solver_initial_condition
      use interpolation
      use linear_algebra
      implicit none
      private
      public ito_proc_params, ito_process
! Integrate an Ito process using Euler-Maruyama method, with path integrals
! involving source and loss terms and reflection term.
! References:
!    [1] L. Zheng (2015), Development and application of stochastic methods for
!        radiation belt simulations, Ph.D. thesis, Rice University
!    [2] B. Oksendal (2000), Stochastic differential equations: An introduction with
!        applications, fifth edition, Springer
!    [3] C. Costantini, B. Pacchiarotti and F. Sartoretto (1998), Numerical
!        approximation for functionals of reflecting diffusion processes, SIAM J.
!        Appl. Math, vol. 58, No. 1, pp. 73-102
!    [4] E. Gobet (2001), Euler schemes and half-space approximation for the
!        simulation of diffusion in a domain, ESAIM: Probability and Statistics,
!        vol. 5, pp. 261-297
!    [5] D. Lepingle (1995), Euler scheme for reflected stochastic differential
!        equations, Mathematics and computers in simulation, vol. 38, pp. 119-126
!    [6] W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery (1992),
!        Numerical recipes in FORTRAN, second edition, Cambridge University Press

      type ito_proc_params
!        critical F value used to control fission
         real(8) :: criticalF
!        total time length of solution
         real(8) :: t_total
!        mean path integrals of from all previous simulations        
         real(8) :: medianU
         real(8) :: medianV
!        see the nomenclatures below 
         real(8) :: pastU
         real(8) :: pastV
      end type ito_proc_params

      type ito_process
!        starting spacetime
         type(spacetime) :: xt_start
!        ending spacetime
         type(spacetime) :: xt_end
!        path integrals, see the nomenclatures below
         real(8) :: Ut, Vt
!        boundary id if the ito process stops on a Dirichlet type boundary,
!        otherwise bid = 0
         integer :: bid
!        stopping time
         real(8) :: tau
!        local time
         real(8) :: kt
!        total number of steps
         integer :: nsteps
!        simulation path length. mean stepsize = pathlen/nsteps
         real(8) :: pathlen
!        maximum stepsize
         real(8) :: maxstep
      contains
         procedure :: walk
      end type ito_process

      interface ito_process
         module procedure start
      end interface ito_process

contains

      function start()
      type(ito_process) :: start

      start%xt_start = spacetime(0.d0, 0.d0, 0.d0, 0.d0)
      start%xt_end = start%xt_start
      start%Ut = 1.d0
      start%Vt = 0.d0
      start%bid = 0
      start%tau = 0.d0
      start%kt = 0.d0
      start%nsteps = 0
      start%pathlen = 0.d0
      start%maxstep = 0.d0

      end function start



      subroutine walk(this, xt_start, t_stop, fissile, params)
!     implement one time-backward ito process segment from xt_start till any of the
!     stopping criteria being met:
!        1) end of time (t = t_stop >= 0)
!        2) hitting a Dirichlet type boundary
!        3) if fissile, the projected functional value of the entire ito process is
!           greater than or equal to a prescribed critical value
!
!     Nomenclature of the path integrals in this routine (T = t_stop + t_total):
!
!        U = exp(Y - Z) = exp( int{c(T-s,Xs)*ds} - int{lambda(T-s,Xs)*dks} ),    (1)
!        V = int{ v(T-s,Xs)*Us*ds }.                                             (2)
!
!     According to the limits of the integral parameter s, U and V are given
!     different names:
!
!       t=T  xt_start    t    xt_end   t=t_stop  t=0
!        +-------+-------+------+--------+--------+---->  time axis
!       s=0     s=s1    s=s'   s=s2    s=t_total
!        | pastU |      Ut      |  futrU |
!        | pastV |      Vt      |  futrV |
!
!     Ut and Vt are integrated in this routine. pastU and pastV are inputs to this
!     routine. futrU and futrV are estimated in this routine using the input medianU
!     and medianV from all previously simulated paths.
      class(ito_process) :: this
      type(spacetime), intent(in) :: xt_start
      real(8), intent(in) :: t_stop
      logical, intent(in) :: fissile
      type(ito_proc_params), intent(in) :: params
!     check the fission condition every Fnsteps steps ("those who run 50 steps laugh
!     at those who run 100 steps.")
      integer, parameter :: Fnsteps = 50
      type(spacetime) :: xt_curr, xt_next
      type(spacetime) :: xt_proj, xt_next_proj, xt_int, xt_map
      type(components) :: ito_cmpnt, ito_cmpnt_proj, ito_cmpnt_next
      type(boundary_condition) :: bcondition
      type(boundary), pointer :: bptr => null()
      integer :: nsteps
      integer :: err
      integer :: bid_curr, bid_next, bid_map
      integer :: marker, marker_map
      real(8) :: rms_ds, step, tp
      real(8) :: dt, dkt, sqrtdt, delt
      real(8) :: W(3), Sigma(3,3)
      real(8) :: vgamma(3), vgamma_dkt(3), dummy_vgamma(3), Fgamma(3)
      real(8) :: lambda, dummy_lambda
      real(8) :: Ut, Vt, Yt, Zt
!     projected functional value of the entire ito process and future path integrals
      real(8) :: projF, futrU, futrV
      logical :: success, end_of_path
      !$OMP threadprivate(bptr)

      if( t_stop<0.d0 ) then
         print *,' mod_solver_ito_process procedure ito_process%walk:'
         print '(A,F10.4,A)','  t_stop =',t_stop,' < 0'
         stop
      end if

      this%xt_start = xt_start
      
!     no need to walk the ito process if t already hits t_stop
      if( (xt_start%t - t_stop)<dtMin ) then
         this%xt_end = xt_start
         return
      end if

!     the ito process cannot start off the domain
      call domain_query(xt_start, marker, bid_curr)
      if( marker==-1 ) then
         print *,' mod_solver_ito_process procedure ito_process%walk:'
         print '(A,I3)','  xt_start is out of domain, bid =',bid_curr
         print '(A,3(F10.4,A))','  xt_start = (',xt_start%x1,',',xt_start%x2,',' &
                                & ,xt_start%x3,')'
         stop
      end if

      Yt = 0.d0
      Zt = 0.d0
      Ut = 1.d0
      Vt = 0.d0
      end_of_path = .false.
      nsteps = 0
      xt_curr = xt_start

!     get the ito components at xt_curr
      ito_cmpnt = equation_cmpnt(xt_curr)

      do while( .not.end_of_path )

!        [0] boundary considerations: root-mean-square spatial stepsize is reduced
!        by half in boundary neighborhood.
         if( marker==0 ) then
            bptr => boundaries(bid_curr)
            rms_ds = 0.5d0*dS
         else
            bptr => null()
            rms_ds = dS
         end if

!        [1] estimate the temporal stepsize dt corresponding to rms_ds
         dt = atss(rms_ds, 3, ito_cmpnt%A, ito_cmpnt%B)
!        temporal stepsize of the last step
         if( (xt_curr%t - dt)<=t_stop ) then
            dt = max(xt_curr%t - t_stop, 0.d0)
            end_of_path = .true.
         end if
         sqrtdt = sqrt(dt)

!        [2] factorize A = Sigma*Sigma^T to obtain the SDE coefficients. N.B., the
!        matrix Sigma is also used in the private function half_space_local_time in
!        step [4].
         call factorization(3, ito_cmpnt%A, Sigma, err)
         if( err==1 ) then
            print *,''
            print *,' mod_solver_ito_process procedure ito_process%walk:'
            print *,'  Indefinite A ='
            call print_matrix(3, ito_cmpnt%A, 3)
            print '(A,4(F10.4,A))','  at xt = (',xt_curr%x1,',',xt_curr%x2,',' &
                               & ,xt_curr%x3,',',xt_curr%t,')'
            print *,''
         end if

200      continue
!        [3] draw three independent standard Gaussian random numbers. N.B., W is
!        also used in the private function half_space_local_time in step [4].
         W(1) = gausrand()
         W(2) = gausrand()
         W(3) = gausrand()

!        [4] if xt_curr is in neighborhood of a Neumann type boundary, find the
!        reflection terms vgamma, dkt, and evaluate lambda on the boundary
         dkt = 0.d0
         vgamma_dkt = 0.d0
         lambda = 0.d0
         if( marker==0 ) then
!           bptr has been associated in step [0]
            if( bptr%b_type==2 ) then
!              get the Neumann type boundary condition at xt_curr
               bcondition = bptr%condition(xt_curr)

!              obtain vgamma at xt_curr
               call third_type_bc(ito_cmpnt, bcondition, vgamma, dummy_lambda)
!              try to find the projection point xt_proj on this boundary
               call projection(xt_curr, bid_curr, vgamma, xt_proj, tp, success)

               if( success ) then
                  dkt = half_space_local_time()
!@                  print '(A,4(F10.4,A))','  xt_proj = (',xt_proj%x1,',' &
!@                            & ,xt_proj%x2,',',xt_proj%x3,',',xt_proj%t,')'
!@                  print *,' dkt =',dkt
!                 re-draw ito component and bcondition at xt_proj
                  ito_cmpnt_proj = equation_cmpnt(xt_proj)
                  bcondition = bptr%condition(xt_proj)
!                 evaluate lambda at xt_proj to be used in step [7]
                  call third_type_bc(ito_cmpnt_proj, bcondition, dummy_vgamma, &
                                   & lambda)
               end if

!              if a projection point was not found, though very unlikely but
!              which could happen if vgamma is almost tangent to the boundary,
!              so that the projection point is very far away from xt_curr, dkt = 0
!              and no reflection occurs. note that, vgamma here is still the gamma
!              vector at xt_curr rather than xt_proj.
               vgamma_dkt = vgamma*dkt
            end if
         end if

!        [5] advance one step using Euler-Maruyama method:
!
!           X(i+1) = X(i) + B(t(i),X(i))*dt(i) + Sigma(t(i),X(i)).dW(i) +
!                    vgamma(X(i))*dkt(i)
!
!        update spatial position
         xt_next%x1 = xt_curr%x1 + ito_cmpnt%B(1) * dt + &
                    & dot_product(Sigma(1,:), W) * sqrtdt + vgamma_dkt(1)
         xt_next%x2 = xt_curr%x2 + ito_cmpnt%B(2) * dt + &
                    & dot_product(Sigma(2,:), W) * sqrtdt + vgamma_dkt(2)
         xt_next%x3 = xt_curr%x3 + ito_cmpnt%B(3) * dt + &
                    & dot_product(Sigma(3,:), W) * sqrtdt + vgamma_dkt(3)
!        without updating time
         xt_next%t = xt_curr%t

100      continue
!        [6] check if xt_next has crossed a boundary
         call domain_query(xt_next, marker, bid_next)
         if( marker==-1 ) then
!           xt_curr must have been in the neighborhood of the same boundary when
!           xt_next crosses the boundary. in other words, a big step larger than
!           the neighborhood width at crossing is not allowed. for a Dirichlet type
!           boundary, this restriction ensures accuracy of the stopping position;
!           for a Neumann type boundary, knowledge of vgamma at xt_curr is required
!           to calculate the projection at xt_next.
            if( bid_next==bid_curr ) then
!              bptr has been pointing to the right boundary
               if( bptr%b_type==1 ) then
!                 find the intersection point if xt_next crosses a Dirichlet type
!                 boundary and assign the boundary id
                  call intersection(xt_curr, xt_next, bid_next, xt_int, tp)
                  xt_next = xt_int
                  dt = tp*dt
                  end_of_path = .true.
                  this%bid = bid_next
               elseif( bptr%b_type==2 ) then
!                 find the projection point if xt_next crosses a Neumann type
!                 boundary. note that vgamma is the gamma vector at xt_curr.
                  call projection(xt_next, bid_next, vgamma, xt_next_proj, tp, &
                                & success)
                  if( success ) then
                     if( tp>=0.d0 ) then
!                       this is the displacement vector from xt_next to xt_next_proj
!                       its norm also contributes to the local time increment (c.f.
!                       Ref. [4])
                        Fgamma = vgamma*tp
                        dkt = dkt + sqrt(dot_product(Fgamma, Fgamma))
                        xt_next = xt_next_proj
!@                        print '(A,4(F10.4,A))','  xt_next_proj = (',xt_next_proj%x1,&
!@                        & ',',xt_next_proj%x2,',',xt_next_proj%x3, &
!@                        & ',',xt_next_proj%t,')'
!@                        print *,' dkt =',dkt
!                       do the boundary check again, because the projected point may
!                       have been out of a Dirichlet type boundary. the case that
!                       the projected points is out of another Neumann type boundary
!                       with a different normal vector will not happen, because the
!                       SDE theory requires vgamma to be C4.
                        goto 100
                     else
!                       tp may not be < 0 because vgamma points inward and xt_next
!                       is out of the boundary.
                        print *,' mod_solver_ito_process procedure ito_process%walk:'
                        print *,'  tp =',tp,' < 0 occurred'
                        print *,'  bid =',bid_next
                        print '(A,3(F10.4,A))','  vgamma = (',vgamma(1),',', &
                              & vgamma(2),',',vgamma(3),')'
                        print '(A,4(F10.4,A))','  xt_next = (',xt_next%x1,',' &
                              & ,xt_next%x2,',',xt_next%x3,',',xt_next%t,')'
                        print '(A,4(F10.4,A))','  xt_next_proj = (', &
                              & xt_next_proj%x1,',',xt_next_proj%x2,',', &
                              & xt_next_proj%x3,',',xt_next_proj%t,')'
                        stop
                     end if
                  else
!                    failed to find a projection point. this means the current
!                    step is too large. reject and re-do the Euler step.
                     goto 200
                  end if
               else
                  print *,' mod_solver_ito_process procedure ito_process%walk:'
                  print *,'  unrecognized b_type =',bptr%b_type
                  stop
               end if
            else
!              a step larger than the neighborhood width occurred when the ito
!              process crosses the boundary. reject and re-do this Euler step.
               goto 200
            end if
         end if

!        update time. note that time runs backwards.
         xt_next%t = xt_curr%t - dt

!        [7] evaluate the path integrals and other accumulants
         ito_cmpnt_next = equation_cmpnt(xt_next)
         Vt = Vt + ito_cmpnt_next%v * Ut * dt
         Yt = Yt + 0.5d0*(ito_cmpnt%c + ito_cmpnt_next%c)*dt
         Zt = Zt + lambda*dkt
         Ut = exp(Yt - Zt)
         
!        accumulate stopping time, local time
         this%tau = this%tau + dt
         this%kt = this%kt + dkt
!        update simulation statistics
         step = dist(xt_curr, xt_next)
         this%pathlen = this%pathlen + step
         this%maxstep = max(this%maxstep, step)

!        [8] end this ito process segment for the upcoming fission (in SOLVER_BATCH)
!        if the projected functional value for the entire path is greater than or
!        equal to the critical value. the projected functional value is estimated by
!
!           projF = pastU*{Ut*[g(xt_map)*futrU + futrV] + Vt} + pastV,           (3)
!
!        where xt_map is the map of xt_next onto t = t_stop plane along the tangent
!        given by the vector ito_cmpnt_next%B, and g is the initial condition.
!        strictly speaking, the mapping should have been done by integrating along
!        the curve whose tangent direction at every point is given by the function
!        B(xt); however, this is too expensive to do. in anticipation that t_total
!        would not be too large, and hence B(xt) does not vary drastically within,
!        B(xt_next) could be a good enough approximation. if xt_map is out of the
!        domain, no fission occurs. also note that, the ito process is allowed to
!        fission right after the first step.
         associate( criticalF=>params%criticalF, pastU=>params%pastU, &
                  & pastV=>params%pastV )
            if( .not.end_of_path .and. fissile .and. criticalF>0.d0 .and. &
              & mod(nsteps, Fnsteps)==0 ) then

               delt = xt_next%t - t_stop
               xt_map%x1 = xt_next%x1 + ito_cmpnt_next%B(1)*delt
               xt_map%x2 = xt_next%x2 + ito_cmpnt_next%B(2)*delt
               xt_map%x3 = xt_next%x3 + ito_cmpnt_next%B(3)*delt
               xt_map%t = t_stop
               call domain_query(xt_map, marker_map, bid_map)

               if( marker_map>=0 ) then
                  call estimate_future_integrals(xt_next%t, futrU, futrV)
                  projF = initial_condition(xt_map)*pastU*Ut*futrU + &
                        & pastU*Ut*futrV + pastU*Vt + pastV
                  if( projF>=criticalF ) end_of_path = .true.
               end if
            end if
         end associate

!        [9] update the iteration items
         xt_curr = xt_next
         bid_curr = bid_next
         ito_cmpnt = ito_cmpnt_next
         nsteps = nsteps + 1

      end do

!     check if the calculation created NaN
      if( isnan(xt_curr%x1) .or. isnan(xt_curr%x2) .or. isnan(xt_curr%x3) .or. &
        & isnan(xt_curr%t) ) then
         print *,' mod_solver_ito_process procedure ito_process%walk:'
         print *,'  NaN occurred:'
         print '(A,4(F10.4,A))','  xt_curr = (',xt_curr%x1,',',xt_curr%x2,',' &
                            & ,xt_curr%x3,',',xt_curr%t,')'
         stop
      end if

      this%xt_end = xt_curr
      this%Ut = Ut
      this%Vt = Vt
      this%nsteps = nsteps

      nullify(bptr)
      contains

         function half_space_local_time() result(dkt)
!        evaluate the half space local time dkt from the information of xt_curr,
!        xt_proj, bcondition%n, W, dt, vgamma, Sigma, and ito_cmpnt%B in the
!        parent routine. the algorithm is adopted from references [4] and [5].
         real(8) :: dkt
         real(8) :: n(3)
         real(8) :: n_dot_gamma
!        displacement from xt_curr to xt_proj
         real(8) :: disp(3)
!        coefficients a and c as appeared on p.281 of [4]
         real(8) :: a(3), a2, c, d
         real(8) :: a_dot_U, V

         disp(1) = xt_proj%x1 - xt_curr%x1
         disp(2) = xt_proj%x2 - xt_curr%x2
         disp(3) = xt_proj%x3 - xt_curr%x3

         n = bcondition%n
         n_dot_gamma = dot_product(n, vgamma)

!        note, coefficient a is a left-product between a vector n and an asymmetric
!        matrix Sigma, therefore their order in matmul matters.
         a = - matmul(n, Sigma)
         c = - dot_product(n, ito_cmpnt%B)
!        note that d < 0 since disp points toward the boundary and n is an inward
!        normal vector
         d = dot_product(n, disp)

!        draw an independent exponential random number V with parameter
!        p = 1/(2*dt), i.e., a random variable V with probability distribution
!        f(V) = p*exp(-p*V), which has <V> = 1/p and <V**2> = 2/p**2. therefore the
!        dimension of V with such a parameter is T.
         V = exprand()
         V = V*2.d0*dt

         a2 = dot_product(a, a)
         a_dot_U = dot_product(a, W)*sqrtdt
         dkt = 0.5d0*(a_dot_U + c*dt + sqrt(a2*V + (a_dot_U + c*dt)**2)) + d

!        theoretically, n_dot_gamma = n.A.n/|n.A| would be positive because A is
!        positive definite. but numerically anything happens.
         if( n_dot_gamma>0.d0 ) then
            dkt = dkt/n_dot_gamma
            dkt = max(dkt, 0.d0)
         else
            dkt = 0.d0
         end if

         end function half_space_local_time


         subroutine estimate_future_integrals(t, futrU, futrV)
!        this routine estimates the future path integrals given a time t between
!        [t_stop, t_stop+t_total], which corresponds to the integral parameter
!        s = s' = t_stop + t_total - t. the estimation is done by assuming that
!
!           <U> = exp(<c>*int{ds}-<lambda>*int{dks}) = exp[(<c>-<lambda*k>)*T'], (4)
!           <V> = <v>int{ exp(<c>*int{dr} - <lambda>*int{dkr})*ds },             (5)
!
!        (T' = t_total), which yield, after some calculations, that if <U> = 1,
!
!           futrU = 1,                                                           (6)
!           futrV = <V>*(1 - s'/T'),                                             (7)
!
!        and otherwise,
!
!           futrU = <U>/exp[(s'/T')*ln(<U>)],                                    (8)
!           futrV = [<V>/(<U> - 1)]*{<U> - exp[(s'/T')*ln(<U>)]}.                (9)
!
!        these estimated future integrals are used to project the funtional value
!        of the entire stochastic path from s=0 till s=T'.
         real(8), intent(in) :: t
         real(8), intent(out) :: futrU, futrV
         real(8), parameter :: Small = 1.d-3
         real(8) :: s_prime, exp_func
         associate( t_total=>params%t_total, medianU=>params%medianU, &
                  & medianV=>params%medianV )

            if( t>=t_stop .and. t<=(t_total + t_stop) ) then
               s_prime = t_total + t_stop - t
            else
               print *,&
               & ' mod_solver_ito_process procedure walk:estimate_future_integrals:'
               print *,'  error: t (= ',t,') < t_stop (=',t_stop,') or'
               print *,'       > t_total + t_stop (=',t_total+t_stop,')'
               stop
            end if
         
            if( abs(medianU - 1.d0)<Small ) then
               futrU = 1.d0
               futrV = medianV*(1.d0 - s_prime/t_total)
            else
               exp_func = exp((s_prime/t_total)*log(medianU))
               futrU = medianU/exp_func
               futrV = medianV*(medianU - exp_func)/(medianU - 1.d0)
            end if

         end associate
         end subroutine estimate_future_integrals

      end subroutine walk



      function dist(r1, r2)
      class(space), intent(in) :: r1, r2
      real(8) :: dist
      dist = sqrt((r1%x1 - r2%x1)**2 + (r1%x2 - r2%x2)**2 + (r1%x3 - r2%x3)**2)
      end function dist



      function atss(dx, na, A, B)
!     adaptive temporal stepsize
      real(8), intent(in) :: dx
      integer, intent(in) :: na
      real(8) :: A(na,na), B(na)
      real(8) :: atss
      real(8) :: normB, trA, dt_adv, dt_dif

      normB = sqrt(dot_product(B, B))
      if( normB>0.d0 ) then
         dt_adv = dx/normB
      else
         dt_adv = dtMax
      end if

      trA = trace(na, A)
      if( trA>0.d0 ) then
         dt_dif = dx*dx/trA
      else
         dt_dif = dtMax
      end if

      atss = min(dt_adv, dt_dif)
      atss = min(atss, dtMax)
      atss = max(atss, dtMin)

      end function atss



      subroutine third_type_bc(ito_cmpnt, condition, vgamma, lambda)
!     this routine returns coefficients of the normalized homogeneous third type
!     boundary condition
!
!        vgamma.grad u - lambda * u = 0
!
!     in which vgamma is a unit conormal vector
!
!        vgamma = n.A/|n.A|
!
!     and
!
!        lambda = 2*(eta + g2)/|n.A|
!
!     where eta (= n.h) and g2 are defined in comments of the source file
!     mod_domain_typedef.F90. the matrix A and vector h are defined as in
!     source file mod_typedef_params.F90.
      type(components), intent(in) :: ito_cmpnt
      type(boundary_condition), intent(in) :: condition
      real(8), intent(out) :: vgamma(3)
      real(8), intent(out) :: lambda
!     the smallest number that the code could distinguish from zero
      real(8), parameter :: eps = epsilon(lambda)
      real(8) :: n_dot_a(3)
      real(8) :: norm_n_dot_a
      associate( g=>condition%g, n=>condition%n, A=>ito_cmpnt%A, h=>ito_cmpnt%h )

      n_dot_a = matmul(n, A)
      norm_n_dot_a = sqrt(dot_product(n_dot_a, n_dot_a))

      if( norm_n_dot_a>eps ) then
         vgamma = n_dot_a/norm_n_dot_a
         lambda = 2.d0*(dot_product(n, h) + g)/norm_n_dot_a
      else
!        |A| is pathetically small here. set lambda = 0 so that no flow is allowed
!        to cross the boundary
         vgamma = n
         lambda = 0.d0
      end if

      end associate
      end subroutine third_type_bc



      subroutine projection(xt_curr, bid, vgamma, xt_proj, tp, success)
!     given a boundary id bid, a spacetime position xt_curr (assumed in vicinity
!     of the boundary), and the local gamma vector vgamma, this routine tries to
!     find the projection point xt_proj of xt_curr along vgamma on this boundary,
!     and returns a logical flag to indicate if successful.
      type(spacetime), intent(in) :: xt_curr
      type(spacetime), intent(out) :: xt_proj
      integer, intent(in) :: bid
!     vgamma*tp gives the displacement vector from xt_curr to the boundary along
!     vgamma. since vgamma is inward pointing, tp > 0 if xt_curr is off the domain
!     and tp <=0 if xt_curr is on the domain.
      real(8), intent(in) :: vgamma(3)
      real(8), intent(out) :: tp
      logical, intent(out) :: success
      real(8) :: tp1, tp2, tpl, tpu
      type(boundary), pointer :: bptr => null()
!     interface func_real is provided in interpolation module
      procedure(func_real), pointer :: dptr => null()
      !$OMP threadprivate(bptr, dptr)

      bptr => boundaries(bid)
      dptr => distance

!     find the root of distance(tp) = 0, where distance is from the point
!     parameterized by tp to the boundary bid. first, bracket the root between
!     lower and upper limits tpl and tpu.
      tpl = -10.d0*Neighborhood
      tpu = 10.d0*Neighborhood
!     tp1 and tp2 are initial guesses of the range -- they will be updated.
!     the following tp1 and tp2 values mean that we assume xt_curr is on the domain
!     and the projection point is within a distance of Neighborhood along vgamma.
!     this is expected to be the most probable case.
      tp1 = -Neighborhood
      tp2 = 0.d0
!     bracket the root of dptr(tp) = 0 inside the updated interval [tp1, tp2]
      call zbrac(dptr, tp1, tp2, tpl, tpu, success)
      if( success ) then
!        root is successfully bracketed. find the root.
!        note that tp has dimension of distance since vgamma is a unit vector. the
!        rigorous rtbis accuracy for tp (tacc) to find a root within a thickness
!        of epsBnd near boundary is tacc = epsBnd/(n.vgamma), where the
!        denominator is a cosine function and is always <= 1. to avoid using the
!        boundary normal vector n, set tacc = epsBnd, and it also gives a more
!        accurate root.
         tp = rtbis(dptr, tp1, tp2, epsBnd)
         call tpoint(tp, xt_proj)
      else
!        failed to find the projection point
         tp = 0.d0
         xt_proj = xt_curr
      end if

      nullify(bptr, dptr)
      contains

         subroutine tpoint(tp, xt_tp)
         real(8), intent(in) :: tp
         type(spacetime), intent(out) :: xt_tp

         xt_tp%t = xt_curr%t
         xt_tp%x1 = xt_curr%x1 + vgamma(1)*tp
         xt_tp%x2 = xt_curr%x2 + vgamma(2)*tp
         xt_tp%x3 = xt_curr%x3 + vgamma(3)*tp

         end subroutine tpoint


         real(8) function distance(tp)
         real(8), intent(in) :: tp
         type(spacetime) :: xt_tp
         call tpoint(tp, xt_tp)
         distance = bptr%equation(xt_tp)
         end function distance

      end subroutine projection



      subroutine intersection(xt_on, xt_off, bid, xt_int, tp)
!     this routine calculates the intersection point xt_int between the boundary
!     bid and a line segment crossing it from xt_on (on the domain) to xt_off (off
!     the domain). tp is a number in [0,1] parameterizing xt_int.
      type(spacetime), intent(in) :: xt_on, xt_off
      integer, intent(in) :: bid
      type(spacetime), intent(out) :: xt_int
!     tp = 0 => xt_int == xt_on
!     tp = 1 => xt_int == xt_off
      real(8), intent(out) :: tp
      real(8) :: equation_x_on, equation_x_off, tacc
      type(boundary), pointer :: bptr => null()
!     interface func_real is provided in interpolation module
      procedure(func_real), pointer :: dptr => null()
      !$OMP threadprivate(bptr, dptr)

      bptr => boundaries(bid)
      dptr => distance

      equation_x_on = bptr%equation(xt_on)
      equation_x_off = bptr%equation(xt_off)

      if( abs(equation_x_on)<epsBnd .or. &
        & abs(equation_x_off - equation_x_on)<epsBnd ) then
!        xt_on is on the boundary, or the line segement is tangent to the boundary.
!        regard xt_on as the intersection point
         tp = 0.d0
         xt_int = xt_on
      elseif( abs(equation_x_off)<epsBnd ) then
!        xt_off is on the boundary. regard xt_off as the intersection point.
         tp = 1.d0
         xt_int = xt_off
      else
!        the demoninator is guaranteed greater than epsBnd by the first if clause
         tacc = epsBnd/abs(dptr(1.d0) - dptr(0.d0))
         tp = rtbis(dptr, 0.d0, 1.d0, tacc)
         call tpoint(tp, xt_int)
      end if

      nullify(bptr, dptr)
      contains

         subroutine tpoint(tp, xt_tp)
         real(8), intent(in) :: tp
         type(spacetime), intent(out) :: xt_tp

         xt_tp%t = xt_on%t
         xt_tp%x1 = (1.d0 - tp)*xt_on%x1 + tp*xt_off%x1
         xt_tp%x2 = (1.d0 - tp)*xt_on%x2 + tp*xt_off%x2
         xt_tp%x3 = (1.d0 - tp)*xt_on%x3 + tp*xt_off%x3

         end subroutine tpoint


         real(8) function distance(tp)
         real(8), intent(in) :: tp
         type(spacetime) :: xt_tp
         call tpoint(tp, xt_tp)
         distance = bptr%equation(xt_tp)
         end function distance

      end subroutine intersection

end module mod_solver_ito_process

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) renamed unwanted module entities to dummy_(original name) in the 'use' clauses
!
! 2) added reference to Oksendal [2000]
!
! 3) made relevant adaptions to the new type definitions, domain, grid and equation
!    modules in this revision
!
! 4) added a path integral about the v term in Eq. (3) in mod_equation_typedef.F90
!
! 5) did not use the staggered-time-and-space sampling of SDE coefficients in the
!    modified Euler scheme of Ref. [3], but used uniform-time-and-space sampling,
!    for the following reasons: (a) the staggered scheme is designed to deal with
!    time-variation irregularity in the coefficients, but does not improve accuracy
!    of the scheme; (b) the staggered scheme complicates the algorithm and slows
!    down calculation by requiring extra evaluations of the coefficients at
!    different times; (c) most of the coefficients met by this program are expected
!    to be time-independent, which make the staggered scheme indifferent from the
!    uniform one, or regularly varied in time in the worst case. 
!
! 6) changed the use of parameter CriticalRatio to an input argument RATIO to 
!    adapt to the improved fission strategy (see mod_solver_batch.F90).
!
! 7) changed type-bound procedure name RUN to WALK to better reflect the meaning
!    of a random walk.
!
! Liheng Zheng, 12/29/2019
!
! 1) corrected the fission criterion and the relevant input/output arguments of
!    procedure WALK per the corrected tree tracing method. previously, partial
!    functional value was used in the fission criterion; the corrected criterion
!    uses the projected functional value of the entire path instead.
!
! 2) modified and moved type stat_moments from mod_solver_batch.F90 to here
!
! Liheng Zheng, 02/25/2020
!
! 1) modified type stat_moments and moved back to mod_solver_batch.F90
!
! 2) change the fission criterion: it no longer uses the statistics of batch means
!    but the statistics of ito prcoesses themselves. modified type ito_proc_params
!    in accordance.
!
! Liheng Zheng, 03/04/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *







