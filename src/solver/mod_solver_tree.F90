module mod_solver_tree
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
      use mod_typedef_params, only: spacetime, nChildren, nGenerations
      implicit none
      private
      public tree
! A tree as exemplified below (nGenerations=3 and nChildren=2)
!
!               +---2a
!               |
!         +---1a+
!         |     |
!         |     +---2b
!         |
!     ---0+           +---3e
!         |           |
!         |     +---2c+
!         |     |     |
!         +---1b+     +---3f
!               |
!               +---2d
!
! is stored in a path table like this
!     
!     +---+---+---+---+
!     | 0 |1a |2a |///|
!     +---+---+---+---+
!     | 0 |1b |2c |3e |
!     +---+---+---+---+
!     | 0 |1a |2b |///|
!     +---+---+---+---+
!     | 0 |1b |2d |///|
!     +---+---+---+---+
!     | 0 |1b |2c |3f |
!     +---+---+---+---+
!     |///|///|///|///|
!     +---+---+---+---+
!     |///|///|///|///|
!     +---+---+---+---+
!     |///|///|///|///|
!     +---+---+---+---+
!
! The table has nChildren**nGenerations rows and nGenerations+1 columns. In this
! case, npaths = 5, and all paths below the 5th are empty.
!
! Each node of the tree corresponds to a segment of a stochastic path, and stores
! the ending spacetime of that path segment (except for the terminal node, see
! below), the aggregated path integrals U and V from the root down to that node
! (c.f. mod_solver_ito_process.F90 for expressions of the path integrals), and two
! auxiliary logicals TERMINAL and EMPTY to indicate the status of the node. For the
! terminal node, the initial or the Dirichlet type boundary condition at the ending
! spacetime is stored in the first entry of the spacetime variable, with the rest
! entries unused. For a node in the j-th (j>0) generation, U and V are aggregated by
!
!    U_j = U_(j-1) * u_j,                                               (1)
!    V_j = V_(j-1) + U_(j-1) * v_j,                                     (2)
!
! where U_j and V_j are the aggregated path integrals stored in that node, and u_j
! and v_j are the path integrals evaluated from the path segment corresponding to
! the node.
!
! The value of the entire tree (its trace) is obtained by a weighted sum over the
! values of tree branches. The value of a branch is the funtional value of the
! stochastic path from root all the way to the terminus, which is calculated from
! the initial/boundary value (g) and the aggragated path integrals by
!
!    F = g*U + V.                                                       (3)

      type node
         real(8) :: U
         real(8) :: V
         type(spacetime) :: xt
         logical :: terminal
         logical :: empty
      end type node

      type tree
         private
         integer :: npaths
         integer :: ipath
         integer :: jgene
         type(node), allocatable :: table(:,:)
      contains
         procedure :: set_node_value
         procedure :: get_node_value
         procedure :: move_to_next_path
         procedure, private :: split
         procedure :: fission
         procedure :: trace
         procedure :: dump
         procedure :: destroy
!         final :: finish
      end type tree

      interface tree
         module procedure start
      end interface tree

contains

      function start()
      type(tree) :: start
      type(node) :: empty_node

      start%npaths = 1
      start%ipath = 1
      start%jgene = 0

      empty_node%U = 1.d0
      empty_node%V = 0.d0
      empty_node%xt = spacetime(0.d0, 0.d0, 0.d0, 0.d0)
      empty_node%terminal = .false.
      empty_node%empty = .true.

      allocate( start%table(1:nChildren**nGenerations, 0:nGenerations) )
      start%table = empty_node

      end function start



      subroutine set_node_value(this, u, v, xt, terminal)
!     set values to the current (ipath, jgene) node
      class(tree) :: this
      real(8), intent(in), optional :: u, v
      type(spacetime), intent(in), optional :: xt
      logical, intent(in), optional :: terminal
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )

         if( i>=1 .and. i<=n .and. j>=0 .and. j<=nGenerations ) then
            if( present(u) ) table(i,j)%U = u
            if( present(v) ) table(i,j)%V = v
            if( present(xt) ) table(i,j)%xt = xt
            if( present(terminal) ) table(i,j)%terminal = terminal
            table(i,j)%empty = .false.
         else
            print *,' mod_solver_tree procedure set_node_value:'
            print *,'  ipath or jgene out of bound'
            print *,'  ipath =',i,', npaths =',n
            print *,'  jgene =',j,', nGenerations =',nGenerations
            stop
         end if

      end associate
      end subroutine set_node_value



      subroutine get_node_value(this, u, v, xt, terminal)
!     get values from the current (ipath, jgene) node
      class(tree) :: this
      real(8), optional :: u, v
      type(spacetime), optional :: xt
      logical, optional :: terminal
      associate( i=>this%ipath, j=>this%jgene, table=>this%table )

         if( .not.table(i,j)%empty ) then
            if( present(u) ) u = table(i,j)%U
            if( present(v) ) v = table(i,j)%V
            if( present(xt) ) xt = table(i,j)%xt
            if( present(terminal) ) terminal = table(i,j)%terminal
         else
            print *,' mod_solver_tree procedure get_node_value:'
            stop '  trying to get values from an empty node'
         end if

      end associate
      end subroutine get_node_value



      function move_to_next_path(this)
!     this function moves the cursor ipath to the next path in the current
!     generation, and returns T if the operation was successful; otherwise, the
!     cursor ipath is moved to the first non-empty path in the current generation
!     and F is returned.
      class(tree) :: this
      logical :: move_to_next_path

      move_to_next_path = .true.
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )
         do while( .true. )
            i = i + 1
            if( i>n ) then
               move_to_next_path = .false.
               exit
            end if
            if( .not.table(i,j)%empty ) exit
         end do
!        if no 'next path' exists, move ipath to the first non-empty path in this
!        generation
         if( .not.move_to_next_path ) then
            i = 1
            do while( i<=n )
               if( .not.table(i,j)%empty ) exit
               i = i + 1
            end do
         end if
      end associate

      end function move_to_next_path



      subroutine split(this)
!     this routine splits one path into nChildren paths at the node indicated by the
!     cursors ipath and jgene. In doing so, the parent node is copied to its first
!     child along the path, and then the path is copied (nChildren - 1) times after
!     the existing paths, so that npaths is updated to (npaths + nChildren - 1)
!     afterwards. note that, after splitting, the cursors remain at the parent node.
      class(tree) :: this
      integer :: k
      type(spacetime) :: xt
      real(8) :: u, v
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )

         if( i>=1 .and. i<=n .and. j<nGenerations ) then
!           copy the content at the current node to its first child
            call this%get_node_value(u=u, v=v, xt=xt)
            j = j + 1
            call this%set_node_value(u=u, v=v, xt=xt)
            j = j - 1
!           copy the current path nChildren-1 times to generate the other children
            do k=1,nChildren-1
               table(n+k,:) = table(i,:)
            end do
!           update npaths
            n = n + nChildren - 1
         else
            print *,' mod_solver_tree procedure split:'
            print *,'  ipath or jgene out of bound'
            print *,'  ipath =',i,', npaths =',n
            print *,'  jgene =',j,', nGenerations =',nGenerations
            stop
         end if

      end associate
      end subroutine split



      function fission(this)
!     this routine splits all fissile paths in the current generation, and moves
!     the cursors to the first path in the latest generation. this function has
!     three integer exit status:
!        0 - the path tree had not grown (no fissile path exists)
!        1 - the path tree had grown but cannot further grow (reached the last
!            generation)
!        2 - the path tree had grown and can still grow
      class(tree) :: this
      integer :: fission
      integer :: old_npaths

      fission = 0
      old_npaths = this%npaths
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )
         if( j<nGenerations ) then
            do i=1,old_npaths
               if( .not.table(i,j)%empty .and. .not.table(i,j)%terminal ) then
                  call split(this)
                  fission = 1
               end if
            end do
         end if
!        if the path tree grew, move jgene to the next generation
         if( fission>0 ) then
            j = j + 1
            if( j<nGenerations ) fission = 2
         end if
!        move ipath to the first non-empty node
         i = 1
         do while( i<=n )
            if( .not.table(i,j)%empty ) exit
            i = i + 1
         end do
      end associate

      end function fission



      subroutine trace(this, trF, trU, trV)
!     this routine calculates the tree traces of the functional value F, the path
!     integral U and the path integral V, by weighted sums over branch values.
      class(tree) :: this
      real(8), intent(out) :: trF, trU, trV
      real(8) :: weight
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )

         trF = 0.d0
         trU = 0.d0
         trV = 0.d0
         do i=1,n
            do j=0,nGenerations
               if( table(i,j)%terminal ) then
                  weight = 1.d0/dble(nChildren**j)
!                 at the terminal node, xt is no longer used to store the ending
!                 spacetime; instead, xt%x1 stores the initial or Dirichlet boundary
!                 condition at the ending spacetime
                  trF = trF + (table(i,j)%xt%x1*table(i,j)%U + table(i,j)%V)*weight
                  trU = trU + table(i,j)%U*weight
                  trV = trV + table(i,j)%V*weight
                  exit
               end if
            end do

            if( j>nGenerations ) then
               print *,' mod_solver_tree procedure trace:'
               print *,'  missing terminal node for branch i =',i
               print *,'  the tree:'
               call this%dump()
               stop
            end if        
         end do

      end associate
      end subroutine trace



      subroutine dump(this)
!     print all paths in the tree on screen
      class(tree) :: this

100   format(A1,2F6.3,A2,3(F6.3,A1),F6.3,A2)
      associate( i=>this%ipath, j=>this%jgene, n=>this%npaths, table=>this%table )
         do i=1,n
            write(6,'(A5,I3,A2)',advance='no') ' Path',i,': '
            do j=0,nGenerations
               write(6,fmt=100,advance='no') '[',table(i,j)%U, table(i,j)%V,',(' &
                                                ,table(i,j)%xt%x1,',' &
                                                ,table(i,j)%xt%x2,',' &
                                                ,table(i,j)%xt%x3,',' &
                                                ,table(i,j)%xt%t,')]'
               if( j==nGenerations .or. table(i,j+1)%empty ) then
                  exit
               else
                  write(6,'(A4)',advance='no') ' -> '
               end if
            end do
            write(6,'(A)') ': root'
         end do
      end associate

      end subroutine dump



      subroutine destroy(this)
      class(tree) :: this
      if( allocated(this%table) ) deallocate( this%table )
      end subroutine destroy

!
!
!      subroutine finish(this)
!      type(tree) :: this
!      if( allocated(this%table) ) deallocate( this%table )
!      end subroutine finish
!
end module mod_solver_tree

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! 1) added a second real field into type node to contain the extra path-integral
!    functional value introduced by the term v in Eq. (2) of mod_equation_
!    typedef.F90, and changed the relevant procedures
!
! 2) renamed the original type-bound procedure FINISH to DESTROY which is used to
!    forcibly destroy a tree object to release memory, and added a final procedure
!    named FINISH. the procedure DESTROY is used in the do-loop of SOLVER_BATCH,
!    where after each iteration the tree object is no longer needed but the compiler
!    would not know.
!
! Liheng Zheng, 11/19/2019
!
! 1) commented out the final subroutine FINISH to avoid duplicated memory-release
!    attempts since DESTROY is used explicitly
!
! Liheng Zheng, 01/29/2020
!
! 1) corrected the tracing method and the underneath meaning of node values U and V.
!    previously, nodal U and V store the segmental path integrals of the
!    corresponding stochastic path segment, and the tracing method was only correct
!    when V vanishes. in this correction, nodal U and V store the aggregated path
!    integrals from the root down to that node. in addition, for the terminal node,
!    the initial/boundary condition value is no longer multiplied to U, but is
!    stored in xt%x1 since xt is useless for the terminal node. this allows TRACE
!    to also calculate the traced U and V values of the tree, which are useful in
!    determining the fission criterion in mod_solver_ito_process.F90.
!
! Liheng Zheng, 02/25/2020
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


