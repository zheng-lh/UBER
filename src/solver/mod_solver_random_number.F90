module mod_solver_random_number
!     this module sets up the pseudo-random number generator for the solver.
      use iso_c_binding
      use omp_lib
      implicit none
      private
      public init_dcmt, free_dcmt, gausrand, exprand
      save

!     thread identity
      integer(c_int) :: tid
      !$OMP threadprivate(tid)

!     procedures in solver_dcmt_omp.c
      external get_mt_parameter_id_st_f
      external sgenrand_mt_f
      external free_mt_struct_f
      real(c_double), external :: gausrand
      real(c_double), external :: exprand

contains

      integer(c_int) function gseed(id)
!     Returns the (C int) seed given an arbitrary (C int) integer id. Note that, on
!     64-bit machines, C int has 32 bits (range: -2^31+1 to +2^31-1).
      implicit none
      integer(c_int), intent(in) :: id
      integer, dimension(8) :: digit
      character(12) :: date, time, zone
      call date_and_time(date,time,zone,digit)
      gseed = digit(8)*10000 + digit(7)*100 + digit(6) + id
      return
      end function gseed
      
      
      
      subroutine init_dcmt()
!     Initialize the DCMT pseudo-random number generator. Two steps: 1. find an
!     mt_struct for each thread; 2. in each thread, generate a seed, and use it to
!     initialize the PRNG.
      implicit none
      integer(c_int) :: seed
      integer :: nt, np
      
      print '(/A)',' Initializing DCMT pseudo-random number generator...'
      
      !$OMP parallel default(private) shared(nt, np)
         !$OMP master
            nt = omp_get_num_threads()
            np = omp_get_num_procs()
            print '(A,I3)',' Number of threads used: ',nt
            print '(A,I3)',' Number of processors used: ',np
         !$OMP end master
         tid = omp_get_thread_num()
         call get_mt_parameter_id_st_f(tid)
         seed = gseed(tid)
         call sgenrand_mt_f(seed)
      !$OMP end parallel

      print *,''      
      return
      end subroutine init_dcmt
      
      
      
      subroutine free_dcmt()
!     Release the memory claimed by the DCMT PRNG.
      implicit none
      !$OMP parallel
         call free_mt_struct_f()
      !$OMP end parallel
      print '(/A/)',' DCMT pseudo-random number generator released.'
      return
      end subroutine free_dcmt

end module mod_solver_random_number




