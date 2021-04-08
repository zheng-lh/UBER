module mod_grid
! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 11/18/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use mod_grid_scalar_field
      use mod_grid_vector_field
      use mod_grid_tensor_field
      use mod_grid_io
      implicit none
      public
!
! This module is the interface between the grid objects and the SDE solver, so that
! the solver needs only 'use mod_grid' rather than looking into the concrete
! mod_grid_* modules. This interface provides access to:
!
!     type grid_scalar_field
!     type grid_vector_field
!     type grid_tensor_field
!     subroutine grid_io_read_data_file
!
end module mod_grid


