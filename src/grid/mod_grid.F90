module mod_grid
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

! * * * * * * * * * * * * * * * * * Editing Log * * * * * * * * * * * * * * * * * *
! created by Liheng Zheng on 11/18/2019
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

