module ice_grid_neighbor

    use ice_kinds_mod
    use ice_blocks, only: nx_block, ny_block
    use ice_exit, only: abort_ice

    implicit none
    private
    public :: grid_neighbor_min, grid_neighbor_max

    contains

!=======================================================================
! Compute the minimum of adjacent values of a field at specific indices,
! depending on the grid location (U, E, N)
!
    real(kind=dbl_kind) function grid_neighbor_min(field, i, j, grid_location) result(mini)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         field    ! field defined at T point

      integer (kind=int_kind), intent(in) :: &
         i, j

      character(len=*), intent(in) :: &
         grid_location ! grid location at which to compute the minumum (U, E, N)

      character(len=*), parameter :: subname = '(grid_neighbor_min)'

      select case (trim(grid_location))
         case('U')
            mini = min(field(i,j), field(i+1,j), field(i,j+1), field(i+1,j+1))
         case('E')
            mini = min(field(i,j), field(i+1,j))
         case('N')
            mini = min(field(i,j), field(i,j+1))
         case default
            call abort_ice(subname // ' unknown grid_location: ' // grid_location, file=__FILE__, line=__LINE__)
      end select

    end function grid_neighbor_min

!=======================================================================
! Compute the maximum of adjacent values of a field at specific indices,
! depending on the grid location (U, E, N)
!
    real(kind=dbl_kind) function grid_neighbor_max(field, i, j, grid_location) result(maxi)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         field    ! field defined at T point

      integer (kind=int_kind), intent(in) :: &
         i, j

      character(len=*), intent(in) :: &
         grid_location ! grid location at which to compute the maximum (U, E, N)


      character(len=*), parameter :: subname = '(grid_neighbor_max)'

      select case (trim(grid_location))
         case('U')
            maxi = max(field(i,j), field(i+1,j), field(i,j+1), field(i+1,j+1))
         case('E')
            maxi = max(field(i,j), field(i+1,j))
         case('N')
            maxi = max(field(i,j), field(i,j+1))
         case default
            call abort_ice(subname // ' unknown grid_location: ' // grid_location, file=__FILE__, line=__LINE__)
      end select

    end function grid_neighbor_max

end module ice_grid_neighbor