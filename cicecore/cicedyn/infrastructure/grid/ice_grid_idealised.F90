module ice_grid_idealised

    use ice_kinds_mod
    use ice_boundary, only: ice_HaloExtrapolate
    use ice_communicate, only: my_task, master_task
    use ice_constants, only: c0, c1, c2, c20, cm_to_m, radius, &
         field_loc_NEcorner, field_type_scalar, field_loc_center
    use ice_fileunits, only: nu_diag
    use ice_exit, only: abort_ice
    use ice_domain_size, only: nx_global, ny_global
    use ice_domain, only: distrb_info, &
        ew_boundary_type, ns_boundary_type, close_boundaries
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_warnings_aborted
    use ice_gather_scatter, only: scatter_global
    use ice_grid_lengths, only: primary_grid_lengths_HTN, primary_grid_lengths_HTE
    
    implicit none
    private
    public :: grid_boxislands_kmt, rectgrid, rectgrid_scale_dxdy

    contains
!=======================================================================

! Regular rectangular grid and mask
!
! author: Elizabeth C. Hunke, LANL

    subroutine rectgrid(save_ghte_ghtn, scale_dxdy, dxrect, dyrect, &
        G_HTE, G_HTN, kmt_type, lonrefrect, latrefrect, hm, kmt, angle, &
        ULON, ULAT, dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN )

        logical (kind=log_kind), intent(in) :: &
            save_ghte_ghtn, & ! flag for saving global hte and htn during initialization
            scale_dxdy        ! flag to apply scale factor to vary dx/dy in rectgrid

        character (len=char_len_long), intent(in)  :: &
            kmt_type  !  options are file, default, boxislands

        real (kind=dbl_kind), intent(in) :: &
            lonrefrect, & ! lower left lon for rectgrid
            latrefrect    ! lower left lat for rectgrid

        ! grid dimensions for rectangular grid
        real (kind=dbl_kind), intent(in)::  &
            dxrect, & !  user_specified spacing (cm) in x-direction (uniform HTN)
            dyrect    !  user_specified spacing (cm) in y-direction (uniform HTE)

        real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            hm     , & ! land/boundary mask, thickness (T-cell)
            kmt      , &   ! ocean topography mask for bathymetry (T-cell)
            angle  ! for conversions between POP grid and lat/lon

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
            ULAT       ! latitude of velocity pts, NE corner of T pts (radians)

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            dxT    , & ! width of T-cell through the middle (m)
            dyT    , & ! height of T-cell through the middle (m)
            dxU    , & ! width of U-cell through the middle (m)
            dyU    , & ! height of U-cell through the middle (m)
            dxN    , & ! width of N-cell through the middle (m)
            dyN    , & ! height of N-cell through the middle (m)
            dxE    , & ! width of E-cell through the middle (m)
            dyE    , & ! height of E-cell through the middle (m)
            HTE    , & ! length of eastern edge of T-cell (m)
            HTN        ! length of northern edge of T-cell (m)
        
        !local vars

        integer (kind=int_kind) :: &
           i, j, &
           imid, jmid
  
        real (kind=dbl_kind) :: &
           length,  &
           rad_to_deg
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g1
  
        character(len=*), parameter :: subname = '(rectgrid)'
  
        !-----------------------------------------------------------------
        ! Calculate various geometric 2d arrays
        !-----------------------------------------------------------------
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        hm (:,:,:) = c0
        kmt(:,:,:) = c0
        angle(:,:,:) = c0   ! "square with the world"
  
        allocate(work_g1(nx_global,ny_global))
  
        if (my_task == master_task) then
            work_g1 = c0
            length = dxrect*cm_to_m/radius*rad_to_deg

            work_g1(1,:) = lonrefrect ! reference lon from namelist

            do j = 1, ny_global
            do i = 2, nx_global
                work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
            enddo
            enddo
            work_g1(:,:) = work_g1(:,:) / rad_to_deg
        endif
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULON, distrb_info, &
                                ew_boundary_type, ns_boundary_type)

        if (my_task == master_task) then
            work_g1 = c0
            length = dyrect*cm_to_m/radius*rad_to_deg

            work_g1(:,1) = latrefrect ! reference latitude from namelist

            do i = 1, nx_global
            do j = 2, ny_global
                work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
            enddo
            enddo
            work_g1(:,:) = work_g1(:,:) / rad_to_deg
        endif
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULAT, distrb_info, &
                                ew_boundary_type, ns_boundary_type)

        if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
                work_g1(i,j) = dxrect             ! HTN
            enddo
            enddo
        endif
        call primary_grid_lengths_HTN(work_g1, save_ghte_ghtn, HTN, &
            dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN)  

        if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
                work_g1(i,j) = dyrect             ! HTE
            enddo
            enddo
        endif
        call primary_grid_lengths_HTE(work_g1, save_ghte_ghtn, HTN, &
            dyT, dyU, dyN, dyE, G_HTE, G_HTN) 
    
        !-----------------------------------------------------------------
        ! Construct T-cell land mask
        ! Keyed on ew_boundary_type; ns_boundary_type should be 'open'.
        !-----------------------------------------------------------------
  
        if (my_task == master_task) then
           work_g1(:,:) = c0      ! initialize hm as land
  
           if (trim(kmt_type) == 'boxislands') then
  
              call grid_boxislands_kmt(work_g1)
  
           elseif (trim(kmt_type) == 'channel') then
  
              do j = 3,ny_global-2     ! closed top and bottom
              do i = 1,nx_global       ! open sides
                 work_g1(i,j) = c1     ! NOTE nx_global > 5
              enddo
              enddo
  
           elseif (trim(kmt_type) == 'channel_oneeast') then
  
              do j = ny_global/2,ny_global/2    ! one channel wide
              do i = 1,nx_global       ! open sides
                 work_g1(i,j) = c1     ! NOTE nx_global > 5
              enddo
              enddo
  
           elseif (trim(kmt_type) == 'channel_onenorth') then
  
              do j = 1,ny_global       ! open sides
              do i = nx_global/2,nx_global/2    ! one channel wide
                 work_g1(i,j) = c1     ! NOTE nx_global > 5
              enddo
              enddo
  
           elseif (trim(kmt_type) == 'wall') then
  
              do j = 1,ny_global       ! open except
              do i = 1,nx_global-2     ! closed east edge
                 work_g1(i,j) = c1
              enddo
              enddo
  
           elseif (trim(kmt_type) == 'default') then
  
              ! land in the upper left and lower right corners,
              ! otherwise open boundaries
              imid = nint(aint(real(nx_global)/c2))
              jmid = nint(aint(real(ny_global)/c2))
  
              do j = 3,ny_global-2
              do i = 3,nx_global-2
                 work_g1(i,j) = c1    ! open central domain
              enddo
              enddo
  
              if (nx_global > 5 .and. ny_global > 5) then
  
                 do j = 1, jmid+2
                 do i = 1, imid+2
                    work_g1(i,j) = c1    ! open lower left corner
                 enddo
                 enddo
  
                 do j = max(jmid-2,1), ny_global
                 do i = max(imid-2,1), nx_global
                    work_g1(i,j) = c1    ! open upper right corner
                 enddo
                 enddo
  
              endif ! > 5x5 grid
  
           else
  
              call abort_ice(subname//' ERROR: unknown kmt_type '//trim(kmt_type), &
                   file=__FILE__, line=__LINE__)
  
           endif ! kmt_type
  
           if (close_boundaries) then
              work_g1(:, 1:2) = c0
              work_g1(:, ny_global-1:ny_global) = c0
              work_g1(1:2, :) = c0
              work_g1(nx_global-1:nx_global, :) = c0
           endif
  
        endif
  
        call scatter_global(hm, work_g1, master_task, distrb_info, &
                            field_loc_center, field_type_scalar)
  
        deallocate(work_g1)
  
    end subroutine rectgrid
  
  !=======================================================================
        ! generate a variable spaced rectangluar grid.
        ! extend spacing from center of grid outward.
        
    subroutine rectgrid_scale_dxdy(save_ghte_ghtn, scale_dxdy, dxrect, dyrect, &
            dxscale, dyscale, lonrefrect, latrefrect, &
            G_HTE, G_HTN, hm, kmt, angle, ULON, ULAT, &
            dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN )

        logical (kind=log_kind), intent(in) :: &
            save_ghte_ghtn, & ! flag for saving global hte and htn during initialization
            scale_dxdy        ! flag to apply scale factor to vary dx/dy in rectgrid

        ! grid dimensions for rectangular grid
        real (kind=dbl_kind), intent(in)::  &
            dxrect, & !  user_specified spacing (cm) in x-direction (uniform HTN)
            dyrect    !  user_specified spacing (cm) in y-direction (uniform HTE)

        ! growth factor for variable spaced grid
        real (kind=dbl_kind), intent(in)::  &
            dxscale, & !  scale factor for grid spacing in x direction (e.g., 1.02)
            dyscale    !  scale factor for gird spacing in y direction (e.g., 1.02)

        real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)

        real (kind=dbl_kind), intent(in) :: &
            lonrefrect, & ! lower left lon for rectgrid
            latrefrect    ! lower left lat for rectgrid

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            hm     , & ! land/boundary mask, thickness (T-cell)
            kmt      , &   ! ocean topography mask for bathymetry (T-cell)
            angle  ! for conversions between POP grid and lat/lon

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
            ULAT       ! latitude of velocity pts, NE corner of T pts (radians)

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            dxT    , & ! width of T-cell through the middle (m)
            dyT    , & ! height of T-cell through the middle (m)
            dxU    , & ! width of U-cell through the middle (m)
            dyU    , & ! height of U-cell through the middle (m)
            dxN    , & ! width of N-cell through the middle (m)
            dyN    , & ! height of N-cell through the middle (m)
            dxE    , & ! width of E-cell through the middle (m)
            dyE    , & ! height of E-cell through the middle (m)
            HTE    , & ! length of eastern edge of T-cell (m)
            HTN        ! length of northern edge of T-cell (m)

        !local vars

        integer (kind=int_kind) :: &
            i, j, &
            center1, center2 ! array centers for expanding dx, dy

        real (kind=dbl_kind) :: &
            length,  &
            rad_to_deg

        real (kind=dbl_kind), dimension(:,:), allocatable :: &
            work_g1

        character(len=*), parameter :: subname = '(rectgrid_scale_dxdy)'

        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        hm (:,:,:) = c0
        kmt(:,:,:) = c0
        angle(:,:,:) = c0   ! "square with the world"
  
        allocate(work_g1(nx_global,ny_global))

        ! determine dx spacing
        ! strategy: initialize with dxrect.
        ! if want to scale the grid, work from center outwards,
        ! multplying neighbor cell by scale factor.
        ! this assumes dx varies in x direction only.
        ! (i.e, dx is the same across same y location)
        if (my_task == master_task) then

            ! initialize with initial dxrect
            work_g1(:,:) = dxrect

            ! check if nx is even or odd
            ! if even, middle 2 columns are center
            ! of odd,  middle 1 column is center
            if (mod(nx_global,2) == 0) then ! nx_global is even

            ! with even number of x locatons,
            ! the center two y columns are center
            center1 = nx_global/2  ! integer math
            center2 = center1 + 1  ! integer math

            else ! nx_global = odd
            ! only one center index. set center2=center1
            center1 = ceiling(real(nx_global/2),int_kind)
            center2 = center1
            endif

            ! note loop over only half the x grid points (center1)-1
            ! working from the center outward.
            do j = 1, ny_global
            do i = 1, center1-1
            ! work from center1 to left
            work_g1(center1-i,j) = dxscale*work_g1(center1-i+1,j)

            ! work from center2 to right
            work_g1(center2+i,j) = dxscale*work_g1(center2+i-1,j)
            enddo ! i
            enddo ! j

        endif       ! my_task == master_task


        ! note work_g1 is converted to meters in primary_grid_lengths_HTN
        call primary_grid_lengths_HTN(work_g1, save_ghte_ghtn, HTN, & 
            dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN)  

        ! make ULON array
        if (my_task == master_task) then

            ! make first column reference lon in radians.
            ! the remaining work_g1 is still dx in meters
            work_g1(1,:) = lonrefrect/rad_to_deg ! radians

            ! loop over remaining points and add spacing to successive
            ! x locations
            do j = 1, ny_global
            do i = 2, nx_global ! start from i=2. i=1 is lonrefrect
            length = work_g1(i,j)/radius             ! grid spacing in radians
            work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
            enddo ! i
            enddo ! j
        endif    ! mytask == master_task
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULON, distrb_info, &
                                ew_boundary_type, ns_boundary_type)

        ! determine dy spacing
        ! strategy: initialize with dyrect.
        ! if want to scale the grid, work from center outwards,
        ! multplying neighbor cell by scale factor.
        ! this assumes dy varies in y direction only.
        ! (i.e, dy is the same across same x location)
        if (my_task == master_task) then

            ! initialize with initial dxrect
            work_g1(:,:) = dyrect

            ! check if ny is even or odd
            ! if even, middle 2 rows are center
            ! of odd,  middle 1 row is center
            if (mod(ny_global,2) == 0) then ! ny_global is even

            ! with even number of x locatons,
            ! the center two y columns are center
            center1 = ny_global/2  ! integer math
            center2 = center1 + 1  ! integer math

            else ! ny_global = odd
            ! only one center index. set center2=center1
            center1 = ceiling(real(ny_global/2),int_kind)
            center2 = center1
            endif

            ! note loop over only half the y grid points (center1)-1
            ! working from the center outward.
            do i = 1, nx_global
            do j = 1, center1-1
            ! work from center1 to bottom
            work_g1(i,center1-j) = dyscale*work_g1(i,center1-j+1)

            ! work from center2 to top
            work_g1(i,center2+j) = dyscale*work_g1(i,center2+j-1)
            enddo ! i
            enddo ! j
        endif    ! mytask == master_task
        ! note work_g1 is converted to meters primary_grid_lengths_HTE
        call primary_grid_lengths_HTE(work_g1, save_ghte_ghtn, HTN, &
            dyT, dyU, dyN, dyE, G_HTE, G_HTN) 

        ! make ULAT array
        if (my_task == master_task) then

            ! make first row reference lat in radians.
            ! the remaining work_g1 is still dy in meters
            work_g1(:,1) = latrefrect/rad_to_deg ! radians


            ! loop over remaining points and add spacing to successive
            ! x locations
            do j = 2, ny_global ! start from j=2. j=1 is latrefrect
            do i = 1, nx_global
            length = work_g1(i,j)/radius             ! grid spacing in radians
            work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
            enddo ! i
            enddo ! j
        endif    ! mytask == master_task
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULAT, distrb_info, &
                                ew_boundary_type, ns_boundary_type)


        deallocate(work_g1)

    end subroutine rectgrid_scale_dxdy
  
!=======================================================================

! Complex land mask for testing box cases
! Requires nx_global, ny_global > 20
! Assumes work array has been initialized to 1 (ocean) and north and
! south land boundaries have been applied (ew_boundary_type='cyclic')

    subroutine grid_boxislands_kmt(work)

        real (kind=dbl_kind), dimension(:,:), intent(inout) :: work
  
        integer (kind=int_kind) :: &
           i, j, k, & ! indices
           nxb, nyb   ! convenient cell-block sizes for building the mask
  
        character(len=*), parameter :: subname = '(grid_boxislands_kmt)'
  
        ! number of cells in 5% of global grid x and y lengths
        nxb = int(real(nx_global, dbl_kind) / c20, int_kind)
        nyb = int(real(ny_global, dbl_kind) / c20, int_kind)
  
        if (nxb < 1 .or. nyb < 1) &
           call abort_ice(subname//' ERROR: requires larger grid size', &
                file=__FILE__, line=__LINE__)
  
        ! initialize work area as all ocean (c1).
        work(:,:) = c1
  
        ! now add land points (c0)
        ! northeast triangle
        k = 0
        do j = ny_global, ny_global-3*nyb, -1
           k = k+1
           do i = nx_global-3*nxb+k, nx_global
              work(i,j) = c0
           enddo
        enddo
  
        ! northwest docks
        do j = ny_global-3*nyb, ny_global
           do i = 1, 1
              work(i,j) = c0
           enddo
        enddo
        do i = 1, 2*nxb
           do j = ny_global-3*nyb, ny_global-nyb-2
              work(i,j) = c0
           enddo
           do j = ny_global-nyb, ny_global-nyb+1
              work(i,j) = c0
           enddo
        enddo
  
        ! southwest docks
        do j = 2*nyb, 3*nyb
           do i = 1, 1
              work(i,j) = c0
           enddo
        enddo
        do j = 1, 2*nyb
           do i = 2, nxb
              work(i,j) = c0
           enddo
           do i = 2*nxb-1, 2*nxb
              work(i,j) = c0
           enddo
           do i = 2*nxb+2,4*nxb
              work(i,j) = c0
           enddo
        enddo
  
        ! tiny island
        do j = 14*nyb, 14*nyb+1
           do i = 14*nxb, 14*nxb+1
              work(i,j) = c0
           enddo
        enddo
  
        ! X islands
        ! left triangle
        k = 0
        do i = 2*nxb, 4*nxb
           k=k+1
           do j = 10*nyb+k, 14*nyb-k
              work(i,j) = c0
           enddo
        enddo
        ! upper triangle
        k = 0
        do j = 14*nyb, 12*nyb, -1
           k=k+1
           do i = 2*nxb+2+k, 6*nxb-2-k
              work(i,j) = c0
           enddo
        enddo
        ! diagonal
        k = 0
        do j = 10*nyb, 14*nyb
           k=k+1
           do i = 2*nxb+4+k, 2*nxb+6+k
              work(i,j) = c0
           enddo
        enddo
        ! lower right triangle
        k = 0
        do j = 12*nyb, 10*nyb, -1
           k=k+1
           do i = 5*nxb+k, 8*nxb
              work(i,j) = c0
           enddo
        enddo
  
        ! bar islands
        do i = 10*nxb, 16*nxb
           do j = 4*nyb, 5*nyb
              work(i,j) = c0
           enddo
           do j = 6*nyb+2, 8*nyb
              work(i,j) = c0
           enddo
           do j = 8*nyb+2, 8*nyb+3
              work(i,j) = c0
           enddo
        enddo
  
    end subroutine grid_boxislands_kmt

end module ice_grid_idealised