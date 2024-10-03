module ice_grid_read_bin

    use ice_kinds_mod
    use ice_broadcast, only: broadcast_scalar
    use ice_boundary, only: ice_HaloExtrapolate
    use ice_constants, only: c0, c1, p5, field_loc_center, field_loc_NEcorner, &
        field_type_angle, field_type_scalar, m_to_cm
    use ice_communicate, only: my_task, master_task
    use ice_blocks, only: block, get_block, nx_block, ny_block
    use ice_domain_size, only: nx_global, ny_global, max_blocks
    use ice_domain, only: blocks_ice, nblocks, distrb_info, &
        ew_boundary_type, ns_boundary_type
    use ice_fileunits, only: nu_diag, nu_grid, nu_kmt
    use ice_gather_scatter, only: gather_global, scatter_global
    use ice_exit, only: abort_ice
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_warnings_aborted
    use ice_grid_lengths, only: primary_grid_lengths_HTN, primary_grid_lengths_HTE
    use ice_gridbox, only: gridbox_verts
    use ice_read_write, only: ice_read, ice_read_global, ice_open

    implicit none
    private
    public :: popgrid, cpomgrid

    contains

!=======================================================================

! POP displaced pole grid and land mask (or tripole).
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)
!
! Land mask record number and field is (1) KMT.
!
! author: Elizabeth C. Hunke, LANL

    subroutine popgrid(save_ghte_ghtn, grid_file, kmt_file, &
        dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
        ULON, ULAT, TLON, TLAT, ANGLE, ANGLET, &
        hm, kmt, lont_bounds, latt_bounds, G_HTE, G_HTN)

        logical (kind=log_kind), intent(in) :: &
         save_ghte_ghtn      ! flag for saving global hte and htn during initialization

        character (len=char_len_long), intent(in) :: &
         grid_file    , & !  input file for POP grid info
         kmt_file      !  input file for POP grid info

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
         HTN    , & ! length of northern edge of T-cell (m)
         ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
         ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
         TLON   , & ! longitude of temp (T) pts (radians)
         TLAT   , & ! latitude of temp (T) pts (radians)
         ANGLE    , & ! for conversions between POP grid and lat/lon
         ANGLET  ! ANGLE converted to T-cells, valid in halo
         
        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)
        
        ! Corners of grid boxes for history output
        real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
            lont_bounds, & ! longitude of gridbox corners for T point
            latt_bounds    ! latitude of gridbox corners for T point

        real (kind=dbl_kind), dimension (:,:), intent(out)  :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)

        !local vars 

        integer (kind=int_kind) :: &
           i, j, iblk, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        logical (kind=log_kind) :: diag
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g1
  
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
           work1
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(popgrid)'
  
        call ice_open(nu_grid,grid_file,64)
        call ice_open(nu_kmt,kmt_file,32)
  
        diag = .true.       ! write diagnostic info
  
        !-----------------------------------------------------------------
        ! topography
        !-----------------------------------------------------------------
  
        call ice_read(nu_kmt,1,work1,'ida4',diag, &
                      field_loc=field_loc_center, &
                      field_type=field_type_scalar)
  
        hm (:,:,:) = c0
        kmt(:,:,:) = c0
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
              kmt(i,j,iblk) = work1(i,j,iblk)
              if (kmt(i,j,iblk) >= p5) hm(i,j,iblk) = c1
           enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
  
        !-----------------------------------------------------------------
        ! lat, lon, angle
        !-----------------------------------------------------------------
  
        allocate(work_g1(nx_global,ny_global))
  
        call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)   ! ULAT
        call gridbox_verts(work_g1,latt_bounds)
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULAT, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
  
        call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULON
        call gridbox_verts(work_g1,lont_bounds)
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULON, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
  
        call ice_read_global(nu_grid,7,work_g1,'rda8',.true.)   ! ANGLE
        call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_angle)
  
        !-----------------------------------------------------------------
        ! cell dimensions
        ! calculate derived quantities from global arrays to preserve
        ! information on boundaries
        !-----------------------------------------------------------------
  
        call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! HTN
        call primary_grid_lengths_HTN(work_g1, save_ghte_ghtn, HTN, &
            dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN) 
  
        call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTE
        call primary_grid_lengths_HTE(work_g1, save_ghte_ghtn, HTN, &
           dyT, dyU, dyN, dyE, G_HTE, G_HTN)                 
  
        deallocate(work_g1)
  
        if (my_task == master_task) then
           close (nu_grid)
           close (nu_kmt)
        endif
  
        end subroutine popgrid

!=======================================================================

! CPOM displaced pole grid and land mask. \\
! Grid record number, field and units are: \\
! (1) ULAT  (degrees)    \\
! (2) ULON  (degrees)    \\
! (3) HTN   (m)          \\
! (4) HTE   (m)          \\
! (7) ANGLE (radians)    \\
!
! Land mask record number and field is (1) KMT.
!
! author: Adrian K. Turner, CPOM, UCL, 09/08/06

    subroutine cpomgrid(save_ghte_ghtn, grid_file, kmt_file, &
        dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
        ULON, ULAT, ANGLE, hm, kmt, G_HTE, G_HTN)

        logical (kind=log_kind), intent(in) :: &
         save_ghte_ghtn      ! flag for saving global hte and htn during initialization

        character (len=char_len_long), intent(in) :: &
         grid_file    , & !  input file for POP grid info
         kmt_file      !  input file for POP grid info

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
         HTN    , & ! length of northern edge of T-cell (m)
         ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
         ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
         ANGLE    ! for conversions between POP grid and lat/lon
         
        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)

        real (kind=dbl_kind), dimension (:,:), intent(out)  :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)

        !local vars 

        integer (kind=int_kind) :: &
             i, j, iblk,           &
             ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        logical (kind=log_kind) :: diag
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g1
  
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
           work1
  
        real (kind=dbl_kind) :: &
           rad_to_deg
  
        type (block) :: &
             this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(cpomgrid)'
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        call ice_open(nu_grid,grid_file,64)
        call ice_open(nu_kmt,kmt_file,32)
  
        diag = .true.       ! write diagnostic info
  
        ! topography
        call ice_read(nu_kmt,1,work1,'ida4',diag)
  
        hm (:,:,:) = c0
        kmt(:,:,:) = c0
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
              kmt(i,j,iblk) = work1(i,j,iblk)
              if (kmt(i,j,iblk) >= c1) hm(i,j,iblk) = c1
           enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
  
        allocate(work_g1(nx_global,ny_global))
  
        ! lat, lon, cell dimensions, angles
        call ice_read_global(nu_grid,1,work_g1, 'rda8',diag)
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
  
        call ice_read_global(nu_grid,2,work_g1, 'rda8',diag)
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
  
        call ice_read_global(nu_grid,3,work_g1,  'rda8',diag)
        work_g1 = work_g1 * m_to_cm
        call primary_grid_lengths_HTN(work_g1, save_ghte_ghtn, HTN, &
           dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN)  
  
        call ice_read_global(nu_grid,4,work_g1,  'rda8',diag)
        work_g1 = work_g1 * m_to_cm
        call primary_grid_lengths_HTE(work_g1, save_ghte_ghtn, HTE, &
           dyT, dyU, dyN, dyE, G_HTE, G_HTN) 
  
        call ice_read_global(nu_grid,7,work_g1,'rda8',diag)
        call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
  
        ! fix units
        ULAT  = ULAT  / rad_to_deg
        ULON  = ULON  / rad_to_deg
  
        deallocate(work_g1)
  
        if (my_task == master_task) then
           close (nu_grid)
           close (nu_kmt)
        endif
  
        write(nu_diag,*) subname," min/max HTN: ", minval(HTN), maxval(HTN)
        write(nu_diag,*) subname," min/max HTE: ", minval(HTE), maxval(HTE)
  
    end subroutine cpomgrid

end module ice_grid_read_bin