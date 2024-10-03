#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================

! Spatial grids, masks, and boundary conditions
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: Option to read from netcdf files (A. Keen, Met Office)
!       Grid reading routines reworked by E. Hunke for boundary values
! 2021: Add N (center of north face) and E (center of east face) grids
!       to support C and CD solvers.  Defining T at center of cells, U at
!       NE corner, N at center of top face, E at center of right face.
!       All cells are quadrilaterals with NE, E, and N associated with
!       directions relative to logical grid.

   module ice_grid

      use ice_kinds_mod
      use ice_broadcast, only: broadcast_array
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c0, c1, c1p5, c2, p5, p25, &
          field_loc_center, field_loc_NEcorner, field_loc_Nface, field_loc_Eface, &
          field_type_scalar, field_type_angle
      use ice_communicate, only: my_task, master_task
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_domain_size, only: nx_global, ny_global, max_blocks
      use ice_domain, only: blocks_ice, nblocks, halo_info, distrb_info, &
          ns_boundary_type, init_domain_distribution
      use ice_fileunits, only: nu_diag, nu_grid, nu_kmt, flush_fileunit
      use ice_gather_scatter, only: scatter_global
      use ice_read_write, only: ice_read_global, &
          ice_read_global_nc, ice_open, ice_open_nc, ice_close_nc
      use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop
      use ice_exit, only: abort_ice
      use ice_global_reductions, only: global_minval, global_maxval
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted, &
          icepack_query_parameters
      
      !components of ice_grid
      use ice_grid_bathy, only: get_bathymetry, get_bathymetry_popfile
      use ice_gridbox, only: gridbox_edges, gridbox_corners
      use ice_grid_neighbor, only: grid_neighbor_min, grid_neighbor_max
      use ice_grid_average, only:  grid_average_X2YS, grid_average_X2YA, &
         grid_average_X2YF, grid_average_X2Y_2
      use ice_grid_latlon, only: Tlatlon, NElatlon
      use ice_grid_idealised, only: rectgrid, rectgrid_scale_dxdy
      use ice_grid_read_bin, only: cpomgrid, popgrid
#ifdef USE_NETCDF
      use ice_grid_read_nc, only: popgrid_nc
#ifdef CESMCOUPLED
      use ice_grid_read_nc, only: latlongrid
#endif
#endif

      implicit none
      private
      public :: init_grid1, init_grid2, grid_average_X2Y, &
                alloc_grid, dealloc_grid, grid_neighbor_min, grid_neighbor_max

      character (len=char_len_long), public :: &
         grid_format  , & ! file format ('bin'=binary or 'nc'=netcdf)
         gridcpl_file , & !  input file for POP coupling grid info
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         kmt_type     , & !  options are file, default, boxislands
         bathymetry_file, & !  input bathymetry for seabed stress
         bathymetry_format, & ! bathymetry file format (default or pop)
         grid_spacing , & !  default of 30.e3m or set by user in namelist
         grid_ice  , & !  Underlying model grid structure (A, B, C, CD)
         grid_ice_thrm, & !  ocean forcing grid for thermo fields (T, U, N, E)
         grid_ice_dynu, & !  ocean forcing grid for dyn U fields  (T, U, N, E)
         grid_ice_dynv, & !  ocean forcing grid for dyn V fields  (T, U, N, E)
         grid_atm     , & !  atmos forcing grid structure (A, B, C, CD)
         grid_atm_thrm, & !  atmos forcing grid for thermo fields (T, U, N, E)
         grid_atm_dynu, & !  atmos forcing grid for dyn U fields  (T, U, N, E)
         grid_atm_dynv, & !  atmos forcing grid for dyn V fields  (T, U, N, E)
         grid_ocn     , & !  ocean forcing grid structure (A B, C, CD)
         grid_ocn_thrm, & !  ocean forcing grid for thermo fields (T, U, N, E)
         grid_ocn_dynu, & !  ocean forcing grid for dyn U fields  (T, U, N, E)
         grid_ocn_dynv, & !  ocean forcing grid for dyn V fields  (T, U, N, E)
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, regional

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
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
         tarea  , & ! area of T-cell (m^2), valid in halo
         uarea  , & ! area of U-cell (m^2), valid in halo
         narea  , & ! area of N-cell (m^2), valid in halo
         earea  , & ! area of E-cell (m^2), valid in halo
         tarear , & ! 1/tarea, valid in halo
         uarear , & ! 1/uarea, valid in halo
         narear , & ! 1/narea, valid in halo
         earear , & ! 1/earea, valid in halo
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
         ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
         TLON   , & ! longitude of temp (T) pts (radians)
         TLAT   , & ! latitude of temp (T) pts (radians)
         NLON   , & ! longitude of center of north face of T pts (radians)
         NLAT   , & ! latitude of center of north face of T pts (radians)
         ELON   , & ! longitude of center of east face of T pts (radians)
         ELAT   , & ! latitude of center of east face of T pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET , & ! ANGLE converted to T-cells, valid in halo
         bathymetry      , & ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac   ! only relevant for lat-lon grids
                             ! gridcell value of [1 - (land fraction)] (T-cell)

      real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         G_HTE  , & ! length of eastern edge of T-cell (global ext.)
         G_HTN      ! length of northern edge of T-cell (global ext.)

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), public ::  &
         dxrect, & !  user_specified spacing (cm) in x-direction (uniform HTN)
         dyrect    !  user_specified spacing (cm) in y-direction (uniform HTE)

      ! growth factor for variable spaced grid
      real (kind=dbl_kind), public ::  &
         dxscale, & !  scale factor for grid spacing in x direction (e.g., 1.02)
         dyscale    !  scale factor for gird spacing in y direction (e.g., 1.02)

      real (kind=dbl_kind), public :: &
         lonrefrect, & ! lower left lon for rectgrid
         latrefrect    ! lower left lat for rectgrid

      ! Corners of grid boxes for history output
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         lont_bounds, & ! longitude of gridbox corners for T point
         latt_bounds, & ! latitude of gridbox corners for T point
         lonu_bounds, & ! longitude of gridbox corners for U point
         latu_bounds, & ! latitude of gridbox corners for U point
         lonn_bounds, & ! longitude of gridbox corners for N point
         latn_bounds, & ! latitude of gridbox corners for N point
         lone_bounds, & ! longitude of gridbox corners for E point
         late_bounds    ! latitude of gridbox corners for E point

      ! masks
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         bm     , & ! task/block id
         uvm    , & ! land/boundary mask (U-cell)
         npm    , & ! land/boundary mask (N-cell)
         epm    , & ! land/boundary mask (E-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)

      logical (kind=log_kind), public :: &
         use_bathymetry, & ! flag for reading in bathymetry_file
         save_ghte_ghtn, & ! flag for saving global hte and htn during initialization
         scale_dxdy        ! flag to apply scale factor to vary dx/dy in rectgrid

      logical (kind=log_kind), dimension (:,:,:), allocatable, public :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask  (U-cell) (1 if all surrounding T cells are ocean)
         umaskCD, & ! land/boundary mask  (U-cell) (1 if at least two surrounding T cells are ocean)
         nmask  , & ! land/boundary mask, (N-cell)
         emask  , & ! land/boundary mask, (E-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         rndex_global       ! global index for local subdomain (dbl)

      interface grid_average_X2Y
      module procedure grid_average_X2Y_base , &
                        grid_average_X2Y_userwghts, &
                        grid_average_X2Y_NEversion
      end interface

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_grid

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(alloc_grid)'

      allocate( &
         dxT      (nx_block,ny_block,max_blocks), & ! width of T-cell through the middle (m)
         dyT      (nx_block,ny_block,max_blocks), & ! height of T-cell through the middle (m)
         dxU      (nx_block,ny_block,max_blocks), & ! width of U-cell through the middle (m)
         dyU      (nx_block,ny_block,max_blocks), & ! height of U-cell through the middle (m)
         dxN      (nx_block,ny_block,max_blocks), & ! width of N-cell through the middle (m)
         dyN      (nx_block,ny_block,max_blocks), & ! height of N-cell through the middle (m)
         dxE      (nx_block,ny_block,max_blocks), & ! width of E-cell through the middle (m)
         dyE      (nx_block,ny_block,max_blocks), & ! height of E-cell through the middle (m)
         HTE      (nx_block,ny_block,max_blocks), & ! length of eastern edge of T-cell (m)
         HTN      (nx_block,ny_block,max_blocks), & ! length of northern edge of T-cell (m)
         tarea    (nx_block,ny_block,max_blocks), & ! area of T-cell (m^2)
         uarea    (nx_block,ny_block,max_blocks), & ! area of U-cell (m^2)
         narea    (nx_block,ny_block,max_blocks), & ! area of N-cell (m^2)
         earea    (nx_block,ny_block,max_blocks), & ! area of E-cell (m^2)
         tarear   (nx_block,ny_block,max_blocks), & ! 1/tarea
         uarear   (nx_block,ny_block,max_blocks), & ! 1/uarea
         narear   (nx_block,ny_block,max_blocks), & ! 1/narea
         earear   (nx_block,ny_block,max_blocks), & ! 1/earea
         tarean   (nx_block,ny_block,max_blocks), & ! area of NH T-cells
         tareas   (nx_block,ny_block,max_blocks), & ! area of SH T-cells
         ULON     (nx_block,ny_block,max_blocks), & ! longitude of U pts, NE corner (radians)
         ULAT     (nx_block,ny_block,max_blocks), & ! latitude of U pts, NE corner (radians)
         TLON     (nx_block,ny_block,max_blocks), & ! longitude of T pts (radians)
         TLAT     (nx_block,ny_block,max_blocks), & ! latitude of T pts (radians)
         NLON     (nx_block,ny_block,max_blocks), & ! longitude of N pts, N face (radians)
         NLAT     (nx_block,ny_block,max_blocks), & ! latitude of N pts, N face (radians)
         ELON     (nx_block,ny_block,max_blocks), & ! longitude of E pts, E face (radians)
         ELAT     (nx_block,ny_block,max_blocks), & ! latitude of E pts, E face (radians)
         ANGLE    (nx_block,ny_block,max_blocks), & ! for conversions between POP grid and lat/lon
         ANGLET   (nx_block,ny_block,max_blocks), & ! ANGLE converted to T-cells
         bathymetry(nx_block,ny_block,max_blocks),& ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac(nx_block,ny_block,max_blocks),& ! only relevant for lat-lon grids
         hm       (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         bm       (nx_block,ny_block,max_blocks), & ! task/block id
         uvm      (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         npm      (nx_block,ny_block,max_blocks), & ! land/boundary mask (N-cell)
         epm      (nx_block,ny_block,max_blocks), & ! land/boundary mask (E-cell)
         kmt      (nx_block,ny_block,max_blocks), & ! ocean topography mask for bathymetry (T-cell)
         tmask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         umask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         umaskCD  (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         nmask    (nx_block,ny_block,max_blocks), & ! land/boundary mask (N-cell)
         emask    (nx_block,ny_block,max_blocks), & ! land/boundary mask (E-cell)
         lmask_n  (nx_block,ny_block,max_blocks), & ! northern hemisphere mask
         lmask_s  (nx_block,ny_block,max_blocks), & ! southern hemisphere mask
         rndex_global(nx_block,ny_block,max_blocks), & ! global index for local subdomain (dbl)
         lont_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for T point
         latt_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for T point
         lonu_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for U point
         latu_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for U point
         lonn_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for N point
         latn_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for N point
         lone_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for E point
         late_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for E point
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory1', file=__FILE__, line=__LINE__)

      if (save_ghte_ghtn) then
         if (my_task == master_task) then
            allocate( &
               G_HTE(nx_global+2*nghost, ny_global+2*nghost), & ! length of eastern edge of T-cell (global ext.)
               G_HTN(nx_global+2*nghost, ny_global+2*nghost), & ! length of northern edge of T-cell (global ext.)
               stat=ierr)
         else
            allocate( &
               G_HTE(1,1), & ! needed for debug checks
               G_HTN(1,1), & ! never used in code
               stat=ierr)
         endif
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory3', file=__FILE__, line=__LINE__)
      endif

      end subroutine alloc_grid

!=======================================================================

!
! DeAllocate space for variables no longer needed after initialization
!
      subroutine dealloc_grid

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(dealloc_grid)'

      if (save_ghte_ghtn) then
         deallocate(G_HTE, G_HTN, stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error1', file=__FILE__, line=__LINE__)
      endif

      end subroutine dealloc_grid

!=======================================================================

! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays.
!
! authors: William Lipscomb and Phil Jones, LANL

      subroutine init_grid1

      integer (kind=int_kind) :: &
         fid_grid, &     ! file id for netCDF grid file
         fid_kmt         ! file id for netCDF kmt file

      character (char_len) :: &
         fieldname       ! field name in netCDF file

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

      integer (kind=int_kind) :: &
         max_blocks_min, & ! min value of max_blocks across procs
         max_blocks_max    ! max value of max_blocks across procs

      real (kind=dbl_kind) :: &
         rad_to_deg

      character(len=*), parameter :: subname = '(init_grid1)'

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      ! check tripole flags here
      ! can't check in init_data because ns_boundary_type is not yet read
      ! can't check in init_domain_blocks because grid_type is not accessible due to circular logic

      if (grid_type == 'tripole' .and. ns_boundary_type /= 'tripole' .and. &
          ns_boundary_type /= 'tripoleT') then
         call abort_ice(subname//' ERROR: grid_type tripole needs tripole ns_boundary_type', &
                        file=__FILE__, line=__LINE__)
      endif

      if (grid_type == 'tripole' .and. (mod(nx_global,2)/=0)) then
         call abort_ice(subname//' ERROR: grid_type tripole requires even nx_global number', &
                        file=__FILE__, line=__LINE__)
      endif

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional'     ) then

         if (trim(grid_format) == 'nc') then

            call ice_open_nc(grid_file,fid_grid)
            call ice_open_nc(kmt_file,fid_kmt)

            fieldname='ulat'
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
            fieldname='kmt'
            call ice_read_global_nc(fid_kmt,1,fieldname,work_g2,.true.)

            if (my_task == master_task) then
               call ice_close_nc(fid_grid)
               call ice_close_nc(fid_kmt)
            endif

         else

            call ice_open(nu_grid,grid_file,64) ! ULAT
            call ice_open(nu_kmt, kmt_file, 32) ! KMT

            call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
            call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

            if (my_task == master_task) then
               close (nu_grid)
               close (nu_kmt)
            endif

         endif

      else   ! rectangular grid

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      !-----------------------------------------------------------------
      ! distribute blocks among processors
      !-----------------------------------------------------------------

      call init_domain_distribution(work_g2, work_g1, grid_ice)  ! KMT, ULAT

      deallocate(work_g1)
      deallocate(work_g2)

      !-----------------------------------------------------------------
      ! write additional domain information
      !-----------------------------------------------------------------

      max_blocks_min = global_minval(max_blocks, distrb_info)
      max_blocks_max = global_maxval(max_blocks, distrb_info)
      if (my_task == master_task) then
        write(nu_diag,*        ) ''
        write(nu_diag,'(2a)'   ) subname,' Block size:'
        write(nu_diag,'(2a,i8)') subname,'   nx_block        = ',nx_block
        write(nu_diag,'(2a,i8)') subname,'   ny_block        = ',ny_block
        write(nu_diag,'(2a,i8)') subname,'   min(max_blocks) = ',max_blocks_min
        write(nu_diag,'(2a,i8)') subname,'   max(max_blocks) = ',max_blocks_max
      endif

      end subroutine init_grid1

!=======================================================================

! Horizontal grid initialization:
!
!     U{LAT,LONG} = true {latitude,longitude} of U points
!     HT{N,E} = cell widths on {N,E} sides of T cell
!     ANGLE = angle between local x direction and true east
!     hm = land mask (c1 for ocean points, c0 for land points)
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!     T-grid and ghost cell values
!     Various grid quantities needed for dynamics and transport
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_grid2

#if defined (_OPENMP)
      use OMP_LIB
#endif

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw, &
         pi, pi2, puny

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      type (block) :: &
         this_block           ! block information for current block

      logical (kind=log_kind) :: &
         l_readCenter ! If anglet exist in grid file read it otherwise calculate it

#if defined (_OPENMP)
      integer(kind=omp_sched_kind) :: ompsk  ! openmp schedule
      integer(kind=int_kind) :: ompcs        ! openmp schedule count
#endif

      character(len=*), parameter :: subname = '(init_grid2)'

      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      l_readCenter = .false.

      call icepack_query_parameters(pi_out=pi, pi2_out=pi2, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional'      ) then
         if (trim(grid_format) == 'nc') then
#ifdef USE_NETCDF
            call popgrid_nc(save_ghte_ghtn, grid_file, kmt_file, &
               dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
               ULON, ULAT, TLON, TLAT, ANGLE, ANGLET, &
               hm, kmt, lont_bounds, latt_bounds, G_HTE, G_HTN)    
               ! read POP grid lengths from nc file
#else
            call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                  file=__FILE__, line=__LINE__)
#endif
         else
            call popgrid(save_ghte_ghtn, grid_file, kmt_file, &
               dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
               ULON, ULAT, TLON, TLAT, ANGLE, ANGLET, &
               hm, kmt, lont_bounds, latt_bounds, G_HTE, G_HTN)        
               ! read POP grid lengths directly
         endif
#ifdef CESMCOUPLED
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid(kmt_file, dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
            tarea, uarea, tarear, uarear, &
            ULON, ULAT, TLON, TLAT, ANGLE, ANGLET, ocn_gridcell_frac, hm, kmt)        
            ! lat lon grid for sequential CESM (CAM mode)
         call makemask
         return
#endif
      elseif (trim(grid_type) == 'cpom_grid') then
         call cpomgrid(save_ghte_ghtn, grid_file, kmt_file, &
            dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
            ULON, ULAT, ANGLE, hm, kmt, G_HTE, G_HTN)
         ! cpom model orca1 type grid
      elseif (scale_dxdy) then
         ! scale grid spacing from center outward.
         ! this different than original method in it
         ! needs to define grid spacing before lat/lon.
         ! original rectgrid defines latlon first
         call rectgrid_scale_dxdy(save_ghte_ghtn, scale_dxdy, dxrect, dyrect, &
            dxscale, dyscale, lonrefrect, latrefrect, &
            G_HTE, G_HTN, hm, kmt, angle, ULON, ULAT, &
            dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN )
      else
         ! regular rectangular grid
         ! original method with addition to use namelist lat/lon reference
         call rectgrid(save_ghte_ghtn, scale_dxdy, dxrect, dyrect, &
            G_HTE, G_HTN, kmt_type, lonrefrect, latrefrect, hm, kmt, angle, &
            ULON, ULAT, dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN )         
      endif

      !-----------------------------------------------------------------
      ! Diagnose OpenMP thread schedule, force order in output
      !-----------------------------------------------------------------

#if defined (_OPENMP)
       !$OMP PARALLEL DO ORDERED PRIVATE(iblk) SCHEDULE(runtime)
       do iblk = 1, nblocks
          if (my_task == master_task) then
             !$OMP ORDERED
             if (iblk == 1) then
                call omp_get_schedule(ompsk,ompcs)
!               write(nu_diag,*) ''
                write(nu_diag,*) subname,' OpenMP runtime thread schedule:'
                write(nu_diag,*) subname,'  omp schedule = ',ompsk,ompcs
             endif
             write(nu_diag,*) subname,' block, thread = ',iblk,OMP_GET_THREAD_NUM()
             !$OMP END ORDERED
          endif
       enddo
       !$OMP END PARALLEL DO
       call flush_fileunit(nu_diag)
#endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      ! Fill halo data locally where possible to avoid missing
      ! data associated with land block elimination
      ! Note: HTN, HTE, dx*, dy* are all defined from global arrays
      ! at halos.
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = 1,ny_block
         do i = 1,nx_block
            tarea(i,j,iblk) = dxT(i,j,iblk)*dyT(i,j,iblk)
            uarea(i,j,iblk) = dxU(i,j,iblk)*dyU(i,j,iblk)
            narea(i,j,iblk) = dxN(i,j,iblk)*dyN(i,j,iblk)
            earea(i,j,iblk) = dxE(i,j,iblk)*dyE(i,j,iblk)

            if (tarea(i,j,iblk) > c0) then
               tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            else
               tarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (uarea(i,j,iblk) > c0) then
               uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            else
               uarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (narea(i,j,iblk) > c0) then
               narear(i,j,iblk) = c1/narea(i,j,iblk)
            else
               narear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (earea(i,j,iblk) > c0) then
               earear(i,j,iblk) = c1/earea(i,j,iblk)
            else
               earear(i,j,iblk) = c0 ! possible on boundaries
            endif

         enddo
         enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ghost cell updates
      ! On the tripole grid, one must be careful with updates of
      !  quantities that involve a difference of cell lengths.
      ! For example, dyhx and dxhy are cell-centered vector components.
      ! Also note that on the tripole grid, cxp and cxm would swap places,
      !  as would cyp and cym.  These quantities are computed only
      !  in north and east ghost cells (above), not south and west.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)

      ! Update just on the tripole seam to ensure bit-for-bit symmetry across seam
      call ice_HaloUpdate (tarea,              halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (uarea,              halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (narea,              halo_info, &
                           field_loc_Nface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (earea,              halo_info, &
                           field_loc_Eface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (tarear,             halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (uarear,             halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (narear,             halo_info, &
                           field_loc_Nface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (earear,             halo_info, &
                           field_loc_Eface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)

      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Calculate ANGLET to be compatible with POP ocean model
      ! First, ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         write(nu_diag,*) subname,' angle = ',minval(ANGLE),maxval(ANGLE),count(out_of_range)
         call abort_ice (subname//' ANGLE out of expected range', &
             file=__FILE__, line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      if (trim(grid_type) == 'cpom_grid') then
         ANGLET(:,:,:) = ANGLE(:,:,:)
      else if (.not. (l_readCenter)) then
         ANGLET = c0

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
         !$OMP                     angle_0,angle_w,angle_s,angle_sw)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               angle_0  = ANGLE(i  ,j  ,iblk) !   w----0
               angle_w  = ANGLE(i-1,j  ,iblk) !   |    |
               angle_s  = ANGLE(i,  j-1,iblk) !   |    |
               angle_sw = ANGLE(i-1,j-1,iblk) !   sw---s
               ANGLET(i,j,iblk) = atan2(p25*(sin(angle_0)+ &
                                             sin(angle_w)+ &
                                             sin(angle_s)+ &
                                             sin(angle_sw)),&
                                        p25*(cos(angle_0)+ &
                                             cos(angle_w)+ &
                                             cos(angle_s)+ &
                                             cos(angle_sw)))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif ! cpom_grid

      if (trim(grid_type) == 'regional' .and. &
          (.not. (l_readCenter))) then
         ! for W boundary extrapolate from interior
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            i = ilo
            if (this_block%i_glob(i) == 1) then
               do j = jlo, jhi
                  ANGLET(i,j,iblk) = c2*ANGLET(i+1,j,iblk)-ANGLET(i+2,j,iblk)
               enddo
            endif
         enddo
         !$OMP END PARALLEL DO
      endif  ! regional

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (ANGLET,           halo_info, &
                           field_loc_center, field_type_angle, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      call makemask          ! velocity mask, hemisphere masks
      if (.not. (l_readCenter)) then
         call Tlatlon(ULAT, ULON, grid_type, TLON, TLAT)  ! get lat, lon on the T grid
      endif
      call NElatlon(ULON, ULAT, TLON, TLAT, & 
         umask, tmask, emask, nmask, grid_type, NLON, NLAT, ELON, ELAT)
         ! get lat, lon on the N, E grid

      !-----------------------------------------------------------------
      ! bathymetry
      !-----------------------------------------------------------------

      if (trim(bathymetry_format) == 'default') then
         call get_bathymetry(use_bathymetry, bathymetry_file, bathymetry, kmt)
      elseif (trim(bathymetry_format) == 'pop') then
         call get_bathymetry_popfile(use_bathymetry, bathymetry_file, bathymetry, kmt)
      else
         call abort_ice(subname//' ERROR: bathymetry_format value must be default or pop', &
            file=__FILE__, line=__LINE__)
      endif

      !----------------------------------------------------------------
      ! Corner coordinates for CF compliant history files
      !----------------------------------------------------------------

      call gridbox_corners(TLAT, TLON, latu_bounds, lonu_bounds, lont_bounds)
      call gridbox_edges(ELAT, ELON, NLAT, NLON,latn_bounds, lonn_bounds, late_bounds, lone_bounds)

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine init_grid2




!=======================================================================

! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask)
! and N and E cells (nmask, emask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! author: Elizabeth C. Hunke, LANL

      subroutine makemask

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         puny

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
            uvmCD

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(makemask)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (kmt,              halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_HaloUpdate (hm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      bm = c0
      allocate(uvmCD(nx_block,ny_block,max_blocks))

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk), &
                                 hm(i,j+1,iblk), hm(i+1,j+1,iblk))
            npm(i,j,iblk) = min (hm(i,j,  iblk), hm(i,j+1,iblk))
            epm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,iblk))
            bm(i,j,iblk) = my_task + iblk/100.0_dbl_kind
            uvmCD(i,j,iblk) = (hm(i,j,  iblk)+hm(i+1,j,  iblk) &
                            +  hm(i,j+1,iblk)+hm(i+1,j+1,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uvm,                halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_HaloUpdate (uvmCD,              halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_HaloUpdate (npm,                halo_info, &
                           field_loc_Nface,    field_type_scalar)
      call ice_HaloUpdate (epm,                halo_info, &
                           field_loc_Eface,    field_type_scalar)
      call ice_HaloUpdate (bm,                 halo_info, &
                           field_loc_center,   field_type_scalar)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! needs to cover halo (no halo update for logicals)
         tmask(:,:,iblk)   = .false.
         umask(:,:,iblk)   = .false.
         umaskCD(:,:,iblk) = .false.
         nmask(:,:,iblk)   = .false.
         emask(:,:,iblk)   = .false.
         do j = jlo-nghost, jhi+nghost
         do i = ilo-nghost, ihi+nghost
            if ( hm(i,j,iblk)   > p5  ) tmask  (i,j,iblk)   = .true.
            if (uvm(i,j,iblk)   > p5  ) umask  (i,j,iblk)   = .true.
            if (uvmCD(i,j,iblk) > c1p5) umaskCD(i,j,iblk)   = .true.
            if (npm(i,j,iblk)   > p5  ) nmask  (i,j,iblk)   = .true.
            if (epm(i,j,iblk)   > p5  ) emask  (i,j,iblk)   = .true.
         enddo
         enddo

      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:,:,iblk) = .false.
         lmask_s(:,:,iblk) = .false.

         tarean(:,:,iblk) = c0
         tareas(:,:,iblk) = c0

         do j = jlo,jhi
         do i = ilo,ihi

            if (ULAT(i,j,iblk) >= -puny) then
               lmask_n(i,j,iblk) = .true. ! N. Hem.
            else
               lmask_s(i,j,iblk) = .true. ! S. Hem.
            endif

            ! N hemisphere area mask (m^2)
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

         enddo
         enddo

      enddo  ! iblk
      !$OMP END PARALLEL DO

      deallocate(uvmCD)

      end subroutine makemask

!=======================================================================

! Shifts quantities from one grid to another
! Constructs the shift based on the grid
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

    subroutine grid_average_X2Y_base(type,work1,grid1,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_base)'

      if (trim(grid1) == trim(grid2)) then
         work2 = work1
      else
         X2Y = trim(grid1)//'2'//trim(grid2)//trim(type)
         call grid_average_X2Y_1(X2Y,work1,work2)
      endif

      end subroutine grid_average_X2Y_base

!=======================================================================

! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_userwghts(type,work1,grid1,wght1,mask1,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         mask1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_userwghts)'

      if (trim(grid1) == trim(grid2)) then
         work2 = work1
      else
         X2Y = trim(grid1)//'2'//trim(grid2)//trim(type)
         call grid_average_X2Y_1f(X2Y,work1,wght1,mask1,work2)
      endif

      end subroutine grid_average_X2Y_userwghts

!=======================================================================

! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_NEversion(type,work1a,grid1a,work1b,grid1b,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1a, grid1b, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1a(:,:,:), work1b(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_NEversion)'

      X2Y = trim(grid1a)//trim(grid1b)//'2'//trim(grid2)//trim(type)

      select case (trim(X2Y))

         ! state masked
         case('NE2US')
            call grid_average_X2Y_2('NE2US',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2US')
            call grid_average_X2Y_2('NE2US',work1b,narea,npm,work1a,earea,epm,work2)
         case('NE2TS')
            call grid_average_X2Y_2('NE2TS',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2TS')
            call grid_average_X2Y_2('NE2TS',work1b,narea,npm,work1a,earea,epm,work2)

         ! state unmasked
         case('NE2UA')
            call grid_average_X2Y_2('NE2UA',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2UA')
            call grid_average_X2Y_2('NE2UA',work1b,narea,npm,work1a,earea,epm,work2)
         case('NE2TA')
            call grid_average_X2Y_2('NE2TA',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2TA')
            call grid_average_X2Y_2('NE2TA',work1b,narea,npm,work1a,earea,epm,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2Y_NEversion

        !=======================================================================
  
  ! Shifts quantities from one grid to another
  ! NOTE: Input array includes ghost cells that must be updated before
  !       calling this routine.
  !
  ! author: T. Craig
  
   subroutine grid_average_X2Y_1(X2Y,work1,work2)
  
      character(len=*) , intent(in) :: &
         X2Y

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=*), parameter :: subname = '(grid_average_X2Y_1)'

      select case (trim(X2Y))

         ! flux unmasked
         case('T2UF')
            call grid_average_X2YF('NE',work1,tarea,work2,uarea)
         case('T2EF')
            call grid_average_X2YF('E' ,work1,tarea,work2,earea)
         case('T2NF')
            call grid_average_X2YF('N' ,work1,tarea,work2,narea)
         case('U2TF')
            call grid_average_X2YF('SW',work1,uarea,work2,tarea)
         case('U2EF')
            call grid_average_X2YF('S' ,work1,uarea,work2,earea)
         case('U2NF')
            call grid_average_X2YF('W' ,work1,uarea,work2,narea)
         case('E2TF')
            call grid_average_X2YF('W' ,work1,earea,work2,tarea)
         case('E2UF')
            call grid_average_X2YF('N' ,work1,earea,work2,uarea)
         case('E2NF')
            call grid_average_X2YF('NW',work1,earea,work2,narea)
         case('N2TF')
            call grid_average_X2YF('S' ,work1,narea,work2,tarea)
         case('N2UF')
            call grid_average_X2YF('E' ,work1,narea,work2,uarea)
         case('N2EF')
            call grid_average_X2YF('SE',work1,narea,work2,earea)

         ! state masked
         case('T2US')
            call grid_average_X2YS('NE',work1,tarea,hm ,work2)
         case('T2ES')
            call grid_average_X2YS('E' ,work1,tarea,hm ,work2)
         case('T2NS')
            call grid_average_X2YS('N' ,work1,tarea,hm ,work2)
         case('U2TS')
            call grid_average_X2YS('SW',work1,uarea,uvm,work2)
         case('U2ES')
            call grid_average_X2YS('S' ,work1,uarea,uvm,work2)
         case('U2NS')
            call grid_average_X2YS('W' ,work1,uarea,uvm,work2)
         case('E2TS')
            call grid_average_X2YS('W' ,work1,earea,epm,work2)
         case('E2US')
            call grid_average_X2YS('N' ,work1,earea,epm,work2)
         case('E2NS')
            call grid_average_X2YS('NW',work1,earea,epm,work2)
         case('N2TS')
            call grid_average_X2YS('S' ,work1,narea,npm,work2)
         case('N2US')
            call grid_average_X2YS('E' ,work1,narea,npm,work2)
         case('N2ES')
            call grid_average_X2YS('SE',work1,narea,npm,work2)

         ! state unmasked
         case('T2UA')
            call grid_average_X2YA('NE',work1,tarea,work2)
         case('T2EA')
            call grid_average_X2YA('E' ,work1,tarea,work2)
         case('T2NA')
            call grid_average_X2YA('N' ,work1,tarea,work2)
         case('U2TA')
            call grid_average_X2YA('SW',work1,uarea,work2)
         case('U2EA')
            call grid_average_X2YA('S' ,work1,uarea,work2)
         case('U2NA')
            call grid_average_X2YA('W' ,work1,uarea,work2)
         case('E2TA')
            call grid_average_X2YA('W' ,work1,earea,work2)
         case('E2UA')
            call grid_average_X2YA('N' ,work1,earea,work2)
         case('E2NA')
            call grid_average_X2YA('NW',work1,earea,work2)
         case('N2TA')
            call grid_average_X2YA('S' ,work1,narea,work2)
         case('N2UA')
            call grid_average_X2YA('E' ,work1,narea,work2)
         case('N2EA')
            call grid_average_X2YA('SE',work1,narea,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select
   
   end subroutine grid_average_X2Y_1

   !=======================================================================
  
  ! Shifts quantities from one grid to another
  ! NOTE: Input array includes ghost cells that must be updated before
  !       calling this routine.
  !
  ! author: T. Craig
  
   subroutine grid_average_X2Y_1f(X2Y,work1,wght1,mask1,work2)
  
      character(len=*) , intent(in) :: &
         X2Y

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         mask1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=*), parameter :: subname = '(grid_average_X2Y_1f)'

      select case (trim(X2Y))

! don't support these for now, requires extra destination wght
!         ! flux unmasked
!         case('T2UF')
!            call grid_average_X2YF('NE',work1,tarea,work2,uarea)
!         case('T2EF')
!            call grid_average_X2YF('E' ,work1,tarea,work2,earea)
!         case('T2NF')
!            call grid_average_X2YF('N' ,work1,tarea,work2,narea)
!         case('U2TF')
!            call grid_average_X2YF('SW',work1,uarea,work2,tarea)
!         case('U2EF')
!            call grid_average_X2YF('S' ,work1,uarea,work2,earea)
!         case('U2NF')
!            call grid_average_X2YF('W' ,work1,uarea,work2,narea)
!         case('E2TF')
!            call grid_average_X2YF('W' ,work1,earea,work2,tarea)
!         case('E2UF')
!            call grid_average_X2YF('N' ,work1,earea,work2,uarea)
!         case('E2NF')
!            call grid_average_X2YF('NW',work1,earea,work2,narea)
!         case('N2TF')
!            call grid_average_X2YF('S' ,work1,narea,work2,tarea)
!         case('N2UF')
!            call grid_average_X2YF('E' ,work1,narea,work2,uarea)
!         case('N2EF')
!            call grid_average_X2YF('SE',work1,narea,work2,earea)

         ! state masked
         case('T2US')
            call grid_average_X2YS('NE',work1,wght1,mask1,work2)
         case('T2ES')
            call grid_average_X2YS('E' ,work1,wght1,mask1,work2)
         case('T2NS')
            call grid_average_X2YS('N' ,work1,wght1,mask1,work2)
         case('U2TS')
            call grid_average_X2YS('SW',work1,wght1,mask1,work2)
         case('U2ES')
            call grid_average_X2YS('S' ,work1,wght1,mask1,work2)
         case('U2NS')
            call grid_average_X2YS('W' ,work1,wght1,mask1,work2)
         case('E2TS')
            call grid_average_X2YS('W' ,work1,wght1,mask1,work2)
         case('E2US')
            call grid_average_X2YS('N' ,work1,wght1,mask1,work2)
         case('E2NS')
            call grid_average_X2YS('NW',work1,wght1,mask1,work2)
         case('N2TS')
            call grid_average_X2YS('S' ,work1,wght1,mask1,work2)
         case('N2US')
            call grid_average_X2YS('E' ,work1,wght1,mask1,work2)
         case('N2ES')
            call grid_average_X2YS('SE',work1,wght1,mask1,work2)

         ! state unmasked
         case('T2UA')
            call grid_average_X2YA('NE',work1,wght1,work2)
         case('T2EA')
            call grid_average_X2YA('E' ,work1,wght1,work2)
         case('T2NA')
            call grid_average_X2YA('N' ,work1,wght1,work2)
         case('U2TA')
            call grid_average_X2YA('SW',work1,wght1,work2)
         case('U2EA')
            call grid_average_X2YA('S' ,work1,wght1,work2)
         case('U2NA')
            call grid_average_X2YA('W' ,work1,wght1,work2)
         case('E2TA')
            call grid_average_X2YA('W' ,work1,wght1,work2)
         case('E2UA')
            call grid_average_X2YA('N' ,work1,wght1,work2)
         case('E2NA')
            call grid_average_X2YA('NW',work1,wght1,work2)
         case('N2TA')
            call grid_average_X2YA('S' ,work1,wght1,work2)
         case('N2UA')
            call grid_average_X2YA('E' ,work1,wght1,work2)
         case('N2EA')
            call grid_average_X2YA('SE',work1,wght1,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select

  end subroutine grid_average_X2Y_1f

!=======================================================================

   end module ice_grid

!=======================================================================
