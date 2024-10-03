module ice_grid_nc

    use ice_kinds_mod
    use ice_broadcast, only: broadcast_scalar
    use ice_boundary, only: ice_HaloExtrapolate
    use ice_constants, only: c0, c1, field_loc_center, field_loc_NEcorner, &
        field_type_angle, field_type_scalar
    use ice_communicate, only: my_task, master_task
    use ice_blocks, only: block, get_block, nx_block, ny_block
    use ice_domain_size, only: nx_global, ny_global, max_blocks
    use ice_domain, only: blocks_ice, nblocks, distrb_info, &
        ew_boundary_type, ns_boundary_type
    use ice_fileunits, only: nu_diag
    use ice_gather_scatter, only: gather_global, scatter_global
    use ice_exit, only: abort_ice
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_warnings_aborted
    use ice_grid_lengths, only: primary_grid_lengths_HTN, primary_grid_lengths_HTE
    use ice_gridbox, only: gridbox_verts
    use ice_read_write, only: ice_read_nc, &
          ice_read_global_nc, ice_open_nc, ice_close_nc, ice_check_nc

    use netcdf, only: nf90_noerr, nf90_inq_varid

    implicit none
    private
    public :: popgrid_nc
#ifdef CESMCOUPLED
    public :: latlongrid
#endif

    contains

!=======================================================================

! POP displaced pole grid and land mask.
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
! Revised for netcdf input: Ann Keen, Met Office, May 2007

    subroutine popgrid_nc(save_ghte_ghtn, grid_file, kmt_file, &
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
            ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
            fid_grid, &            ! file id for netCDF grid file
            fid_kmt                ! file id for netCDF kmt file

        logical (kind=log_kind) :: diag

        character (char_len) :: &
            fieldname              ! field name in netCDF file

        real (kind=dbl_kind) :: &
            pi

        real (kind=dbl_kind), dimension(:,:), allocatable :: &
            work_g1

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
            work1

        type (block) :: &
            this_block           ! block information for current block

        integer(kind=int_kind) :: &
            varid
        integer (kind=int_kind) :: &
            status                ! status flag

        logical (kind=log_kind) :: &
            l_readCenter ! If anglet exist in grid file read it otherwise calculate it

        character(len=*), parameter :: subname = '(popgrid_nc)'
        
        call icepack_query_parameters(pi_out=pi)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

        call ice_open_nc(grid_file,fid_grid)
        call ice_open_nc(kmt_file,fid_kmt)

        diag = .true.       ! write diagnostic info
        !-----------------------------------------------------------------
        ! topography
        !-----------------------------------------------------------------

        fieldname='kmt'
        call ice_read_nc(fid_kmt,1,fieldname,work1,diag, &
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
            if (kmt(i,j,iblk) >= c1) hm(i,j,iblk) = c1
            enddo
            enddo
        enddo
        !$OMP END PARALLEL DO

        !-----------------------------------------------------------------
        ! lat, lon, angle
        !-----------------------------------------------------------------

        allocate(work_g1(nx_global,ny_global))

        fieldname='ulat'
        call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULAT
        call gridbox_verts(work_g1,latt_bounds)
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULAT, distrb_info, &
                                ew_boundary_type, ns_boundary_type)

        fieldname='ulon'
        call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULON
        call gridbox_verts(work_g1,lont_bounds)
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULON, distrb_info, &
                                ew_boundary_type, ns_boundary_type)

        fieldname='angle'
        call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ANGLE
        call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_angle)
        ! fix ANGLE: roundoff error due to single precision
        where (ANGLE >  pi) ANGLE =  pi
        where (ANGLE < -pi) ANGLE = -pi

        ! if grid file includes anglet then read instead
        fieldname='anglet'
        if (my_task == master_task) then
            status = nf90_inq_varid(fid_grid, trim(fieldname) , varid)
            if (status /= nf90_noerr) then
            write(nu_diag,*) subname//' CICE will calculate angleT, TLON and TLAT'
            else
            write(nu_diag,*) subname//' angleT, TLON and TLAT is read from grid file'
            l_readCenter = .true.
            endif
        endif
        call broadcast_scalar(l_readCenter,master_task)
        if (l_readCenter) then
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
            call scatter_global(ANGLET, work_g1, master_task, distrb_info, &
                                field_loc_center, field_type_angle)
            where (ANGLET >  pi) ANGLET =  pi
            where (ANGLET < -pi) ANGLET = -pi
            fieldname="tlon"
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
            call scatter_global(TLON, work_g1, master_task, distrb_info, &
                                field_loc_center, field_type_scalar)
            fieldname="tlat"
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
            call scatter_global(TLAT, work_g1, master_task, distrb_info, &
                                field_loc_center, field_type_scalar)
        endif
        !-----------------------------------------------------------------
        ! cell dimensions
        ! calculate derived quantities from global arrays to preserve
        ! information on boundaries
        !-----------------------------------------------------------------

        fieldname='htn'
        call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTN
        call primary_grid_lengths_HTN(work_g1, save_ghte_ghtn, HTN, &
            dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN)                  
        fieldname='hte'
        call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTE
        call primary_grid_lengths_HTE(work_g1, save_ghte_ghtn, HTE, &
            dyT, dyU, dyN, dyE, G_HTE, G_HTN)                 

        deallocate(work_g1)

        if (my_task == master_task) then
            call ice_close_nc(fid_grid)
            call ice_close_nc(fid_kmt)
        endif
        
    end subroutine popgrid_nc

#ifdef CESMCOUPLED
!=======================================================================

! Read in kmt file that matches CAM lat-lon grid and has single column
! functionality
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90 and cice ncdf calls

    subroutine latlongrid(kmt_file, dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE, HTE, HTN, &
        tarea, uarea, tarear, uarear, &
        ULON, ULAT, TLON, TLAT, ANGLE, ANGLET, ocn_gridcell_frac, hm, kmt)

        use ice_scam, only: scmlat, scmlon, single_column
        use netcdf
        use ice_constants, only: radius

        character (len=char_len_long), intent(in) :: &
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
         tarea  , & ! area of T-cell (m^2), valid in halo
         uarea  , & ! area of U-cell (m^2), valid in halo
         tarear , & ! 1/tarea, valid in halo
         uarear , & ! 1/uarea, valid in halo
         ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
         ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
         TLON   , & ! longitude of temp (T) pts (radians)
         TLAT   , & ! latitude of temp (T) pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET , &  ! ANGLE converted to T-cells, valid in halo
         ocn_gridcell_frac
         
        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)

        
        !local vars 

      integer (kind=int_kind) :: &
         i, j, iblk

      integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           closelat, &        ! Single-column latitude value
           closelon, &        ! Single-column longitude value
           closelatidx, &     ! Single-column latitude index to retrieve
           closelonidx        ! Single-column longitude index to retrieve

      integer (kind=int_kind) :: &
           start(2), &        ! Start index to read in
           count(2)           ! Number of points to read in

      integer (kind=int_kind) :: &
           start3(3), &        ! Start index to read in
           count3(3)           ! Number of points to read in

      integer (kind=int_kind) :: &
        status                ! status flag

      real (kind=dbl_kind), allocatable :: &
           lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries

      real (kind=dbl_kind) :: &
         pos_scmlon,&         ! temporary
         pi, &
         puny, &
         scamdata             ! temporary

      character(len=*), parameter :: subname = '(lonlatgrid)'

      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      call icepack_query_parameters(pi_out=pi, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! Determine dimension of domain file and check for consistency

      if (my_task == master_task) then
         call ice_open_nc(kmt_file, ncid)

         status = nf90_inq_dimid (ncid, 'ni', dimid)
         call ice_check_nc(status, subname//' ERROR: inq_dimid ni', file=__FILE__, line=__LINE__)
         status = nf90_inquire_dimension(ncid, dimid, len=ni)
         call ice_check_nc(status, subname//' ERROR: inq dim ni', file=__FILE__, line=__LINE__)
         status = nf90_inq_dimid (ncid, 'nj', dimid)
         call ice_check_nc(status, subname//' ERROR: inq_dimid nj', file=__FILE__, line=__LINE__)
         status = nf90_inquire_dimension(ncid, dimid, len=nj)
         call ice_check_nc(status, subname//' ERROR: inq dim nj', file=__FILE__, line=__LINE__)
      end if

      ! Determine start/count to read in for either single column or global lat-lon grid
      ! If single_column, then assume that only master_task is used since there is only one task

      if (single_column) then
         ! Check for consistency
         if (my_task == master_task) then
            if ((nx_global /= 1).or. (ny_global /= 1)) then
               write(nu_diag,*) 'Because you have selected the column model flag'
               write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
               write(nu_diag,*) 'ice_domain_size.F and recompile'
               call abort_ice (subname//' ERROR: check nx_global, ny_global', file=__FILE__, line=__LINE__)
            endif
         end if

         ! Read in domain file for single column
         allocate(lats(nj))
         allocate(lons(ni))
         allocate(pos_lons(ni))
         allocate(glob_grid(ni,nj))

         start3=(/1,1,1/)
         count3=(/ni,nj,1/)
         status = nf90_inq_varid(ncid, 'xc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid xc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         call ice_check_nc(status, subname//' ERROR: get_var xc', file=__FILE__, line=__LINE__)
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do

         status = nf90_inq_varid(ncid, 'yc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid yc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         call ice_check_nc(status, subname//' ERROR: get_var yc', file=__FILE__, line=__LINE__)
         do j = 1,nj
            lats(j) = glob_grid(1,j)
         end do

         ! convert lons array and scmlon to 0,360 and find index of value closest to 0
         ! and obtain single-column longitude/latitude indices to retrieve

         pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
         pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)
         start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
         start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

         deallocate(lats)
         deallocate(lons)
         deallocate(pos_lons)
         deallocate(glob_grid)

         status = nf90_inq_varid(ncid, 'xc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid xc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var xc', file=__FILE__, line=__LINE__)
         TLON = scamdata
         status = nf90_inq_varid(ncid, 'yc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid yc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var yc', file=__FILE__, line=__LINE__)
         TLAT = scamdata
         status = nf90_inq_varid(ncid, 'area' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid area', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var are', file=__FILE__, line=__LINE__)
         tarea = scamdata
         status = nf90_inq_varid(ncid, 'mask' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid mask', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var mask', file=__FILE__, line=__LINE__)
         hm = scamdata
         status = nf90_inq_varid(ncid, 'frac' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid frac', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var frac', file=__FILE__, line=__LINE__)
         ocn_gridcell_frac = scamdata
      else
         ! Check for consistency
         if (my_task == master_task) then
            if (nx_global /= ni .and. ny_global /= nj) then
              write(nu_diag,*) 'latlongrid: ni,nj = ',ni,nj
              write(nu_diag,*) 'latlongrid: nx_g,ny_g = ',nx_global, ny_global
              call abort_ice (subname//' ERROR: ni,nj not equal to nx_global,ny_global', &
                              file=__FILE__, line=__LINE__)
            end if
         end if

         ! Read in domain file for global lat-lon grid
         call ice_read_nc(ncid, 1, 'xc'  , TLON             , diag=.true.)
         call ice_read_nc(ncid, 1, 'yc'  , TLAT             , diag=.true.)
         call ice_read_nc(ncid, 1, 'area', tarea            , diag=.true., &
            field_loc=field_loc_center,field_type=field_type_scalar)
         call ice_read_nc(ncid, 1, 'mask', hm               , diag=.true.)
         call ice_read_nc(ncid, 1, 'frac', ocn_gridcell_frac, diag=.true.)
      end if

      if (my_task == master_task) then
         call ice_close_nc(ncid)
      end if

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            ! Convert from degrees to radians
            TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

            ! Convert from degrees to radians
            TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

            ! Convert from radians^2 to m^2
            ! (area in domain file is in radians^2 and tarea is in m^2)
            tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
         end do
         end do
      end do
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      ! The U grid (velocity) is not used when run with sequential CAM
      ! because we only use thermodynamic sea ice.  However, ULAT is used
      ! in the default initialization of CICE so we calculate it here as
      ! a "dummy" so that CICE will initialize with ice.  If a no ice
      ! initialization is OK (or desired) this can be commented out and
      ! ULAT will remain 0 as specified above.  ULAT is located at the
      ! NE corner of the grid cell, TLAT at the center, so here ULAT is
      ! hacked by adding half the latitudinal spacing (in radians) to
      ! TLAT.
      !-----------------------------------------------------------------

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            if (ny_global == 1) then
               uarea(i,j,iblk)  = tarea(i,j,  iblk)
            else
               uarea(i,j,iblk)  = p25*  &
                                 (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                                + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
            endif
            tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
            uarear(i,j,iblk)   = c1/uarea(i,j,iblk)

            if (single_column) then
               ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)
            else
               if (ny_global == 1) then
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)
               else
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)
               endif
            endif
            ULON  (i,j,iblk) = c0
            NLON  (i,j,iblk) = c0
            NLAT  (i,j,iblk) = c0
            ELON  (i,j,iblk) = c0
            ELAT  (i,j,iblk) = c0
            ANGLE (i,j,iblk) = c0

            ANGLET(i,j,iblk) = c0
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxT   (i,j,iblk) = 1.e36_dbl_kind
            dyT   (i,j,iblk) = 1.e36_dbl_kind
            dxU   (i,j,iblk) = 1.e36_dbl_kind
            dyU   (i,j,iblk) = 1.e36_dbl_kind
            dxN   (i,j,iblk) = 1.e36_dbl_kind
            dyN   (i,j,iblk) = 1.e36_dbl_kind
            dxE   (i,j,iblk) = 1.e36_dbl_kind
            dyE   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine latlongrid
#endif

end module ice_grid_nc