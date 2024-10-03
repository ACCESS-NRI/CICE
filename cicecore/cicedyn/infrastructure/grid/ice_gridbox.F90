module ice_gridbox

    use ice_kinds_mod
    use ice_constants, only: c0, c2, c360, field_loc_NEcorner, field_type_scalar
    use ice_communicate, only: my_task, master_task
    use ice_fileunits, only: nu_diag
    use ice_exit, only: abort_ice
    use ice_domain, only: blocks_ice, nblocks, distrb_info
    use ice_blocks, only: block, get_block, nx_block, ny_block
    use ice_domain_size, only: nx_global, ny_global, max_blocks
    use icepack_intfc, only: icepack_warnings_flush, icepack_query_parameters, &
        icepack_warnings_aborted
    use ice_gather_scatter, only: scatter_global, gather_global
    
    implicit none
    private
    public :: gridbox_verts, gridbox_edges, gridbox_corners

    contains

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================

! These fields are only used for netcdf history output, and the
! ghost cell values are not needed.
! NOTE:  Extrapolations were used: these fields are approximate!
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL

    subroutine gridbox_corners(TLAT, TLON, latu_bounds, lonu_bounds, lont_bounds)

        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
            TLAT, TLON

        real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
        lonu_bounds, & ! longitude of gridbox corners for U point
        latu_bounds    ! latitude of gridbox corners for U point

        real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: lont_bounds

        !local vars

        integer (kind=int_kind) :: &
            i,j,iblk,icorner,& ! index counters
            ilo,ihi,jlo,jhi    ! beginning and end of physical domain
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g2
  
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
           work1
  
        real (kind=dbl_kind) :: &
           rad_to_deg
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(gridbox_corners)'
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        !-------------------------------------------------------------
        ! Get coordinates of grid boxes for each block as follows:
        ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
        !-------------------------------------------------------------
  
        latu_bounds(:,:,:,:) = c0
        lonu_bounds(:,:,:,:) = c0
  
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
  
              latu_bounds(1,i,j,iblk)=TLAT(i  ,j  ,iblk)*rad_to_deg
              latu_bounds(2,i,j,iblk)=TLAT(i+1,j  ,iblk)*rad_to_deg
              latu_bounds(3,i,j,iblk)=TLAT(i+1,j+1,iblk)*rad_to_deg
              latu_bounds(4,i,j,iblk)=TLAT(i  ,j+1,iblk)*rad_to_deg
  
              lonu_bounds(1,i,j,iblk)=TLON(i  ,j  ,iblk)*rad_to_deg
              lonu_bounds(2,i,j,iblk)=TLON(i+1,j  ,iblk)*rad_to_deg
              lonu_bounds(3,i,j,iblk)=TLON(i+1,j+1,iblk)*rad_to_deg
              lonu_bounds(4,i,j,iblk)=TLON(i  ,j+1,iblk)*rad_to_deg
  
           enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
  
        !----------------------------------------------------------------
        ! extrapolate on global grid to get edge values
        !----------------------------------------------------------------
  
        if (my_task == master_task) then
           allocate(work_g2(nx_global,ny_global))
        else
           allocate(work_g2(1,1))
        endif
  
        work1(:,:,:) = latu_bounds(2,:,:,:)
  !     work_g2 = c0
  
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latu_bounds(2,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = latu_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latu_bounds(3,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = latu_bounds(4,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latu_bounds(4,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lonu_bounds(2,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonu_bounds(2,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lonu_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonu_bounds(3,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lonu_bounds(4,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonu_bounds(4,:,:,:) = work1(:,:,:)
  
        deallocate(work_g2)
  
        !----------------------------------------------------------------
        ! Convert longitude to Degrees East >0 for history output
        !----------------------------------------------------------------
  
        allocate(work_g2(nx_block,ny_block))  ! not used as global here
        !OMP fails in this loop
        do iblk = 1, nblocks
           do icorner = 1, 4
              work_g2(:,:) = lont_bounds(icorner,:,:,iblk) + c360
              where (work_g2 > c360) work_g2 = work_g2 - c360
              where (work_g2 < c0 )  work_g2 = work_g2 + c360
              lont_bounds(icorner,:,:,iblk) = work_g2(:,:)
              work_g2(:,:) = lonu_bounds(icorner,:,:,iblk) + c360
              where (work_g2 > c360) work_g2 = work_g2 - c360
              where (work_g2 < c0 )  work_g2 = work_g2 + c360
              lonu_bounds(icorner,:,:,iblk) = work_g2(:,:)
           enddo
        enddo
        deallocate(work_g2)
  
        end subroutine gridbox_corners

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================

! These fields are only used for netcdf history output, and the
! ghost cell values are not needed.
! NOTE:  Extrapolations were used: these fields are approximate!
!

    subroutine gridbox_edges( ELAT, ELON, NLAT, NLON, latn_bounds, lonn_bounds, late_bounds, lone_bounds)

        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
            ELAT, ELON, NLAT, NLON

        ! Corners of grid boxes for history output
        real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
        lonn_bounds, & ! longitude of gridbox corners for N point
        latn_bounds, & ! latitude of gridbox corners for N point
        lone_bounds, & ! longitude of gridbox corners for E point
        late_bounds    ! latitude of gridbox corners for E point

        !local vars

        integer (kind=int_kind) :: &
            i,j,iblk,icorner,& ! index counters
            ilo,ihi,jlo,jhi    ! beginning and end of physical domain
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g2
  
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
           work1
  
        real (kind=dbl_kind) :: &
           rad_to_deg
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(gridbox_edges)'
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        !-------------------------------------------------------------
        ! Get coordinates of grid boxes for each block as follows:
        ! for N pt: (1) W edge, (2) E edge, (3) E edge j+1, (4) W edge j+1
        ! for E pt: (1) S edge, (2) S edge i+1, (3) N edge, i+1 (4) N edge
        !-------------------------------------------------------------
  
        latn_bounds(:,:,:,:) = c0
        lonn_bounds(:,:,:,:) = c0
        late_bounds(:,:,:,:) = c0
        lone_bounds(:,:,:,:) = c0
  
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
  
              latn_bounds(1,i,j,iblk)=ELAT(i-1,j  ,iblk)*rad_to_deg
              latn_bounds(2,i,j,iblk)=ELAT(i  ,j  ,iblk)*rad_to_deg
              latn_bounds(3,i,j,iblk)=ELAT(i  ,j+1,iblk)*rad_to_deg
              latn_bounds(4,i,j,iblk)=ELAT(i-1,j+1,iblk)*rad_to_deg
  
              lonn_bounds(1,i,j,iblk)=ELON(i-1,j  ,iblk)*rad_to_deg
              lonn_bounds(2,i,j,iblk)=ELON(i  ,j  ,iblk)*rad_to_deg
              lonn_bounds(3,i,j,iblk)=ELON(i  ,j+1,iblk)*rad_to_deg
              lonn_bounds(4,i,j,iblk)=ELON(i-1,j+1,iblk)*rad_to_deg
  
              late_bounds(1,i,j,iblk)=NLAT(i  ,j-1,iblk)*rad_to_deg
              late_bounds(2,i,j,iblk)=NLAT(i+1,j-1,iblk)*rad_to_deg
              late_bounds(3,i,j,iblk)=NLAT(i+1,j  ,iblk)*rad_to_deg
              late_bounds(4,i,j,iblk)=NLAT(i  ,j  ,iblk)*rad_to_deg
  
              lone_bounds(1,i,j,iblk)=NLON(i  ,j-1,iblk)*rad_to_deg
              lone_bounds(2,i,j,iblk)=NLON(i+1,j-1,iblk)*rad_to_deg
              lone_bounds(3,i,j,iblk)=NLON(i+1,j  ,iblk)*rad_to_deg
              lone_bounds(4,i,j,iblk)=NLON(i  ,j  ,iblk)*rad_to_deg
  
           enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
  
        !----------------------------------------------------------------
        ! extrapolate on global grid to get edge values
        !----------------------------------------------------------------
  
        if (my_task == master_task) then
           allocate(work_g2(nx_global,ny_global))
        else
           allocate(work_g2(1,1))
        endif
  
        ! latn_bounds
  
        work1(:,:,:) = latn_bounds(1,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) &
                              - work_g2(3,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latn_bounds(1,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = latn_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latn_bounds(3,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = latn_bounds(4,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) &
                              - work_g2(3,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        latn_bounds(4,:,:,:) = work1(:,:,:)
  
        ! lonn_bounds
  
        work1(:,:,:) = lonn_bounds(1,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) &
                              - work_g2(3,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonn_bounds(1,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lonn_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonn_bounds(3,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lonn_bounds(4,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                      - work_g2(i,ny_global-2)
           enddo
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) &
                              - work_g2(3,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lonn_bounds(4,:,:,:) = work1(:,:,:)
  
        ! late_bounds
  
        work1(:,:,:) = late_bounds(1,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g2(i,2) &
                              - work_g2(i,3)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        late_bounds(1,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = late_bounds(2,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g2(i,2) &
                              - work_g2(i,3)
           enddo
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        late_bounds(2,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = late_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        late_bounds(3,:,:,:) = work1(:,:,:)
  
        ! lone_bounds
  
        work1(:,:,:) = lone_bounds(1,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g2(i,2) &
                              - work_g2(i,3)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lone_bounds(1,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lone_bounds(2,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g2(i,2) &
                              - work_g2(i,3)
           enddo
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lone_bounds(2,:,:,:) = work1(:,:,:)
  
        work1(:,:,:) = lone_bounds(3,:,:,:)
        call gather_global(work_g2, work1, master_task, distrb_info)
        if (my_task == master_task) then
           do j = 1, ny_global
              work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                      - work_g2(nx_global-2,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        lone_bounds(3,:,:,:) = work1(:,:,:)
  
        deallocate(work_g2)
  
        !----------------------------------------------------------------
        ! Convert longitude to Degrees East >0 for history output
        !----------------------------------------------------------------
  
        allocate(work_g2(nx_block,ny_block))  ! not used as global here
        !OMP fails in this loop
        do iblk = 1, nblocks
           do icorner = 1, 4
              work_g2(:,:) = lonn_bounds(icorner,:,:,iblk) + c360
              where (work_g2 > c360) work_g2 = work_g2 - c360
              where (work_g2 < c0 )  work_g2 = work_g2 + c360
              lonn_bounds(icorner,:,:,iblk) = work_g2(:,:)
              work_g2(:,:) = lone_bounds(icorner,:,:,iblk) + c360
              where (work_g2 > c360) work_g2 = work_g2 - c360
              where (work_g2 < c0 )  work_g2 = work_g2 + c360
              lone_bounds(icorner,:,:,iblk) = work_g2(:,:)
           enddo
        enddo
        deallocate(work_g2)
  
        end subroutine gridbox_edges

!=======================================================================

! NOTE:  Boundary conditions for fields on NW, SW, SE corners
!        have not been implemented; using NE corner location for all.
!        Extrapolations are also used: these fields are approximate!
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL

    subroutine gridbox_verts(work_g,vbounds)

        real (kind=dbl_kind), dimension(:,:), intent(in) :: &
            work_g
  
        real (kind=dbl_kind), dimension(4,nx_block,ny_block,max_blocks), intent(out) :: &
            vbounds
  
        integer (kind=int_kind) :: &
            i,j                 ! index counters
  
        real (kind=dbl_kind) :: &
            rad_to_deg
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g2
  
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
           work1
  
        character(len=*), parameter :: subname = '(gridbox_verts)'
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        if (my_task == master_task) then
           allocate(work_g2(nx_global,ny_global))
        else
           allocate(work_g2(1,1))
        endif
  
        !-------------------------------------------------------------
        ! Get coordinates of grid boxes for each block as follows:
        ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
        !-------------------------------------------------------------
  
        work_g2(:,:) = c0
        if (my_task == master_task) then
           do j = 2, ny_global
           do i = 2, nx_global
              work_g2(i,j) = work_g(i-1,j-1) * rad_to_deg
           enddo
           enddo
           ! extrapolate
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
           enddo
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g2(i,2) - work_g2(i,3)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        vbounds(1,:,:,:) = work1(:,:,:)
  
        work_g2(:,:) = c0
        if (my_task == master_task) then
           do j = 2, ny_global
           do i = 1, nx_global
              work_g2(i,j) = work_g(i,j-1) * rad_to_deg
           enddo
           enddo
           ! extrapolate
           do i = 1, nx_global
              work_g2(i,1) = (c2*work_g2(i,2) - work_g2(i,3))
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        vbounds(2,:,:,:) = work1(:,:,:)
  
        work_g2(:,:) = c0
        if (my_task == master_task) then
           do j = 1, ny_global
           do i = 1, nx_global
              work_g2(i,j) = work_g(i,j) * rad_to_deg
           enddo
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        vbounds(3,:,:,:) = work1(:,:,:)
  
        work_g2(:,:) = c0
        if (my_task == master_task) then
           do j = 1, ny_global
           do i = 2, nx_global
              work_g2(i,j) = work_g(i-1,j  ) * rad_to_deg
           enddo
           enddo
           ! extrapolate
           do j = 1, ny_global
              work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
           enddo
        endif
        call scatter_global(work1, work_g2, &
                            master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        vbounds(4,:,:,:) = work1(:,:,:)
  
        deallocate (work_g2)
  
    end subroutine gridbox_verts

end module ice_gridbox