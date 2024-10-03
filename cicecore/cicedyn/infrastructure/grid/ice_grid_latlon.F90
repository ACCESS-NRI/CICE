module ice_grid_latlon

    use ice_kinds_mod
    use ice_blocks, only: block, get_block
    use ice_exit, only: abort_ice
    use ice_communicate, only: my_task, master_task
    use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop
    use ice_boundary, only: ice_HaloUpdate, ice_HaloExtrapolate
    use ice_domain, only: blocks_ice, nblocks, distrb_info, halo_info, &
        ew_boundary_type, ns_boundary_type
    use ice_fileunits, only: nu_diag
    use ice_constants, only: c0, c1, c1p5, c2, c4, p5,  &
        field_loc_center, field_loc_Nface, field_loc_Eface, &
        field_type_scalar
    use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted, &
        icepack_query_parameters, icepack_init_parameters


    implicit none
    private
    public :: Tlatlon, NElatlon

    contains

!=======================================================================

! Initializes latitude and longitude on T grid
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! generation routine

    subroutine Tlatlon(ULAT, ULON, grid_type, TLON, TLAT)


        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
        ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
        ULAT       ! latitude of velocity pts, NE corner of T pts (radians)

        character (len=char_len_long), intent(in)  :: &
        grid_type

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
        ! ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
        ! ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
        TLON   , & ! longitude of temp (T) pts (radians)
        TLAT    ! latitude of temp (T) pts (radians)
        ! NLON   , & ! longitude of center of north face of T pts (radians)
        ! ! NLAT   , & ! latitude of center of north face of T pts (radians)
        ! ELON   , & ! longitude of center of east face of T pts (radians)
        ! ELAT   , & ! latitude of center of east face of T pts (radians)

        integer (kind=int_kind) :: &
             i, j, iblk       , & ! horizontal indices
             ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        real (kind=dbl_kind) :: &
             z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da, &
             rad_to_deg
  
        type (block) :: &
             this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(Tlatlon)'
  
        if (my_task==master_task) then
           write(nu_diag,*) subname,' called'
        endif
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        TLAT(:,:,:) = c0
        TLON(:,:,:) = c0
  
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
        !$OMP                     x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
        !$OMP                     tx,ty,tz,da)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
  
              z1 = cos(ULAT(i-1,j-1,iblk))
              x1 = cos(ULON(i-1,j-1,iblk))*z1
              y1 = sin(ULON(i-1,j-1,iblk))*z1
              z1 = sin(ULAT(i-1,j-1,iblk))
  
              z2 = cos(ULAT(i,j-1,iblk))
              x2 = cos(ULON(i,j-1,iblk))*z2
              y2 = sin(ULON(i,j-1,iblk))*z2
              z2 = sin(ULAT(i,j-1,iblk))
  
              z3 = cos(ULAT(i-1,j,iblk))
              x3 = cos(ULON(i-1,j,iblk))*z3
              y3 = sin(ULON(i-1,j,iblk))*z3
              z3 = sin(ULAT(i-1,j,iblk))
  
              z4 = cos(ULAT(i,j,iblk))
              x4 = cos(ULON(i,j,iblk))*z4
              y4 = sin(ULON(i,j,iblk))*z4
              z4 = sin(ULAT(i,j,iblk))
  
              ! ---------
              ! TLON/TLAT 4 pt computation (pts 1, 2, 3, 4)
              ! ---------
  
              tx = (x1+x2+x3+x4)/c4
              ty = (y1+y2+y3+y4)/c4
              tz = (z1+z2+z3+z4)/c4
              da = sqrt(tx**2+ty**2+tz**2)
  
              tz = tz/da
  
              ! TLON in radians East
              TLON(i,j,iblk) = c0
              if (tx /= c0 .or. ty /= c0) TLON(i,j,iblk) = atan2(ty,tx)
  
              ! TLAT in radians North
              TLAT(i,j,iblk) = asin(tz)
  
           enddo                  ! i
           enddo                  ! j
        enddo                     ! iblk
        !$OMP END PARALLEL DO
  
        if (trim(grid_type) == 'regional') then
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
                    TLON(i,j,iblk) = c2*TLON(i+1,j,iblk) - &
                                        TLON(i+2,j,iblk)
                    TLAT(i,j,iblk) = c2*TLAT(i+1,j,iblk) - &
                                        TLAT(i+2,j,iblk)
                 enddo
              endif
           enddo
           !$OMP END PARALLEL DO
        endif   ! regional
  
        call ice_timer_start(timer_bound)
        call ice_HaloUpdate (TLON,             halo_info, &
                             field_loc_center, field_type_scalar, &
                             fillValue=c1)
        call ice_HaloUpdate (TLAT,             halo_info, &
                             field_loc_center, field_type_scalar, &
                             fillValue=c1)
        call ice_HaloExtrapolate(TLON, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
        call ice_HaloExtrapolate(TLAT, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
  
    end subroutine Tlatlon
  
  !=======================================================================
  
  ! Initializes latitude and longitude on N, E grid from Tlatlon
  !
  ! author: T. Craig 
  
    subroutine NElatlon(ULON, ULAT, TLON, TLAT, umask, tmask, emask, nmask, grid_type,  NLON, NLAT, ELON, ELAT)

        use ice_global_reductions, only: global_minval, global_maxval

        real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
            ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
            ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
            TLON   , & ! longitude of temp (T) pts (radians)
            TLAT    ! latitude of temp (T) pts (radians)

        logical (kind=log_kind), dimension (:,:,:), intent(in) :: &
            tmask  , & ! land/boundary mask, thickness (T-cell)
            umask  , & ! land/boundary mask  (U-cell) (1 if all surrounding T cells are ocean)
            ! umaskCD, & ! land/boundary mask  (U-cell) (1 if at least two surrounding T cells are ocean)
            nmask  , & ! land/boundary mask, (N-cell)
            emask   ! land/boundary mask, (E-cell)

        character (len=char_len_long), intent(in)  :: &
            grid_type

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            NLON   , & ! longitude of center of north face of T pts (radians)
            NLAT   , & ! latitude of center of north face of T pts (radians)
            ELON   , & ! longitude of center of east face of T pts (radians)
            ELAT       ! latitude of center of east face of T pts (radians)

        ! local vars
        integer (kind=int_kind) :: &
             i, j, iblk       , & ! horizontal indices
             ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        real (kind=dbl_kind) :: &
             z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da, &
             rad_to_deg
  
        type (block) :: &
             this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(NElatlon)'
  
        if (my_task==master_task) then
           write(nu_diag,*) subname,' called'
        endif
  
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
           file=__FILE__, line=__LINE__)
  
        NLAT(:,:,:) = c0
        NLON(:,:,:) = c0
        ELAT(:,:,:) = c0
        ELON(:,:,:) = c0
  
        !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
        !$OMP                     x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
        !$OMP                     tx,ty,tz,da)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi
  
           do j = jlo, jhi
           do i = ilo, ihi
  
              z1 = cos(ULAT(i-1,j-1,iblk))
              x1 = cos(ULON(i-1,j-1,iblk))*z1
              y1 = sin(ULON(i-1,j-1,iblk))*z1
              z1 = sin(ULAT(i-1,j-1,iblk))
  
              z2 = cos(ULAT(i,j-1,iblk))
              x2 = cos(ULON(i,j-1,iblk))*z2
              y2 = sin(ULON(i,j-1,iblk))*z2
              z2 = sin(ULAT(i,j-1,iblk))
  
              z3 = cos(ULAT(i-1,j,iblk))
              x3 = cos(ULON(i-1,j,iblk))*z3
              y3 = sin(ULON(i-1,j,iblk))*z3
              z3 = sin(ULAT(i-1,j,iblk))
  
              z4 = cos(ULAT(i,j,iblk))
              x4 = cos(ULON(i,j,iblk))*z4
              y4 = sin(ULON(i,j,iblk))*z4
              z4 = sin(ULAT(i,j,iblk))
  
              ! ---------
              ! NLON/NLAT 2 pt computation (pts 3, 4)
              ! ---------
  
              tx = (x3+x4)/c2
              ty = (y3+y4)/c2
              tz = (z3+z4)/c2
              da = sqrt(tx**2+ty**2+tz**2)
  
              tz = tz/da
  
              ! NLON in radians East
              NLON(i,j,iblk) = c0
              if (tx /= c0 .or. ty /= c0) NLON(i,j,iblk) = atan2(ty,tx)
  
              ! NLAT in radians North
              NLAT(i,j,iblk) = asin(tz)
  
              ! ---------
              ! ELON/ELAT 2 pt computation (pts 2, 4)
              ! ---------
  
              tx = (x2+x4)/c2
              ty = (y2+y4)/c2
              tz = (z2+z4)/c2
              da = sqrt(tx**2+ty**2+tz**2)
  
              tz = tz/da
  
              ! ELON in radians East
              ELON(i,j,iblk) = c0
              if (tx /= c0 .or. ty /= c0) ELON(i,j,iblk) = atan2(ty,tx)
  
              ! ELAT in radians North
              ELAT(i,j,iblk) = asin(tz)
  
           enddo                  ! i
           enddo                  ! j
        enddo                     ! iblk
        !$OMP END PARALLEL DO
  
        if (trim(grid_type) == 'regional') then
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
                    NLON(i,j,iblk) = c1p5*TLON(i+1,j,iblk) - &
                                       p5*TLON(i+2,j,iblk)
                    NLAT(i,j,iblk) = c1p5*TLAT(i+1,j,iblk) - &
                                       p5*TLAT(i+2,j,iblk)
                 enddo
              endif
           enddo
           !$OMP END PARALLEL DO
        endif   ! regional
  
        call ice_timer_start(timer_bound)
        call ice_HaloUpdate (NLON,             halo_info, &
                             field_loc_Nface,  field_type_scalar, &
                             fillValue=c1)
        call ice_HaloUpdate (NLAT,             halo_info, &
                             field_loc_Nface,  field_type_scalar, &
                             fillValue=c1)
        call ice_HaloUpdate (ELON,             halo_info, &
                             field_loc_Eface,  field_type_scalar, &
                             fillValue=c1)
        call ice_HaloUpdate (ELAT,             halo_info, &
                             field_loc_Eface,  field_type_scalar, &
                             fillValue=c1)
        call ice_HaloExtrapolate(NLON, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
        call ice_HaloExtrapolate(NLAT, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
        call ice_HaloExtrapolate(ELON, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
        call ice_HaloExtrapolate(ELAT, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)
        call ice_timer_stop(timer_bound)
  
        x1 = global_minval(TLON, distrb_info, tmask)
        x2 = global_maxval(TLON, distrb_info, tmask)
        x3 = global_minval(TLAT, distrb_info, tmask)
        x4 = global_maxval(TLAT, distrb_info, tmask)
  
        y1 = global_minval(ULON, distrb_info, umask)
        y2 = global_maxval(ULON, distrb_info, umask)
        y3 = global_minval(ULAT, distrb_info, umask)
        y4 = global_maxval(ULAT, distrb_info, umask)
  
        if (my_task==master_task) then
           write(nu_diag,*) ' '
           write(nu_diag,*) subname,' min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
           write(nu_diag,*) subname,' min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
           write(nu_diag,*) subname,' min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
           write(nu_diag,*) subname,' min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
        endif                     ! my_task
  
        x1 = global_minval(NLON, distrb_info, nmask)
        x2 = global_maxval(NLON, distrb_info, nmask)
        x3 = global_minval(NLAT, distrb_info, nmask)
        x4 = global_maxval(NLAT, distrb_info, nmask)
  
        y1 = global_minval(ELON, distrb_info, emask)
        y2 = global_maxval(ELON, distrb_info, emask)
        y3 = global_minval(ELAT, distrb_info, emask)
        y4 = global_maxval(ELAT, distrb_info, emask)
  
        if (my_task==master_task) then
           write(nu_diag,*) ' '
           write(nu_diag,*) subname,' min/max NLON:', x1*rad_to_deg, x2*rad_to_deg
           write(nu_diag,*) subname,' min/max NLAT:', x3*rad_to_deg, x4*rad_to_deg
           write(nu_diag,*) subname,' min/max ELON:', y1*rad_to_deg, y2*rad_to_deg
           write(nu_diag,*) subname,' min/max ELAT:', y3*rad_to_deg, y4*rad_to_deg
        endif                     ! my_task
  
    end subroutine NElatlon

end module ice_grid_latlon