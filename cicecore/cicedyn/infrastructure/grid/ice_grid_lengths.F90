module ice_grid_lengths

    use ice_kinds_mod
    use ice_constants, only: c0, c2, p5, p25, cm_to_m, &
        field_loc_center, field_loc_NEcorner, field_loc_Nface, field_loc_Eface, &
        field_type_scalar
    use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
    use ice_domain_size, only: nx_global, ny_global
    use ice_domain, only: blocks_ice, nblocks, halo_info, distrb_info
    use ice_communicate, only: my_task, master_task
    use ice_gather_scatter, only: gather_global, scatter_global

    implicit none
    private
    public :: primary_grid_lengths_HTN, primary_grid_lengths_HTE

    contains

!=======================================================================

! Calculate dxU and dxT from HTN on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dxU, dxT and HTN to all processors.
!
! author: Elizabeth C. Hunke, LANL

    subroutine primary_grid_lengths_HTN(work_g, save_ghte_ghtn, HTN, dxT, dxU, dyU, dxN, dxE, G_HTE, G_HTN)

        real (kind=dbl_kind), dimension(:,:), intent(inout)  :: work_g ! global array holding HTN

        logical (kind=log_kind), intent(in) :: save_ghte_ghtn

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            HTN    , & ! length of northern edge of T-cell (m)
            dxT    , & ! width of T-cell through the middle (m)
            dxU    , & ! width of U-cell through the middle (m)
            dyU    , & ! height of U-cell through the middle (m)
            dxN    , &     ! width of N-cell through the middle (m)
            dxE    ! width of E-cell through the middle (m)

        real (kind=dbl_kind), dimension (:,:), intent(out)  :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)
  
        ! local variables
  
        integer (kind=int_kind) :: &
           i, j, &
           ip1     ! i+1
  
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
           work_g2
  
        character(len=*), parameter :: subname = '(primary_grid_lengths_HTN)'
  
        if (my_task == master_task) then
           allocate(work_g2(nx_global,ny_global))
        else
           allocate(work_g2(1,1))
        endif
  
        ! HTN, dxU = average of 2 neighbor HTNs in i
  
        if (my_task == master_task) then
           do j = 1, ny_global
           do i = 1, nx_global
              work_g(i,j) = work_g(i,j) * cm_to_m                ! HTN
           enddo
           enddo
           do j = 1, ny_global
           do i = 1, nx_global
              ! assume cyclic; noncyclic will be handled during scatter
              ip1 = i+1
              if (i == nx_global) ip1 = 1
              work_g2(i,j) = p5*(work_g(i,j) + work_g(ip1,j))    ! dxU
           enddo
           enddo
           if (save_ghte_ghtn) then
              do j = 1, ny_global
              do i = 1,nx_global
                 G_HTN(i+nghost,j+nghost) = work_g(i,j)
              enddo
              enddo
              call global_ext_halo(G_HTN)
           endif
        endif
        call scatter_global(HTN, work_g, master_task, distrb_info, &
                            field_loc_Nface, field_type_scalar)
        call scatter_global(dxU, work_g2, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
  
        ! dxT = average of 2 neighbor HTNs in j
  
        if (my_task == master_task) then
           do j = 2, ny_global
           do i = 1, nx_global
              work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j-1)) ! dxT
           enddo
           enddo
           ! extrapolate to obtain dxT along j=1
           do i = 1, nx_global
              work_g2(i,1) = c2*work_g(i,2) - work_g(i,3) ! dxT
           enddo
        endif
        call scatter_global(dxT, work_g2, master_task, distrb_info, &
                            field_loc_center, field_type_scalar)
  
        ! dxN = HTN
  
        dxN(:,:,:) = HTN(:,:,:)   ! dxN
  
        ! dxE = average of 4 surrounding HTNs
  
        if (my_task == master_task) then
           do j = 2, ny_global
           do i = 1, nx_global
              ! assume cyclic; noncyclic will be handled during scatter
              ip1 = i+1
              if (i == nx_global) ip1 = 1
              work_g2(i,j) = p25*(work_g(i,j)+work_g(ip1,j)+work_g(i,j-1)+work_g(ip1,j-1))   ! dxE
           enddo
           enddo
           ! extrapolate to obtain dxT along j=1
           do i = 1, nx_global
              ! assume cyclic; noncyclic will be handled during scatter
              ip1 = i+1
              if (i == nx_global) ip1 = 1
              work_g2(i,1) = p5*(c2*work_g(i  ,2) - work_g(i  ,3) + &
                                 c2*work_g(ip1,2) - work_g(ip1,3))      ! dxE
           enddo
        endif
        call scatter_global(dxE, work_g2, master_task, distrb_info, &
                            field_loc_center, field_type_scalar)
  
        deallocate(work_g2)
  
    end subroutine primary_grid_lengths_HTN


!=======================================================================
! Calculate dyU and dyT from HTE on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dyU, dyT and HTE to all processors.
!
! author: Elizabeth C. Hunke, LANL

    subroutine primary_grid_lengths_HTE(work_g, save_ghte_ghtn, HTE, dyT, dyU, dyN, dyE, G_HTE, G_HTN)

        real (kind=dbl_kind), dimension(:,:), intent(inout)  :: work_g ! global array holding HTE

        logical (kind=log_kind), intent(in) :: save_ghte_ghtn

        real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
            HTE    , & ! length of northern edge of T-cell (m)
            dyT    , & ! height of T-cell through the middle (m)
            dyU    , & ! height of U-cell through the middle (m)
            dyN    , & ! height of N-cell through the middle (m)
            dyE        ! height of E-cell through the middle (m)

        real (kind=dbl_kind), dimension (:,:), intent(out)  :: &
            G_HTE  , & ! length of eastern edge of T-cell (global ext.)
            G_HTN      ! length of northern edge of T-cell (global ext.)
     
        ! local variables
     
        integer (kind=int_kind) :: &
            i, j, &
            im1     ! i-1
     
        real (kind=dbl_kind), dimension(:,:), allocatable :: &
            work_g2
     
        character(len=*), parameter :: subname = '(primary_grid_lengths_HTE)'
     
        if (my_task == master_task) then
            allocate(work_g2(nx_global,ny_global))
        else
            allocate(work_g2(1,1))
        endif
     
        ! HTE, dyU = average of 2 neighbor HTE in j
     
        if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
                work_g(i,j) = work_g(i,j) * cm_to_m                ! HTE
            enddo
            enddo
            do j = 1, ny_global-1
            do i = 1, nx_global
                work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j+1)) ! dyU
            enddo
            enddo
            ! extrapolate to obtain dyU along j=ny_global
            if (ny_global > 1) then
                do i = 1, nx_global
                    work_g2(i,ny_global) = c2*work_g(i,ny_global-1) - work_g(i,ny_global-2)  ! dyU
                enddo
            endif
            if (save_ghte_ghtn) then
                do j = 1, ny_global
                do i = 1, nx_global
                    G_HTE(i+nghost,j+nghost) = work_g(i,j)
                enddo
                enddo
                call global_ext_halo(G_HTE)
            endif
        endif
        call scatter_global(HTE, work_g, master_task, distrb_info, &
                            field_loc_Eface, field_type_scalar)
        call scatter_global(dyU, work_g2, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
     
        ! dyT = average of 2 neighbor HTE in i
     
        if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
                ! assume cyclic; noncyclic will be handled during scatter
                im1 = i-1
                if (i == 1) im1 = nx_global
                work_g2(i,j) = p5*(work_g(i,j) + work_g(im1,j))    ! dyT
            enddo
            enddo
        endif
        call scatter_global(dyT, work_g2, master_task, distrb_info, &
                            field_loc_center, field_type_scalar)
     
        ! dyN = average of 4 neighbor HTEs
     
        if (my_task == master_task) then
            do j = 1, ny_global-1
            do i = 1, nx_global
                ! assume cyclic; noncyclic will be handled during scatter
                im1 = i-1
                if (i == 1) im1 = nx_global
                work_g2(i,j) = p25*(work_g(i,j) + work_g(im1,j) + work_g(i,j+1) + work_g(im1,j+1))   ! dyN
            enddo
            enddo
            ! extrapolate to obtain dyN along j=ny_global
            if (ny_global > 1) then
                do i = 1, nx_global
                    ! assume cyclic; noncyclic will be handled during scatter
                    im1 = i-1
                    if (i == 1) im1 = nx_global
                    work_g2(i,ny_global) = p5*(c2*work_g(i  ,ny_global-1) - work_g(i  ,ny_global-2) + &
                                            c2*work_g(im1,ny_global-1) - work_g(im1,ny_global-2))     ! dyN
                enddo
            endif
        endif
        call scatter_global(dyN, work_g2, master_task, distrb_info, &
                            field_loc_center, field_type_scalar)
     
        ! dyE = HTE
     
        dyE(:,:,:) = HTE(:,:,:)
     
        deallocate(work_g2)
     
     end subroutine primary_grid_lengths_HTE

!=======================================================================

!  This subroutine fills ghost cells in global extended grid

    subroutine global_ext_halo(array)

        use ice_domain, only: ew_boundary_type, ns_boundary_type

        real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
           array   ! extended global grid size nx+2*nghost, ny+2*nghost
                   ! nghost+1:nghost+nx_global and nghost+1:nghost+ny_global filled on entry
  
        integer (kind=int_kind) :: n
  
        character(len=*), parameter :: subname = '(global_ext_halo)'
  
        do n = 1,nghost
           if (ns_boundary_type =='cyclic') then
              array(:,n)                  = array(:,ny_global+n)
              array(:,ny_global+nghost+n) = array(:,nghost+n)
           elseif (ns_boundary_type == 'open') then
              array(:,n)                  = array(:,nghost+1)
              array(:,ny_global+nghost+n) = array(:,ny_global+nghost)
           else
              array(:,n)                  = c0
              array(:,ny_global+nghost+n) = c0
           endif
        enddo
  
        do n = 1,nghost
           if (ew_boundary_type =='cyclic') then
              array(n                 ,:) = array(nx_global+n,:)
              array(nx_global+nghost+n,:) = array(nghost+n   ,:)
           elseif (ew_boundary_type == 'open') then
              array(n                 ,:) = array(nghost+1        ,:)
              array(nx_global+nghost+n,:) = array(nx_global+nghost,:)
           else
              array(n                 ,:) = c0
              array(nx_global+nghost+n,:) = c0
           endif
        enddo
  
        end subroutine global_ext_halo
  
  !=======================================================================

      
end module ice_grid_lengths