
      module gathscatchk_data

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind, log_kind
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot, nghost, &
          i_global, j_global, nblocks_x, nblocks_y, &
          ew_boundary_type, ns_boundary_type
      use ice_boundary, only: ice_HaloUpdate, ice_HaloUpdate_stress
      use ice_gather_scatter
      use ice_constants, only: c0, c1, c2, p5, spval_dbl, spval_int, &
          field_loc_unknown, field_loc_noupdate, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, &
          field_type_unknown, field_type_noupdate, &
          field_type_scalar, field_type_vector, field_type_angle
      use ice_communicate, only: my_task, master_task, get_num_procs, MPI_COMM_ICE
      use ice_distribution, only: ice_distributionGetBlockID, ice_distributionGet, &
          ice_distributionGetBlockLoc
      use ice_domain_size, only: nx_global, ny_global, &
          block_size_x, block_size_y, max_blocks
      use ice_domain, only: distrb_info, halo_info, nblocks
      use ice_exit, only: abort_ice, end_run
      use ice_global_reductions, only: global_minval, global_maxval, global_sum

      implicit none

      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1, &
         skipflag = 2

      end module gathscatchk_data

!=======================================================================

      program gathscatchk

      ! This tests the CICE halo update methods by
      ! using CICE_InitMod (from the standalone model) to read/initialize
      ! a CICE grid/configuration.

      use gathscatchk_data

      implicit none

      integer(int_kind) :: nn, nl, nt, nf, n
      integer(int_kind) :: i, j, ii, jj, i1, j1, i2, j2, i3, j3
      integer(int_kind) :: ig, jg, ilo, ihi, jlo, jhi
      integer(int_kind) :: nxg, nyg, nxgx, nygx, nxb, nyb
      integer(int_kind) :: iblock, itrip, ioffset, joffset
      integer(int_kind) :: blockID, numBlocks, jtrip, processor
      type (block) :: this_block

      ! temporary fields sent for computation
      real(dbl_kind)   , allocatable :: dg1(:,:), dg1x(:,:), da1(:,:,:)
      real(dbl_kind)   , allocatable :: dg2(:,:), dg2x(:,:), da2(:,:,:)
      real(real_kind)  , allocatable :: rg1(:,:), rg1x(:,:), ra1(:,:,:)
      integer(int_kind), allocatable :: ig1(:,:), ig1x(:,:), ia1(:,:,:)
      logical(log_kind), allocatable :: lg1(:,:), lg1x(:,:), la1(:,:,:)

      ! baseline array values
      real(dbl_kind), allocatable :: bg1x(:,:),ba1(:,:,:)

      ! expected results
      real(dbl_kind), allocatable :: bgxa(:,:),baa(:,:,:)  ! extended
      real(dbl_kind), allocatable :: bgxb(:,:),bab(:,:,:)  ! land block elim
      real(dbl_kind), allocatable :: bgxc(:,:),bac(:,:,:)  ! fully corrected

      integer(int_kind), parameter :: maxtests = 5
      integer(int_kind), parameter :: maxtypes = 0
      integer(int_kind), parameter :: maxlocs = 0
      integer(int_kind), parameter :: maxfills = 0
      integer(int_kind), parameter :: nz1 = 0
      integer(int_kind), parameter :: nz2 = 0
      real(dbl_kind)    :: aichk,ajchk,cichk,cjchk,rival,rjval,rsign
      real(dbl_kind)    :: fillexpected,cichk1,cichk2,cjchk1,cjchk2,wgt1,wgt2
      character(len=16) :: locs_name(maxlocs), types_name(maxtypes), fill_name(maxfills)
      integer(int_kind) :: field_loc(maxlocs), field_type(maxtypes)
      logical :: halofill, found
      integer(int_kind) :: npes, ierr, ntask, testcnt, tottest, tpcnt, tfcnt, tscnt
      integer(int_kind) :: errorflag0, gflag, k1m, k2m, ptcntsum, failcntsum
      integer(int_kind) :: nblocks_used
      integer(int_kind), allocatable :: errorflag(:)
      integer(int_kind), allocatable :: ptcnt(:), failcnt(:)
      character(len=128), allocatable :: teststring(:)
      character(len=32) :: halofld
      logical :: lbe  ! land block elimination
      logical :: tripole_average, tripole_pole
      logical :: first_call = .true.

      real(dbl_kind)   , parameter :: dgsspval = -9.99e25_dbl_kind
      real(dbl_kind)   , parameter :: fillval = -88888.0_dbl_kind
      character(len=*) , parameter :: subname='(gathscatchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize
      npes = get_num_procs()
      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      nblocks_used = global_sum(nblocks, distrb_info)
      if (nblocks_used < nblocks_tot) then
         lbe = .true.
      else
         lbe = .false.
      endif

      tottest = maxtests
      allocate(errorflag(tottest))
      allocate(teststring(tottest))
      allocate(ptcnt(tottest))
      allocate(failcnt(tottest))
      ptcnt(:) = 0
      failcnt(:) = 0

      !-----------------------------------------------------------------
      ! Testing
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest GATHSCATCHK'
         write(6,*) ' '
         write(6,*) ' npes         = ',npes
         write(6,*) ' my_task      = ',my_task
         write(6,*) ' nx_global    = ',nx_global
         write(6,*) ' ny_global    = ',ny_global
         write(6,*) ' block_size_x = ',block_size_x
         write(6,*) ' block_size_y = ',block_size_y
         write(6,*) ' nblocks_used = ',nblocks_used
         write(6,*) ' nblocks_tot  = ',nblocks_tot
         write(6,*) ' tottest      = ',tottest
         if (lbe) then
            write(6,*) ' NOTE: lbe is TRUE, haloUpdate check will be skipped '
         endif
         write(6,*) ' '
      endif

      teststring(:) = ' '

      !-----------------------------------------------------------------
      ! Test gathscat
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         nxg = nx_global
         nyg = ny_global
         nxgx = nx_global+2*nghost
         nygx = ny_global+2*nghost
      else
         nxg = 1
         nyg = 1
         nxgx = 1
         nygx = 1
      endif
      nxb = nx_block
      nyb = ny_block
      allocate(dg1(nxg,nyg), dg1x(nxgx,nygx), da1(nx_block,ny_block,max_blocks))
      allocate(dg2(nxg,nyg), dg2x(nxgx,nygx), da2(nx_block,ny_block,max_blocks))
      allocate(rg1(nxg,nyg), rg1x(nxgx,nygx), ra1(nx_block,ny_block,max_blocks))
      allocate(ig1(nxg,nyg), ig1x(nxgx,nygx), ia1(nx_block,ny_block,max_blocks))
      allocate(lg1(nxg,nyg), lg1x(nxgx,nygx), la1(nx_block,ny_block,max_blocks))

      allocate(              bg1x(nxgx,nygx), ba1(nx_block,ny_block,max_blocks))

      allocate(              bgxa(nxgx,nygx), baa(nx_block,ny_block,max_blocks))
      allocate(              bgxb(nxgx,nygx), bab(nx_block,ny_block,max_blocks))
      allocate(              bgxc(nxgx,nygx), bac(nx_block,ny_block,max_blocks))

      bgxa(:,:) = spval_dbl ! expected result (iglobal,jglobal)
      bgxb(:,:) = spval_dbl ! expected result (iglobal,jglobal)
      bgxc(:,:) = spval_dbl ! expected result corrected for bc
      bg1x(:,:) = spval_dbl ! baseline test field (1:nxg,1:nyg)

   if (my_task == master_task) then
      ! i_global is 1-nghost:n*_global+nghost
      ! local extended grids are 1:n*_global+2*nghost

      do j = 1,nygx
      do i = 1,nxgx
         ! bg1x is global by actual i,j index, no messing around
         ii = i-nghost
         jj = j-nghost
         bg1x(i,j) = real(ii*1000+jj,dbl_kind)

         ! bgxa is bg1x but using i_global, j_global index and land
         ! block elimination
         ! abs is used here because j_global = -j_global for tripole
         ! tcx -j_global should go away soon
         ii = i_global(i-nghost)
         jj = j_global(j-nghost)
         if (jj < 1-nghost) jj = abs(jj)
         bgxa(i,j) = real(ii*1000+jj,dbl_kind)

         ! bgxb is bg1x with spval in eliminated blocks
         bgxb(i,j) = bg1x(i,j)

         ! bgxc is bgxa + full halo calc, fill, including tripole
         bgxc(i,j) = bgxa(i,j)
         ! undefined
         if (ii < 1 .or. ii > nx_global) bgxc(i,j) = spval_dbl
         if (jj < 1 .or. jj > ny_global) bgxc(i,j) = spval_dbl
         ! tripole
         if (j > ny_global .and. j <= ny_global+nghost .and. &
            (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then
            if (ns_boundary_type == 'tripole') then
               ii = nx_global - ii + 1
               jj = 2*ny_global - jj + 1
            elseif (ns_boundary_type == 'tripoleT') then
               ii = nx_global - ii + 2
               jj = 2*ny_global - jj
            endif
            if (ii < 1)         ii = ii + nx_global
            if (ii > nx_global) ii = ii - nx_global
            bgxc(i,j) = real(ii*1000+jj,dbl_kind)
         elseif (j == ny_global+nghost .and. ns_boundary_type == 'tripoleT') then
            ! average current value and shifted value on boundary
            ii = nx_global - ii + 2
            jj = 2*ny_global - jj
            if (ii < 1)         ii = ii + nx_global
            if (ii > nx_global) ii = ii - nx_global
            bgxc(i,j) = (bgxc(i,j) + real(ii*1000+jj,dbl_kind))/c2
         elseif (i == 1 .and. &
            (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap')) then
            ! inside neighbor
            ii = i_global(i+1-nghost)
            if (j == 1) then
               jj = j_global(j+1-nghost)
            elseif (j == nygx) then
               jj = j_global(j-1-nghost)
            else
               jj = j_global(j-nghost)
            endif
            bgxc(i,j) = real(ii*1000+jj,dbl_kind)
            ! combine with second neighbor
            if (ew_boundary_type == 'linear_extrap') then
               ii = i_global(i+2-nghost)
               if (j == 1) then
                  jj = j_global(j+2-nghost)
               elseif (j == nygx) then
                  jj = j_global(j-2-nghost)
               else
                  jj = j_global(j-nghost)
               endif
               bgxc(i,j) = c2*bgxc(i,j)-real(ii*1000+jj,dbl_kind)
            endif
         elseif (i == nxgx .and. &
            (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap')) then
            ! inside neighbor
            ii = i_global(i-1-nghost)
            if (j == 1) then
               jj = j_global(j+1-nghost)
            elseif (j == nygx) then
               jj = j_global(j-1-nghost)
            else
               jj = j_global(j-nghost)
            endif
            bgxc(i,j) = real(ii*1000+jj,dbl_kind)
            ! combine with second neighbor
            if (ew_boundary_type == 'linear_extrap') then
               ii = i_global(i-2-nghost)
               if (j == 1) then
                  jj = j_global(j+2-nghost)
               elseif (j == nygx) then
                  jj = j_global(j-2-nghost)
               else
                  jj = j_global(j-nghost)
               endif
               bgxc(i,j) = c2*bgxc(i,j)-real(ii*1000+jj,dbl_kind)
            endif
         elseif (j == 1 .and. &
            (ns_boundary_type == 'zero_gradient' .or. ns_boundary_type == 'linear_extrap')) then
            ! inside neighbor
            ii = i_global(i-nghost)
            jj = j_global(j+1-nghost)
            bgxc(i,j) = real(ii*1000+jj,dbl_kind)
            ! combine with second neighbor
            if (ns_boundary_type == 'linear_extrap') then
               ii = i_global(i-nghost)
               jj = j_global(j+2-nghost)
               bgxc(i,j) = c2*bgxc(i,j)-real(ii*1000+jj,dbl_kind)
            endif
         elseif (j == nygx .and. &
            (ns_boundary_type == 'zero_gradient' .or. ns_boundary_type == 'linear_extrap')) then
            ! inside neighbor
            ii = i_global(i-nghost)
            jj = j_global(j-1-nghost)
            bgxc(i,j) = real(ii*1000+jj,dbl_kind)
            ! combine with second neighbor
            if (ns_boundary_type == 'linear_extrap') then
               ii = i_global(i-nghost)
               jj = j_global(j-2-nghost)
               bgxc(i,j) = c2*bgxc(i,j)-real(ii*1000+jj,dbl_kind)
            endif
         endif
      enddo
      enddo

      ! find eliminated blocks and set those gridcells to spval including on outer halo
      do blockID = 1,nblocks_tot
         call ice_distributionGetBlockLoc(distrb_info, blockID, processor, iblock)
         if (iblock == 0) then
!            write(6,*) 'Eliminate Block ',blockID
            this_block = get_block(blockID, blockID)
            i1 = this_block%ilo
            i2 = this_block%ihi
            j1 = this_block%jlo
            j2 = this_block%jhi
            do j = j1, j2
            do i = i1, i2
               ! find global indices of eliminated blocks
               ii = this_block%i_glob(i)
               jj = this_block%j_glob(j)
               ! then look for those indices in the global grid (need to consider halo)
               do j3 = 1,nygx
               do i3 = 1,nxgx
                  if (i_global(i3-nghost) == ii .and. j_global(j3-nghost) == jj) then
!                     write(6,*) 'Elimination gridpointa ',i3-nghost,j3-nghost
                     bgxa(i3,j3) = spval_dbl
                     bgxc(i3,j3) = spval_dbl
                  endif
               enddo
               enddo
            enddo
            enddo

            i1 = this_block%ilo
            i2 = this_block%ihi
            j1 = this_block%jlo
            j2 = this_block%jhi
            if (this_block%iblock == 1        ) i1 = this_block%ilo-nghost
            if (this_block%iblock == nblocks_x) i2 = this_block%ihi+nghost
            if (this_block%jblock == 1        ) j1 = this_block%jlo-nghost
            if (this_block%jblock == nblocks_y) j2 = this_block%jhi+nghost
            do j = j1, j2
            do i = i1, i2
!               ii = this_block%i_glob(this_block%ilo) + i - nghost - 1
!               jj = this_block%j_glob(this_block%jlo) + j - nghost - 1
               ii = this_block%i_glob(this_block%ilo) + i - this_block%ilo
               jj = this_block%j_glob(this_block%jlo) + j - this_block%jlo
!               write(6,'(a,6i6)') 'Elimination gridpointb ',i,j,ii,jj,ii+nghost,jj+nghost
               ! for bgxb, just eliminate block and outer halo
               bgxb(ii+nghost,jj+nghost) = spval_dbl
            enddo
            enddo
         endif
      enddo

   endif

      baa(:,:,:) = spval_dbl ! expect result (i_glob,j_glob)
      bac(:,:,:) = spval_dbl ! expect result corrected for bc
      ba1(:,:,:) = spval_dbl ! baseline test field (1:nxb,1:nyb)
      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         do j = 1,ny_block
         do i = 1,nx_block
            ii = this_block%i_glob(this_block%ilo) + i - nghost - 1
            jj = this_block%j_glob(this_block%jlo) + j - nghost - 1
            ba1(i,j,iblock) = real(ii*1000+jj,dbl_kind)
            ! padding
            if (ii < 1-nghost .or. ii > nx_global+nghost) ba1(i,j,iblock) = spval_dbl
            if (jj < 1-nghost .or. jj > ny_global+nghost) ba1(i,j,iblock) = spval_dbl

            ! abs is used here because j_glob = -j_glob for tripole
            ! tcx this should go away
            ii = this_block%i_glob(i)
            jj = this_block%j_glob(j)
            if (jj < 1-nghost) jj = abs(jj)
            baa(i,j,iblock) = real(ii*1000+jj,dbl_kind)

            bac(i,j,iblock) = baa(i,j,iblock)
            ! undefined
            if (ii < 1 .or. ii > nx_global) bac(i,j,iblock) = spval_dbl
            if (jj < 1 .or. jj > ny_global) bac(i,j,iblock) = spval_dbl

            if (this_block%i_glob(i) == spval_int .or. this_block%j_glob(j) == spval_int) then
                ! skip fill, padded blocks
            elseif (this_block%jblock == nblocks_y .and. &
                j > this_block%jhi .and. j <= this_block%jhi+nghost .and. &
                (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then
               if (ns_boundary_type == 'tripole') then
                  ii = nx_global - ii + 1
                  jj = 2*ny_global - jj + 1
               elseif (ns_boundary_type == 'tripoleT') then
                  ii = nx_global - ii + 2
                  jj = 2*ny_global - jj
               endif
               if (ii < 1)         ii = ii + nx_global
               if (ii > nx_global) ii = ii - nx_global
               bac(i,j,iblock) = real(ii*1000+jj,dbl_kind)
            elseif (this_block%jblock == nblocks_y .and. &
                    j == this_block%jhi .and. ns_boundary_type == 'tripoleT') then
               ! average current value and shifted value on boundary
               ii = nx_global - ii + 2
               jj = 2*ny_global - jj
               if (ii < 1)         ii = ii + nx_global
               if (ii > nx_global) ii = ii - nx_global
               bac(i,j,iblock) = (bac(i,j,iblock) + real(ii*1000+jj,dbl_kind))/c2
!tcx, this seems to have a problem with nghost>1
            elseif (this_block%iblock == 1 .and. i < this_block%ilo .and. &
               (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap')) then
               ! inside neighbor
               ii = this_block%i_glob(i+1)
               if (this_block%jblock == 1 .and. j < this_block%jlo) then
                  jj = this_block%j_glob(j+1)
               elseif (this_block%jblock == nblocks_y .and. &
                  j > this_block%jhi .and. j <= this_block%jhi+nghost) then
                  jj = this_block%j_glob(j-1)
               else
                  jj = this_block%j_glob(j)
               endif
               bac(i,j,iblock) = real(ii*1000+jj,dbl_kind)
               ! combine with second neighbor
               if (ew_boundary_type == 'linear_extrap') then
                  ii = this_block%i_glob(i+2)
                  if (this_block%jblock == 1 .and. j < this_block%jlo) then
                     jj = this_block%j_glob(j+2)
                  elseif (this_block%jblock == nblocks_y .and. &
                     j > this_block%jhi .and. j <= this_block%jhi+nghost) then
                     jj = this_block%j_glob(j-2)
                  else
                     jj = this_block%j_glob(j)
                  endif
                  bac(i,j,iblock) = c2*bac(i,j,iblock)-real(ii*1000+jj,dbl_kind)
               endif
            elseif (this_block%iblock == nblocks_x .and. &
               i > this_block%ihi .and. i <= this_block%ihi+nghost .and. &
               (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap')) then
               ! inside neighbor
               ii = this_block%i_glob(i-1)
               if (this_block%jblock == 1 .and. j < this_block%jlo) then
                  jj = this_block%j_glob(j+1)
               elseif (this_block%jblock == nblocks_y .and. &
                  j > this_block%jhi .and. j <= this_block%jhi+nghost) then
                  jj = this_block%j_glob(j-1)
               else
                  jj = this_block%j_glob(j)
               endif
               bac(i,j,iblock) = real(ii*1000+jj,dbl_kind)
               ! combine with second neighbor
               if (ew_boundary_type == 'linear_extrap') then
                  ii = this_block%i_glob(i-2)
                  if (this_block%jblock == 1 .and. j < this_block%jlo) then
                     jj = this_block%j_glob(j+2)
                  elseif (this_block%jblock == nblocks_y .and. &
                     j > this_block%jhi .and. j <= this_block%jhi+nghost) then
                     jj = this_block%j_glob(j-2)
                  else
                     jj = this_block%j_glob(j)
                  endif
                  bac(i,j,iblock) = c2*bac(i,j,iblock)-real(ii*1000+jj,dbl_kind)
               endif
            elseif (this_block%jblock == 1 .and. j < this_block%jlo .and. &
               (ns_boundary_type == 'zero_gradient' .or. ns_boundary_type == 'linear_extrap')) then
               ! inside neighbor
               ii = this_block%i_glob(i)
               jj = this_block%j_glob(j+1)
               bac(i,j,iblock) = real(ii*1000+jj,dbl_kind)
               ! combine with second neighbor
               if (ns_boundary_type == 'linear_extrap') then
                  ii = this_block%i_glob(i)
                  jj = this_block%j_glob(j+2)
                  bac(i,j,iblock) = c2*bac(i,j,iblock)-real(ii*1000+jj,dbl_kind)
               endif
            elseif (this_block%jblock == nblocks_y .and. &
               j > this_block%jhi .and. j <= this_block%jhi+nghost .and. &
               (ns_boundary_type == 'zero_gradient' .or. ns_boundary_type == 'linear_extrap')) then
               ! inside neighbor
               ii = this_block%i_glob(i)
               jj = this_block%j_glob(j-1)
               bac(i,j,iblock) = real(ii*1000+jj,dbl_kind)
               ! combine with second neighbor
               if (ns_boundary_type == 'linear_extrap') then
                  ii = this_block%i_glob(i)
                  jj = this_block%j_glob(j-2)
                  bac(i,j,iblock) = c2*bac(i,j,iblock)-real(ii*1000+jj,dbl_kind)
               endif
            endif
         enddo
         enddo
      enddo

      errorflag(:)  = passflag
      do n = 1,maxtests

         if (n == 1) then
            teststring(n) = 'GATH_EXT'
            da1(:,:,:) = ba1(:,:,:)
            dg1x(:,:) = spval_dbl
            call gather_global(dg1x,da1,master_task,distrb_info,grid_ext=.true.)
            if (my_task == master_task) then
               do j = 1,nygx
               do i = 1,nxgx
                  ptcnt(n) = ptcnt(n) + 1
                  if (dg1x(i,j) /= bgxb(i,j)) then
                     errorflag(n) = failflag
                     failcnt(n) = failcnt(n) + 1
                     write(100+my_task,'(a,3i6,2g16.6)') 'e1 ',n,i,j,dg1x(i,j),bgxb(i,j)
                  endif
               enddo
               enddo
            endif

         elseif (n == 2) then
            teststring(n) = 'GATH_STD'
            da1(:,:,:) = ba1(:,:,:)
            dg1(:,:) = spval_dbl
            call gather_global(dg1,da1,master_task,distrb_info)
            if (my_task == master_task) then
               do j = 1,nyg
               do i = 1,nxg
                  ptcnt(n) = ptcnt(n) + 1
                  if (dg1(i,j) /= bgxa(i+nghost,j+nghost)) then
                     errorflag(n) = failflag
                     failcnt(n) = failcnt(n) + 1
                     write(100+my_task,'(a,3i6,2g16.6)') 'e2 ',n,i,j,dg1(i,j),bgxa(i+nghost,j+nghost)
                  endif
               enddo
               enddo
            endif

         elseif (n == 3) then
            teststring(n) = 'SCAT_EXT'
            dg1x(:,:) = bg1x(:,:)
            da1(:,:,:) = spval_dbl
            call scatter_global(da1,dg1x,master_task,distrb_info,fillValue=spval_dbl,grid_ext=.true.)
            do iblock = 1,numBlocks
               do j = 1,ny_block
               do i = 1,nx_block
                  ptcnt(n) = ptcnt(n) + 1
                  if (da1(i,j,iblock) /= ba1(i,j,iblock)) then
                     errorflag(n) = failflag
                     failcnt(n) = failcnt(n) + 1
                     write(100+my_task,'(a,4i6,2g16.6)') 'e3 ',n,i,j,iblock,da1(i,j,iblock),ba1(i,j,iblock)
                  endif
               enddo
               enddo
            enddo

         elseif (n == 4) then
            teststring(n) = 'SCAT_STD'
            if (my_task == master_task) then
               dg1(1:nx_global,1:ny_global) = bg1x(1+nghost:nx_global+nghost,1+nghost:ny_global+nghost)
            endif
            da1(:,:,:) = spval_dbl
            call scatter_global(da1,dg1,master_task,distrb_info,field_loc_center,field_type_scalar,fillValue=spval_dbl)
            do iblock = 1,numBlocks
               do j = 1,ny_block
               do i = 1,nx_block
                  ptcnt(n) = ptcnt(n) + 1
                  if (da1(i,j,iblock) /= bac(i,j,iblock)) then
                     errorflag(n) = failflag
                     failcnt(n) = failcnt(n) + 1
                     write(100+my_task,'(a,4i6,2g16.6)') 'e4 ',n,i,j,iblock,da1(i,j,iblock),bac(i,j,iblock)
                  endif
               enddo
               enddo
            enddo

         elseif (n == 5) then
            teststring(n) = 'GATHHALO'
            da2(:,:,:) = spval_dbl
            do iblock = 1,numBlocks
               call ice_distributionGetBlockID(distrb_info, iblock, blockID)
               this_block = get_block(blockID, blockID)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo,jhi
               do i = ilo,ihi
                  da2(i,j,iblock) = da1(i,j,iblock)
               enddo
               enddo
            enddo
            call ice_haloUpdate(da2, halo_info, field_loc_center, field_type_scalar, fillvalue=spval_dbl)
            do iblock = 1,numBlocks
               do j = 1,ny_block
               do i = 1,nx_block
                  ptcnt(n) = ptcnt(n) + 1
                  if (da2(i,j,iblock) /= da1(i,j,iblock)) then
                     if (lbe) then
                        errorflag(n) = skipflag
                     else
                        errorflag(n) = failflag
                     endif
                     failcnt(n) = failcnt(n) + 1
                     write(100+my_task,'(a,4i6,2g16.6)') 'e5 ',n,i,j,iblock,da2(i,j,iblock),da1(i,j,iblock)
                  endif
               enddo
               enddo
            enddo

         endif

      enddo

      ! ---------------------------
      ! SUMMARY
      ! ---------------------------

      errorflag0 = passflag
      do n = 1,tottest
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
         ptcntsum = global_sum(ptcnt(n),distrb_info)
         ptcnt(n) = ptcntsum
         failcntsum = global_sum(failcnt(n),distrb_info)
         failcnt(n) = failcntsum
         if (errorflag(n) == failflag) errorflag0 = failflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GATHSCATCHK COMPLETED SUCCESSFULLY'
         write(6,*) ' '
         tpcnt = 0
         tfcnt = 0
         tscnt = 0
         do n = 1,tottest
            if (errorflag(n) == passflag) then
               tpcnt = tpcnt + 1
               write(6,'(2a,2i9)') 'PASS ',trim(teststring(n)),ptcnt(n),failcnt(n)
            elseif (errorflag(n) == skipflag) then
               tscnt = tscnt + 1
               write(6,'(2a,2i9,a)') 'SKIP ',trim(teststring(n)),ptcnt(n),failcnt(n),' (lbe active)'
            else
               tfcnt = tfcnt + 1
               write(6,'(2a,2i9)') 'FAIL ',trim(teststring(n)),ptcnt(n),failcnt(n)
            endif
         enddo
         write(6,*) ' '
         write(6,*) ' total skip = ',tscnt
         write(6,*) ' total pass = ',tpcnt
         write(6,*) ' total fail = ',tfcnt
         write(6,*) ' '
         if (errorflag0 == passflag) then
            write(6,*) 'GATHSCATCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'GATHSCATCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif

      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program gathscatchk

!=======================================================================
