
      module gathscatchk_data

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind, real_kind, log_kind
      use ice_blocks, only: block, get_block, nx_block, ny_block, nblocks_tot, nghost, &
          i_global, j_global, nblocks_x, nblocks_y
      use ice_boundary, only: ice_HaloUpdate, ice_HaloUpdate_stress
      use ice_gather_scatter
      use ice_constants, only: c0, c1, c2, p5, spval_dbl, &
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
      use ice_domain, only: distrb_info, halo_info, nblocks, &
          ew_boundary_type, ns_boundary_type
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

         ! bgxa is is bg1x but using i_global, j_global index and land
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
            ! tripole
            if (this_block%jblock == nblocks_y .and. &
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
#if (1 == 0)
      locs_name (:) = 'unknown'
      types_name(:) = 'unknown'
      fill_name (:) = 'unknown'
      field_type(:) = field_type_unknown
      field_loc (:) = field_loc_unknown

      types_name(1) = 'scalar'
      field_type(1) = field_type_scalar
      types_name(2) = 'vector'
      field_type(2) = field_type_vector
      types_name(3) = 'angle'
      field_type(3) = field_type_angle
      types_name(4) = 'none'
      field_type(4) = field_type_noupdate
!      types_name(5) = 'unknown'
!      field_type(5) = field_type_unknown  ! aborts in CICE, as expected

      locs_name (1) = 'center'
      field_loc (1)  = field_loc_center
      locs_name (2) = 'NEcorn'
      field_loc (2)  = field_loc_NEcorner
      locs_name (3) = 'Nface'
      field_loc (3)  = field_loc_Nface
      locs_name (4) = 'Eface'
      field_loc (4)  = field_loc_Eface
      locs_name (5) = 'none'
      field_loc (5)  = field_loc_noupdate
!      locs_name (6) = 'unknown'
!      field_loc (6)  = field_loc_unknown  ! aborts in CICE, as expected

      fill_name (1) = 'fill'
      fill_name (2) = 'nofill'

      tottest = maxtests * maxlocs * maxtypes * maxfills
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
         write(6,*) ' nblocks_tot  = ',nblocks_tot
         write(6,*) ' tottest      = ',tottest
         write(6,*) ' '
      endif

      errorflag0    = passflag
      errorflag(:)  = passflag
      teststring(:) = ' '

      ! ---------------------------
      ! TEST HALO UPDATE
      ! ---------------------------

      if (my_task == master_task) write(6,*) ' '

      allocate(darrayi1   (nx_block,ny_block,max_blocks))
      allocate(darrayj1   (nx_block,ny_block,max_blocks))
      allocate(darrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(darrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(darrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(darrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(rarrayi1   (nx_block,ny_block,max_blocks))
      allocate(rarrayj1   (nx_block,ny_block,max_blocks))
      allocate(rarrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(rarrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(rarrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(rarrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(iarrayi1   (nx_block,ny_block,max_blocks))
      allocate(iarrayj1   (nx_block,ny_block,max_blocks))
      allocate(iarrayi2   (nx_block,ny_block,nz1,max_blocks))
      allocate(iarrayj2   (nx_block,ny_block,nz1,max_blocks))
      allocate(iarrayi3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(iarrayj3   (nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(larrayi1   (nx_block,ny_block,max_blocks))
      allocate(larrayj1   (nx_block,ny_block,max_blocks))
      allocate(darrayi1str(nx_block,ny_block,max_blocks))
      allocate(darrayj1str(nx_block,ny_block,max_blocks))
      allocate(darrayi10  (nx_block,ny_block,max_blocks))
      allocate(darrayj10  (nx_block,ny_block,max_blocks))

      allocate(cidata_bas(nx_block,ny_block,nz1,nz2,max_blocks))
      allocate(cjdata_bas(nx_block,ny_block,nz1,nz2,max_blocks))

      darrayi1 = fillval
      darrayj1 = fillval
      darrayi2 = fillval
      darrayj2 = fillval
      darrayi3 = fillval
      darrayj3 = fillval
      rarrayi1 = fillval
      rarrayj1 = fillval
      rarrayi2 = fillval
      rarrayj2 = fillval
      rarrayi3 = fillval
      rarrayj3 = fillval
      iarrayi1 = fillval
      iarrayj1 = fillval
      iarrayi2 = fillval
      iarrayj2 = fillval
      iarrayi3 = fillval
      iarrayj3 = fillval
      larrayi1 = .false.
      larrayj1 = .true.
      darrayi1str = fillval
      darrayj1str = fillval
      darrayi10  = fillval
      darrayj10  = fillval
      cidata_bas = fillval
      cjdata_bas = fillval

      call ice_distributionGet(distrb_info, numLocalBlocks = numBlocks)

      !--- baseline data ---
      ! set to the global index
      ! i/j valid everywhere for "cyclic"
      ! i/j valid for "open" with extrapolation on outer boundary
      ! i/j zero on outer boundary for "closed"

      do iblock = 1,numBlocks
         call ice_distributionGetBlockID(distrb_info, iblock, blockID)
         this_block = get_block(blockID, blockID)
         do k2 = 1,nz2
         do k1 = 1,nz1
         do j = 1,ny_block
         do i = 1,nx_block
            cidata_bas(i,j,k1,k2,iblock) = real(this_block%i_glob(i),kind=dbl_kind) + &
                   real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind
            cjdata_bas(i,j,k1,k2,iblock) = real(this_block%j_glob(j),kind=dbl_kind) + &
                   real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind
         enddo
         enddo
         enddo
         enddo
      enddo

      !---------------------------------------------------------------

      testcnt = 0
      do nn = 1, maxtests
      do nl = 1, maxlocs
      do nt = 1, maxtypes
      do nf = 1, maxfills

         !--- setup test ---
         first_call = .true.
         testcnt = testcnt + 1
         if (nf == 1) then
            halofill = .true.
            fillexpected = dhalofillval
         elseif (nf == 2) then
            halofill = .false.
            fillexpected = fillval
         else
            write(6,*) subname,' nf = ',nf
            if (my_task == master_task) then
               write(6,*) ' '
               write(6,*) 'GATHSCATCHK FAILED'
               write(6,*) ' '
            endif
            call abort_ice(subname//' invalid value of nf',file=__FILE__,line=__LINE__)
         endif
         if (testcnt > tottest) then
            if (my_task == master_task) then
               write(6,*) ' '
               write(6,*) 'GATHSCATCHK FAILED'
               write(6,*) ' '
            endif
            call abort_ice(subname//' testcnt > tottest',file=__FILE__,line=__LINE__)
         endif

         !--- fill arrays ---
         darrayi1(:,:,:) = fillval
         darrayj1(:,:,:) = fillval
         darrayi2(:,:,:,:) = fillval
         darrayj2(:,:,:,:) = fillval
         darrayi3(:,:,:,:,:) = fillval
         darrayj3(:,:,:,:,:) = fillval
         darrayi1str(:,:,:) = fillval
         darrayj1str(:,:,:) = fillval
         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            do j = jb,je
               do i = ib,ie
                  k1 = 1
                  k2 = 1
                  darrayi1(i,j,iblock) = cidata_bas(i,j,k1,k2,iblock)
                  darrayj1(i,j,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                  do k1 = 1,nz1
                     k2 = 1
                     darrayi2(i,j,k1,iblock) = cidata_bas(i,j,k1,k2,iblock)
                     darrayj2(i,j,k1,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                     do k2 = 1,nz2
                        darrayi3(i,j,k1,k2,iblock) = cidata_bas(i,j,k1,k2,iblock)
                        darrayj3(i,j,k1,k2,iblock) = cjdata_bas(i,j,k1,k2,iblock)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ! copy original darray1 for "stress" compare
         darrayi10 = darrayi1
         darrayj10 = darrayj1

         !--- halo update ---

         if (nn == 1) then
            k1m = 1
            k2m = 1
            halofld = '2DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 2) then
            k1m = nz1
            k2m = 1
            halofld = '3DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 3) then
            k1m = nz1
            k2m = nz2
            halofld = '4DR8'
            if (halofill) then
               call ice_haloUpdate(darrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate(darrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate(darrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(darrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
         elseif (nn == 4) then
            k1m = 1
            k2m = 1
            halofld = '2DR4'
            rarrayi1 = real(darrayi1,kind=real_kind)
            rarrayj1 = real(darrayj1,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = real(rarrayi1,kind=dbl_kind)
            darrayj1 = real(rarrayj1,kind=dbl_kind)
         elseif (nn == 5) then
            k1m = nz1
            k2m = 1
            halofld = '3DR4'
            rarrayi2 = real(darrayi2,kind=real_kind)
            rarrayj2 = real(darrayj2,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi2 = real(rarrayi2,kind=dbl_kind)
            darrayj2 = real(rarrayj2,kind=dbl_kind)
         elseif (nn == 6) then
            k1m = nz1
            k2m = nz2
            halofld = '4DR4'
            rarrayi3 = real(darrayi3,kind=real_kind)
            rarrayj3 = real(darrayj3,kind=real_kind)
            if (halofill) then
               call ice_haloUpdate(rarrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
               call ice_haloUpdate(rarrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=rhalofillval)
            else
               call ice_haloUpdate(rarrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(rarrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi3 = real(rarrayi3,kind=dbl_kind)
            darrayj3 = real(rarrayj3,kind=dbl_kind)
         elseif (nn == 7) then
            k1m = 1
            k2m = 1
            halofld = '2DI4'
            iarrayi1 = nint(darrayi1)
            iarrayj1 = nint(darrayj1)
            if (halofill) then
               call ice_haloUpdate(iarrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = real(iarrayi1,kind=dbl_kind)
            darrayj1 = real(iarrayj1,kind=dbl_kind)
         elseif (nn == 8) then
            k1m = nz1
            k2m = 1
            halofld = '3DI4'
            iarrayi2 = nint(darrayi2)
            iarrayj2 = nint(darrayj2)
            if (halofill) then
               call ice_haloUpdate(iarrayi2, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj2, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi2, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj2, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi2 = real(iarrayi2,kind=dbl_kind)
            darrayj2 = real(iarrayj2,kind=dbl_kind)
         elseif (nn == 9) then
            k1m = nz1
            k2m = nz2
            halofld = '4DI4'
            iarrayi3 = nint(darrayi3)
            iarrayj3 = nint(darrayj3)
            if (halofill) then
               call ice_haloUpdate(iarrayi3, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
               call ice_haloUpdate(iarrayj3, halo_info, field_loc(nl), field_type(nt), fillvalue=ihalofillval)
            else
               call ice_haloUpdate(iarrayi3, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(iarrayj3, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi3 = real(iarrayi3,kind=dbl_kind)
            darrayj3 = real(iarrayj3,kind=dbl_kind)
         elseif (nn == 10) then
            k1m = 1
            k2m = 1
            halofld = '2DL1'
            where (darrayi1 == fillval)
               larrayi1 = .false.
            elsewhere
               larrayi1 = (mod(nint(darrayi1),2) == 1)
            endwhere
            where (darrayj1 == fillval)
               larrayj1 = .true.
            elsewhere
               larrayj1 = (mod(nint(darrayj1),2) == 1)
            endwhere
            if (halofill) then
               call ice_haloUpdate(larrayi1, halo_info, field_loc(nl), field_type(nt), fillvalue=.false.)
               call ice_haloUpdate(larrayj1, halo_info, field_loc(nl), field_type(nt), fillvalue=.true.)
            else
               call ice_haloUpdate(larrayi1, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate(larrayj1, halo_info, field_loc(nl), field_type(nt))
            endif
            darrayi1 = c0
            where (larrayi1) darrayi1 = c1
            darrayj1 = c0
            where (larrayj1) darrayj1 = c1
         elseif (nn == 11) then
            k1m = 1
            k2m = 1
            halofld = 'STRESS'
            darrayi1str = -darrayi1  ! flip sign for testing
            darrayj1str = -darrayj1
            if (halofill) then
               call ice_haloUpdate_stress(darrayi1, darrayi1str, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
               call ice_haloUpdate_stress(darrayj1, darrayj1str, halo_info, field_loc(nl), field_type(nt), fillvalue=dhalofillval)
            else
               call ice_haloUpdate_stress(darrayi1, darrayi1str, halo_info, field_loc(nl), field_type(nt))
               call ice_haloUpdate_stress(darrayj1, darrayj1str, halo_info, field_loc(nl), field_type(nt))
            endif
         endif

         write(teststring(testcnt),'(4a8,2a12)') trim(halofld),trim(locs_name(nl)),trim(types_name(nt)),trim(fill_name(nf)), &
                      trim(ew_boundary_type),trim(ns_boundary_type)

         do iblock = 1,numBlocks
            call ice_distributionGetBlockID(distrb_info, iblock, blockID)
            this_block = get_block(blockID, blockID)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            ! just check non-padded gridcells
            do j = jb-nghost, je+nghost
            do i = ib-nghost, ie+nghost
            do k1 = 1,k1m
            do k2 = 1,k2m
               tripole_average = .false.
               tripole_pole = .false.
               if (index(halofld,'2D') > 0) then
                  aichk = darrayi1(i,j,iblock)
                  ajchk = darrayj1(i,j,iblock)
               elseif (index(halofld,'STRESS') > 0) then
                  aichk = darrayi1(i,j,iblock)
                  ajchk = darrayj1(i,j,iblock)
               elseif (index(halofld,'3D') > 0) then
                  aichk = darrayi2(i,j,k1,iblock)
                  ajchk = darrayj2(i,j,k1,iblock)
               elseif (index(halofld,'4D') > 0) then
                  aichk = darrayi3(i,j,k1,k2,iblock)
                  ajchk = darrayj3(i,j,k1,k2,iblock)
               else
                  if (my_task == master_task) then
                     write(6,*) ' '
                     write(6,*) 'GATHSCATCHK FAILED'
                     write(6,*) ' '
                  endif
                  call abort_ice(subname//' halofld not matched '//trim(halofld),file=__FILE__,line=__LINE__)
               endif

               cichk = cidata_bas(i,j,k1,k2,iblock)
               cjchk = cjdata_bas(i,j,k1,k2,iblock)

               ! halo special cases

               if (field_loc (nl) == field_loc_noupdate .or. &
                   field_type(nt) == field_type_noupdate) then
                  if (i < ib .or. j < jb .or. i > ie .or. j > je) then
                     ! no halo update anywhere, doesn't even see fillvalue passed in
                     cichk = fillval
                     cjchk = fillval
                  endif

               else

                  ! if boundary_type is not cyclic set outer boundary to fill, other special cases below
                  if (ew_boundary_type /= 'cyclic' .and. &
                      ((this_block%i_glob(ib) == 1         .and. i < ib) .or. &  ! west outer face
                       (this_block%i_glob(ie) == nx_global .and. i > ie))) then  ! east outer face
                     cichk = fillexpected
                     cjchk = fillexpected
                  endif

                  ! if boundary_type is not cyclic set outer boundary to fill, other special cases below
                  ! - tripole north edge will be haloed and is updated below, default to fill value for now
                  ! - tripole south edge will be set to the fillvalue or to haloupdate internal default (c0)
                  !   tripole basically assumes south edge is land or always ice free in CICE
                  if (ns_boundary_type /= 'cyclic' .and. &
                      ((this_block%j_glob(jb) == 1         .and. j < jb) .or. &  ! south outer face
                       (this_block%j_glob(je) == ny_global .and. j > je))) then  ! north outer face
                     if ((ns_boundary_type == 'tripole' .or. &
                          ns_boundary_type == 'tripoleT') .and. &
                         .not. halofill) then
                        cichk = c0
                        cjchk = c0
                     else
                        cichk = fillexpected
                        cjchk = fillexpected
                     endif
                  endif

                  ! zero_gradient and linear_extrap edges then corners
                  if (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap') then
                     wgt1 = c1  ! zero_gradient
                     wgt2 = c0
                     if (this_block%i_glob(ib) == 1         .and. i < ib) then  ! West
                        if (ew_boundary_type == 'linear_extrap') then
                           wgt1 = real(ib-i+1,dbl_kind)
                           wgt2 = real(ib-i  ,dbl_kind)   ! wgt1 - 1
                        endif
                        cichk = wgt1*cidata_bas(ib,j,k1,k2,iblock) - wgt2*cidata_bas(ib+1,j,k1,k2,iblock)
                        cjchk = wgt1*cjdata_bas(ib,j,k1,k2,iblock) - wgt2*cjdata_bas(ib+1,j,k1,k2,iblock)
                     elseif (this_block%i_glob(ie) == nx_global .and. i > ie) then  ! East
                        if (ew_boundary_type == 'linear_extrap') then
                           wgt1 = real(i-ie+1,dbl_kind)
                           wgt2 = real(i-ie  ,dbl_kind)   ! wgt1 - 1
                        endif
                        cichk = wgt1*cidata_bas(ie,j,k1,k2,iblock) - wgt2*cidata_bas(ie-1,j,k1,k2,iblock)
                        cjchk = wgt1*cjdata_bas(ie,j,k1,k2,iblock) - wgt2*cjdata_bas(ie-1,j,k1,k2,iblock)
                     endif
                  endif

                  if (ns_boundary_type == 'zero_gradient' .or. ns_boundary_type == 'linear_extrap') then
                     wgt1 = c1   ! zero_gradient
                     wgt2 = c0
                     if (this_block%j_glob(jb) == 1         .and. j < jb) then  ! South
                        if (ns_boundary_type == 'linear_extrap') then
                           wgt1 = real(jb-j+1,dbl_kind)
                           wgt2 = real(jb-j  ,dbl_kind)   ! wgt1 - 1
                        endif
                        cichk = wgt1*cidata_bas(i,jb,k1,k2,iblock) - wgt2*cidata_bas(i,jb+1,k1,k2,iblock)
                        cjchk = wgt1*cjdata_bas(i,jb,k1,k2,iblock) - wgt2*cjdata_bas(i,jb+1,k1,k2,iblock)
                     elseif (this_block%j_glob(je) == ny_global .and. j > je) then  ! North
                        if (ns_boundary_type == 'linear_extrap') then
                           wgt1 = real(j-je+1,dbl_kind)
                           wgt2 = real(j-je  ,dbl_kind)   ! wgt1 - 1
                        endif
                        cichk = wgt1*cidata_bas(i,je,k1,k2,iblock) - wgt2*cidata_bas(i,je-1,k1,k2,iblock)
                        cjchk = wgt1*cjdata_bas(i,je,k1,k2,iblock) - wgt2*cjdata_bas(i,je-1,k1,k2,iblock)
                     endif

                     ! Boundary Corners, can come at it either direction, do ns then ew
                     if (ew_boundary_type == 'zero_gradient' .or. ew_boundary_type == 'linear_extrap') then
                        wgt1 = c1   ! zero_gradient
                        wgt2 = c0
                        found = .false.
                        if (this_block%i_glob(ib) == 1         .and. i < ib) then
                           if (this_block%j_glob(jb) == 1         .and. j < jb) then
                              found = .true.  ! Southwest
                              if (ns_boundary_type == 'linear_extrap') then
                                 wgt1 = real(jb-j+1,dbl_kind)
                                 wgt2 = real(jb-j  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk1 = wgt1*cidata_bas(ib  ,jb,k1,k2,iblock) - wgt2*cidata_bas(ib  ,jb+1,k1,k2,iblock)
                              cichk2 = wgt1*cidata_bas(ib+1,jb,k1,k2,iblock) - wgt2*cidata_bas(ib+1,jb+1,k1,k2,iblock)
                              cjchk1 = wgt1*cjdata_bas(ib  ,jb,k1,k2,iblock) - wgt2*cjdata_bas(ib  ,jb+1,k1,k2,iblock)
                              cjchk2 = wgt1*cjdata_bas(ib+1,jb,k1,k2,iblock) - wgt2*cjdata_bas(ib+1,jb+1,k1,k2,iblock)
                           elseif (this_block%j_glob(je) == ny_global .and. j > je) then
                              found = .true.  ! Northwest
                              if (ns_boundary_type == 'linear_extrap') then
                                 wgt1 = real(j-je+1,dbl_kind)
                                 wgt2 = real(j-je  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk1 = wgt1*cidata_bas(ib  ,je,k1,k2,iblock) - wgt2*cidata_bas(ib  ,je-1,k1,k2,iblock)
                              cichk2 = wgt1*cidata_bas(ib+1,je,k1,k2,iblock) - wgt2*cidata_bas(ib+1,je-1,k1,k2,iblock)
                              cjchk1 = wgt1*cjdata_bas(ib  ,je,k1,k2,iblock) - wgt2*cjdata_bas(ib  ,je-1,k1,k2,iblock)
                              cjchk2 = wgt1*cjdata_bas(ib+1,je,k1,k2,iblock) - wgt2*cjdata_bas(ib+1,je-1,k1,k2,iblock)
                           endif
                           if (found) then
                              wgt1 = c1   ! zero_gradient
                              wgt2 = c0
                              if (ew_boundary_type == 'linear_extrap') then
                                 wgt1 = real(ib-i+1,dbl_kind)
                                 wgt2 = real(ib-i  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk = wgt1*cichk1 - wgt2*cichk2
                              cjchk = wgt1*cjchk1 - wgt2*cjchk2
                           endif
                        elseif (this_block%i_glob(ie) == nx_global .and. i > ie) then
                           if (this_block%j_glob(jb) == 1         .and. j < jb) then
                              found = .true.  ! Southeast
                              if (ns_boundary_type == 'linear_extrap') then
                                 wgt1 = real(jb-j+1,dbl_kind)
                                 wgt2 = real(jb-j  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk1 = wgt1*cidata_bas(ie  ,jb,k1,k2,iblock) - wgt2*cidata_bas(ie  ,jb+1,k1,k2,iblock)
                              cichk2 = wgt1*cidata_bas(ie-1,jb,k1,k2,iblock) - wgt2*cidata_bas(ie-1,jb+1,k1,k2,iblock)
                              cjchk1 = wgt1*cjdata_bas(ie  ,jb,k1,k2,iblock) - wgt2*cjdata_bas(ie  ,jb+1,k1,k2,iblock)
                              cjchk2 = wgt1*cjdata_bas(ie-1,jb,k1,k2,iblock) - wgt2*cjdata_bas(ie-1,jb+1,k1,k2,iblock)
                           elseif (this_block%j_glob(je) == ny_global .and. j > je) then
                              found = .true.  ! Northeast
                              if (ns_boundary_type == 'linear_extrap') then
                                 wgt1 = real(j-je+1,dbl_kind)
                                 wgt2 = real(j-je  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk1 = wgt1*cidata_bas(ie  ,je,k1,k2,iblock) - wgt2*cidata_bas(ie  ,je-1,k1,k2,iblock)
                              cichk2 = wgt1*cidata_bas(ie-1,je,k1,k2,iblock) - wgt2*cidata_bas(ie-1,je-1,k1,k2,iblock)
                              cjchk1 = wgt1*cjdata_bas(ie  ,je,k1,k2,iblock) - wgt2*cjdata_bas(ie  ,je-1,k1,k2,iblock)
                              cjchk2 = wgt1*cjdata_bas(ie-1,je,k1,k2,iblock) - wgt2*cjdata_bas(ie-1,je-1,k1,k2,iblock)
                           endif
                           if (found) then
                              wgt1 = c1   ! zero_gradient
                              wgt2 = c0
                              if (ew_boundary_type == 'linear_extrap') then
                                 wgt1 = real(i-ie+1,dbl_kind)
                                 wgt2 = real(i-ie  ,dbl_kind)   ! wgt1 - 1
                              endif
                              cichk = wgt1*cichk1 - wgt2*cichk2
                              cjchk = wgt1*cjchk1 - wgt2*cjchk2
                           endif
                        endif
                     endif
                  endif  ! zero_gradient, linear_extrap

                  if (index(halofld,'STRESS') > 0) then
                     ! only updates on tripole zipper for tripole grids
                     ! darrayi10 is copy of darrayi1 before halo call
                     cichk = darrayi10(i,j,iblock)
                     cjchk = darrayj10(i,j,iblock)
                  endif

                  !--- tripole on north boundary, need to hardcode ---
                  !--- tripole and tripoleT slightly different     ---
                  !--- establish special set of points here        ---
                  if ((this_block%j_glob(je) == ny_global) .and. &
                     ((ns_boundary_type == 'tripole'  .and. &
                       (j > je .or. &
                        (j == je .and. (field_loc(nl) == field_loc_Nface .or. field_loc(nl) == field_loc_NEcorner)))) .or. &
                      (ns_boundary_type == 'tripoleT' .and. &
                       (j >= je)))) then

                     ! flip sign for vector/angle except for logical halo updates
                     rsign = c1
                     if ((field_type(nt) == field_type_vector .or. field_type(nt) == field_type_angle) .and. &
                         .not. (index(halofld,'L1') > 0)) then
                        rsign = -c1
                     endif

                     ! for tripole
                     if (ns_boundary_type == 'tripole') then

                        ! compute itrip and jtrip, these are the location where the halo values are defined for i,j
                        ! for j=je averaging, itrip and jtrip are the 2nd gridpoint associated with averaging

                        ! standard center tripole u-fold
                        itrip = nx_global-this_block%i_glob(i)+1
                        jtrip = max(je - (j-je) + 1 , je)
                        ioffset = 0
                        joffset = 0

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Nface) then
                           ! need j offset
                           joffset = -1
                           if (j == je) then
                              tripole_average = .true.
                           endif
                        endif

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Eface) then
                           ! fold plus cell offset
                           ioffset = -1
                           ! CICE treats j=ny_global tripole edge points incorrectly
                           ! should do edge wraparound and average
                           ! CICE does not update those points, assumes it's "land"
                           if (j == je) then
                              if (this_block%i_glob(i) == nx_global/2) then
                                  tripole_pole = .true.
                              elseif (this_block%i_glob(i) == nx_global  ) then
                                  tripole_pole = .true.
                              endif
                           endif
                        endif

                     ! for tripoleT
                     elseif (ns_boundary_type == 'tripoleT') then

                        ! compute itrip and jtrip, these are the location where the halo values are defined for i,j
                        ! for j=je averaging, itrip and jtrip are the 2nd gridpoint associated with averaging

                        ! standard center tripoleT t-fold
                        itrip = nx_global-this_block%i_glob(i)+2
                        jtrip = je - (j-je)
                        ioffset = 0
                        joffset = 0

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Eface) then
                           ! fold plus cell offset
                           ioffset = -1
                        endif

                        if (field_loc(nl) == field_loc_NEcorner .or. field_loc(nl) == field_loc_Nface) then
                           ! need j offset
                           joffset = -1
                        endif

                        if (field_loc(nl) == field_loc_Center .or. field_loc(nl) == field_loc_Eface) then
                           if (j == je) then
                              tripole_average = .true.
                           endif
                        endif

                        ! center point poles need to be treated special
                        if (field_loc(nl) == field_loc_Center) then
                           if (j == je .and. &
                              (this_block%i_glob(i) == 1 .or. this_block%i_glob(i) == nx_global/2+1)) then
                              tripole_pole = .true.
                           endif
                        endif

                     endif

                     itrip = mod(itrip + ioffset + nx_global-1,nx_global)+1
                     jtrip = jtrip + joffset

                     rival = (real(itrip,kind=dbl_kind) + &
                              real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind)
                     rjval = (real(this_block%j_glob(jtrip),kind=dbl_kind) + &
                              real(k1,kind=dbl_kind)*1000._dbl_kind + real(k2,kind=dbl_kind)*10000._dbl_kind)

                     if (index(halofld,'STRESS') > 0) then
                        ! only updates on tripole zipper for tripole grids, not tripoleT
                        ! note: L1 and STRESS never overlap so don't worry about L1 here
                        if (tripole_pole) then
                           ! flip sign due to sign of darrayi1str
                           ! ends of tripole seam not averaged in CICE
                           cichk = -rsign * cidata_bas(i,j,k1,k2,iblock)
                           cjchk = -rsign * cjdata_bas(i,j,k1,k2,iblock)
                        else
                           cichk = -rsign * rival
                           cjchk = -rsign * rjval
                        endif

                     elseif (tripole_pole) then
                        ! ends of tripole seam not averaged in CICE
                        cichk = rsign * cidata_bas(i,j,k1,k2,iblock)
                        cjchk = rsign * cjdata_bas(i,j,k1,k2,iblock)

                     elseif (tripole_average) then
                        if (index(halofld,'L1') > 0) then
                           ! logical math doesn't work this way, force to correct answer
                           cichk = aichk ! p5 * (mod(nint(cidata_bas(i,j,k1,k2,iblock)),2) + rsign * mod(nint(rival),2))
                           cjchk = ajchk ! p5 * (mod(nint(cidata_bas(i,j,k1,k2,iblock)),2) + rsign * mod(nint(rjval),2))
                        else
                           cichk = p5 * (cidata_bas(i,j,k1,k2,iblock) + rsign * rival)
                           cjchk = p5 * (cjdata_bas(i,j,k1,k2,iblock) + rsign * rjval)
                        endif

                     else
                        ! standard tripole fold
                        cichk = rsign * rival
                        cjchk = rsign * rjval
                     endif

                  endif  ! tripole or tripoleT

               endif

               if (index(halofld,'I4') > 0) then
                  cichk = real(nint(cichk),kind=dbl_kind)
                  cjchk = real(nint(cjchk),kind=dbl_kind)
               endif

               if (index(halofld,'L1') > 0) then
                  if (cichk == dhalofillval .or. cichk == fillval) then
                     cichk = c0
                  else
                     cichk = mod(nint(cichk),2)
                  endif
                  if (cjchk == dhalofillval .or. cjchk == fillval) then
                     cjchk = c1
                  else
                     cjchk = mod(nint(cjchk),2)
                  endif
               endif

               ptcnt(testcnt) = ptcnt(testcnt) + 1
               call chkresults(aichk,cichk,errorflag(testcnt),testcnt,failcnt(testcnt), &
                    i,j,k1,k2,iblock,first_call,teststring(testcnt),trim(halofld)//'_I')
               call chkresults(ajchk,cjchk,errorflag(testcnt),testcnt,failcnt(testcnt), &
                    i,j,k1,k2,iblock,first_call,teststring(testcnt),trim(halofld)//'_J')
            enddo  ! k2
            enddo  ! k1
            enddo  ! i
            enddo  ! j
         enddo  ! iblock

      enddo  ! maxfills
      enddo  ! maxtypes
      enddo  ! maxlocs
      enddo  ! maxtests

      ! ---------------------------
      ! SUMMARY
      ! ---------------------------

      do n = 1,tottest
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
         ptcntsum = global_sum(ptcnt(n),distrb_info)
         ptcnt(n) = ptcntsum
         failcntsum = global_sum(failcnt(n),distrb_info)
         failcnt(n) = failcntsum
      enddo
      errorflag0 = maxval(errorflag(:))

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GATHSCATCHK COMPLETED SUCCESSFULLY'
         write(6,*) ' '
         tpcnt = 0
         tfcnt = 0
         do n = 1,tottest
            if (errorflag(n) == passflag) then
               tpcnt = tpcnt + 1
               write(6,'(2a,2i9)') 'PASS ',trim(teststring(n)),ptcnt(n),failcnt(n)
            else
               tfcnt = tfcnt + 1
               write(6,'(2a,2i9)') 'FAIL ',trim(teststring(n)),ptcnt(n),failcnt(n)
            endif
         enddo
         write(6,*) ' '
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

      subroutine chkresults(a1,r1,errorflag,testcnt,failcnt,i,j,k1,k2,iblock,first_call,teststring,halofld)

      use gathscatchk_data

      implicit none

      real(dbl_kind)   , intent(in)    :: a1,r1
      integer(int_kind), intent(inout) :: errorflag, failcnt
      integer(int_kind), intent(in)    :: i,j,k1,k2,iblock,testcnt
      logical          , intent(inout) :: first_call
      character(len=*) , intent(in)    :: teststring,halofld

      logical,parameter :: print_always = .false.
      character(len=*) , parameter :: subname='(chkresults)'

      if (a1 /= r1 .or. print_always) then
         if (a1 /= r1) then
            errorflag = failflag
            failcnt = failcnt + 1
         endif
         if (first_call) then
            write(100+my_task,*) ' '
            write(100+my_task,'(a,i4,2a)') '------- TEST = ',testcnt,' ',trim(teststring)
            write(100+my_task,*) ' '
            write(100+my_task,'(a)') '           test  task    i     j    k1    k2  iblock  expected   halocomp       diff'
            first_call = .false.
         endif
         write(100+my_task,1001) trim(halofld),testcnt,my_task,i,j,k1,k2,iblock,r1,a1,r1-a1
      endif

 1001 format(a8,7i6,3f12.3)

      end subroutine chkresults
!=======================================================================
#endif

