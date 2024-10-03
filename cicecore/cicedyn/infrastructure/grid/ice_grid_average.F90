module ice_grid_average

    use ice_kinds_mod
    use ice_constants, only: c0, p5, p25
    use ice_domain, only: nblocks, blocks_ice
    use ice_blocks, only: block, get_block
    use ice_exit, only: abort_ice

    implicit none
    private
    public :: grid_average_X2YS, grid_average_X2YA, grid_average_X2YF, grid_average_X2Y_2


    contains
  !=======================================================================
  ! Shifts quantities from one grid to another
  ! State masked version, simple area weighted averager
  ! NOTE: Input array includes ghost cells that must be updated before
  !       calling this routine.
  !
  ! author: T. Craig
  
    subroutine grid_average_X2YS(dir,work1,wght1,mask1,work2)
  
        character(len=*) , intent(in) :: &
           dir
  
        real (kind=dbl_kind), intent(in) :: &
           work1(:,:,:), &
           wght1(:,:,:), &
           mask1(:,:,:)
  
        real (kind=dbl_kind), intent(out) :: &
           work2(:,:,:)
  
        ! local variables
  
        integer (kind=int_kind) :: &
           i, j, iblk, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        real (kind=dbl_kind) :: &
           wtmp
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(grid_average_X2YS)'
  
        work2(:,:,:) = c0
  
        select case (trim(dir))
  
           case('NE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                          + mask1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                          + mask1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                          + mask1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + mask1(i+1,j  ,iblk)*work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                     + mask1(i  ,j+1,iblk)*work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                     + mask1(i+1,j+1,iblk)*work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                          + mask1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                          + mask1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                          + mask1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + mask1(i-1,j  ,iblk)*work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                     + mask1(i  ,j-1,iblk)*work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                     + mask1(i-1,j-1,iblk)*work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                          + mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                          + mask1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                          + mask1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i-1,j  ,iblk)*work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                     + mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + mask1(i-1,j+1,iblk)*work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                     + mask1(i  ,j+1,iblk)*work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                          + mask1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                          + mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                          + mask1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i  ,j-1,iblk)*work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                     + mask1(i+1,j-1,iblk)*work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                     + mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + mask1(i+1,j  ,iblk)*work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('E')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                          + mask1(i+1,j,iblk)*wght1(i+1,j,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i  ,j,iblk)*work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                     + mask1(i+1,j,iblk)*work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('W')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                          + mask1(i  ,j,iblk)*wght1(i  ,j,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i-1,j,iblk)*work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                     + mask1(i  ,j,iblk)*work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('N')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                          + mask1(i,j+1,iblk)*wght1(i,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i,j  ,iblk)*work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                     + mask1(i,j+1,iblk)*work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('S')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                          + mask1(i,j  ,iblk)*wght1(i,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1(i,j-1,iblk)*work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                     + mask1(i,j  ,iblk)*work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case default
              call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
           end select
  
    end subroutine grid_average_X2YS

!=======================================================================
! Shifts quantities from one grid to another
! State unmasked version, simple weighted averager
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

    subroutine grid_average_X2YA(dir,work1,wght1,work2)

        character(len=*) , intent(in) :: &
           dir
  
        real (kind=dbl_kind), intent(in) :: &
           work1(:,:,:), &
           wght1(:,:,:)
  
        real (kind=dbl_kind), intent(out) :: &
           work2(:,:,:)
  
        ! local variables
  
        integer (kind=int_kind) :: &
           i, j, iblk, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        real (kind=dbl_kind) :: &
           wtmp
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(grid_average_X2YA)'
  
        work2(:,:,:) = c0
  
        select case (trim(dir))
  
           case('NE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i  ,j  ,iblk)  &
                          + wght1(i+1,j  ,iblk)  &
                          + wght1(i  ,j+1,iblk)  &
                          + wght1(i+1,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                     + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                     + work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i  ,j  ,iblk)  &
                          + wght1(i-1,j  ,iblk)  &
                          + wght1(i  ,j-1,iblk)  &
                          + wght1(i-1,j-1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                     + work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                     + work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i-1,j  ,iblk)  &
                          + wght1(i  ,j  ,iblk)  &
                          + wght1(i-1,j+1,iblk)  &
                          + wght1(i  ,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                     + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                     + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i  ,j-1,iblk)  &
                          + wght1(i+1,j-1,iblk)  &
                          + wght1(i  ,j  ,iblk)  &
                          + wght1(i+1,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                     + work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                     + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('E')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i  ,j,iblk)  &
                          + wght1(i+1,j,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                     + work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('W')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i-1,j,iblk)  &
                          + wght1(i  ,j,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                     + work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('N')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i,j  ,iblk)  &
                          + wght1(i,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                     + work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('S')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1(i,j-1,iblk)  &
                          + wght1(i,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                     + work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case default
              call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
           end select
  
    end subroutine grid_average_X2YA
  

!=======================================================================
! Shifts quantities from one grid to another
! Flux masked, original implementation based on earlier t2u and u2t versions
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

    subroutine grid_average_X2YF(dir,work1,wght1,work2,wght2)

        character(len=*) , intent(in) :: &
           dir
  
        real (kind=dbl_kind), intent(in) :: &
           work1(:,:,:), &
           wght1(:,:,:), &
           wght2(:,:,:)
  
        real (kind=dbl_kind), intent(out) :: &
           work2(:,:,:)
  
        ! local variables
  
        integer (kind=int_kind) :: &
           i, j, iblk, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(grid_average_X2YF)'
  
        work2(:,:,:) = c0
  
        select case (trim(dir))
  
           case('NE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p25 * &
                                      (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                     + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                     + work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                     / wght2(i  ,j  ,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p25 *  &
                                     (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                    + work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                    + work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                    + work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                    / wght2(i  ,j  ,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NW')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p25 * &
                                      (work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                     + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                     + work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                     + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                     / wght2(i  ,j  ,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('SE')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p25 *  &
                                     (work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                    + work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                    + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                    + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                    / wght2(i  ,j  ,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('E')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p5 * &
                                      (work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                     + work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                     / wght2(i  ,j,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('W')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p5 * &
                                      (work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                     + work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                     / wght2(i  ,j,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('N')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p5 * &
                                      (work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                     + work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                     / wght2(i  ,j,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('S')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    work2(i,j,iblk) = p5 * &
                                      (work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                     + work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                     / wght2(i  ,j,iblk)
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case default
              call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
        
        end select
  
    end subroutine grid_average_X2YF

!=======================================================================
! Shifts quantities from one grid to another
! State masked version, simple weighted averager
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

    subroutine grid_average_X2Y_2(dir,work1a,wght1a,mask1a,work1b,wght1b,mask1b,work2)

        character(len=*) , intent(in) :: &
           dir
  
        real (kind=dbl_kind), intent(in) :: &
           work1a(:,:,:), work1b(:,:,:), &
           wght1a(:,:,:), wght1b(:,:,:), &
           mask1a(:,:,:), mask1b(:,:,:)
  
        real (kind=dbl_kind), intent(out) :: &
           work2(:,:,:)
  
        ! local variables
  
        integer (kind=int_kind) :: &
           i, j, iblk, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain
  
        real (kind=dbl_kind) :: &
           wtmp
  
        type (block) :: &
           this_block           ! block information for current block
  
        character(len=*), parameter :: subname = '(grid_average_X2Y_2)'
  
        work2(:,:,:) = c0
  
        select case (trim(dir))
  
           case('NE2US')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                          + mask1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                          + mask1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                          + mask1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1a(i  ,j  ,iblk)*work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                     + mask1a(i+1,j  ,iblk)*work1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                                     + mask1b(i  ,j  ,iblk)*work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                                     + mask1b(i  ,j+1,iblk)*work1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NE2TS')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (mask1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                          + mask1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                          + mask1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                          + mask1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (mask1a(i  ,j-1,iblk)*work1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                                     + mask1a(i  ,j  ,iblk)*work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                     + mask1b(i-1,j  ,iblk)*work1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                                     + mask1b(i  ,j  ,iblk)*work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NE2UA')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1a(i  ,j  ,iblk)  &
                          + wght1a(i+1,j  ,iblk)  &
                          + wght1b(i  ,j  ,iblk)  &
                          + wght1b(i  ,j+1,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                     + work1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                                     + work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                                     + work1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case('NE2TA')
              !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
              do iblk = 1, nblocks
                 this_block = get_block(blocks_ice(iblk),iblk)
                 ilo = this_block%ilo
                 ihi = this_block%ihi
                 jlo = this_block%jlo
                 jhi = this_block%jhi
                 do j = jlo, jhi
                 do i = ilo, ihi
                    wtmp = (wght1a(i  ,j-1,iblk)  &
                          + wght1a(i  ,j  ,iblk)  &
                          + wght1b(i-1,j  ,iblk)  &
                          + wght1b(i  ,j  ,iblk))
                    if (wtmp /= c0) &
                    work2(i,j,iblk) = (work1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                                     + work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                     + work1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                                     + work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)) &
                                     / wtmp
                 enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
  
           case default
              call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
           end select
  
    end subroutine grid_average_X2Y_2

end module ice_grid_average