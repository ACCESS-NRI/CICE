
subroutine handle_nc_err(err_code, string)
   use ice_exit, only: abort_ice
   use netCDF, only: nf90_strerror
   integer, intent(in) :: err_code
   character (len=*), intent(in) :: string

   if (my_task == master_task) then
      write(6,*)  string
      write(6,*)  nf90_strerror(err_code)
   endif
   call abort_ice
end subroutine handle_nc_err

module test_grid_var
   use ice_kinds_mod, only: int_kind, dbl_kind
   use ice_communicate, only: my_task, master_task
   use ice_gather_scatter, only: gather_global

   implicit none

   type :: t_grid_var
      real (kind=dbl_kind), dimension (:,:,:), allocatable :: gv
      character(len=32) :: fn
      logical :: optional = .False.
   end type

   type :: t_area_var
      real (kind=dbl_kind), dimension (:,:,:), allocatable :: xlen,ylen,area
      character(len=32) :: name
   end type

   logical :: debug = .False.
   real (kind = dbl_kind), parameter :: PI = 4*atan(1.0_16)

   real (kind=dbl_kind), dimension(:,:), public,  allocatable :: &
      grid_var_glob, test_var, calc_var, eps_var, xlen, ylen
   real (kind=dbl_kind), dimension(:,:,:), public,  allocatable :: grid_var
   character(len=*), public,  parameter ::  &
      passflag = 'PASS', &
      failflag = 'FAIL', &
      skipflag = 'SKIP'
   real (kind=dbl_kind), parameter :: eps = 1.e-6

   contains

   subroutine test_var_size(fid, varid, var, errorflag) 

      ! Get a var from the original netcdf
      ! Confirm it exists and is of the right size 
      ! Set errorflag according to the result

      use netcdf, only: nf90_inquire_variable,nf90_inquire_dimension, NF90_NOERR, nf90_max_var_dims      
      use ice_domain_size, only: nx_global, ny_global

      integer, intent(in) :: fid, varid
      integer :: status, test_nx, test_ny
      integer, dimension(nf90_max_var_dims) :: dimids
      character(len=8), intent(out)  :: errorflag
      type(t_grid_var), intent(in) :: var
      
      status = nf90_inquire_variable(fid, varid, dimids = dimids)
      if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error ")
      status = nf90_inquire_dimension(fid, dimids(1), len = test_nx)
      if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error  ")
      status = nf90_inquire_dimension(fid, dimids(2), len = test_ny)
      if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error ")

      if (debug) write (6,*) "LOG: checking size of "//trim(var%fn)
      if (test_nx/=nx_global .or. test_ny/=ny_global) then
         if (my_task == master_task) then
            errorflag = failflag
            write (6,*) "Error in ", var%fn, " nx_global or ny_global from ice_in does not match grid file"
         endif
      endif

   end subroutine test_var_size

   subroutine test_var_unchanged(fid, varid, var, errorflag) 

      ! Get a var from the original netcdf
      ! Get the var from cice through the "global_gather"
      ! Confirm they are the same 
      ! Set errorflag according to the result

      use netcdf, only: nf90_get_var, NF90_NOERR
      use ice_domain, only: distrb_info

      integer, intent(in) :: fid, varid
      integer :: status
      character(len=8), intent(out)  :: errorflag
      type(t_grid_var), intent(in) :: var
      
      if (debug) write (6,*) "LOG: attempting global gather "//trim(var%fn)
      call gather_global(grid_var_glob, var%gv, master_task, distrb_info)

      if (my_task == master_task) then
         status = nf90_get_var(fid, varid, test_var)
         if (status==NF90_NOERR) then
            if (var%fn=='htn' .or. var%fn=='hte') test_var=test_var/100
            call assert_array_equal(errorflag, var%fn)
         else
            call handle_nc_err(status, \
               "ERR: Could not load variable '"//trim(var%fn)//"' from grid file")
            errorflag = failflag
         endif             
      endif

   end subroutine test_var_unchanged

   subroutine assert_array_equal(errorflag, varname)
      use ice_domain_size, only: nx_global, ny_global

      character(len=8), intent(out)  :: errorflag
      character(len=*), intent(in) :: varname
      integer(kind=int_kind) :: i,j
      

      eps_var = abs(grid_var_glob-test_var)
      if (all(eps_var<eps)) then
         write (6,*) "PASS: '"//trim(varname)//"' in cice matches"
         errorflag = passflag
      else
         write (6,*) "FAIL: '"//trim(varname)//"' in cice does not match"
         write (6,*) "--------CICE Var - Expected (Test) Var--------"
         do j=1,ny_global
            do i=1,nx_global 
               if (eps_var(i,j)>=eps) then
                  write (6,*) "Fail at i = ",i," j = ",j, " err = ", eps_var(i,j)
               endif
            enddo
         enddo
         errorflag = failflag
      endif
   end subroutine assert_array_equal

   subroutine test_interp(errorflag, varname)
      use ice_domain_size, only: nx_global, ny_global
      use ice_domain, only: distrb_info
      use ice_grid
      use ice_exit, only: abort_ice


      character(len=8), intent(out)  :: errorflag
      character(len=*), intent(in) :: varname
      integer(kind=int_kind) :: i,j
      real (kind=dbl_kind) :: angle_0, angle_w, angle_s, angle_sw

      select case (varname)
         case ("angleT")
            call gather_global(calc_var, angle, master_task, distrb_info)
            if (my_task==master_task) then
               do j=2,ny_global
                  do i=2,nx_global
                     angle_0  = calc_var(i  ,j)   !   w----0
                     angle_w  = calc_var(i-1,j)   !   |    |
                     angle_s  = calc_var(i,  j-1) !   |    |
                     angle_sw = calc_var(i-1,j-1) !   sw---s
                     test_var(i,j) = atan2(0.25*(sin(angle_0)+ & 
                                                   sin(angle_w)+ &
                                                   sin(angle_s)+ &
                                                   sin(angle_sw)),&
                                             0.25*(cos(angle_0)+ &
                                                   cos(angle_w)+ &
                                                   cos(angle_s)+ &
                                                   cos(angle_sw)))                           ! test_var(i,j) = 0.5*calc_var(i+1,j) + 0.5*c
                  enddo
               enddo
               ! j=1
               do i=2,nx_global
                  angle_0  = calc_var(i  ,1)
                  angle_w  = calc_var(i-1,1) 
                  angle_s  = calc_var(i,  1) 
                  angle_sw = calc_var(i-1,1) 
                  test_var(i,1) = atan2(0.25*(sin(angle_0)+ & 
                                                sin(angle_w)+ &
                                                sin(angle_s)+ &
                                                sin(angle_sw)),&
                                          0.25*(cos(angle_0)+ &
                                                cos(angle_w)+ &
                                                cos(angle_s)+ &
                                                cos(angle_sw)))                           ! test_var(i,j) = 0.5*calc_var(i+1,j) + 0.5*c
               enddo
               ! i=1
               do j=2,ny_global
                  angle_0  = calc_var(1  ,j)
                  angle_w  = calc_var(nx_global,j)
                  angle_s  = calc_var(1,  j-1) 
                  angle_sw = calc_var(nx_global,j-1) 
                  test_var(1,j) = atan2(0.25*(sin(angle_0)+ & 
                                                sin(angle_w)+ &
                                                sin(angle_s)+ &
                                                sin(angle_sw)),&
                                          0.25*(cos(angle_0)+ &
                                                cos(angle_w)+ &
                                                cos(angle_s)+ &
                                                cos(angle_sw)))                           ! test_var(i,j) = 0.5*calc_var(i+1,j) + 0.5*c
               enddo                     
            endif
         case default
            call abort_ice(varname//" not tested, consider adding")
      end select
      if (my_task == master_task) call assert_array_equal(errorflag, varname)

   end subroutine test_interp


end module test_grid_var


      
program gridloadchk

   ! This tests the CICE grid netcdf routines by setting up cice to the point the grid is loaded
   ! and verifies results from hardwired inputs with known outputs

   ! Test was written with these options:
   ! ./cice.setup -m gadi -e intel --test unittest -s gridloadchk,iohdf5,gx3nc --testid unittest 
   ! ./cice.build gridloadchk

   use ice_kinds_mod, only: int_kind, dbl_kind
   use ice_grid
   use ice_domain, only: distrb_info
   use ice_domain_size, only: nx_global, ny_global
   use ice_communicate, only: my_task, master_task
   use ice_exit, only: abort_ice, end_run
   use ice_timers, only: timer_total, init_ice_timers, ice_timer_start
   use ice_gather_scatter, only: gather_global

   use test_grid_var
   
   use CICE_InitGrid
   use netcdf, only: nf90_open, nf90_inq_varid, NF90_NOERR, NF90_NOWRITE, nf90_get_var
   implicit none

   integer(kind=int_kind) :: n,m,i,j
   integer(kind=int_kind), parameter :: ntests = 25 , ngridvar = 8 , nareavar = 4
   character(len=8)  :: errorflag0,errorflag(1:ntests),errorflagtmp
   ! character(len=32) :: testname(ntests)
   integer :: fid , varid, status , itest
   type(t_grid_var) :: grid_vars(ngridvar)
   type(t_area_var) :: area_vars(nareavar)

   call cice_init

   errorflag0   = passflag
   errorflag(:) = "NOTRUN"
   ! testname(:) = ''
   ! testname(1) = 'compute_elapsed_days'

   if (my_task /= master_task) then
      ! only debug on master task
      debug = .False.
   endif

   call init_grid1 ! domain distribution
   call alloc_grid  
   call init_ice_timers      ! initialize all timers, this is needed by init_grid2
   call init_grid2           ! loads grid variables

#ifndef USE_NETCDF
   call abort_ice("Test only supports netcdf grid files (ERROR: USE_NETCDF cpp not defined)")
#else
   if (trim(grid_format) /= 'nc') then
      call abort_ice("Test only supports netcdf grid files, set grid_format to 'nc' in namslist")
   endif

   if (my_task == master_task) then
      write(6,*) '-------------------------'
      write(6,*) 'RunningUnitTest GRIDLOADCHK'
      write(6,*) '-------------------------'
   endif

   ! --------------------
   ! This section confirms that 
   ! if we gather the grid from each task, this is unchanged from the original grid
   ! ---------------------

   ! These are the fields in popgrid_nc subroutine. except htn/hte - these are tested later
   grid_vars(1)%fn = "ulat"
   grid_vars(1)%gv =  ULAT
   grid_vars(2)%fn = "ulon"
   grid_vars(2)%gv =  ULON
   grid_vars(3)%fn = "angle"
   grid_vars(3)%gv =  ANGLE
   grid_vars(4)%fn = "angleT"
   grid_vars(4)%gv =  ANGLET
   grid_vars(4)%optional = .True.
   grid_vars(5)%fn = "tlon"
   grid_vars(5)%gv =  TLON
   grid_vars(5)%optional = .True.
   grid_vars(6)%fn = "tlat"
   grid_vars(6)%gv =  TLAT
   grid_vars(6)%optional = .True.
   grid_vars(7)%fn = "hte"
   grid_vars(7)%gv =  HTE
   grid_vars(8)%fn = "htn"
   grid_vars(8)%gv =  HTN

   if (my_task == master_task) then
      allocate(test_var(nx_global, ny_global))
      allocate(grid_var_glob(nx_global, ny_global))
      allocate(calc_var(nx_global, ny_global))
      allocate(eps_var(nx_global, ny_global))
   else
      allocate(test_var(1,1))
      allocate(grid_var_glob(1,1))
      allocate(calc_var(1,1))
      allocate(eps_var(1,1))
   endif

   ! loop each expected variable in the grid file      

   if (debug) write (6,*) "LOG: opening test instance of grid file"
   status = nf90_open(grid_file, NF90_NOWRITE, fid)
   if (status/=NF90_NOERR) call handle_nc_err(status, "Could not open "//grid_file//" for testing")

   do itest = 1 , ngridvar

      if (debug) write (6,*) "LOG: finding '"//trim(grid_vars(iTest)%fn)//"' var"
      status = nf90_inq_varid(fid, grid_vars(iTest)%fn, varid) 
      if (status/=NF90_NOERR) then
         if (grid_vars(iTest)%optional) then
            ! confirm CICE has internally calculated good values
            call gather_global(grid_var_glob, grid_vars(iTest)%gv, master_task, distrb_info)
            call test_interp(errorflag(itest), grid_vars(iTest)%fn)
         else
            call handle_nc_err(status, "Could not find variable '"//trim(grid_vars(iTest)%fn)//"' in "//grid_file )
         endif
      else              
         call test_var_size(fid, varid, grid_vars(iTest), errorflag(itest))
         if (errorflag(itest) /= failflag) call test_var_unchanged(fid, varid, grid_vars(iTest), errorflag(itest)) 
      endif
   enddo

   ! --------------------
   ! This section tests the averaging for cell sidelengths 
   ! Source file contains dxN, dxE

   ! dXN should equal htn
   ! dxT should equal htn averaged from j & j-1 cell
   ! dxE should equal htn averaged from j & j-1 of i & i+1 cell
   ! dxU should equal htn averaged from i & i+1 cell
   ! dYE should equal hte
   ! dyT should equal hte averaged from i-1 & i cell
   ! dyN should equal hte averaged from j & j+1 of i-1 & i cell
   ! dyU should equal hte averaged from j & j+1 
   ! ---------------------

   ! dxN
   itest=ngridvar+1
   call gather_global(grid_var_glob, dxN, master_task, distrb_info)
   if (my_task == master_task) then
      status = nf90_inq_varid(fid, "htn", varid) 
      if (status/=NF90_NOERR) call handle_nc_err(status, "Could not find variable htn in grid file")
      status = nf90_get_var(fid, varid, test_var)
      if (status/=NF90_NOERR) call handle_nc_err(status, "Could not load variable htn from grid file")
      test_var = test_var/100 !cm to m
      call assert_array_equal(errorflag(itest), "dxN")
   endif
   calc_var = test_var

   ! dxT
   itest=itest+1
   call gather_global(grid_var_glob, dxT, master_task, distrb_info)
   if (my_task == master_task) then
      do j=2,ny_global
         do i = 1, nx_global
            test_var(i,j)=0.5*(calc_var(i,j-1)+calc_var(i,j))
         enddo
      enddo
      ! j=1
      do i = 1, nx_global
         ! this is a strange interpolation
         ! surely it should be 2*calc_var(i,1)-calc_var(i,2)
         test_var(i,1)=2*calc_var(i,2)-calc_var(i,3)
      enddo
      call assert_array_equal(errorflag(itest), "dxT")
   endif

   ! dxE
   itest=itest+1
   call gather_global(grid_var_glob, dxE, master_task, distrb_info)
   if (my_task == master_task) then
      do j=2,ny_global
         do i = 1, (nx_global-1)
            test_var(i,j)=0.25*(calc_var(i,j-1)+calc_var(i,j)+calc_var(i+1,j-1)+calc_var(i+1,j))
         enddo
      enddo
      !j=1
      do i = 1, nx_global-1
         test_var(i,1)=calc_var(i,2)-0.5*calc_var(i,3)+calc_var(i+1,2)-0.5*calc_var(i+1,3)
      enddo
      !i=nx_global, cycle around
      do j=2,ny_global
         test_var(nx_global,j)=0.25*(calc_var(nx_global,j-1)+calc_var(nx_global,j)+calc_var(1,j-1)+calc_var(1,j))
      enddo
      !i=nx_global, j=1
      test_var(nx_global,1)=calc_var(nx_global,2)-0.5*calc_var(nx_global,3)+calc_var(1,2)-0.5*calc_var(1,3)
      call assert_array_equal(errorflag(itest), "dxE")
   endif

   ! dxU
   itest=itest+1
   call gather_global(grid_var_glob, dxU, master_task, distrb_info)
   if (my_task == master_task) then
      do j=1,ny_global
         do i = 1, nx_global-1
            test_var(i,j)=0.5*(calc_var(i,j)+calc_var(i+1,j))
         enddo
      enddo
      do j = 1, ny_global
         test_var(nx_global,j)=0.5*(calc_var(nx_global,j)+calc_var(1,j))
      enddo
      call assert_array_equal(errorflag(itest), "dxU")
   endif

   ! dyE
   itest=itest+1
   call gather_global(grid_var_glob, dyE, master_task, distrb_info)
   if (my_task == master_task) then
      status = nf90_inq_varid(fid, "hte", varid) 
      if (status/=NF90_NOERR) call handle_nc_err(status, "Could not find variable hte in grid file")
      status = nf90_get_var(fid, varid, test_var)
      if (status/=NF90_NOERR) call handle_nc_err(status, "Could not load variable hte from grid file")
      test_var = test_var/100 !cm to m
      call assert_array_equal(errorflag(itest), "dyE")
   endif
   calc_var = test_var

   ! dyT
   itest=itest+1
   call gather_global(grid_var_glob, dyT, master_task, distrb_info)
   if (my_task == master_task) then
      do j=1,ny_global
         do i = 2, nx_global
            test_var(i,j)=0.5*(calc_var(i-1,j)+calc_var(i,j))
         enddo
      enddo
      do j=1,ny_global
         test_var(1,j)=0.5*(calc_var(nx_global,j)+calc_var(1,j))
      enddo
      call assert_array_equal(errorflag(itest), "dyT")
   endif

   ! dyN
   itest=itest+1
   call gather_global(grid_var_glob, dyN, master_task, distrb_info)
   if (my_task == master_task) then
      do j=1,ny_global-1
         do i = 2, nx_global
            test_var(i,j)=0.25*(calc_var(i,j)+calc_var(i,j+1)+calc_var(i-1,j)+calc_var(i-1,j+1))
         enddo
      enddo
      !i=1
      do j=1,ny_global
         test_var(1,j)=0.25*(calc_var(1,j)+calc_var(1,j+1)+calc_var(nx_global,j)+calc_var(nx_global,j+1))
      enddo
      !j=ny_global
      do i = 2, nx_global
         ! this approximation would be improved by cycling around for ns_boundary_type != open but does depend on the boundary type
         test_var(i,ny_global)=calc_var(i, ny_global-1)-0.5*calc_var(i, ny_global-2)+calc_var(i-1, ny_global-1)-0.5*calc_var(i-1, ny_global-2)
      enddo
      !i=1, j=ny_global
      test_var(1,ny_global)=calc_var(1, ny_global-1)-0.5*calc_var(1, ny_global-2)+calc_var(nx_global, ny_global-1)-0.5*calc_var(nx_global, ny_global-2)
      call assert_array_equal(errorflag(itest), "dyN")
   endif

   ! dyU
   itest=itest+1
   call gather_global(grid_var_glob, dyU, master_task, distrb_info)
   if (my_task == master_task) then
      do j=1,ny_global-1
         do i = 1, nx_global
            test_var(i,j)=0.5*(calc_var(i,j)+calc_var(i,j+1))
         enddo
      enddo
      ! j=ny_global
      do i = 1, nx_global
         ! there would be a better assumption for this line by closing the tripole instead
         test_var(i,ny_global)=2*calc_var(i,ny_global-1)-calc_var(i,ny_global-2)
      enddo
      call assert_array_equal(errorflag(itest), "dyU")
   endif


   ! --------------------
   ! This section tests the internally calculated cice cell areas
   ! ---------------------

   ! These are the fields in calculated by CICE in init_grid2 subroutine
   ! tarea, uarea, narea, earea

   area_vars(1)%xlen=dxT
   area_vars(1)%ylen=dyT
   area_vars(1)%area=tarea
   area_vars(1)%name='tarea'
   area_vars(2)%xlen=dxE
   area_vars(2)%ylen=dyE
   area_vars(2)%area=earea
   area_vars(2)%name='earea'
   area_vars(3)%xlen=dxN
   area_vars(3)%ylen=dyN
   area_vars(3)%area=narea
   area_vars(3)%name='narea'
   area_vars(4)%xlen=dxU
   area_vars(4)%ylen=dyU
   area_vars(4)%area=uarea
   area_vars(4)%name='uarea'


   if (my_task == master_task) then
      allocate(xlen(nx_global, ny_global))
      allocate(ylen(nx_global, ny_global))
   else
      allocate(xlen(1,1))
      allocate(ylen(1,1))
   endif

   do n = 1, nareavar
      itest=itest+1
      call gather_global(xlen, area_vars(n)%xlen, master_task, distrb_info)
      call gather_global(ylen, area_vars(n)%ylen, master_task, distrb_info)
      call gather_global(grid_var_glob, area_vars(n)%area, master_task, distrb_info)
      if (my_task == master_task) then
         test_var = xlen * ylen
         call assert_array_equal(errorflag(itest), area_vars(n)%name)
      endif
   enddo

   deallocate(xlen)
   deallocate(ylen)

   ! -------------------
   ! This section sanity tests the internally calculated lat/lons
   ! To-do: test the values cell by cell
   ! -------------------

   ! ELON / ELAT / NLON / NLAT
   itest=itest+1
   call gather_global(grid_var_glob, ELON, master_task, distrb_info)
   if (my_task==master_task) then
      if ((maxval(grid_var_glob)-minval(grid_var_glob))>2*PI) then
         write(6,*) "ELON out of range"
         errorflag(itest) = failflag
      endif
   endif

   itest=itest+1
   call gather_global(grid_var_glob, ELAT, master_task, distrb_info)
   if (my_task==master_task) then
      if ((maxval(grid_var_glob)-minval(grid_var_glob))>PI) then
         write(6,*) "ELAT out of range"
         errorflag(itest) = failflag
      endif
   endif

   itest=itest+1
   call gather_global(grid_var_glob, NLON, master_task, distrb_info)
   if (my_task==master_task) then
      if ((maxval(grid_var_glob)-minval(grid_var_glob))>2*PI) then
         write(6,*) "NLON out of range"
         errorflag(itest) = failflag
      endif
   endif

   itest=itest+1
   call gather_global(grid_var_glob, NLAT, master_task, distrb_info)
   if (my_task==master_task) then
      if ((maxval(grid_var_glob)-minval(grid_var_glob))>PI) then
         write(6,*) "NLAT out of range"
         errorflag(itest) = failflag
      endif
   endif

   ! ---------------------
   ! Add test for grid-cell corners?
   ! ---------------------

   deallocate(calc_var)
   deallocate(grid_var_glob)
   deallocate(test_var)
   deallocate(eps_var)

   if (my_task == master_task) then
      write(6,*) ' '
      ! write(6,*) 'GRIDLOADCHK COMPLETED SUCCESSFULLY'
      do n = 1,ntests
         if (errorflag(n)==failflag) errorflag0 = failflag
      enddo
      if (errorflag0 == passflag) then
         write(6,*) 'GRIDLOADCHK TEST COMPLETED SUCCESSFULLY'
      else
         write(6,*) 'GRIDLOADCHK TEST FAILED (AND COMPLETED SUCCESSFULLY)'
      endif
      write(6,*) ' '
      write(6,*) '=========================================================='
      write(6,*) ' '
      
   endif
#endif

   !-----------------------------------------------------------------
   ! Gracefully end
   !-----------------------------------------------------------------

   call dealloc_grid 
   call end_run()

end program

