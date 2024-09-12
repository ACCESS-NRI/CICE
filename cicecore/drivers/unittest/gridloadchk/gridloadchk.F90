
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

         private
         public :: test_var_unchanged, t_grid_var

         real (kind=dbl_kind), dimension(:,:), public,  allocatable :: grid_var_glob, test_var
         real (kind=dbl_kind), dimension(:,:,:), public,  allocatable :: grid_var
         character(len=*), public,  parameter ::  &
            passflag = 'PASS', &
            failflag = 'FAIL'

         contains

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
            
            if (my_task == master_task) write (6,*) "LOG: attempting global gather "//trim(var%fn)
            call gather_global(grid_var_glob, var%gv, master_task, distrb_info)

            if (my_task == master_task) then
               if ( .not. var%optional) then 
                  status = nf90_get_var(fid, varid, test_var)
                  if (status/=NF90_NOERR) call handle_nc_err(status, "Could not load variable '"//trim(var%fn)//"' from grid file")

                  if (all(grid_var_glob==test_var)) then
                     write (6,*) "LOG: "//trim(var%fn)//" in cice matches original grid file"
                     errorflag = passflag
                  else
                     errorflag = failflag
                  endif
               endif
            endif

         end subroutine test_var_unchanged

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
      use test_grid_var

      use CICE_InitGrid
      use netcdf
      implicit none

      integer(kind=int_kind) :: n,m !,ny,nm,nd,nf1,nf2,xadd,nfa,nfb,nfc,ns1,ns2
      character(len=32) :: calstr !,unitstr,signstr

      integer(kind=int_kind), parameter :: ntests = 9 , ngridvar = 10
      character(len=8)  :: errorflag0,errorflag(1:ntests),errorflagtmp
      character(len=32) :: testname(ntests)
      integer :: fid , varid, status , test_nx, test_ny, itest
      integer, dimension(nf90_max_var_dims) :: dimids

      

      type(t_grid_var) :: grid_vars(ngridvar)


      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'RunningUnitTest GRIDLOADCHK'
         write(6,*) ' '
      endif
      write(6,*) ' my task: ', my_task !Why are these all 0 ?
      write(6,*) ' master task: ', master_task

      errorflag0   = passflag
      errorflag(:) = passflag
      testname(:) = ''
      ! testname(1) = 'compute_elapsed_days'

      call cice_init
      call init_grid1 ! domain distribution

      ! -------------------
      ! init_grid1 creates a distribution of blocks on each task
      ! this needs its own test
      ! -------------------
      call alloc_grid  
      call init_ice_timers      ! initialize all timers, this is needed by init_grid2
      ! call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! loads grid variables

! #ifndef USE_NETCDF

      ! --------------------
      ! This test confirms that 
      ! if we gather the grid from each task, this is unchanged from the original grid?
      ! ---------------------

      itest = 1 

      grid_vars(1)%fn = "ulat"
      grid_vars(1)%gv =  ULAT
      grid_vars(2)%fn = "ulon"
      grid_vars(2)%gv =  ULON
      grid_vars(3)%fn = "angle"
      grid_vars(3)%gv =  ANGLE
      grid_vars(3)%fn = "angleT"
      grid_vars(3)%gv =  ANGLET
      grid_vars(3)%optional = .True.
      

      if (my_task == master_task) then
         allocate(test_var(nx_global, ny_global))
         allocate(grid_var_glob(nx_global, ny_global))
      else
         allocate(test_var(1,1))
         allocate(grid_var_glob(1,1))
      endif

      ! add support for binary format?
      if (trim(grid_format) == 'nc') then

         if (my_task == master_task) write (6,*) "LOG: opening test instance of grid file"
         status = nf90_open(grid_file, NF90_NOWRITE, fid)
         if (status/=NF90_NOERR) call handle_nc_err(status, "Could not open "//grid_file//" for testing")

         do itest = 1 , 3

            if (my_task == master_task) write (6,*) "LOG: finding "//trim(grid_vars(iTest)%fn)//" var"
            status = nf90_inq_varid(fid, grid_vars(iTest)%fn, varid) 
            if (status/=NF90_NOERR) call handle_nc_err(status, "Could not find variable '"//trim(grid_vars(iTest)%fn)//"' in "//grid_file )

            status = nf90_inquire_variable(fid, varid, dimids = dimids)
            if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error ")
            status = nf90_inquire_dimension(fid, dimids(1), len = test_nx)
            if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error  ")
            status = nf90_inquire_dimension(fid, dimids(2), len = test_ny)
            if(status /= NF90_NOERR) call handle_nc_err(status, "Inq var error ")

            if (my_task == master_task) write (6,*) "LOG: checking size of "//trim(grid_vars(iTest)%fn)
            if (test_nx/=nx_global .or. test_ny/=ny_global) then
               if (my_task == master_task) then
                  errorflag(itest) = failflag
                  write (6,*) "Error in test ",itest," nx_global or ny_global from ice_in does not match grid file"
               endif
            else

               call test_var_unchanged(fid, varid, grid_vars(iTest), errorflag(itest)) 

            
            endif

         enddo

      endif


      deallocate(grid_var_glob)
      deallocate(test_var)

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'GRIDLOADCHK COMPLETED SUCCESSFULLY'
         do n = 1,ntests
            if (errorflag(n) == passflag) then
               write(6,*) 'PASS '!,trim(stringflag(n))
            else
               write(6,*) 'FAIL '!,trim(stringflag(n))
            endif
         enddo
         if (errorflag0 == passflag) then
            write(6,*) 'GRIDLOADCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'GRIDLOADCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif

      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      end program

