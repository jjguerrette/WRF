program gen_be_ep2
!
!---------------------------------------------------------------------- 
!  Purpose : To convert WRF ensemble to format required for use as 
!  flow-dependent perturbations in WRF-Var (alpha control variable, 
!  alphacv_method = 2).
!
!  Dale Barker (NCAR/MMM)      January 2007
!  Arthur P. Mizzi (NCAR/MMM)  February 2011  Modified to use .vari extension for
!                                             ensemble variance file output from
!                                             gen_be_ensmean.f90
!
!  JJ Guerrette (NCAR/MMM)     March 2019     Refactored for generalization across 
!                                             different variables, enabled (again?) 
!                                             usage of ensemble mean/variance files
!                                             output from gen_be_ensmean.f90, added
!                                             MPI distribution of grid
!----------------------------------------------------------------------

#ifdef crayx1
#define iargc ipxfargc
#endif

   use da_control, only : stderr, stdout, filename_len, root
   use da_gen_be, only : da_stage0_initialize
!#ifdef DM_PARALLEL
!   use module_dm, only : MPASPECT
!#endif

   implicit none

#ifdef DM_PARALLEL
   include 'mpif.h'
#endif
   type field_type
      real, allocatable :: data(:,:,:)
   end type field_type

   character (len=filename_len)   :: directory        ! General filename stub.
   character (len=filename_len)   :: filestub         ! General filename stub.
   character (len=filename_len)   :: filebase   ! General filename base.
   character (len=filename_len)   :: input_file       ! Input file.
   character (len=10)    :: date                      ! Character date.
   character (len=10)    :: varname                   ! Variable to search for.
   character (len=10), allocatable  :: varnames(:)               ! Variables to work on

   character (len=3)     :: cne                       ! Ensemble size.
   character (len=3)     :: pe                    ! Member/processor indices -> character.

   integer, external     :: iargc
   integer               :: numarg
   integer               :: ne                        ! Ensemble size
   integer               :: ivar                      ! Loop counter.
   integer               :: ni, nj, nk                ! Global dimensions of grid (T points).
   integer               :: ips, ipe, jps, jpe        ! Local patch dimensions
   integer               :: mp_physics                ! microphysics option
   real                  :: ds                        ! Grid resolution.
   logical               :: remove_mean               ! Remove mean from standard fields.

   integer               :: nvars
   character (len=10), allocatable :: select_vars(:)
   integer, allocatable            :: var_nk(:)

   integer :: cdfid            ! NETCDF file id.
   integer :: rcode            ! Return code (0=ok).

   integer              :: num_procs, myproc, ierr, ntasks_x, ntasks_y
   logical              :: log_proc
   integer, allocatable :: nip(:), njp(:), &
                           ips_a(:), ipe_a(:), &
                           jps_a(:), jpe_a(:)

   stderr = 0
   stdout = 6

#ifdef DM_PARALLEL
   call mpi_init(ierr)
   call mpi_comm_size(mpi_comm_world,num_procs,ierr)
   call mpi_comm_rank(mpi_comm_world,myproc,ierr)
#else
   num_procs = 1
   myproc = 0
#endif
   write(pe,'(i3.3)')myproc

   log_proc = ( myproc == root )
   if (log_proc) write(stdout,'(A,I0)') 'num_procs = ', num_procs

!---------------------------------------------------------------------------------------------
   if (log_proc) write(stdout,'(/a)')' [A] Initialize information.'
!---------------------------------------------------------------------------------------------

   remove_mean = .true.

   numarg = iargc()
   if ( numarg /= 4 )then
      if (log_proc) write(stdout,'(a)') &
        "Usage: gen_be_ep2 date ne <directory> <filename> Stop"
      stop
   end if

   ! Initialse to stop Cray compiler complaining
   date=""
   cne=""
   directory=""
   filestub=""

   call getarg( 1, date )
   call getarg( 2, cne )
   read(cne,'(i3)')ne
   call getarg( 3, directory )
   call getarg( 4, filestub )

   if ( remove_mean ) then
      if (log_proc) write(stdout,'(a,a)') &
      ' Computing gen_be ensemble perturbation files for date ', date
   else
      if (log_proc) write(stdout,'(a,a)') &
      ' Computing gen_be ensemble forecast files for date ', date
   end if
   if (log_proc) write(stdout,'(a)') &
      ' Perturbations are in MODEL SPACE (u, v, t, q, ps)'
   if (log_proc) write(stdout,'(a,i4)') &
      ' Ensemble Size = ', ne
   if (log_proc) write(stdout,'(a,a)') &
      ' Directory = ', trim(directory)
   if (log_proc) write(stdout,'(a,a)') &
      ' Filename = ', trim(filestub)

!---------------------------------------------------------------------------------------------
   if (log_proc) write(stdout,'(/a)') &
      ' [B] Set up data dimensions and allocate data arrays:' 
!---------------------------------------------------------------------------------------------

   ! select_vars could come from namelist or external caller. Hard-wired for now.
   nvars = 10
   allocate(select_vars(nvars))
   select_vars(1:nvars) = &
      [ character(len=10) :: "U", "V", "T", "QVAPOR", "PSFC", & 
                   "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP" ]

   ! Get domain grid dimensions from first T field:
   varname = "T"
   filebase = trim(directory)//'/'//trim(filestub)
   input_file = trim(filebase)//'.e001'
   call da_stage0_initialize( input_file, varname, ni, nj, nk, ds, mp_physics )

   allocate(var_nk(nvars))
   var_nk = (/nk, nk, nk, nk, 1, &
              nk, nk, nk, nk, nk/)

   ! Assign local patch grid dimensions using num_procs, ni, nj
   ! Note: the patch and domain limits could alternatively come from WRF framework
   !       in an on-line WRFDA calc. for perturbations, would not need to write 
   !       pert files
#ifdef DM_PARALLEL
   allocate( nip(num_procs)  , njp(num_procs)   )
   allocate( ips_a(num_procs), ipe_a(num_procs) )
   allocate( jps_a(num_procs), jpe_a(num_procs) )

   call MPASPECT   ( num_procs, ntasks_x, ntasks_y, 1, 1)
   call split_grid ( num_procs, ntasks_x, ntasks_y, ni, nj, nip, njp, ips_a, jps_a)
   ipe_a = ips_a + nip - 1
   jpe_a = jps_a + njp - 1

   if (log_proc) write(stdout,'(A,3(I8))') ' ni, nj, nk ', ni, nj, nk
   if (log_proc) write(stdout,'(A,*(I8))') ' ips_a = ', ips_a
   if (log_proc) write(stdout,'(A,*(I8))') ' ipe_a = ', ipe_a
   if (log_proc) write(stdout,'(A,*(I8))') ' jps_a = ', jps_a
   if (log_proc) write(stdout,'(A,*(I8))') ' jpe_a = ', jpe_a

   ips = ips_a(myproc+1)
   ipe = ipe_a(myproc+1)
   jps = jps_a(myproc+1)
   jpe = jpe_a(myproc+1)
   deallocate(nip, njp)
   deallocate(ips_a, ipe_a)
   deallocate(jps_a, jpe_a)
#else
   ips = 1
   ipe = ni
   jps = 1
   jpe = nj
#endif

!---------------------------------------------------------------------------------------------
   if (log_proc) write(stdout,'(/a)') &
      ' [C] Loop over variables:' 
!---------------------------------------------------------------------------------------------

   call gen_be_ep( filebase, select_vars, nvars, remove_mean, .true., &
                   ni, nj, var_nk, ne, ips, ipe, jps, jpe, &
                   num_procs, myproc )

!   allocate(varnames(1))
!   do ivar = 1, nvars
!      varnames(1) = trim(select_vars(ivar))
!
!Alternatively loops over variables for writing
!      call gen_be_ep( filebase, varnames, 1, remove_mean, .true., &
!                      ni, nj, var_nk(ivar:ivar), ne, ips, ipe, jps, jpe, &
!                      num_procs, myproc )
!
!OR read them directly into WRFDA (da_setup_flow_predictors.inc)
!
!      !Note: need to update filebase if time stamp changes (4DEnVar)
!
!      allocate(dummy(ips:ipe, jps:jpe, 1:var_nk(ivar), ne)
!      call gen_be_ep( filebase, varnames, 1, .true., .false., &
!                      ni, nj, var_nk(ivar:ivar), ne, ips, ipe, jps, jpe, &
!                      num_procs, myproc, model_field = dummy(:,:,:,1:ne)  )
!
!      do iens = 1, ne
!         if (varnames(1) == "U") then
!            ep % v1(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!         end if
!         if (varnames(1) == "V") then
!            ep % v2(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!         end if
!         if (varnames(1) == "T") then
!            ep % v3(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!         end if
!         if (varnames(1) == "QVAPOR") then
!            ep % v4(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!         end if
!         if (varnames(1) == "PSFC") then
!            ep % v5(ips:ipe,jps:jpe,   1,ie) = dummy(:,:,1,ie)
!         end if
!
!         !  Optional include hydrometeors: 
!         !  (could adjust varnames based on these options instead of if statement)
!         if ( alphacv_method == alphacv_method_xa .and. alpha_hydrometeors ) then  ! xa space
!            if (varnames(1) == "QCLOUD") then
!               ep % cl(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!            end if
!            if (varnames(1) == "QRAIN") then
!               ep % rn(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!            end if
!            if (varnames(1) == "QICE") then
!               ep % ci(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!            end if
!            if (varnames(1) == "QSNOW") then
!               ep % sn(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!            end if
!            if (varnames(1) == "QGRAUP") then
!               ep % gr(ips:ipe,jps:jpe,1:nk,ie) = dummy(:,:,:,ie)
!            end if
!         end if
!      end do
!
!      deallocate(dummy)
!
!   end do
!   deallocate(varnames)


   deallocate ( select_vars )
   deallocate ( var_nk )


#ifdef DM_PARALLEL
   call mpi_finalize(ierr)
#endif

contains

!-----------------------------------------------------------------------------

#ifdef crayx1
   subroutine getarg(i, harg)
     implicit none
     character(len=*) :: harg
     integer :: ierr, ilen, i

     call pxfgetarg(i, harg, ilen, ierr)
     return
   end subroutine getarg
#endif

!-----------------------------------------------------------------------------

#ifdef DM_PARALLEL
   !copied from module_dm.f90/da_module_dm due to linking issue with module_dm.o
   subroutine mpaspect( p, minm, minn, procmin_m, procmin_n )
      implicit none
      integer p, m, n, mini, minm, minn, procmin_m, procmin_n
      mini = 2*p
      minm = 1
      minn = p
      do m = 1, p
        if ( mod( p, m ) .eq. 0 ) then
          n = p / m
          if ( abs(m-n) .lt. mini                &
               .and. m .ge. procmin_m            &
               .and. n .ge. procmin_n            &
             ) then
            mini = abs(m-n)
            minm = m
            minn = n
          end if
        end if
      end do
!      if ( minm .lt. procmin_m .or. minn .lt. procmin_n ) then
!        write( wrf_err_message , * )'mpaspect: unable to generate processor mesh.  stopping.'
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        write( wrf_err_message , * )' procmin_m ', procmin_m
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        write( wrf_err_message , * )' procmin_n ', procmin_n
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        write( wrf_err_message , * )' p         ', p
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        write( wrf_err_message , * )' minm      ', minm
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        write( wrf_err_message , * )' minn      ', minn
!        call da_wrf_message ( trim ( wrf_err_message ) )
!        call da_wrf_error_fatal3("<stdin>",127,&
!'da_module_dm: mpaspect' )
!      end if
   return
   end subroutine mpaspect

   subroutine split_grid( nprocs, ntx, nty, &
                          nx_global, ny_global, &
                          nx_grid,   ny_grid, &
                          xs_grid,   ys_grid )
      implicit none

      integer, intent(in)  :: nprocs                               ! Total number of processors
      integer, intent(in)  :: ntx, nty                             ! Number of tasks in x-y decomposition
      integer, intent(in)  :: nx_global, ny_global                 ! Number of grid points in domain
      integer, intent(out) :: nx_grid(nprocs), ny_grid(nprocs), &  ! Number of grid points in x-y patches
                              xs_grid(nprocs), ys_grid(nprocs)     ! Starting grid points in x-y patches

      integer, target   :: nx_vec(ntx), xs_vec(ntx)
      integer, target   :: ny_vec(nty), ys_vec(nty)
      integer, pointer  :: nvec(:), svec(:)
      integer           :: mm, i, j, ii, iproc, igrid, ntasks, nglobal, fact

      do igrid = 1, 2
         if (igrid.eq.1) then
            nvec => ny_vec
            svec => ys_vec
            ntasks = nty
            nglobal = ny_global
         else if (igrid.eq.2) then
            nvec => nx_vec
            svec => xs_vec
            ntasks = ntx
            nglobal = nx_global
         end if

         nvec = nglobal / ntasks
         mm = mod( nglobal , ntasks )
         do j = 1, ntasks
            if ( mm .eq. 0 ) exit
            nvec(j) = nvec(j) + 1
            mm = mm - 1
         end do

         svec(1) = 1
         do j = 1, ntasks
            if (j .lt. ntasks) then
               svec(j+1) = svec(j) + nvec(j)
            end if
         end do
      end do

      iproc = 0
      do j = 1, nty
      do i = 1, ntx
         iproc = iproc + 1
         ny_grid(iproc) = ny_vec(j)
         ys_grid(iproc) = ys_vec(j)
         nx_grid(iproc) = nx_vec(i)
         xs_grid(iproc) = xs_vec(i)
      end do
      end do

   end subroutine split_grid
#endif

!-----------------------------------------------------------------------------

   subroutine gen_be_ep( filebase, varnames, nvars, remove_mean, write_binary, &
                         ni, nj, nks, ne, ips, ipe, jps, jpe, &
                         num_procs, myproc, &
                         model_field, model_field_mean, model_field_stdv )

      !EXTERNAL DEPENDENCIES:  root, stdout

      use da_reporting, only : da_error,message

      implicit none

      include 'netcdf.inc'

      character (len=filename_len), intent(in) :: filebase    ! General filename base.
      character (len=10), intent(in)           :: varnames(nvars) ! Variable to search for.
      logical, intent(in) :: remove_mean        ! Remove mean from standard fields.
      logical, intent(in) :: write_binary       ! Write binary output for all fields.
      integer, intent(in) :: ne, nvars          ! Ensemble/variable size
      integer, intent(in) :: ni, nj             ! Global dimensions of grid (T points).
      integer, intent(in) :: nks(nvars)         ! Array of # levels
      integer, intent(in) :: ips, ipe, jps, jpe ! Local patch dimensions
      integer, intent(in) :: num_procs, myproc
      real, optional, intent(out) :: model_field      ( ips:ipe, jps:jpe, nks(1), ne, nvars )
      real, optional, intent(out) :: model_field_mean ( ips:ipe, jps:jpe, nks(1), nvars )
      real, optional, intent(out) :: model_field_stdv ( ips:ipe, jps:jpe, nks(1), nvars )

      character (len=filename_len)   :: input_file       ! Input file.
      character (len=3)     :: ce                        ! Member index -> character.
      integer               :: i, j, k, imem, ivar, jvar ! Loop counters.
      real                  :: imem_inv                  ! 1 / imem.

      integer, parameter    :: nvars_max = 10
      character (len=10)    :: var_shortnames(nvars_max), var_longnames(nvars_max)
      logical               :: var_is_meteor(nvars_max)

      integer               :: this_var

      integer               :: ndims(nvars)
      character (len=10)    :: shortnames(nvars)
      character (len=10)    :: longnames(nvars)
      type(field_type)      :: member_fields(nvars)
      type(field_type)      :: mean_fields(nvars)
      type(field_type)      :: stdv_fields(nvars)

      integer :: cdfid            ! NETCDF file id.
      integer :: id_var
      integer :: rcode, rcode_     ! Return code (0=ok).
      logical :: mean_vari_present
      logical :: log_proc

      log_proc = ( myproc == root )

      ! long variable names
      var_longnames(1:nvars_max) = &
         [ character(len=10) :: "U", "V", "T", "QVAPOR", "PSFC", & 
                      "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP" ]

      ! short variable names
      var_shortnames(1:nvars_max) = &
         [ character(len=10) :: "u", "v", "t", "q", "ps", &
                      "qcloud", "qrain", "qice", "qsnow", "qgraup" ]

      var_is_meteor = (/.false., .false., .false., .false., .false., &
                        .true. , .true. , .true. , .true. , .true. /)

      this_var = 0

      do ivar = 1, nvars
         do jvar = 1, nvars_max
            if (trim(varnames(ivar)) == trim(var_longnames(jvar))) then
               this_var = jvar
               exit
            end if
         end do
         if (this_var < 1) then
            write(message(1),'(A,A)') &
            varnames(ivar), ' variable is not an option in gen_be_ep2'
            if (log_proc) write(stdout,'(A)') message(1)
            cycle
         end if

         ! Check variable is in e001 file:
         input_file = trim(filebase)//'.e001'
         rcode_ = nf_open(input_file(1:len_trim(input_file)), NF_NOwrite, cdfid)
         rcode = nf_inq_varid ( cdfid, varnames(ivar), id_var )
         rcode_ = nf_close(cdfid)
         if ( rcode /= 0 ) then
            write(message(1),'(A,A)') &
               varnames(ivar), ' variable is not in input file'
            if ( var_is_meteor(this_var) ) then
               if (log_proc) write(stdout,'(A)') message(1)
               cycle
            else
               call da_error(__FILE__,__LINE__,message(1:1))
            end if
         end if

!---------------------------------------------------------------------------------------------
         if (log_proc) write(stdout,'(/2A)') &
            ' [0] Allocating fields for variable ',trim(varnames(ivar))
!---------------------------------------------------------------------------------------------
         ! cloud variable determination 
         if ( var_is_meteor(this_var) .and. log_proc) write(stdout,'(A)')&
            '          [cloud variable]'

         shortnames(ivar) = var_shortnames(this_var)
         longnames(ivar)  = var_longnames(this_var)
         ndims(ivar) = 3
         if (nks(ivar) == 1) ndims(ivar) = 2

         ! Allocate data arrays in output fields:
         allocate( member_fields(ivar) % data (ips:ipe, jps:jpe, 1:nks(ivar) ) )
         allocate( mean_fields(ivar) % data   (ips:ipe, jps:jpe, 1:nks(ivar) ) )
         allocate( stdv_fields(ivar) % data   (ips:ipe, jps:jpe, 1:nks(ivar) ) )
         member_fields(ivar) % data = 0.0
         mean_fields(ivar) % data = 0.0
         stdv_fields(ivar) % data = 0.0
      end do !ivar loop

!---------------------------------------------------------------------------------------------
      if (log_proc) write(stdout,'(/a)')&
         ' [1] Generate mean and variance (read or compute)'
!---------------------------------------------------------------------------------------------

      !If mean/vari files exist, use those instead of calculating mean and stddev online
      mean_vari_present = .false.

      input_file = trim(filebase)//'.mean'
      rcode = nf_open(input_file(1:len_trim(input_file)), NF_NOwrite, cdfid)
      if (rcode /= 0) then
         if (log_proc) write(stdout,'(A)') &
            ' Mean file does not exist, computing mean/variance.'
      else
         rcode = nf_close(cdfid)
         input_file = trim(filebase)//'.vari'
         rcode = nf_open(input_file(1:len_trim(input_file)), NF_NOwrite, cdfid)
         if (rcode /= 0) then
            if (log_proc) write(stdout,'(A)') &
               ' Variance file does not exist, computing mean/variance.'
         else
            rcode = nf_close(cdfid)
            mean_vari_present = .true.
         end if
      end if


      if (mean_vari_present) then
         if (log_proc) write(stdout,'(A)') &
            '         Reading vars from mean file...'
         input_file = trim(filebase)//'.mean'
         do ivar = 1, nvars
            call read_field_from_nc( input_file, mean_fields(ivar) % data, longnames(ivar), &
                                     ips, ipe, jps, jpe, nks(ivar), ndims(ivar) )
         end do

         if (log_proc) write(stdout,'(A)') &
            '         Reading vars from variance file...'
         input_file = trim(filebase)//'.vari'
         do ivar = 1, nvars
            call read_field_from_nc( input_file, stdv_fields(ivar) % data, longnames(ivar), &
                                     ips, ipe, jps, jpe, nks(ivar), ndims(ivar) )
            stdv_fields(ivar) % data = sqrt( stdv_fields(ivar) % data )
         end do
      else
         do imem = 1, ne
            write(ce,'(i3.3)')imem
            if (log_proc) write(stdout,'(2A)') &
               '         Working on ensemble member ',trim(ce)

            input_file = trim(filebase)//'.e'//trim(ce)
            do ivar = 1, nvars
               call read_field_from_nc( input_file, member_fields(ivar) % data, longnames(ivar), &
                                        ips, ipe, jps, jpe, nks(ivar), ndims(ivar) )
               ! Sum values and square values:
               mean_fields(ivar) % data = &
                    mean_fields(ivar) % data + member_fields(ivar) % data

               stdv_fields(ivar) % data = &
                    stdv_fields(ivar) % data + member_fields(ivar) % data * member_fields(ivar) % data
            end do
         end do

         ! Finalize mean and stdv
         imem_inv = 1.0 / real ( ne )
         do ivar = 1, nvars
            mean_fields(ivar) % data = mean_fields(ivar) % data * imem_inv
            stdv_fields(ivar) % data = stdv_fields(ivar) % data * imem_inv

            stdv_fields(ivar) % data = &
                sqrt( stdv_fields(ivar) % data - mean_fields(ivar) % data * mean_fields(ivar) % data )
         end do
      end if ! mean_vari_present

!---------------------------------------------------------------------------------------------
      if (log_proc) write(stdout,'(/a)') &
         ' [2] Compute perturbations and output' 
!---------------------------------------------------------------------------------------------

      if ( remove_mean ) then
         if (log_proc) write(stdout,'(a)') &
      "  Calculate ensemble perturbations"
      else
         if (log_proc) write(stdout,'(a)') &
      "  WARNING: Not removing ensemble mean (outputs are full-fields)"
      end if

      do imem = 1, ne
         write(ce,'(i3.3)')imem
         if (log_proc) write(stdout,'(2A)') &
            '         Working on ensemble member ',trim(ce)
         if (log_proc) write(stdout,'(A)')  &
            '          Reading member field...'

         input_file = trim(filebase)//'.e'//trim(ce)
         do ivar = 1, nvars
            call read_field_from_nc( input_file, member_fields(ivar) % data, longnames(ivar), &
                                     ips, ipe, jps, jpe, nks(ivar), ndims(ivar) )

            if ( remove_mean ) then
               member_fields(ivar) % data = member_fields(ivar) % data - mean_fields(ivar) % data
            end if

            if ( present(model_field) ) then
               model_field(:,:,1:nks(ivar),imem,ivar) = member_fields(ivar) % data
            end if

!            if ( write_binary ) then
!               call write_fields_to_bin( member_fields(ivar:ivar), shortnames(ivar:ivar), &
!                                        '.e'//trim(ce),&
!                                        num_procs, myproc, 1, &
!                                        ni, nj, nks(ivar:ivar), ips, ipe, jps, jpe )
!            end if
         end do
         if ( write_binary ) then
            if (log_proc) write(stdout,'(A)') &
               '         Writing fields...'
            call write_fields_to_bin( member_fields, shortnames, &
                                     '.e'//trim(ce),&
                                     num_procs, myproc, nvars, &
                                     ni, nj, nks, ips, ipe, jps, jpe )
         end if
      end do

      do ivar = 1, nvars
         if ( present(model_field_mean) ) then
            model_field_mean(:,:,1:nks(ivar),ivar) = mean_fields(ivar) % data
         end if
         if ( present(model_field_stdv) ) then
            model_field_stdv(:,:,1:nks(ivar),ivar) = stdv_fields(ivar) % data
         end if
!         if ( write_binary ) then
!            ! Write out mean/stdv field (stdv stored in stdv data arrays):
!            call write_fields_to_bin( mean_fields(ivar:ivar), shortnames(ivar:ivar), &
!                                     '.mean', &
!                                     num_procs, myproc, 1, &
!                                     ni, nj, nks(ivar:ivar), ips, ipe, jps, jpe )
!
!            call write_fields_to_bin( stdv_fields(ivar:ivar), shortnames(ivar:ivar), &
!                                     '.stdv', &
!                                     num_procs, myproc, 1, &
!                                     ni, nj, nks(ivar:ivar), ips, ipe, jps, jpe )
!         end if

      end do
      if ( write_binary ) then
         ! Write out mean/stdv field (stdv stored in stdv data arrays):
         if (log_proc) write(stdout,'(A)') &
            '         Writing mean/stdv field...'
         call write_fields_to_bin( mean_fields, shortnames, &
                                  '.mean', &
                                  num_procs, myproc, nvars, &
                                  ni, nj, nks, ips, ipe, jps, jpe )

         call write_fields_to_bin( stdv_fields, shortnames, &
                                  '.stdv', &
                                  num_procs, myproc, nvars, &
                                  ni, nj, nks, ips, ipe, jps, jpe )
      end if



      ! Deallocate data arrays
      do ivar = 1, nvars
         deallocate( member_fields(ivar) % data )
         deallocate( mean_fields(ivar) % data   )
         deallocate( stdv_fields(ivar) % data   )
      end do

   end subroutine gen_be_ep

!-----------------------------------------------------------------------------

   subroutine read_field_from_nc( read_file, field, varname, &
                                  ips, ipe, jps, jpe, nk, ndims )

      use da_gen_be, only : da_get_field, da_get_trh

      implicit none

      character (len=filename_len), intent(in) :: read_file
      character (len=10), intent(in)   :: varname
      integer, intent(in)              :: ips, ipe, jps, jpe, nk, ndims
      real, intent(inout)  :: field( ips:ipe, jps:jpe, 1:nk )

      real, allocatable  :: ttmp(:,:)                 ! temperature.
      real, allocatable  :: dummy(:,:)                ! dummy.
      integer :: i, j, k, ierr


      do k = 1, nk !Loop over this field''s 3rd dimension
         if ( trim(varname) == "U" ) then
            ! Read and interpolate u to mass pts:
            if (k==1) allocate( dummy ( ips:ipe+1, jps:jpe)   ) ! u on Arakawa C-grid.
            call da_get_field( read_file, varname, ndims, &
                               ips, ipe+1, jps, jpe, k, dummy )
            do j = jps, jpe
               do i = ips, ipe
                  field (i,j,k) = 0.5 * ( dummy(i,j) + dummy(i+1,j) )
               end do
            end do
         else if ( trim(varname) == "V" ) then
            ! Read and interpolate v to mass pts:
            if (k==1) allocate( dummy ( ips:ipe,   jps:jpe+1) ) ! v on Arakawa C-grid.
            call da_get_field( read_file, varname, ndims, &
                               ips, ipe, jps, jpe+1, k, dummy )
            do j = jps, jpe
               do i = ips, ipe
                  field (i,j,k) = 0.5 * ( dummy(i,j) + dummy(i,j+1) )
               end do
            end do
         else if ( trim(varname) == "T" ) then
            ! Read theta, and convert to temperature:
            if (k==1) allocate( ttmp ( ips:ipe,   jps:jpe)   )
            if (k==1) allocate( dummy( ips:ipe,   jps:jpe)   )
            call da_get_trh( read_file, ips, ipe, jps, jpe, k, ttmp, dummy)
            field (:,:,k) = ttmp(:,:)
         else
            ! Read all other variables
            if (k==1) allocate( dummy( ips:ipe,   jps:jpe)   )
            call da_get_field( read_file, varname, ndims, &
                               ips, ipe, jps, jpe, k, dummy )
            if ( trim(varname) == "QVAPOR" ) then
               ! Convert mixing ratio to specific humidity:
               field (:,:,k) = dummy(:,:) / ( 1.0 + dummy(:,:) )
            else
               ! 1-to-1 correspondance between netcdf field and analyzed variable
               field (:,:,k) = dummy(:,:)
            end if
         end if
      end do

      ! Deallocate temporary arrays:
      if (allocated(ttmp)) deallocate( ttmp )
      if (allocated(dummy)) deallocate( dummy )

!#ifdef DM_PARALLEL
!      call mpi_barrier(MPI_COMM_WORLD,ierr)
!#endif

   end subroutine read_field_from_nc

!-----------------------------------------------------------------------------

   subroutine write_fields_to_bin( fields_local, varnames, suf, &
                                  num_procs, myproc, nvars, &
                                  ni, nj, nks, ips, ipe, jps, jpe )

      use da_tools_serial, only : da_get_unit, da_free_unit
      use da_par_util1, only: true_mpi_real

      implicit none

      type(field_type), intent(in)  :: fields_local(nvars)
      character (len=*), intent(in) :: varnames(nvars), suf
      integer, intent(in)           :: num_procs, myproc, nvars
      integer, intent(in)           :: ni, nj, nks(nvars), ips, ipe, jps, jpe

      integer :: write_proc
      integer            :: write_unit, k !, j
      real, allocatable  :: globbuf(:,:,:), locbuf(:)

#ifdef DM_PARALLEL
      integer :: ips_, ipe_, jps_, jpe_
      integer :: icount, jcount, tag, iproc, ierr
      integer :: status(MPI_STATUS_SIZE)
#endif

      do ivar = 1, nvars
         write_proc = mod(ivar - 1, num_procs)

         if (myproc == write_proc) then
            ! Create file and allocate global buffer
            call da_get_unit(write_unit)

            open (write_unit, &
                  file = trim(varnames(ivar))//trim(suf), &
                  form='unformatted')
            allocate( globbuf(1:ni, 1:nj, 1:nks(ivar)) )
         end if

#ifdef DM_PARALLEL
         if (num_procs > 1) then
            !Communicate fields_local segments to globbuf
            do iproc = 0, num_procs - 1
               ips_ = ips
               ipe_ = ipe
               jps_ = jps
               jpe_ = jpe
               call mpi_bcast(ips_, 1, mpi_integer, iproc, MPI_COMM_WORLD, ierr)
               call mpi_bcast(ipe_, 1, mpi_integer, iproc, MPI_COMM_WORLD, ierr)
               call mpi_bcast(jps_, 1, mpi_integer, iproc, MPI_COMM_WORLD, ierr)
               call mpi_bcast(jpe_, 1, mpi_integer, iproc, MPI_COMM_WORLD, ierr)

               icount = ipe_ - ips_ + 1
               jcount = jpe_ - jps_ + 1
!!!               allocate( locbuf(1:icount) )
               allocate( locbuf(1:icount*jcount) )
               do k = 1, nks(ivar)
!!!                  do j = jps, jpe
!!!                     tag = 20*j+300000*k
                     tag = iproc+10000*k
                     if (myproc == iproc) then
                        if (myproc == write_proc) then
                           globbuf( ips_:ipe_, jps_:jpe_, k ) = fields_local(ivar) % data(:,:, k )
                        else
                   
!!!                        locbuf = fields_local(ivar) % data( ips:ipe, j, k )
!!!                        call MPI_SEND( locbuf, icount, true_mpi_real, write_proc, tag, MPI_COMM_WORLD, ierr)
                           locbuf = reshape( fields_local(ivar) % data(:,:, k ), (/ icount * jcount /) )
                           call MPI_SEND( locbuf, icount * jcount, true_mpi_real, &
                                          write_proc, tag, MPI_COMM_WORLD, ierr)
                        end if
                     else if (myproc == write_proc) then
!!!                        call MPI_RECV( locbuf, icount, true_mpi_real, iproc,      tag, MPI_COMM_WORLD, status, ierr)
!!!                        globbuf( ips:ipe, j, k ) = locbuf
                        call MPI_RECV( locbuf, icount * jcount, true_mpi_real, &
                                       iproc,      tag, MPI_COMM_WORLD, status, ierr)
                        globbuf( ips_:ipe_, jps_:jpe_, k ) = reshape( locbuf, (/icount, jcount/) )
                     end if
!!!                  end do

               end do
               deallocate( locbuf )
            end do

!            call mpi_barrier(MPI_COMM_WORLD,ierr)

         else
#endif
            globbuf = fields_local(ivar) % data
#ifdef DM_PARALLEL
         end if
#endif
         if (myproc == write_proc) then
            ! Write to file, close file, and deallocate buffer
            write(write_unit) ni, nj, nks(ivar)
            write(write_unit) globbuf

            close(write_unit)

            deallocate( globbuf )

            call da_free_unit(write_unit)
         end if
      end do

   end subroutine write_fields_to_bin

!-----------------------------------------------------------------------------

!   subroutine logger(mess, fmt, proc)
!
!      implicit none
!      character(len=*) :: mess
!      character(len=*) :: fmt
!      integer, parameter :: log_proc = root
!
!      if (proc == log_proc) write(stdout,fmt) mess
!
!   end subroutine

end program gen_be_ep2
