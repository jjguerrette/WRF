subroutine da_get_field( input_file, var, field_dims, i_start, i_end, j_start, j_end, k,field)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

#include "netcdf.inc"
   
   character(len=200), intent(in)  :: input_file       ! NETCDF file nane.
   character(len=10),  intent(in)  :: var              ! Variable to search for.
   integer,            intent(in)  :: field_dims       ! # Dimensions of field. 
   integer,            intent(in)  :: i_start, i_end   ! Dimension 1 of field. 
   integer,            intent(in)  :: j_start, j_end   ! Dimension 2 of field. 
   integer,            intent(in)  :: k                ! Vertical index.
   real,               intent(out) :: field(i_start:i_end,j_start:j_end) ! Output field

   integer                   :: cdfid            ! NETCDF file id.
   integer                   :: rcode            ! Return code(0=ok).
   integer                   :: length           ! Length of filename.
   integer                   :: id_var           ! NETCDF variable ID. 

   integer                   :: istart(4)        ! Start value of arrays.
   integer                   :: icount(4)        ! Length of arrays.
!   integer                   :: iend(4)          ! End value of arrays.
   real(kind=4), allocatable :: field1d(:)       ! Used if 1D field read. 
   real(kind=4), allocatable :: field2d(:,:)     ! Used if 2D field read. 
   real(kind=4), allocatable :: field3d(:,:,:)   ! Used if 3D field read. 

   if (trace_use_dull) call da_trace_entry("da_get_field")

   length = len_trim(input_file)
   rcode = nf_open( input_file(1:length), NF_NOwrite, cdfid)
   if (rcode /= 0) then
      write(message(1),'(3a,i0)')' nf_open(',input_file(1:length),') returned ',rcode
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   !  Check variable is in file:
   rcode = nf_inq_varid( cdfid, var, id_var)
   if (rcode /= 0) then
      write(message(1),'(3a)')var,' variable is not in input file ',input_file(1:length)
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   istart = 1
   istart(1) = i_start
   icount(1) = i_end - i_start + 1
   istart(2) = j_start
   icount(2) = j_end - j_start + 1
   icount(4) = 1          ! Single time assumed.

   if (field_dims == 1) then
      icount(2) = 1
      icount(3) = 1
!      iend = istart + icount - 1
      allocate( field1d(1:icount(1)))
      call ncvgt( cdfid, id_var, istart, icount, field1d, rcode)
      field(:,1) = field1d(:)
      rcode = nf_close( cdfid)
      deallocate( field1d)
   else if (field_dims == 2) then
      icount(3) = 1
!      iend = istart + icount - 1
      allocate( field2d(1:icount(1),1:icount(2)))
      call ncvgt( cdfid, id_var, istart, icount, field2d, rcode)
      field(:,:) = field2d(:,:)
      rcode = nf_close( cdfid)
      deallocate( field2d)
   else if (field_dims == 3) then
      istart(3)=k
      icount(3)=1
!      iend = istart + icount - 1
      allocate( field3d(1:icount(1),1:icount(2),1:1))
      call ncvgt( cdfid, id_var, istart, icount, field3d(:,:,1:1), rcode)

      field(:,:) = field3d(:,:,1)
      deallocate( field3d)
   end if

   rcode = nf_close( cdfid)

   if (trace_use_dull) call da_trace_exit("da_get_field")

end subroutine da_get_field


