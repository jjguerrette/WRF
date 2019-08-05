subroutine da_perform_2drf(ni, nj, num_passes, rf_scale, field)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ni               ! Array dimension 1.
   integer, intent(in)    :: nj               ! Array dimension 2.
   integer, intent(in)    :: num_passes       ! Number of passes of RF.
   real,    intent(in)    :: rf_scale         ! Recursive filter scaling parameter.
   real*8,  intent(inout) :: field(1:ni,1:nj) ! Field to be filtered.

   integer               :: i, j, pass         ! Loop counters.
   real*8                :: e, alpha           ! Recursive filter parameters.
   real                  :: mean_field         ! Mean field.
   real*8                :: field_t(1:nj,1:ni) ! Transpose of field to be filtered.

   if (trace_use) call da_trace_entry("da_perform_2drf")

   e = 0.25 * num_passes / (rf_scale * rf_scale)
   alpha = 1 + e - sqrt(e * (e + 2.0))

   mean_field = sum(field(1:ni,1:nj)) / real(ni*nj)

   ! Perform filter in I-direction:
   call da_recursive_filter_1d(num_passes, alpha, field, ni, nj)

   field_t = TRANSPOSE(field)

   ! Perform filter in J-direction:
   call da_recursive_filter_1d(num_passes, alpha, field_t, nj, ni)

   field = TRANSPOSE(field_t)

   if (trace_use) call da_trace_exit("da_perform_2drf")

end subroutine da_perform_2drf


