subroutine da_recursive_filter_1d_adj(pass, alpha, field, n)

   !---------------------------------------------------------------------------
   ! Purpose: Perform one pass of recursive filter on 1D array - adjoint.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: pass           ! Current pass of filter.
   real*8, intent(in)     :: alpha          ! Alpha coefficient for RF.
   real*8, intent(inout)  :: field(1:n)     ! Array to be filtered.
   integer, intent(in)    :: n              ! Size of field array.
   
   integer                :: j              ! Loop counter.
   real                   :: one_alpha      ! 1 - alpha.
   real                   :: a(1:n)         ! Input field.
   real                   :: b(1:n)         ! Field after left-right pass.
   real                   :: c(1:n)         ! Field after right-left pass.

!   character(len=10) :: this_time
!   character(len=8) :: this_date
!
!   integer :: dy_i, hr_i, mn_i, sc_i, msc_i, &
!              dy_f, hr_f, mn_f, sc_f, msc_f
!   real    :: time_sec_i, time_sec_f

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d_adj")

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_i
!read(this_time(1:2),'(I2.2)') hr_i
!read(this_time(3:4),'(I2.2)') mn_i
!read(this_time(5:6),'(I2.2)') sc_i
!read(this_time(8:10),'(I3.3)') msc_i
!time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   !-------------------------------------------------------------------------
   ! [4.0] Copy and tidy up:
   !-------------------------------------------------------------------------

   c(1:n) = field(1:n)

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_f
!read(this_time(1:2),'(I2.2)') hr_f
!read(this_time(3:4),'(I2.2)') mn_f
!read(this_time(5:6),'(I2.2)') sc_f
!read(this_time(8:10),'(I3.3)') msc_f
!time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0
!
!write(unit=message(1),fmt='(A,F12.4,A)') &
!   '   inside da_recursive_filter_1d_adj: Initialise: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
!call da_message(message(1:1))

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_i
!read(this_time(1:2),'(I2.2)') hr_i
!read(this_time(3:4),'(I2.2)') mn_i
!read(this_time(5:6),'(I2.2)') sc_i
!read(this_time(8:10),'(I3.3)') msc_i
!time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0

   
   !-------------------------------------------------------------------------
   ! [3.0] Perform left-moving filter:
   !-------------------------------------------------------------------------

   ! [3.2] Perform pass left to right:

   b(1:n) = 0.0   

   do j = 1, n-1
      c(j+1) = c(j+1) + alpha * c(j)
      b(j) = one_alpha * c(j)
   end do

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      b(n) = b(n) + c(n) / (1.0 + alpha)
   else
      b(n) = b(n) + one_alpha * c(n) / (1.0 - alpha**2)**2
      b(n-1) = b(n-1) - one_alpha * alpha**3 * c(n) / (1.0 - alpha**2)**2
   end if

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_f
!read(this_time(1:2),'(I2.2)') hr_f
!read(this_time(3:4),'(I2.2)') mn_f
!read(this_time(5:6),'(I2.2)') sc_f
!read(this_time(8:10),'(I3.3)') msc_f
!time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0
!
!write(unit=message(1),fmt='(A,F12.4,A)') &
!   '   inside da_recursive_filter_1d_adj: left-moving: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
!call da_message(message(1:1))

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_i
!read(this_time(1:2),'(I2.2)') hr_i
!read(this_time(3:4),'(I2.2)') mn_i
!read(this_time(5:6),'(I2.2)') sc_i
!read(this_time(8:10),'(I3.3)') msc_i
!time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0


   !-------------------------------------------------------------------------
   ! [2.0] Perform right-moving filter:
   !-------------------------------------------------------------------------

   a(1:n) = 0.0

   ! [2.2] Perform pass left to right:

   do j = n, 2, -1
      b(j-1) = b(j-1) + alpha * b(j)
      a(j) = a(j) + one_alpha * b(j)
   end do

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      a(1) = a(1) + one_alpha * b(1)
   else if (pass == 2) then
      a(1) = a(1) + b(1) / (1.0 + alpha)
   else
      a(1) = a(1) + one_alpha * b(1) / (1.0 - alpha**2)**2
      a(2) = a(2) - one_alpha * alpha**3 * b(1) / (1.0 - alpha**2)**2
   end if

   field(1:n) = a(1:n)

!call date_and_time(date=this_date, time=this_time)
!read(this_date(7:8),'(I2.2)') dy_f
!read(this_time(1:2),'(I2.2)') hr_f
!read(this_time(3:4),'(I2.2)') mn_f
!read(this_time(5:6),'(I2.2)') sc_f
!read(this_time(8:10),'(I3.3)') msc_f
!time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0

!write(unit=message(1),fmt='(A,F12.4,A)') &
!   '   inside da_recursive_filter_1d_adj: right-moving: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
!call da_message(message(1:1))

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d_adj")
   
end subroutine da_recursive_filter_1d_adj


