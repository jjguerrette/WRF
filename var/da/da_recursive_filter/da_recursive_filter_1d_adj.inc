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

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d_adj")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   !-------------------------------------------------------------------------
   ! [4.0] Copy and tidy up:
   !-------------------------------------------------------------------------

   c(1:n) = field(1:n)
  
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

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d_adj")
   
end subroutine da_recursive_filter_1d_adj


subroutine da_recursive_filter_vec_1d_adj(npass, alpha, field2d, n1, n2)

   !--------------------------------------------------------------------------------------------
   ! Purpose: Perform npass passes of recursive filter across 1D segments of 2D array - adjoint.
   !--------------------------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: npass          ! Number of filter passes.
   real*8,  intent(in)    :: alpha          ! Alpha coefficient for RF.
   real*8,  intent(inout) :: field(n1,n2)   ! Array to be filtered.
   integer, intent(in)    :: n1,n2          ! Size of field2d array.

   integer :: i1, i2, pass   ! Loop counters.
   real    :: one_alpha      ! 1 - alpha.
   real    :: field_rl(n1)   ! Field after right-left pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_vec_1d_adj")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   do pass = npass, 1, -1
     do i2 = 1, n2

         !-------------------------------------------------------------------------
         ! [3.0] Perform left-moving filter:
         !-------------------------------------------------------------------------

         ! [3.2] Perform pass left to right:

         field_rl = 0.0   

         do i1 = 1, n1-1
            field2d(i1+1,i2) = field2d(i1+1,i2) + alpha * field2d(i1,i2)
            field_rl(i1) = one_alpha * field2d(i1,i2)
         end do

         ! use turning conditions as in the appendix of Hayden & Purser (1995):

         if (pass == 1) then
            field_rl(n1) = field_rl(n1) + field2d(n1,i2) / (1.0 + alpha)
         else
            field_rl(n1) = field_rl(n1) + one_alpha * field2d(n1,i2) / (1.0 - alpha**2)**2
            field_rl(n1-1) = field_rl(n1-1) - one_alpha * alpha**3 * field2d(n1,i2) / (1.0 - alpha**2)**2
         end if

         !-------------------------------------------------------------------------
         ! [2.0] Perform right-moving filter:
         !-------------------------------------------------------------------------

         field2d(:,i2) = 0.0

         ! [2.2] Perform pass left to right:

         do i1 = n1, 2, -1
            field_rl(i1-1) = field_rl(i1-1) + alpha * field_rl(i1)
            field2d(i1,i2) = field2d(i1,i2) + one_alpha * field_rl(i1)
         end do

         ! use turning conditions as in the appendix of Hayden & Purser (1995):

         if (pass == 1) then
            field2d(1,i2) = field2d(1,i2) + one_alpha * field_rl(1)
         else if (pass == 2) then
            field2d(1,i2) = field2d(1,i2) + field_rl(1) / (1.0 + alpha)
         else
            field2d(1,i2) = field2d(1,i2) + one_alpha * field_rl(1) / (1.0 - alpha**2)**2
            field2d(2,i2) = field2d(2,i2) - one_alpha * alpha**3 * field_rl(1) / (1.0 - alpha**2)**2
         end if

      end do !i2 loop
   end do !pass loop


   if (trace_use_dull) call da_trace_exit("da_recursive_filter_vec_1d_adj")
   
end subroutine da_recursive_filter_vec_1d_adj

