subroutine da_recursive_filter_1d(pass, alpha, field, n)

   !---------------------------------------------------------------------------
   ! Purpose: Perform one pass of recursive filter on 1D array.
   !
   ! Method:  Perform right-moving filter followed by left-moving filter.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: pass           ! Current pass of filter.
   real*8,  intent(in)    :: alpha          ! Alpha coefficient for RF.
   real*8,  intent(inout) :: field(1:n)     ! Array to be filtered.
   integer, intent(in)    :: n              ! Size of field array.

   integer :: j              ! Loop counter.
   real    :: one_alpha      ! 1 - alpha.
   real    :: a(1:n)         ! Input field.
   real    :: b(1:n)         ! Field after left-right pass.
   real    :: c(1:n)         ! Field after right-left pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d")
   
   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha
   
   a(1:n) = field(1:n)

   !-------------------------------------------------------------------------
   ! [2.0] Perform right-moving filter:
   !-------------------------------------------------------------------------

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      b(1) = one_alpha * a(1)
   else if (pass == 2) then
      b(1) = a(1) / (1.0 + alpha)
   else
      b(1) = one_alpha * (a(1) - alpha**3 * a(2)) / (1.0 - alpha**2)**2
   end if

   ! [2.2] Perform pass left to right:

   do j = 2, n
      b(j) = alpha * b(j-1) + one_alpha * a(j)
   end do

   !-------------------------------------------------------------------------
   ! [3.0] Perform left-moving filter:
   !-------------------------------------------------------------------------

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      c(n) = b(n) / (1.0 + alpha)
   else
      c(n) = one_alpha * (b(n) - alpha**3 * b(n-1)) / (1.0 - alpha**2)**2
   end if

   ! [3.2] Perform pass left to right:

   do j = n-1, 1, -1
      c(j) = alpha * c(j+1) + one_alpha * b(j)
   end do
        
   field(1:n) = c(1:n)

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d")
   
end subroutine da_recursive_filter_1d


subroutine da_recursive_filter_vec_1d(npass, alpha, field2d, n1, n2)

   !----------------------------------------------------------------------------------
   ! Purpose: Perform npass passes of recursive filter across 1D segments of 2D array
   !
   ! Method:  Perform right-moving filter followed by left-moving filter.
   !----------------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: npass          ! Number of filter passes.
   real*8,  intent(in)    :: alpha          ! Alpha coefficient for RF.
   real*8,  intent(inout) :: field2d(n1,n2)     ! Array to be filtered.
   integer, intent(in)    :: n1,n2              ! Size of field2d array.

   integer :: i1, i2, pass   ! Loop counters.
   real    :: one_alpha      ! 1 - alpha.
   real    :: field_lr(n1)   ! Field after left-right pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_vec_1d")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   do pass = 1, npass
      do i2 = 1, n2

         !-------------------------------------------------------------------------
         ! [2.0] Perform right-moving filter:
         !-------------------------------------------------------------------------

         ! use turning conditions as in the appendix of Hayden & Purser (1995):

         if (pass == 1) then
            field_lr(1) = one_alpha * field2d(1,i2)
         else if (pass == 2) then
            field_lr(1) = field2d(1,i2) / (1.0 + alpha)
         else
            field_lr(1) = one_alpha * (field2d(1,i2) - alpha**3 * field2d(2,i2)) / (1.0 - alpha**2)**2
         end if

         ! [2.2] Perform pass left to right:

         do i1 = 2, n1
            field_lr(i1) = alpha * field_lr(i1-1) + one_alpha * field2d(i1,i2)
         end do

         !-------------------------------------------------------------------------
         ! [3.0] Perform left-moving filter:
         !-------------------------------------------------------------------------

         ! use turning conditions as in the appendix of Hayden & Purser (1995):

         if (pass == 1) then
            field2d(n1,i2) = field_lr(n1) / (1.0 + alpha)
         else
            field2d(n1,i2) = one_alpha * (field_lr(n1) - alpha**3 * field_lr(n1-1)) / (1.0 - alpha**2)**2
         end if

         ! [3.2] Perform pass left to right:

         do i1 = n1-1, 1, -1
            field2d(i1,i2) = alpha * field2d(i1+1,i2) + one_alpha * field_lr(i1)
         end do

      end do !i2 loop
   end do !pass loop

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_vec_1d")
   
end subroutine da_recursive_filter_vec_1d

