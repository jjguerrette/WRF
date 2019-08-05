subroutine da_recursive_filter_1d_adj(npass, alpha, field2d, n1, n2)

   !--------------------------------------------------------------------------------------------
   ! Purpose: Perform npass passes of recursive filter across 1D segments of 2D array - adjoint.
   !--------------------------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: npass          ! Number of filter passes.
   real*8,  intent(in)    :: alpha          ! Alpha coefficient for RF.
   real*8,  intent(inout) :: field2d(n1,n2) ! Array to be filtered.
   integer, intent(in)    :: n1,n2          ! Size of field2d array.

   integer                :: i1, i2, pass   ! Loop counters.
   real                   :: one_alpha      ! 1 - alpha.
   real                   :: a(n1)          ! Input field.
   real                   :: b(n1)          ! Field after left-right pass.
   real                   :: c(n1)          ! Field after right-left pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d_adj")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   do pass = npass, 1, -1

      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( i1, i2, a, b, c )
      do i2 = 1, n2

         c = field2d(:,i2)

         !-------------------------------------------------------------------------
         ! [3.0] Perform left-moving filter:
         !-------------------------------------------------------------------------

         ! [3.2] Perform pass right to left (adjoint):

         b = 0.0   

         do i1 = 1, n1-1
            c(i1+1) = c(i1+1) + alpha * c(i1)
            b(i1) = one_alpha * c(i1)
         end do

         ! use turning conditions as in the appendix of Hayden & Purser (1995):

         if (pass == 1) then
            b(n1) = b(n1) + c(n1) / (1.0 + alpha)
         else
            b(n1) = b(n1) + one_alpha * c(n1) / (1.0 - alpha**2)**2
            b(n1-1) = b(n1-1) - one_alpha * alpha**3 * c(n1) / (1.0 - alpha**2)**2
         end if

         !-------------------------------------------------------------------------
         ! [2.0] Perform right-moving filter:
         !-------------------------------------------------------------------------

         a(:) = 0.0

         ! [2.2] Perform pass left to right (adjoint):

         do i1 = n1, 2, -1
            b(i1-1) = b(i1-1) + alpha * b(i1)
            a(i1) = a(i1) + one_alpha * b(i1)
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

         field2d(:,i2) = a

      end do !i2 loop
      !$OMP END PARALLEL DO

   end do !pass loop


   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d_adj")
   
end subroutine da_recursive_filter_1d_adj

