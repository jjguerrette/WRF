subroutine da_get_height( input_file, i_end, j_end, k_end, height)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates T, RH from input WRF file.
   !---------------------------------------------------------------------------
   
   implicit none
   
   character(len=200), intent(in)  :: input_file       ! NETCDF file nane.
   integer,            intent(in)  :: i_end, j_end, k_end          ! Dimensions.
   real,               intent(out) :: height(1:i_end,1:j_end,1:k_end) ! Height.

   character(len=10) :: var                            ! Variable to search for.
   integer           :: k                              ! Loop counter.
   real              :: gravity_inv                    ! 1/gravity.
   real              :: phb(1:i_end,1:j_end)             ! Base state geopotential.
   real              :: ph(1:i_end,1:j_end)              ! Perturbation geopotential.
   real              :: phfull(1:i_end,1:j_end,1:k_end+1) ! Geopotential.

   if (trace_use) call da_trace_entry("da_get_height")

   gravity_inv = 1.0 / gravity

   do k = 1, k_end+1
      var = "PHB"
      call da_get_field( input_file, var, 3, 1, i_end, 1, j_end, k, phb)
      var = "PH"
      call da_get_field( input_file, var, 3, 1, i_end, 1, j_end, k, ph)

      phfull(:,:,k) = phb + ph ! Calculate full geopotential on full(w) model levels:
   end do

   do k = 1, k_end
      height(:,:,k) = 0.5 *( phfull(:,:,k+1) + phfull(:,:,k)) * gravity_inv
   end do

   if (trace_use) call da_trace_exit("da_get_height")

end subroutine da_get_height


