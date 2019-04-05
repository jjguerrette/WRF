subroutine da_get_trh( input_file, i_start, i_end, j_start, j_end, k, temp, rh )

   !---------------------------------------------------------------------------
   ! Purpose: Calculates T, RH from input WRF file.
   !---------------------------------------------------------------------------

   implicit none

   character(len=200), intent(in)  :: input_file       ! NETCDF file name.
   integer,            intent(in)  :: i_start, i_end   ! Dimension 1 of field. 
   integer,            intent(in)  :: j_start, j_end   ! Dimension 2 of field. 
   integer,            intent(in)  :: k                              ! Model level.
   real,               intent(out) :: temp(i_start:i_end,j_start:j_end)          ! Temperature.
   real,               intent(out) :: rh(i_start:i_end,j_start:j_end)            ! Relative humidity.

   character(len=10) :: var                       ! Variable to search for. var = "T"
   integer           :: i, j                      ! Loop counters.

   real              :: thetap(i_start:i_end,j_start:j_end)     ! Perturbation potential temperature.
   real              :: pb(i_start:i_end,j_start:j_end)         ! Base state pressure.
   real              :: pp(i_start:i_end,j_start:j_end)         ! Pressure perturbation.
   real              :: x(i_start:i_end,j_start:j_end)          ! Vapor mixing ratio.

   real              :: theta                     ! Potential temperature.
   real              :: p                         ! Pressure.
   real              :: q                         ! Specific humidity.
   real              :: t_c                       ! Temp(Celsius).
   real              :: es                        ! Saturation vapor pressure.
   real              :: qs                        ! Saturation specific humidity.

   if (trace_use) call da_trace_entry("da_get_trh")

   var = "T" ! Perturbation potential temperature in WRF.
   call da_get_field( input_file, var, 3, i_start, i_end, j_start, j_end, k, thetap)

   var = "PB"  ! Base state pressure in WRF.
   call da_get_field( input_file, var, 3, i_start, i_end, j_start, j_end, k, pb)

   var = "P" ! Perturbation pressure in WRF.
   call da_get_field( input_file, var, 3, i_start, i_end, j_start, j_end, k, pp)

   var = "QVAPOR"  ! Water vapor mixing ratio.
   call da_get_field( input_file, var, 3, i_start, i_end, j_start, j_end, k, x)

   do j = j_start, j_end
      do i = i_start, i_end
         ! Convert p', theta' to T:
         theta = t0 + thetap(i,j)                 ! Theta = Theta0 + Thetap
         p = pb(i,j) + pp(i,j)                     ! p = base p + pert p.
         temp(i,j) = theta *( p/base_pres)**kappa ! Theta to T.

         ! Convert to specific humidity.
         q = x(i,j) /( 1.0 + x(i,j))

         ! Calculate relative humidity:
         t_c = temp(i,j) - t_kelvin
         es = es_alpha * exp( es_beta * t_c /( t_c + es_gamma))
         qs = rd_over_rv * es /( p - rd_over_rv1 * es)
         rh(i,j) = 100.0 * q / qs
      end do
   end do

   if (trace_use) call da_trace_exit("da_get_trh")

end subroutine da_get_trh


