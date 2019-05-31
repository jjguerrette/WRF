subroutine da_transform_vtox(grid, cv_size, xbx, be, ep, cv, vv, vp)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform x' = Uv. 
   !--------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   type(ep_type),  intent(in)    :: ep   ! Ensemble perturbations.
   real,           intent(in)    :: cv(1:cv_size)   ! control variables.
   type(vp_type),  intent(out)   :: vv   ! grdipt/eof cv (local).
   type(vp_type),  intent(inout) :: vp   ! grdipt/level cv (local).

   character(len=10) :: this_time
   character(len=8) :: this_date

   integer :: dy_i, hr_i, mn_i, sc_i, msc_i, &
              dy_f, hr_f, mn_f, sc_f, msc_f
   real    :: time_sec_i, time_sec_f

   if (trace_use) call da_trace_entry("da_transform_vtox")

   call da_zero_x (grid%xa)
   
   if (.not. use_background_errors) then
      if (trace_use) call da_trace_exit("da_transform_vtox")
      return
   end if

   !----------------------------------------------------------------------
   ! [1.0]: Perform vv = u_h cv transform:
   !----------------------------------------------------------------------

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_i
read(this_time(1:2),'(I2.2)') hr_i
read(this_time(3:4),'(I2.2)') mn_i
read(this_time(5:6),'(I2.2)') sc_i
read(this_time(8:10),'(I3.3)') msc_i
time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0

   if (global) then
      call da_transform_vtovv_global(cv_size, xbx, be, cv, vv)
   else if ( (fg_format == fg_format_wrf_arw_regional .or.   &
              fg_format == fg_format_wrf_nmm_regional) .and. &
              (.not. cv_options == 3) )then
      call da_transform_vtovv(grid, cv_size, be, cv, vv)
   end if

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_f
read(this_time(1:2),'(I2.2)') hr_f
read(this_time(3:4),'(I2.2)') mn_f
read(this_time(5:6),'(I2.2)') sc_f
read(this_time(8:10),'(I3.3)') msc_f
time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0

write(unit=message(1),fmt='(A,F12.4,A)') &
   '   inside da_transform_vtox: da_transform_vtovv duration: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
call da_message(message(1:1))

  
   !----------------------------------------------------------------------
   ! [2.0]: Perform vp = u_v vv transform:
   !----------------------------------------------------------------------

   if (  cv_options == 3 ) then

      call da_apply_be( be, cv, vp, grid)
      call da_transform_bal( vp, be, grid)
   
   else

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_i
read(this_time(1:2),'(I2.2)') hr_i
read(this_time(3:4),'(I2.2)') mn_i
read(this_time(5:6),'(I2.2)') sc_i
read(this_time(8:10),'(I3.3)') msc_i
time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0

   if (vert_corr == vert_corr_2) then      
      call da_vertical_transform(grid, 'u', be, grid%xb % vertical_inner_product, vv, vp)
   else
      vp % v1(its:ite,jts:jte,kts:kte) = vv % v1(its:ite,jts:jte,kts:kte)
      vp % v2(its:ite,jts:jte,kts:kte) = vv % v2(its:ite,jts:jte,kts:kte)
      vp % v3(its:ite,jts:jte,kts:kte) = vv % v3(its:ite,jts:jte,kts:kte)
      vp % v4(its:ite,jts:jte,kts:kte) = vv % v4(its:ite,jts:jte,kts:kte)
      vp % v5(its:ite,jts:jte,kts:kte) = vv % v5(its:ite,jts:jte,kts:kte)

      if ( cloud_cv_options >= 2 ) then
         vp % v6(its:ite,jts:jte,kts:kte)  = vv % v6(its:ite,jts:jte,kts:kte)
         vp % v7(its:ite,jts:jte,kts:kte)  = vv % v7(its:ite,jts:jte,kts:kte)
         vp % v8(its:ite,jts:jte,kts:kte)  = vv % v8(its:ite,jts:jte,kts:kte)
         vp % v9(its:ite,jts:jte,kts:kte)  = vv % v9(its:ite,jts:jte,kts:kte)
         vp % v10(its:ite,jts:jte,kts:kte) = vv % v10(its:ite,jts:jte,kts:kte)
      end if
      if ( use_cv_w ) then
         vp % v11(its:ite,jts:jte,kts:kte) = vv % v11(its:ite,jts:jte,kts:kte)
      end if

      if (be % ne > 0) then
!        vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = vv%alpha(its:ite,jts:jte,kts:kte,1:be%ne)
         vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) =  &
             vv%alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)
      end if
   end if

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_f
read(this_time(1:2),'(I2.2)') hr_f
read(this_time(3:4),'(I2.2)') mn_f
read(this_time(5:6),'(I2.2)') sc_f
read(this_time(8:10),'(I3.3)') msc_f
time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0

write(unit=message(1),fmt='(A,F12.4,A)') &
   '   inside da_transform_vtox: vert_corr duration: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
call da_message(message(1:1))

   !----------------------------------------------------------------------  
   ! [3.0]: Perform x = u_p vp transform::
   !----------------------------------------------------------------------

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_i
read(this_time(1:2),'(I2.2)') hr_i
read(this_time(3:4),'(I2.2)') mn_i
read(this_time(5:6),'(I2.2)') sc_i
read(this_time(8:10),'(I3.3)') msc_i
time_sec_i = real( ( ( hr_i * 60 + mn_i ) * 60 + sc_i) * 1000 + msc_i ) / 1000.0

   call da_transform_vptox(grid, vp, be, ep)

call date_and_time(date=this_date, time=this_time)
read(this_date(7:8),'(I2.2)') dy_f
read(this_time(1:2),'(I2.2)') hr_f
read(this_time(3:4),'(I2.2)') mn_f
read(this_time(5:6),'(I2.2)') sc_f
read(this_time(8:10),'(I3.3)') msc_f
time_sec_f = real( ( ( hr_f * 60 + mn_f ) * 60 + sc_f) * 1000 + msc_f ) / 1000.0

write(unit=message(1),fmt='(A,F12.4,A)') &
   '   inside da_transform_vtox: da_transform_vptox duration: ', time_sec_f - time_sec_i + real(dy_f - dy_i) * 86400.0, ' sec.'
call da_message(message(1:1))

   end if
   
   if (trace_use) call da_trace_exit("da_transform_vtox")

end subroutine da_transform_vtox

