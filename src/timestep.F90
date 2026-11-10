subroutine timestep

  use blmod, only: cs2, r, vel, dtime, dtime_p, nt, time, delta_time, &
                   shockpos, shockpos_stop, rho, temp, p, p_rad, lambda, &
                   cfl_limit_zone, cfl_raw_min_dt, mass, delta_cmass, Q, cr, &
                   do_piston, do_bomb, gravity_switch, expanded_input_detected
  use parameters
  use physical_constants
  implicit none

  integer :: i
  integer, parameter :: LEGACY_SHOCK_WINDOW = 10
  integer, parameter :: CFL_RAMP_STEPS = 200

  real*8 :: sound
  real*8 :: signal, eul_signal, lag_signal
  real*8 :: compr
  real*8 :: dttrans
  real*8 :: raw_min_dt, after_cfl_dt, before_bounds_dt
  real*8 :: limit_dr, limit_sound, limit_signal, limit_lag_signal
  real*8 :: limit_eul_signal, limit_vel_i, limit_vel_ip1
  real*8 :: limit_rho, limit_temp, limit_p_matter, limit_p_rad_eq
  real*8 :: limit_cs2
  integer :: limit_i
  logical :: growth_limited, dttrans_limited, dtmin_limited, dtmax_limited
  character(len=16) :: current_mode, limit_mode

  integer(int64), save :: nt_shock_stopped = -1
  integer, save :: shockpos_at_stop = -1
  integer, save :: cfl_diag_unit = -1
  logical, save :: cfl_diag_disabled = .false.
  real*8 :: alpha_eul, ramp_xi
  integer :: anchor, shock_window_active
  logical :: in_ramp_window, use_eulerian_here

  real*8, parameter :: CSM_GEOM_MIN_WIDTH_FRAC = 5.0d-1
  real*8 :: accel_pred(imax), vel_trial(imax)
  real*8 :: geom_candidate_dt, geom_lo, geom_hi, geom_mid, geom_dtv
  real*8 :: accel_dt, accel_scale, accel_dp_matter, accel_dp_rad, accel_lambda
  real*8 :: geom_ratio_min, geom_dr, geom_dr_trial
  integer :: geom_iter, geom_limit_zone
  logical :: geom_active, geom_valid


  shock_window_active = LEGACY_SHOCK_WINDOW

  ! start from the configured upper timestep bound
  dttrans = dtmax

  if(nt.le.1 .and. .not.(expanded_input_detected .and. shockpos_stop.eq.1)) then
    dttrans = 1.0d-8
    dtime = 1.0d-8
  endif

  dtime_p = dtime

  if (shockpos_stop .eq. 0 .and. nt_shock_stopped .ge. 0) then
     nt_shock_stopped = -1
     shockpos_at_stop = -1
  endif

  if (shockpos_stop .eq. 1 .and. nt_shock_stopped .lt. 0) then
     nt_shock_stopped = nt
     shockpos_at_stop = shockpos
  endif

  if (shockpos_stop .eq. 0) then
     alpha_eul = 1.0d0
     anchor = shockpos
     in_ramp_window = (shockpos .gt. 0) .and. (shockpos .lt. imax)
  else if (expanded_input_detected) then
     alpha_eul = 0.0d0
     anchor = 0
     in_ramp_window = .false.
  else if (nt - nt_shock_stopped .le. CFL_RAMP_STEPS) then
     ramp_xi = dble(nt - nt_shock_stopped) / dble(CFL_RAMP_STEPS)
     alpha_eul = 0.5d0 * (1.0d0 + cos(pi * ramp_xi))
     anchor = shockpos_at_stop
     in_ramp_window = (anchor .gt. 0) .and. (anchor .lt. imax)
  else
     alpha_eul = 0.0d0
     anchor = 0
     in_ramp_window = .false.
  endif

  dtime = dtmax
  raw_min_dt = huge(1.0d0)
  limit_i = 0
  limit_dr = 0.0d0
  limit_sound = 0.0d0
  limit_signal = 0.0d0
  limit_lag_signal = 0.0d0
  limit_eul_signal = 0.0d0
  limit_vel_i = 0.0d0
  limit_vel_ip1 = 0.0d0
  limit_rho = 0.0d0
  limit_temp = 0.0d0
  limit_p_matter = 0.0d0
  limit_p_rad_eq = 0.0d0
  limit_cs2 = 0.0d0
  limit_mode = 'none'
  ! limit each interval by its acoustic and compression crossing time
  do i=1,imax-1
    sound = sqrt(max(cs2(i), 0.0d0))
    compr = max(0.0d0, vel(i) - vel(i+1))
    lag_signal = sound + compr

    use_eulerian_here = in_ramp_window .and. abs(i-anchor).le.shock_window_active
    if (use_eulerian_here) then
      eul_signal = max(abs(vel(i)+sound), abs(vel(i)-sound))
      signal = alpha_eul * eul_signal + (1.0d0 - alpha_eul) * lag_signal
      if (alpha_eul.ge.1.0d0 - 1.0d-12) then
         current_mode = 'eulerian'
      else if (alpha_eul.le.1.0d-12) then
         current_mode = 'lagrangian'
      else
         current_mode = 'blended'
      endif
    else
      eul_signal = 0.0d0
      signal = lag_signal
      current_mode = 'lagrangian'
    endif

    delta_time(i) = (r(i+1) - r(i)) / signal
    dtime = min(dtime, delta_time(i))
    if (delta_time(i).lt.raw_min_dt) then
       raw_min_dt = delta_time(i)
       limit_i = i
       limit_dr = r(i+1) - r(i)
       limit_sound = sound
       limit_signal = signal
       limit_lag_signal = lag_signal
       limit_eul_signal = eul_signal
       limit_vel_i = vel(i)
       limit_vel_ip1 = vel(i+1)
       limit_rho = rho(i)
       limit_temp = temp(i)
       limit_p_rad_eq = p_rad(i)
       limit_p_matter = max(p(i) - p_rad(i), 0.0d0)
       limit_cs2 = cs2(i)
       limit_mode = current_mode
    endif
  enddo

  cfl_limit_zone = limit_i
  cfl_raw_min_dt = raw_min_dt

  ! apply the requested courant factor to the smallest cell limit
  dtime = dt_cfl_factor * dtime
  after_cfl_dt = dtime

  ! limit timestep growth between accepted updates
  growth_limited = (dt_growth_factor*dtime_p.lt.dtime)
  if (expanded_input_detected .and. shockpos_stop.eq.1) then
     growth_limited = .false.
  else
     dtime = min(dt_growth_factor*dtime_p,dtime)
  endif

  dttrans_limited = (dttrans.lt.0.75d0*dtime)
  if(dttrans.lt.0.75d0*dtime) then
    dtime = dttrans
  endif

  geom_active = trim(adjustl(csm)).ne.'none' .and. &
       (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
        trim(adjustl(mass_scheme)).eq.'constant_mass_loss' .or. &
        trim(adjustl(mass_scheme)).eq.'constant_density')
  geom_limit_zone = 0
  ! check that predicted acceleration cannot cross cell boundaries
  if (geom_active .and. dtime.gt.0.0d0) then
     accel_pred(:) = 0.0d0
     do i=2,imax
        if (radiation) then
           if (trim(adjustl(csm)).eq.'active' .and. &
                trim(adjustl(mass_scheme)).eq.'constant_density') then
              accel_dp_matter=p(i)-p(i-1)
              accel_dp_rad=0.0d0
              accel_lambda=0.0d0
           else
              accel_dp_matter=(p(i)-p_rad(i))-(p(i-1)-p_rad(i-1))
              accel_dp_rad=p_rad(i)-p_rad(i-1)
              accel_lambda=max(0.0d0,min(1.0d0,lambda(i)))
           endif
        else
           accel_dp_matter=p(i)-p(i-1)
           accel_dp_rad=0.0d0
           accel_lambda=0.0d0
        endif
        accel_pred(i) = &
             - ggrav*mass(i) / r(i)**2 * gravity_switch &
             - 4.0d0*pi*r(i)**2 * &
               (accel_dp_matter+accel_lambda*accel_dp_rad) / delta_cmass(i-1) &
             - 4.0d0*pi * (cr(i)**2*Q(i)-cr(i-1)**2*Q(i-1)) / delta_cmass(i-1)
     enddo
     if (do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
        accel_pred(1:2) = 0.0d0
     endif

     accel_dt=dtime
     do i=2,imax
        accel_scale=max(1.0d7,max(abs(vel(i)),sqrt(max(cs2(max(1,i-1)),0.0d0))))
        if (abs(accel_pred(i)).gt.1.0d-99) then
           accel_dt=min(accel_dt,hydro_accel_cfl*accel_scale/abs(accel_pred(i)))
        endif
     enddo
     dtime=min(dtime,accel_dt)

     geom_candidate_dt = dtime
     call geometry_trial_ok(geom_candidate_dt, geom_valid, geom_ratio_min, geom_limit_zone)
     if (.not.geom_valid) then
        geom_lo = 0.0d0
        geom_hi = geom_candidate_dt
        do geom_iter=1,40
           geom_mid = 0.5d0*(geom_lo+geom_hi)
           call geometry_trial_ok(geom_mid, geom_valid, geom_ratio_min, geom_limit_zone)
           if (geom_valid) then
              geom_lo = geom_mid
           else
              geom_hi = geom_mid
           endif
        enddo
        dtime = max(0.98d0*geom_lo, 0.0d0)
     endif
  endif

  ! enforce the configured minimum and maximum timestep
  before_bounds_dt = dtime
  dtmin_limited = (before_bounds_dt.lt.dtmin)
  dtmax_limited = (max(before_bounds_dt,dtmin).gt.dtmax)
  dtime = min(max(dtime,dtmin),dtmax)

  if (grid_debug.ge.1) call write_timestep_cfl_diag

contains

  ! predict cell geometry before accepting the trial timestep
  subroutine geometry_trial_ok(dt_try, ok, ratio_min, limit_zone)
    implicit none
    real*8, intent(in) :: dt_try
    logical, intent(out) :: ok
    real*8, intent(out) :: ratio_min
    integer, intent(out) :: limit_zone
    integer :: j

    geom_dtv = 0.5d0*(dt_try+dtime_p)
    vel_trial(:) = vel(:) + geom_dtv*accel_pred(:)
    if (do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
       vel_trial(1) = piston_vel
       vel_trial(2) = piston_vel
    endif
    if (do_piston.and.time.ge.piston_tend) then
       vel_trial(1)=0.0d0
    else if (do_bomb) then
       vel_trial(1)=0.0d0
    endif

    ratio_min = huge(1.0d0)
    limit_zone = 0
    do j=1,imax-1
       geom_dr = r(j+1)-r(j)
       if (geom_dr.le.0.0d0) then
          ratio_min = -huge(1.0d0)
          limit_zone = j
          exit
       endif
       geom_dr_trial = geom_dr + dt_try*(vel_trial(j+1)-vel_trial(j))
       if (geom_dr_trial/geom_dr.lt.ratio_min) then
          ratio_min = geom_dr_trial/geom_dr
          limit_zone = j
       endif
    enddo
    ok = ratio_min.ge.CSM_GEOM_MIN_WIDTH_FRAC
  end subroutine geometry_trial_ok

  subroutine write_timestep_cfl_diag
    implicit none

    integer :: ios
    character(len=512) :: diag_file

    if (cfl_diag_disabled) return

    if (cfl_diag_unit.lt.0) then
       cfl_diag_unit = 221
       diag_file = trim(adjustl(outdir))//'/timestep_cfl_diag.dat'
       open(unit=cfl_diag_unit, file=trim(diag_file), status='replace', &
            action='write', iostat=ios)
       if (ios.ne.0) then
          write(*,*) 'WARNING: could not open timestep CFL diagnostic: ', trim(diag_file)
          cfl_diag_disabled = .true.
          cfl_diag_unit = -1
          return
       endif
       write(cfl_diag_unit,'(A)') &
            '# nt time_s dtime_s raw_min_dt_s after_cfl_dt_s limiting_zone ' // &
            'dr_cm sound_cm_s vel_i_cm_s vel_ip1_cm_s lag_signal eul_signal ' // &
            'chosen_signal rho_g_cm3 temp_K p_matter_dyn_cm2 p_rad_eq_dyn_cm2 ' // &
            'cs2_cm2_s2 mode shockpos shockpos_stop anchor alpha_eul ' // &
            'shock_half_window growth_cap dttrans_cap dtmin_cap dtmax_cap'
    endif

    write(cfl_diag_unit,'(I8,1X,4ES16.7,1X,I6,1X,12ES16.7,1X,A16,1X,3I7,1X,ES12.4,1X,I4,1X,4L2)') &
         nt, time, dtime, raw_min_dt, after_cfl_dt, limit_i, &
         limit_dr, limit_sound, limit_vel_i, limit_vel_ip1, &
         limit_lag_signal, limit_eul_signal, limit_signal, &
         limit_rho, limit_temp, limit_p_matter, limit_p_rad_eq, &
         limit_cs2, &
         trim(adjustl(limit_mode)), shockpos, shockpos_stop, anchor, &
         alpha_eul, shock_window_active, growth_limited, dttrans_limited, dtmin_limited, &
         dtmax_limited
    flush(cfl_diag_unit)

  end subroutine write_timestep_cfl_diag

end subroutine timestep
