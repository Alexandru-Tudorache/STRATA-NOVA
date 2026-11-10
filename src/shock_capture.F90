subroutine shock_capture

  use blmod, only: shockpos_stop, shockpos, shockpos_prev, rho, mass, vel, &
            tau, r, bomb_total_energy, time, time_prev, radshock, &
            radshock_prev, breakoutflag
  use parameters
  use physical_constants
  use csm_switch, only: csm_generated_active, csm_generated_surface_radius_cm
  implicit none

  integer :: i, trace_idx

  real*8 :: veldiff(imax)
  real*8 :: veldiffmin
  real*8 :: rhoshock, mshock, taushock
  real*8, save :: velshock = 1.0d0
  real*8 :: velshock_analyt, velbreakout_safe, trace_radius

  character(len=1024) :: filename


  call get_trace_surface(trace_radius, trace_idx)
  if (breakoutflag.eq.0 .and. time_prev.le.0.0d0 .and. &
       radshock_prev.le.0.0d0 .and. already_expanded(trace_radius, trace_idx)) then
    shockpos = trace_idx
    radshock = r(shockpos)
    taushock = tau(shockpos)
    velshock = 0.0d0
    velbreakout_safe = safe_breakout_velocity(velshock, &
         expanded_trace_speed(trace_radius, trace_idx), shockpos)
    call finish_shock_trace('expanded_profile_already_broken_out', &
         velshock, velbreakout_safe, taushock)
    return
  endif

  veldiffmin = 1.0d0
  do i=shockpos_prev,imax-1
    veldiff(i) = (vel(i+1) - vel(i))
    if(veldiff(i).lt.veldiffmin) then
        veldiffmin = veldiff(i)
        shockpos = i
    endif
  enddo

  if(veldiffmin.eq.1.0d0) then
    shockpos = imax
    radshock = r(shockpos)
    taushock = tau(shockpos)

    if(time.gt.time_prev) then
      velshock = (radshock - radshock_prev)/(time-time_prev)
    endif
    velbreakout_safe = safe_breakout_velocity(velshock, 0.0d0, shockpos)

    if(breakoutflag.eq.0) then
      if(radshock_prev.gt.0.0d0 .and. time_prev.gt.0.0d0 .and. &
           breakout_condition(taushock, velbreakout_safe)) then
        call finish_shock_trace('optical_depth_breakout', velshock, &
             velbreakout_safe, taushock)

      else if (shockpos .ge. imax) then

        call finish_shock_trace('outer_boundary_trace_stop', velshock, &
             velbreakout_safe, taushock)

      endif
    endif

    return
  endif

  rhoshock = rho(shockpos)
  mshock = mass(shockpos)
  radshock = r(shockpos)
  taushock = tau(shockpos)
  !equation (19) from Matzner & McKee, ApJ 510:379-403
  if(rhoshock.gt.0.0d0 .and. radshock.gt.0.0d0 .and. mshock.gt.0.0d0) then
    velshock_analyt = 0.794*(bomb_total_energy/mshock)**(0.5d0) &
        *(mshock/(rhoshock*radshock*radshock*radshock))**0.19d0
  else
    velshock_analyt = 0.0d0
  endif


  if(shockpos.ne.shockpos_prev) then

    if(time.gt.time_prev) then
      velshock = (radshock - radshock_prev)/(time-time_prev)
    endif
    time_prev = time
    radshock_prev = radshock

    if(shockpos.ne.imax .and. output_mode.eq.0) then
        filename = trim(adjustl(outdir))//"/velshock_index.dat"
        open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
        write(666,"(I5.4, 5E15.6)") shockpos,time,radshock,velshock, &
                                            velshock_analyt,taushock
        close(666)
    endif

  endif
  velbreakout_safe = safe_breakout_velocity(velshock, velshock_analyt, shockpos)

  if(breakoutflag.eq.0 .and. breakout_condition(taushock, velbreakout_safe) &
       .and. (shockpos.gt.(imax/2) .or. radshock.gt.(0.8d0*r(imax)))) then

    call finish_shock_trace('optical_depth_breakout', velshock, &
         velbreakout_safe, taushock)
    return

  end if

  shockpos_prev = shockpos

contains

  ! classify an already expanding profile without introducing a new explosion
  subroutine get_trace_surface(trace_radius_out, trace_idx_out)
    real*8, intent(out) :: trace_radius_out
    integer, intent(out) :: trace_idx_out

    trace_radius_out = r(imax)
    if (csm_generated_active .and. csm_generated_surface_radius_cm.gt.0.0d0) then
       trace_radius_out = csm_generated_surface_radius_cm
    endif
    trace_idx_out = nearest_radius_index(trace_radius_out)
  end subroutine get_trace_surface

  integer function nearest_radius_index(target_radius) result(idx)
    real*8, intent(in) :: target_radius
    real*8 :: best_dist, dist
    integer :: j

    idx = 1
    best_dist = huge(1.0d0)
    do j = 1, imax
       if (r(j).gt.0.0d0 .and. target_radius.gt.0.0d0) then
          dist = abs(log(r(j) / target_radius))
       else
          dist = abs(r(j) - target_radius)
       endif
       if (dist.lt.best_dist) then
          best_dist = dist
          idx = j
       endif
    enddo
  end function nearest_radius_index

  logical function already_expanded(trace_radius_in, trace_idx_in)
    real*8, intent(in) :: trace_radius_in
    integer, intent(in) :: trace_idx_in
    real*8 :: v_trace, expansion_age

    v_trace = expanded_trace_speed(trace_radius_in, trace_idx_in)

    already_expanded = .false.
    if (trace_radius_in.le.1.0d10) return
    if (.not.usable_speed(v_trace)) return
    if (v_trace.le.1.0d7) return
    expansion_age = trace_radius_in / v_trace
    if (expansion_age.gt.1.0d0 .and. expansion_age.lt.1.0d8) then
       already_expanded = .true.
    endif
  end function already_expanded

  real*8 function expanded_trace_speed(trace_radius_in, trace_idx_in) result(v_trace)
    real*8, intent(in) :: trace_radius_in
    integer, intent(in) :: trace_idx_in
    integer :: j, j_lo, j_hi

    v_trace = 0.0d0
    j_hi = max(1,min(imax,trace_idx_in))
    j_lo = max(1,j_hi-4)
    do j=j_lo,j_hi
       if (r(j).le.1.01d0*trace_radius_in .and. usable_speed(abs(vel(j)))) then
          v_trace = max(v_trace,abs(vel(j)))
       endif
    enddo
  end function expanded_trace_speed

  logical function usable_speed(speed)
    real*8, intent(in) :: speed
    usable_speed = speed.eq.speed .and. speed.gt.0.0d0 .and. &
         abs(speed).lt.huge(1.0d0)*1.0d-6
  end function usable_speed

  real*8 function safe_breakout_velocity(raw_speed, fallback_speed, idx) result(speed)
    real*8, intent(in) :: raw_speed, fallback_speed
    integer, intent(in) :: idx
    integer :: safe_idx

    safe_idx = max(1, min(imax, idx))
    if (usable_speed(raw_speed) .and. raw_speed.lt.0.95d0*clite) then
       speed = raw_speed
    else
       speed = max(abs(fallback_speed), abs(vel(safe_idx)))
       if (.not.usable_speed(speed)) then
          speed = 0.0d0
       else
          speed = min(speed, 0.95d0*clite)
       endif
    endif
  end function safe_breakout_velocity

  logical function breakout_condition(tau_value, speed)
    real*8, intent(in) :: tau_value, speed
    breakout_condition = usable_speed(speed) .and. tau_value.ge.0.0d0 .and. &
         tau_value.lt.(clite / speed)
  end function breakout_condition

  subroutine finish_shock_trace(reason, raw_speed, safe_speed, tau_value)
    character(len=*), intent(in) :: reason
    real*8, intent(in) :: raw_speed, safe_speed, tau_value

    if (trim(reason).eq.'optical_depth_breakout') then
       write(*,*) "Time of breakout is", time
    else if (trim(reason).eq.'expanded_profile_already_broken_out') then
       write(*,*) "Expanded profile is already post-breakout! stopping shock trace"
    else
       write(*,*) "Shock trace reached end of outer boundary at", time
    endif

    open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
         status="unknown",form='formatted',position="append")
    write(666,*) 'Shock trace stop reason = ', trim(reason)
    if (trim(reason).eq.'expanded_profile_already_broken_out') then
       write(666,*) 'Time of breakout = before simulation start'
    else if (trim(reason).eq.'optical_depth_breakout') then
       write(666,*) 'Time of breakout = ', time, 'seconds'
    else
       write(666,*) 'Shock trace stopped = ', time, 'seconds'
    endif
    write(666,*) 'Gridpoint of trace stop = ', shockpos
    write(666,*) 'Shock radius at trace stop = ', radshock, 'cm'
    write(666,*) 'Optical depth at trace stop = ', tau_value
    write(666,*) 'Raw shock speed = ', raw_speed, 'cm/s'
    write(666,*) 'Safe breakout speed = ', safe_speed, 'cm/s'
    close(666)

    breakoutflag = 1
    shockpos_stop = 1
    radshock_prev = radshock
    time_prev = time
    shockpos_prev = shockpos
    write(*,*) 'Finished tracing the shock position!'
  end subroutine finish_shock_trace

end subroutine shock_capture
