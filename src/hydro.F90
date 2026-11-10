subroutine hydro

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,ie,j,m
  integer :: keytemp,keyerr
  real*8 :: dtv
  logical, save :: eos_debug_dumped = .false.
  real*8 :: denom, shell_volume, rho_candidate

  real*8 :: function_star
  real*8 :: function_deriv
  real*8 :: delta_temp(imax)

  real*8 :: delta_max
  integer :: location_max

  real*8 :: p_temp(imax), eps_temp(imax), temp_temp(imax)

  real*8, parameter :: RHO_FLOOR        = 1.0d-30
  real*8, parameter :: MAX_REL_TEMP_STEP = 2.5d-1


  ! follow dt convention of MB93:
  ! use the mean of the current and previous timestep for nodal velocity
  dtv = 0.5d0 * (dtime + dtime_p)
  ! save the old state for this trial update
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)

  if(do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
        vel(1) = piston_vel
        vel(2) = piston_vel
  endif

  ! update velocity from gravity, pressure and artificial viscosity
  do i=2,imax
    vel(i) = vel_p(i) &

     - dtv * ggrav*mass(i) / r(i)**2 *gravity_switch   &

     - dtv * 4.0d0*pi*r(i)**2 * (p(i) - p(i-1)) / delta_cmass(i-1)   &

     - dtv * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)

    if (vel(i).ne.vel(i)) then
       scratch_step=.true.
       diag_scratch_reason=3
       diag_scratch_status=1
       diag_scratch_location=i
       return
    endif
  enddo

  if(do_piston.and.time.ge.piston_tend) then
      vel(1) = 0.0d0
  else if(do_bomb) then
      vel(1) = 0.0d0
  endif

  ! move the cell boundaries and reject crossed or invalid radii
  do i=1,imax
   r(i) = r_p(i) + dtime * vel(i)
   if (r(i).ne.r(i) .or. r(i).le.0.0d0) then
       scratch_step=.true.
       diag_scratch_reason=4
       diag_scratch_status=1
       diag_scratch_location=i
       return
   endif
   ! keep this guard separate because fortran may evaluate both sides of .and.
   if (i.gt.1) then
      if (r(i).le.r(i-1)) then
         scratch_step=.true.
         diag_scratch_reason=4
         diag_scratch_status=2
         diag_scratch_location=i
         return
      endif
   endif
  enddo


  ! keep cell mass fixed while updating its spherical volume
  do i=1,imax-1
     shell_volume = (4.0d0*pi/3.0d0) * (r(i+1)**3 - r(i)**3)
     if (shell_volume.le.0.0d0 .or. shell_volume.ne.shell_volume) then
        scratch_step=.true.
        diag_scratch_reason=5
        diag_scratch_status=1
        diag_scratch_location=i
        return
     endif
     rho_candidate = delta_mass(i) / shell_volume
     if (rho_candidate.ne.rho_candidate .or. rho_candidate.lt.RHO_FLOOR) then
        scratch_step=.true.
        diag_scratch_reason=5
        diag_scratch_status=2
        diag_scratch_location=i
        return
     endif
     rho(i) = rho_candidate
  enddo
  rho(imax) = RHO_FLOOR


  ! place cell centres halfway through the spherical volume
  do i=1,imax-1
     cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo
  cr(imax) = r(imax) + (r(imax) - cr(imax-1))


    call artificial_viscosity


  ! refresh radioactive deposition when its scheduled time is reached
  if(time.ge.time_Ni) then
      if (time - time_Ni_last_fired .gt. &
           0.5d0 * grid_update_interval_days * 86400.0d0) then
         call nickel_heating
         time_Ni_last_fired = time
         diag_ni_eval_total = diag_ni_eval_total + 1
      endif
      time_Ni = time + Ni_period
  endif

  ! apply bomb heating only inside its injection window
  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
      call bomb_pattern
  else
      bomb_heating(:) = 0.0d0
  endif

  temp_temp(1:imax) = temp(1:imax)
  delta_temp(1:imax) = temp_temp(1:imax)
  ie=0
  delta_max = 1

  ! iterate on temperature until the energy residual converges
  do while(delta_max > epstol_hydro)

    keytemp = 1
    call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), &
             abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
             cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), &
             entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)

    delta_max = 0

    do i=1,imax-1

        function_star = eps_temp(i) - eps(i) + 0.5d0*(p_temp(i) + p(i)) &
            * (1.0d0/rho(i)-1.0d0/rho_p(i)) - bomb_heating(i)*dtime &
            - Ni_heating(i)*dtime &
            - Qterm(i)

        function_deriv = dedt(i) + 0.5d0 * dpdt(i)*(1.0d0/rho(i)-1.0d0/rho_p(i))
        delta_temp(i) = -function_star/function_deriv

        ! limit each temperature correction and enforce the evolution floor
        denom = max(temp_temp(i), evolution_temperature_floor_k)
        if(abs(delta_temp(i)).gt.MAX_REL_TEMP_STEP * denom) then
            delta_temp(i) = sign(MAX_REL_TEMP_STEP * denom, delta_temp(i))
        endif
        if(temp_temp(i) + delta_temp(i) .lt. evolution_temperature_floor_k) then
            delta_temp(i) = evolution_temperature_floor_k - temp_temp(i)
        endif

        temp_temp(i) = temp_temp(i) + delta_temp(i)

        temp(i) = temp_temp(i)

        if(abs(delta_temp(i)/max(temp(i), evolution_temperature_floor_k)).gt.delta_max) then
            delta_max = abs(delta_temp(i)/max(temp(i), evolution_temperature_floor_k))
            location_max = i
        end if

    enddo

    ie=ie+1

    if(ie.gt.itmax_hydro) then
        write(6,*) "EOS problem", delta_max, location_max
        if(.not.eos_debug_dumped) then
            eos_debug_dumped = .true.
            write(6,*) "EOS DEBUG START (hydro)"
            write(6,*) "EOS DEBUG nt,time,dtime:", nt, time, dtime
            write(6,*) "EOS DEBUG delta_max,location_max:", delta_max, location_max
            if(location_max.ge.1 .and. location_max.le.imax-1) then
                i = location_max
                write(6,*) "EOS DEBUG mass,delta_mass,r,cr:", mass(i), delta_mass(i), r(i), cr(i)
                write(6,*) "EOS DEBUG rho,rho_p:", rho(i), rho_p(i)
                write(6,*) "EOS DEBUG temp,temp_p,temp_temp:", temp(i), temp_p(i), temp_temp(i)
                write(6,*) "EOS DEBUG p,p_p,p_temp:", p(i), p_p(i), p_temp(i)
                write(6,*) "EOS DEBUG eps,eps_p,eps_temp:", eps(i), eps_p(i), eps_temp(i)
                write(6,*) "EOS DEBUG dpdt,dedt,cs2:", dpdt(i), dedt(i), cs2(i)
                write(6,*) "EOS DEBUG Q,Qterm,delta_temp(i):", Q(i), Qterm(i), delta_temp(i)
                write(6,*) "EOS DEBUG ye,abar,free_electron_frac:", ye(i), abar(i), free_electron_frac(i)
                write(6,*) "EOS DEBUG Ni_heating,bomb_heating:", Ni_heating(i), bomb_heating(i)
                write(6,*) "EOS DEBUG local zone stencil (zone,mass,rho,temp,ye,abar):"
                do j=max(1,i-2),min(imax-1,i+2)
                    write(6,*) "EOS DEBUG ZONE", j, mass(j), rho(j), temp_temp(j), ye(j), abar(j)
                enddo
                write(6,*) "EOS DEBUG composition (j,comp(i,j)):"
                do m=1,ncomps
                    write(6,*) "EOS DEBUG COMP", m, comp(i,m)
                enddo
            endif
            write(6,*) "EOS DEBUG END (hydro)"
        endif
        scratch_step = .true.
        exit
    endif

  enddo

  keytemp = 1
  call eos(rho(1:imax-1),temp(1:imax-1),ye(1:imax-1), &
           abar(1:imax-1),p(1:imax-1),eps(1:imax-1), &
           cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), &
           entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)

  ! set the passive exterior state and zero outer pressure
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  p(imax) = 0.0d0

  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))

  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))


end subroutine hydro
