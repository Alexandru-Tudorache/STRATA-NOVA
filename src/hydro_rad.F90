subroutine hydro_rad

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,k,location_max,rad_attempt
  integer :: keytemp,keyerr
  real*8 :: denom, applied_step, shell_volume, rho_candidate
  real*8 :: dp_matter, dp_radiation, lambda_momentum
  real*8 :: mechanical_work(imax), power_face(imax), adv_flux(imax), vbar(imax)
  real*8 :: bigr, eddington, eface

  external :: dgtsv
  integer :: info
  real*8 :: dl(imax-2), d(imax-1), du(imax-2), b(imax-1)

  real*8 :: delta_max, last_delta_max
  real*8 :: rad_resid_signed, rad_resid_l1, rad_energy_scale
  real*8 :: rad_resid_signed_rel, rad_resid_l1_rel
  real*8 :: rad_abs_energy_tol
  logical :: rad_signed_converged

  real*8 :: Aarray(imax-1), Barray(imax-1), Carray(imax-1), Darray(imax-1)
  real*8 :: p_temp(imax), eps_temp(imax), lum_temp(imax), temp_temp(imax)
  real*8 :: lambda_temp(imax), lum_old(imax)
  real*8 :: zone_energy_residual, divL_old, divL_new

  real*8, parameter :: RHO_FLOOR        = 1.0d-30
  real*8, parameter :: MAX_REL_TEMP_STEP = 2.5d-1
  real*8, parameter :: RAD_RESID_SIGNED_TOL = 1.0d-6
  real*8, parameter :: RAD_RESID_L1_TOL = 1.0d-4
  real*8, parameter :: RAD_RESID_ABS_TOTAL_FRAC = 2.0d-8


  diag_rad_energy_residual_erg = 0.0d0
  diag_rad_energy_residual_max = 0.0d0
  diag_temp_floor_zones = 0
  diag_radius_clamps = 0
  diag_velocity_caps = 0
  lum_old(:) = lum(:)

  theta = 0.5d0
  rad_attempt = 1

  ! save the old state before moving the grid
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)
  p_p(:) = p(:)
  temp_p(:) = temp(:)
  kappa_p(:) = kappa(:)

  if(do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
        vel(1) = piston_vel
        vel(2) = piston_vel
  endif

  ! update nodal momentum before solving the energy equation
  do i=2,imax
    dp_matter = (p(i)-p_rad(i)) - (p(i-1)-p_rad(i-1))
    dp_radiation = p_rad(i)-p_rad(i-1)
    lambda_momentum = max(0.0d0,min(1.0d0,lambda(i)))

    vel(i) = vel_p(i) &

        - dtime * ggrav*mass(i) / r(i)**2 *gravity_switch   &

        - dtime * 4.0d0*pi*r(i)**2 * &
          (dp_matter + lambda_momentum*dp_radiation) / delta_cmass(i-1)   &

        - dtime * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)

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

  ! Sum the gas+radiation mixed-frame balances (Krumholz et al. 2007,
  ! eqs. 40-42). Besides diffusion, radiation transports g*E*v relative
  ! to a material cell, g=(1-f)/2, and exchanges the work of its force.
  ! The equilibrium one-temperature and non-relativistic approximations remain.
  ! Work from the same discrete forces and kick as the momentum update.
  ! Midpoint velocity makes impulse times velocity equal the nodal kinetic
  ! energy change, including artificial viscosity. No global energy rescale.
  vbar=0.5d0*(vel_p+vel)
  power_face=0.0d0
  adv_flux=0.0d0
  do i=2,imax
     bigr=8.0d0*pi*r(i)**2*abs(temp(i)**4-temp(i-1)**4)*inv_kappa(i) / &
          (delta_cmass(i-1)*max(temp(i)**4+temp(i-1)**4,1.0d-99))
     eddington=lambda(i)/3.0d0+(lambda(i)*bigr/3.0d0)**2
     eddington=max(1.0d0/3.0d0,min(1.0d0,eddington))
     eface=1.5d0*(p_rad(i)+p_rad(i-1))
     adv_flux(i)=4.0d0*pi*r(i)**2*vbar(i)*0.5d0*(1.0d0-eddington)*eface
     power_face(i)=4.0d0*pi*r(i)**2*vbar(i)*lambda(i)*(p_rad(i)-p_rad(i-1))
  enddo
  mechanical_work=0.0d0
  do i=1,imax-1
     mechanical_work(i)=(-dtime*4.0d0*pi*( &
          (p(i)-p_rad(i))*(r(i+1)**2*vbar(i+1)-r(i)**2*vbar(i)) + &
          Q(i)*cr(i)**2*(vbar(i+1)-vbar(i))) + &
          dtime*0.5d0*(power_face(i)+power_face(i+1)) - &
          dtime*(adv_flux(i+1)-adv_flux(i)))/delta_mass(i)
  enddo
  mechanical_work(imax-1)=mechanical_work(imax-1)+ &
       0.5d0*dtime*power_face(imax)/delta_mass(imax-1)
  radiation_work_boundary_power=adv_flux(imax)-adv_flux(1)

  ! reject invalid or crossed radial boundaries
  do i=1,imax
   r(i) = r_p(i) + dtime * vbar(i)
   if (r(i).ne.r(i) .or. r(i).le.0.0d0) then
       scratch_step=.true.
       diag_scratch_reason=4
       diag_scratch_status=1
       diag_scratch_location=i
       return
   endif
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

  do i=1,imax-1
     ! derive density from the unchanged cell mass and new volume
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


  do i=1,imax-1
    cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo
  cr(imax) = r(imax) + (r(imax) - cr(imax-1))


  call artificial_viscosity


  ! refresh radioactive deposition on its own schedule
  if(time.ge.time_Ni) then
      if (time - time_Ni_last_fired .gt. &
           0.5d0 * grid_update_interval_days * 86400.0d0) then
         call nickel_heating
         time_Ni_last_fired = time
         diag_ni_eval_total = diag_ni_eval_total + 1
      endif
      time_Ni = time + Ni_period
  endif

  ! apply bomb heating only during the requested interval
  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
      call bomb_pattern
  else
      bomb_heating(:) = 0.0d0
  endif

200 continue
  p_temp(1:imax) = p(1:imax)
  eps_temp(1:imax) = eps(1:imax)
  temp_temp(1:imax) = temp(1:imax)
  last_delta_max = huge(1.0d0)
  location_max = 0

  do k=1, itmax_rad

    keytemp = 1
    call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), &
         abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
         cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), entropy(1:imax-1), &
         p_rad(1:imax-1),keyerr,keytemp,eoskey)


    call luminosity(r(:), temp_temp(:), kappa_p(:), &
                    lambda_temp(:), inv_kappa(:), lum_temp(:))

    ! assemble the tridiagonal temperature correction
    call matrix_arrays(temp_temp(:), lambda_temp(:), inv_kappa(:), &
        eps_temp(:), lum_temp(:), mechanical_work, &
        Aarray(:), Barray(:), Carray(:), Darray(:))

    ! require energy residuals and temperature corrections to converge
    rad_resid_signed = 0.0d0
    rad_resid_l1 = 0.0d0
    rad_energy_scale = 0.0d0
    do i=1,imax-1
       rad_resid_signed = rad_resid_signed + Darray(i)*delta_mass(i)
       rad_resid_l1 = rad_resid_l1 + abs(Darray(i))*delta_mass(i)
       rad_energy_scale = rad_energy_scale + &
            max(max(abs(eps_p(i)),abs(eps_temp(i))),1.0d12)*delta_mass(i)
    enddo
    rad_resid_signed_rel = abs(rad_resid_signed)/max(rad_energy_scale,1.0d0)
    rad_resid_l1_rel = rad_resid_l1/max(rad_energy_scale,1.0d0)
    rad_abs_energy_tol = RAD_RESID_ABS_TOTAL_FRAC * &
         max(abs(total_initial_energy),rad_energy_scale,1.0d0)
    rad_signed_converged = (rad_resid_signed_rel.le.RAD_RESID_SIGNED_TOL) .or. &
         (abs(rad_resid_signed).le.rad_abs_energy_tol)
    if (k.gt.1 .and. last_delta_max.le.epstol_rad .and. &
         rad_signed_converged .and. rad_resid_l1_rel.le.RAD_RESID_L1_TOL) then
       diag_hydro_rad_iters = k-1
       goto 101
    endif

    if (rad_attempt.eq.1 .and. k.ge.8 .and. last_delta_max.le.epstol_rad .and. &
         ((.not.rad_signed_converged) .or. rad_resid_l1_rel.gt.RAD_RESID_L1_TOL)) then
       goto 150
    endif

    du(1:imax-2) = Aarray(1:imax-2)
    d(1:imax-1)  = Barray(1:imax-1)
    dl(1:imax-2) = Carray(2:imax-1)

    b(1:imax-1) = Darray(1:imax-1)

    info = 0
    ! solve for the coupled temperature corrections
    call dgtsv(imax-1,1,dl,d,du,b,imax-1,info)

    if(info.ne.0) then
       open(unit=666, &
           file=trim(adjustl(trim(adjustl(outdir))//"/failed_matrix.dat")), &
           status="unknown",form='formatted',position="append")
       do i=1,imax-1
           write(666,*) Aarray(i), Barray(i), Carray(i), Darray(i)
       enddo
       close(666)
       error stop "problem in the matrix inversion (see Data/failed_matrix.dat)"
    endif

    delta_max = 0.0d0
    do i=1,imax-1
        denom = max(temp_temp(i), evolution_temperature_floor_k)
        applied_step = b(i)
        if(abs(applied_step).gt.MAX_REL_TEMP_STEP * denom) then
            applied_step = sign(MAX_REL_TEMP_STEP * denom, applied_step)
        endif
        if(temp_temp(i) + applied_step .lt. evolution_temperature_floor_k) then
            applied_step = evolution_temperature_floor_k - temp_temp(i)
        endif
        b(i) = applied_step
        if(abs(applied_step/denom).gt.delta_max) then
            delta_max = abs(applied_step/denom)
            location_max = i
        endif
    enddo
    do i=1, imax-1
        temp_temp(i) = temp_temp(i) + b(i)
        if(temp_temp(i).lt.evolution_temperature_floor_k) then
            temp_temp(i) = evolution_temperature_floor_k
        end if
    end do
    last_delta_max = delta_max

  enddo

150 continue
  ! retry failed crank-nicolson diffusion once with backward euler
  if (rad_attempt.eq.1) then
     rad_attempt = 2
     theta = 1.0d0
     open(unit=995,file=trim(adjustl(outdir))//'/radiation_stiff_fallback.dat', &
          status='unknown',position='append')
     write(995,'(I10,1X,I6,1X,4ES16.7)') nt,k,time,dtime,rad_resid_signed_rel,rad_resid_l1_rel
     close(995)
     goto 200
  endif

  diag_hydro_rad_iters = itmax_rad
  diag_scratch_reason = 1
  diag_scratch_status = 0
  if (last_delta_max.gt.epstol_rad) diag_scratch_status = diag_scratch_status + 1
  if (.not.rad_signed_converged) diag_scratch_status = diag_scratch_status + 2
  if (rad_resid_l1_rel.gt.RAD_RESID_L1_TOL) diag_scratch_status = diag_scratch_status + 4
  diag_scratch_location = location_max
  open(unit=997,file=trim(adjustl(outdir))//'/hydro_rad_failure.dat', &
       status='unknown',position='append')
  write(997,'(I10,1X,7ES16.7,I8)') nt,time,dtime,last_delta_max, &
       rad_resid_signed_rel,rad_resid_l1_rel,theta,rad_resid_signed,location_max
  close(997)
  scratch_step = .true.

  101 continue


  eps(1:imax-1) = eps_temp(1:imax-1)
  p(1:imax-1)   = p_temp(1:imax-1)
  temp(1:imax-1)  = temp_temp(1:imax-1)


  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  p(imax) = 0.0d0


  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))
  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))

  ! measure the accepted discrete energy residual
  diag_temp_floor_zones = count(temp(1:imax-1) .le. &
       evolution_temperature_floor_k*(1.0d0+1.0d-10))
  do i=1,imax-1
     divL_old = (lum_old(i+1)-lum_old(i))/delta_mass(i)
     divL_new = (lum(i+1)-lum(i))/delta_mass(i)
     zone_energy_residual = eps(i)-eps_p(i)-mechanical_work(i) + &
          dtime*((1.0d0-theta)*divL_old + theta*divL_new) - &
          dtime*(bomb_heating(i)+Ni_heating(i))
     diag_rad_energy_residual_erg = diag_rad_energy_residual_erg + &
          zone_energy_residual*delta_mass(i)
     diag_rad_energy_residual_max = max(diag_rad_energy_residual_max,abs(zone_energy_residual))
  enddo

end subroutine hydro_rad
