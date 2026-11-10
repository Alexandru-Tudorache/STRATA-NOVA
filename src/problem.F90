subroutine problem

  use blmod, only: ntstart, tstart, ncomps, mass, r, Rstar, &
            opacity_floor, metallicity, envelope_metallicity, do_piston, &
            do_bomb, total_initial_energy, bomb_total_energy, bomb_spread, &
            rho, vel, kappa, kappa_table, dkappadt, tau, temp, p, delta_mass, &
            lambda, inv_kappa, lum, eos_gamma1, expanded_input_detected
  use parameters
  use eosmodule
  use physical_constants
  use startup_log, only: print_startup_separator
  use csm_switch, only: csm_generated_active, csm_progenitor_surface_radius_cm
  use staged_state, only: restore_staged_state
  implicit none


  logical :: coast_restored
  integer :: i
  character(len=256) :: filename

  real*4 :: opact,dopact,dopacr,dopactd
  common/e/ opact,dopact,dopacr,dopactd


  ntstart       = 0
  tstart        = 0.0d0


  if(eoskey.eq.1) then
      write(6,'(A)') ' EOS: ideal gas'
      write(6,'(A)') ' WARNING: ideal EOS does not return radiation pressure.'
      eos_gamma1 = 1.4d0
  else if(eoskey.eq.2) then
      write(6,'(A)') ' EOS: Paczynski'
  else
     error stop "Choice of EOS not available, check the parameter eoskey."
  endif
  call print_startup_separator


  call get_ncomps_from_profile(composition_profile_name,ncomps)

  call allocate_vars

  call initialize_vars


  call get_inner_outer_mass_from_profile(profile_name,mass(1),mass(imax))

  if(mass_excision) then
      mass(1) = mass_excised*msun
  endif

  call get_inner_outer_radius_from_profile(profile_name,mass(1),r(1),r(imax))

  Rstar=r(imax)

  open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
      status="unknown",form='formatted',position="append")
  write(666,*) 'Mass of the model = ', mass(imax)/msun, 'solar masses'
  write(666,*) 'Initial radius = ', Rstar/rsun, 'solar radii'
  close(666)

  call grid


  write(*,'(A,A)') ' Profile: ',trim(profile_name)
  ! read the initial hydrodynamic and composition state
  call read_profile(profile_name)
  call restore_staged_state(coast_restored)
  call print_startup_separator

  if (coast_restored) then
     Rstar=r(imax)
  else if (csm_generated_active) then
     Rstar=r(imax)
  else
     call integrate_radius_initial
  endif

  expanded_input_detected = detect_expanded_input_profile()


  do i=1, imax
    opacity_floor(i) = (envelope_metallicity*of_core - of_env      &
                + metallicity(i)*(of_env - of_core))/(envelope_metallicity - 1)
  end do

  if (output_mode .eq. 0) then
     filename = trim(adjustl(outdir))//"/opacity_floor.dat"
     call output_screenshot(opacity_floor,filename,imax)
  endif


  ! measure the initial energy before setting an explosion source
  call conservation_compute_energies

  if(initial_data.eq."Piston_Explosion") then

    if (expanded_input_detected) then
       do_piston = .false.
       write(6,*) 'Expanded input profile detected: piston disabled.'
       open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
            status="unknown",form='formatted',position="append")
       write(666,*) 'Piston suppressed = expanded input profile already exploding'
       close(666)
    else
       do_piston = .true.
    endif

  ! choose the total thermal energy to deposit
  else if(initial_data.eq."Thermal_Bomb") then

    do_bomb = .true.

    if(bomb_mode.eq.1) then
       bomb_total_energy = final_energy - total_initial_energy
    else if (bomb_mode.eq.2) then
       bomb_total_energy = final_energy
    else
       write(6,"(A10,I3,A12)") "bomb_mode ",bomb_mode," not defined!"
       stop 1
    endif

    call print_startup_separator
    write(6,'(A,I0)') ' Thermal bomb mode: ',bomb_mode

    if (bomb_start_point.lt.1 .or. bomb_start_point.ge.imax) then
       write(6,*) 'bomb_start_point must be in the active grid range 1..imax-1.'
       write(6,*) 'bomb_start_point = ', bomb_start_point, ' imax = ', imax
       stop 1
    endif
    if (bomb_mass_spread.le.0.0d0) then
       write(6,*) 'bomb_mass_spread must be > 0.'
       stop 1
    endif

    bomb_spread = 0
    do i=bomb_start_point+1, imax-1
        if(mass(i).ge.(mass(bomb_start_point)+bomb_mass_spread*msun)) then
            bomb_spread = i - bomb_start_point
            open(unit=666, &
                file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
                status="unknown",form='formatted',position="append")
            write(666,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
            write(6,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
            close(666)
            exit
        end if
    end do

    if (bomb_spread.le.0) then
       write(6,*) 'Bomb mass window extends beyond the active grid.'
       write(6,*) 'bomb_start_point = ', bomb_start_point
       write(6,*) 'bomb_mass_spread (Msun) = ', bomb_mass_spread
       write(6,*) 'available active mass (Msun) = ', &
            (mass(imax-1) - mass(bomb_start_point)) / msun
       stop 1
    endif


    open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
        status="unknown",form='formatted',position="append")
    write(666,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
    write(666,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
    write(6,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
    write(6,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
    close(666)

  else
    error stop "Wrong type of explosion, check the parameter 'initial_data'"
  endif

  call compose_opacity_tables_OPAL

  call opacity(rho(:),temp(:),kappa(:),kappa_table(:),dkappadt(:))
  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))


  call read_BolCorr


  if (output_mode .eq. 0) then

  filename = trim(adjustl(outdir))//"/rho_initial.dat"
  call output_screenshot(rho,filename,imax)

  filename = trim(adjustl(outdir))//"/rad_initial.dat"
  call output_screenshot(r,filename,imax)

  filename = trim(adjustl(outdir))//"/mass_initial.dat"
  call output_screenshot(mass,filename,imax)

  filename = trim(adjustl(outdir))//"/delta_mass_initial.dat"
  call output_screenshot(delta_mass,filename,imax)

  filename = trim(adjustl(outdir))//"/press_initial.dat"
  call output_screenshot(p,filename,imax)

  filename = trim(adjustl(outdir))//"/density_profile.dat"
  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
          form='formatted',position="append")
  do i=1, imax-1
      write(666,"(1P20E19.10E3)") r(imax)-r(imax-i), &
          rho(imax-i),mass(imax)-mass(imax-i)
  enddo
  close(666)

  endif

contains

  logical function detect_expanded_input_profile() result(is_expanded)
    integer :: j
    real*8 :: ejecta_surface_radius, v_outer, expansion_age

    is_expanded = .false.
    ejecta_surface_radius = Rstar
    if (csm_generated_active .and. csm_progenitor_surface_radius_cm.gt.0.0d0) then
       ejecta_surface_radius = csm_progenitor_surface_radius_cm
    endif

    if (ejecta_surface_radius.le.1.0d10) return
    v_outer = 0.0d0
    do j=1,imax
       if (r(j).le.1.01d0*ejecta_surface_radius) then
          if (vel(j).eq.vel(j)) v_outer = max(v_outer,abs(vel(j)))
       endif
    enddo
    if (v_outer.le.1.0d7) return

    expansion_age = ejecta_surface_radius/v_outer
    if (expansion_age.gt.1.0d0 .and. expansion_age.lt.1.0d8) then
       is_expanded = .true.
    endif
  end function detect_expanded_input_profile

end subroutine problem
