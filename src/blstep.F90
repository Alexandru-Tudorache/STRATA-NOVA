subroutine blstep

  use blmod, only: rho, temp, ye, abar, eps, p, cs2, vel, r, cr, scratch_step, &
                    Q, Qterm, kappa, kappa_table, dkappadt, lambda, inv_kappa, lum, &
                    p_rad, dpdt, dedt, entropy, Ni_heating, Ni_deposit_function, &
                    bomb_heating, Ni_total_luminosity, Ni_energy_rate, &
                    time, dtime, shockpos_stop, nt, tdump_scalar, &
                    expanded_input_detected, &
                    time_Ni, time_Ni_last_fired, &
                    diag_scratch_count, diag_scratch_reason, &
                    diag_scratch_status, diag_scratch_location, diag_analysis_called
  use parameters
  use physical_constants
  implicit none
  logical :: scalar_analysis_needed

  ! save variables
  real*8 :: rho_save(imax),temp_save(imax),vel_save(imax)
  real*8 :: p_save(imax),ye_save(imax),cs2_save(imax), abar_save(imax)
  real*8 :: r_save(imax), cr_save(imax), eps_save(imax)
  real*8 :: Q_save(imax), Qterm_save(imax)
  real*8 :: kappa_save(imax), kappa_table_save(imax), dkappadt_save(imax)
  real*8 :: lambda_save(imax), inv_kappa_save(imax), lum_save(imax)
  real*8 :: p_rad_save(imax), dpdt_save(imax), dedt_save(imax), entropy_save(imax)
  real*8 :: Ni_heating_save(imax), Ni_deposit_function_save(imax)
  real*8 :: bomb_heating_save(imax)
  real*8 :: time_Ni_save, time_Ni_last_fired_save
  real*8 :: Ni_total_luminosity_save, Ni_energy_rate_save

  integer :: iterations
  logical :: scratch_file_exists
  character(len=512) :: scratch_diag_path

!------------------------------------------------------------------------------

  ! save old values in case we have to redo the step
  rho_save   = rho
  temp_save  = temp
  ye_save    = ye
  abar_save    = abar
  eps_save    = eps
  p_save     = p
  cs2_save   = cs2
  vel_save  = vel
  r_save     = r
  cr_save    = cr
  Q_save = Q
  Qterm_save = Qterm
  kappa_save = kappa
  kappa_table_save = kappa_table
  dkappadt_save = dkappadt
  lambda_save = lambda
  inv_kappa_save = inv_kappa
  lum_save = lum
  p_rad_save = p_rad
  dpdt_save = dpdt
  dedt_save = dedt
  entropy_save = entropy
  Ni_heating_save = Ni_heating
  Ni_deposit_function_save = Ni_deposit_function
  bomb_heating_save = bomb_heating
  time_Ni_save = time_Ni
  time_Ni_last_fired_save = time_Ni_last_fired
  Ni_total_luminosity_save = Ni_total_luminosity
  Ni_energy_rate_save = Ni_energy_rate

  iterations = 0
  diag_scratch_count = 0
  diag_scratch_reason = 0
  diag_scratch_status = 0
  diag_scratch_location = 0
  diag_analysis_called = .false.

  do while(scratch_step.or.(iterations.eq.0))

     iterations = iterations + 1

     if (scratch_step) then
        ! redo step with half the timestep; this is 
        ! sometimes necessary if convergence failes
        if (iterations.gt.10) then
           write(6,*) '******* Time step failed after 10 scratch retries.'
           write(6,*) '******* reason/status/zone = ',diag_scratch_reason, &
                diag_scratch_status, diag_scratch_location
           error stop 1
        endif
        write(6,*) "Scratching entire step...", nt
        write(6,"(A25,1P10E15.6)") "time,dt,dt_new: ",time, dtime, dtime/2.0d0
        scratch_diag_path=trim(adjustl(outdir))//'/scratch_diagnostics.dat'
        inquire(file=trim(scratch_diag_path),exist=scratch_file_exists)
        open(unit=998, file=trim(scratch_diag_path), status='unknown', position='append')
        if (.not.scratch_file_exists) then
           write(998,'(A)') '# nt time_s old_dt new_dt retry reason status zone'
        endif
        write(998,'(I8,3ES16.7,I5,3I6)') nt, time, dtime, &
             dtime/2.0d0, iterations-1, diag_scratch_reason, &
             diag_scratch_status, diag_scratch_location
        close(998)

        dtime = dtime / 2.0d0
        rho = rho_save
        temp = temp_save
        ye = ye_save
        abar = abar_save
        eps = eps_save
        p = p_save
        cs2 = cs2_save
        vel = vel_save
        r = r_save
        cr = cr_save

        ! hydro_rad.F90 updates several derived fields before it checks for convergence
        ! older 'saved' terms MUST be fed into the next evolution loop if a failure occurs
        ! if not, this can poison the next step and destroy the output's physicality

        Q = Q_save
        Qterm = Qterm_save
        kappa = kappa_save
        kappa_table = kappa_table_save
        dkappadt = dkappadt_save
        lambda = lambda_save
        inv_kappa = inv_kappa_save
        lum = lum_save
        p_rad = p_rad_save
        dpdt = dpdt_save
        dedt = dedt_save
        entropy = entropy_save
        Ni_heating = Ni_heating_save
        Ni_deposit_function = Ni_deposit_function_save
        bomb_heating = bomb_heating_save
        time_Ni = time_Ni_save
        time_Ni_last_fired = time_Ni_last_fired_save
        Ni_total_luminosity = Ni_total_luminosity_save
        Ni_energy_rate = Ni_energy_rate_save
        diag_scratch_count = diag_scratch_count + 1
        scratch_step = .false.

     endif

     ! select between pure hydro or
     ! radiation+hydro solver

     if(radiation) then
        call hydro_rad
     else
        call hydro
     end if

     ! restore state at start of loop 
     if (scratch_step) cycle

     ! expansion checker: some ejecta supplied by the user can already be expanded
     ! it doesn't make sense for the code to set up a piston once again
     ! this checks for an expanded profile and calls shock_capture to stop tracing the shock
     ! lets ejecta homologously expand after detection
     if ((time.ge.bomb_tend .or. expanded_input_detected) .and. &
          shockpos_stop.eq.0) then
        call shock_capture
     endif
     
     ! compute luminosity
     ! time is the "old" time here (same one used in snec.F90 for tdump_scalar)
     if (.not.sedov) then
        scalar_analysis_needed = time.eq.0.0d0 .or. time.ge.tdump_scalar .or. &
             time+dtime.ge.tend .or. nt+1.ge.ntmax
        if (ntout_scalar.gt.0) scalar_analysis_needed = scalar_analysis_needed .or. &
             mod(nt+1,int(ntout_scalar,int64)).eq.0
        if (scalar_analysis_needed) then
           call analysis
           diag_analysis_called = .true.
        endif
     endif

     ! energy conservation: we already did an initial call in problem.F90
     ! here we just sample it on the scalar dump cadence
     if (time.ge.tdump_scalar) then
        call conservation_compute_energies
     endif

  enddo

end subroutine blstep
