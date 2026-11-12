! write the products selected by output mode and cadence
subroutine output_all(modeflag)

  use blmod
  use parameters
  use physical_constants
  implicit none

  character(len=1024) :: filename
  integer :: modeflag


  if(modeflag.eq.0) then


  else if(modeflag.eq.1) then

    if (output_mode .eq. 0 .or. output_mode .eq. 2) then

    filename = trim(adjustl(outdir))//"/vel.xg"
    call output_single_mass(vel,filename)

    filename = trim(adjustl(outdir))//"/rho.xg"
    call output_single_mass(rho,filename)

    if (output_mode .eq. 0) then
       filename = trim(adjustl(outdir))//"/ye.xg"
       call output_single_mass(ye,filename)
    endif

    filename = trim(adjustl(outdir))//"/press.xg"
    call output_single_mass(p,filename)

    filename = trim(adjustl(outdir))//"/cs2.xg"
    call output_single_mass(cs2,filename)

    filename = trim(adjustl(outdir))//"/Q.xg"
    call output_single_mass(Q,filename)

    filename = trim(adjustl(outdir))//"/eps.xg"
    call output_single_mass(eps,filename)

    filename = trim(adjustl(outdir))//"/mass.xg"
    call output_single_radius(mass,filename)

    filename = trim(adjustl(outdir))//"/temp.xg"
    call output_single_mass(temp,filename)

    filename = trim(adjustl(outdir))//"/lum.xg"
    call output_single_mass(lum,filename)

    filename = trim(adjustl(outdir))//"/tau.xg"
    call output_single_mass(tau,filename)

    filename = trim(adjustl(outdir))//"/delta_time.xg"
    call output_single_mass(delta_time,filename)

    filename = trim(adjustl(outdir))//"/radius.xg"
    call output_single_mass(r,filename)

    filename = trim(adjustl(outdir))//"/kappa.xg"
    call output_single_mass(kappa,filename)

    if (output_mode .eq. 0) then
       filename = trim(adjustl(outdir))//"/kappa_table.xg"
       call output_single_mass(kappa_table,filename)

       filename = trim(adjustl(outdir))//"/logR_op.xg"
       call output_single_mass(logR_op,filename)

       filename = trim(adjustl(outdir))//"/logT.xg"
       call output_single_mass(logT,filename)
    endif

    filename = trim(adjustl(outdir))//"/p_rad.xg"
    call output_single_mass(p_rad,filename)

    if (output_mode .eq. 0) then

    filename = trim(adjustl(outdir))//"/Ni_deposit_function.xg"
    call output_single_mass(Ni_deposit_function,filename)

    filename = trim(adjustl(outdir))//"/He_1.xg"
    call output_single_mass(ion_fractions(He_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/He_2.xg"
    call output_single_mass(ion_fractions(He_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/He_3.xg"
    call output_single_mass(ion_fractions(He_number,3,:),filename)

    filename = trim(adjustl(outdir))//"/H_1.xg"
    call output_single_mass(ion_fractions(H_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/H_2.xg"
    call output_single_mass(ion_fractions(H_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/free_electron_frac.xg"
    call output_single_mass(free_electron_frac,filename)

    filename = trim(adjustl(outdir))//"/E_shell.xg"
    call output_single_mass(E_shell,filename)

    filename = trim(adjustl(outdir))//"/time_diff.xg"
    call output_single_mass(time_diff,filename)

    filename = trim(adjustl(outdir))//"/time_exp.xg"
    call output_single_mass(time_exp,filename)

    filename = trim(adjustl(outdir))//"/photosphere_tracer.xg"
    call output_single_mass(photosphere_tracer,filename)

    endif
    endif

  else if(modeflag.eq.2) then

     filename = trim(adjustl(outdir))//"/lum_observed.dat"
     call output_scalar(lum_public,filename)

     if (.not.sedov .and. time.gt.0.0d0 .and. lum_public.gt.0.0d0) then
        open(666,file=trim(adjustl(outdir))//'/magnitudes.dat', &
             status='unknown',position='append')
        write(666,'(13E18.9)') time,T_color,magnitudes(1:11)
        close(666)
     endif

     call output_hydro_runtime_diagnostics()

     filename = trim(adjustl(outdir))//"/vel_photo.dat"
     call output_scalar(vel_photo,filename)

     if (output_mode .eq. 0 .or. output_mode .eq. 2) then

     filename = trim(adjustl(outdir))//"/T_eff.dat"
     call output_scalar(T_eff,filename)

     filename = trim(adjustl(outdir))//"/index_photo.dat"
     call output_integer(index_photo,filename)

     filename = trim(adjustl(outdir))//"/photometry_valid.dat"
     call output_integer(photometry_valid,filename)

     filename = trim(adjustl(outdir))//"/temperature_floor_fraction.dat"
     call output_scalar(temperature_floor_fraction,filename)

     filename = trim(adjustl(outdir))//"/mass_photo.dat"
     call output_scalar(mass_photo,filename)

     filename = trim(adjustl(outdir))//"/opacity_corrupted.dat"
     call output_integer(opacity_corrupted,filename)

     filename = trim(adjustl(outdir))//"/index_lumshell.dat"
     call output_integer(index_lumshell,filename)

     filename = trim(adjustl(outdir))//"/mass_lumshell.dat"
     call output_scalar(mass_lumshell,filename)

     endif

     filename = trim(adjustl(outdir))//"/rad_photo.dat"
     call output_scalar(rad_photo,filename)

  endif


end subroutine output_all


subroutine output_single_mass(var,filename)

  use blmod, only: mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(len=*) filename
  integer i


  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E29.20E3)") mass(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)


end subroutine output_single_mass


subroutine output_single_radius(var,filename)

  use blmod, only: r,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E19.10E3)") r(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)

end subroutine output_single_radius


subroutine output_scalar(var,filename)

  use blmod, only: time

  implicit none
  real*8 var
  character(len=*) filename

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var

  close(666)

end subroutine output_scalar

subroutine output_integer(var,filename)

  use blmod, only: time

  implicit none
  integer var
  character(len=*) filename

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(E19.10E3, I5.4)") time,var

  close(666)

end subroutine output_integer


subroutine output_screenshot(var,filename,imaximum)

  implicit none
  real*8 var(*)
  character(len=*) filename
  integer i
  integer imaximum

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
       form='formatted',position="append")

  do i=1, imaximum
      write(666,"(I5.4, E25.16E3)") i,var(i)
  enddo

  close(666)

end subroutine output_screenshot


subroutine output_hydro_runtime_diagnostics()

  use blmod, only: time, rho, r, cr, vel, cs2, shockpos, &
       shockpos_stop, radshock, cfl_limit_zone, cfl_raw_min_dt
  use parameters, only: imax, outdir
  use outinfomod, only: locate_live_csm_interaction
  implicit none

  integer :: u, ios, trace_idx, interaction_idx, maxrho_idx, cfl_idx
  real*8 :: qmax, qratio, trace_radius, interaction_radius, maxrho
  real*8 :: cfl_rho, cfl_dr, cfl_cs, cfl_compression
  logical :: exists
  character(len=1024) :: path

  call locate_live_csm_interaction(interaction_idx,qmax,qratio)

  trace_idx=max(1,min(imax,shockpos))
  interaction_idx=max(1,min(imax-1,interaction_idx))
  cfl_idx=max(1,min(imax-1,cfl_limit_zone))
  maxrho_idx=maxloc(rho(1:imax-1),dim=1)

  trace_radius=cr(trace_idx)
  interaction_radius=cr(interaction_idx)
  maxrho=rho(maxrho_idx)
  cfl_rho=rho(cfl_idx)
  cfl_dr=r(cfl_idx+1)-r(cfl_idx)
  cfl_cs=sqrt(max(cs2(cfl_idx),0.0d0))
  cfl_compression=max(0.0d0,vel(cfl_idx)-vel(cfl_idx+1))

  if(trace_radius.ne.trace_radius) trace_radius=0.0d0
  if(interaction_radius.ne.interaction_radius) interaction_radius=0.0d0
  if(maxrho.ne.maxrho) maxrho=0.0d0
  if(cfl_rho.ne.cfl_rho) cfl_rho=0.0d0
  if(cfl_dr.ne.cfl_dr) cfl_dr=0.0d0
  if(cfl_cs.ne.cfl_cs) cfl_cs=0.0d0
  if(cfl_compression.ne.cfl_compression) cfl_compression=0.0d0

  path=trim(adjustl(outdir))//'/hydro_runtime_diagnostics.dat'
  inquire(file=trim(path),exist=exists)
  open(newunit=u,file=trim(path),status='unknown',position='append',iostat=ios)
  if(ios.ne.0) return
  if(.not.exists) then
     write(u,'(A)') '# time_s shock_trace_zone shock_trace_stopped shock_trace_radius_cm ' // &
          'stored_radshock_cm interaction_zone interaction_radius_cm interaction_qmax ' // &
          'interaction_q_over_p maxrho_g_cm3 maxrho_zone cfl_zone cfl_rho_g_cm3 ' // &
          'cfl_dr_cm cfl_cs_cm_s cfl_compression_cm_s cfl_raw_dt_s'
  endif
  write(u,'(1P17E22.12)') time,dble(trace_idx),dble(shockpos_stop),trace_radius,radshock, &
       dble(interaction_idx),interaction_radius,qmax,qratio,maxrho,dble(maxrho_idx), &
       dble(cfl_idx),cfl_rho,cfl_dr,cfl_cs,cfl_compression,cfl_raw_min_dt
  close(u)

end subroutine output_hydro_runtime_diagnostics

subroutine write_run_summary(nsteps)

  use blmod, only: time, photosphere_fell_on_the_center, &
       diag_ni_eval_total, diag_ni_eval_remap
  use parameters
  use csm_switch, only: csm_generated_active, csm_base_profile_path, &
       csm_base_composition_path, csm_generated_profile_path, &
       csm_generated_composition_path, csm_mode_active, csm_preset_active, &
       csm_mix_summary, csm_generated_mass_msun
  implicit none

  integer(int64), intent(in) :: nsteps
  integer :: u
  character(len=1024) :: filepath

  u = 199
  filepath = trim(adjustl(outdir))//"/run_summary.dat"

  open(unit=u, file=trim(adjustl(filepath)), status='replace', form='formatted')

  write(u,'(A)') '# STRATA-NOVA grey-FLD run summary'
  write(u,'(A)') '# key = value'
  write(u,'(A)') '#'

  write(u,'(A)') '# --- Model identity ---'
  write(u,'(A,A)')     'profile_name = ', trim(adjustl(profile_name))
  if (csm_generated_active) then
     write(u,'(A,A)')  'base_profile_name = ', trim(adjustl(csm_base_profile_path))
     write(u,'(A,A)')  'base_comp_profile_name = ', trim(adjustl(csm_base_composition_path))
     write(u,'(A,A)')  'generated_profile_name = ', trim(adjustl(csm_generated_profile_path))
     write(u,'(A,A)')  'generated_comp_profile_name = ', trim(adjustl(csm_generated_composition_path))
     write(u,'(A,A)')  'csm = ', trim(adjustl(csm_mode_active))
     write(u,'(A,A)')  'csm_composition_preset = ', trim(adjustl(csm_preset_active))
     write(u,'(A,A)')  'csm_mix = ', trim(adjustl(csm_mix_summary))
     write(u,'(A,ES14.6)') 'csm_generated_mass_msun = ', csm_generated_mass_msun
     write(u,'(A,A)') 'mass_scheme = ', trim(adjustl(mass_scheme))
     write(u,'(A,ES14.6)') 'csm_outer_radius_cm = ', csm_outer_radius_cm
     if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss') then
        write(u,'(A,ES14.6)') 'csm_wind_velocity_kms = ', csm_wind_velocity_kms
     endif
  endif
  write(u,'(A,A)')     'initial_data = ', trim(adjustl(initial_data))


  write(u,'(A)') '# --- Physical parameters ---'
  write(u,'(A,ES14.6)') 'Ni_mass = ', Ni_mass
  write(u,'(A,F10.6)')  'Ni_mix_fraction = ', Ni_mix_fraction
  write(u,'(A,I4)')     'Ni_mix_kernel = ', Ni_mix_kernel
  write(u,'(A,I4)')     'boxcar_smooth_ni = ', boxcar_smooth_ni
  write(u,'(A,F10.6)')  'boxcar_nominal_mass_msun = ', boxcar_nominal_mass_msun
  write(u,'(A,I4)')     'boxcar_number_iterations = ', boxcar_number_iterations
  if (initial_data .eq. 'Thermal_Bomb') then
     write(u,'(A,ES14.6)') 'final_energy = ', final_energy
     write(u,'(A,ES14.6)') 'bomb_mass_spread = ', bomb_mass_spread
     write(u,'(A,I6)')     'bomb_start_point = ', bomb_start_point
  endif

  write(u,'(A)') '# --- Grid configuration ---'
  write(u,'(A,I6)')    'imax = ', imax
  write(u,'(A,A)')     'grid_mode = ', trim(adjustl(grid_mode))

  write(u,'(A)') '# --- Speed parameters ---'
  write(u,'(A,ES14.6)') 'Ni_period = ', Ni_period
  write(u,'(A,ES14.6)') 'Ni_period_max = ', Ni_period_max_target
  write(u,'(A,F10.6)')  'Ni_fractional_change = ', Ni_fractional_change_target
  write(u,'(A,I6)')     'Ni_quad_npoints = ', Ni_quad_npoints
  write(u,'(A,ES10.2)') 'epstol_hydro = ', epstol_hydro
  write(u,'(A,ES10.2)') 'epstol_rad = ', epstol_rad
  write(u,'(A,I6)')     'itmax_hydro = ', itmax_hydro
  write(u,'(A,I6)')     'itmax_rad = ', itmax_rad
  write(u,'(A,ES14.6)') 'dtmax = ', dtmax
  write(u,'(A,ES14.6)') 'dt_cfl_factor = ', dt_cfl_factor
  write(u,'(A,ES14.6)') 'dt_growth_factor = ', dt_growth_factor
  write(u,'(A,ES14.6)') 'evolution_temperature_floor_k = ', evolution_temperature_floor_k
  write(u,'(A,I6)')     'omp_num_threads = ', omp_num_threads

  write(u,'(A)') '# --- Run statistics ---'
  write(u,'(A,I0)')     'timesteps = ', nsteps
  write(u,'(A,I10)')    'Ni_evals = ', diag_ni_eval_total
  write(u,'(A,I10)')    'Ni_evals_remap = ', diag_ni_eval_remap
  write(u,'(A,ES14.6)') 'tend_actual = ', time
  write(u,'(A,I2)')     'photosphere_fell = ', photosphere_fell_on_the_center

  write(u,'(A)') '# --- Output configuration ---'
  write(u,'(A,I2)')  'output_mode = ', output_mode
  write(u,'(A)') 'luminosity_public_file = lum_observed.dat'
  write(u,'(A)') 'lum_observed_definition = all-frequency selected grey luminosity; ' // &
       'staged final applies the configured observer-time integration'
  ! keep the old summary keys for existing readers
  write(u,'(A,I2)') 'csm_luminosity_handover_done = ', 0
  write(u,'(A,ES14.6)') 'csm_luminosity_handover_time_s = ', -1.0d0
  write(u,'(A,I2)')  'eoskey = ', eoskey

  close(u)

end subroutine write_run_summary
