subroutine input_parser
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use blmod, only: gravity_switch, wipe_outdir
  use parameters
  use filesystem_utils, only: shell_quote, reset_output_directory
  use startup_log, only: print_startup_separator
  use omp_lib, only: omp_get_max_threads, omp_set_num_threads
  implicit none

  character(len=2048) :: cpstring
  logical :: opt
  logical :: outdirthere
  integer :: i_mode, c_mode
  integer :: legacy_boxcar_mode3_number_iterations
  real*8 :: legacy_boxcar_mode3_nominal_mass_msun
  character(len=32) :: csm_norm
  character(len=1024) :: protected_inputs(4)
  character(len=1024) :: parameter_file_path
  integer :: env_status, cmd_status


  opt = .false.

  parameter_file_path = 'parameters'
  call get_environment_variable('SNEC_PARAMETER_FILE', parameter_file_path, status=env_status)
  if (env_status.eq.-1) error stop 'SNEC_PARAMETER_FILE is too long'
  if (env_status.ne.0 .or. len_trim(parameter_file_path).eq.0) parameter_file_path='parameters'


  ! read named settings from the selected parameter file
  call get_string_parameter('outdir',outdir,opt)


  call get_string_parameter('profile_name',profile_name,opt)
  call get_string_parameter('comp_profile_name',composition_profile_name,opt)
  call get_string_parameter('restart_state_file',restart_state_file,.true.)
  if (trim(restart_state_file).eq.'NOTTHERE') restart_state_file=''


  call get_string_parameter('csm',csm,.true.)
  if (len_trim(csm).eq.0 .or. trim(adjustl(csm)).eq.'NOTTHERE') then
     csm = 'none'
  endif
  csm_norm = trim(adjustl(csm))
  do i_mode=1,len_trim(csm_norm)
     c_mode = iachar(csm_norm(i_mode:i_mode))
     if (c_mode.ge.iachar('A') .and. c_mode.le.iachar('Z')) then
        csm_norm(i_mode:i_mode) = achar(c_mode+32)
     endif
  enddo
  csm = trim(adjustl(csm_norm))
  if (trim(csm).ne.'none' .and. trim(csm).ne.'active') then
     write(*,*) '******* Unsupported csm: ',trim(csm)
     write(*,*) '******* Supported states: none, active'
     stop 1
  endif

  if (trim(csm).ne.'none') then
     call get_string_parameter('mass_scheme',mass_scheme,.true.)
     if (len_trim(mass_scheme).eq.0 .or. trim(adjustl(mass_scheme)).eq.'NOTTHERE') then
        mass_scheme = 'continuous_mass'
     endif
     do i_mode=1,len_trim(mass_scheme)
        c_mode = iachar(mass_scheme(i_mode:i_mode))
        if (c_mode.ge.iachar('A') .and. c_mode.le.iachar('Z')) then
           mass_scheme(i_mode:i_mode) = achar(c_mode+32)
        endif
     enddo
     if (trim(adjustl(mass_scheme)).eq.'continous_mass') mass_scheme='continuous_mass'
     mass_scheme=trim(adjustl(mass_scheme))
     if (trim(mass_scheme).ne.'continuous_mass' .and. &
          trim(mass_scheme).ne.'constant_mass_loss' .and. &
          trim(mass_scheme).ne.'non_constant_mass_loss' .and. &
          trim(mass_scheme).ne.'constant_density') then
        write(*,*) '******* Unsupported mass_scheme: ',trim(mass_scheme)
        write(*,*) '******* Supported schemes: continuous_mass, constant_mass_loss, non_constant_mass_loss, constant_density'
        stop 1
     endif

     if (trim(mass_scheme).eq.'continuous_mass') then
        call get_double_parameter('csm_mass_msun',csm_mass_msun,.true.)

     else if (trim(mass_scheme).eq.'constant_density') then
        call get_double_parameter('csm_mass_msun',csm_mass_msun,.true.)
        call get_double_parameter('csm_wind_velocity_kms',csm_wind_velocity_kms,.true.)

     else if (trim(mass_scheme).eq.'constant_mass_loss') then
        call get_double_parameter('csm_wind_velocity_kms',csm_wind_velocity_kms,.true.)
        call get_double_parameter('csm_mdot_msun_per_yr',csm_mdot_msun_per_yr,.true.)

     else if (trim(mass_scheme).eq.'non_constant_mass_loss') then
        call get_double_parameter('csm_wind_velocity_kms',csm_wind_velocity_kms,.true.)
     endif

     call get_double_parameter('csm_outer_radius_cm',csm_outer_radius_cm,.true.)
     call get_double_parameter('csm_temperature_k',csm_temperature_k,.true.)
     call get_integer_parameter('csm_extension_zones',csm_extension_zones,.true.)
     if (csm_extension_zones.eq.-666) csm_extension_zones=48
     call get_integer_parameter('csm_blend_zones',csm_blend_zones,.true.)
     if (csm_blend_zones.eq.-666) csm_blend_zones=8
     call get_string_parameter('csm_composition_preset',csm_composition_preset,.true.)
     if (len_trim(csm_composition_preset).eq.0 .or. &
          trim(adjustl(csm_composition_preset)).eq.'NOTTHERE') then
        csm_composition_preset='solar'
     endif
     call get_integer_parameter('csm_surface_mix_zones',csm_surface_mix_zones,.true.)
     if (csm_surface_mix_zones.eq.-666) csm_surface_mix_zones=10

     if (csm_outer_radius_cm.le.0.0d0) then
        write(*,*) '******* csm_outer_radius_cm must be > 0'
        stop 1
     endif
     if (csm_extension_zones.lt.2) then
        write(*,*) '******* csm_extension_zones must be >= 2'
        stop 1
     endif
     if (csm_blend_zones.lt.1 .or. csm_blend_zones.ge.csm_extension_zones) then
        write(*,*) '******* csm_blend_zones must lie in [1,csm_extension_zones-1]'
        stop 1
     endif
     if ((trim(mass_scheme).eq.'continuous_mass' .or. trim(mass_scheme).eq.'constant_density') .and. &
          csm_mass_msun.le.0.0d0) then
        write(*,*) '******* csm_mass_msun must be > 0 for direct-mass schemes'
        stop 1
     endif
     if (trim(mass_scheme).eq.'constant_density') then
        if (csm_wind_velocity_kms.lt.0.0d0 .or. csm_wind_velocity_kms.ne.csm_wind_velocity_kms) then
           write(*,*) '******* constant_density requires csm_wind_velocity_kms >= 0'
           stop 1
        endif
     endif
     if (trim(mass_scheme).eq.'constant_mass_loss') then
        if (csm_wind_velocity_kms.le.0.0d0 .or. &
             csm_wind_velocity_kms.ne.csm_wind_velocity_kms) then
           write(*,*) '******* constant_mass_loss requires csm_wind_velocity_kms > 0'
           stop 1
        endif
        if (csm_mdot_msun_per_yr.le.0.0d0 .or. &
             csm_mdot_msun_per_yr.ne.csm_mdot_msun_per_yr) then
           write(*,*) '******* constant_mass_loss requires csm_mdot_msun_per_yr > 0'
           stop 1
        endif
     endif
     if (trim(mass_scheme).eq.'non_constant_mass_loss') then
        if (csm_wind_velocity_kms.le.0.0d0 .or. &
             csm_wind_velocity_kms.ne.csm_wind_velocity_kms) then
           write(*,*) '******* non_constant_mass_loss requires csm_wind_velocity_kms > 0'
           stop 1
        endif
     endif
  endif


  call get_string_parameter('initial_data',initial_data,opt)

  if(initial_data.eq."Piston_Explosion") then
     call get_double_parameter('piston_vel',piston_vel,opt)
     call get_double_parameter('piston_tstart',piston_tstart,opt)
     call get_double_parameter('piston_tend',piston_tend,opt)
  endif

  if(initial_data.eq."Thermal_Bomb") then
     call get_double_parameter('final_energy',final_energy,opt)
     call get_double_parameter('bomb_tstart',bomb_tstart,opt)
     call get_double_parameter('bomb_tend',bomb_tend,opt)
     call get_double_parameter('bomb_mass_spread',bomb_mass_spread,opt)
     call get_integer_parameter('bomb_start_point',bomb_start_point,opt)
     call get_integer_parameter('bomb_mode',bomb_mode,.true.)
     if(bomb_mode.eq.-666) then
        bomb_mode = 1
     endif
  endif


   call get_integer_parameter('imax',imax,opt)
   call get_double_parameter('grid_update_interval_days',grid_update_interval_days,.true.)
   call get_string_parameter('grid_mode',grid_mode,.true.)
   do i_mode=1,len_trim(grid_mode)
      c_mode = iachar(grid_mode(i_mode:i_mode))
      if (c_mode.ge.iachar('A') .and. c_mode.le.iachar('Z')) then
         grid_mode(i_mode:i_mode) = achar(c_mode + 32)
      endif
   enddo
   if (len_trim(grid_mode).eq.0 .or. trim(adjustl(grid_mode)).eq.'notthere') grid_mode='legacy_native'
   if (trim(adjustl(grid_mode)).eq.'legacy_pattern') then
      write(*,*) ' legacy_pattern is deprecated; using Fortran-native legacy grid.'
      grid_mode='legacy_native'
   endif
   if (trim(adjustl(grid_mode)).ne.'legacy_native') then
      write(*,*) '******* Unsupported grid_mode: ',trim(grid_mode)
      write(*,*) '******* Final STRATA-NOVA production mode: legacy_native'
      stop 1
   endif
   grid_debug=0
   call get_logical_parameter('mass_excision',mass_excision,opt)
  if(mass_excision) then
     call get_double_parameter('mass_excised',mass_excised,opt)
  end if


  call get_logical_parameter('radiation',radiation,opt)
  call get_integer_parameter('eoskey',eoskey,opt)
  call get_integer_parameter('Ni_switch',Ni_switch,opt)
  call get_double_parameter('initial_model_age_s',initial_model_age_s,.true.)
  if (initial_model_age_s.eq.-666.0d0) initial_model_age_s=0.0d0
  if (initial_model_age_s.lt.0.0d0) then
     write(*,*) '******* initial_model_age_s must be >= 0'
     stop 1
  endif
  call get_double_parameter('Ni_mass',Ni_mass,opt)
  call get_double_parameter('Ni_mix_fraction',Ni_mix_fraction,opt)
  call get_integer_parameter('Ni_mix_kernel',Ni_mix_kernel,.true.)
  if (Ni_mix_kernel.eq.-666) Ni_mix_kernel = 2
  call get_double_parameter('Ni_mix_component2_fraction',Ni_mix_component2_fraction,.true.)
  call get_double_parameter('Ni_mix_component2_extent',Ni_mix_component2_extent,.true.)
  if (Ni_mix_component2_extent.lt.0.0d0) Ni_mix_component2_extent = Ni_mix_fraction
  if (Ni_mix_component2_fraction.lt.0.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_fraction < 0; clamping to 0.'
     Ni_mix_component2_fraction = 0.0d0
  else if (Ni_mix_component2_fraction.gt.1.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_fraction > 1; clamping to 1.'
     Ni_mix_component2_fraction = 1.0d0
  endif
  if (Ni_mix_component2_extent.lt.0.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_extent < 0; clamping to 0.'
     Ni_mix_component2_extent = 0.0d0
  else if (Ni_mix_component2_extent.gt.1.0d0) then
     write(*,*) 'WARNING: Ni_mix_component2_extent > 1; clamping to 1.'
     Ni_mix_component2_extent = 1.0d0
  endif
  call get_double_parameter('Ni_period',Ni_period,opt)
  Ni_period_min_target = Ni_period
  call get_double_parameter('Ni_period_max',Ni_period_max_target,opt)
  if (Ni_period_max_target.le.0.0d0) then
     Ni_period_max_target = 10.0d0 * Ni_period_min_target
  endif
  call get_double_parameter('Ni_fractional_change',Ni_fractional_change_target,opt)
  if (Ni_fractional_change_target.le.0.0d0) then
     Ni_fractional_change_target = 0.25d0
  endif
  call get_integer_parameter('Ni_quad_npoints',Ni_quad_npoints,opt)
  if (Ni_quad_npoints.le.0) Ni_quad_npoints = 150
  call get_integer_parameter('Ni_by_hand',Ni_by_hand,.true.)
  if (Ni_by_hand.eq.-666) Ni_by_hand = 1

  call get_integer_parameter('saha_ncomps',saha_ncomps,opt)

  call get_logical_parameter('boxcar_smoothing',boxcar_smoothing,opt)
  call get_integer_parameter('boxcar_smooth_ni',boxcar_smooth_ni,opt)
  call get_double_parameter('boxcar_nominal_mass_msun',boxcar_nominal_mass_msun,.true.)
  if (boxcar_nominal_mass_msun.le.0.0d0) then
     write(*,*) 'WARNING: boxcar_nominal_mass_msun <= 0; using default 0.4'
     boxcar_nominal_mass_msun = 0.4d0
  endif
  call get_integer_parameter('boxcar_number_iterations',boxcar_number_iterations,.true.)
  if (boxcar_number_iterations.eq.-666) then
     boxcar_number_iterations = 4
  else if (boxcar_number_iterations.le.0) then
     write(*,*) 'WARNING: boxcar_number_iterations <= 0; using default 4'
     boxcar_number_iterations = 4
  endif
  call get_double_parameter('boxcar_ni_nominal_mass_msun',boxcar_ni_nominal_mass_msun,.true.)
  legacy_boxcar_mode3_nominal_mass_msun = -1.23456789d307
  call get_double_parameter('boxcar_mode3_nominal_mass_msun',legacy_boxcar_mode3_nominal_mass_msun,.true.)
  if (boxcar_ni_nominal_mass_msun.eq.-1.0d0 .and. &
       legacy_boxcar_mode3_nominal_mass_msun.ne.-1.23456789d307) then
     boxcar_ni_nominal_mass_msun = legacy_boxcar_mode3_nominal_mass_msun
     write(*,*) 'WARNING: boxcar_mode3_nominal_mass_msun is deprecated; use boxcar_ni_nominal_mass_msun.'
  endif
  if (boxcar_ni_nominal_mass_msun.eq.0.0d0) then
     write(*,*) 'WARNING: boxcar_ni_nominal_mass_msun = 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_nominal_mass_msun = -1.0d0
  else if (boxcar_ni_nominal_mass_msun.lt.0.0d0 .and. &
       boxcar_ni_nominal_mass_msun.ne.-1.0d0) then
     write(*,*) 'WARNING: boxcar_ni_nominal_mass_msun < 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_nominal_mass_msun = -1.0d0
  endif

  call get_integer_parameter('boxcar_ni_number_iterations',boxcar_ni_number_iterations,.true.)
  if (boxcar_ni_number_iterations.eq.-666) then
     boxcar_ni_number_iterations = -1
  endif
  legacy_boxcar_mode3_number_iterations = -666
  call get_integer_parameter('boxcar_mode3_number_iterations',legacy_boxcar_mode3_number_iterations,.true.)
  if (boxcar_ni_number_iterations.eq.-1 .and. legacy_boxcar_mode3_number_iterations.ne.-666) then
     boxcar_ni_number_iterations = legacy_boxcar_mode3_number_iterations
     write(*,*) 'WARNING: boxcar_mode3_number_iterations is deprecated; use boxcar_ni_number_iterations.'
  endif
  if (boxcar_ni_number_iterations.eq.0) then
     write(*,*) 'WARNING: boxcar_ni_number_iterations = 0; disabling secondary Ni-smoothing override.'
     boxcar_ni_number_iterations = -1
  else if (boxcar_ni_number_iterations.lt.-1) then
     write(*,*) 'WARNING: boxcar_ni_number_iterations < -1; disabling secondary Ni-smoothing override.'
     boxcar_ni_number_iterations = -1
  endif

  call get_double_parameter('opacity_floor_envelope',of_env,opt)
  call get_double_parameter('opacity_floor_core',of_core,opt)


  call get_int64_parameter('ntmax',ntmax,opt)
  if (ntmax.le.0) error stop 'ntmax must be positive'
  call get_double_parameter('tend',tend,opt)
  if (tend.le.0.0d0) error stop 'tend must be positive'
  call get_double_parameter('dtout',dtout,opt)
  dtout_late_target = dtout
  call get_double_parameter('dtout_mid',dtout_mid_target,opt)
  call get_double_parameter('dtout_fast',dtout_fast_target,opt)
  call get_double_parameter('dtout_scalar',dtout_scalar,opt)
  dtout_scalar_late_target = dtout_scalar
  call get_double_parameter('dtout_scalar_mid',dtout_scalar_mid_target,opt)
  call get_double_parameter('dtout_scalar_fast',dtout_scalar_fast_target,opt)
  call get_double_parameter('output_mid_transition_days',output_mid_transition_days,opt)
  call get_double_parameter('output_late_transition_days',output_late_transition_days,opt)
  call get_double_parameter('dtout_check',dtout_check,opt)
  call get_integer_parameter('ntout',ntout,opt)
  call get_integer_parameter('ntout_scalar',ntout_scalar,opt)
  call get_integer_parameter('ntout_check',ntout_check,opt)
  call get_integer_parameter('ntinfo',ntinfo,opt)
  call get_double_parameter('dtmin',dtmin,opt)
  call get_double_parameter('dtmax',dtmax,opt)
  call get_double_parameter('dt_cfl_factor',dt_cfl_factor,.true.)
  if (dt_cfl_factor.eq.-666.0d0) dt_cfl_factor = 0.95d0
  dt_cfl_factor = max(1.0d-3, min(1.0d0, dt_cfl_factor))
  call get_double_parameter('dt_growth_factor',dt_growth_factor,.true.)
  if (dt_growth_factor.eq.-666.0d0) dt_growth_factor = 1.025d0
  dt_growth_factor = max(1.0d0, dt_growth_factor)
  call get_double_parameter('evolution_temperature_floor_k', &
       evolution_temperature_floor_k,.true.)
  if (evolution_temperature_floor_k.eq.-666.0d0) evolution_temperature_floor_k = 1.00d2
  if (evolution_temperature_floor_k.le.0.0d0) then
     write(*,*) 'WARNING: evolution_temperature_floor_k must be positive; using 100 K'
     evolution_temperature_floor_k = 1.00d2
  endif
  call get_double_parameter('hydro_accel_cfl',hydro_accel_cfl,.true.)
  if (hydro_accel_cfl.eq.-666.0d0) hydro_accel_cfl=0.25d0
  hydro_accel_cfl=max(0.02d0,min(1.0d0,hydro_accel_cfl))
  call get_integer_parameter('omp_num_threads', omp_num_threads, .true.)
  if (omp_num_threads.eq.-666) omp_num_threads = 0
  if (omp_num_threads.gt.0) then
     call omp_set_num_threads(omp_num_threads)
     write(*,'(A,I4,A,I4)') ' OpenMP threads requested: ', &
          omp_num_threads, '  runtime max now: ', omp_get_max_threads()
  else
     omp_num_threads = 0
  endif


  call get_integer_parameter('smooth_ni_luminosity',smooth_ni_luminosity,.true.)
  if (smooth_ni_luminosity.eq.-666) smooth_ni_luminosity = 1


  call get_double_parameter('epstol_hydro',epstol_hydro,.true.)
  if (epstol_hydro .le. 0.0d0) epstol_hydro = 1.0d-4
  call get_double_parameter('epstol_rad',epstol_rad,.true.)
  if (epstol_rad .le. 0.0d0) epstol_rad = 1.0d-4
  call get_integer_parameter('itmax_hydro',itmax_hydro,.true.)
  if (itmax_hydro .le. 0) itmax_hydro = 100
  call get_integer_parameter('itmax_rad',itmax_rad,.true.)
  if (itmax_rad .le. 0) itmax_rad = 300
  write(*,'(A,ES9.2,A,ES9.2,A,I4,A,I4)') &
       ' Solver: hydro=', epstol_hydro, '  rad=', epstol_rad, &
       '  iterations=', itmax_hydro, '/', itmax_rad


  call get_integer_parameter('ni_raytrace_opt',ni_raytrace_opt,.true.)
  if (ni_raytrace_opt.eq.-666) ni_raytrace_opt = 1

  call get_integer_parameter('ni_ray_interp',ni_ray_interp,.true.)
  if (ni_ray_interp.eq.-666) ni_ray_interp = 1


  call get_integer_parameter('output_mode',output_mode,.true.)
  if (output_mode.eq.-666) output_mode = 0


  call get_logical_parameter('sedov',sedov,opt)
  if(sedov) then
     gravity_switch = 0
  else
     gravity_switch = 1
  end if


#if __INTEL_COMPILER
  inquire(directory=trim(adjustl(outdir)),exist=outdirthere)
#else
  inquire(file=trim(adjustl(outdir)),exist=outdirthere)
#endif
  if(.not.outdirthere) then
     write(6,*) "*** Output directory does not exist."
     write(6,*) "Please create the output directory: ", trim(adjustl(outdir))
     stop 1
  endif

  call print_startup_separator
  if(wipe_outdir) then
     protected_inputs(1)=parameter_file_path
     protected_inputs(2)=profile_name
     protected_inputs(3)=composition_profile_name
     protected_inputs(4)=restart_state_file
     call reset_output_directory(trim(outdir),protected_inputs)
     write(*,'(A,A,A)') ' Output directory: ',trim(outdir),'  [reset]'
  else
     write(*,'(A,A,A)') ' Output directory: ',trim(outdir),'  [preserved]'
  endif
  call print_startup_separator

  cpstring='cp -- '//trim(shell_quote(trim(adjustl(parameter_file_path))))//' '// &
       trim(shell_quote(trim(adjustl(outdir))//'/parameters'))
  call execute_command_line(trim(cpstring),wait=.true.,exitstat=cmd_status)
  if (cmd_status.ne.0) then
     write(*,*) 'ERROR: could not copy parameters into output directory.'
     stop 1
  endif


contains


    subroutine get_string_parameter(parname,par,opt)
        character(len=*), intent(in) :: parname
        character(len=*), intent(out) :: par
        logical, intent(in) :: opt
        character(len=4096) :: line, value, key
        character :: quote
        integer :: u, ios, i, eq, n

        open(newunit=u,file=trim(parameter_file_path),status='old',action='read',iostat=ios)
        if (ios.ne.0) error stop 'Cannot open parameter file'
        do
           read(u,'(A)',iostat=ios) line
           if (ios.lt.0) exit
           if (ios.ne.0) error stop 'Cannot read parameter file'
           if (len_trim(line).eq.len(line)) error stop 'Parameter line is too long'
           line=adjustl(line)
           if (line(1:1).eq.'#' .or. line(1:1).eq.'!') cycle
           eq=index(line,'=')
           if (eq.le.1) cycle
           key=line(1:eq-1)
           do i=1,len_trim(key)
              if (iachar(key(i:i)).eq.9) key(i:i)=' '
           enddo
           if (trim(key).ne.parname) cycle
           value=adjustl(line(eq+1:))
           quote=' '
           do i=1,len_trim(value)
              if (quote.ne.' ') then
                 if (value(i:i).eq.quote) quote=' '
              else if (value(i:i).eq.'"' .or. value(i:i).eq."'") then
                 quote=value(i:i)
              else if (value(i:i).eq.'#' .or. value(i:i).eq.'!') then
                 value(i:)=' '
                 exit
              endif
           enddo
           if (quote.ne.' ') error stop 'Unclosed quote in parameter file'
           n=len_trim(value)
           if (n.gt.0) then
              quote=value(1:1)
              if (quote.eq.'"' .or. quote.eq."'") then
                 if (n.lt.2) error stop 'Invalid quoted parameter'
                 if (value(n:n).ne.quote) error stop 'Text after quoted parameter'
                 line=' '
                 if (n.gt.2) line(1:n-2)=value(2:n-1)
                 value=line
              endif
           endif
           if (len_trim(value).gt.len(par)) then
              write(*,*) 'ERROR: parameter value is too long: ',parname
              stop 1
           endif
           par=trim(value)
           close(u)
           return
        enddo
        close(u)
        if (.not.opt) then
           write(*,*) 'ERROR: missing parameter: ',parname
           stop 1
        endif
        par='NOTTHERE'
    end subroutine get_string_parameter

  ! reject malformed values before converting the numeric token
    subroutine check_numeric_token(parname,value,alphabet)
        character(len=*), intent(in) :: parname,value,alphabet
        if (len_trim(value).eq.0 .or. verify(trim(value),alphabet).ne.0) then
           write(*,*) 'ERROR: invalid numeric parameter: ',parname,' = ',trim(value)
           stop 1
        endif
    end subroutine check_numeric_token

  ! read one finite real value using the full parameter key
    subroutine get_double_parameter(parname,par,opt)
        character(len=*), intent(in) :: parname
        real*8, intent(inout) :: par
        logical, intent(in) :: opt
        character(len=256) :: value
        integer :: ios
        call get_string_parameter(parname,value,opt)
        if (opt .and. value.eq.'NOTTHERE') return
        call check_numeric_token(parname,value,'0123456789.+-dDeE')
        read(value,*,iostat=ios) par
        if (ios.ne.0) then
           write(*,*) 'ERROR: invalid real parameter: ',parname
           stop 1
        endif
        if (.not.ieee_is_finite(par)) error stop 'Non-finite real parameter'
    end subroutine get_double_parameter

  ! check the integer range before accepting a timestep or count limit
    subroutine get_int64_parameter(parname,par,opt)
        character(len=*), intent(in) :: parname
        integer(int64), intent(out) :: par
        logical, intent(in) :: opt
        character(len=256) :: value
        integer :: ios
        call get_string_parameter(parname,value,opt)
        if (opt .and. value.eq.'NOTTHERE') then
           par=-666
           return
        endif
        call check_numeric_token(parname,value,'0123456789+-')
        read(value,*,iostat=ios) par
        if (ios.ne.0) then
           write(*,*) 'ERROR: invalid or overflowing integer parameter: ',parname
           stop 1
        endif
    end subroutine get_int64_parameter

    subroutine get_integer_parameter(parname,par,opt)
        character(len=*), intent(in) :: parname
        integer, intent(out) :: par
        logical, intent(in) :: opt
        integer(int64) :: value
        call get_int64_parameter(parname,value,opt)
        if (value.gt.int(huge(par),int64) .or. value.lt.-int(huge(par),int64)-1) &
             error stop 'Integer parameter is out of range'
        par=int(value)
    end subroutine get_integer_parameter

    subroutine get_logical_parameter(parname,par,opt)
        character(len=*), intent(in) :: parname
        logical, intent(inout) :: par
        logical, intent(in) :: opt
        integer :: value
        character(len=256) :: text
        call get_string_parameter(parname,text,opt)
        if (opt .and. text.eq.'NOTTHERE') return
        call get_integer_parameter(parname,value,.false.)
        par=value.ne.0
    end subroutine get_logical_parameter
end subroutine input_parser
