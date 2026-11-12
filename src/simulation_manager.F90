module simulation_manager

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use staged_state, only: write_staged_state
  use observer_projection, only: equal_arrival_time_curve, linear_interp
  use filesystem_utils, only: shell_quote, check_output_directory, reset_output_directory
  use startup_log, only: print_startup_separator, print_stage_banner, is_stage_child
  implicit none

  private
  public :: simulation_manager_dispatch
  public :: simulation_manager_write_child_handoff

  real*8, parameter :: day_seconds = 86400.0d0
  real*8, parameter :: sigma_sb = 5.670374419d-5
  real*8, parameter :: c_light = 2.99792458d10
  real*8, parameter :: pi_mgr = 3.1415926535897932384626433832795d0
  integer, parameter :: path_len = 1024
  integer, parameter :: line_len = 4096
  integer, parameter :: max_overrides = 96
  integer, parameter :: staged_continuous_shell_evolution_fixed = 1

  type :: stage_state
     real*8 :: elapsed_time_s = -1.0d0
     character(len=32) :: termination_reason = ''
     real*8 :: mass_g = 0.0d0
     real*8 :: radius_cm = 0.0d0
     real*8 :: rho_outer = 0.0d0
     real*8 :: temp_outer = 0.0d0
     real*8 :: lum_grey_selected = 0.0d0
     real*8 :: max_abs_velocity_cms = 0.0d0
  end type stage_state

contains

  ! route staged runs through fresh solver processes
  subroutine simulation_manager_dispatch(handled)
    logical, intent(out) :: handled

    character(len=64) :: simulation_mode, csm_setting
    character(len=path_len) :: manager_parameter_path
    integer :: stat

    handled = .false.
    if (is_stage_child()) return

    manager_parameter_path='parameters'
    call get_environment_variable('SNEC_PARAMETER_FILE',manager_parameter_path,status=stat)
    if (stat.eq.-1) error stop 'SNEC_PARAMETER_FILE is too long'
    if (stat.ne.0 .or. len_trim(manager_parameter_path).eq.0) manager_parameter_path='parameters'
    call manager_get_string(trim(manager_parameter_path), 'simulation_mode', simulation_mode, 'single')
    call lowercase_inplace(simulation_mode)
    call manager_get_string(trim(manager_parameter_path), 'csm', csm_setting, 'none')
    call lowercase_inplace(csm_setting)
    if (trim(adjustl(csm_setting)).ne.'none' .and. trim(adjustl(csm_setting)).ne.'active') then
       write(*,*) '******* Unsupported csm setting: ',trim(csm_setting)
       write(*,*) '******* Supported states: none, active'
       stop 1
    endif

    select case (trim(adjustl(simulation_mode)))
    case ('single')
       call validate_mode_parameter_shapes(trim(manager_parameter_path),'single')
       return
    case ('staged_pipeline')
       if (trim(adjustl(csm_setting)).ne.'active') then
          write(*,*) '******* staged_pipeline requires csm = "active".'
          stop 1
       endif
       call validate_mode_parameter_shapes(trim(manager_parameter_path),'staged_pipeline')
       handled = .true.
       call run_staged_pipeline(trim(manager_parameter_path))
    case default
       write(*,*) '******* Unsupported simulation_mode: ', trim(simulation_mode)
       write(*,*) '******* Supported modes: single, staged_pipeline'
       stop 1
    end select
  end subroutine simulation_manager_dispatch


  ! export the final live state for the next encounter
  subroutine simulation_manager_write_child_handoff()
    use blmod, only: mass, r, temp, rho, vel, comp, comp_details, ncomps, ye, &
         lum_public, time
    use parameters, only: imax

    character(len=path_len) :: profile_path, comp_path, meta_path, handoff_dir
    integer :: stat, i, j, up, uc, um
    real*8 :: sumx, xnorm, ye_local

    profile_path = ''
    comp_path = ''
    meta_path = ''
    call get_environment_variable('STRATA_HANDOFF_PROFILE', profile_path, status=stat)
    if (stat.ne.0 .or. len_trim(profile_path).eq.0) return
    call get_environment_variable('STRATA_HANDOFF_COMP', comp_path, status=stat)
    if (stat.ne.0 .or. len_trim(comp_path).eq.0) return
    call get_environment_variable('STRATA_HANDOFF_META', meta_path, status=stat)
    if (stat.ne.0 .or. len_trim(meta_path).eq.0) return

    call parent_directory(profile_path,handoff_dir)
    call mkdir_p(trim(handoff_dir))

    open(newunit=up, file=trim(profile_path), status='replace', action='write')
    write(up,'(I10)') imax
    do i=1,imax
       ye_local = 0.0d0
       sumx = sum(max(comp(i,1:ncomps),0.0d0))
       if (sumx.le.0.0d0) then
          ye_local = ye(i)
       else
          do j=1,ncomps
             xnorm = max(comp(i,j),0.0d0)/sumx
             ye_local = ye_local + xnorm*comp_details(j,2)/max(comp_details(j,1),1.0d0)
          enddo
       endif
       if (i.eq.imax) then
          write(up,'(I8,1X,7ES26.16E3)') i, mass(i), r(i), temp(imax-1), rho(imax-1), &
               vel(i), ye_local, 0.0d0
       else
          write(up,'(I8,1X,7ES26.16E3)') i, mass(i), r(i), temp(i), rho(i), &
               vel(i), ye_local, 0.0d0
       endif
    enddo
    close(up)

    open(newunit=uc, file=trim(comp_path), status='replace', action='write')
    write(uc,'(I10,1X,I10)') imax, ncomps
    write(uc,'(*(ES18.10E3,1X))') comp_details(1:ncomps,1)
    write(uc,'(*(ES18.10E3,1X))') comp_details(1:ncomps,2)
    do i=1,imax
       sumx = sum(max(comp(i,1:ncomps),0.0d0))
       write(uc,'(2ES26.16E3,1X)',advance='no') mass(i), r(i)
       do j=1,ncomps
          if (sumx.gt.0.0d0) then
             xnorm = max(comp(i,j),0.0d0)/sumx
          else
             xnorm = comp(i,j)
          endif
          write(uc,'(ES26.16E3,1X)',advance='no') xnorm
       enddo
       write(uc,*)
    enddo
    close(uc)

    open(newunit=um, file=trim(meta_path), status='replace', action='write')
    write(um,'(A,ES26.16E3)') 'elapsed_time_s = ',time
    write(um,'(A)') 'termination_reason = tend'
    write(um,'(A,ES26.16E3)') 'mass_g = ', mass(imax)
    write(um,'(A,ES26.16E3)') 'radius_cm = ', r(imax)
    write(um,'(A,ES26.16E3)') 'rho_outer = ', rho(max(1,imax-1))
    write(um,'(A,ES26.16E3)') 'temp_outer = ', temp(max(1,imax-1))
    write(um,'(A,ES26.16E3)') 'lum_grey_selected = ', max(lum_public,0.0d0)
    write(um,'(A,ES26.16E3)') 'max_abs_velocity_cms = ', maxval(abs(vel(1:imax)))
    close(um)
    call write_staged_state(trim(profile_path)//'.state')
  end subroutine simulation_manager_write_child_handoff


  ! evolve the ejecta through the configured physical shell radii
  subroutine run_staged_pipeline(parameter_path)
    character(len=*), intent(in) :: parameter_path

    integer :: n_shells, i, nstages, maxstages, exitstat
    integer :: inherited_imax, inherited_ni_quad, inherited_output_mode
    integer :: inherited_extension_zones, inherited_blend_zones
    integer :: staged_equal_arrival, staged_observer_nmu
    real*8 :: staged_final_day, inherited_dtmax, inherited_csm_mdot
    real*8 :: inherited_csm_velocity, inherited_csm_temperature
    real*8 :: initial_age_s, final_s, current_abs_s, shell_start_s
    real*8 :: shell_end_s, radius_attach, shell_temp, initial_radius, initial_velocity
    real*8, allocatable :: encounter_days(:), shell_inner_radii(:), shell_mass(:), shell_width(:)
    real*8, allocatable :: shell_velocity(:)
    real*8, allocatable :: shell_temperature(:)
    integer, allocatable :: shell_extension_zones(:), shell_blend_zones(:)
    real*8, allocatable :: attachment_radius(:), used_shell_temperature(:)
    real*8, allocatable :: attached_shell_mass(:)
    character(len=path_len) :: workdir, base_profile, base_comp
    character(len=path_len) :: protected_inputs(3)
    character(len=32) :: staged_mass_scheme
    character(len=path_len) :: exe_path, scratch_dir, final_dir, archived_log
    character(len=path_len) :: current_profile, current_comp
    character(len=path_len) :: stage_dir, stage_param, stage_log
    character(len=path_len) :: handoff_profile, handoff_comp, handoff_meta
    character(len=128) :: stage_label
    character(len=path_len), allocatable :: stage_dirs(:)
    real*8, allocatable :: stage_offsets(:)
    type(stage_state) :: state
    logical :: final_stage, direct_mass_scheme, wind_mass_scheme

    call manager_get_integer(parameter_path,'staged_n_shells',n_shells,-1)
    if (n_shells.lt.1) then
       write(*,*) '******* staged_pipeline requires staged_n_shells >= 1'
       stop 1
    endif

    allocate(encounter_days(n_shells), shell_inner_radii(n_shells), shell_mass(n_shells), shell_width(n_shells))
    encounter_days(:)=0.0d0
    allocate(shell_velocity(n_shells))
    allocate(shell_temperature(n_shells), shell_extension_zones(n_shells))
    allocate(shell_blend_zones(n_shells), attachment_radius(n_shells))
    allocate(used_shell_temperature(n_shells), attached_shell_mass(n_shells))
    attached_shell_mass(:)=0.0d0

    call manager_get_string(parameter_path,'mass_scheme',staged_mass_scheme,'constant_density')
    inherited_csm_mdot=-1.0d0
    direct_mass_scheme=is_direct_mass_scheme(staged_mass_scheme)
    wind_mass_scheme=is_wind_mass_scheme(staged_mass_scheme)
    if (.not.direct_mass_scheme .and. .not.wind_mass_scheme) then
       write(*,*) '******* Unsupported staged mass_scheme: ',trim(staged_mass_scheme)
       stop 1
    endif
    if (trim(adjustl(staged_mass_scheme)).eq.'constant_mass_loss') then
       call manager_get_double(parameter_path,'csm_mdot_msun_per_yr',inherited_csm_mdot,-1.0d0)
       if (inherited_csm_mdot.le.0.0d0) then
          write(*,*) '******* constant_mass_loss requires csm_mdot_msun_per_yr > 0'
          stop 1
       endif
    endif

    call manager_get_real_array(parameter_path,'staged_shell_inner_radii_cm',n_shells,shell_inner_radii)
    if (direct_mass_scheme) then
       call manager_get_real_array(parameter_path,'staged_shell_masses_msun',n_shells,shell_mass)
    else
       shell_mass(:)=0.0d0
    endif
    call manager_get_real_array(parameter_path,'staged_shell_widths_cm',n_shells,shell_width)

    call manager_get_double(parameter_path,'csm_wind_velocity_kms',inherited_csm_velocity,100.0d0)
    call manager_get_double(parameter_path,'csm_temperature_k',inherited_csm_temperature,-1.0d0)
    shell_velocity(:)=inherited_csm_velocity
    shell_temperature(:)=inherited_csm_temperature
    call manager_get_string(parameter_path,'staged_workdir',workdir,'staged_out')
    call manager_get_string(parameter_path,'profile_name',base_profile,'NOTTHERE')
    call manager_get_string(parameter_path,'comp_profile_name',base_comp,'NOTTHERE')
    if (trim(base_profile).eq.'NOTTHERE' .or. trim(base_comp).eq.'NOTTHERE') then
       write(*,*) '******* staged_pipeline uses profile_name and comp_profile_name as the base model.'
       stop 1
    endif

    call manager_get_double(parameter_path,'staged_final_day',staged_final_day,-1.0d0)
    if (staged_final_day.le.0.0d0) then
       write(*,*) '******* staged_final_day must be > 0'
       stop 1
    endif
    call manager_get_integer(parameter_path,'imax',inherited_imax,100)
    call manager_get_integer(parameter_path,'Ni_quad_npoints',inherited_ni_quad,150)
    call manager_get_integer(parameter_path,'output_mode',inherited_output_mode,0)
    call manager_get_double(parameter_path,'dtmax',inherited_dtmax,300.0d0)
    call manager_get_integer(parameter_path,'csm_extension_zones',inherited_extension_zones,48)
    call manager_get_integer(parameter_path,'csm_blend_zones',inherited_blend_zones,8)
    call manager_get_integer(parameter_path,'staged_equal_arrival_time',staged_equal_arrival,1)
    call manager_get_integer(parameter_path,'staged_observer_nmu',staged_observer_nmu,32)

    shell_extension_zones(:)=inherited_extension_zones
    shell_blend_zones(:)=inherited_blend_zones

    call validate_staged_inputs(n_shells,staged_mass_scheme,shell_inner_radii,shell_mass,shell_width, &
         shell_velocity,shell_extension_zones,shell_blend_zones, &
         inherited_imax,inherited_ni_quad,staged_observer_nmu, &
         inherited_output_mode)

    call infer_homologous_age(trim(base_profile), initial_age_s)
    final_s = staged_final_day*day_seconds
    call read_profile_surface(trim(base_profile),initial_radius,initial_velocity)
    if (shell_inner_radii(1).le.initial_radius) then
       write(*,*) '******* First shell radius must exceed the starting ejecta radius'
       stop 1
    endif

    call get_command_argument(0,exe_path)
    if (len_trim(exe_path).eq.0) exe_path='./snec'

    protected_inputs(1)=parameter_path
    protected_inputs(2)=base_profile
    protected_inputs(3)=base_comp
    call prepare_workdir(trim(workdir),protected_inputs)
    final_dir = trim(workdir)//'/final'
    scratch_dir = trim(workdir)//'/.manager_tmp'
    call mkdir_p(trim(scratch_dir))

    maxstages = n_shells + 1
    allocate(stage_dirs(maxstages),stage_offsets(maxstages))
    stage_dirs(:)=''
    stage_offsets(:)=0.0d0
    nstages=0

    call print_startup_separator
    write(*,'(A)') ' Native staged simulation manager'
    write(*,'(A,I0,A,F9.4,A,F9.4,A)') ' Shells: ',n_shells, &
         '   Base age: ',initial_age_s/day_seconds,' d   Final age: ',staged_final_day,' d'
    write(*,'(A,A)') ' Staged CSM mass scheme: ',trim(adjustl(staged_mass_scheme))
    if (wind_mass_scheme) then
       if (trim(adjustl(staged_mass_scheme)).eq.'constant_mass_loss') then
          write(*,'(A,ES12.4,A)') ' CSM Mdot: ',inherited_csm_mdot,' Msun/yr'
       else
          write(*,'(A)') ' Badenes CSM mass follows the pre-explosion wind history'
       endif
       write(*,'(A)') ' staged_shell_masses_msun is not read for this mass scheme.'
    endif
    write(*,'(A)') ' Shell encounters occur at their configured physical radii'
    write(*,'(A)') ' Parent csm_outer_radius_cm is ignored by staged_pipeline.'
    call print_startup_separator

    stage_dir = trim(workdir)//'/stage000_base'
    stage_param = trim(scratch_dir)//'/stage000_base.parameters'
    stage_log = trim(scratch_dir)//'/stage000_base.log'
    call stage_handoff_paths(stage_dir,handoff_profile,handoff_comp,handoff_meta)
    call run_stage_to_radius(parameter_path,exe_path,stage_param,stage_log,stage_dir, &
         base_profile,base_comp,handoff_profile,handoff_comp,handoff_meta, &
         initial_age_s,final_s,shell_inner_radii(1),initial_radius,initial_velocity, &
         inherited_imax,inherited_ni_quad,inherited_output_mode,inherited_dtmax, &
         .false.,.false.,direct_mass_scheme,0.0d0,0.0d0,0.0d0,1000.0d0,48,3, &
         current_abs_s,state)
    call archive_stage_control_files(stage_param,stage_log,stage_dir,archived_log)
    write(*,'(A)') ' Stage 0 reached the first shell radius'
    call append_stage(stage_dirs,stage_offsets,nstages,stage_dir,initial_age_s)
    call require_stage_end(state,current_abs_s-initial_age_s)
    call validate_nonrelativistic_stage(state,'stage 0')
    current_profile=handoff_profile
    current_comp=handoff_comp

    do i=1,n_shells
       shell_start_s = current_abs_s
       encounter_days(i) = shell_start_s/day_seconds

       radius_attach=state%radius_cm
       attachment_radius(i)=radius_attach
       if (shell_temperature(i).gt.0.0d0) then
          shell_temp=shell_temperature(i)
       else
          shell_temp=equilibrium_shell_temperature(state%lum_grey_selected,radius_attach)
       endif
       used_shell_temperature(i)=shell_temp

       if (i.lt.n_shells) then
          if (shell_inner_radii(i+1).le.radius_attach+shell_width(i)) then
             write(*,*) '******* Shell radii overlap after the preceding shell is attached'
             stop 1
          endif
       endif
       final_stage = i.eq.n_shells
       if (final_stage) then
          stage_dir = trim(final_dir)
       else
          write(stage_dir,'(A,"/stage",I3.3,"_shell")') trim(workdir),i
       endif
       write(stage_param,'(A,"/stage",I3.3,"_shell.parameters")') trim(scratch_dir),i
       write(stage_log,'(A,"/stage",I3.3,"_shell.log")') trim(scratch_dir),i
       call stage_handoff_paths(stage_dir,handoff_profile,handoff_comp,handoff_meta)

       write(stage_label,'(A,I0,A)') 'SHELL ',i,' INTERACTION'
       write(*,'(A,A)') ' Starting physical interaction: ',trim(stage_label)
       if (direct_mass_scheme) then
          write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)') ' Shell parameters: M=',shell_mass(i), &
               ' Msun   R=',radius_attach,' cm   width=',shell_width(i)
       else
          write(*,'(A,A,A,ES11.3,A,ES11.3)') ' Shell parameters: M=derived (', &
               trim(adjustl(staged_mass_scheme)),')   R=',radius_attach,' cm   width=',shell_width(i)
       endif

       if (final_stage) then
          shell_end_s=final_s
          call write_child_parameters(parameter_path,stage_param,stage_dir,current_profile,current_comp, &
               shell_start_s,shell_end_s,inherited_imax,inherited_ni_quad,inherited_output_mode, &
               inherited_dtmax,.true.,.true.,direct_mass_scheme, &
               shell_mass(i),radius_attach+shell_width(i),shell_velocity(i),shell_temp, &
               shell_extension_zones(i),shell_blend_zones(i))
          call run_child(exe_path,stage_param,stage_log,handoff_profile,handoff_comp,handoff_meta,exitstat)
          call require_child_success(exitstat,'shell stage',stage_log)
          call read_stage_state(handoff_meta,state)
       else
          call run_stage_to_radius(parameter_path,exe_path,stage_param,stage_log,stage_dir, &
               current_profile,current_comp,handoff_profile,handoff_comp,handoff_meta, &
               shell_start_s,final_s,shell_inner_radii(i+1),radius_attach+shell_width(i), &
               max(state%max_abs_velocity_cms,1.0d7),inherited_imax,inherited_ni_quad, &
               inherited_output_mode,inherited_dtmax,.true.,.true.,direct_mass_scheme, &
               shell_mass(i),radius_attach+shell_width(i),shell_velocity(i),shell_temp, &
               shell_extension_zones(i),shell_blend_zones(i),shell_end_s,state)
       endif
       call archive_stage_control_files(stage_param,stage_log,stage_dir,archived_log)
       write(*,'(A,I0,A)') ' Stage ',i,' completed.'
       call manager_get_required_double(trim(stage_dir)//'/generated_csm/boundary_diagnostics.txt', &
            'csm_mass_msun',attached_shell_mass(i))
       write(*,'(A,ES11.3,A)') ' Attached CSM mass: ',attached_shell_mass(i),' Msun'
       call append_stage(stage_dirs,stage_offsets,nstages,stage_dir,shell_start_s)
       call require_stage_end(state,shell_end_s-shell_start_s)
       write(stage_label,'(A,I0)') 'stage ',i
       call validate_nonrelativistic_stage(state,trim(stage_label))
       current_profile=handoff_profile
       current_comp=handoff_comp
       current_abs_s=shell_end_s

       if (current_abs_s.ge.final_s-1.0d-6) exit
    enddo

    call stitch_stage_outputs(trim(final_dir),stage_dirs(1:nstages),stage_offsets(1:nstages), &
         staged_observer_nmu,staged_equal_arrival.ne.0)
    call stitch_stage_science_outputs(trim(final_dir),stage_dirs(1:nstages),stage_offsets(1:nstages), &
         staged_observer_nmu,staged_equal_arrival.ne.0)
    call write_staged_manifest(parameter_path,trim(final_dir),n_shells,encounter_days,shell_inner_radii,shell_mass, &
         attached_shell_mass,shell_width,shell_velocity,used_shell_temperature,attachment_radius, &
         shell_extension_zones,shell_blend_zones,initial_age_s,staged_final_day,nstages, &
         stage_dirs(1:nstages),stage_offsets(1:nstages),staged_continuous_shell_evolution_fixed, &
         staged_equal_arrival,staged_observer_nmu,inherited_imax,inherited_ni_quad,inherited_csm_mdot)
    call write_staged_metrics(trim(final_dir),n_shells,encounter_days,staged_final_day)
    call execute_command_line('rm -rf '//shell_quote(trim(scratch_dir)),wait=.true.)

    call print_startup_separator
    write(*,'(A)') ' Simulation completed: staged pipeline finished.'
    write(*,'(A,A)') ' Final products: ',trim(final_dir)
    write(*,'(A,A)') ' Observer bolometric luminosity: ',trim(final_dir)//'/lum_observed.dat'
    write(*,'(A,A)') ' Manifest: ',trim(final_dir)//'/manifest.json'
    call print_startup_separator
  end subroutine run_staged_pipeline


  ! check physical ranges and ordering before accepting a handoff
  subroutine validate_staged_inputs(n_shells,mass_scheme_name,times,masses,widths,velocities, &
       extension_zones,blend_zones,staged_imax,ni_quad,nmu,output_mode)
    integer, intent(in) :: n_shells, extension_zones(n_shells), blend_zones(n_shells)
    character(len=*), intent(in) :: mass_scheme_name
    integer, intent(in) :: staged_imax, ni_quad, nmu, output_mode
    real*8, intent(in) :: times(n_shells), masses(n_shells), widths(n_shells)
    real*8, intent(in) :: velocities(n_shells)
    integer :: i
    logical :: direct_mass

    direct_mass=is_direct_mass_scheme(mass_scheme_name)

    if (staged_imax.lt.20) then
       write(*,*) '******* imax must be >= 20 for staged_pipeline'
       stop 1
    endif
    if (output_mode.lt.0 .or. output_mode.gt.2) then
       write(*,*) '******* output_mode must be 0, 1, or 2'
       stop 1
    endif
    if (ni_quad.lt.1 .or. nmu.lt.2) then
       write(*,*) '******* Ni_quad_npoints >= 1 and staged_observer_nmu >= 2 are required'
       stop 1
    endif
    do i=1,n_shells
       if (times(i).le.0.0d0 .or. widths(i).le.0.0d0) then
          write(*,*) '******* Shell radius and width must be positive. shell=',i
          stop 1
       endif
       if (direct_mass .and. masses(i).le.0.0d0) then
          write(*,*) '******* staged_shell_masses_msun must be positive for direct-mass schemes. shell=',i
          stop 1
       endif
       if (is_wind_mass_scheme(mass_scheme_name)) then
          if (velocities(i).le.0.0d0) then
             write(*,*) '******* Wind-shell velocity must be >0. shell=',i
             stop 1
          endif
       else
          if (velocities(i).lt.0.0d0) then
             write(*,*) '******* Shell velocity must be >=0. shell=',i
             stop 1
          endif
       endif
       if (extension_zones(i).lt.2 .or. blend_zones(i).lt.1 .or. &
            blend_zones(i).ge.extension_zones(i)) then
          write(*,*) '******* Invalid extension/blend zones for shell ',i
          stop 1
       endif
    enddo
    do i=2,n_shells
       if (times(i).le.times(i-1)) then
          write(*,*) '******* staged_shell_inner_radii_cm must be strictly increasing'
          stop 1
       endif
    enddo

  end subroutine validate_staged_inputs


  ! pass the live profile and absolute age to a fresh solver process
  subroutine write_child_parameters(template_path,output_path,stage_outdir,profile_path,comp_path, &
       abs_start_s,abs_end_s,staged_imax,ni_quad,staged_output_mode,dtmax_s, &
       with_shell,restart_from_handoff,direct_mass_control,shell_mass,shell_outer,shell_velocity,shell_temp,extension_zones, &
       blend_zones)
    character(len=*), intent(in) :: template_path,output_path,stage_outdir,profile_path,comp_path
    real*8, intent(in) :: abs_start_s,abs_end_s,dtmax_s
    logical, intent(in) :: with_shell,restart_from_handoff,direct_mass_control
    real*8, intent(in) :: shell_mass,shell_outer,shell_velocity,shell_temp
    integer, intent(in) :: staged_imax,ni_quad,staged_output_mode,extension_zones,blend_zones

    character(len=64) :: keys(max_overrides)
    character(len=path_len) :: vals(max_overrides)
    integer :: n

    call mkdir_p(trim(stage_outdir))

    keys(:)=''; vals(:)=''; n=0
    call add_override(keys,vals,n,'simulation_mode',quoted('single'))
    call add_override(keys,vals,n,'outdir',quoted(trim(stage_outdir)))
    call add_override(keys,vals,n,'profile_name',quoted(trim(profile_path)))
    call add_override(keys,vals,n,'comp_profile_name',quoted(trim(comp_path)))
    call add_override(keys,vals,n,'imax',int_string(staged_imax))
    call add_override(keys,vals,n,'Ni_quad_npoints',int_string(ni_quad))
    call add_override(keys,vals,n,'initial_model_age_s',real_string(abs_start_s))
    call add_override(keys,vals,n,'tend',real_string(abs_end_s-abs_start_s))
    call add_override(keys,vals,n,'dtmax',real_string(dtmax_s))
    call add_override(keys,vals,n,'output_mode',int_string(staged_output_mode))
    if (restart_from_handoff .and. .not.with_shell) then
       call add_override(keys,vals,n,'restart_state_file',quoted(trim(profile_path)//'.state'))
    else
       call add_override(keys,vals,n,'restart_state_file',quoted(''))
    endif
    if (restart_from_handoff) then
       call add_override(keys,vals,n,'Ni_by_hand','0')
       call add_override(keys,vals,n,'boxcar_smoothing','0')
    endif
    if (with_shell) then
       call add_override(keys,vals,n,'csm',quoted('active'))
       if (direct_mass_control) then
          call add_override(keys,vals,n,'csm_mass_msun',real_string(shell_mass))
       else
          call add_override(keys,vals,n,'csm_mass_msun','0.0d0')
       endif
       call add_override(keys,vals,n,'csm_outer_radius_cm',real_string(shell_outer))
       call add_override(keys,vals,n,'csm_wind_velocity_kms',real_string(shell_velocity))
       call add_override(keys,vals,n,'csm_temperature_k',real_string(shell_temp))
       call add_override(keys,vals,n,'csm_extension_zones',int_string(extension_zones))
       call add_override(keys,vals,n,'csm_blend_zones',int_string(blend_zones))
    else
       call add_override(keys,vals,n,'csm',quoted('none'))
       call add_override(keys,vals,n,'mass_scheme',quoted('continuous_mass'))
    endif

    call rewrite_parameter_file(template_path,output_path,keys,vals,n)
  end subroutine write_child_parameters


  subroutine read_profile_surface(profile_path,radius_cm,velocity_cms)
    character(len=*), intent(in) :: profile_path
    real*8, intent(out) :: radius_cm,velocity_cms
    integer :: u,i,n,ios,k
    real*8 :: mass_g,temperature,density,ye,omega
    open(newunit=u,file=trim(profile_path),status='old',action='read',iostat=ios)
    if (ios.ne.0) stop 1
    read(u,*,iostat=ios) n
    if (ios.ne.0.or.n.lt.1) stop 1
    do i=1,n
       read(u,*,iostat=ios) k,mass_g,radius_cm,temperature,density,velocity_cms,ye,omega
       if (ios.ne.0) stop 1
    enddo
    close(u)
  end subroutine read_profile_surface


  ! bracket and bisect the duration needed to reach the next shell
  subroutine run_stage_to_radius(template_path,exe_path,param_path,log_path,stage_dir, &
       profile_path,comp_path,handoff_profile,handoff_comp,handoff_meta, &
       start_s,limit_s,target_radius,start_radius,velocity_guess, &
       imax_value,ni_quad,output_mode,dtmax_s,with_shell,from_handoff,direct_mass, &
       shell_mass,shell_outer,shell_velocity,shell_temp,extension_zones,blend_zones, &
       reached_s,state)
    character(len=*), intent(in) :: template_path,exe_path,param_path,log_path,stage_dir
    character(len=*), intent(in) :: profile_path,comp_path,handoff_profile,handoff_comp,handoff_meta
    real*8, intent(in) :: start_s,limit_s,target_radius,start_radius,velocity_guess
    real*8, intent(in) :: dtmax_s,shell_mass,shell_outer,shell_velocity,shell_temp
    integer, intent(in) :: imax_value,ni_quad,output_mode,extension_zones,blend_zones
    logical, intent(in) :: with_shell,from_handoff,direct_mass
    real*8, intent(out) :: reached_s
    type(stage_state), intent(out) :: state
    real*8 :: low_time,high_time,low_radius,high_radius,candidate,actual_time,observed_speed
    integer :: attempt,exitstat

    if (target_radius.le.start_radius) then
       write(*,*) '******* Shell radius must exceed live exterior radius'
       stop 1
    endif
    low_time=start_s
    low_radius=start_radius
    high_time=-1.0d0
    high_radius=0.0d0
    candidate=min(limit_s,start_s+max(2.0d0,0.95d0*(target_radius-start_radius)/ &
         max(velocity_guess,1.0d7)))
    do attempt=1,18
       if (candidate.le.start_s) then
          write(*,*) '******* No positive time available for shell encounter'
          stop 1
       endif
       call write_child_parameters(template_path,param_path,stage_dir,profile_path,comp_path, &
            start_s,candidate,imax_value,ni_quad,output_mode,dtmax_s, &
            with_shell,from_handoff,direct_mass,shell_mass,shell_outer,shell_velocity, &
            shell_temp,extension_zones,blend_zones)
       call run_child(exe_path,param_path,log_path,handoff_profile,handoff_comp,handoff_meta,exitstat)
       call require_child_success(exitstat,'physical radius encounter',log_path)
       call read_stage_state(handoff_meta,state)
       actual_time=start_s+state%elapsed_time_s
       if (state%radius_cm.ge.target_radius) then
          high_time=actual_time
          high_radius=state%radius_cm
          if (abs(high_radius-target_radius)/target_radius.le.1.0d-3) exit
       else
          low_time=actual_time
          low_radius=state%radius_cm
          if (actual_time.ge.limit_s-1.0d-4) then
             write(*,*) '******* Ejecta did not reach the next physical shell before final age'
             stop 1
          endif
       endif
       if (attempt.eq.18) then
          if (high_time.lt.0.0d0) then
             write(*,*) '******* Shell radius not reached with numerical bracketing'
             stop 1
          endif
          call write_child_parameters(template_path,param_path,stage_dir,profile_path,comp_path, &
               start_s,high_time,imax_value,ni_quad,output_mode,dtmax_s, &
               with_shell,from_handoff,direct_mass,shell_mass,shell_outer,shell_velocity, &
               shell_temp,extension_zones,blend_zones)
          call run_child(exe_path,param_path,log_path,handoff_profile,handoff_comp,handoff_meta,exitstat)
          call require_child_success(exitstat,'physical radius encounter',log_path)
          call read_stage_state(handoff_meta,state)
          exit
       endif
       if (high_time.gt.0.0d0) then
          candidate=0.5d0*(low_time+high_time)
       else
          observed_speed=max((state%radius_cm-start_radius)/max(actual_time-start_s,1.0d-12),1.0d6)
          candidate=min(limit_s,actual_time+max(2.0d0,1.05d0*(target_radius-state%radius_cm)/observed_speed))
       endif
    enddo
    reached_s=start_s+state%elapsed_time_s
    if (reached_s.ge.limit_s-1.0d-5) then
       write(*,*) '******* Shell encounter did not occur before the final simulation age'
       stop 1
    endif
    if (state%radius_cm.lt.target_radius) then
       write(*,*) '******* Physical shell radius not reached'
       stop 1
    endif
    if (abs(state%radius_cm-target_radius)/target_radius.gt.1.0d-3) then
       write(*,'(A,ES11.3)') ' ******* WARNING encounter radius tolerance exceeded: ', &
            abs(state%radius_cm-target_radius)/target_radius
    endif
    write(*,'(A,ES12.4,A,F10.5,A)') ' Encounter radius: ',state%radius_cm, &
         ' cm, simulation age ',reached_s/day_seconds,' d'
  end subroutine run_stage_to_radius


  ! retain the solver exit status while recording its output
  subroutine run_child(exe_path,param_path,log_path,handoff_profile,handoff_comp,handoff_meta,exitstat)
    character(len=*), intent(in) :: exe_path,param_path,log_path,handoff_profile,handoff_comp,handoff_meta
    integer, intent(out) :: exitstat
    character(len=line_len*4) :: cmd
    character(len=line_len*3) :: inner_cmd
    logical :: exists_profile, exists_comp, exists_meta

    inner_cmd = 'SNEC_PARAMETER_FILE='//shell_quote(trim(param_path))//' '// &
         'STRATA_STAGE_CHILD=1 GFORTRAN_UNBUFFERED_ALL=1 '// &
         'STRATA_HANDOFF_PROFILE='//shell_quote(trim(handoff_profile))//' '// &
         'STRATA_HANDOFF_COMP='//shell_quote(trim(handoff_comp))//' '// &
         'STRATA_HANDOFF_META='//shell_quote(trim(handoff_meta))//' '// &
         shell_quote(trim(exe_path))//' 2>&1 | tee '//shell_quote(trim(log_path))
    cmd='bash -o pipefail -c '//shell_quote(trim(inner_cmd))

    call execute_command_line(trim(cmd),wait=.true.,exitstat=exitstat)
    inquire(file=trim(handoff_profile),exist=exists_profile)
    inquire(file=trim(handoff_comp),exist=exists_comp)
    inquire(file=trim(handoff_meta),exist=exists_meta)
    if (exitstat.eq.0 .and. (.not.exists_profile .or. .not.exists_comp .or. .not.exists_meta)) then
       exitstat=999
    endif
  end subroutine run_child


  subroutine archive_stage_control_files(temp_param,temp_log,stage_dir,archived_log)
    character(len=*), intent(in) :: temp_param,temp_log,stage_dir
    character(len=*), intent(out) :: archived_log
    integer :: stat
    logical :: exists

    archived_log = trim(stage_dir)//'/stage.log'
    call mkdir_p(trim(stage_dir))

    inquire(file=trim(temp_log),exist=exists)
    if (exists) then
       call execute_command_line('mv -f '//shell_quote(trim(temp_log))//' '// &
            shell_quote(trim(archived_log)),wait=.true.,exitstat=stat)
       if (stat.ne.0) then
          write(*,*) '******* Could not archive staged log: ',trim(temp_log)
          stop 1
       endif
    endif

    inquire(file=trim(temp_param),exist=exists)
    if (exists) then
       call execute_command_line('rm -f '//shell_quote(trim(temp_param)),wait=.true.)
    endif
  end subroutine archive_stage_control_files


  subroutine require_child_success(exitstat,label,log_path)
    integer, intent(in) :: exitstat
    character(len=*), intent(in) :: label,log_path
    if (exitstat.ne.0) then
       write(*,*) '******* Staged child failed: ',trim(label)
       write(*,*) '******* See log: ',trim(log_path)
       stop 1
    endif
  end subroutine require_child_success


  subroutine append_stage(stage_dirs,stage_offsets,nstages,stage_dir,offset)
    character(len=path_len), intent(inout) :: stage_dirs(:)
    real*8, intent(inout) :: stage_offsets(:)
    integer, intent(inout) :: nstages
    character(len=*), intent(in) :: stage_dir
    real*8, intent(in) :: offset
    nstages=nstages+1
    if (nstages.gt.size(stage_dirs)) then
       write(*,*) '******* simulation_manager: stage array overflow'
       stop 1
    endif
    stage_dirs(nstages)=trim(stage_dir)
    stage_offsets(nstages)=offset
  end subroutine append_stage


  subroutine parent_directory(path,parent)
    character(len=*), intent(in) :: path
    character(len=*), intent(out) :: parent
    integer :: i
    parent='.'
    do i=len_trim(path),1,-1
       if (path(i:i).eq.'/') then
          if (i.gt.1) parent=path(1:i-1)
          return
       endif
    enddo
  end subroutine parent_directory


  subroutine stage_handoff_paths(stage_dir,profile_path,comp_path,meta_path)
    character(len=*), intent(in) :: stage_dir
    character(len=*), intent(out) :: profile_path,comp_path,meta_path
    character(len=path_len) :: local_dir

    local_dir=trim(stage_dir)//'/profiles'
    profile_path=trim(local_dir)//'/handoff_profile.dat'
    comp_path=trim(local_dir)//'/handoff_comp.dat'
    meta_path=trim(local_dir)//'/handoff_state.dat'
  end subroutine stage_handoff_paths


  ! read the recorded end state before accepting another encounter
  subroutine read_stage_state(path,state)
    character(len=*), intent(in) :: path
    type(stage_state), intent(out) :: state
    character(len=line_len) :: line,key,value
    integer :: u,ios,p

    state=stage_state()
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write(*,*) '******* Cannot read staged handoff metadata: ',trim(path)
       stop 1
    endif
    do
       read(u,'(A)',iostat=ios) line
       if (ios.ne.0) exit
       p=index(line,'=')
       if (p.le.0) cycle
       key=trim(adjustl(line(1:p-1)))
       value=trim(adjustl(line(p+1:)))
       select case(trim(key))
       case('elapsed_time_s'); read(value,*) state%elapsed_time_s
       case('termination_reason'); state%termination_reason=trim(value)
       case('mass_g'); read(value,*) state%mass_g
       case('radius_cm'); read(value,*) state%radius_cm
       case('rho_outer'); read(value,*) state%rho_outer
       case('temp_outer'); read(value,*) state%temp_outer
       case('lum_grey_selected'); read(value,*) state%lum_grey_selected
       case('max_abs_velocity_cms'); read(value,*) state%max_abs_velocity_cms
       end select
    enddo
    close(u)
    if (state%radius_cm.le.0.0d0) then
       write(*,*) '******* Invalid handoff radius in ',trim(path)
       stop 1
    endif
  end subroutine read_stage_state


  subroutine require_stage_end(state,expected_elapsed_s)
    type(stage_state), intent(in) :: state
    real*8, intent(in) :: expected_elapsed_s
    if (.not.ieee_is_finite(state%elapsed_time_s)) error stop 'Non-finite stage end time'
    if (state%termination_reason.ne.'tend' .or. &
         abs(state%elapsed_time_s-expected_elapsed_s).gt.max(1.0d-6,1.0d-12*expected_elapsed_s)) &
         error stop 'Stage did not reach the requested encounter age'
  end subroutine require_stage_end


  ! stop the handoff if material reaches the configured newtonian speed limit
  subroutine validate_nonrelativistic_stage(state,label)
    type(stage_state), intent(in) :: state
    character(len=*), intent(in) :: label
    real*8 :: beta

    if (.not.ieee_is_finite(state%max_abs_velocity_cms) .or. &
         state%max_abs_velocity_cms.lt.0.0d0) then
       write(*,*) '******* Invalid staged maximum velocity in ',trim(label)
       stop 1
    endif

    beta=state%max_abs_velocity_cms/c_light
    if (beta.ge.0.30d0) then
       write(*,'(A,A,A,F7.3,A)') ' ******* ',trim(label), &
            ' reaches max |v|/c = ',beta,'.'
       write(*,'(A)') '         STRATA-NOVA/SNEC is non-relativistic; this staged state is not valid.'
       write(*,'(A)') '         Increase the physical CSM mass loading/width or choose a denser wind before continuing.'
       stop 1
    else if (beta.ge.0.20d0) then
       write(*,'(A,A,A,F7.3,A)') ' WARNING: ',trim(label), &
            ' reaches max |v|/c = ',beta,'. Keep this run out of the final science set unless convergence is demonstrated.'
    endif
  end subroutine validate_nonrelativistic_stage


  ! estimate the shell temperature from the incoming luminosity and radius
  real*8 function equilibrium_shell_temperature(luminosity,radius) result(temp_shell)
    real*8, intent(in) :: luminosity,radius
    if (luminosity.gt.0.0d0 .and. radius.gt.0.0d0) then
       temp_shell=(luminosity/(4.0d0*pi_mgr*sigma_sb*radius*radius))**0.25d0
       temp_shell=max(1000.0d0,temp_shell)
    else
       temp_shell=1000.0d0
    endif
  end function equilibrium_shell_temperature


  ! join stage-local luminosities on the absolute time axis
  subroutine stitch_stage_outputs(workdir,stage_dirs,stage_offsets,nmu,apply_eat)
    character(len=*), intent(in) :: workdir
    character(len=path_len), intent(in) :: stage_dirs(:)
    real*8, intent(in) :: stage_offsets(:)
    integer, intent(in) :: nmu
    logical, intent(in) :: apply_eat

    integer :: j,n,total_l,total_r,pos_l,pos_r,nloc
    real*8, allocatable :: t(:),lum(:),tr(:),rad(:),local_t(:),local_y(:)
    real*8, allocatable :: rr(:),observer(:)
    character(len=path_len) :: path

    total_l=0; total_r=0
    do j=1,size(stage_dirs)
       path=trim(stage_dirs(j))//'/lum_observed.dat'
       call count_two_column_rows(path,nloc)
       total_l=total_l+nloc
       path=trim(stage_dirs(j))//'/rad_photo.dat'
       call count_two_column_rows(path,nloc)
       total_r=total_r+nloc
    enddo
    if (total_l.lt.2 .or. total_r.lt.2) then
       write(*,*) '******* Not enough scalar rows to stitch staged calculation.'
       stop 1
    endif
    allocate(t(total_l),lum(total_l),tr(total_r),rad(total_r))
    pos_l=0; pos_r=0

    do j=1,size(stage_dirs)
       path=trim(stage_dirs(j))//'/lum_observed.dat'
       call read_two_column_file(path,local_t,local_y)
       call append_scalar_stage(local_t,local_y,stage_offsets(j),j.gt.1,t,lum,pos_l)
       deallocate(local_t,local_y)
       path=trim(stage_dirs(j))//'/rad_photo.dat'
       call read_two_column_file(path,local_t,local_y)
       call append_scalar_stage(local_t,local_y,stage_offsets(j),j.gt.1,tr,rad,pos_r)
       deallocate(local_t,local_y)
    enddo

    call compact_monotonic_series(t,lum,pos_l,n)
    call compact_monotonic_series(tr,rad,pos_r,nloc)
    allocate(rr(n),observer(n))
    do j=1,n
       rr(j)=linear_interp(tr(1:nloc),rad(1:nloc),t(j))
    enddo
    if (apply_eat) then
       call equal_arrival_time_curve(t(1:n),lum(1:n),rr,n,nmu,observer)
    else
       observer(:)=lum(1:n)
    endif

    call write_two_column(trim(workdir)//'/rad_photo_combined.dat',t(1:n),rr,n)
    call write_two_column(trim(workdir)//'/rad_photo.dat',t(1:n),rr,n)
    call write_two_column(trim(workdir)//'/lum_observed.dat',t(1:n),observer,n)

    deallocate(t,lum,tr,rad,rr,observer)
  end subroutine stitch_stage_outputs


  ! join the available scalar and radial outputs across handoffs
  subroutine stitch_stage_science_outputs(final_dir,stage_dirs,stage_offsets,nmu,apply_eat)
    character(len=*), intent(in) :: final_dir
    character(len=path_len), intent(in) :: stage_dirs(:)
    real*8, intent(in) :: stage_offsets(:)
    integer, intent(in) :: nmu
    logical, intent(in) :: apply_eat

    integer :: i
    character(len=64), parameter :: scalar_files(10) = [character(len=64) :: &
         'vel_photo.dat', 'T_eff.dat', 'index_photo.dat', &
         'photometry_valid.dat', 'temperature_floor_fraction.dat', &
         'mass_photo.dat', 'opacity_corrupted.dat', 'index_lumshell.dat', &
         'mass_lumshell.dat', 'hydro_runtime_diagnostics.dat']
    character(len=64), parameter :: xg_files(27) = [character(len=64) :: &
         'vel.xg', 'rho.xg', 'ye.xg', 'press.xg', 'cs2.xg', 'Q.xg', &
         'eps.xg', 'mass.xg', 'temp.xg', 'lum.xg', 'tau.xg', &
         'delta_time.xg', 'radius.xg', 'kappa.xg', 'kappa_table.xg', &
         'logR_op.xg', 'logT.xg', 'p_rad.xg', 'Ni_deposit_function.xg', &
         'He_1.xg', 'He_2.xg', 'He_3.xg', 'H_1.xg', 'H_2.xg', &
         'free_electron_frac.xg', 'E_shell.xg', 'photosphere_tracer.xg']

    do i=1,size(scalar_files)
       call stitch_generic_timeseries(trim(final_dir),trim(scalar_files(i)), &
            stage_dirs,stage_offsets)
    enddo
    do i=1,size(xg_files)
       call stitch_xg_profiles(trim(final_dir),trim(xg_files(i)),stage_dirs,stage_offsets)
    enddo
    call stitch_stage_magnitudes(trim(final_dir),stage_dirs,stage_offsets,nmu,apply_eat)
  end subroutine stitch_stage_science_outputs


  ! apply observer delays to linear passband flux before restoring magnitudes
  subroutine stitch_stage_magnitudes(final_dir,stage_dirs,stage_offsets,nmu,apply_eat)
    character(len=*), intent(in) :: final_dir
    character(len=path_len), intent(in) :: stage_dirs(:)
    real*8, intent(in) :: stage_offsets(:)
    integer, intent(in) :: nmu
    logical, intent(in) :: apply_eat

    integer :: j,i,k,u,ios,nrows,n,pos
    real*8 :: vals(13),abs_t
    real*8, allocatable :: t(:),teff(:),mags(:,:),rad_t(:),rad_y(:),rr(:)
    real*8, allocatable :: flux(:),flux_obs(:)
    character(len=path_len) :: path,outpath
    logical :: exists

    nrows=0
    do j=1,size(stage_dirs)
       path=trim(stage_dirs(j))//'/magnitudes.dat'
       inquire(file=trim(path),exist=exists)
       if (.not.exists) cycle
       open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
       if (ios.ne.0) cycle
       do
          read(u,*,iostat=ios) vals
          if (ios.ne.0) exit
          if (j.gt.1 .and. vals(1).le.1.0d-9) cycle
          nrows=nrows+1
       enddo
       close(u)
    enddo
    if (nrows.lt.2) return

    allocate(t(nrows),teff(nrows),mags(nrows,11)); pos=0
    do j=1,size(stage_dirs)
       path=trim(stage_dirs(j))//'/magnitudes.dat'
       inquire(file=trim(path),exist=exists)
       if (.not.exists) cycle
       open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
       if (ios.ne.0) cycle
       do
          read(u,*,iostat=ios) vals
          if (ios.ne.0) exit
          if (j.gt.1 .and. vals(1).le.1.0d-9) cycle
          abs_t=vals(1)+stage_offsets(j)
          if (pos.gt.0) then
             if (abs_t.le.t(pos)+1.0d-9) cycle
          endif
          pos=pos+1
          t(pos)=abs_t; teff(pos)=vals(2); mags(pos,1:11)=vals(3:13)
       enddo
       close(u)
    enddo
    n=pos
    if (n.lt.2) then
       deallocate(t,teff,mags); return
    endif

    call read_two_column_file(trim(final_dir)//'/rad_photo.dat',rad_t,rad_y)
    allocate(rr(n),flux(n),flux_obs(n))
    do i=1,n
       rr(i)=linear_interp(rad_t,rad_y,t(i))
    enddo

    if (apply_eat) then
       do k=1,11
          do i=1,n
             if (isfinite_real(mags(i,k)) .and. abs(mags(i,k)).lt.100.0d0) then
                flux(i)=10.0d0**(-0.4d0*mags(i,k))
             else
                flux(i)=1.0d-99
             endif
          enddo
          call equal_arrival_time_curve(t(1:n),flux,rr,n,nmu,flux_obs)
          do i=1,n
             mags(i,k)=-2.5d0*log10(max(flux_obs(i),1.0d-99))
          enddo
       enddo
    endif

    outpath=trim(final_dir)//'/magnitudes.dat'
    open(newunit=u,file=trim(outpath),status='replace',action='write')
    do i=1,n
       write(u,'(13E18.9)') t(i),teff(i),mags(i,1:11)
    enddo
    close(u)

    deallocate(t,teff,mags,rad_t,rad_y,rr,flux,flux_obs)
  end subroutine stitch_stage_magnitudes


  subroutine stitch_generic_timeseries(final_dir,filename,stage_dirs,stage_offsets)
    character(len=*), intent(in) :: final_dir,filename
    character(len=path_len), intent(in) :: stage_dirs(:)
    real*8, intent(in) :: stage_offsets(:)

    integer :: j,u_in,u_out,ios,p0,p1,stat
    real*8 :: local_t,abs_t
    character(len=line_len) :: line,remainder
    character(len=path_len) :: input_path,output_path,temp_path
    logical :: exists,any_input

    any_input=.false.
    do j=1,size(stage_dirs)
       input_path=trim(stage_dirs(j))//'/'//trim(filename)
       inquire(file=trim(input_path),exist=exists)
       if (exists) any_input=.true.
    enddo
    if (.not.any_input) return

    output_path=trim(final_dir)//'/'//trim(filename)
    temp_path=trim(final_dir)//'/.combined_'//trim(filename)
    open(newunit=u_out,file=trim(temp_path),status='replace',action='write')

    do j=1,size(stage_dirs)
       input_path=trim(stage_dirs(j))//'/'//trim(filename)
       inquire(file=trim(input_path),exist=exists)
       if (.not.exists) cycle
       open(newunit=u_in,file=trim(input_path),status='old',action='read',iostat=ios)
       if (ios.ne.0) cycle
       do
          read(u_in,'(A)',iostat=ios) line
          if (ios.ne.0) exit
          if (len_trim(line).eq.0) cycle
          read(line,*,iostat=ios) local_t
          if (ios.ne.0) then
             if (j.eq.1) write(u_out,'(A)') trim(line)
             cycle
          endif
          if (j.gt.1 .and. local_t.le.1.0d-9) cycle
          abs_t=local_t+stage_offsets(j)
          call first_token_bounds(line,p0,p1)
          remainder=''
          if (p1.lt.len_trim(line)) remainder=line(p1+1:len_trim(line))
          if (len_trim(remainder).gt.0) then
             write(u_out,'(ES24.14E3,1X,A)') abs_t,trim(adjustl(remainder))
          else
             write(u_out,'(ES24.14E3)') abs_t
          endif
       enddo
       close(u_in)
    enddo
    close(u_out)
    call execute_command_line('mv -f '//shell_quote(trim(temp_path))//' '// &
         shell_quote(trim(output_path)),wait=.true.,exitstat=stat)
    if (stat.ne.0) then
       write(*,*) '******* Could not finalize staged scalar product: ',trim(output_path)
       stop 1
    endif
  end subroutine stitch_generic_timeseries


  subroutine stitch_xg_profiles(final_dir,filename,stage_dirs,stage_offsets)
    character(len=*), intent(in) :: final_dir,filename
    character(len=path_len), intent(in) :: stage_dirs(:)
    real*8, intent(in) :: stage_offsets(:)

    integer :: j,u_in,u_out,ios,p,stat
    real*8 :: local_t,abs_t
    character(len=line_len) :: line
    character(len=path_len) :: input_path,output_path,temp_path
    logical :: exists,any_input,skip_block

    any_input=.false.
    do j=1,size(stage_dirs)
       input_path=trim(stage_dirs(j))//'/'//trim(filename)
       inquire(file=trim(input_path),exist=exists)
       if (exists) any_input=.true.
    enddo
    if (.not.any_input) return

    output_path=trim(final_dir)//'/'//trim(filename)
    temp_path=trim(final_dir)//'/.combined_'//trim(filename)
    open(newunit=u_out,file=trim(temp_path),status='replace',action='write')

    do j=1,size(stage_dirs)
       input_path=trim(stage_dirs(j))//'/'//trim(filename)
       inquire(file=trim(input_path),exist=exists)
       if (.not.exists) cycle
       open(newunit=u_in,file=trim(input_path),status='old',action='read',iostat=ios)
       if (ios.ne.0) cycle
       skip_block=.false.
       do
          read(u_in,'(A)',iostat=ios) line
          if (ios.ne.0) exit
          p=index(line,'Time =')
          if (p.gt.0) then
             read(line(p+6:),*,iostat=ios) local_t
             if (ios.ne.0) cycle
             skip_block=(j.gt.1 .and. local_t.le.1.0d-9)
             if (.not.skip_block) then
                abs_t=local_t+stage_offsets(j)
                write(u_out,'(A,ES24.14E3)') '"Time = ',abs_t
             endif
          else if (.not.skip_block) then
             write(u_out,'(A)') trim(line)
          endif
       enddo
       close(u_in)
    enddo
    close(u_out)
    call execute_command_line('mv -f '//shell_quote(trim(temp_path))//' '// &
         shell_quote(trim(output_path)),wait=.true.,exitstat=stat)
    if (stat.ne.0) then
       write(*,*) '******* Could not finalize staged XG product: ',trim(output_path)
       stop 1
    endif
  end subroutine stitch_xg_profiles


  subroutine first_token_bounds(line,p0,p1)
    character(len=*), intent(in) :: line
    integer, intent(out) :: p0,p1
    integer :: i,n
    n=len_trim(line); p0=1; p1=n
    do i=1,n
       if (line(i:i).ne.' ' .and. iachar(line(i:i)).ne.9) then
          p0=i; exit
       endif
    enddo
    p1=p0
    do i=p0,n
       if (line(i:i).eq.' ' .or. iachar(line(i:i)).eq.9) then
          p1=i-1; return
       endif
       p1=i
    enddo
  end subroutine first_token_bounds


  subroutine append_scalar_stage(local_t,local_y,offset,restarted,t,y,pos)
    real*8, intent(in) :: local_t(:),local_y(:),offset
    logical, intent(in) :: restarted
    real*8, intent(inout) :: t(:),y(:)
    integer, intent(inout) :: pos
    integer :: i
    do i=1,size(local_t)
       if (restarted .and. local_t(i).le.1.0d-9) cycle
       if (.not.isfinite_real(local_t(i)) .or. .not.isfinite_real(local_y(i))) cycle
       pos=pos+1
       if (pos.gt.size(t)) then
          write(*,*) '******* simulation_manager: scalar stitch allocation overflow'
          stop 1
       endif
       t(pos)=local_t(i)+offset
       y(pos)=local_y(i)
    enddo
  end subroutine append_scalar_stage


  ! remove duplicate times while preserving the preceding stage at a handoff
  subroutine compact_monotonic_series(t,y,nin,nout)
    real*8, intent(inout) :: t(:),y(:)
    integer, intent(in) :: nin
    integer, intent(out) :: nout
    integer :: i
    nout=0
    do i=1,nin
       if (nout.eq.0) then
          nout=1; t(nout)=t(i); y(nout)=y(i)
       else if (t(i).gt.t(nout)+1.0d-9) then
          nout=nout+1; t(nout)=t(i); y(nout)=y(i)
       else if (abs(t(i)-t(nout)).le.1.0d-9) then
          cycle
       else
          write(*,*) '******* Non-monotonic staged scalar times encountered.'
          stop 1
       endif
    enddo
  end subroutine compact_monotonic_series


  ! measure discrete extrema in the configured encounter windows
  subroutine write_staged_metrics(workdir,n_shells,encounter_days,final_day)
    character(len=*), intent(in) :: workdir
    integer, intent(in) :: n_shells
    real*8, intent(in) :: encounter_days(n_shells),final_day
    real*8, allocatable :: t(:),l(:)
    real*8 :: lo,hi,peak_t,peak_l,valley_t,valley_l,ratio,main_t,main_l
    integer :: i,u

    call read_two_column_file(trim(workdir)//'/lum_observed.dat',t,l)
    open(newunit=u,file=trim(workdir)//'/staged_metrics.dat',status='replace',action='write')
    write(u,'(A)') '# interaction_shell peak_time_day peak_lum_erg_s valley_time_day valley_lum_erg_s peak_to_valley'
    do i=1,n_shells
       lo=encounter_days(i)*day_seconds
       if (i.lt.n_shells) then
          hi=min((encounter_days(i)+0.75d0)*day_seconds,encounter_days(i+1)*day_seconds)
       else
          hi=min((encounter_days(i)+0.75d0)*day_seconds,final_day*day_seconds)
       endif
       call window_extremum(t,l,lo,hi,.true.,peak_t,peak_l)
       if (i.lt.n_shells) then
          hi=min((encounter_days(i)+2.0d0)*day_seconds,encounter_days(i+1)*day_seconds)
       else
          hi=min((encounter_days(i)+2.0d0)*day_seconds,final_day*day_seconds)
       endif
       call window_extremum(t,l,peak_t,hi,.false.,valley_t,valley_l)
       if (valley_l.gt.0.0d0) then
          ratio=peak_l/valley_l
       else
          ratio=0.0d0
       endif
       write(u,'(I6,1X,5ES20.10E3)') i,peak_t/day_seconds,peak_l, &
            valley_t/day_seconds,valley_l,ratio
    enddo
    lo=min((encounter_days(n_shells)+2.0d0)*day_seconds,final_day*day_seconds)
    hi=final_day*day_seconds
    call window_extremum(t,l,lo,hi,.true.,main_t,main_l)
    write(u,'(A,1X,2ES20.10E3)') '# radioactive_maximum',main_t/day_seconds,main_l
    close(u)
    deallocate(t,l)
  end subroutine write_staged_metrics


  subroutine window_extremum(t,l,lo,hi,want_max,te,le)
    real*8, intent(in) :: t(:),l(:),lo,hi
    logical, intent(in) :: want_max
    real*8, intent(out) :: te,le
    integer :: i
    logical :: found
    found=.false.; te=lo; le=0.0d0
    do i=1,size(t)
       if (t(i).lt.lo .or. t(i).gt.hi) cycle
       if (.not.found) then
          te=t(i); le=l(i); found=.true.
       else if ((want_max .and. l(i).gt.le) .or. ((.not.want_max) .and. l(i).lt.le)) then
          te=t(i); le=l(i)
       endif
    enddo
  end subroutine window_extremum


  subroutine write_staged_manifest(parameter_path,workdir,n_shells,times,inner_radii,masses,attached_masses,widths,velocities, &
       temps,radii,extension_zones,blend_zones,initial_age_s,final_day,nstages, &
       stage_dirs,stage_offsets,continuous_shell_evolution,eat_enabled,nmu,staged_imax,ni_quad,inherited_csm_mdot)
    character(len=*), intent(in) :: parameter_path,workdir
    integer, intent(in) :: n_shells,nstages,extension_zones(n_shells),blend_zones(n_shells)
    integer, intent(in) :: continuous_shell_evolution,eat_enabled,nmu,staged_imax,ni_quad
    real*8, intent(in) :: times(n_shells),inner_radii(n_shells),masses(n_shells),attached_masses(n_shells)
    real*8, intent(in) :: widths(n_shells),velocities(n_shells)
    real*8, intent(in) :: temps(n_shells),radii(n_shells)
    real*8, intent(in) :: initial_age_s,final_day,stage_offsets(nstages),inherited_csm_mdot
    character(len=path_len), intent(in) :: stage_dirs(nstages)
    integer :: u,i
    character(len=32) :: manifest_mass_scheme

    call manager_get_string(parameter_path,'mass_scheme',manifest_mass_scheme,'constant_density')
    open(newunit=u,file=trim(workdir)//'/manifest.json',status='replace',action='write')
    write(u,'(A)') '{'
    write(u,'(A,ES24.14E3,A)') '  "initial_age_s": ',initial_age_s,','
    write(u,'(A,ES24.14E3,A)') '  "final_day": ',final_day,','
    write(u,'(A,I0,A)') '  "n_shells": ',n_shells,','
    write(u,'(A,A,A)') '  "mass_scheme": "',trim(json_escape(manifest_mass_scheme)),'",'
    if (is_direct_mass_scheme(manifest_mass_scheme)) then
       write(u,'(A)') '  "mass_control": "staged_shell_masses_msun",'
    else if (trim(adjustl(manifest_mass_scheme)).eq.'constant_mass_loss') then
       write(u,'(A)') '  "mass_control": "csm_mdot_msun_per_yr",'
    else
       write(u,'(A)') '  "mass_control": "non_constant_mass_loss",'
    endif
    if (trim(adjustl(manifest_mass_scheme)).eq.'constant_mass_loss') then
       write(u,'(A,ES24.14E3,A)') '  "csm_mdot_msun_per_yr": ',inherited_csm_mdot,','
    endif
    write(u,'(A,I0,A)') '  "staged_imax": ',staged_imax,','
    write(u,'(A,I0,A)') '  "ni_quad_npoints": ',ni_quad,','
    write(u,'(A,I0,A)') '  "observer_nmu": ',nmu,','
    write(u,'(A,I0,A)') '  "equal_arrival_time": ',eat_enabled,','
    write(u,'(A)') '  "luminosity_default": "grey_bolometric",'
    write(u,'(A)') '  "luminosity_file": "lum_observed.dat",'
    write(u,'(A)') '  "equal_arrival_time_channels": ["lum_observed.dat", "magnitudes.dat"],'
    write(u,'(A,I0,A)') '  "continuous_shell_evolution": ',continuous_shell_evolution,','
    write(u,'(A)') '  "shells": ['
    do i=1,n_shells
       write(u,'(A)') '    {'
       write(u,'(A,I0,A)') '      "index": ',i,','
       write(u,'(A,ES24.14E3,A)') '      "encounter_time_day": ',times(i),','
       write(u,'(A,ES24.14E3,A)') '      "configured_inner_radius_cm": ',inner_radii(i),','
       write(u,'(A,ES24.14E3,A)') '      "mass_msun": ',attached_masses(i),','
       if (is_direct_mass_scheme(manifest_mass_scheme)) then
          write(u,'(A,ES24.14E3,A)') '      "configured_mass_msun": ',masses(i),','
       endif
       write(u,'(A,ES24.14E3,A)') '      "width_cm": ',widths(i),','
       write(u,'(A,ES24.14E3,A)') '      "velocity_kms": ',velocities(i),','
       write(u,'(A,ES24.14E3,A)') '      "temperature_k": ',temps(i),','
       write(u,'(A,ES24.14E3,A)') '      "attachment_radius_cm": ',radii(i),','
       write(u,'(A,I0,A)') '      "extension_zones": ',extension_zones(i),','
       write(u,'(A,I0)') '      "blend_zones": ',blend_zones(i)
       if (i.lt.n_shells) then
          write(u,'(A)') '    },'
       else
          write(u,'(A)') '    }'
       endif
    enddo
    write(u,'(A)') '  ],'
    write(u,'(A)') '  "stages": ['
    do i=1,nstages
       write(u,'(A,A,A,ES24.14E3,A)',advance='no') '    {"path": "', &
            trim(json_escape(stage_dirs(i))), '", "offset_s": ',stage_offsets(i),'}'
       if (i.lt.nstages) then
          write(u,'(A)') ','
       else
          write(u,*)
       endif
    enddo
    write(u,'(A)') '  ]'
    write(u,'(A)') '}'
    close(u)

  end subroutine write_staged_manifest


  ! use the median valid radius-to-velocity ratio as the initial age
  subroutine infer_homologous_age(profile_path,age_s)
    character(len=*), intent(in) :: profile_path
    real*8, intent(out) :: age_s
    integer :: u,ios,n,i,nvalid,idx
    character(len=line_len) :: line
    real*8 :: m,r,t,rho,v,ye,omega
    real*8, allocatable :: ages(:)

    open(newunit=u,file=trim(profile_path),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write(*,*) '******* Cannot open base profile for age inference: ',trim(profile_path)
       stop 1
    endif
    read(u,*,iostat=ios) n
    if (ios.ne.0 .or. n.lt.1) then
       write(*,*) '******* simulation_manager: invalid base profile header'
       stop 1
    endif
    allocate(ages(n)); nvalid=0
    do i=1,n
       read(u,'(A)',iostat=ios) line
       if (ios.ne.0) exit
       read(line,*,iostat=ios) idx,m,r,t,rho,v,ye,omega
       if (ios.ne.0) cycle
       if (v.gt.1.0d6 .and. r.gt.0.0d0) then
          nvalid=nvalid+1
          ages(nvalid)=r/v
       endif
    enddo
    close(u)
    if (nvalid.lt.1) then
       age_s=0.0d0
       deallocate(ages)
       return
    endif
    call quicksort_real(ages,1,nvalid)
    if (mod(nvalid,2).eq.1) then
       age_s=ages((nvalid+1)/2)
    else
       age_s=0.5d0*(ages(nvalid/2)+ages(nvalid/2+1))
    endif
    deallocate(ages)
  end subroutine infer_homologous_age


  recursive subroutine quicksort_real(a,left,right)
    real*8, intent(inout) :: a(:)
    integer, intent(in) :: left,right
    integer :: i,j
    real*8 :: pivot,tmp
    if (left.ge.right) return
    i=left; j=right; pivot=a((left+right)/2)
    do
       do while(a(i).lt.pivot); i=i+1; enddo
       do while(a(j).gt.pivot); j=j-1; enddo
       if (i.le.j) then
          tmp=a(i); a(i)=a(j); a(j)=tmp
          i=i+1; j=j-1
       endif
       if (i.gt.j) exit
    enddo
    if (left.lt.j) call quicksort_real(a,left,j)
    if (i.lt.right) call quicksort_real(a,i,right)
  end subroutine quicksort_real


  logical function is_direct_mass_scheme(name) result(is_direct)
    character(len=*), intent(in) :: name
    character(len=32) :: normalized
    normalized=trim(adjustl(name))
    is_direct=(normalized.eq.'continuous_mass' .or. normalized.eq.'constant_density')
  end function is_direct_mass_scheme


  logical function is_wind_mass_scheme(name) result(is_wind)
    character(len=*), intent(in) :: name
    character(len=32) :: normalized
    normalized=trim(adjustl(name))
    is_wind=(normalized.eq.'non_constant_mass_loss' .or. normalized.eq.'constant_mass_loss')
  end function is_wind_mass_scheme


  subroutine prepare_workdir(workdir,protected_paths)
    character(len=*), intent(in) :: workdir,protected_paths(:)
    call check_output_directory(workdir,protected_paths)
    call mkdir_p(workdir)
    call reset_output_directory(workdir,protected_paths)
  end subroutine prepare_workdir


  subroutine mkdir_p(path)
    character(len=*), intent(in) :: path
    integer :: stat
    call execute_command_line('mkdir -p -- '//shell_quote(trim(path)),wait=.true.,exitstat=stat)
    if (stat.ne.0) then
       write(*,*) '******* Could not create directory: ',trim(path)
       stop 1
    endif
  end subroutine mkdir_p


  ! replace exact parameter keys while retaining the remaining settings
  subroutine rewrite_parameter_file(template_path,output_path,keys,vals,nkeys)
    character(len=*), intent(in) :: template_path,output_path
    character(len=64), intent(in) :: keys(:)
    character(len=path_len), intent(in) :: vals(:)
    integer, intent(in) :: nkeys
    logical :: seen(max_overrides)
    character(len=line_len) :: line,key
    integer :: ui,uo,ios,p,j

    seen=.false.
    open(newunit=ui,file=trim(template_path),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write(*,*) '******* simulation_manager: cannot open parameter template: ',trim(template_path)
       stop 1
    endif
    open(newunit=uo,file=trim(output_path),status='replace',action='write')
    do
       read(ui,'(A)',iostat=ios) line
       if (ios.ne.0) exit
       key=''; p=index(line,'=')
       if (p.gt.0) then
          if (first_nonblank_char(line).ne.'#' .and. first_nonblank_char(line).ne.'!') then
             key=trim(adjustl(line(1:p-1)))
          endif
       endif
       j=find_override(keys,nkeys,trim(key))
       if (j.gt.0) then
          if (.not.seen(j)) then
             write(uo,'(A,A,A)') trim(keys(j)),' = ',trim(vals(j))
             seen(j)=.true.
          endif
       else
          write(uo,'(A)') trim(line)
       endif
    enddo
    do j=1,nkeys
       if (.not.seen(j)) write(uo,'(A,A,A)') trim(keys(j)),' = ',trim(vals(j))
    enddo
    close(ui); close(uo)
  end subroutine rewrite_parameter_file


  integer function find_override(keys,nkeys,key) result(idx)
    character(len=64), intent(in) :: keys(:)
    integer, intent(in) :: nkeys
    character(len=*), intent(in) :: key
    integer :: i
    idx=0
    if (len_trim(key).eq.0) return
    do i=1,nkeys
       if (trim(keys(i)).eq.trim(key)) then
          idx=i; return
       endif
    enddo
  end function find_override


  subroutine add_override(keys,vals,n,key,val)
    character(len=64), intent(inout) :: keys(:)
    character(len=path_len), intent(inout) :: vals(:)
    integer, intent(inout) :: n
    character(len=*), intent(in) :: key,val
    n=n+1
    if (n.gt.size(keys)) then
       write(*,*) '******* simulation_manager: too many parameter overrides'
       stop 1
    endif
    keys(n)=trim(key); vals(n)=trim(val)
  end subroutine add_override


  function first_nonblank_char(line) result(ch)
    character(len=*), intent(in) :: line
    character(len=1) :: ch
    integer :: i
    ch=' '
    do i=1,len_trim(line)
       if (line(i:i).ne.' ' .and. iachar(line(i:i)).ne.9) then
          ch=line(i:i); return
       endif
    enddo
  end function first_nonblank_char


  ! check that each control has the shape required by its simulation mode
  subroutine validate_mode_parameter_shapes(path,mode)
    character(len=*), intent(in) :: path,mode
    integer :: n_shells,nvals
    character(len=64) :: scheme
    logical :: found

    call require_scalar_if_present(path,'csm')
    call require_scalar_if_present(path,'csm_wind_velocity_kms')
    call require_scalar_if_present(path,'csm_temperature_k')

    select case (trim(mode))
    case ('single')
       call require_at_most_one_if_present(path,'staged_shell_inner_radii_cm')
       call require_at_most_one_if_present(path,'staged_shell_masses_msun')
       call require_at_most_one_if_present(path,'staged_shell_widths_cm')
    case ('staged_pipeline')
       call manager_get_integer(path,'staged_n_shells',n_shells,-1)
       if (n_shells.lt.1) then
          write(*,*) '******* staged_pipeline requires staged_n_shells >= 1'
          stop 1
       endif
       call require_exact_array_count(path,'staged_shell_inner_radii_cm',n_shells,.true.)
       call require_exact_array_count(path,'staged_shell_widths_cm',n_shells,.true.)
       call manager_get_string(path,'mass_scheme',scheme,'constant_density')
       call lowercase_inplace(scheme)
       if (is_direct_mass_scheme(scheme)) then
          call require_exact_array_count(path,'staged_shell_masses_msun',n_shells,.true.)
       else
          call count_parameter_values(path,'staged_shell_masses_msun',nvals,found)
          if (found .and. nvals.ne.n_shells) then
             call report_array_count_error('staged_shell_masses_msun',n_shells,nvals)
          endif
       endif
    case default
    end select
  end subroutine validate_mode_parameter_shapes


  subroutine require_scalar_if_present(path,key)
    character(len=*), intent(in) :: path,key
    integer :: nvals
    logical :: found
    call count_parameter_values(path,key,nvals,found)
    if (found .and. nvals.ne.1) then
       write(*,'(A,A,A,I0,A)') ' ******* ',trim(key),' must contain exactly 1 value; found ',nvals,'.'
       stop 1
    endif
  end subroutine require_scalar_if_present


  subroutine require_at_most_one_if_present(path,key)
    character(len=*), intent(in) :: path,key
    integer :: nvals
    logical :: found
    call count_parameter_values(path,key,nvals,found)
    if (found .and. nvals.gt.1) then
       write(*,'(A,A,A,I0,A)') ' ******* single mode uses 1 shell: ',trim(key), &
            ' must contain at most 1 value; found ',nvals,'.'
       stop 1
    endif
  end subroutine require_at_most_one_if_present


  subroutine require_exact_array_count(path,key,nrequired,required)
    character(len=*), intent(in) :: path,key
    integer, intent(in) :: nrequired
    logical, intent(in) :: required
    integer :: nvals
    logical :: found
    call count_parameter_values(path,key,nvals,found)
    if (.not.found) then
       if (required) then
          write(*,'(A,A,A,I0,A)') ' ******* ',trim(key),' requires ',nrequired, &
               ' values for the requested staged shells; found 0.'
          stop 1
       endif
       return
    endif
    if (nvals.ne.nrequired) call report_array_count_error(key,nrequired,nvals)
  end subroutine require_exact_array_count


  subroutine report_array_count_error(key,nrequired,nfound)
    character(len=*), intent(in) :: key
    integer, intent(in) :: nrequired,nfound
    write(*,'(A,A,A,I0,A,I0,A)') ' ******* ',trim(key),' requires ',nrequired, &
         ' values for ',nrequired,' staged shells; found '//trim(int_string(nfound))//'.'
    stop 1
  end subroutine report_array_count_error


  subroutine count_parameter_values(path,key,nvals,found)
    character(len=*), intent(in) :: path,key
    integer, intent(out) :: nvals
    logical, intent(out) :: found
    character(len=line_len) :: raw
    integer :: i
    logical :: in_token

    call manager_get_raw(path,key,raw,found)
    nvals=0
    if (.not.found) return
    in_token=.false.
    do i=1,len_trim(raw)
       if (raw(i:i).eq.' ' .or. iachar(raw(i:i)).eq.9 .or. raw(i:i).eq.',') then
          in_token=.false.
       else
          if (.not.in_token) nvals=nvals+1
          in_token=.true.
       endif
    enddo
  end subroutine count_parameter_values


  ! match a full parameter name before extracting its value
  subroutine manager_get_raw(path,key,value,found)
    character(len=*), intent(in) :: path,key
    character(len=*), intent(out) :: value
    logical, intent(out) :: found
    character(len=line_len) :: line,lkey
    integer :: u,ios,p,c,i
    character :: quote

    found=.false.; value=''
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if (ios.ne.0) return
    do
       read(u,'(A)',iostat=ios) line
       if (ios.ne.0) exit
       if (first_nonblank_char(line).eq.'#' .or. first_nonblank_char(line).eq.'!') cycle
       if (len_trim(line).eq.len(line)) error stop 'Parameter line is too long'
       p=index(line,'=')
       if (p.le.0) cycle
       lkey=trim(adjustl(line(1:p-1)))
       if (trim(lkey).ne.trim(key)) cycle
       quote=' '; c=len_trim(line)+1
       do i=p+1,len_trim(line)
          if (quote.ne.' ') then
             if (line(i:i).eq.quote) quote=' '
          else if (line(i:i).eq.'"' .or. line(i:i).eq."'") then
             quote=line(i:i)
          else if (line(i:i).eq.'#' .or. line(i:i).eq.'!') then
             c=i
             exit
          endif
       enddo
       if (quote.ne.' ') error stop 'Unclosed quote in parameter file'
       lkey=trim(adjustl(line(p+1:c-1)))
       if (len_trim(lkey).gt.len(value)) error stop 'Parameter value is too long'
       value=trim(lkey)
       call strip_quotes(value)
       found=.true.; exit
    enddo
    close(u)
  end subroutine manager_get_raw


  subroutine manager_get_string(path,key,value,default)
    character(len=*), intent(in) :: path,key,default
    character(len=*), intent(out) :: value
    logical :: found
    call manager_get_raw(path,key,value,found)
    if (.not.found) value=default
  end subroutine manager_get_string


  subroutine manager_get_double(path,key,value,default)
    character(len=*), intent(in) :: path,key
    real*8, intent(out) :: value
    real*8, intent(in) :: default
    character(len=line_len) :: raw
    logical :: found
    integer :: ios
    call manager_get_raw(path,key,raw,found)
    if (.not.found) then
       value=default; return
    endif
    read(raw,*,iostat=ios) value
    if (ios.ne.0) then
       write(*,*) '******* Could not parse staged real parameter: ',trim(key)
       stop 1
    endif
    if (.not.ieee_is_finite(value)) error stop 'Non-finite staged real parameter'
  end subroutine manager_get_double


  subroutine manager_get_required_double(path,key,value)
    character(len=*), intent(in) :: path,key
    real*8, intent(out) :: value
    character(len=line_len) :: raw
    logical :: found
    integer :: ios
    call manager_get_raw(path,key,raw,found)
    if (.not.found) then
       write(*,*) '******* Required staged diagnostic is missing: ',trim(key)
       write(*,*) '        File: ',trim(path)
       stop 1
    endif
    read(raw,*,iostat=ios) value
    if (ios.ne.0 .or. value.ne.value .or. value.le.0.0d0) then
       write(*,*) '******* Invalid staged diagnostic value: ',trim(key)
       write(*,*) '        File: ',trim(path)
       stop 1
    endif
    if (.not.ieee_is_finite(value)) error stop 'Non-finite handoff value'
  end subroutine manager_get_required_double


  subroutine manager_get_integer(path,key,value,default)
    character(len=*), intent(in) :: path,key
    integer, intent(out) :: value
    integer, intent(in) :: default
    character(len=line_len) :: raw
    logical :: found
    integer :: ios
    call manager_get_raw(path,key,raw,found)
    if (.not.found) then
       value=default; return
    endif
    read(raw,*,iostat=ios) value
    if (ios.ne.0) then
       write(*,*) '******* Could not parse staged integer parameter: ',trim(key)
       stop 1
    endif
  end subroutine manager_get_integer


  subroutine manager_get_real_array(path,key,n,a)
    character(len=*), intent(in) :: path,key
    integer, intent(in) :: n
    real*8, intent(out) :: a(n)
    character(len=line_len) :: raw
    logical :: found
    integer :: ios
    call manager_get_raw(path,key,raw,found)
    if (.not.found) then
       write(*,*) '******* Missing staged array parameter: ',trim(key)
       stop 1
    endif
    call commas_to_spaces(raw)
    read(raw,*,iostat=ios) a
    if (ios.ne.0) then
       write(*,*) '******* Could not parse ',n,' values for ',trim(key)
       stop 1
    endif
    if (any(.not.ieee_is_finite(a))) error stop 'Non-finite staged array'
  end subroutine manager_get_real_array


  subroutine commas_to_spaces(s)
    character(len=*), intent(inout) :: s
    integer :: i
    do i=1,len(s)
       if (s(i:i).eq.',') s(i:i)=' '
    enddo
  end subroutine commas_to_spaces


  subroutine strip_quotes(s)
    character(len=*), intent(inout) :: s
    integer :: n
    s=trim(adjustl(s)); n=len_trim(s)
    if (n.ge.2) then
       if ((s(1:1).eq.'"' .and. s(n:n).eq.'"') .or. &
            (s(1:1).eq."'" .and. s(n:n).eq."'")) then
          s=s(2:n-1)
       endif
    endif
  end subroutine strip_quotes


  subroutine lowercase_inplace(s)
    character(len=*), intent(inout) :: s
    integer :: i,c
    do i=1,len_trim(s)
       c=iachar(s(i:i))
       if (c.ge.iachar('A') .and. c.le.iachar('Z')) s(i:i)=achar(c+32)
    enddo
  end subroutine lowercase_inplace


  function quoted(s) result(out)
    character(len=*), intent(in) :: s
    character(len=path_len) :: out
    out='"'//trim(s)//'"'
  end function quoted


  function real_string(x) result(s)
    real*8, intent(in) :: x
    character(len=64) :: s
    write(s,'(ES26.16E3)') x
  end function real_string


  function int_string(i) result(s)
    integer, intent(in) :: i
    character(len=64) :: s
    write(s,'(I0)') i
  end function int_string


  function json_escape(s) result(out)
    character(len=*), intent(in) :: s
    character(len=path_len*2) :: out
    integer :: i,p
    out=''; p=0
    do i=1,len_trim(s)
       if (s(i:i).eq.'"' .or. s(i:i).eq.'\\') then
          p=p+1; out(p:p)='\\'
       endif
       p=p+1; out(p:p)=s(i:i)
    enddo
  end function json_escape


  logical function isfinite_real(x) result(ok)
    real*8, intent(in) :: x
    ok=ieee_is_finite(x)
  end function isfinite_real


  subroutine count_two_column_rows(path,n)
    character(len=*), intent(in) :: path
    integer, intent(out) :: n
    integer :: u,ios
    real*8 :: a,b
    n=0
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write(*,*) '******* Missing staged scalar file: ',trim(path)
       stop 1
    endif
    do
       read(u,*,iostat=ios) a,b
       if (ios.ne.0) exit
       n=n+1
    enddo
    close(u)
  end subroutine count_two_column_rows


  subroutine read_two_column_file(path,x,y)
    character(len=*), intent(in) :: path
    real*8, allocatable, intent(out) :: x(:),y(:)
    integer :: n,u,ios,i
    call count_two_column_rows(path,n)
    allocate(x(n),y(n))
    open(newunit=u,file=trim(path),status='old',action='read')
    do i=1,n
       read(u,*,iostat=ios) x(i),y(i)
       if (ios.ne.0) then
          write(*,*) '******* simulation_manager: malformed two-column file: ',trim(path)
          stop 1
       endif
    enddo
    close(u)
  end subroutine read_two_column_file


  subroutine write_two_column(path,x,y,n)
    character(len=*), intent(in) :: path
    real*8, intent(in) :: x(n),y(n)
    integer, intent(in) :: n
    integer :: u,i
    open(newunit=u,file=trim(path),status='replace',action='write')
    do i=1,n
       write(u,'(2ES26.16E3)') x(i),y(i)
    enddo
    close(u)
  end subroutine write_two_column


end module simulation_manager
