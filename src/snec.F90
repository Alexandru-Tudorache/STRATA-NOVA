program snec

  use blmod, only: dtime, dtime_p, time, nt, ntstart, tstart,   &
     tdump, tdump_scalar, tdump_check, time_Ni, &
     expanded_input_detected, &
     diag_hydro_rad_iters, diag_scratch_count, &
     diag_analysis_called, diag_ni_eval_total, diag_ni_eval_remap, &
     diag_scratch_reason, diag_scratch_status, diag_scratch_location
  use parameters
  use cadence_control, only: compute_adaptive_periods
  use csm_switch, only: prepare_csm_inputs
  use startup_log, only: print_title_card, print_startup_separator, is_stage_child
  use simulation_manager, only: simulation_manager_dispatch, &
       simulation_manager_write_child_handoff
  implicit none

  logical :: OutputFlag = .false.
  logical :: OutputFlagScalar = .false.
  logical :: OutputFlagCheck = .false.
  integer :: timing_unit = 198
  integer(8) :: t_count_start, t_count_end, t_count_rate
  real*8 :: blstep_walltime
  logical :: manager_handled
  logical :: record_step_timing


  call print_title_card

  ! dispatch staged runs before the ordinary solver setup
  call simulation_manager_dispatch(manager_handled)
  if (manager_handled) stop

  call input_parser
  call prepare_csm_inputs

  call problem

  call artificial_viscosity

  call output_all(0)
  call output_all(1)
  ! skip the uninitialized observer sample for imported ejecta
  if (.not.expanded_input_detected) call output_all(2)

  call timestep

  time = tstart
  nt = ntstart

  call compute_adaptive_periods(time, dtime, dtout, dtout_scalar, Ni_period)

  tdump_check = tstart+dtout_check
  tdump_scalar = time+dtout_scalar
  tdump = time+dtout

  record_step_timing = output_mode .eq. 0 .and. grid_debug .ge. 1
  if (record_step_timing) then
     call system_clock(count_rate=t_count_rate)
     open(unit=timing_unit, file=trim(adjustl(outdir))//'/blstep_timing.dat', &
          status='replace')
     write(timing_unit,'(A)') '# nt time_s time_d dtime blstep_ms iters scratches analysis ' // &
          'scratch_reason scratch_status scratch_zone'
  endif


  ! advance one accepted hydro step at a time
  IntegrationLoop: do

     dtime_p = dtime
     call timestep

     call compute_adaptive_periods(time, dtime, dtout, dtout_scalar, Ni_period)

     if (time.lt.tdump) then
        tdump = min(tdump, time + dtout)
     endif
     if (time.lt.tdump_scalar) then
        tdump_scalar = min(tdump_scalar, time + dtout_scalar)
     endif
     if (time.lt.time_Ni) then
        time_Ni = min(time_Ni, time + Ni_period)
     endif


     if(ntinfo.gt.0) then
        if(mod(nt,int(ntinfo,int64)).eq.0) then
           call outinfo
        endif
     endif

     if((time+dtime).gt.tend) dtime = tend-time

     if (record_step_timing) call system_clock(t_count_start)
     ! retry rejected trials inside blstep before accumulating sources
     call blstep
     call conservation_update_sources(dtime)
     if (record_step_timing) then
        call system_clock(t_count_end)
        blstep_walltime = dble(t_count_end - t_count_start) / dble(t_count_rate) * 1.0d3
        write(timing_unit,'(I8,3ES16.7,F12.4,2I5,L3,3I6)') &
            nt, time, time/86400.0d0, dtime, blstep_walltime, &
            diag_hydro_rad_iters, diag_scratch_count, &
            diag_analysis_called, diag_scratch_reason, &
            diag_scratch_status, diag_scratch_location
     endif

     ! advance counters only after the accepted update
     nt = nt + 1

     if (ntout.gt.0) then
        if ( mod(nt,int(ntout,int64)) .eq. 0) OutputFlag = .true.
     endif

     if (ntout_scalar.gt.0) then
        if ( mod(nt,int(ntout_scalar,int64)) .eq. 0 ) OutputFlagScalar = .true.
     endif

     if (ntout_check.gt.0) then
        if ( mod(nt,int(ntout_check,int64)) .eq. 0 ) OutputFlagCheck = .true.
     endif

     if ( time.ge.tdump) then
        tdump=time+dtout
        OutputFlag = .true.
     endif

     if ( time.ge.tdump_scalar) then
        tdump_scalar=time+dtout_scalar
        OutputFlagScalar = .true.
     endif

     if ( time.ge.tdump_check) then
        tdump_check=tdump_check+dtout_check
        OutputFlagCheck = .true.
     endif

     time = time+dtime

     if (OutputFlag) then
        call output_all(0)
        call output_all(1)
        OutputFlag = .false.
     endif

     if (OutputFlagScalar) then
        call output_all(2)
        OutputFlagScalar = .false.
     endif

     if (OutputFlagCheck) then
        OutputFlagCheck = .false.
     endif

     if (time.eq.tend) then
        if (.not.is_stage_child()) then
           call print_startup_separator
           write(*,'(A,I8,A)') " Simulation completed: tend reached after ", nt, " timesteps"
           write(*,'(A,I6,A,I6,A)') " Ni evals: ", diag_ni_eval_total, &
                "  (remap-forced resets: ", diag_ni_eval_remap, ")"
        endif
        call output_all(0)
        call output_all(1)
        call output_all(2)
        call write_run_summary(nt)
        exit
     endif

     if (nt.ge.ntmax) then
        if (.not.is_stage_child()) then
           call print_startup_separator
           write(*,*) "ERROR: simulation incomplete; ntmax reached before tend"
           write(*,'(A,I6,A,I6,A)') " Ni evals: ", diag_ni_eval_total, &
                "  (remap-forced resets: ", diag_ni_eval_remap, ")"
        endif
        call output_all(0)
        call output_all(1)
        call output_all(2)
        call write_run_summary(nt)
        stop 1
     endif

  enddo IntegrationLoop

  if (record_step_timing) close(timing_unit)

  ! export the accepted live state when this process is a stage child
  call simulation_manager_write_child_handoff()

end program snec
