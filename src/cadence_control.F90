! control the pace of data output 

module cadence_control

  use parameters, only: dtout_fast_target, dtout_mid_target, dtout_late_target, &
       dtout_scalar_fast_target, dtout_scalar_mid_target, dtout_scalar_late_target, &
       output_mid_transition_days, output_late_transition_days, &
       Ni_period_min_target, Ni_period_max_target, Ni_fractional_change_target, initial_model_age_s
  implicit none

  real*8, parameter :: day_seconds = 86400.0d0

contains

  ! logic block for switching output cadence when the code transitions from early to mid to late simulation time
  real*8 function select_cadence(current_time, fast_value, mid_value, slow_value) result(cadence)
    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(in) :: fast_value, mid_value, slow_value
    real*8 :: time_days
    real*8 :: mid_transition, late_transition

    time_days = max(0.0d0, current_time) / day_seconds
    mid_transition = max(0.0d0, output_mid_transition_days)
    late_transition = max(mid_transition, output_late_transition_days)

    ! slow, medium and fast data output
    if (time_days < mid_transition) then
       cadence = fast_value ! early cadence
    else if (time_days < late_transition) then
       cadence = mid_value ! later cadence
    else
       cadence = slow_value ! latest cadence with slowest data output
    endif
  end function select_cadence


  ! rounds to nearest whole number of needed timesteps
  real*8 function quantize_to_step(target, step) result(period)

    implicit none
    real*8, intent(in) :: target
    real*8, intent(in) :: step
    integer(8) :: nsteps

    if (step <= 0.0d0) then
       period = target
       return
    endif

    if (target <= 0.0d0) then
       period = step
       return
    endif

    nsteps = max(1_8, nint(target/step, 8))
    period = dble(nsteps) * step

  end function quantize_to_step

  ! computes spatial profiles at the requested cadence
  real*8 function compute_profile_period(current_time, step) result(period)

    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(in) :: step
    real*8 :: target

    target = select_cadence(current_time, dtout_fast_target, dtout_mid_target, dtout_late_target)
    period = quantize_to_step(target, step)

  end function compute_profile_period

  ! similar concept here; computes the interval for scalar outputs (like luminosity.dat)
  real*8 function compute_scalar_period(current_time, step) result(period)

    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(in) :: step
    real*8 :: target

    target = select_cadence(current_time, dtout_scalar_fast_target, dtout_scalar_mid_target, &
         dtout_scalar_late_target)
    period = quantize_to_step(target, step)

  end function compute_scalar_period


  ! estimating derivative of nickel heating (in time)
  ! focusing on one component of the ni-56 decay chain, being ni-56 to co-56
  real*8 function compute_radioactive_heating_timescale(current_time) result(t_heat)

    implicit none
    real*8, intent(in) :: current_time
    real*8, parameter :: tau_Ni_days = 8.8d0
    real*8, parameter :: tau_Co_days = 111.3d0
    real*8, parameter :: epsilon_Ni = 3.9d10
    real*8, parameter :: epsilon_Co = 6.78d9
    real*8 :: tau_Ni, tau_Co
    real*8 :: exp_Ni, exp_Co
    real*8 :: heating, dheating_dt

    tau_Ni = tau_Ni_days * day_seconds
    tau_Co = tau_Co_days * day_seconds

    exp_Ni = exp(-max(0.0d0, current_time+initial_model_age_s)/tau_Ni)
    exp_Co = exp(-max(0.0d0, current_time+initial_model_age_s)/tau_Co)

    heating = epsilon_Ni*exp_Ni + epsilon_Co*(exp_Co - exp_Ni)
    dheating_dt = -epsilon_Ni*exp_Ni/tau_Ni + epsilon_Co*(-exp_Co/tau_Co + exp_Ni/tau_Ni)

    if (heating <= 0.0d0 .or. abs(dheating_dt) <= tiny(1.0d0)) then
       t_heat = Ni_period_max_target
    else
       t_heat = min(Ni_period_max_target, abs(heating/dheating_dt))
    endif

  end function compute_radioactive_heating_timescale


  ! interval between ni heating calculations; akin to quantize_step
  real*8 function compute_Ni_period(current_time, step) result(period)
    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(in) :: step
    real*8 :: t_heat, target
    real*8 :: min_period

    t_heat = compute_radioactive_heating_timescale(current_time)
    target = Ni_fractional_change_target * t_heat
    min_period = Ni_period_min_target
    target = max(min_period, min(target, Ni_period_max_target))
    period = quantize_to_step(target, step)

  end function compute_Ni_period


  ! final calculation function making use of quantize_to_step for integer intervals
  subroutine compute_adaptive_periods(current_time, step, profile_period, scalar_period, Ni_period_out)

    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(in) :: step
    real*8, intent(out) :: profile_period
    real*8, intent(out) :: scalar_period
    real*8, intent(out) :: Ni_period_out

    profile_period = compute_profile_period(current_time, step)
    scalar_period = compute_scalar_period(current_time, step)
    Ni_period_out = compute_Ni_period(current_time, step)

  end subroutine compute_adaptive_periods


end module cadence_control
