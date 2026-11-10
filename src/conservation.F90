subroutine conservation_update_sources(step_dt)

  use blmod, only: Ni_energy_rate, Ni_deposit_function, delta_mass, lum, &
       cumulative_ni_deposited_energy, cumulative_boundary_radiated_energy, &
       energy_audit_prev_ni_power, energy_audit_prev_rad_power, &
       energy_audit_source_initialized, diag_rad_energy_residual_erg, &
       cumulative_rad_equation_residual
  use blmod, only: radiation_work_boundary_power, cumulative_boundary_work_energy
  use parameters, only: imax, Ni_switch
  implicit none

  real*8, intent(in) :: step_dt
  real*8 :: ni_power, rad_power

  ni_power = 0.0d0
  if (Ni_switch.eq.1) then
     ni_power = Ni_energy_rate * sum(Ni_deposit_function(1:imax)*delta_mass(1:imax))
  endif
  rad_power = max(lum(imax),0.0d0)

  if (.not.energy_audit_source_initialized) then
     energy_audit_prev_ni_power = ni_power
     energy_audit_prev_rad_power = rad_power
     energy_audit_source_initialized = .true.
  endif

  ! accumulate ni deposition and escaping radiation with trapezoidal time weights
  cumulative_ni_deposited_energy = cumulative_ni_deposited_energy + &
       0.5d0*(energy_audit_prev_ni_power + ni_power)*step_dt
  cumulative_boundary_radiated_energy = cumulative_boundary_radiated_energy + &
       0.5d0*(energy_audit_prev_rad_power + rad_power)*step_dt

  cumulative_boundary_work_energy = cumulative_boundary_work_energy + radiation_work_boundary_power*step_dt
  energy_audit_prev_ni_power = ni_power
  energy_audit_prev_rad_power = rad_power
  cumulative_rad_equation_residual = cumulative_rad_equation_residual + &
       diag_rad_energy_residual_erg

end subroutine conservation_update_sources

subroutine conservation_compute_energies

  use blmod, only: delta_mass, delta_cmass, cmass, cr, gravity_switch, eps, vel, &
                    total_initial_energy, time, tdump_scalar, &
                    cumulative_ni_deposited_energy, cumulative_boundary_radiated_energy, &
                    energy_budget_reference, energy_budget_reference_set
  use blmod, only: cumulative_boundary_work_energy
  use parameters
  use physical_constants
  implicit none

  integer :: i

  real*8 :: egrav, eint, ekin
  real*8 :: legacy_egrav, legacy_eint, legacy_ekin
  real*8 :: source_corrected, budget_residual, audit_start_time


  egrav = 0.0d0
  eint = 0.0d0
  ekin = 0.0d0

  ! sum material energy over physical cells, excluding the outer sentinel
  do i=1,imax-1
     egrav = egrav - ggrav*delta_mass(i)*cmass(i)/cr(i) *gravity_switch
     eint = eint + eps(i)*delta_mass(i)
  enddo

  ! use the same nodal masses as the momentum update
  ekin = 0.25d0*delta_mass(1)*vel(1)**2 + &
       0.5d0*sum(delta_cmass(1:imax-1)*vel(2:imax)**2)

  if(time.eq.0.0d0) then
      ! keep the inherited initial-energy scale used by the radiation solve
      legacy_egrav = 0.0d0
      legacy_eint = 0.0d0
      legacy_ekin = 0.0d0
      do i=1,imax
         legacy_egrav = legacy_egrav - ggrav*delta_mass(i)*cmass(i)/cr(i) *gravity_switch
         legacy_eint = legacy_eint + eps(i)*delta_mass(i)
      enddo
      do i=1,imax-1
         legacy_ekin = legacy_ekin + 0.5d0*(0.50d0*(vel(i+1)+vel(i)))**2 * delta_mass(i)
      enddo
      legacy_ekin = legacy_ekin + 0.5d0*vel(imax)**2 * delta_mass(imax)
      total_initial_energy = legacy_egrav+legacy_eint+legacy_ekin
  endif

  audit_start_time = max(piston_tend,bomb_tend)
  ! remove deposited ni energy and add escaped radiation and boundary work
  source_corrected = egrav+eint+ekin-total_initial_energy &
       - cumulative_ni_deposited_energy + cumulative_boundary_radiated_energy + cumulative_boundary_work_energy
  if (.not.energy_budget_reference_set .and. time.ge.audit_start_time) then
     energy_budget_reference = source_corrected
     energy_budget_reference_set = .true.
  endif
  if (energy_budget_reference_set) then
     budget_residual = source_corrected-energy_budget_reference
  else
     budget_residual = 0.0d0
  endif

  if(output_mode .eq. 0 .and. (time.eq.0.0d0.or.time.ge.tdump_scalar)) then
     open(666,file=trim(adjustl(outdir))//"/conservation.dat",&
          status='unknown',position='append')
     write(666,"(1P10E18.9)") time,egrav,eint,ekin,egrav+eint+ekin, &
         egrav+eint+ekin-total_initial_energy
     close(666)
     open(667,file=trim(adjustl(outdir))//"/conservation_budget.dat", &
          status='unknown',position='append')
     write(667,"(1P8E20.10)") time,egrav+eint+ekin, &
          egrav+eint+ekin-total_initial_energy, cumulative_ni_deposited_energy, &
          cumulative_boundary_radiated_energy, source_corrected, budget_residual, cumulative_boundary_work_energy
     close(667)
  endif

end subroutine conservation_compute_energies
