module csm_shared

  use filesystem_utils, only: shell_quote
  use parameters, only: csm_surface_mix_zones
  implicit none

  real*8, parameter :: msun_cgs = 1.98847d33
  real*8, parameter :: year_cgs = 3.15576d7
  real*8, parameter :: pi_cgs = 3.14159265358979323846d0
  real*8, parameter :: floor_abundance = 1.0d-99

  type :: csm_profile_data
     integer :: n_zones = 0
     integer, allocatable :: idx(:)
     real*8, allocatable :: mass(:)
     real*8, allocatable :: radius(:)
     real*8, allocatable :: temp(:)
     real*8, allocatable :: rho(:)
     real*8, allocatable :: vel(:)
     real*8, allocatable :: ye(:)
     real*8, allocatable :: omega(:)
  end type csm_profile_data

  type :: csm_composition_data
     integer :: n_zones = 0
     integer :: n_isotopes = 0
     real*8, allocatable :: mass_numbers(:)
     real*8, allocatable :: charge_numbers(:)
     character(len=16), allocatable :: species(:)
     real*8, allocatable :: mass(:)
     real*8, allocatable :: radius(:)
     real*8, allocatable :: fractions(:,:)
  end type csm_composition_data

  type :: csm_boundary_diagnostics
     real*8 :: mass_loss_rate = 0.0d0
     integer :: extension_zones = 0
     real*8 :: attachment_radius_cm = 0.0d0
     real*8 :: outer_radius_cm = 0.0d0
     real*8 :: csm_mass_msun = 0.0d0
     real*8 :: surface_mass_g = 0.0d0
     real*8 :: surface_radius_cm = 0.0d0
     real*8 :: surface_density = 0.0d0
     real*8 :: first_csm_density = 0.0d0
     real*8 :: analytic_first_csm_density = 0.0d0
     real*8 :: density_jump_ratio = 1.0d0
     real*8 :: analytic_density_jump_ratio = 1.0d0
     integer :: transition_anchor_index = 0
     real*8 :: transition_anchor_radius_cm = 0.0d0
     real*8 :: transition_anchor_density = 0.0d0
     real*8 :: bridged_surface_density = 0.0d0
     real*8 :: bridge_density_jump_ratio = 1.0d0
     character(len=32) :: density_smoothing_mode = 'none'
     integer :: density_transition_zones = 0
     real*8 :: density_inner_log_slope = 0.0d0
     real*8 :: density_outer_log_slope = 0.0d0
     real*8 :: surface_temperature = 0.0d0
     real*8 :: first_csm_temperature = 0.0d0
     real*8 :: temperature_jump_ratio = 1.0d0
     real*8 :: surface_velocity_cgs = 0.0d0
     real*8 :: first_csm_velocity_cgs = 0.0d0
     real*8 :: outer_csm_velocity_cgs = 0.0d0
     real*8 :: requested_wind_velocity_kms = 0.0d0
     integer :: velocity_transition_zones = 0
     character(len=32) :: velocity_smoothing_mode = 'none'
     real*8 :: wind_ye = 0.0d0
     real*8 :: invariant_min = 0.0d0
     real*8 :: invariant_max = 0.0d0
     character(len=32) :: mass_scheme = 'continuous_mass'
     real*8 :: mass_loss_duration_yr = 0.0d0
     real*8 :: mass_loss_rate_start = 0.0d0
     real*8 :: mass_loss_rate_burst = 0.0d0
     real*8 :: outer_density_normalization = 0.0d0
     logical :: terminal_surface_repaired = .false.
     integer :: terminal_surface_anchor_index = 0
     real*8 :: terminal_density_original = 0.0d0
     real*8 :: terminal_temperature_original = 0.0d0
  end type csm_boundary_diagnostics

contains

  subroutine read_profile_file(path, profile)
    character(len=*), intent(in) :: path
    type(csm_profile_data), intent(out) :: profile
    integer :: i, u

    u = 901
    open(unit=u, file=trim(path), status='old', form='formatted', action='read')
    read(u,*) profile%n_zones
    allocate(profile%idx(profile%n_zones))
    allocate(profile%mass(profile%n_zones))
    allocate(profile%radius(profile%n_zones))
    allocate(profile%temp(profile%n_zones))
    allocate(profile%rho(profile%n_zones))
    allocate(profile%vel(profile%n_zones))
    allocate(profile%ye(profile%n_zones))
    allocate(profile%omega(profile%n_zones))
    do i = 1, profile%n_zones
       read(u,*) profile%idx(i), profile%mass(i), profile%radius(i), &
            profile%temp(i), profile%rho(i), profile%vel(i), &
            profile%ye(i), profile%omega(i)
    enddo
    close(u)
  end subroutine read_profile_file

  subroutine read_composition_file(path, composition)
    character(len=*), intent(in) :: path
    type(csm_composition_data), intent(out) :: composition
    integer :: i, u

    u = 902
    open(unit=u, file=trim(path), status='old', form='formatted', action='read')
    read(u,*) composition%n_zones, composition%n_isotopes
    allocate(composition%mass_numbers(composition%n_isotopes))
    allocate(composition%charge_numbers(composition%n_isotopes))
    allocate(composition%species(composition%n_isotopes))
    allocate(composition%mass(composition%n_zones))
    allocate(composition%radius(composition%n_zones))
    allocate(composition%fractions(composition%n_zones, composition%n_isotopes))
    read(u,*) composition%mass_numbers(1:composition%n_isotopes)
    read(u,*) composition%charge_numbers(1:composition%n_isotopes)
    do i = 1, composition%n_isotopes
       composition%species(i) = species_name(int(nint(composition%mass_numbers(i))), &
            int(nint(composition%charge_numbers(i))))
    enddo
    do i = 1, composition%n_zones
       read(u,*) composition%mass(i), composition%radius(i), &
            composition%fractions(i,1:composition%n_isotopes)
    enddo
    close(u)
  end subroutine read_composition_file

  subroutine validate_alignment(profile, composition)
    type(csm_profile_data), intent(in) :: profile
    type(csm_composition_data), intent(in) :: composition
    integer :: i

    if (profile%n_zones.ne.composition%n_zones) then
       write(*,*) '******* Profile and composition zone counts do not match'
       stop 1
    endif
    do i = 1, profile%n_zones
       if (.not.values_close(profile%mass(i), composition%mass(i), &
            1.0d-6, 1.0d20)) then
          write(*,*) '******* Profile and composition mass coordinates do not match'
          stop 1
       endif
       if (.not.values_close(profile%radius(i), composition%radius(i), &
            1.0d-6, 1.0d2)) then
          write(*,*) '******* Profile and composition radius coordinates do not match'
          stop 1
       endif
    enddo
  end subroutine validate_alignment

  subroutine write_profile_file(path, profile)
    character(len=*), intent(in) :: path
    type(csm_profile_data), intent(in) :: profile
    integer :: i, u

    u = 903
    open(unit=u, file=trim(path), status='replace', form='formatted', &
         action='write')
    write(u,'(I8)') profile%n_zones
    do i = 1, profile%n_zones
       write(u,'(I8,1X,7(ES24.16E3,1X))') i, profile%mass(i), profile%radius(i), &
            profile%temp(i), profile%rho(i), profile%vel(i), profile%ye(i), &
            profile%omega(i)
    enddo
    close(u)
  end subroutine write_profile_file

  subroutine write_composition_file(path, composition)
    character(len=*), intent(in) :: path
    type(csm_composition_data), intent(in) :: composition
    integer :: i, u

    u = 904
    open(unit=u, file=trim(path), status='replace', form='formatted', &
         action='write')
    write(u,'(I8,1X,I8)') composition%n_zones, composition%n_isotopes
    write(u,'(*(ES14.6E3,1X))') composition%mass_numbers(:)
    write(u,'(*(ES14.6E3,1X))') composition%charge_numbers(:)
    do i = 1, composition%n_zones
       write(u,'(2ES18.10E3,1X,*(ES14.6E3,1X))') composition%mass(i), &
            composition%radius(i), composition%fractions(i,1:composition%n_isotopes)
    enddo
    close(u)
  end subroutine write_composition_file

  ! append the chosen csm mixture without changing the ejecta species order
  subroutine extend_composition(base_composition, new_profile, mix_aligned, &
       new_composition)
    type(csm_composition_data), intent(in) :: base_composition
    type(csm_profile_data), intent(in) :: new_profile
    real*8, intent(in) :: mix_aligned(:)
    type(csm_composition_data), intent(out) :: new_composition
    integer :: i, j, copied_zones, extension_zones, max_copy

    max_copy = min(base_composition%n_zones, new_profile%n_zones)
    copied_zones = 0
    do i = 1, max_copy
       if (.not.values_close(new_profile%mass(i), base_composition%mass(i), &
            1.0d-6, 1.0d20)) exit
       if (.not.values_close(new_profile%radius(i), base_composition%radius(i), &
            1.0d-6, 1.0d2)) exit
       copied_zones = i
    enddo
    if (copied_zones.le.0) then
       write(*,*) '******* No aligned base composition zones found for CSM extension'
       stop 1
    endif

    extension_zones = new_profile%n_zones - copied_zones
    new_composition%n_zones = new_profile%n_zones
    new_composition%n_isotopes = base_composition%n_isotopes
    allocate(new_composition%mass_numbers(base_composition%n_isotopes))
    allocate(new_composition%charge_numbers(base_composition%n_isotopes))
    allocate(new_composition%species(base_composition%n_isotopes))
    allocate(new_composition%mass(new_profile%n_zones))
    allocate(new_composition%radius(new_profile%n_zones))
    allocate(new_composition%fractions(new_profile%n_zones, &
         base_composition%n_isotopes))
    new_composition%mass_numbers = base_composition%mass_numbers
    new_composition%charge_numbers = base_composition%charge_numbers
    new_composition%species = base_composition%species
    new_composition%mass = new_profile%mass
    new_composition%radius = new_profile%radius
    new_composition%fractions(1:copied_zones, :) = &
         base_composition%fractions(1:copied_zones, :)
    do i = 1, extension_zones
       new_composition%fractions(copied_zones+i, :) = &
            max(mix_aligned(:), floor_abundance)
       do j = 1, base_composition%n_isotopes
          if (trim(base_composition%species(j)).eq.'ni56') then
             new_composition%fractions(copied_zones+i, j) = 0.0d0
          endif
       enddo
    enddo
  end subroutine extend_composition

  subroutine get_preset_mix(composition, preset_norm, mix_aligned)
    type(csm_composition_data), intent(in) :: composition
    character(len=*), intent(in) :: preset_norm
    real*8, intent(out) :: mix_aligned(composition%n_isotopes)
    integer :: i, zone_count
    real*8 :: total

    mix_aligned(:) = 0.0d0
    select case (trim(preset_norm))
    case ('solar')
       call assign_species_fraction(composition, mix_aligned, 'h1', 0.73d0)
       call assign_species_fraction(composition, mix_aligned, 'he4', 0.25d0)
       call assign_species_fraction(composition, mix_aligned, 'c12', 0.004d0)
       call assign_species_fraction(composition, mix_aligned, 'o16', 0.010d0)
       call assign_species_fraction(composition, mix_aligned, 'ne20', 0.004d0)
       call assign_species_fraction(composition, mix_aligned, 'fe52', 0.002d0)
    case ('wd_surface_mean')
       zone_count = min(max(csm_surface_mix_zones, 1), composition%n_zones)
       do i = 1, composition%n_isotopes
          mix_aligned(i) = sum(composition%fractions( &
               composition%n_zones-zone_count+1:composition%n_zones, i)) / &
               dble(zone_count)
       enddo
    case default
       write(*,*) '******* Unsupported csm_composition_preset: ', trim(preset_norm)
       write(*,*) '******* Supported presets: solar, wd_surface_mean'
       stop 1
    end select

    do i = 1, composition%n_isotopes
       if (trim(composition%species(i)).eq.'ni56') mix_aligned(i) = 0.0d0
    enddo

    total = sum(mix_aligned)
    if (total.le.0.0d0) then
       write(*,*) '******* Selected CSM composition preset does not map to available isotopes'
       stop 1
    endif
    mix_aligned(:) = mix_aligned(:) / total
  end subroutine get_preset_mix

  subroutine assign_species_fraction(composition, mix_aligned, target_species, value)
    type(csm_composition_data), intent(in) :: composition
    real*8, intent(inout) :: mix_aligned(composition%n_isotopes)
    character(len=*), intent(in) :: target_species
    real*8, intent(in) :: value
    integer :: i

    do i = 1, composition%n_isotopes
       if (trim(composition%species(i)).eq.trim(target_species)) then
          mix_aligned(i) = value
          return
       endif
    enddo
  end subroutine assign_species_fraction

  ! derive composition scalars from normalized isotope mass fractions
  subroutine compute_mix_scalars(mix_aligned, mass_numbers, charge_numbers, &
       ye_value, abar_value)
    real*8, intent(in) :: mix_aligned(:), mass_numbers(:), charge_numbers(:)
    real*8, intent(out) :: ye_value, abar_value
    real*8 :: inv_abar
    integer :: i

    ye_value = 0.0d0
    inv_abar = 0.0d0
    do i = 1, size(mix_aligned)
       ye_value = ye_value + mix_aligned(i) * charge_numbers(i) / &
            max(mass_numbers(i), 1.0d0)
       inv_abar = inv_abar + mix_aligned(i) / max(mass_numbers(i), 1.0d0)
    enddo
    if (ye_value.le.0.0d0) ye_value = 0.5d0
    if (inv_abar.gt.0.0d0) then
       abar_value = 1.0d0 / inv_abar
    else
       abar_value = 1.0d0
    endif
  end subroutine compute_mix_scalars

  function electron_fraction_from_mix(mix_aligned, mass_numbers, charge_numbers) &
       result(ye_value)
    real*8, intent(in) :: mix_aligned(:), mass_numbers(:), charge_numbers(:)
    real*8 :: ye_value, abar_value
    call compute_mix_scalars(mix_aligned, mass_numbers, charge_numbers, &
         ye_value, abar_value)
  end function electron_fraction_from_mix

  function summarize_mix(species, mix_aligned) result(summary)
    character(len=16), intent(in) :: species(:)
    real*8, intent(in) :: mix_aligned(:)
    character(len=256) :: summary
    character(len=48) :: token
    integer :: i

    summary = ''
    do i = 1, size(species)
       if (mix_aligned(i).le.floor_abundance) cycle
       write(token,'(A,A,F7.4)') trim(species(i)), '=', mix_aligned(i)
       if (len_trim(summary).eq.0) then
          summary = trim(token)
       else
          summary = trim(summary)//', '//trim(token)
       endif
    enddo
  end function summarize_mix

  subroutine write_boundary_diagnostics(path, diag, mode_norm, preset_norm, &
       mix_summary)
    character(len=*), intent(in) :: path, mode_norm, preset_norm, mix_summary
    type(csm_boundary_diagnostics), intent(in) :: diag
    integer :: u

    u = 905
    open(unit=u, file=trim(path), status='replace', form='formatted', &
         action='write')
    write(u,'(A,A)') 'csm = ', trim(mode_norm)
    write(u,'(A,A)') 'csm_composition_preset = ', trim(preset_norm)
    write(u,'(A,A)') 'csm_mix = ', trim(mix_summary)
    write(u,'(A,ES14.6)') 'mass_loss_rate_msun_per_yr = ', diag%mass_loss_rate
    write(u,'(A,I8)') 'extension_zones = ', diag%extension_zones
    write(u,'(A,ES14.6)') 'attachment_radius_cm = ', diag%attachment_radius_cm
    write(u,'(A,ES14.6)') 'outer_radius_cm = ', diag%outer_radius_cm
    write(u,'(A,ES14.6)') 'csm_mass_msun = ', diag%csm_mass_msun
    write(u,'(A,ES14.6)') 'surface_density = ', diag%surface_density
    write(u,'(A,ES14.6)') 'first_csm_density = ', diag%first_csm_density
    write(u,'(A,ES14.6)') 'analytic_first_csm_density = ', &
         diag%analytic_first_csm_density
    write(u,'(A,ES14.6)') 'density_jump_ratio = ', diag%density_jump_ratio
    write(u,'(A,ES14.6)') 'analytic_density_jump_ratio = ', &
         diag%analytic_density_jump_ratio
    write(u,'(A,I8)') 'transition_anchor_index = ', diag%transition_anchor_index
    write(u,'(A,ES14.6)') 'transition_anchor_radius_cm = ', &
         diag%transition_anchor_radius_cm
    write(u,'(A,ES14.6)') 'transition_anchor_density = ', &
         diag%transition_anchor_density
    write(u,'(A,ES14.6)') 'bridged_surface_density = ', diag%bridged_surface_density
    write(u,'(A,ES14.6)') 'bridge_density_jump_ratio = ', &
         diag%bridge_density_jump_ratio
    write(u,'(A,A)') 'density_smoothing_mode = ', &
         trim(diag%density_smoothing_mode)
    write(u,'(A,I8)') 'density_transition_zones = ', &
         diag%density_transition_zones
    write(u,'(A,ES14.6)') 'density_inner_log_slope = ', &
         diag%density_inner_log_slope
    write(u,'(A,ES14.6)') 'density_outer_log_slope = ', &
         diag%density_outer_log_slope
    write(u,'(A,ES14.6)') 'surface_temperature = ', diag%surface_temperature
    write(u,'(A,ES14.6)') 'first_csm_temperature = ', diag%first_csm_temperature
    write(u,'(A,ES14.6)') 'temperature_jump_ratio = ', diag%temperature_jump_ratio
    write(u,'(A,ES14.6)') 'surface_velocity_cgs = ', diag%surface_velocity_cgs
    write(u,'(A,ES14.6)') 'first_csm_velocity_cgs = ', diag%first_csm_velocity_cgs
    write(u,'(A,ES14.6)') 'outer_csm_velocity_cgs = ', diag%outer_csm_velocity_cgs
    write(u,'(A,ES14.6)') 'requested_wind_velocity_kms = ', &
         diag%requested_wind_velocity_kms
    write(u,'(A,I8)') 'velocity_transition_zones = ', diag%velocity_transition_zones
    write(u,'(A,A)') 'velocity_smoothing_mode = ', &
         trim(diag%velocity_smoothing_mode)
    write(u,'(A,ES14.6)') 'wind_ye = ', diag%wind_ye
    write(u,'(A,ES14.6)') 'invariant_min = ', diag%invariant_min
    write(u,'(A,ES14.6)') 'invariant_max = ', diag%invariant_max
    write(u,'(A,A)') 'mass_scheme = ', trim(diag%mass_scheme)
    write(u,'(A,ES14.6)') 'mass_loss_duration_yr = ', diag%mass_loss_duration_yr
    write(u,'(A,ES14.6)') 'mass_loss_rate_start_msun_per_yr = ', diag%mass_loss_rate_start
    write(u,'(A,ES14.6)') 'mass_loss_rate_burst_msun_per_yr = ', diag%mass_loss_rate_burst
    write(u,'(A,ES14.6)') 'outer_density_normalization = ', diag%outer_density_normalization
    write(u,'(A,L1)') 'terminal_surface_repaired = ', diag%terminal_surface_repaired
    write(u,'(A,I8)') 'terminal_surface_anchor_index = ', diag%terminal_surface_anchor_index
    write(u,'(A,ES14.6)') 'terminal_density_original = ', diag%terminal_density_original
    write(u,'(A,ES14.6)') 'terminal_temperature_original = ', diag%terminal_temperature_original
    close(u)
  end subroutine write_boundary_diagnostics

  subroutine ensure_directory(path)
    character(len=*), intent(in) :: path
    character(len=2048) :: cmd
    cmd='mkdir -p -- '//shell_quote(trim(path))
    call system(trim(cmd))
  end subroutine ensure_directory

  ! integrate each finite spherical shell using its mean density
  subroutine compute_shell_masses(inner_radius, outer_radii, densities, shell_masses)
    real*8, intent(in) :: inner_radius
    real*8, intent(in) :: outer_radii(:), densities(:)
    real*8, intent(out) :: shell_masses(size(outer_radii))
    real*8 :: r_inner
    integer :: i

    r_inner = inner_radius
    do i = 1, size(outer_radii)
       shell_masses(i) = densities(i) * (4.0d0/3.0d0) * pi_cgs * &
            (outer_radii(i)**3 - r_inner**3)
       r_inner = outer_radii(i)
    enddo
  end subroutine compute_shell_masses

  subroutine cumulative_from_surface(surface_mass, shell_masses, new_mass)
    real*8, intent(in) :: surface_mass
    real*8, intent(in) :: shell_masses(:)
    real*8, intent(out) :: new_mass(size(shell_masses))
    integer :: i

    new_mass(1) = surface_mass + shell_masses(1)
    do i = 2, size(shell_masses)
       new_mass(i) = new_mass(i-1) + shell_masses(i)
    enddo
  end subroutine cumulative_from_surface

  subroutine geomspace_positive(start_val, end_val, n_terms, values)
    real*8, intent(in) :: start_val, end_val
    integer, intent(in) :: n_terms
    real*8, intent(out) :: values(n_terms)
    integer :: i
    real*8 :: log_start, log_end, frac

    if (n_terms.le.0) return
    if (n_terms.eq.1) then
       values(1) = end_val
       return
    endif
    log_start = log(max(start_val, floor_abundance))
    log_end = log(max(end_val, floor_abundance))
    do i = 1, n_terms
       frac = dble(i-1) / dble(n_terms-1)
       values(i) = exp(log_start + frac * (log_end - log_start))
    enddo
    values(n_terms) = end_val
  end subroutine geomspace_positive

  function lower_trimmed(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len_trim(adjustl(text))) :: out
    integer :: i, code

    out = trim(adjustl(text))
    do i = 1, len(out)
       code = iachar(out(i:i))
       if (code.ge.iachar('A') .and. code.le.iachar('Z')) out(i:i) = &
            achar(code + 32)
    enddo
  end function lower_trimmed

  logical function values_close(a, b, rtol, atol)
    real*8, intent(in) :: a, b, rtol, atol
    values_close = abs(a-b).le.max(atol, rtol*max(abs(a), abs(b)))
  end function values_close

  function safe_ratio(num, den) result(value)
    real*8, intent(in) :: num, den
    real*8 :: value
    if (abs(den).le.floor_abundance) then
       value = 0.0d0
    else
       value = num / den
    endif
  end function safe_ratio

  function species_name(mass_num, charge_num) result(name)
    integer, intent(in) :: mass_num, charge_num
    character(len=16) :: name

    name = 'unknown'
    select case (mass_num)
    case (1)
       if (charge_num.eq.0) name = 'n'
       if (charge_num.eq.1) name = 'h1'
    case (4)
       if (charge_num.eq.2) name = 'he4'
    case (12)
       if (charge_num.eq.6) name = 'c12'
    case (16)
       if (charge_num.eq.8) name = 'o16'
    case (20)
       if (charge_num.eq.10) name = 'ne20'
    case (24)
       if (charge_num.eq.12) name = 'mg24'
    case (28)
       if (charge_num.eq.14) name = 'si28'
       if (charge_num.eq.28) name = 'ni56'
    case (32)
       if (charge_num.eq.16) name = 's32'
    case (36)
       if (charge_num.eq.18) name = 'ar36'
    case (40)
       if (charge_num.eq.20) name = 'ca40'
    case (44)
       if (charge_num.eq.22) name = 'ti44'
    case (48)
       if (charge_num.eq.24) name = 'cr48'
    case (52)
       if (charge_num.eq.26) name = 'fe52'
    case (56)
       if (charge_num.eq.28) name = 'ni56'
    end select
  end function species_name

end module csm_shared

module csm_profile_builders

  use csm_shared
  use parameters, only: mass_scheme, csm_mass_msun, csm_outer_radius_cm, &
       csm_wind_velocity_kms, csm_temperature_k, csm_mdot_msun_per_yr, &
       csm_extension_zones, csm_blend_zones
  implicit none

  real*8, parameter :: wind_history_rise_end_yr=5.0d4
  real*8, parameter :: wind_history_plateau_end_yr=4.0d5
  real*8, parameter :: wind_history_fall_end_yr=8.0d5
  real*8, parameter :: wind_history_tsn_yr=6.0d5
  real*8, parameter :: wind_history_total_mass_msun=5.0d-1
  real*8, parameter :: wind_history_area_to_tsn_yr=5.25d5
  real*8, parameter :: wind_history_mdot_peak=wind_history_total_mass_msun/wind_history_area_to_tsn_yr

contains

  subroutine build_single_attached(base_profile, base_composition, mix_aligned, &
       new_profile, diag)
    type(csm_profile_data), intent(in) :: base_profile
    type(csm_composition_data), intent(in) :: base_composition
    real*8, intent(in) :: mix_aligned(:)
    type(csm_profile_data), intent(out) :: new_profile
    type(csm_boundary_diagnostics), intent(out) :: diag

    type(csm_profile_data) :: base_work
    integer :: n_old, n_new, density_transition_zones, velocity_transition_zones
    integer :: surface_anchor_index
    logical :: terminal_surface_repaired
    real*8 :: surface_mass, surface_radius, surface_density
    real*8 :: surface_temp, surface_vel, surface_omega, surface_step, wind_ye
    real*8 :: inner_slope, outer_slope, target_mass_msun, target_mass_cgs
    real*8 :: delta_t_yr, mdot_start, mdot_burst, outer_norm
    real*8 :: terminal_density_original, terminal_temperature_original
    character(len=32) :: scheme_norm
    real*8, allocatable :: new_radius(:), rho_new(:), rho_analytic(:)
    real*8, allocatable :: temp_new(:), vel_new(:), new_mass(:), invariant(:)
    real*8, allocatable :: shell_masses(:)

    ! repair an exported boundary sentinel before attaching material
    base_work = base_profile
    call repair_terminal_surface_state(base_work, terminal_surface_repaired, &
         surface_anchor_index, terminal_density_original, terminal_temperature_original)

    n_old = base_work%n_zones
    n_new = csm_extension_zones
    if (n_new.lt.2) then
       write(*,*) '******* csm_extension_zones must be >= 2 for attached CSM'
       stop 1
    endif

    surface_mass = base_work%mass(n_old)
    surface_radius = base_work%radius(n_old)
    surface_density = base_work%rho(n_old)
    surface_temp = base_work%temp(n_old)
    surface_vel = base_work%vel(n_old)
    surface_omega = base_work%omega(n_old)
    if (csm_outer_radius_cm.le.surface_radius) then
       write(*,*) '******* csm_outer_radius_cm must exceed the progenitor surface radius'
       stop 1
    endif

    surface_step = max(surface_radius-base_work%radius(max(1,n_old-1)), &
         surface_radius*1.0d-6)
    wind_ye = electron_fraction_from_mix(mix_aligned, base_composition%mass_numbers, &
         base_composition%charge_numbers)
    scheme_norm = lower_trimmed(mass_scheme)
    if (trim(scheme_norm).eq.'continous_mass') scheme_norm = 'continuous_mass'

    ! direct-mass and wind schemes use different mass normalizations
    call derive_target_csm_mass(surface_radius, scheme_norm, target_mass_msun, &
         delta_t_yr, mdot_start, mdot_burst)
    target_mass_cgs = target_mass_msun*msun_cgs

    allocate(new_radius(n_new), rho_new(n_new), rho_analytic(n_new))
    allocate(temp_new(n_new), vel_new(n_new), new_mass(n_new), invariant(n_new))
    allocate(shell_masses(n_new))

    call build_wind_radii(surface_radius, csm_outer_radius_cm, n_new, &
         surface_step, new_radius)
    ! build the density profile before assigning shell masses
    call solve_attached_density(base_work, new_radius, target_mass_cgs, &
         rho_new, rho_analytic, outer_norm, density_transition_zones, &
         inner_slope, outer_slope)

    call compute_shell_masses(surface_radius, new_radius, rho_new, shell_masses)
    call cumulative_from_surface(surface_mass, shell_masses, new_mass)
    call build_temperature_bridge(surface_radius, surface_temp, new_radius, temp_new)
    call build_velocity_bridge(surface_vel, scheme_norm, new_radius, vel_new, &
         velocity_transition_zones)

    invariant(:) = 4.0d0*pi_cgs*new_radius(:)**3*rho_new(:)

    call copy_profile_with_append(base_work, base_work%mass, base_work%temp, &
         base_work%rho, new_mass, new_radius, temp_new, rho_new, vel_new, &
         wind_ye, surface_omega, new_profile)

    if (trim(scheme_norm).eq.'non_constant_mass_loss' .or. trim(scheme_norm).eq.'constant_mass_loss') then
       diag%mass_loss_rate = mdot_burst
    else
       diag%mass_loss_rate = 0.0d0
    endif
    diag%extension_zones = n_new
    diag%attachment_radius_cm = surface_radius
    diag%outer_radius_cm = csm_outer_radius_cm
    diag%csm_mass_msun = (new_mass(n_new)-surface_mass)/msun_cgs
    diag%surface_mass_g = surface_mass
    diag%surface_radius_cm = surface_radius
    diag%surface_density = surface_density
    diag%first_csm_density = rho_new(1)
    diag%analytic_first_csm_density = rho_analytic(1)
    diag%density_jump_ratio = safe_ratio(rho_new(1),surface_density)
    diag%analytic_density_jump_ratio = safe_ratio(rho_analytic(1),surface_density)
    diag%transition_anchor_index = n_old
    diag%transition_anchor_radius_cm = surface_radius
    diag%transition_anchor_density = surface_density
    diag%bridged_surface_density = surface_density
    diag%bridge_density_jump_ratio = safe_ratio(rho_new(1),surface_density)
    if (trim(scheme_norm).eq.'constant_density') then
       diag%density_smoothing_mode = 'constant_density_mass_exact'
    else if (trim(scheme_norm).eq.'non_constant_mass_loss' .or. trim(scheme_norm).eq.'constant_mass_loss') then
       diag%density_smoothing_mode = 'physical_discontinuity_r2_wind'
    else
       diag%density_smoothing_mode = 'attached_c1_mass_exact'
    endif
    diag%density_transition_zones = density_transition_zones
    diag%density_inner_log_slope = inner_slope
    diag%density_outer_log_slope = outer_slope
    diag%surface_temperature = surface_temp
    diag%first_csm_temperature = temp_new(1)
    diag%temperature_jump_ratio = safe_ratio(temp_new(1),surface_temp)
    diag%surface_velocity_cgs = surface_vel
    diag%first_csm_velocity_cgs = vel_new(1)
    diag%outer_csm_velocity_cgs = vel_new(n_new)
    diag%requested_wind_velocity_kms = csm_wind_velocity_kms
    diag%velocity_transition_zones = velocity_transition_zones
    if (trim(scheme_norm).eq.'non_constant_mass_loss' .or. trim(scheme_norm).eq.'constant_mass_loss') then
       diag%velocity_smoothing_mode = 'physical_ejecta_wind_jump'
    else if (trim(scheme_norm).eq.'constant_density') then
       diag%velocity_smoothing_mode = 'log_radius_c2_smootherstep'
    else
       diag%velocity_smoothing_mode = 'surface_comoving_unchanged'
    endif
    diag%wind_ye = wind_ye
    diag%invariant_min = minval(invariant)
    diag%invariant_max = maxval(invariant)
    diag%mass_scheme = scheme_norm
    diag%mass_loss_duration_yr = delta_t_yr
    diag%mass_loss_rate_start = mdot_start
    diag%mass_loss_rate_burst = mdot_burst
    diag%outer_density_normalization = outer_norm
    diag%terminal_surface_repaired = terminal_surface_repaired
    diag%terminal_surface_anchor_index = surface_anchor_index
    diag%terminal_density_original = terminal_density_original
    diag%terminal_temperature_original = terminal_temperature_original

    deallocate(new_radius,rho_new,rho_analytic,temp_new,vel_new,new_mass,invariant, &
         shell_masses)
  end subroutine build_single_attached

  ! integrate the requested wind history or use the direct shell mass
  subroutine derive_target_csm_mass(surface_radius, scheme_norm, target_mass_msun, &
       delta_t_yr, mdot_start, mdot_burst)
    real*8, intent(in) :: surface_radius
    character(len=*), intent(in) :: scheme_norm
    real*8, intent(out) :: target_mass_msun, delta_t_yr, mdot_start, mdot_burst
    real*8 :: wind_cgs, launch_start

    target_mass_msun = 0.0d0
    delta_t_yr = 0.0d0
    mdot_start = 0.0d0
    mdot_burst = 0.0d0

    if (trim(scheme_norm).eq.'continuous_mass' .or. trim(scheme_norm).eq.'constant_density') then
       if (csm_mass_msun.le.0.0d0) then
          write(*,*) '******* csm_mass_msun must be > 0 for direct-mass schemes'
          stop 1
       endif
       target_mass_msun = csm_mass_msun
       return
    endif

    if (trim(scheme_norm).ne.'constant_mass_loss' .and. trim(scheme_norm).ne.'non_constant_mass_loss') then
       write(*,*) '******* Unsupported mass_scheme: ',trim(scheme_norm)
       write(*,*) '******* Supported schemes: continuous_mass, constant_mass_loss, non_constant_mass_loss, constant_density'
       stop 1
    endif

    wind_cgs=csm_wind_velocity_kms*1.0d5
    delta_t_yr=(csm_outer_radius_cm-surface_radius)/wind_cgs/year_cgs

    if (trim(scheme_norm).eq.'constant_mass_loss') then
       ! Moriya et al. (2023): steady spherical wind,
       ! rho = Mdot/(4*pi*r^2*vwind).  The finite-domain mass follows directly
       ! from the wind crossing time, Delta M = Mdot*Delta r/vwind
       mdot_start=csm_mdot_msun_per_yr
       mdot_burst=csm_mdot_msun_per_yr
       target_mass_msun=csm_mdot_msun_per_yr*delta_t_yr
    else
       if (delta_t_yr.ge.wind_history_tsn_yr) then
          write(*,*) '******* Badenes CSM radius reaches before the onset of its wind history.'
          write(*,*) '        Reduce csm_outer_radius_cm or increase csm_wind_velocity_kms.'
          stop 1
       endif
       launch_start=wind_history_tsn_yr-delta_t_yr
       mdot_start=wind_history_mdot(launch_start)
       mdot_burst=wind_history_mdot(wind_history_tsn_yr)
       target_mass_msun=wind_history_cumulative_mass(wind_history_tsn_yr)-wind_history_cumulative_mass(launch_start)
    endif

    if (target_mass_msun.le.0.0d0 .or. target_mass_msun.ne.target_mass_msun) then
       write(*,*) '******* Wind integration produced non-positive or invalid CSM mass'
       stop 1
    endif
  end subroutine derive_target_csm_mass


  ! evaluate the prescribed rise, plateau and decline in mass-loss rate
  real*8 function wind_history_mdot(time_yr) result(mdot)
    real*8, intent(in) :: time_yr
    real*8 :: x

    if (time_yr.le.0.0d0 .or. time_yr.ge.wind_history_fall_end_yr) then
       mdot=0.0d0
    else if (time_yr.lt.wind_history_rise_end_yr) then
       mdot=wind_history_mdot_peak*time_yr/wind_history_rise_end_yr
    else if (time_yr.le.wind_history_plateau_end_yr) then
       mdot=wind_history_mdot_peak
    else
       x=(time_yr-wind_history_plateau_end_yr)/(wind_history_fall_end_yr-wind_history_plateau_end_yr)
       mdot=wind_history_mdot_peak*(1.0d0-x)
    endif
  end function wind_history_mdot


  ! integrate the piecewise wind history analytically
  real*8 function wind_history_cumulative_mass(time_yr) result(mass_msun)
    real*8, intent(in) :: time_yr
    real*8 :: t, x, mass_rise, mass_plateau

    t=max(0.0d0,min(time_yr,wind_history_fall_end_yr))
    mass_rise=0.5d0*wind_history_mdot_peak*wind_history_rise_end_yr
    mass_plateau=wind_history_mdot_peak*(wind_history_plateau_end_yr-wind_history_rise_end_yr)

    if (t.le.wind_history_rise_end_yr) then
       mass_msun=0.5d0*wind_history_mdot_peak*t*t/wind_history_rise_end_yr
    else if (t.le.wind_history_plateau_end_yr) then
       mass_msun=mass_rise+wind_history_mdot_peak*(t-wind_history_rise_end_yr)
    else
       x=t-wind_history_plateau_end_yr
       mass_msun=mass_rise+mass_plateau+wind_history_mdot_peak*(x- &
            0.5d0*x*x/(wind_history_fall_end_yr-wind_history_plateau_end_yr))
    endif
  end function wind_history_cumulative_mass


  ! replace an exported passive sentinel with an admissible surface state
  subroutine repair_terminal_surface_state(profile, repaired, anchor_index, &
       rho_original, temp_original)
    type(csm_profile_data), intent(inout) :: profile
    logical, intent(out) :: repaired
    integer, intent(out) :: anchor_index
    real*8, intent(out) :: rho_original, temp_original

    integer :: n, j, i0
    real*8 :: slope, dx, rho_ratio, temp_anchor
    logical :: sentinel

    n = profile%n_zones
    repaired = .false.
    anchor_index = max(1,n)
    rho_original = 0.0d0
    temp_original = 0.0d0
    if (n.lt.2) return

    rho_original = profile%rho(n)
    temp_original = profile%temp(n)
    rho_ratio = profile%rho(n)/max(profile%rho(n-1),1.0d-99)

    sentinel = (profile%temp(n).le.0.0d0 .and. profile%temp(n-1).gt.0.0d0) .or. &
         (rho_ratio.lt.1.0d-10)
    if (.not.sentinel) return

    anchor_index = n-1
    do while (anchor_index.gt.1)
       if (profile%temp(anchor_index).gt.0.0d0 .and. &
            profile%rho(anchor_index).gt.1.0d-90) exit
       anchor_index = anchor_index-1
    enddo
    if (anchor_index.lt.1) then
       write(*,*) '******* Could not identify a physical ejecta surface below terminal sentinel.'
       stop 1
    endif

    i0 = max(1,anchor_index-1)
    if (anchor_index.gt.i0) then
       dx = log(max(profile%radius(anchor_index),1.0d-99)) - &
            log(max(profile%radius(i0),1.0d-99))
       if (abs(dx).gt.1.0d-12) then
          slope = (log(max(profile%rho(anchor_index),1.0d-99)) - &
               log(max(profile%rho(i0),1.0d-99)))/dx
       else
          slope = -3.0d0
       endif
    else
       slope = -3.0d0
    endif
    slope = min(-5.0d-2,max(-1.2d1,slope))
    temp_anchor = max(profile%temp(anchor_index),csm_temperature_k)

    do j=anchor_index+1,n
       profile%rho(j) = max(1.0d-99,profile%rho(anchor_index) * &
            (profile%radius(j)/profile%radius(anchor_index))**slope)
       profile%temp(j) = temp_anchor
    enddo

    repaired = .true.
    write(*,*) ' Reconstructed terminal ejecta boundary state for attached CSM.'
    write(*,*) '  anchor_index, radius_cm = ',anchor_index,profile%radius(anchor_index)
    write(*,*) '  original terminal rho,T = ',rho_original,temp_original
    write(*,*) '  repaired terminal rho,T = ',profile%rho(n),profile%temp(n)
  end subroutine repair_terminal_surface_state

  ! normalize the finite-volume envelope and bridge to the requested mass
  subroutine solve_attached_density(base_profile, radii, target_mass_cgs, &
       rho_final, rho_analytic, outer_norm, transition_zones, inner_slope, outer_slope)
    type(csm_profile_data), intent(in) :: base_profile
    real*8, intent(in) :: radii(:), target_mass_cgs
    real*8, intent(out) :: rho_final(size(radii)), rho_analytic(size(radii))
    real*8, intent(out) :: outer_norm, inner_slope, outer_slope
    integer, intent(out) :: transition_zones

    integer :: n_blend, iter, i
    real*8 :: surface_radius, surface_rho, log_lo, log_hi, log_mid
    real*8 :: mass_lo, mass_hi, mass_mid, relerr, bridge_dex, shell_inner_radius
    real*8, parameter :: max_bridge_dex_per_zone = 2.5d-1
    real*8, allocatable :: trial(:), shell_masses(:)

    surface_radius = base_profile%radius(base_profile%n_zones)
    surface_rho = base_profile%rho(base_profile%n_zones)
    inner_slope = clipped_density_slope(base_profile%radius,base_profile%rho, &
         base_profile%n_zones)

    ! a uniform shell keeps its physical density jump at attachment
    if (trim(lower_trimmed(mass_scheme)).eq.'constant_density') then
       outer_norm = target_mass_cgs / max((4.0d0*pi_cgs/3.0d0) * &
            (radii(size(radii))**3-surface_radius**3),1.0d-99)
       rho_final(:)=max(outer_norm,1.0d-99)
       rho_analytic(:)=rho_final(:)
       transition_zones=0
       inner_slope=0.0d0
       outer_slope=0.0d0
       return
    endif

    if (trim(lower_trimmed(mass_scheme)).eq.'constant_mass_loss' .or. &
         trim(lower_trimmed(mass_scheme)).eq.'non_constant_mass_loss') then
       call fill_wind_density(surface_radius,radii,rho_final,rho_analytic,outer_norm)
       transition_zones=0
       inner_slope=-2.0d0
       outer_slope=-2.0d0
       return
    endif

    ! only the attached continuous-mass envelope uses the density bridge
    n_blend = min(max(csm_blend_zones,0),max(0,size(radii)-1))
    shell_inner_radius = surface_radius
    if (n_blend.gt.0) shell_inner_radius = radii(n_blend)
    outer_slope = attached_background_log_slope()

    allocate(trial(size(radii)),shell_masses(size(radii)))

    do
       log_lo = log(max(surface_rho,1.0d-99))-40.0d0
       log_hi = log(max(surface_rho,1.0d-99))+4.0d0

       call fill_attached_candidate(base_profile,radii,exp(log_lo), &
            n_blend,trial,rho_analytic)
       call compute_shell_masses(surface_radius,radii,trial,shell_masses)
       mass_lo = sum(shell_masses)

       call fill_attached_candidate(base_profile,radii,exp(log_hi), &
            n_blend,trial,rho_analytic)
       call compute_shell_masses(surface_radius,radii,trial,shell_masses)
       mass_hi = sum(shell_masses)

       do while (target_mass_cgs.gt.mass_hi .and. log_hi.lt.650.0d0)
          log_hi = log_hi + 4.0d0
          call fill_attached_candidate(base_profile,radii,exp(log_hi), &
               n_blend,trial,rho_analytic)
          call compute_shell_masses(surface_radius,radii,trial,shell_masses)
          mass_hi = sum(shell_masses)
       enddo

       if (target_mass_cgs.ge.mass_lo .and. target_mass_cgs.le.mass_hi) exit
       if (target_mass_cgs.lt.mass_lo .and. n_blend.gt.0) then
          n_blend = n_blend-1
          cycle
       endif
       write(*,*) '******* Could not bracket the requested attached-CSM mass.'
       write(*,*) '        target_mass_g, minimum_bridge_mass_g, bracket_mass_g = ', &
            target_mass_cgs,mass_lo,mass_hi
       write(*,*) '        Increase csm_extension_zones or reduce csm_blend_zones.'
       stop 1
    enddo

    ! bisect the density normalization until the integrated mass matches
    do iter=1,180
       log_mid = 0.5d0*(log_lo+log_hi)
       call fill_attached_candidate(base_profile,radii,exp(log_mid), &
            n_blend,trial,rho_analytic)
       call compute_shell_masses(surface_radius,radii,trial,shell_masses)
       mass_mid = sum(shell_masses)
       if (mass_mid.lt.target_mass_cgs) then
          log_lo = log_mid
       else
          log_hi = log_mid
       endif
       relerr = abs(mass_mid-target_mass_cgs)/max(target_mass_cgs,1.0d-99)
       if (relerr.lt.1.0d-12) exit
    enddo

    outer_norm = exp(0.5d0*(log_lo+log_hi))
    call fill_attached_candidate(base_profile,radii,outer_norm,n_blend, &
         rho_final,rho_analytic)
    transition_zones = n_blend
    shell_inner_radius = surface_radius
    if (n_blend.gt.0) shell_inner_radius = radii(n_blend)
    outer_slope = attached_background_log_slope()

    do i=1,size(rho_final)
       if (rho_final(i).ne.rho_final(i) .or. rho_final(i).le.0.0d0) then
          write(*,*) '******* Attached CSM density bridge produced an invalid state at zone ',i
          stop 1
       endif
    enddo

    call compute_shell_masses(surface_radius,radii,rho_final,shell_masses)
    relerr = abs(sum(shell_masses)-target_mass_cgs)/max(target_mass_cgs,1.0d-99)
    if (relerr.gt.1.0d-9) then
       write(*,*) '******* Attached CSM mass conservation failed after bridge solve.'
       write(*,*) '        fractional mass error = ',relerr
       stop 1
    endif

    bridge_dex = attached_bridge_max_dex(surface_rho,rho_final,n_blend)
    if (bridge_dex.gt.max_bridge_dex_per_zone) then
       write(*,*) '******* WARNING: attached CSM bridge is under-resolved at the requested width.'
       write(*,*) '        max |Delta log10 rho| per bridge zone = ',bridge_dex
       write(*,*) '        Increase csm_blend_zones or csm_extension_zones if needed.'
    endif

    deallocate(trial,shell_masses)
  end subroutine solve_attached_density

  real*8 function attached_bridge_max_dex(surface_rho,rho_values,n_blend) result(max_dex)
    real*8, intent(in) :: surface_rho,rho_values(:)
    integer, intent(in) :: n_blend
    integer :: i,ncheck
    real*8 :: prev_rho,dex

    max_dex = 0.0d0
    if (n_blend.le.0) return
    prev_rho = max(surface_rho,1.0d-99)
    ncheck = min(n_blend,size(rho_values))
    do i=1,ncheck
       dex = abs(log10(max(rho_values(i),1.0d-99)/prev_rho))
       max_dex = max(max_dex,dex)
       prev_rho = max(rho_values(i),1.0d-99)
    enddo
  end function attached_bridge_max_dex

  ! construct a trial log-density bridge for the mass solve
  subroutine fill_attached_candidate(base_profile,radii,outer_norm, &
       n_blend,rho_candidate,rho_analytic)
    type(csm_profile_data), intent(in) :: base_profile
    real*8, intent(in) :: radii(:), outer_norm
    integer, intent(in) :: n_blend
    real*8, intent(out) :: rho_candidate(size(radii)),rho_analytic(size(radii))

    integer :: i
    real*8 :: surface_radius,surface_rho,inner_slope,x0,x1,x,u,dx,shell_inner_radius
    real*8 :: y0,y1,m0,m1,h00,h10,h01,h11,shape,secant

    surface_radius = base_profile%radius(base_profile%n_zones)
    surface_rho = base_profile%rho(base_profile%n_zones)
    inner_slope = clipped_density_slope(base_profile%radius,base_profile%rho, &
         base_profile%n_zones)
    x0 = log(max(surface_radius,1.0d-99))
    shell_inner_radius = surface_radius
    if (n_blend.gt.0) shell_inner_radius = radii(n_blend)

    do i=1,size(radii)
       shape = attached_background_shape(radii(i),surface_radius)
       rho_analytic(i) = max(1.0d-99,outer_norm*shape)
       rho_candidate(i) = rho_analytic(i)
    enddo

    if (n_blend.le.0) return

    ! join the density endpoints with a cubic in log radius and log density
    x1 = log(max(radii(n_blend),surface_radius*(1.0d0+1.0d-12)))
    y0 = log(max(surface_rho,1.0d-99))
    y1 = log(max(rho_analytic(n_blend),1.0d-99))
    dx = max(x1-x0,1.0d-99)
    secant = (y1-y0)/dx
    m0 = inner_slope
    m1 = attached_background_log_slope()
    call limit_monotone_hermite_slopes(secant,m0,m1)

    do i=1,n_blend
          x = log(max(radii(i),1.0d-99))
          u = max(0.0d0,min(1.0d0,(x-x0)/dx))
          h00 = 2.0d0*u**3-3.0d0*u**2+1.0d0
          h10 = u**3-2.0d0*u**2+u
          h01 = -2.0d0*u**3+3.0d0*u**2
          h11 = u**3-u**2
          rho_candidate(i) = exp(h00*y0+h10*dx*m0+h01*y1+h11*dx*m1)
          rho_candidate(i) = max(rho_candidate(i),1.0d-99)
    enddo
    rho_candidate(n_blend) = rho_analytic(n_blend)
  end subroutine fill_attached_candidate

  ! keep the log-density bridge monotone
  subroutine limit_monotone_hermite_slopes(secant,slope_left,slope_right)
    real*8, intent(in) :: secant
    real*8, intent(inout) :: slope_left,slope_right
    real*8 :: alpha,beta,scale,norm2

    if (abs(secant).le.1.0d-14) then
       slope_left = 0.0d0
       slope_right = 0.0d0
       return
    endif
    if (slope_left*secant.le.0.0d0) slope_left = 0.0d0
    if (slope_right*secant.le.0.0d0) slope_right = 0.0d0
    ! limit the endpoint slopes to avoid overshooting the log-density bridge
    alpha = slope_left/secant
    beta = slope_right/secant
    norm2 = alpha*alpha+beta*beta
    if (norm2.gt.9.0d0) then
       scale = 3.0d0/sqrt(norm2)
       slope_left = scale*alpha*secant
       slope_right = scale*beta*secant
    endif
  end subroutine limit_monotone_hermite_slopes

  ! derive each wind density from integrated mass divided by shell volume
  subroutine fill_wind_density(surface_radius,radii,rho_final,rho_analytic,outer_norm)
    real*8, intent(in) :: surface_radius,radii(:)
    real*8, intent(out) :: rho_final(size(radii)),rho_analytic(size(radii)),outer_norm
    integer :: i
    real*8 :: rin,rout,rmid,volume,wind_cgs,lookback_in,lookback_out
    real*8 :: launch_lo,launch_hi,shell_mass_msun,mdot_mid,mdot_cgs
    character(len=32) :: scheme_norm

    scheme_norm=lower_trimmed(mass_scheme)
    wind_cgs=csm_wind_velocity_kms*1.0d5
    rin=surface_radius

    ! integrate wind mass over each radial interval, then divide by its volume
    do i=1,size(radii)
       rout=radii(i)
       volume=(4.0d0*pi_cgs/3.0d0)*(rout**3-rin**3)
       if (trim(scheme_norm).eq.'constant_mass_loss') then
          shell_mass_msun=csm_mdot_msun_per_yr*(rout-rin)/wind_cgs/year_cgs
          mdot_mid=csm_mdot_msun_per_yr
       else
          lookback_in=(rin-surface_radius)/wind_cgs/year_cgs
          lookback_out=(rout-surface_radius)/wind_cgs/year_cgs
          launch_lo=wind_history_tsn_yr-lookback_out
          launch_hi=wind_history_tsn_yr-lookback_in
          shell_mass_msun=wind_history_cumulative_mass(launch_hi)-wind_history_cumulative_mass(launch_lo)
          rmid=0.5d0*(rin+rout)
          mdot_mid=wind_history_mdot(wind_history_tsn_yr-(rmid-surface_radius)/wind_cgs/year_cgs)
       endif
       rho_final(i)=max(shell_mass_msun*msun_cgs/max(volume,1.0d-99),1.0d-99)
       rmid=0.5d0*(rin+rout)
       mdot_cgs=mdot_mid*msun_cgs/year_cgs
       rho_analytic(i)=max(mdot_cgs/(4.0d0*pi_cgs*wind_cgs*rmid**2),1.0d-99)
       rin=rout
    enddo

    outer_norm=rho_analytic(size(radii))* &
         (0.5d0*(radii(max(1,size(radii)-1))+radii(size(radii)))/surface_radius)**2
  end subroutine fill_wind_density


  real*8 function attached_background_shape(radius,surface_radius) result(shape)
    real*8, intent(in) :: radius,surface_radius
    character(len=32) :: scheme_norm
    scheme_norm = lower_trimmed(mass_scheme)
    if (trim(scheme_norm).eq.'constant_density') then
       shape = 1.0d0
    else if (trim(scheme_norm).eq.'continuous_mass') then
       shape = (radius/surface_radius)**(-3.0d0)
    else
       shape = (radius/surface_radius)**(-2.0d0)
    endif
  end function attached_background_shape

  real*8 function attached_background_log_slope() result(slope)
    character(len=32) :: scheme_norm
    scheme_norm = lower_trimmed(mass_scheme)
    if (trim(scheme_norm).eq.'constant_density') then
       slope = 0.0d0
    else if (trim(scheme_norm).eq.'continuous_mass') then
       slope = -3.0d0
    else
       slope = -2.0d0
    endif
  end function attached_background_log_slope


  real*8 function smoothstep01(value) result(weight)
    real*8, intent(in) :: value
    real*8 :: x
    x = max(0.0d0,min(1.0d0,value))
    weight = x*x*(3.0d0-2.0d0*x)
  end function smoothstep01

  real*8 function smootherstep01(value) result(weight)
    real*8, intent(in) :: value
    real*8 :: x
    x = max(0.0d0,min(1.0d0,value))
    weight = x*x*x*(10.0d0+x*(-15.0d0+6.0d0*x))
  end function smootherstep01

  ! assign comoving material, free wind or a smooth uniform-shell transition
  subroutine build_velocity_bridge(surface_vel,scheme_norm,radii,vel_new, &
       transition_zones)
    real*8, intent(in) :: surface_vel,radii(:)
    character(len=*), intent(in) :: scheme_norm
    real*8, intent(out) :: vel_new(size(radii))
    integer, intent(out) :: transition_zones
    integer :: i
    real*8 :: x0,x1,x,weight,wind_vel

    if (trim(scheme_norm).ne.'non_constant_mass_loss' .and. trim(scheme_norm).ne.'constant_mass_loss' .and. &
         trim(scheme_norm).ne.'constant_density') then
       vel_new(:) = surface_vel
       transition_zones = 0
       return
    endif

    ! convert the requested wind speed to cgs units
    wind_vel = csm_wind_velocity_kms*1.0d5

    if (trim(scheme_norm).eq.'non_constant_mass_loss' .or. trim(scheme_norm).eq.'constant_mass_loss') then
       vel_new(:) = wind_vel
       transition_zones = 0
       return
    endif

    ! smooth the uniform-shell velocity transition in log radius
    transition_zones = min(size(radii),max(2,csm_blend_zones))
    x0 = log(max(radii(1),1.0d-99))
    x1 = log(max(radii(transition_zones),radii(1)*(1.0d0+1.0d-12)))

    do i=1,size(radii)
       if (i.le.transition_zones) then
          x = (log(max(radii(i),1.0d-99))-x0)/max(x1-x0,1.0d-99)
          weight = smootherstep01(x)
          vel_new(i) = surface_vel+weight*(wind_vel-surface_vel)
       else
          vel_new(i) = wind_vel
       endif
    enddo

    vel_new(1) = surface_vel
    vel_new(transition_zones:size(radii)) = wind_vel
  end subroutine build_velocity_bridge

  ! join positive endpoint temperatures in log radius and log temperature
  subroutine build_temperature_bridge(surface_radius,surface_temp,radii,temp_new)
    real*8, intent(in) :: surface_radius,surface_temp,radii(:)
    real*8, intent(out) :: temp_new(size(radii))
    integer :: i,n_blend
    real*8 :: x0,x1,x,w,log_left,log_right

    ! blend positive endpoint temperatures in log space
    n_blend = min(max(csm_blend_zones,1),max(1,size(radii)-1))
    x0 = log(max(surface_radius,1.0d-99))
    x1 = log(max(radii(n_blend),surface_radius*(1.0d0+1.0d-12)))
    do i=1,size(radii)
       if (i.le.n_blend .and. surface_temp.gt.0.0d0 .and. csm_temperature_k.gt.0.0d0) then
          x = log(max(radii(i),1.0d-99))
          w = smoothstep01((x-x0)/max(x1-x0,1.0d-99))
          log_left = log(surface_temp)
          log_right = log(csm_temperature_k)
          temp_new(i) = exp((1.0d0-w)*log_left+w*log_right)
       else
          temp_new(i) = csm_temperature_k
       endif
    enddo
  end subroutine build_temperature_bridge

  ! resolve the attachment before expanding the radial spacing
  subroutine build_wind_radii(surface_radius,outer_radius,n_terms,first_step,radii)
    real*8, intent(in) :: surface_radius,outer_radius,first_step
    integer, intent(in) :: n_terms
    real*8, intent(out) :: radii(n_terms)
    real*8 :: total_span,first_step_eff
    integer :: i

    total_span = outer_radius-surface_radius
    if (total_span.le.0.0d0) then
       write(*,*) '******* Attached-CSM radial span must be positive'
       stop 1
    endif

    first_step_eff = min(first_step,0.5d0*total_span)
    first_step_eff = max(first_step_eff,total_span*1.0d-12)
    call fill_geometric_radii(first_step_eff,radii)

    if (.not.radii_strictly_increasing(radii)) then
       first_step_eff = min(first_step,4.0d0*total_span/dble(max(n_terms,1)))
       first_step_eff = max(first_step_eff,total_span*1.0d-12)
       call fill_geometric_radii(first_step_eff,radii)
    endif

    do i=2,n_terms
       if (radii(i).le.radii(i-1)) then
          write(*,*) '******* Attached-CSM radial grid is not strictly increasing at zone ',i
          stop 1
       endif
    enddo

  contains

    subroutine fill_geometric_radii(step0,rout)
      real*8, intent(in) :: step0
      real*8, intent(out) :: rout(n_terms)
      real*8 :: local_growth,local_width
      integer :: j
      local_growth = solve_growth_factor(step0,n_terms,total_span)
      local_width = step0
      rout(1)=surface_radius+local_width
      do j=2,n_terms
         local_width=local_width*local_growth
         rout(j)=rout(j-1)+local_width
      enddo
      rout(n_terms)=outer_radius
    end subroutine fill_geometric_radii

    logical function radii_strictly_increasing(rin) result(ok)
      real*8, intent(in) :: rin(n_terms)
      integer :: j
      ok=.true.
      do j=2,n_terms
         if (rin(j).le.rin(j-1)) then
            ok=.false.
            return
         endif
      enddo
    end function radii_strictly_increasing
  end subroutine build_wind_radii

  subroutine copy_profile_with_append(base_profile,base_mass,base_temp,base_rho, &
       append_mass,append_radius,append_temp,append_rho,append_vel,append_ye, &
       append_omega,new_profile)
    type(csm_profile_data), intent(in) :: base_profile
    real*8, intent(in) :: base_mass(:),base_temp(:),base_rho(:)
    real*8, intent(in) :: append_mass(:),append_radius(:),append_temp(:)
    real*8, intent(in) :: append_rho(:),append_vel(:),append_ye,append_omega
    type(csm_profile_data), intent(out) :: new_profile
    integer :: i,n_old,n_new,n_append

    n_old = size(base_mass)
    n_append = size(append_mass)
    n_new = n_old+n_append
    new_profile%n_zones = n_new
    allocate(new_profile%idx(n_new),new_profile%mass(n_new), &
         new_profile%radius(n_new),new_profile%temp(n_new),new_profile%rho(n_new), &
         new_profile%vel(n_new),new_profile%ye(n_new),new_profile%omega(n_new))

    do i=1,n_old
       new_profile%idx(i)=i
       new_profile%mass(i)=base_mass(i)
       new_profile%radius(i)=base_profile%radius(i)
       new_profile%temp(i)=base_temp(i)
       new_profile%rho(i)=base_rho(i)
       new_profile%vel(i)=base_profile%vel(i)
       new_profile%ye(i)=base_profile%ye(i)
       new_profile%omega(i)=base_profile%omega(i)
    enddo
    do i=1,n_append
       new_profile%idx(n_old+i)=n_old+i
       new_profile%mass(n_old+i)=append_mass(i)
       new_profile%radius(n_old+i)=append_radius(i)
       new_profile%temp(n_old+i)=append_temp(i)
       new_profile%rho(n_old+i)=append_rho(i)
       new_profile%vel(n_old+i)=append_vel(i)
       new_profile%ye(n_old+i)=append_ye
       new_profile%omega(n_old+i)=append_omega
    enddo
  end subroutine copy_profile_with_append

  real*8 function clipped_density_slope(radius,rho,i_anchor) result(slope)
    real*8, intent(in) :: radius(:),rho(:)
    integer, intent(in) :: i_anchor
    integer :: i0,i1
    real*8 :: dx

    i0=max(1,i_anchor-1)
    i1=min(size(radius),i_anchor+1)
    if (i1.le.i0) then
       slope=-3.0d0
       return
    endif
    dx=log(max(radius(i1),1.0d-99))-log(max(radius(i0),1.0d-99))
    if (abs(dx).le.1.0d-99) then
       slope=-3.0d0
    else
       slope=(log(max(rho(i1),1.0d-99))-log(max(rho(i0),1.0d-99)))/dx
    endif
    slope=min(-0.05d0,max(-12.0d0,slope))
  end function clipped_density_slope

  real*8 function geometric_series_sum(first_step,growth,n_terms) result(series_sum)
    real*8, intent(in) :: first_step,growth
    integer, intent(in) :: n_terms
    if (abs(growth-1.0d0).le.1.0d-12) then
       series_sum=first_step*dble(n_terms)
    else
       series_sum=first_step*(growth**n_terms-1.0d0)/(growth-1.0d0)
    endif
  end function geometric_series_sum

  ! bisect the geometric spacing needed to reach the requested radius
  real*8 function solve_growth_factor(first_step,n_terms,total_span) result(growth)
    real*8, intent(in) :: first_step,total_span
    integer, intent(in) :: n_terms
    real*8 :: lower,upper,middle,span,uniform_span
    integer :: iter

    if (n_terms.le.1) then
       growth=1.0d0
       return
    endif
    if (total_span.lt.first_step*(1.0d0-1.0d-12)) then
       write(*,*) '******* First attached-CSM radial step exceeds the total CSM span'
       stop 1
    endif

    uniform_span = first_step*dble(n_terms)
    if (abs(total_span-uniform_span).le.1.0d-12*max(total_span,1.0d0)) then
       growth=1.0d0
       return
    endif

    if (total_span.gt.uniform_span) then
       lower=1.0d0
       upper=1.02d0
       do while (geometric_series_sum(first_step,upper,n_terms).lt.total_span)
          upper=upper*1.1d0
          if (upper.gt.10.0d0) then
             write(*,*) '******* Could not solve for attached-CSM radial growth factor'
             stop 1
          endif
       enddo
    else
       lower=1.0d-8
       upper=1.0d0
       if (geometric_series_sum(first_step,lower,n_terms).gt.total_span) then
          write(*,*) '******* Could not fit compact attached-CSM radial grid'
          stop 1
       endif
    endif

    do iter=1,120
       middle=0.5d0*(lower+upper)
       span=geometric_series_sum(first_step,middle,n_terms)
       if (span.lt.total_span) then
          lower=middle
       else
          upper=middle
       endif
    enddo
    growth=0.5d0*(lower+upper)
  end function solve_growth_factor

end module csm_profile_builders

module csm_switch

  use parameters, only: outdir, profile_name, composition_profile_name, &
       csm, mass_scheme, csm_extension_zones, csm_blend_zones, csm_composition_preset, &
       csm_surface_mix_zones, csm_wind_velocity_kms, csm_mdot_msun_per_yr
  use csm_shared, only: csm_profile_data, csm_composition_data, &
       csm_boundary_diagnostics, read_profile_file, read_composition_file, &
       validate_alignment, write_profile_file, write_composition_file, &
       write_boundary_diagnostics, extend_composition, get_preset_mix, &
       summarize_mix, ensure_directory, lower_trimmed
  use csm_profile_builders, only: build_single_attached
  implicit none

  logical :: csm_generated_active = .false.
  character(len=512) :: csm_base_profile_path = ''
  character(len=512) :: csm_base_composition_path = ''
  character(len=512) :: csm_generated_profile_path = ''
  character(len=512) :: csm_generated_composition_path = ''
  character(len=512) :: csm_generated_directory = ''
  character(len=32) :: csm_mode_active = 'none'
  character(len=32) :: csm_preset_active = 'solar'
  character(len=256) :: csm_mix_summary = ''
  real*8 :: csm_generated_mass_msun = 0.0d0
  real*8 :: csm_generated_surface_mass_g = 0.0d0
  real*8 :: csm_generated_surface_radius_cm = 0.0d0
  real*8 :: csm_generated_outer_radius_cm = 0.0d0
  real*8 :: csm_progenitor_surface_radius_cm = 0.0d0
  integer :: csm_generated_density_transition_zones = 0

contains

  ! prepare generated inputs once before the solver reads its grid
  subroutine prepare_csm_inputs()
    type(csm_profile_data) :: base_profile, new_profile
    type(csm_composition_data) :: base_composition, new_composition
    type(csm_boundary_diagnostics) :: diag
    real*8, allocatable :: mix_aligned(:)
    character(len=32) :: mode_norm, preset_norm

    call reset_csm_state
    csm_base_profile_path = trim(adjustl(profile_name))
    csm_base_composition_path = trim(adjustl(composition_profile_name))

    mode_norm = lower_trimmed(csm)
    if (len_trim(mode_norm).eq.0) mode_norm = 'none'
    preset_norm = lower_trimmed(csm_composition_preset)
    if (len_trim(preset_norm).eq.0) preset_norm = 'solar'
    csm_mode_active = mode_norm
    csm_preset_active = preset_norm

    if (mode_norm.eq.'none') return
    if (mode_norm.ne.'active') then
       write(*,*) '******* Unsupported csm: ', trim(mode_norm)
       write(*,*) '******* Supported states: none, active'
       stop 1
    endif
    if (csm_extension_zones.le.0) then
       write(*,*) '******* csm_extension_zones must be > 0'
       stop 1
    endif
    if (csm_blend_zones.lt.0) then
       write(*,*) '******* csm_blend_zones must be >= 0'
       stop 1
    endif
    if (csm_surface_mix_zones.le.0) then
       write(*,*) '******* csm_surface_mix_zones must be > 0'
       stop 1
    endif

    call read_profile_file(csm_base_profile_path, base_profile)
    call read_composition_file(csm_base_composition_path, base_composition)
    call validate_alignment(base_profile, base_composition)
    csm_progenitor_surface_radius_cm = base_profile%radius(base_profile%n_zones)

    allocate(mix_aligned(base_composition%n_isotopes))
    call get_preset_mix(base_composition, preset_norm, mix_aligned)
    csm_mix_summary = summarize_mix(base_composition%species, mix_aligned)

    call build_single_attached(base_profile, base_composition, mix_aligned, &
         new_profile, diag)
    call extend_composition(base_composition, new_profile, mix_aligned, &
         new_composition)

    csm_generated_directory = trim(adjustl(outdir))//'/generated_csm'
    call ensure_directory(csm_generated_directory)
    csm_generated_profile_path = trim(csm_generated_directory)//'/wdprofile_generated.dat'
    csm_generated_composition_path = trim(csm_generated_directory)//'/wdcomp_generated.dat'

    call write_profile_file(csm_generated_profile_path, new_profile)
    call write_composition_file(csm_generated_composition_path, new_composition)
    call write_boundary_diagnostics(trim(csm_generated_directory)// &
         '/boundary_diagnostics.txt', diag, mode_norm, preset_norm, csm_mix_summary)
    call append_info_log(diag)

    profile_name = trim(csm_generated_profile_path)
    composition_profile_name = trim(csm_generated_composition_path)
    csm_generated_mass_msun = diag%csm_mass_msun
    csm_generated_surface_mass_g = diag%surface_mass_g
    csm_generated_surface_radius_cm = diag%surface_radius_cm
    csm_generated_outer_radius_cm = diag%outer_radius_cm
    csm_generated_density_transition_zones = diag%density_transition_zones
    csm_generated_active = .true.
  end subroutine prepare_csm_inputs

  subroutine reset_csm_state
    csm_generated_active = .false.
    csm_generated_profile_path = ''
    csm_generated_composition_path = ''
    csm_generated_directory = ''
    csm_mode_active = 'none'
    csm_preset_active = 'solar'
    csm_mix_summary = ''
    csm_generated_mass_msun = 0.0d0
    csm_generated_surface_mass_g = 0.0d0
    csm_generated_surface_radius_cm = 0.0d0
    csm_generated_outer_radius_cm = 0.0d0
    csm_progenitor_surface_radius_cm = 0.0d0
    csm_generated_density_transition_zones = 0
  end subroutine reset_csm_state

  subroutine append_info_log(diag)
    type(csm_boundary_diagnostics), intent(in) :: diag
    character(len=1024) :: filepath
    integer :: u

    filepath = trim(adjustl(outdir))//'/info.dat'
    u = 906
    open(unit=u, file=trim(filepath), status='unknown', form='formatted', &
         position='append')
    write(u,'(A)') '--- CSM generated inputs ---'
    write(u,'(A,A)') 'Base profile = ', trim(csm_base_profile_path)
    write(u,'(A,A)') 'Base composition = ', trim(csm_base_composition_path)
    write(u,'(A,A)') 'Generated profile = ', trim(csm_generated_profile_path)
    write(u,'(A,A)') 'Generated composition = ', trim(csm_generated_composition_path)
    write(u,'(A,A)') 'CSM = ', trim(csm_mode_active)
    write(u,'(A,A)') 'Mass scheme = ', trim(adjustl(mass_scheme))
    write(u,'(A,A)') 'CSM composition preset = ', trim(csm_preset_active)
    write(u,'(A,A)') 'CSM mix = ', trim(csm_mix_summary)
    write(u,'(A,ES14.6)') 'CSM mass (Msun) = ', diag%csm_mass_msun
    write(u,'(A,ES14.6)') 'Surface mass (g) = ', diag%surface_mass_g
    write(u,'(A,A)') 'Density smoothing = ', trim(diag%density_smoothing_mode)
    if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
         trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
       write(u,'(A,ES14.6)') 'CSM wind velocity (km/s) = ', csm_wind_velocity_kms
       if (trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
          write(u,'(A,ES14.6)') 'CSM Mdot (Msun/yr) = ', csm_mdot_msun_per_yr
       else
          write(u,'(A)') 'CSM Mdot = Badenes trapezoidal history'
       endif
       write(u,'(A,A)') 'Velocity smoothing = ', trim(diag%velocity_smoothing_mode)
       write(u,'(A,I8)') 'Velocity transition zones = ', diag%velocity_transition_zones
    endif
    close(u)
  end subroutine append_info_log

end module csm_switch
