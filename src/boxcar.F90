! boxcar smoothing with kernel rehaul from SuperSNEC

!-------------------------------------------------------------------------------

! compute boxcar profile
subroutine boxcar_profile(mass_coords, abundances, zones, nspecies, &
     boxcar_nominal_mass, number_iterations)

  use mixing_utils, only: build_effective_delta_mass
  implicit none

  integer, intent(in) :: zones, nspecies
  real*8, intent(in) :: mass_coords(zones)
  real*8, intent(inout) :: abundances(zones, nspecies)
  real*8, intent(in) :: boxcar_nominal_mass
  integer, intent(in) :: number_iterations
  real*8, allocatable :: local_delta(:)

  if (zones <= 1 .or. nspecies <= 0) return ! stop code if not enough zones or composition species

  allocate(local_delta(zones))
  call build_effective_delta_mass(zones, mass_coords, local_delta)
  call boxcar_kernel(abundances, local_delta, zones, nspecies, &
       boxcar_nominal_mass, number_iterations)
  deallocate(local_delta)

end subroutine boxcar_profile

! compute mixing and conserve mass fractions
!-------------------------------------------------------------------------------
subroutine boxcar_kernel(abundances, delta_mass, zones, nspecies, &
     boxcar_nominal_mass, number_iterations)
  implicit none

  integer, intent(in) :: zones, nspecies
  real*8, intent(inout) :: abundances(zones, nspecies)
  real*8, intent(in) :: delta_mass(zones)
  real*8, intent(in) :: boxcar_nominal_mass
  integer, intent(in) :: number_iterations

  real*8 :: boxcar_actual_mass
  real*8 :: mass_element
  integer :: i, l, l_max, k, n
  integer :: el

  ! same invalid setup protection logic
  if (zones <= 0 .or. nspecies <= 0) return
  if (boxcar_nominal_mass <= 0.0d0 .or. number_iterations <= 0) return

  do n=1, number_iterations ! for each boxcar iteration

      do i=1, zones ! for each cell

          l_max = zones
          do l=i, zones
          ! sum of cells containing greater mass than boxcar mass  
              if( sum(delta_mass(i:l)).gt.boxcar_nominal_mass ) then
                l_max = l
                exit
              endif
              if( l.eq.zones ) then
                l_max = zones ! outer boundary
                exit
              endif
          enddo

          boxcar_actual_mass = sum(delta_mass(i:l_max))
          ! repeat until no longer less than or equal to zero; protection against invalid averaging
          if (boxcar_actual_mass <= 0.0d0) cycle

          ! mass fraction conservation
          do el=1, nspecies
              mass_element = sum(delta_mass(i:l_max)*abundances(i:l_max,el))
              do k = i, l_max
                  abundances(k,el) = mass_element/boxcar_actual_mass
              enddo
          enddo

          if( boxcar_actual_mass .lt. boxcar_nominal_mass ) exit

      enddo

  enddo

end subroutine boxcar_kernel
