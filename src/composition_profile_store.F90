
module composition_profile_store
  use mixing_utils, only: build_effective_delta_mass
  implicit none
  real*8, parameter :: tiny_mass_interval = 1.0d-30
  integer :: stored_profile_zones = 0
  integer :: stored_profile_species = 0
  real*8, allocatable :: stored_profile_mass(:)
  real*8, allocatable :: stored_profile_comp(:,:)
  real*8, allocatable :: stored_profile_delta(:)
  real*8, allocatable :: stored_profile_edges(:)
contains

  subroutine store_processed_profile(pmass, pcomp, zones, species, left_edge_semantics)

    implicit none
    integer, intent(in) :: zones, species
    real*8, intent(in) :: pmass(zones)
    real*8, intent(in) :: pcomp(zones,species)
    logical, intent(in), optional :: left_edge_semantics
    logical :: use_left_edges

    ! replace the previous in-memory composition
    call release_stored_profile()

    if (zones <= 0 .or. species <= 0) then
       stored_profile_zones = 0
       stored_profile_species = 0
       return
    endif

    allocate(stored_profile_mass(zones))
    allocate(stored_profile_comp(zones,species))
    allocate(stored_profile_delta(zones))
    allocate(stored_profile_edges(zones+1))
    stored_profile_mass(1:zones) = pmass(1:zones)
    stored_profile_comp(1:zones,1:species) = pcomp(1:zones,1:species)

    ! restart rows may describe intervals beginning at pmass(i)
    use_left_edges = .false.
    if (present(left_edge_semantics)) use_left_edges = left_edge_semantics
    if (use_left_edges) then
       call build_left_edge_profile(zones)
    else
       call build_effective_delta_mass(zones, pmass, stored_profile_delta)
       call build_profile_edges(zones)
    endif
    stored_profile_zones = zones
    stored_profile_species = species

  end subroutine store_processed_profile


  subroutine release_stored_profile()

    if (allocated(stored_profile_mass)) deallocate(stored_profile_mass)
    if (allocated(stored_profile_comp)) deallocate(stored_profile_comp)
    if (allocated(stored_profile_delta)) deallocate(stored_profile_delta)
    if (allocated(stored_profile_edges)) deallocate(stored_profile_edges)
    stored_profile_zones = 0
    stored_profile_species = 0

  end subroutine release_stored_profile


  ! check that the stored arrays exist
  logical function has_stored_profile()

    has_stored_profile = stored_profile_zones > 0 .and. stored_profile_species > 0 &
         & .and. allocated(stored_profile_mass) .and. allocated(stored_profile_comp) &
         & .and. allocated(stored_profile_edges)

  end function has_stored_profile


  ! map cell averages by their overlap in mass coordinate
  subroutine map_stored_profile_to_grid(target_edges, target_points, nspecies, comp_out)

    implicit none
    integer, intent(in) :: target_points, nspecies
    real*8, intent(in) :: target_edges(target_points)
    real*8, intent(out) :: comp_out(target_points,nspecies)
    integer :: species, new_idx
    integer :: old_idx
    real*8 :: left_edge, right_edge, denom, accumulator
    real*8 :: segment_end, overlap
    real*8 :: clip_left, clip_right
    real*8 :: outer_edge

    if (.not.has_stored_profile()) return
    if (nspecies > stored_profile_species) then
       write(*,*) 'ERROR: requested species exceed stored profile content'
       stop 1
    endif

    if (target_points <= 0) return

    outer_edge = stored_profile_edges(stored_profile_zones+1)

    ! average each species over the mass shared by old and new cells
    do species=1,nspecies
       old_idx = 1
       do new_idx=1,target_points-1
          left_edge = target_edges(new_idx)
          right_edge = target_edges(new_idx+1)

          ! use the nearest composition outside the stored mass range
          if (right_edge <= stored_profile_edges(1)) then
             comp_out(new_idx,species) = stored_profile_comp(1,species)
             cycle
          endif
          if (left_edge >= outer_edge) then
             comp_out(new_idx,species) = stored_profile_comp(stored_profile_zones,species)
             cycle
          endif

          ! restrict the average to the covered part of this cell
          clip_left = max(left_edge, stored_profile_edges(1))
          clip_right = min(right_edge, outer_edge)
          denom = max(tiny_mass_interval, clip_right - clip_left)
          accumulator = 0.0d0
          left_edge = clip_left
          right_edge = clip_right

          ! advance to the first old cell overlapping this new interval
          do while (old_idx < stored_profile_zones .and. &
               stored_profile_edges(old_idx+1) <= left_edge + tiny_mass_interval)
             old_idx = old_idx + 1
          enddo

          ! accumulate abundance times overlap mass
          do
             if (left_edge >= right_edge - tiny_mass_interval) exit
             segment_end = min(right_edge, stored_profile_edges(old_idx+1))
             overlap = max(0.0d0, segment_end - left_edge)
             accumulator = accumulator + stored_profile_comp(old_idx,species) * overlap
             left_edge = segment_end
             if (left_edge >= right_edge - tiny_mass_interval) exit
             if (old_idx < stored_profile_zones) old_idx = old_idx + 1
          enddo

          comp_out(new_idx,species) = accumulator/denom
       enddo

       ! copy the last physical composition into the boundary row
       if (target_points .ge. 2) then
          comp_out(target_points,species) = comp_out(target_points-1,species)
       else
          comp_out(target_points,species) = stored_profile_comp(1,species)
       endif
    enddo

  end subroutine map_stored_profile_to_grid


  ! treat each mass coordinate as the beginning of its abundance interval
  subroutine build_left_edge_profile(zones)

    implicit none
    integer, intent(in) :: zones
    integer :: i
    real*8 :: last_width

    if (zones <= 0) return
    ! use the supplied left edges and extend the final interval
    do i=1,zones
       stored_profile_edges(i)=stored_profile_mass(i)
    enddo
    if (zones > 1) then
       last_width=max(tiny_mass_interval,stored_profile_mass(zones)-stored_profile_mass(zones-1))
    else
       last_width=max(tiny_mass_interval,abs(stored_profile_mass(1)))
    endif
    stored_profile_edges(zones+1)=stored_profile_mass(zones)+last_width
    do i=1,zones-1
       stored_profile_delta(i)=max(tiny_mass_interval,stored_profile_mass(i+1)-stored_profile_mass(i))
    enddo
    stored_profile_delta(zones)=last_width

  end subroutine build_left_edge_profile


  ! infer interval boundaries for centre-labelled source profiles
  subroutine build_profile_edges(zones)
    implicit none
    integer, intent(in) :: zones
    integer :: i

    if (zones <= 0) return

    ! construct mass intervals from the centre-coordinate convention
    stored_profile_edges(1) = stored_profile_mass(1)
    if (zones >= 1) then
       stored_profile_edges(1) = stored_profile_mass(1) - 0.5d0*stored_profile_delta(1)
    endif

    do i=1,zones
       stored_profile_edges(i+1) = stored_profile_edges(i) + stored_profile_delta(i)
    enddo

    stored_profile_edges(1) = min(stored_profile_edges(1), stored_profile_mass(1))
    stored_profile_edges(zones+1) = max(stored_profile_edges(zones+1), stored_profile_mass(zones))

    do i=2,zones+1
       if (stored_profile_edges(i) <= stored_profile_edges(i-1)) then
          stored_profile_edges(i) = stored_profile_edges(i-1) + tiny_mass_interval
       endif
    enddo

  end subroutine build_profile_edges

end module composition_profile_store
