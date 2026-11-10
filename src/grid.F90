subroutine grid

  use blmod, only: mass, cmass, delta_mass, delta_cmass
  use parameters, only: imax, grid_mode
  use legacy_grid_pattern, only: build_legacy_mass_grid
  implicit none

  character(len=32) :: mode
  real*8 :: mass_start, mass_end

  if (imax.le.1) return

  mode = normalize_mode(grid_mode)
  mass_start = mass(1)
  mass_end = mass(imax)

  if (mode.eq.'legacy_native' .or. mode.eq.'legacy_pattern') then
     call build_legacy_mass_grid(mass,imax,mass_start,mass_end)
  else
     write(*,*) '******* Unsupported grid_mode: ',trim(grid_mode)
     write(*,*) '******* Production mode is legacy_native (legacy_pattern is an alias).'
     stop 1
  endif

  call finalize_differences()

contains

  function normalize_mode(raw) result(clean)
    character(len=*), intent(in) :: raw
    character(len=len(raw)) :: clean
    integer :: i, code
    clean = trim(adjustl(raw))
    do i=1,len(clean)
       code=iachar(clean(i:i))
       if (code.ge.iachar('A') .and. code.le.iachar('Z')) clean(i:i)=achar(code+32)
    enddo
  end function normalize_mode

  ! derive cell masses and staggered mass intervals from grid boundaries
  subroutine finalize_differences()
    integer :: i
    do i=1,imax-1
       cmass(i)=mass(i)+0.5d0*(mass(i+1)-mass(i))
    enddo
    cmass(imax)=mass(imax)+0.5d0*(mass(imax)-mass(imax-1))
    do i=1,imax-1
       delta_mass(i)=mass(i+1)-mass(i)
       delta_cmass(i)=cmass(i+1)-cmass(i)
    enddo
    delta_mass(imax)=delta_mass(imax-1)
    delta_cmass(imax)=delta_cmass(imax-1)
  end subroutine finalize_differences

end subroutine grid
