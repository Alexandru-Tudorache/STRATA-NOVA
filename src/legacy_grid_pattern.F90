module legacy_grid_pattern
  implicit none
  private
  public :: build_legacy_mass_grid, build_legacy_unit_pattern

contains

  ! retain the inherited nonuniform mass-grid spacing
  subroutine build_legacy_unit_pattern(n, pattern)
    integer, intent(in) :: n
    real*8, intent(out) :: pattern(n)

    integer :: i, itran
    real*8 :: ratio1, ratio2, f1, f2, delta_tran
    real*8, allocatable :: delta(:)

    if (n.lt.5) then
       write(*,*) '******* legacy_native grid requires imax >= 5'
       stop 1
    endif

    if (n.ge.22) then
       itran = 20
    else
       itran = max(3,min(n-2,n/3+1))
    endif

    ratio1 = 0.1d0
    ratio2 = 0.01d0
    f1 = ratio1**(1.0d0/dble(itran-2))
    f2 = ratio2**(1.0d0/dble(n-itran-1))

    delta_tran = ((1.0d0-f2)*(1.0d0-f1)) / &
         ((1.0d0-f2)*(1.0d0-f1**(itran-1)) + &
          (1.0d0-f1)*(1.0d0-f2**(n-itran)))

    allocate(delta(n-1))
    do i=1,itran-1
       delta(i) = delta_tran*f1**(itran-1-i)
    enddo
    do i=itran,n-1
       delta(i) = delta_tran*f2**(i-itran)
    enddo

    pattern(1)=0.0d0
    do i=2,n
       pattern(i)=pattern(i-1)+delta(i-1)
    enddo
    pattern(n)=1.0d0

    do i=1,n-1
       if (pattern(i+1).le.pattern(i)) then
          write(*,*) '******* legacy_native pattern is non-monotone at ',i
          stop 1
       endif
    enddo

    deallocate(delta)
  end subroutine build_legacy_unit_pattern

  subroutine build_legacy_mass_grid(mass,n,mstart,mend)
    integer, intent(in) :: n
    real*8, intent(in) :: mstart,mend
    real*8, intent(out) :: mass(n)
    real*8, allocatable :: pattern(:)
    integer :: i

    if (mend.le.mstart) then
       write(*,*) '******* legacy_native grid requires mend > mstart'
       stop 1
    endif

    allocate(pattern(n))
    call build_legacy_unit_pattern(n,pattern)
    do i=1,n
       mass(i)=mstart+pattern(i)*(mend-mstart)
    enddo
    mass(1)=mstart
    mass(n)=mend
    deallocate(pattern)
  end subroutine build_legacy_mass_grid

end module legacy_grid_pattern
