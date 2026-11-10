module outinfomod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: out_info_string_every = 200
  integer :: outinfo_count

contains

  subroutine locate_live_csm_interaction(front_idx,qmax_out,qratio_out)
    use blmod, only: Q,p,mass
    use parameters, only: imax,csm_blend_zones
    use csm_switch, only: csm_generated_active,csm_generated_surface_mass_g
    implicit none

    integer, intent(out) :: front_idx
    real*8, intent(out) :: qmax_out,qratio_out
    integer :: i,i_surface,i0,strongest_idx
    real*8 :: qval,qratio

    i_surface=1
    if (csm_generated_active .and. csm_generated_surface_mass_g.gt.0.0d0) then
       do i=1,imax
          if (mass(i).ge.csm_generated_surface_mass_g) then
             i_surface=i
             exit
          endif
       enddo
    else
       i_surface=max(1,imax/2)
    endif

    i0=max(1,min(imax-2,i_surface-max(8,csm_blend_zones)))
    front_idx=i_surface
    strongest_idx=i_surface
    qmax_out=0.0d0
    qratio_out=0.0d0

    do i=i0,imax-2
       qval=abs(Q(i))
       if (qval.gt.qmax_out) then
          qmax_out=qval
          strongest_idx=i
       endif
    enddo

    if (qmax_out.gt.1.0d-99) then
       front_idx=strongest_idx
       qratio_out=abs(Q(strongest_idx))/ &
            (abs(Q(strongest_idx))+abs(p(strongest_idx))+1.0d-99)
       do i=strongest_idx,imax-2
          qval=abs(Q(i))
          qratio=qval/(qval+abs(p(i))+1.0d-99)
          if (qval.ge.1.0d-3*qmax_out .and. qratio.ge.5.0d-2) then
             front_idx=i
             qratio_out=qratio
          endif
       enddo
    endif

    front_idx=max(1,min(imax-1,front_idx))
  end subroutine locate_live_csm_interaction


  subroutine locate_live_shock(front_idx,front_radius,shock_strength)
    use blmod, only: cr
    use parameters, only: imax
    implicit none

    integer, intent(out) :: front_idx
    real*8, intent(out) :: front_radius,shock_strength
    real*8 :: qmax,qratio

    call locate_live_csm_interaction(front_idx,qmax,qratio)
    shock_strength=qratio
    if (qmax.le.1.0d-99 .or. qratio.lt.1.0d-3) then
       front_idx=max(1,imax-1)
       shock_strength=0.0d0
    endif
    front_idx=max(1,min(imax-1,front_idx))

    front_radius=cr(front_idx)
    if (.not.ieee_is_finite(front_radius) .or. front_radius.lt.0.0d0) front_radius=0.0d0
    if (.not.ieee_is_finite(shock_strength)) shock_strength=0.0d0
    shock_strength=max(0.0d0,min(1.0d0,shock_strength))
  end subroutine locate_live_shock

end module outinfomod


subroutine outinfo
  use blmod, only: nt,time,dtime
  use outinfomod, only: outinfo_count,out_info_string_every,locate_live_shock
  implicit none

  real*8 :: shock_radius,shock_strength
  integer :: shock_idx

  ! report the current run state
  call locate_live_shock(shock_idx,shock_radius,shock_strength)

  if (mod(outinfo_count,out_info_string_every).eq.0) then
     outinfo_count=0
     write(6,'(A8,2A14,2A16)') 'nt','time','dtime','shock_r(cm)','shock_strength'
  endif
  write(*,'(I8,1P4E16.6)') nt,time,dtime,shock_radius,shock_strength

  outinfo_count=outinfo_count+1
end subroutine outinfo
