subroutine nickel_heating

  use blmod, only: ye, comp, Ni_heating, Ni_total_luminosity, Ni_number, &
                   r, cr, rho, time, Ni_deposit_function, Ni_energy_rate, delta_mass
  use parameters
  use physical_constants
  implicit none


  real*8 :: r_Ni
  real*8 :: delta_r
  real*8 :: r_x, r_max
  real*8 :: r_j, comp_Ni_j
  real*8 :: delta_tau_j

  real*8 :: th, delta_th
  real*8 :: th_min(imax)

  real*8 :: I_prime, I_prime_av

  integer :: i, i_Ni, i_theta, i_radius

  integer :: index, index_hint

  real*8 :: kappa_gamma(imax)

  real*8 :: interp_frac
  integer :: idx_lo, idx_hi

  real*8, parameter :: nimin = 1.0d-5
  real*8, parameter :: th_max = pi

!------------------------------------------------------------------------------
!based on the work of Swartz et al., ApJ 446:766 (1995)


  if(Ni_switch.eq.0) then

     Ni_heating(:) = 0.0d0
     Ni_deposit_function(:) = 0.0d0
     Ni_total_luminosity = 0.0d0

  else

     ! find the outermost cell containing appreciable ni
     i_Ni = 0
     do i=imax, 1, -1
        if(comp(i,Ni_number).gt.nimin) then
           r_Ni = r(i)
           i_Ni = i
           exit
        endif
     enddo

     ! restrict rays to directions that intersect the ni-bearing region
     if(i_Ni.eq.0) then
        write(*,*) 'mass fraction of Ni is lower than NIMIN in every grid point'
        write(*,*) 'try reducing NIMIN'
        stop 1
     else if(i_Ni.eq.1) then
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     else if(i_Ni.eq.imax) then
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
     else
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     end if


     ! cache gamma-ray extinction for the optimized tracing path
     if (ni_raytrace_opt .ge. 1) then
        do i=1, imax
           kappa_gamma(i) = ye(i) * 0.06d0 * rho(i)
        end do
     endif

     ! average the attenuated ni source over rays at each grid point
     do i=1, imax
        I_prime_av = 0
        delta_th = (th_max-th_min(i))/Ni_quad_npoints
        index_hint = 1

        do i_theta=0,Ni_quad_npoints-1
           th = th_max-dble(i_theta)*delta_th
           r_max = -r(i)*cos(th) + sqrt((r(i)*cos(th))**2-(r(i)**2-r_Ni**2))
           delta_r = r_max/Ni_quad_npoints
           I_prime = 0

           do i_radius=Ni_quad_npoints,1,-1
              r_x = dble(i_radius)*delta_r
              r_j = sqrt(r(i)*r(i) + r_x*r_x + 2.0d0*r(i)*r_x*cos(th))

              if(r_j.le.r(1)) then
                 delta_tau_j = 0.0d0
                 comp_Ni_j = 0.0d0
              else if(r_j.ge.r(imax-1)) then
                 if (ni_raytrace_opt .ge. 1) then
                    delta_tau_j = delta_r * kappa_gamma(imax-1)
                 else
                    delta_tau_j = delta_r * ye(imax-1) * 0.06d0 * rho(imax-1)
                 endif
                 comp_Ni_j = comp(imax-1,Ni_number)
              else
                 if (ni_raytrace_opt .ge. 1) then
                    do while (index_hint > 1 .and. r_j .le. r(index_hint))
                       index_hint = index_hint - 1
                    end do
                    do while (index_hint < imax-1 .and. r_j .gt. r(index_hint+1))
                       index_hint = index_hint + 1
                    end do
                    index = index_hint
                    if (ni_ray_interp .ge. 1) then
                       ! interpolate using cell centres rather than radial boundaries
                       if (r_j .lt. cr(index) .and. index .gt. 1) then
                          idx_lo = index - 1
                          idx_hi = index
                       else if (index .lt. imax-1) then
                          idx_lo = index
                          idx_hi = index + 1
                       else
                          idx_lo = index
                          idx_hi = index
                       endif
                       if (idx_lo .eq. idx_hi) then
                          delta_tau_j = delta_r * kappa_gamma(index)
                          comp_Ni_j = comp(index,Ni_number)
                       else
                          interp_frac = (r_j - cr(idx_lo)) / (cr(idx_hi) - cr(idx_lo))
                          interp_frac = max(0.0d0, min(1.0d0, interp_frac))
                          delta_tau_j = delta_r * ((1.0d0-interp_frac)*kappa_gamma(idx_lo) &
                               + interp_frac*kappa_gamma(idx_hi))
                          comp_Ni_j = (1.0d0-interp_frac)*comp(idx_lo,Ni_number) &
                               + interp_frac*comp(idx_hi,Ni_number)
                       endif
                    else
                       delta_tau_j = delta_r * kappa_gamma(index)
                       comp_Ni_j = comp(index,Ni_number)
                    endif
                 else
                    index = 1
                    do while (index < imax-1 .and. r_j .gt. r(index+1))
                       index = index + 1
                    end do
                    if (ni_ray_interp .ge. 1) then
                       if (r_j .lt. cr(index) .and. index .gt. 1) then
                          idx_lo = index - 1
                          idx_hi = index
                       else if (index .lt. imax-1) then
                          idx_lo = index
                          idx_hi = index + 1
                       else
                          idx_lo = index
                          idx_hi = index
                       endif
                       if (idx_lo .eq. idx_hi) then
                          delta_tau_j = delta_r * ye(index)*0.06d0*rho(index)
                          comp_Ni_j = comp(index,Ni_number)
                       else
                          interp_frac = (r_j - cr(idx_lo)) / (cr(idx_hi) - cr(idx_lo))
                          interp_frac = max(0.0d0, min(1.0d0, interp_frac))
                          delta_tau_j = delta_r * ( &
                               (1.0d0-interp_frac)*ye(idx_lo)*0.06d0*rho(idx_lo) &
                               + interp_frac*ye(idx_hi)*0.06d0*rho(idx_hi))
                          comp_Ni_j = (1.0d0-interp_frac)*comp(idx_lo,Ni_number) &
                               + interp_frac*comp(idx_hi,Ni_number)
                       endif
                    else
                       delta_tau_j = delta_r * ye(index) * 0.06d0 * rho(index)
                       comp_Ni_j = comp(index,Ni_number)
                    endif
                 endif
              end if

              ! advance the source through one attenuating ray segment
              I_prime = (I_prime-comp_Ni_j)*exp(-min(max(delta_tau_j,0.0d0),700.0d0)) + comp_Ni_j
           end do

           I_prime_av = I_prime_av + I_prime*sin(th)*delta_th*0.5d0
        end do

        Ni_deposit_function(i) = I_prime_av
     end do

     ! evaluate decay power using the absolute model age
     Ni_energy_rate = &
          3.24d10*exp(-(time+initial_model_age_s)*overtau_Ni) + &
          7.29d9*exp(-(time+initial_model_age_s)*overtau_Co)

     Ni_heating(1:imax) = Ni_energy_rate*Ni_deposit_function(1:imax)

     Ni_total_luminosity = &
          Ni_energy_rate*sum(Ni_deposit_function(1:imax)*delta_mass(1:imax))

  end if

end subroutine nickel_heating
