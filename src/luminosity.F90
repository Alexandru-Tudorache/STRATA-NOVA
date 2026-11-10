
subroutine luminosity(r_x,temp_x,kappa_x,lambda_x,inv_kappa_x,lum_x)

  use blmod, only: delta_cmass
  use parameters
  use physical_constants
  implicit none

  real*8 :: r_x(imax)
  real*8 :: temp_x(imax)
  real*8 :: kappa_x(imax)

  real*8 :: lambda_x(imax)
  real*8 :: inv_kappa_x(imax)
  real*8 :: lum_x(imax)

  integer :: i
  real*8 :: bigr_x(imax)


  ! close the inner radiative boundary
  inv_kappa_x(1) = 1.0d0
  bigr_x(1) = 0.0d0
  lambda_x(1) = 1.0d0
  lum_x(1) = 0.0d0

  ! construct face opacities and the flux limiter from neighbouring cells
  do i=2, imax

    inv_kappa_x(i) = (temp_x(i)**4/kappa_x(i) + temp_x(i-1)**4/kappa_x(i-1)) &
        / (temp_x(i)**4 + temp_x(i-1)**4)

    bigr_x(i) = 8.0d0*pi*r_x(i)**2 * abs(temp_x(i)**4-(temp_x(i-1))**4) &
         * inv_kappa_x(i) / (delta_cmass(i-1)*(temp_x(i)**4+temp_x(i-1)**4))

    lambda_x(i) = (6.0d0 + 3.0d0*bigr_x(i)) &
        /(6.0d0 + 3.0d0*bigr_x(i) + bigr_x(i)*bigr_x(i))

  end do

  do i=2, imax

    ! convert the limited diffusion gradient into spherical luminosity
    lum_x(i) = -( 4.0d0*pi*r_x(i)**2 )**2  &
        * a_rad*clite*lambda_x(i)*inv_kappa_x(i)/3.0d0  &
        * (temp_x(i)**4-temp_x(i-1)**4)/delta_cmass(i-1)

  end do


end subroutine luminosity


subroutine optical_depth(rho, r, kappa, tau)

  use parameters
  implicit none

  real*8 rho(imax)
  real*8 r(imax)
  real*8 kappa(imax)

  real*8 tau(imax)

  integer i


  ! integrate optical depth inward from the vacuum exterior
  tau(imax) = 0

  do i=imax-1, 1, -1
      tau(i) = tau(i+1) + (r(i+1) - r(i))*rho(i)*(kappa(i+1) + kappa(i))/2.0d0
  enddo

end subroutine optical_depth


