subroutine artificial_viscosity

  use blmod, only: vel, rho, Q, cr, cr_p, delta_mass, dtime, Qterm, cs2, mass
  use parameters
  use physical_constants
  use csm_switch, only: csm_generated_active, csm_generated_surface_mass_g
  implicit none

  integer :: i
  real*8, parameter :: cvisc = 2.0d0
  real*8, parameter :: cvisc_linear_constant_shell = 4.0d0
  real*8, parameter :: cvisc_linear_wind = 8.0d0
  real*8 :: dv_comp, sound

! Richtmeyer-Von Neuman artificial viscosity

  do i=1,imax-1
      if(vel(i+1).lt.vel(i)) then
          dv_comp = vel(i) - vel(i+1)
          Q(i) = cvisc*rho(i)*dv_comp**2

          ! three different viscous terms are used for no-csm runs and csm runs 
          ! this helps with stability when the csm medium is compressed by a violent shockwave
          ! these are Nuemann-Richtmyer terms (see Neumann+Richtmeyer, 1950)
          ! addtional testing of arbitrary choices for csm-specific constants still needed

          ! csm-specific selector
          ! continous_mass gets cvisc = 2.0
          ! constant_density gets cvisc = 4.0
          ! wind mass-loss schemes get cvisc = 8.0
          
          if (csm_generated_active .and. &
               (trim(adjustl(mass_scheme)).eq.'constant_density' .or. &
                trim(adjustl(mass_scheme)).eq.'constant_mass_loss' .or. &
                trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss') .and. &
               mass(i+1).gt.csm_generated_surface_mass_g) then
             sound = sqrt(max(cs2(i),0.0d0)) ! +ve sound speed only
             if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
                  trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
                Q(i) = Q(i) + cvisc_linear_wind*rho(i)*sound*dv_comp
             else
                Q(i) = Q(i) + cvisc_linear_constant_shell*rho(i)*sound*dv_comp 
             endif
          endif
      else
          Q(i) = 0
      endif
  enddo

  ! outer boundary
  Q(imax) = 0.0d0

  ! the term with viscosity at the right hand side of the energy equation
  do i=1,imax-1
      Qterm(i) = - 4.0d0*pi*dtime* (0.5d0 * (cr(i) + cr_p(i)))**2 &
          * Q(i) *(vel(i+1) - vel(i))/delta_mass(i)
  enddo

  ! once again, a hard zero boundary
  Qterm(imax) = 0.0d0

end subroutine artificial_viscosity
