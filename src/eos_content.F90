subroutine ideal_eos(rhox,tempx,px,ex,cs2x,dpdtx,dedtx,np)

  use blmod, only: eos_gamma1
  use physical_constants
  implicit none

  integer :: np
  real*8 :: rhox(np),tempx(np)

  real*8 :: px(np),ex(np),cs2x(np),dpdtx(np),dedtx(np)

  integer :: i


  do i=1, np

    px(i) = (kboltz*avo_real)*rhox(i)*tempx(i)

    ex(i) = (kboltz*avo_real)/(eos_gamma1-1.0d0)*tempx(i)

    cs2x(i) = eos_gamma1*(kboltz*avo_real)*tempx(i)

    dpdtx(i) = (kboltz*avo_real)*rhox(i)

    dedtx(i) = (kboltz*avo_real)/(eos_gamma1-1.0d0)

  end do

end subroutine ideal_eos


subroutine paczynski_eos(rhox,tempx,yex,abarx,px,ex,cs2x,dpdtx,dedtx, &
    pradx,k,np)

  use blmod, only: free_electron_frac, ion_fractions, comp_details, comp, ncomps
  use parameters
  use physical_constants
  use eosmodule, only: xxip, ionization_energy
  implicit none

  integer :: np
  integer :: k(np)
  real*8 :: rhox(np),tempx(np),yex(np),abarx(np)

  real*8 :: px(np),ex(np),cs2x(np),dpdtx(np),dedtx(np),pradx(np)

  real*8 :: N
  real*8 :: ybar
  real*8 :: pion, pe, ped, pend, pednr, pedr
  real*8 :: ne

  real*8 :: invrhox

  real*8 :: e_ioncorr
  real*8 :: ffactor

  real*8 :: sum1, sum2, sum3, sum_pot_frac
  real*8 :: ymax
  real*8 :: zion
  real*8 :: y_r, chi_r, chi_T, chi_rho, c_V
  real*8 :: nu_j
  real*8 :: Gamma1
  real*8 :: p_matter, e_matter, dpdt_matter, dedt_matter

  integer :: i,j
  integer :: l,lmax


!------------------------------------------------------------------------------
!variables ybar, N, y_r, chi_r, chi_T, chi_rho, c_V, nu_j, Gamma1 are named as 
!the analogous quantities in paragraph 9.18 of "Cox & Giuli's principles of
!stellar structure" by Weiss et al. (2004). However, the quantities presented
!here are obtained without assuming that only one element at a time changes
!ionization state.
!The degenerate electron gas is treated as 
!in Paczynski B., ApJ, 267:315 (1983).

  do i = 1, np

    ! solve the ion populations before assembling the pressure and energy
    call simple_saha(saha_ncomps,k(i),tempx(i),rhox(i),ne)


    N = 1.0d0 / (mproton * abarx(i))

    invrhox = 1.0d0/rhox(i)
    ybar = ne / N * invrhox

    sum1 = 0.0d0
    sum2 = 0.0d0
    sum3 = 0.0d0
    e_ioncorr = 0.0d0

    do j=1, ncomps

      zion = comp_details(j,2)

      if(zion.lt.(1.0d0-1.0d-6)) cycle

      sum_pot_frac = 0.0d0

      ! start from the neutral state so the ionization lookup stays in bounds
      lmax = 1
      ymax = ion_fractions(j,1,k(i))
      do l=2, int(zion) + 1
        if(ion_fractions(j,l,k(i)).gt.ymax) then
            ymax = ion_fractions(j,l,k(i))
            lmax = l
        end if
        sum_pot_frac = sum_pot_frac + ion_fractions(j,l,k(i)) * ionization_energy(int(zion),l-1)
      enddo

      if(lmax.eq.1) then
        y_r = ion_fractions(j,lmax+1,k(i))
        chi_r = xxip(int(zion),lmax)
      else if(lmax.eq.(int(zion) + 1)) then
        y_r = ion_fractions(j,lmax,k(i))
        chi_r = xxip(int(zion),lmax-1)
      else
        if(ion_fractions(j,lmax+1,k(i)).gt.ion_fractions(j,lmax-1,k(i))) then
            y_r = ion_fractions(j,lmax + 1,k(i))
            chi_r = xxip(int(zion),lmax)
        else
            y_r = ion_fractions(j,lmax,k(i))
            chi_r = xxip(int(zion),lmax-1)
        end if
      end if

      ! convert each mass fraction into a number abundance
      nu_j = abarx(i) * comp(k(i),j) / comp_details(j,1)

      sum1 = sum1 + nu_j * y_r * (1 - y_r)

      sum2 = sum2 + nu_j * chi_r * y_r * (1 - y_r)

      sum3 = sum3 + nu_j * chi_r * chi_r * y_r * (1 - y_r)

      e_ioncorr = e_ioncorr + N * nu_j * sum_pot_frac

    enddo

    ! combine ion pressure with nondegenerate and degenerate electron terms
    pion = N * rhox(i) * kboltz * tempx(i)

    pend = ybar * pion
    pednr = kdnr * (yex(i)*rhox(i))**(fivethirds)
    pedr = kdr * (yex(i)*rhox(i))**(fourthirds)
    ped = (1.0d0/(pednr*pednr)+1.0d0/(pedr*pedr))**(-0.5d0)
    pe = (pend*pend + ped*ped)**(0.5d0)

    pradx(i) = a_rad * tempx(i)**4 * overthree

    ffactor = (fivethirds)*ped*ped/(pednr*pednr) &
              +(fourthirds)*ped*ped/(pedr*pedr)

    p_matter = pion + pe

    e_matter = 1.5d0 * pion * invrhox &
              + pe / (rhox(i) * (ffactor - 1.0d0)) + e_ioncorr

    dpdt_matter = pion / tempx(i) &
          + pend*pend / (pe * tempx(i)) * (1.0d0 &
              + ( 1.5d0*sum1 + sum2/(kboltz*tempx(i)) )/(ybar + sum1) )

    dedt_matter = 1.5d0 * N * kboltz &
          + pend*pend / (pe * tempx(i) * rhox(i) * (ffactor - 1.0d0)) * (1.0d0 &
              + ( 1.5d0*sum1 + sum2/(kboltz*tempx(i)) )/(ybar + sum1) )        &
          + N * ( sum2 * ( 1.5d0*ybar - sum2/(kboltz*tempx(i)) )/(ybar + sum1) &
              + sum3/(kboltz*tempx(i)) ) / tempx(i)

    ! include lte radiation in the total pressure and specific energy
    px(i) = p_matter + pradx(i)
    ex(i) = e_matter + 3.0d0 * pradx(i) * invrhox
    dpdtx(i) = dpdt_matter + 4.0d0 * pradx(i) / tempx(i)
    dedtx(i) = dedt_matter + 4.0d0 * a_rad * tempx(i)**3 * invrhox

    chi_T = tempx(i) * dpdtx(i) / px(i)

    chi_rho = rhox(i) / px(i) * ( N * kboltz * tempx(i)   &
       + (ybar * pend*pend /(ybar + sum1) + ffactor * ped*ped )/(pe * rhox(i)) )

    c_V = dedtx(i)

    ! derive the adiabatic sound speed from the eos derivatives
    Gamma1 = chi_T*chi_T * px(i) / (c_V * rhox(i) * tempx(i)) + chi_rho

    cs2x(i) = Gamma1 * px(i) * invrhox

    free_electron_frac(k(i)) = ne * mproton / (yex(i) * rhox(i))

  end do

end subroutine paczynski_eos
