  ! iterate the ion populations and the free-electron density together
subroutine simple_saha(ncomps_consider,k,tempx,rhox,ne)

  use blmod, only: comp, comp_details, ncomps, zav, ion_fractions
  use parameters
  use physical_constants
  use eosmodule, only: xxip, stat_weight_p1_ratio
  implicit none

  integer :: k
  integer :: ncomps_consider
  real*8  :: tempx, rhox

  real*8 :: ne

  integer, parameter :: maxstates_cap = 30
  integer, parameter :: itmax = 100
  real*8,  parameter :: zavtol = 1.0d-10
  real*8,  parameter :: zavtolinv = 1.0d10

  real*8 :: y(maxstates_cap)
  real*8 :: f(maxstates_cap)

  real*8 :: zion
  real*8 :: coef, saha_prefac
  real*8 :: n_k, zav_star, delta_zav
  real*8 :: function_x, function_deriv
  real*8 :: ktinv
  real*8 :: sum0, sum1, sum2, sum3
  real*8 :: productx, denom
  real*8 :: zavstar_x_nk_inv

  integer :: i,j,ie,izion
  integer :: min_state
  integer :: max_state
  integer :: number_states
  integer :: maxstates
  integer :: nsolve

!------------------------------------------------------------------------------
! based on the work of Zaghloul et al., J. Phys. D: Appl. Phys. 33:977 (2000)
! here the state number 1 is neutral and the total number of states is Z+1

  if(tempx.lt.0.0d0) then
     write(*,*) 'temperature is negative', k
     stop 1
  end if

  ne = 0.0d0

  nsolve = max(0, min(ncomps_consider, ncomps))
  maxstates = 2
  if(nsolve.gt.0) then
     maxstates = int(maxval(comp_details(1:nsolve,2))) + 1
  end if
  maxstates = max(2, min(maxstates_cap, maxstates))

  ktinv = 1.0d0/(tempx*kboltz)
  coef = saha_coeff * tempx**(1.5d0)
  saha_prefac = coef * exp(10.0d0*rhox)

  ! solve each selected element using its current mean-charge estimate
  do j=1, ncomps

     zion = comp_details(j,2)
     izion = int(zion)
     number_states = izion + 1
     if(number_states.gt.maxstates) then
        number_states = maxstates
     end if
     n_k = rhox*comp(k,j)/(comp_details(j,1)*mproton)

     if(j.le.nsolve.and.comp_details(j,2).gt.0.0d0) then

        zav_star = zav(j,k)

        min_state = 1
        max_state = number_states

        zavstar_x_nk_inv = 1.0d0/(zav_star*n_k)
        do i = 1, number_states - 1

            f(i) = stat_weight_p1_ratio(j,i) * saha_prefac &
                                * exp(max(-700.0d0,-xxip(izion,i)*ktinv))

            if( (f(i)*zavstar_x_nk_inv) .gt. zavtolinv ) then
                y(i) = 0.0d0
                min_state = i + 1
            end if

            if( (f(i)*zavstar_x_nk_inv) .lt. zavtol ) then
                max_state = i
                y(i+1:number_states) = 0.0d0
                exit
            end if
        end do

        if(max_state.eq.1) then
            y(1) = 1.0d0
            zav_star = 1.0d-6
        else if(min_state.eq.number_states) then
            y(number_states) = 1.0d0
            zav_star = zion
        else if(min_state.eq.max_state) then
            y(min_state) = 1.0d0
            zav_star = dble(min_state) - 1.0d0
        else

            delta_zav = zav_star
            ie = 0

            do while(abs(delta_zav/zav_star) > zavtol)

                ie = ie+1

                productx = 1.0d0

                sum0 = 0.0d0
                sum1 = 0.0d0
                sum2 = 0.0d0
                sum3 = 0.0d0

                do i = min_state, max_state - 1
                    productx = productx*f(i)/(zav_star*n_k)

                    sum0 = sum0 + productx
                    sum1 = sum1 + productx * i
                    sum2 = sum2 + productx * (i - dble(min_state) + 1.0d0)
                    sum3 = sum3 + productx * i * (i - dble(min_state) + 1.0d0)

                end do

                denom = 1.0d0/(dble(min_state) - 1.0d0 + sum1)

                function_x = 1.0d0 - zav_star * (1.0d0 + sum0) * denom

                function_deriv = (sum2 - (1+sum0) * (1+sum3*denom)) * denom

                delta_zav = -function_x/function_deriv
                zav_star = zav_star + delta_zav

                if(zav_star.gt.zion) then
                    zav_star = 0.5d0*zion
                end if

                if(ie.gt.ITMAX) then
                    write(6,*) "convergence problem in saha solver", k
                    exit
                endif

            end do

            y(min_state) = zav_star * denom
            do i=min_state+1, max_state
                y(i) = y(i-1)*f(i-1)/(zav_star*n_k)
            end do

        end if
        do i = 1, number_states
            ion_fractions(j,i,k) = y(i)
        end do
    else
        y(1:number_states-1) = 0.0d0
        y(number_states) = 1.0d0
        zav_star = zion
        if(number_states.gt.1) then
            ion_fractions(j,1:number_states-1,k) = 0.0d0
        end if
        ion_fractions(j,number_states,k) = 1.0d0
    end if

    ne = ne + zav_star*n_k

    zav(j,k) = zav_star

  end do


end subroutine simple_saha
