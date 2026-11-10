subroutine matrix_arrays(temp_x, lambda_x, inv_kappa_x, eps_x, lum_x, &
     mechanical_work, t1_array_x, t2_array_x, t3_array_x, t4_array_x)

  use blmod, only: delta_mass, delta_cmass, dtime, theta, eps_p, lum, dedt, r, &
       bomb_heating, Ni_heating
  use parameters, only: imax
  use physical_constants, only: pi, a_rad, clite
  implicit none
  real*8, intent(in) :: temp_x(imax), lambda_x(imax), inv_kappa_x(imax)
  real*8, intent(in) :: eps_x(imax), lum_x(imax), mechanical_work(imax)
  real*8, intent(out) :: t1_array_x(imax-1), t2_array_x(imax-1), t3_array_x(imax-1), t4_array_x(imax-1)
  real*8 :: const, left, right
  integer :: i

  ! the basis of the triadiagonal solver
  ! temperature correction is of a 3-way split -> composed of 4 components: term 1, term 2, term 3 and term 4
  ! t1_i delta T_{i+1} + t2_i delta T_{i} + t3_i delta T_{i-1} = t4_i

  ! be careful if refining how many cells a target cell has an effect on
  ! for example, instead of the target cell influencing the temperature of its two neighbour cells,
  ! the code could be changed to extent two or even three neighbour cells, so sharp discontinuities are problematic

  ! constant arises from differenating lum. with respect to temperature
  
  const=theta*dtime*64.0d0*pi*pi*a_rad*clite/3.0d0
  do i=1,imax-1
     right=const/delta_mass(i)*r(i+1)**4*lambda_x(i+1)*inv_kappa_x(i+1)/delta_cmass(i)
     left=0.0d0
     t3_array_x(i)=0.0d0
     if (i.gt.1) then
        left=const/delta_mass(i)*r(i)**4*lambda_x(i)*inv_kappa_x(i)/delta_cmass(i-1)
        t3_array_x(i)=-left*temp_x(i-1)**3
     endif
     t2_array_x(i)=dedt(i)+(left+right)*temp_x(i)**3
     t1_array_x(i)=0.0d0
     if (i.lt.imax-1) t1_array_x(i)=-right*temp_x(i+1)**3

     ! outer boundary: t4 terms consitute 'free' terms; not affected by temperature
     t4_array_x(i)=eps_p(i)+mechanical_work(i)-eps_x(i) &
          -dtime/delta_mass(i)*((1.0d0-theta)*(lum(i+1)-lum(i)) &
          +theta*(lum_x(i+1)-lum_x(i))) &
          +dtime*(bomb_heating(i)+Ni_heating(i))
  enddo

end subroutine matrix_arrays
