subroutine opacity(rho_x,temp_x,kappa_x,kappa_table_x,dkappadt_x)

  use blmod, only: logT, logR_op, opacity_floor
  use parameters
  use physical_constants
  implicit none

  real*8 rho_x(imax)
  real*8 temp_x(imax)

  real*8 kappa_x(imax)
  real*8 kappa_table_x(imax)
  real*8 dkappadt_x(imax)

  integer i
  real*8 t6
  real*8 R_op
  real*8 log_kappa
  real*8 temp_for_opacity
  real*8, parameter :: opacity_temp_min_k = 10.0d0**2.7d0

  integer no_tabular_value


  i = 1
  temp_for_opacity = max(temp_x(i), opacity_temp_min_k)
  t6 = temp_for_opacity*1.0d-6
  logT(i) = log10(temp_for_opacity)
  R_op = rho_x(i)/(t6**3)
  logR_op(i) = log10(R_op)
  ! interpolate the tabulated opacity in log temperature and log r
  call bicubic_interpolation(i,logR_op(i),logT(i),log_kappa,no_tabular_value)
  if (no_tabular_value.eq.0) then
     kappa_x(i) = 10.0d0**log_kappa
  else
     kappa_x(i) = opacity_floor(i)
  endif
  dkappadt_x(i) = 0.0d0

  do i=2, imax - 1
     temp_for_opacity = max(temp_x(i), opacity_temp_min_k)
     t6 = temp_for_opacity*1.0d-6
     logT(i) = log10(temp_for_opacity)
     R_op = rho_x(i)/(t6**3)
     logR_op(i) = log10(R_op)
     call bicubic_interpolation(i,logR_op(i),logT(i),log_kappa,no_tabular_value)
     if (no_tabular_value.eq.0) then
        kappa_x(i) = 10.0d0**log_kappa
     else
        kappa_x(i) = kappa_x(i-1)
     endif
     dkappadt_x(i) = 0.0d0
  enddo

  kappa_table_x(1:imax-1) = kappa_x(1:imax-1)

  do i=1, imax-1
    if( kappa_x(i) .lt. opacity_floor(i)) then
        kappa_x(i)   =  opacity_floor(i)
    end if
  end do

  kappa_x(imax) = kappa_x(imax-1)
  kappa_table_x(imax) = kappa_table_x(imax-1)
  dkappadt_x(imax) = dkappadt_x(imax-1)

end subroutine opacity


subroutine compose_opacity_tables_OPAL

  use blmod, only: logR_array, logT_array, opacity_tables, rpoints, tpoints, &
      op_rows, op_cols, envelope_metallicity, xxo, xxc, H_number, C_number, &
      O_number, comp
  use parameters
  use physical_constants
  use startup_log, only: print_startup_separator
  implicit none

  integer :: i,j,k
  real*8 :: R_op, t6

  real*4 opact,dopact,dopacr,dopactd
  common/e/ opact,dopact,dopacr,dopactd


  call print_startup_separator
  write(*,'(A)') ' Composing opacity tables from OPAL...'

	  tpoints = 134
	  rpoints = 19

	  if (allocated(opacity_tables)) deallocate(opacity_tables)
	  if (allocated(logT_array)) deallocate(logT_array)
	  if (allocated(logR_array)) deallocate(logR_array)
	  if (allocated(op_rows)) deallocate(op_rows)
	  if (allocated(op_cols)) deallocate(op_cols)

	  allocate(opacity_tables(imax,tpoints,rpoints))


  allocate(logT_array(tpoints))
  allocate(logR_array(rpoints))

  logT_array(1) = 2.7d0
  do i=2, 5
    logT_array(i) = logT_array(i-1)+0.05d0
  end do
  logT_array(6) = logT_array(5)+0.01d0
  do i=7, 11
    logT_array(i) = logT_array(i-1)+0.02d0
  end do
  do i=12, 60
    logT_array(i) = logT_array(i-1)+0.01d0
  end do
  do i=61, 110
    logT_array(i) = logT_array(i-1)+0.05d0
  end do
  do i=111, 131
    logT_array(i) = logT_array(i-1)+0.1d0
  end do
  do i=132,134
    logT_array(i) = logT_array(i-1)+0.2d0
  end do

  do j=1, rpoints
    logR_array(j) = -8.0d0 + (j-1)*0.5d0
  end do


  allocate(op_rows(tpoints,2))
  allocate(op_cols(rpoints,2))

  do i=1, tpoints
    op_rows(i,1) = 1
  end do
  do i=1, tpoints-13
    op_rows(i,2) = 19
  end do
  op_rows(tpoints-12,2) = 18
  op_rows(tpoints-11,2) = 17
  op_rows(tpoints-10,2) = 17
  op_rows(tpoints- 9,2) = 16
  op_rows(tpoints- 8,2) = 16
  op_rows(tpoints- 7,2) = 16
  op_rows(tpoints- 6,2) = 16
  op_rows(tpoints- 5,2) = 15
  op_rows(tpoints- 4,2) = 15
  op_rows(tpoints- 3,2) = 15
  op_rows(tpoints- 2,2) = 15
  op_rows(tpoints- 1,2) = 15
  op_rows(tpoints,2)    = 14

  do j=1, rpoints
    op_cols(j,1) = 1
  end do
  do j=1, 14
    op_cols(j,2) = 134
  end do
  op_cols(15,2) = 133
  op_cols(16,2) = 128
  op_cols(17,2) = 124
  op_cols(18,2) = 122
  op_cols(19,2) = 121

  do k=1, imax

    xxc(k) = MAX( 0.0d0, comp(k,C_number) - envelope_metallicity*C_frac_sol )

    xxo(k) = MAX( 0.0d0, comp(k,O_number) - envelope_metallicity*O_frac_sol )

    if (k.gt.1) then
      if (comp(k,H_number).eq.comp(k-1,H_number) .and. &
          xxc(k).eq.xxc(k-1) .and. xxo(k).eq.xxo(k-1)) then
        opacity_tables(k,:,:) = opacity_tables(k-1,:,:)
        cycle
      endif
    endif

    do i = 1, tpoints
        do j = 1, rpoints
            if (i.lt.op_cols(j,1) .or. i.gt.op_cols(j,2) &
                .or. j.lt.op_rows(i,1) .or. j.gt.op_rows(i,2)) then

                opacity_tables(k,i,j) = 0.0d0

            else

                if (i.eq.op_cols(j,1)) then
                    t6 = (10.0d0**(logT_array(i)+1.0d-5))*1.0d-6
                else if (i.eq.op_cols(j,2)) then
                    t6 = (10.0d0**(logT_array(i)-1.0d-5))*1.0d-6
                else
                    t6 = (10.0d0**logT_array(i))*1.0d-6
                end if
                if (j.eq.op_rows(i,1)) then
                    R_op = 10.0d0**(logR_array(j)+1.0d-5)
                else if (j.eq.op_rows(i,2)) then
                    R_op = 10.0d0**(logR_array(j)-1.0d-5)
                else
                    R_op = 10.0d0**logR_array(j)
                end if

                call opac( real(envelope_metallicity), real(comp(k,H_number)), &
                    real(xxc(k)), real(xxo(k)), real(t6), real(R_op) )

                opacity_tables(k,i,j) = dble(opact)

            end if
        end do
    end do

  end do

  write(*,*) "Finished composing opacity tables!"

end subroutine compose_opacity_tables_OPAL
