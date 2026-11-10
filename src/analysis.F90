subroutine analysis

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i

  real*8 :: vel_delta, lum_safe, lum_for_mag, rad_for_mag, analysis_time
  real*8 :: lum_radiative_observed
  real*8 :: bol_corr_used(11), band_rel_instant(11)
  real*8, parameter :: T_eff_min = 5000.0d0

  ! Fractional zone contribution for smooth Ni heating sum
  real*8 :: frac_ni_zone, dm_span, ni_zone_step, step_frac

!------------------------------------------------------------------------------
!---------- Calculating optical depth and tracing the photosphere -------------

  ! analysis runs inside blstep after the hydrodynamic/radiation update but
  ! before the main loop advances time, so its state belongs to time+dtime.
  analysis_time = time + dtime
  ! kappa_table (without the opacity floor) is used to trace the photosphere
  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))
  do i=1,imax
     if(tau(i).ne.tau(i) .or. tau(i).lt.0.0d0) tau(i) = 0.0d0
  enddo

  ! find the grid point, where the photosphere is located
  do i=imax-1, 1, -1
     if(tau(i).gt.0.66d0) then
        index_photo = i + 1
        exit
     endif
  enddo

  ! fix the moment, when the photosphere reaches the inner boundary, if it does
  if(tau(1).lt.0.66d0) then
     index_photo = 1
     if(photosphere_fell_on_the_center.eq.0) then
        photosphere_fell_on_the_center = 1
        open(unit=666, &
             file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
             status="unknown",form='formatted',position="append")
        write(666,*) 'Photosphere reached the center at ', analysis_time, 'seconds'
        close(666)
     end if
  end if

  !find the values of some variables at the photosphere
  if(photosphere_fell_on_the_center.eq.0) then
     if(index_photo.ge.imax) then
        lum_photo  = lum(imax)
        mass_photo = mass(imax)
        vel_photo  = vel(imax)
        rad_photo  = r(imax)
     else if(index_photo.le.1) then
        lum_photo  = lum(1)
        mass_photo = mass(1)
        vel_photo  = vel(1)
        rad_photo  = r(1)
     else
        call map_map(lum_photo,  0.66d0, lum(imax:1:-1),  tau(imax:1:-1),imax)
        call map_map(mass_photo, 0.66d0, mass(imax:1:-1), tau(imax:1:-1),imax)
        call map_map(vel_photo,  0.66d0, vel(imax:1:-1),  tau(imax:1:-1),imax)
        call map_map(rad_photo,  0.66d0, r(imax:1:-1),    tau(imax:1:-1),imax)
     endif
  else
     lum_photo   = lum(1)
     mass_photo  = mass(1)
     vel_photo   = vel(1)
     rad_photo   = r(1)
  end if

  ! calculates valid photometry (bounded to be physical), and calculates temperature floor
  photometry_valid = 1 - min(1,max(0,photosphere_fell_on_the_center))
  temperature_floor_fraction = dble(count(temp(1:imax-1).le. &
       evolution_temperature_floor_k*(1.0d0+1.0d-10)))/dble(max(1,imax-1))

  if(index_photo.lt.1) index_photo = 1
  if(index_photo.gt.imax) index_photo = imax
  if(lum_photo.ne.lum_photo .or. lum_photo.lt.0.0d0) then
     if(lum(index_photo).eq.lum(index_photo) .and. lum(index_photo).gt.0.0d0) then
        lum_photo = lum(index_photo)
     else
        lum_photo = 0.0d0
     endif
  endif
  if(mass_photo.ne.mass_photo) then
     if(mass(index_photo).eq.mass(index_photo)) then
        mass_photo = mass(index_photo)
     else
        mass_photo = 0.0d0
     endif
  endif
  if(vel_photo.ne.vel_photo) then
     if(vel(index_photo).eq.vel(index_photo)) then
        vel_photo = vel(index_photo)
     else
        vel_photo = 0.0d0
     endif
  endif
  if(rad_photo.ne.rad_photo .or. rad_photo.le.0.0d0) then
     if(r(index_photo).eq.r(index_photo) .and. r(index_photo).gt.0.0d0) then
        rad_photo = r(index_photo)
     else
        rad_photo = 1.0d0
     endif
  endif
  if(rad_photo.lt.r(1)) rad_photo = r(1)
  if(rad_photo.gt.r(imax)) rad_photo = r(imax)
  if(mass_photo.lt.mass(1)) mass_photo = mass(1)
  if(mass_photo.gt.mass(imax)) mass_photo = mass(imax)

  ! Local LTE temperature at the same tau=2/3 surface.  This provides an
  ! independent photospheric diagnostic while the selected grey luminosity may
  ! temporarily hand over between the photosphere and the outer FLD boundary
  if(index_photo.ge.imax) then
     T_photo_local=max(temp(max(1,imax-1)),0.0d0)
  else if(index_photo.le.1) then
     T_photo_local=max(temp(1),0.0d0)
  else
     call map_map(T_photo_local,0.66d0,temp(imax:1:-1),tau(imax:1:-1),imax)
     T_photo_local=max(T_photo_local,0.0d0)
  endif
  if(T_photo_local.ne.T_photo_local) T_photo_local=0.0d0
  lum_radiative_observed=lum_photo

  !photosphere tracer is equal to 1 at the photosphere, and 0 everywhere else
  !used for visualization
  photosphere_tracer(:) = 0.0d0
  photosphere_tracer(index_photo) = 1.0d0

  ! check if the photosphere moves through the regions with wrong opacity
  ! log10(T) = 3.75 is the lower boundary of the OPAL opacity tables
  if(metallicity(index_photo).gt.envelope_metallicity .and. &
       log10(temp(index_photo)).lt.3.75d0 .and. &
       kappa_table(index_photo).gt.kappa(index_photo) &
       .and. index_photo.gt.1 ) then
     opacity_corrupted = 1
  else
     opacity_corrupted = 0
  end if

!------------------------ Tracing the luminosity shell ------------------------
  index_lumshell = imax
  do i=1, imax - 1
     if(vel(imax-i).gt.0.0d0 .and. vel(imax-i).eq.vel(imax-i) .and. &
          tau(imax-i).gt.(clite/vel(imax-i))) then
        index_lumshell = imax - i + 1
        exit
     end if
  end do
  if (vel(2).gt.0.0d0 .and. vel(2).eq.vel(2) .and. tau(2).le.(clite/vel(2))) then
     index_lumshell = 1
  end if
  mass_lumshell = mass(index_lumshell)

  ! characteristic diffusion and expansion times for different shells
  do i=1, imax-1
     time_diff(i) = kappa(i)*rho(i)*(r(imax)-r(i))**2.0/clite
     if(time_diff(i).ne.time_diff(i) .or. time_diff(i).lt.0.0d0) time_diff(i) = 0.0d0
     vel_delta = vel(imax)-vel(i)
     if(abs(vel_delta).gt.1.0d-99 .and. vel_delta.eq.vel_delta) then
        time_exp(i) = (r(imax)-r(i))/vel_delta
     else
        time_exp(i) = 0.0d0
     endif
     if(time_exp(i).ne.time_exp(i)) time_exp(i) = 0.0d0
  end do

  ! internal energies of shells from the given radius out to the surface
  do i=1, imax
     E_shell(i) = sum(eps(i:imax)*delta_mass(i:imax))
  end do

!------------------- Calculate the observed luminosity ------------------------
  if(photosphere_fell_on_the_center.eq.0) then
     lum_observed = lum_radiative_observed + sum(Ni_energy_rate* &
          Ni_deposit_function(index_photo:imax)*delta_mass(index_photo:imax))

     if (smooth_ni_luminosity .ge. 1 .and. index_photo .gt. 1) then
        ni_zone_step = Ni_energy_rate * &
             Ni_deposit_function(index_photo-1) * delta_mass(index_photo-1) ! ni smoothing
        step_frac = ni_zone_step / max(lum_observed, 1.0d30)
        if (step_frac .gt. 0.02d0) then
           dm_span = mass(index_photo) - mass(index_photo-1)
           if (dm_span .gt. 1.0d-20) then
              frac_ni_zone = (mass(index_photo) - mass_photo) / dm_span
              frac_ni_zone = max(0.0d0, min(1.0d0, frac_ni_zone))
           else
              frac_ni_zone = 0.0d0
           endif
           lum_observed = lum_observed + frac_ni_zone * ni_zone_step
        endif
     endif
  else
     lum_observed = lum(1) + &
          sum(Ni_energy_rate*Ni_deposit_function(1:imax)*delta_mass(1:imax)) ! calculating total contribution
  end if
  if(lum_observed.ne.lum_observed .or. lum_observed.lt.0.0d0) lum_observed = max(lum_photo, 0.0d0)

  ! grey emergent luminosity: use the radiative flux that has actually
  ! propagated through the hydrorad solve and leaves the outer
  ! boundary. In multiple-shell CSM this captures
  ! radiation generated exterior to the legacy tau=2/3 photosphere

  lum_emergent = max(lum(imax), 0.0d0)
  if(lum_emergent.ne.lum_emergent) lum_emergent = 0.0d0
  
  ! publish outer-boundary flux. a regular 2/3 tau photosphere estimator can be near-zero
  ! even when a physical flux exists. An extended CSM envelope can bring about a delayed numerical feature
  ! that the usual observable (lum. above photosphere + gamma rays) cannot capture for such thin media. As such, we track boundary flux 
  ! (that actually takes into consideration the emerging FLD radiation)
  ! when the two channels cross, we revert back to the usual SNEC observable
  ! the two channels meet at the same point. we essentially swap observable curves
  ! to account for the effect mentioned above, then switch back when appropriate

  ! some artifical components of luminosity can be created when switching channels. here,
  ! we select which observable is higher rather than adding them together, since that would be redundant
  if (trim(adjustl(csm)).eq.'active' .and. &
       trim(adjustl(mass_scheme)).eq.'constant_density') then
     lum_csm_filtered=max(lum_observed,lum_emergent)
  else if (trim(adjustl(csm)).eq.'active' .and. &
       (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss')) then
     ! if using multiple stage mode, this flag controls the luminosity estimator in this evolution loop   
     if (csm_luminosity_handover_done.eq.0 .and. analysis_time.gt.0.25d0*86400.0d0 .and. &
          lum_emergent.gt.0.0d0 .and. lum_observed.ge.lum_emergent) then
        csm_luminosity_handover_done = 1
        csm_luminosity_handover_time = analysis_time
     endif
     if (csm_luminosity_handover_done.eq.1) then
        lum_csm_filtered = lum_observed
     else
        lum_csm_filtered = lum_emergent
     endif
  else
     lum_csm_filtered = lum_observed
  endif

  ! bound to be above zero
  lum_public = max(lum_csm_filtered,0.0d0)

  ! write down the time when the contribution of the Ni above the
  ! photosphere to the luminosity is greater than 5%
  if(shockpos_stop.eq.1 .and. lum_photo.gt.0.0d0 .and. &
       abs(lum_observed-lum_photo).gt.0.05d0*max(abs(lum_photo),1.0d-99) .and. &
       Ni_contributes_five_percents.eq.0) then

     Ni_contributes_five_percents = 1
     open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
          status="unknown",form='formatted',position="append")
     write(666,*) 'Ni contribution to the luminosity is 5% at ', &
          analysis_time, 'seconds'
     close(666)

  end if

!-------- Calculation of color magnitudes using bolometric corrections --------

  ! Broadband magnitudes use the selected grey luminosity and the same
  ! inherited colour-temperature and bolometric-correction tables.
  rad_for_mag=max(rad_photo,1.0d0)
  lum_for_mag=lum_public
  if (lum_for_mag.gt.0.0d0) then
     T_eff=(lum_for_mag/(4.0d0*pi*sigma_SB*rad_for_mag**2))**0.25d0
  else
     T_eff=0.0d0
  endif
  if(T_eff.ne.T_eff .or. T_eff.lt.0.0d0) T_eff=0.0d0
  T_eff_raw=T_eff

  ! adhere to 5000K lower limit
  T_color = max(T_eff,T_eff_min)
  T_color = min(T_color,temp_bol_corr(nlines_bol_corr))

  lum_safe = max(lum_for_mag/sun_lum, 1.0d-99)
  do i=1, 11
     call map_map(bol_corr_used(i),T_color,bol_corr(1:nlines_bol_corr,i), &
          temp_bol_corr,nlines_bol_corr)
     !  M_band = M_sun - 2.5 log10[(L/L_sun) 10^(0.4 BC)].
     band_rel_instant(i)=lum_safe*10.0d0**(0.4d0*bol_corr_used(i))
     magnitudes(i)=sun_mag-2.5d0*log10(max(band_rel_instant(i),1.0d-99))
  enddo

end subroutine analysis
