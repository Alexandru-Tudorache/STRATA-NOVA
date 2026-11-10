subroutine read_profile_compositions(prof_name)

  use blmod, only: comp, H_number, He_number, C_number, &
                    O_number, Ni_number, mass, ye, abar, metallicity, &
                    comp_details, ncomps, delta_mass
  use composition_profile_helpers, only: load_raw_composition_profile, &
       identify_major_species, normalize_raw_composition, &
       verify_mass_fraction_closure
  use composition_profile_store, only: store_processed_profile, map_stored_profile_to_grid
  use parameters
  use physical_constants
  use mixing_utils, only: build_effective_delta_mass
  use startup_log, only: print_startup_separator
  implicit none

  character(*) :: prof_name

  real*8,allocatable :: pmass(:),pradius(:),pcomp(:,:)

  real*8 :: Ni_mass_g, m_cut, m_ej, q_mix, delta_q, q, q0
  real*8 :: sum_fdm, scaling, Ni_uniform_frac
  real*8, allocatable :: fNi(:)
  real*8, allocatable :: profile_delta_mass(:)
  real*8, allocatable :: saved_csm_comp(:,:), csm_grid_mix(:)
  real*8 :: pmass_inner, pmass_outer
  real*8 :: ejecta_inner_mass, ejecta_outer_mass
  real*8 :: mix_outer_mass, mix_outer_mass2
  real*8 :: q_mix2, f_outer
  real*8 :: Ni_core_mass, Ni_outer_mass
  real*8 :: support_core_mass, support_outer_mass
  real*8 :: core_fraction, outer_fraction
  real*8 :: diag_ni_total, diag_ni_above, diag_ni_below
  real*8 :: diag_ni_grid, diag_ejecta_grid
  real*8 :: smooth_nominal_mass_msun
  real*8 :: smooth_ni_nominal_mass_msun
  integer :: smooth_number_iterations
  integer :: smooth_ni_number_iterations
  logical :: ni_override_active, preserve_ni_during_general
  integer :: i_excise, smooth_zones, ejecta_last_zone, csm_first_zone, csm_zones

  real*8, allocatable :: saved_ni(:)
  real*8 :: sum_other, target_frac, scale_frac, csm_surface_mass
  integer :: j

  integer :: profile_zones
  integer :: i,l

  character(len=1024) :: filename


  call load_raw_composition_profile(prof_name, profile_zones, pmass, pradius, pcomp)
  write(*,'(A,I0)') ' Composition zones: ',profile_zones

  call identify_major_species(H_number, He_number, C_number, O_number, Ni_number)

  ejecta_last_zone = profile_zones
  csm_first_zone = profile_zones + 1
  csm_zones = 0
  csm_surface_mass = 0.0d0
  if (trim(adjustl(csm)).ne.'none') then
     csm_zones = csm_extension_zones
     ejecta_last_zone = profile_zones - csm_zones
     if (csm_zones.le.0 .or. ejecta_last_zone.lt.1 .or. &
          ejecta_last_zone.ge.profile_zones) then
        write(*,*) 'ERROR: invalid ejecta/CSM zone boundary in generated composition profile.'
        stop 1
     endif
     csm_first_zone = ejecta_last_zone + 1
     csm_surface_mass = pmass(ejecta_last_zone)
     allocate(saved_csm_comp(csm_zones,ncomps), csm_grid_mix(ncomps))
     saved_csm_comp(:,:) = pcomp(csm_first_zone:profile_zones,1:ncomps)
     if (Ni_number.ne.0) saved_csm_comp(:,Ni_number) = 0.0d0
     do i = 1, csm_zones
        sum_other = sum(saved_csm_comp(i,1:ncomps))
        if (sum_other.le.0.0d0) then
           write(*,*) 'ERROR: generated CSM composition has zero total abundance.'
           stop 1
        endif
        saved_csm_comp(i,1:ncomps) = saved_csm_comp(i,1:ncomps)/sum_other
     enddo
     csm_grid_mix(:) = saved_csm_comp(1,1:ncomps)
     pcomp(csm_first_zone:profile_zones,1:ncomps) = saved_csm_comp(:,:)
  endif

  if(Ni_switch.eq.1 .and. Ni_number.eq.0) then
      write(*,*) 'radioactive Ni is absent in the composition profile'
      write(*,*) 'please, add a column for it (see documentation of the code)'
      write(*,*) 'or put Ni_switch = 0'
      stop 1
  endif


  ! measure the ejecta mass range before placing any manual ni
  pmass_inner = minval(pmass(1:ejecta_last_zone))
  pmass_outer = maxval(pmass(1:ejecta_last_zone))
  ejecta_inner_mass = max(pmass_inner, mass(1))
  ejecta_outer_mass = max(pmass_outer, ejecta_inner_mass)

  if (Ni_by_hand.ne.0 .and. Ni_number.ne.0) then
     do i = 1, ejecta_last_zone
        pcomp(i, Ni_number) = 0.d0
     enddo
  endif

  if (Ni_by_hand.ne.0 .and. Ni_number.ne.0 .and. Ni_mass.gt.0.d0) then

     Ni_mass_g = Ni_mass*msun
     m_cut = ejecta_inner_mass
     m_ej  = ejecta_outer_mass - m_cut
     if (m_ej.le.0.d0) then
        write(*,*) 'ERROR in read_profile_compositions: m_ej <= 0'
        stop 1
     endif

     q_mix = Ni_mix_fraction
     if (q_mix.le.0.d0) then
        write(*,*) 'WARNING: Ni_mix_fraction <= 0; skipping Ni_by_hand remapping.'
        goto 100
     endif
     if (q_mix.gt.1.d0) q_mix = 1.0d0

     mix_outer_mass = m_cut + q_mix * m_ej

     if (Ni_mix_kernel.ne.1 .and. Ni_mix_kernel.ne.2 .and. Ni_mix_kernel.ne.3) then
        write(*,*) 'WARNING: unsupported Ni_mix_kernel =', Ni_mix_kernel, &
                   '; using kernel=2 (exponential).'
        Ni_mix_kernel = 2
     endif

     call print_startup_separator

     ! kernel 1 places a uniform ni fraction inside the selected mass interval
     if (Ni_mix_kernel .eq. 1) then
        Ni_uniform_frac = Ni_mass_g / (q_mix * m_ej)
        if (Ni_uniform_frac .gt. 1.0d0) then
           q0 = Ni_mass_g / m_ej
           q0 = min(q0, 1.0d0)
           write(*,'(A,F8.4,A,F8.4,A,F8.4)') &
                '  WARNING: Ni_mass cannot fit in q_mix=', q_mix, &
                ' (would need X_Ni=', Ni_uniform_frac, &
                '); expanding to q_mix=', q0
           q_mix = q0
           mix_outer_mass = m_cut + q_mix * m_ej
           Ni_uniform_frac = Ni_mass_g / (q_mix * m_ej)
        endif
        write(*,*) 'Ni mixing (step): kernel=1, q_mix =', q_mix, &
                   ' mix_outer_mass =', mix_outer_mass/msun, ' Msun'
        write(*,'(A,F8.5)') '    X_Ni = ', Ni_uniform_frac
        do i = 1, ejecta_last_zone
           if (pmass(i) .ge. m_cut .and. pmass(i) .le. mix_outer_mass) then
              pcomp(i,Ni_number) = Ni_uniform_frac
           else
              pcomp(i,Ni_number) = 0.d0
           endif
        enddo

     ! kernel 2 uses a flat core followed by an exponential tail
     else if (Ni_mix_kernel .eq. 2) then
        allocate(profile_delta_mass(ejecta_last_zone))
        call build_effective_delta_mass(ejecta_last_zone, pmass(1:ejecta_last_zone), profile_delta_mass)
        allocate(fNi(ejecta_last_zone))

        write(*,*) 'Ni mixing (exponential): kernel=2, q_mix =', q_mix, &
                   ' mix_outer_mass =', mix_outer_mass/msun, ' Msun'

        q0 = 0.5d0*q_mix
        delta_q = (q_mix - q0)/3.d0
        if (delta_q.le.1.d-6) delta_q = 1.d-6

        sum_fdm = 0.d0
        do i = 1, ejecta_last_zone
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass) then
               fNi(i) = 0.d0
               cycle
           endif
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0

           if (q.le.q0) then
              fNi(i) = 1.d0
           else if (q.le.q_mix) then
              fNi(i) = exp( -(q - q0)/delta_q )
           else
              fNi(i) = 0.d0
           endif

           sum_fdm = sum_fdm + fNi(i)*profile_delta_mass(i)
        enddo

        if (sum_fdm.le.0.d0) then
           write(*,*) 'WARNING: sum_fdm <= 0, skipping Ni_by_hand remapping.'
        else
           scaling = Ni_mass_g/sum_fdm
           do i = 1, ejecta_last_zone
              if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass) cycle
              pcomp(i,Ni_number) = scaling*fNi(i)
           enddo
        endif

        deallocate(fNi)
        deallocate(profile_delta_mass)

     ! kernel 3 divides ni between an inner core and an outer shell
     else if (Ni_mix_kernel .eq. 3) then
        allocate(profile_delta_mass(ejecta_last_zone))
        call build_effective_delta_mass(ejecta_last_zone, pmass(1:ejecta_last_zone), profile_delta_mass)

        q_mix2 = Ni_mix_component2_extent
        if (q_mix2.lt.0.d0) q_mix2 = q_mix
        if (q_mix2.gt.1.d0) q_mix2 = 1.0d0
        f_outer = Ni_mix_component2_fraction
        if (f_outer.lt.0.d0) f_outer = 0.d0
        if (f_outer.gt.1.d0) f_outer = 1.d0

        if (q_mix2.lt.q_mix) then
           write(*,*) 'WARNING: kernel=3 requires component2 extent >= Ni_mix_fraction.'
           write(*,*) '         Clamping Ni_mix_component2_extent to', q_mix
           q_mix2 = q_mix
        endif

        q0 = (1.d0 - f_outer) * Ni_mass_g / m_ej
        if (q0 .gt. q_mix) then
           write(*,'(A,F8.4,A,F8.4)') &
                '  WARNING: core Ni cannot fit in q_mix=', q_mix, &
                '; expanding to q_mix=', min(q0, 1.0d0)
           q_mix = min(q0, 1.0d0)
           mix_outer_mass = m_cut + q_mix * m_ej
        endif
        if (f_outer .gt. 0.d0) then
           q0 = q_mix + f_outer * Ni_mass_g / m_ej
           if (q0 .gt. q_mix2) then
              write(*,'(A,F8.4,A,F8.4)') &
                   '  WARNING: outer Ni cannot fit in q_mix2=', q_mix2, &
                   '; expanding to q_mix2=', min(q0, 1.0d0)
              q_mix2 = min(q0, 1.0d0)
           endif
        endif

        mix_outer_mass2 = m_cut + q_mix2 * m_ej
        Ni_outer_mass = f_outer * Ni_mass_g
        Ni_core_mass = Ni_mass_g - Ni_outer_mass
        if (Ni_core_mass.lt.0.d0) Ni_core_mass = 0.d0

        support_core_mass = 0.d0
        support_outer_mass = 0.d0
        do i = 1, ejecta_last_zone
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass2) cycle
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0
           if (q.le.q_mix) then
              support_core_mass = support_core_mass + profile_delta_mass(i)
           else if (q.le.q_mix2) then
              support_outer_mass = support_outer_mass + profile_delta_mass(i)
           endif
        enddo

        if (Ni_core_mass.gt.0.d0 .and. support_core_mass.le.0.d0) then
           write(*,*) 'ERROR: kernel=3 core support mass is zero.'
           write(*,*) '       Increase Ni_mix_fraction.'
           stop 1
        endif

        if (Ni_outer_mass.gt.0.d0 .and. support_outer_mass.le.0.d0) then
           write(*,*) 'WARNING: kernel=3 outer shell has zero support mass.'
           write(*,*) '         Moving component-2 Ni mass back into core component.'
           Ni_core_mass = Ni_core_mass + Ni_outer_mass
           Ni_outer_mass = 0.d0
        endif

        core_fraction = 0.d0
        outer_fraction = 0.d0
        if (Ni_core_mass.gt.0.d0) core_fraction = Ni_core_mass / support_core_mass
        if (Ni_outer_mass.gt.0.d0) outer_fraction = Ni_outer_mass / support_outer_mass

        write(*,*) 'Ni mixing (two-component): kernel=3, q_core =', q_mix, &
                   ' q_outer =', q_mix2, ' f_outer =', f_outer
        write(*,*) '    mix_outer_mass(core)  =', mix_outer_mass/msun, ' Msun'
        write(*,*) '    mix_outer_mass(shell) =', mix_outer_mass2/msun, ' Msun'
        write(*,*) '    M_core =', Ni_core_mass/msun, ' Msun, M_outer =', &
                   Ni_outer_mass/msun, ' Msun'

        do i = 1, ejecta_last_zone
           if (pmass(i) .lt. m_cut .or. pmass(i) .gt. mix_outer_mass2) then
              pcomp(i,Ni_number) = 0.d0
              cycle
           endif
           q = (pmass(i) - m_cut)/m_ej
           if (q.lt.0.d0) q = 0.d0
           if (q.gt.1.d0) q = 1.d0

           if (q.le.q_mix) then
              pcomp(i,Ni_number) = core_fraction
           else if (q.le.q_mix2) then
              pcomp(i,Ni_number) = outer_fraction
           else
              pcomp(i,Ni_number) = 0.d0
           endif
        enddo

        if (maxval(pcomp(:,Ni_number)) .gt. 1.0d0) then
           write(*,'(A,F8.4)') '  WARNING: kernel=3 peak Ni fraction exceeds 1.0 (peak=', &
                maxval(pcomp(:,Ni_number))
        endif

        deallocate(profile_delta_mass)

     endif

  endif
100 continue

  call normalize_raw_composition(pcomp, profile_zones, ncomps, Ni_number)
  ! keep the attached csm composition separate from ejecta mixing
  if (csm_zones.gt.0) then
     pcomp(csm_first_zone:profile_zones,1:ncomps) = saved_csm_comp(:,:)
  endif

  if (Ni_number.ne.0) then
     allocate(profile_delta_mass(profile_zones))
     call build_effective_delta_mass(profile_zones, pmass, profile_delta_mass)
     diag_ni_total = 0.d0
     diag_ni_above = 0.d0
     diag_ni_below = 0.d0
     do i = 1, profile_zones
        diag_ni_total = diag_ni_total + pcomp(i,Ni_number)*profile_delta_mass(i)
        if (pmass(i) .ge. mass(1)) then
           diag_ni_above = diag_ni_above + pcomp(i,Ni_number)*profile_delta_mass(i)
        else
           diag_ni_below = diag_ni_below + pcomp(i,Ni_number)*profile_delta_mass(i)
        endif
     enddo
  endif

  ! smooth only the ejecta above the excision boundary
  if(boxcar_smoothing) then
      i_excise = 1
      do i = 1, ejecta_last_zone
         if (pmass(i) .ge. mass(1)) then
            i_excise = i
            exit
         endif
      enddo
      smooth_zones = ejecta_last_zone - i_excise + 1

      smooth_nominal_mass_msun = boxcar_nominal_mass_msun
      smooth_number_iterations = boxcar_number_iterations
      smooth_ni_nominal_mass_msun = smooth_nominal_mass_msun
      smooth_ni_number_iterations = smooth_number_iterations
      ni_override_active = .false.
      if (boxcar_smooth_ni .ne. 0 .and. Ni_number .ne. 0) then
         if (boxcar_ni_nominal_mass_msun.gt.0.d0) then
            smooth_ni_nominal_mass_msun = boxcar_ni_nominal_mass_msun
            ni_override_active = .true.
         endif
         if (boxcar_ni_number_iterations.gt.0) then
            smooth_ni_number_iterations = boxcar_ni_number_iterations
            ni_override_active = .true.
         endif
      endif
      preserve_ni_during_general = (boxcar_smooth_ni.eq.0 .or. ni_override_active) .and. &
           Ni_number.ne.0


      ! save ni before the general composition smoothing pass
      if (preserve_ni_during_general) then
         allocate(saved_ni(profile_zones))
         saved_ni(:) = pcomp(:, Ni_number)
      endif

      ! apply the mass-weighted boxcar to the processed profile
      call boxcar_profile(pmass(i_excise:ejecta_last_zone), &
                          pcomp(i_excise:ejecta_last_zone,:), &
                          smooth_zones, ncomps, &
                          smooth_nominal_mass_msun*msun, &
                          smooth_number_iterations)

      if (preserve_ni_during_general) then
         if (maxval(saved_ni(1:ejecta_last_zone)) .gt. 1.0d0) then
            write(*,'(A,F8.4)') '  WARNING: Ni fraction capped at 1.0 ' // &
                 '(peak was', maxval(saved_ni(1:ejecta_last_zone))
            write(*,'(A)')      '           Ni_mix_fraction is too small ' // &
                 'for given Ni_mass. Increase Ni_mix_fraction.'
         endif
         pcomp(1:ejecta_last_zone, Ni_number) = &
              min(saved_ni(1:ejecta_last_zone), 1.0d0)
         do i = i_excise, ejecta_last_zone
            sum_other = 0.0d0
            do j = 1, ncomps
               if (j /= Ni_number) sum_other = sum_other + pcomp(i,j)
            enddo
            target_frac = max(0.0d0, 1.0d0 - pcomp(i, Ni_number))
            if (sum_other > 0.0d0) then
               scale_frac = target_frac / sum_other
               do j = 1, ncomps
                  if (j /= Ni_number) pcomp(i,j) = pcomp(i,j) * scale_frac
               enddo
            endif
         enddo
      endif

      ! optionally smooth ni again with its own settings
      if (ni_override_active) then
         call boxcar_profile(pmass(i_excise:ejecta_last_zone), &
                             pcomp(i_excise:ejecta_last_zone,Ni_number:Ni_number), &
                             smooth_zones, 1, &
                             smooth_ni_nominal_mass_msun*msun, &
                             smooth_ni_number_iterations)
         if (maxval(pcomp(1:ejecta_last_zone,Ni_number)) .gt. 1.0d0) then
            write(*,'(A,F8.4)') '  WARNING: Ni fraction capped at 1.0 after Ni-only smoothing ' // &
                 '(peak was', maxval(pcomp(1:ejecta_last_zone,Ni_number))
         endif
         pcomp(1:ejecta_last_zone, Ni_number) = &
              min(pcomp(1:ejecta_last_zone, Ni_number), 1.0d0)
         do i = i_excise, ejecta_last_zone
            sum_other = 0.0d0
            do j = 1, ncomps
               if (j /= Ni_number) sum_other = sum_other + pcomp(i,j)
            enddo
            target_frac = max(0.0d0, 1.0d0 - pcomp(i, Ni_number))
            if (sum_other > 0.0d0) then
               scale_frac = target_frac / sum_other
               do j = 1, ncomps
                  if (j /= Ni_number) pcomp(i,j) = pcomp(i,j) * scale_frac
               enddo
            endif
         enddo
      endif

      if (preserve_ni_during_general) then
         deallocate(saved_ni)
      endif
  endif

  if (csm_zones.gt.0) then
     pcomp(csm_first_zone:profile_zones,1:ncomps) = saved_csm_comp(:,:)
  endif

  if (Ni_number.ne.0) then
     diag_ni_total = 0.d0
     diag_ni_above = 0.d0
     diag_ni_below = 0.d0
     do i = 1, profile_zones
        diag_ni_total = diag_ni_total + pcomp(i,Ni_number)*profile_delta_mass(i)
        if (pmass(i) .ge. mass(1)) then
           diag_ni_above = diag_ni_above + pcomp(i,Ni_number)*profile_delta_mass(i)
        else
           diag_ni_below = diag_ni_below + pcomp(i,Ni_number)*profile_delta_mass(i)
        endif
     enddo
     deallocate(profile_delta_mass)
  endif

  ! restart abundances describe intervals beginning at each mass coordinate
  if (initial_model_age_s.gt.0.0d0 .and. Ni_by_hand.eq.0) then
     call store_processed_profile(pmass, pcomp, profile_zones, ncomps, .true.)
  else
     call store_processed_profile(pmass, pcomp, profile_zones, ncomps)
  endif
  ! average the processed composition onto the hydro mass intervals
  call map_stored_profile_to_grid(mass, imax, ncomps, comp)

  if (csm_zones.gt.0) then
     do i = 1, imax
        if (mass(i).ge.csm_surface_mass) comp(i,1:ncomps) = csm_grid_mix(1:ncomps)
     enddo
  endif

  deallocate(pmass)
  deallocate(pradius)
  deallocate(pcomp)
  if (allocated(saved_csm_comp)) deallocate(saved_csm_comp)
  if (allocated(csm_grid_mix)) deallocate(csm_grid_mix)


  call verify_mass_fraction_closure(comp, imax, ncomps, 'composition mapping (custom)')

  if (Ni_number.ne.0) then
     diag_ni_grid = 0.d0
     diag_ejecta_grid = 0.d0
     do i = 1, imax
        diag_ni_grid = diag_ni_grid + comp(i,Ni_number)*delta_mass(i)
        diag_ejecta_grid = diag_ejecta_grid + delta_mass(i)
     enddo
     write(*,'(A,ES11.4,A,ES11.4,A)') ' Ni grid mass: ',diag_ni_grid/msun, &
          ' Msun   ejecta: ',diag_ejecta_grid/msun,' Msun'
     if (Ni_by_hand.ne.0 .and. Ni_mass.gt.0.0d0) then
        if (abs(diag_ni_grid/(Ni_mass*msun)-1.0d0).gt.5.0d-2) then
           write(*,'(A,F7.2,A)') ' WARNING: mapped Ni differs from Ni_mass by ', &
                100.0d0*(diag_ni_grid/(Ni_mass*msun)-1.0d0),' %'
        endif
     endif
  endif

  do i=1,imax
      ye(i) = sum(comp(i,1:ncomps)* &
          (comp_details(1:ncomps,2)/comp_details(1:ncomps,1)))
      abar(i) = 1.0d0/sum(comp(i,1:ncomps)/comp_details(1:ncomps,1))
  enddo


  if (output_mode .eq. 0) then
  do l=1,ncomps
  if(l.eq.H_number) then
      filename = trim(adjustl(outdir))//"/H_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.He_number) then
      filename = trim(adjustl(outdir))//"/He_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.C_number) then
      filename = trim(adjustl(outdir))//"/C_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.O_number) then
      filename = trim(adjustl(outdir))//"/O_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.Ni_number) then
      filename = trim(adjustl(outdir))//"/Ni_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  end if
  enddo


  filename = trim(adjustl(outdir))//"/comp_init_all.dat"
  open(unit=668,file=trim(adjustl(filename)),status='unknown',form='formatted')
  write(668,'(A,I4)') '# ncomps = ', ncomps
  do i=1,imax
     write(668,'(I5.4)',advance='no') i
     do l=1,ncomps
        write(668,'(E25.16E3)',advance='no') comp(i,l)
     enddo
     write(668,*)
  enddo
  close(668)
  endif

  if (output_mode .eq. 2) then
     filename = trim(adjustl(outdir))//"/comp_init_all.dat"
     open(unit=668,file=trim(adjustl(filename)),status='unknown',form='formatted')
     write(668,'(A,I4)') '# ncomps = ', ncomps
     do i=1,imax
        write(668,'(I5.4)',advance='no') i
        do l=1,ncomps
           write(668,'(E25.16E3)',advance='no') comp(i,l)
        enddo
        write(668,*)
     enddo
     close(668)
  endif

  do i=1,imax
      if(H_number.ne.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,H_number) - comp(i,He_number)
      else if(H_number.ne.0 .and. He_number.eq.0) then
          metallicity(i) = 1 - comp(i,H_number)
      else if(H_number.eq.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,He_number)
      else
          metallicity(i) = 1
      end if
  enddo

  if (output_mode .eq. 0) then
     filename = trim(adjustl(outdir))//"/metallicity_init.dat"
     call output_screenshot(metallicity,filename,imax)
  endif


end subroutine read_profile_compositions
