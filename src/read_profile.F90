subroutine get_inner_outer_mass_from_profile(prof_name,inner_mass,outer_mass)

  implicit none

  character(*)         :: prof_name

  real*8               :: inner_mass
  real*8               :: outer_mass

  integer              :: profile_zones
  integer              :: i,ibuffer
  real*8, allocatable  :: pmass(:)


  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')

  read(666,*) profile_zones
  allocate(pmass(profile_zones))

  do i=1,profile_zones
     read(666,*) ibuffer,pmass(i)
  enddo
  close(666)

  outer_mass = pmass(profile_zones)
  inner_mass = pmass(1)

  deallocate(pmass)

end subroutine get_inner_outer_mass_from_profile


subroutine get_inner_outer_radius_from_profile(prof_name,inner_mass, &
        inner_radius,outer_radius)

  implicit none

  character(*)         :: prof_name
  real*8               :: inner_mass

  real*8               :: inner_radius
  real*8               :: outer_radius

  integer              :: profile_zones
  integer              :: i,ibuffer
  real*8, allocatable  :: pradius(:)
  real*8, allocatable  :: pmass(:)


  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')

  read(666,*) profile_zones
  allocate(pradius(profile_zones))
  allocate(pmass(profile_zones))

  do i=1,profile_zones
     read(666,*) ibuffer, pmass(i), pradius(i)
  enddo
  close(666)

  call map_map(inner_radius,inner_mass,pradius,pmass,profile_zones)

  outer_radius = pradius(profile_zones)

  deallocate(pradius)
  deallocate(pmass)

end subroutine get_inner_outer_radius_from_profile


subroutine get_ncomps_from_profile(prof_name,xncomps)

  implicit none

  character(*)         :: prof_name

  integer              :: xncomps

  integer              :: profile_zones


  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')
  read(666,*) profile_zones,xncomps
  close(666)

end subroutine get_ncomps_from_profile


subroutine read_profile(prof_name)

  use blmod, only: mass, cmass, vel, rho, temp, ncomps, ye, abar, comp_details,&
                    eps, p, cs2, dedt, dpdt, entropy, zav, p_rad
  use parameters
  use physical_constants
  use eosmodule, only: init_ionpot
  use csm_switch, only: csm_generated_active
  implicit none

  character(*) :: prof_name
  integer :: profile_zones
  integer :: i,l
  integer :: ibuffer
  integer :: keytemp, keyerr
  logical :: preserve_csm_geometry

  real*8,allocatable :: pmass(:), pradius(:), ptemp(:), prho(:), pvel(:)


  open(666,file=trim(prof_name),status='unknown',form='formatted',action='read')

  read(666,*) profile_zones
  write(*,'(A,I0)') ' Profile zones: ',profile_zones

  allocate(pmass(profile_zones))
  allocate(pradius(profile_zones))
  allocate(ptemp(profile_zones))
  allocate(prho(profile_zones))
  allocate(pvel(profile_zones))

  do i=1,profile_zones
     read(666,*) ibuffer, pmass(i), pradius(i), ptemp(i), prho(i), pvel(i)
  enddo

  preserve_csm_geometry = csm_generated_active

  if (preserve_csm_geometry) then
     call redistribute_wind_outer_grid(profile_zones,pmass)
  endif

  do i=1,imax
      call map_map(vel(i), mass(i),pvel, pmass,profile_zones)
  enddo

  do i=1,imax-1
     call map_map(temp(i),cmass(i),ptemp,pmass,profile_zones)
  enddo

  if (preserve_csm_geometry) then
     call map_wind_geometry(profile_zones,pmass,pradius,pvel)
  else
     do i=1,imax-1
        call map_map(rho(i),cmass(i),prho,pmass,profile_zones)
     enddo
  endif

  rho(imax) = 0.0d0
  temp(imax) = 0.0d0


  deallocate(pmass)
  deallocate(pradius)
  deallocate(ptemp)
  deallocate(prho)
  deallocate(pvel)

  if(ncomps.gt.0) then
     call read_profile_compositions(composition_profile_name)
  endif

  if(eoskey.eq.2) then
     call init_ionpot
  endif


  !initialize zav (for the Saha solver, Paczynski EOS), assuming full ionization
  do l=1, ncomps
    do i=1, imax
        zav(l,i) = comp_details(l,2)
    end do
  end do

  keytemp = 1
  call eos(rho(1:imax-1),temp(1:imax-1),ye(1:imax-1), &
         abar(1:imax-1),p(1:imax-1),eps(1:imax-1), &
         cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), &
         entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)

  eps(imax) = 0.0d0
  p(imax) = 0.0d0

end subroutine read_profile


subroutine map_linterp(x1,x2,y1,y2,x,y)

  implicit none

  real*8 :: slope,x1,x2,y1,y2,x,y

    if (x2.lt.x1) then
       error stop "Error in linterp!"
    endif

    if (x2.ne.x1) then
        slope = (y2 - y1) / (x2 - x1)
    else
        slope = 0
    endif

    y = slope*(x-x1) + y1

end subroutine  map_linterp

subroutine map_find_index(zones,array,goal,upper_index,lower_index)

  implicit none

  integer :: zones
  real*8 :: array(*)
  real*8 :: goal
  integer :: middle_index,upper_index,lower_index

    if(goal.lt.array(1) .or. goal.gt.array(zones)) then
        write(*,*) 'value passed to map_find_index is out of the array'
        stop 1
    end if

    lower_index = 1
    upper_index = zones

    do while ( (upper_index - lower_index) .gt. 1 )
       middle_index = (lower_index + upper_index) * 0.5d0
       if ( goal .le. array(middle_index) ) then
            upper_index = middle_index
       else
            lower_index = middle_index
       end if
    enddo

end subroutine map_find_index


  ! locate a coordinate bracket before mapping the tabulated quantity
subroutine map_map(point_value,point_radius,parray,pradius,zones)

  implicit none

  real*8 :: point_value, point_radius
  real*8 :: pradius(*), parray(*)
  integer :: zones
  integer :: upper_index, lower_index

    if (point_radius .ge. pradius(1) .and. &
       point_radius .le. pradius(zones) )  then

     call map_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)

     call map_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  &
          point_radius, point_value )

    else if (point_radius .lt. pradius(1)) then
     call map_linterp(pradius(1),pradius(2), &
          parray(1),parray(2),point_radius,point_value)

    else if (point_radius .gt. pradius(zones)) then
     call map_linterp(pradius(zones-1),pradius(zones), &
          parray(zones-1),parray(zones),point_radius,point_value)
    endif


end subroutine map_map

  ! reserve hydro mass intervals for the attached csm and its transition
subroutine redistribute_wind_outer_grid(profile_zones,pmass)

  use blmod, only: mass, cmass, delta_mass, delta_cmass
  use parameters, only: imax, csm_blend_zones, mass_scheme, csm, csm_mass_msun
  use csm_switch, only: csm_generated_surface_mass_g, csm_generated_density_transition_zones
  use legacy_grid_pattern, only: build_legacy_mass_grid
  implicit none

  integer, intent(in) :: profile_zones
  real*8, intent(in) :: pmass(profile_zones)
  integer :: i, j, p_surface, p_bridge_end
  integer :: i_surface, n_csm, n_bridge, n_outer
  real*8 :: mstart, mend, msurf, mbridge, frac
  logical :: attached_mass_loss

  mstart = mass(1)
  mend = pmass(profile_zones)
  msurf = csm_generated_surface_mass_g
  attached_mass_loss = (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_density') .and. &
       trim(adjustl(csm)).eq.'active'

  if (msurf.le.mstart .or. msurf.ge.mend) then
     write(*,*) '******* Invalid generated CSM surface mass for wind-grid allocation.'
     write(*,*) '******* mstart, msurf, mend = ',mstart,msurf,mend
     stop 1
  endif

  n_csm = max(csm_blend_zones+8, max(24, imax/8))
  n_csm = min(64,n_csm)
  if (trim(adjustl(mass_scheme)).eq.'continuous_mass') then
     n_csm=max(csm_generated_density_transition_zones+8,max(24,imax/8))
  endif
  if (trim(adjustl(mass_scheme)).eq.'constant_density') then
     if (csm_mass_msun.le.1.0d-5) then
        n_csm=16
     else if (csm_mass_msun.le.3.0d-5) then
        n_csm=20
     else
        n_csm=24
     endif
  endif
  n_csm = min(imax-25,n_csm)
  if (trim(adjustl(mass_scheme)).ne.'constant_density') then
     if (n_csm.lt.csm_blend_zones+2) then
        write(*,*) '******* imax is too small to resolve the requested CSM bridge.'
        write(*,*) '******* imax, blend_zones = ',imax,csm_blend_zones
        stop 1
     endif
  endif
  i_surface = imax-n_csm

  ! construct the fixed native mass grid
  call build_legacy_mass_grid(mass(1:i_surface),i_surface,mstart,msurf)
  mass(i_surface)=msurf

  p_surface=1
  do i=1,profile_zones
     if (pmass(i).le.msurf*(1.0d0+1.0d-14)) p_surface=i
  enddo
  p_bridge_end=min(profile_zones,max(p_surface+1, &
       p_surface+csm_generated_density_transition_zones))
  mbridge=max(msurf,min(pmass(p_bridge_end),mend))

  n_bridge=min(csm_blend_zones,n_csm-2)
  if (trim(adjustl(mass_scheme)).eq.'continuous_mass') &
       n_bridge=min(csm_generated_density_transition_zones,n_csm-2)
  if (trim(adjustl(mass_scheme)).eq.'constant_density') then
     n_bridge=0
  else if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
     n_bridge=0
  else if (mbridge.le.msurf*(1.0d0+1.0d-14) .or. &
           mbridge.ge.mend*(1.0d0-1.0d-14)) then
     n_bridge=0
  else if (attached_mass_loss) then
     n_bridge=min(2,n_csm-2)
  endif

  if (n_bridge.gt.0) then
     do j=1,n_bridge
        frac=dble(j)/dble(n_bridge)
        mass(i_surface+j)=msurf+frac*(mbridge-msurf)
     enddo
     n_outer=n_csm-n_bridge
     do j=1,n_outer
        frac=dble(j)/dble(n_outer)
        mass(i_surface+n_bridge+j)=mbridge+frac*(mend-mbridge)
     enddo
  else
     do j=1,n_csm
        frac=dble(j)/dble(n_csm)
        mass(i_surface+j)=msurf+frac*(mend-msurf)
     enddo
  endif
  mass(imax)=mend

  do i=1,imax-1
     if (mass(i+1).le.mass(i)) then
        write(*,*) '******* Non-monotone wind mass grid after CSM allocation.'
        write(*,*) '******* zone, m_i, m_ip1 = ',i,mass(i),mass(i+1)
        stop 1
     endif
     cmass(i)=mass(i)+0.5d0*(mass(i+1)-mass(i))
     ! derive physical cell masses from adjacent mass boundaries
     delta_mass(i)=mass(i+1)-mass(i)
  enddo
  cmass(imax)=mass(imax)+0.5d0*(mass(imax)-mass(imax-1))
  delta_mass(imax)=delta_mass(imax-1)
  do i=1,imax-1
     delta_cmass(i)=cmass(i+1)-cmass(i)
  enddo
  delta_cmass(imax)=delta_cmass(imax-1)

  write(*,*) ' Fixed legacy wind-grid allocation:'
  write(*,*) '  ejecta/CSM boundary index = ',i_surface
  write(*,*) '  CSM hydro intervals       = ',n_csm
  write(*,*) '  density-bridge intervals  = ',n_bridge

end subroutine redistribute_wind_outer_grid

  ! retain the generated radius-mass relation and recompute cell densities
subroutine map_wind_geometry(profile_zones,pmass,pradius,pvel)

  use blmod, only: mass, r, cr, rho, vel, delta_mass
  use parameters, only: imax, csm_blend_zones, csm_wind_velocity_kms, outdir, &
       mass_scheme, csm
  use csm_switch, only: csm_generated_surface_mass_g, csm_generated_outer_radius_cm
  use physical_constants, only: pi
  implicit none

  integer, intent(in) :: profile_zones
  real*8, intent(in) :: pmass(profile_zones), pradius(profile_zones), pvel(profile_zones)

  integer :: i, i_surface, i_end, u
  integer :: i_density_end, i_velocity_start, n_velocity_transition
  real*8 :: shell_volume, surface_velocity, wind_velocity
  real*8 :: x, rel_outer_error, pre_anchor_rel_outer_error, max_dv, dv
  real*8 :: outer_abs_error
  real*8 :: r_velocity_start, r_velocity_end
  logical :: attached_mass_loss
  character(len=512) :: diag_path

  attached_mass_loss = (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_density') .and. &
       trim(adjustl(csm)).eq.'active'

  do i=1,imax
     call map_map(r(i),mass(i),pradius,pmass,profile_zones)
  enddo

  outer_abs_error = abs(r(imax)-csm_generated_outer_radius_cm)
  rel_outer_error = outer_abs_error / &
       max(abs(csm_generated_outer_radius_cm),1.0d0)
  pre_anchor_rel_outer_error = rel_outer_error

  if (r(imax-1).ge.csm_generated_outer_radius_cm) then
     write(*,*) '******* CSM hydro mapping has crossed the requested outer radius.'
     write(*,*) '******* requested/generated = ',csm_generated_outer_radius_cm
     write(*,*) '******* penultimate radius  = ',r(imax-1)
     stop 1
  endif
  r(imax) = csm_generated_outer_radius_cm

  do i=1,imax-1
     if (r(i+1).le.r(i)) then
        write(*,*) '******* Non-monotone radius after CSM mass-grid mapping.'
        write(*,*) '******* zone, r_i, r_ip1 = ',i,r(i),r(i+1)
        stop 1
     endif
     ! map density using finite spherical shell volumes
     shell_volume = (4.0d0*pi/3.0d0)*(r(i+1)**3-r(i)**3)
     if (shell_volume.le.0.0d0 .or. delta_mass(i).le.0.0d0) then
        write(*,*) '******* Invalid shell geometry in mass-loss CSM mapping.'
        stop 1
     endif
     rho(i) = delta_mass(i)/shell_volume
     cr(i) = ((r(i)**3+r(i+1)**3)/2.0d0)**(1.0d0/3.0d0)
  enddo
  cr(imax) = r(imax) + (r(imax)-cr(imax-1))

  rel_outer_error = abs(r(imax)-csm_generated_outer_radius_cm) / &
       max(abs(csm_generated_outer_radius_cm),1.0d0)
  if (rel_outer_error.gt.1.0d-14) then
     write(*,*) '******* Internal error: snapped CSM outer radius is inconsistent.'
     stop 1
  endif

  i_surface = imax
  do i=1,imax
     if (mass(i).ge.csm_generated_surface_mass_g) then
        i_surface = i
        exit
     endif
  enddo

  if (i_surface.le.1 .or. i_surface.ge.imax) then
     write(*,*) '******* Could not resolve the ejecta/CSM surface on the hydro grid.'
     write(*,*) '******* surface mass, grid range = ',csm_generated_surface_mass_g, &
          mass(1),mass(imax)
     stop 1
  endif

  call map_map(surface_velocity,csm_generated_surface_mass_g,pvel,pmass,profile_zones)
  wind_velocity = csm_wind_velocity_kms*1.0d5
  vel(i_surface) = surface_velocity
  i_density_end = i_surface
  i_velocity_start = i_surface

  if (trim(adjustl(mass_scheme)).eq.'continuous_mass') then
     vel(i_surface:imax)=surface_velocity
     wind_velocity=surface_velocity
     i_end=i_surface
  else if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
     i_end = i_surface + csm_blend_zones
     if (i_end.gt.imax) then
        write(*,*) '******* Insufficient hydro cells for wind velocity transition.'
        write(*,*) '******* first CSM boundary = ',i_surface
        write(*,*) '******* requested transition cells = ',csm_blend_zones
        write(*,*) '******* available cells = ',imax-i_surface
        stop 1
     endif
     do i=i_surface+1,i_end-1
        x = dble(i-i_surface)/dble(csm_blend_zones)
        vel(i) = surface_velocity + x*(wind_velocity-surface_velocity)
     enddo
     vel(i_end) = wind_velocity
     do i=i_end+1,imax
        vel(i) = wind_velocity
     enddo
  else if (attached_mass_loss) then
     if (trim(adjustl(mass_scheme)).eq.'constant_density') then
        i_density_end=i_surface
     else
        i_density_end=min(i_surface+2,imax-1)
     endif
     i_velocity_start=i_density_end
     n_velocity_transition=min(csm_blend_zones,max(2,imax-i_velocity_start-2))
     i_end=i_velocity_start+n_velocity_transition
     if (i_end.gt.imax) then
        write(*,*) '******* Insufficient hydro cells for separated attached-CSM velocity bridge.'
        write(*,*) '******* density bridge end = ',i_density_end
        write(*,*) '******* effective velocity cells = ',n_velocity_transition
        write(*,*) '******* available cells = ',imax-i_density_end
        stop 1
     endif

     do i=i_surface+1,i_velocity_start
        vel(i)=surface_velocity
     enddo
        r_velocity_start=r(i_velocity_start)
        r_velocity_end=r(i_end)
        if (r_velocity_end.le.r_velocity_start) then
           write(*,*) '******* Invalid separated attached-CSM velocity-bridge radius span.'
           stop 1
        endif
        do i=i_velocity_start+1,i_end-1
           x=(r(i)-r_velocity_start)/(r_velocity_end-r_velocity_start)
           x=max(0.0d0,min(1.0d0,x))
           vel(i)=surface_velocity+x*(wind_velocity-surface_velocity)
        enddo
     vel(i_end)=wind_velocity
     do i=i_end+1,imax
        vel(i)=wind_velocity
     enddo
  else
     i_end = i_surface + csm_blend_zones
     if (i_end.gt.imax) then
        write(*,*) '******* Insufficient hydro cells to resolve csm_blend_zones.'
        write(*,*) '******* first CSM boundary = ',i_surface
        write(*,*) '******* requested blend cells = ',csm_blend_zones
        write(*,*) '******* available cells = ',imax-i_surface
        write(*,*) '******* Increase imax/surface concentration or reduce csm_blend_zones.'
        stop 1
     endif

     if (csm_blend_zones.gt.0) then
        do i=i_surface+1,i_end-1
           x = dble(i-i_surface)/dble(csm_blend_zones)
           vel(i) = surface_velocity + x*(wind_velocity-surface_velocity)
        enddo
     endif
     vel(i_end) = wind_velocity
     do i=i_end+1,imax
        vel(i) = wind_velocity
     enddo
  endif

  max_dv = 0.0d0
  do i=i_surface,imax-1
     dv = abs(vel(i+1)-vel(i))
     max_dv = max(max_dv,dv)
  enddo

  diag_path = trim(adjustl(outdir))//'/generated_csm/hydro_mapping_diagnostics.dat'
  u = 919
  open(unit=u,file=trim(diag_path),status='replace',form='formatted')
  write(u,'(A)') '# Initialization-only CSM hydro-grid mapping diagnostics'
  write(u,'(A,ES24.15)') 'generated_outer_radius_cm = ',csm_generated_outer_radius_cm
  write(u,'(A,ES24.15)') 'mapped_outer_radius_cm = ',r(imax)
  write(u,'(A,ES24.15)') 'outer_radius_relative_error = ',rel_outer_error
  write(u,'(A,ES24.15)') 'pre_anchor_endpoint_interpolation_error = ',pre_anchor_rel_outer_error
  write(u,'(A,I3)') 'outer_radius_endpoint_anchored = ',1
  write(u,'(A,ES24.15)') 'surface_mass_g = ',csm_generated_surface_mass_g
  write(u,'(A,I10)') 'first_csm_boundary = ',i_surface
  write(u,'(A,I10)') 'csm_hydro_intervals = ',imax-i_surface
  write(u,'(A,ES24.15)') 'mapped_csm_mass_g = ',mass(imax)-csm_generated_surface_mass_g
  write(u,'(A,I10)') 'density_attachment_end_boundary = ',i_density_end
  write(u,'(A,I10)') 'velocity_transition_start_boundary = ',i_velocity_start
  write(u,'(A,I10)') 'transition_end_boundary = ',i_end
  if (trim(adjustl(mass_scheme)).eq.'continuous_mass') then
     write(u,'(A,I10)') 'velocity_transition_hydro_cells = ',0
     write(u,'(A)') 'velocity_transition_scheme = surface_comoving_unchanged'
  else if (trim(adjustl(mass_scheme)).eq.'non_constant_mass_loss' .or. &
       trim(adjustl(mass_scheme)).eq.'constant_mass_loss') then
     write(u,'(A,I10)') 'velocity_transition_hydro_cells = ',csm_blend_zones
     write(u,'(A)') 'velocity_transition_scheme = hydro_cell_linear_minimax'
  else
     if (attached_mass_loss) then
        write(u,'(A,I10)') 'velocity_transition_hydro_cells = ',n_velocity_transition
        write(u,'(A,I10)') 'requested_blend_zones = ',csm_blend_zones
        write(u,'(A)') 'velocity_transition_scheme = separated_physical_radius_linear'
     else
        write(u,'(A,I10)') 'velocity_transition_hydro_cells = ',csm_blend_zones
        write(u,'(A)') 'velocity_transition_scheme = hydro_cell_linear_minimax'
     endif
  endif
  write(u,'(A,ES24.15)') 'surface_velocity_cms = ',surface_velocity
  write(u,'(A,ES24.15)') 'wind_velocity_cms = ',wind_velocity
  write(u,'(A,ES24.15)') 'max_adjacent_velocity_delta_cms = ',max_dv
  close(u)

end subroutine map_wind_geometry

  ! recover radial boundaries from fixed cell masses and densities
subroutine integrate_radius_initial

  use blmod, only: r, rho, cr, delta_mass
  use parameters
  use physical_constants
  implicit none

  integer :: i


  do i=1,imax-1
   r(i+1) = ( 3.0d0/(4.0d0*pi) * delta_mass(i)/rho(i) + r(i)**3 )**(1.0d0/3.0d0)

   if(rho(i).le.0.0d0) then
      error stop "Negative density. Profile extends not to large enough radii."
   endif
  enddo


  do i=1,imax-1
    cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo

  cr(imax) = r(imax) + (r(imax) - cr(imax-1))


end subroutine integrate_radius_initial
