! same-executable handoff of the live grid and eos state
module staged_state
  use blmod, only: mass, cmass, delta_mass, delta_cmass, r, cr, rho, temp, &
       vel, eps, p, cs2, dpdt, dedt, entropy, p_rad, comp, comp_details, &
       ye, abar, metallicity, zav, ncomps, ion_fractions, free_electron_frac
  use parameters, only: imax, eoskey, csm, restart_state_file
  implicit none
  private
  public :: write_staged_state, restore_staged_state
contains
  subroutine write_staged_state(path)
    character(len=*), intent(in) :: path
    integer :: u
    open(newunit=u,file=path,status='replace',access='stream',form='unformatted')
    write(u) 'STRATA_STATE_002', imax, ncomps, eoskey
    write(u) comp_details
    write(u) mass, cmass, delta_mass, delta_cmass, r, cr, rho, temp, vel, &
         eps, p, cs2, dpdt, dedt, entropy, p_rad, comp, ye, abar, metallicity, zav, &
         ion_fractions, free_electron_frac
    close(u)
  end subroutine

  subroutine restore_staged_state(restored)
    logical, intent(out) :: restored
    character(len=16) :: magic
    integer :: u, stat, nz, nc, ek
    real*8, allocatable :: species(:,:)
    restored=.false.
    if (len_trim(restart_state_file).eq.0) return
    if (trim(csm).ne.'none') error stop 'Exact staged state requires a coast'
    open(newunit=u,file=trim(restart_state_file),status='old',access='stream',form='unformatted',iostat=stat)
    if (stat.ne.0) error stop 'Cannot open staged state'
    read(u,iostat=stat) magic, nz, nc, ek
    if (stat.ne.0) error stop 'Cannot read staged state header'
    if (magic.ne.'STRATA_STATE_002' .or. nz.ne.imax .or. nc.ne.ncomps .or. ek.ne.eoskey) &
         error stop 'Incompatible staged state'
    allocate(species(size(comp_details,1),size(comp_details,2)))
    read(u,iostat=stat) species
    if (stat.ne.0) error stop 'Cannot read staged species'
    if (any(species.ne.comp_details)) error stop 'Staged species mismatch'
    read(u,iostat=stat) mass, cmass, delta_mass, delta_cmass, r, cr, rho, temp, vel, &
         eps, p, cs2, dpdt, dedt, entropy, p_rad, comp, ye, abar, metallicity, zav, &
         ion_fractions, free_electron_frac
    if (stat.ne.0) error stop 'Incomplete staged state'
    close(u)
    restored=.true.
    write(*,'(A)') ' Preserving exact coast grid and thermodynamic state.'
  end subroutine
end module staged_state
