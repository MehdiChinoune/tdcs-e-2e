submodule(input) input
  implicit none

contains

  module subroutine read_input(in_unit, Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange)
    integer, intent(in) :: in_unit
    character(len=2), intent(out) :: Atom_name, Orbit_name
    real(wp), intent(out) :: Ei, Es, Ee, thetas
    integer, intent(out) :: step(3)
    integer, intent(out), optional :: exchange
    !
    read( in_unit, * ) Atom_name
    read( in_unit, * ) Orbit_name
    read( in_unit, * ) Ei, Es, Ee
    read( in_unit, * ) thetas
    read( in_unit, * ) step
    if( present(exchange) ) then
      read( in_unit, * ) exchange
    end if
    ! Rewind, to use it more than once.
    rewind( in_unit )
    !
  end subroutine read_input

  module subroutine read_orbit(orbit_file, orbit_target )
    use constants, only: ev
    character(len=5), intent(in)  :: orbit_file
    type(orbit), intent(out) :: orbit_target
    !
    real(wp) :: Ie
    integer :: INPUT

    open( newunit = INPUT, file = 'data/'//orbit_file//'.dat', status = 'old' &
      , action = 'read')

    read( INPUT, * ) orbit_target%nelec
    read( INPUT, * ) orbit_target%l
    read( INPUT, * ) orbit_target%nf
    allocate ( orbit_target%a(orbit_target%nf), orbit_target%e(orbit_target%nf), orbit_target%n(orbit_target%nf) )
    read( INPUT, * ) orbit_target%n
    read( INPUT, * ) orbit_target%a
    read( INPUT, * ) orbit_target%e
    read( INPUT, * ) Ie
    orbit_target%Ie = Ie*eV

    close(INPUT)

  end subroutine read_orbit

end submodule input
