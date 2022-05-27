submodule(input) input
  implicit none

contains

  module subroutine read_input(in_unit, Ei, Es, Ee, thetas, step, Atom, Orbit, exchange)
    integer, intent(in) :: in_unit
    character(len=2), intent(out) :: Atom, Orbit
    real(wp), intent(out) :: Ei, Es, Ee, thetas
    integer, intent(out) :: step(3)
    integer, intent(out), optional :: exchange
    !
    read( in_unit, * ) Atom
    read( in_unit, * ) Orbit
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

  module subroutine read_orbit(orbit_file, nelec, lo, no, n, a, e )
    character(len=5), intent(in)  :: orbit_file
    integer, intent(out) :: nelec, lo, no
    integer, allocatable, intent(out) :: n(:)
    real(wp), allocatable, intent(out) :: a(:), e(:)

    integer :: INPUT

    open( newunit = INPUT, file = 'data/'//orbit_file//'.dat', status = 'old' &
      , action = 'read')

    read( INPUT, * ) nelec
    read( INPUT, * ) lo
    read( INPUT, * ) no
    allocate ( a(no), e(no), n(no) )
    read( INPUT, * ) n
    read( INPUT, * ) a
    read( INPUT, * ) e

    close(INPUT)

  end subroutine read_orbit

end submodule input
