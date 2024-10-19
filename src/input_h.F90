module input
  use constants, only : wp
  implicit none

  interface

    module subroutine read_input(in_unit, Ei, Es, Ee, thetas, step, Atom, Orbit, exchange)
      integer, intent(in) :: in_unit
      real(wp), intent(out) :: Ei, Es, Ee, thetas
      integer, intent(out) :: step(3)
      integer, optional, intent(out) :: exchange
      character(len=2), intent(out) :: Atom, Orbit
    end subroutine read_input

    module subroutine read_orbit(orbit_file, Ie, nelec, lo, no, n, a, e )
      character(len=5), intent(in)  :: orbit_file
      real(wp), intent(out) :: Ie
      integer, intent(out) :: nelec, lo, no
      integer, allocatable, intent(out) :: n(:)
      real(wp), allocatable, intent(out) :: a(:), e(:)
    end subroutine read_orbit

  end interface

end module input
