module input
  use constants, only : wp
  use types, only: orbit
  implicit none

  interface

    module subroutine read_fdcs_input(in_unit, Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange)
      integer, intent(in) :: in_unit
      real(wp), intent(out) :: Ei, Es, Ee, thetas
      integer, intent(out) :: step(3)
      integer, optional, intent(out) :: exchange
      character(len=2), intent(out) :: Atom_name, Orbit_name
    end subroutine read_fdcs_input

    module subroutine read_tcs_input(in_unit, Atom_name, Orbit_name, exchange)
      integer, intent(in) :: in_unit
      integer, optional, intent(out) :: exchange
      character(len=2), intent(out) :: Atom_name, Orbit_name
    end subroutine read_tcs_input

    module subroutine read_orbit(orbit_file, orbit_target )
      character(*), intent(in)  :: orbit_file
      type(orbit), intent(out) :: orbit_target
    end subroutine read_orbit

  end interface

end module input
