MODULE input
  USE constants ,only: RP
  IMPLICIT NONE

  INTERFACE

    MODULE SUBROUTINE read_input(in_unit, Ei, Es, Ee, thetas, step, Atom, Orbit, exchange)
      INTEGER, INTENT(IN) :: in_unit
      REAL(KIND=RP)   , INTENT(OUT) :: Ei, Es, Ee, thetas
      INTEGER         , INTENT(OUT) :: step(3)
      INTEGER, OPTIONAL, INTENT(OUT) :: exchange
      CHARACTER(LEN=2), INTENT(OUT) :: Atom, Orbit
    END SUBROUTINE read_input

    MODULE SUBROUTINE read_orbit(orbit_file, nelec, lo, no, n, a, e )
      CHARACTER(LEN=5), INTENT(IN)  :: orbit_file
      INTEGER         , INTENT(OUT) :: nelec, lo, no
      INTEGER, ALLOCATABLE, INTENT(OUT) :: n(:)
      REAL(KIND=RP), ALLOCATABLE, INTENT(OUT) :: a(:), e(:)
    END SUBROUTINE read_orbit

  END INTERFACE
END MODULE input
