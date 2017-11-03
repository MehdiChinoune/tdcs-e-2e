SUBMODULE (input) input
  implicit none

CONTAINS

  MODULE SUBROUTINE read_input(in_unit, Ei, Es, Ee, thetas, step, Atom, Orbit, exchange)
    INTEGER, INTENT(IN) :: in_unit
    REAL(KIND=RP)   , INTENT(OUT) :: Ei, Es, Ee, thetas
    INTEGER         , INTENT(OUT) :: step(3)
    INTEGER, OPTIONAL, INTENT(OUT) :: exchange
    CHARACTER(LEN=2), INTENT(OUT) :: Atom, Orbit

    READ( in_unit, * ) Atom
    READ( in_unit, * ) Orbit
    READ( in_unit, * ) Ei, Es, Ee
    READ( in_unit, * ) thetas
    READ( in_unit, * ) step
    READ( in_unit, * ) exchange

  END SUBROUTINE read_input

  MODULE SUBROUTINE read_orbit(orbit_file, nelec, lo, no, n, a, e )
    CHARACTER(LEN=5), INTENT(IN)  :: orbit_file
    INTEGER         , INTENT(OUT) :: nelec, lo, no
    INTEGER, ALLOCATABLE, INTENT(OUT) :: n(:)
    REAL(KIND=RP), ALLOCATABLE, INTENT(OUT) :: a(:), e(:)

    INTEGER :: IN

    OPEN( newunit=IN, FILE='Data/'//orbit_file//'.dat', STATUS='old', ACTION='read')

    READ( IN, * ) nelec
    READ( IN, * ) lo
    READ( IN, * ) no
    ALLOCATE ( a(no), e(no), n(no) )
    READ( IN, * ) n
    READ( IN, * ) a
    READ( IN, * ) e
    CLOSE(IN)

  END SUBROUTINE read_orbit

END SUBMODULE input
