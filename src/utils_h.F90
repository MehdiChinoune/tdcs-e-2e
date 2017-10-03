MODULE utils
  USE constants, ONLY : RP
  IMPLICIT NONE
  ! Common Variables
  !REAL(KIND=RP), PROTECTED :: fac(0:34), lnfac(0:400)
  REAL(KIND=RP) :: fac(0:34), lnfac(0:400), fak(0:400)
  logical :: fac_called = .false.

  INTERFACE

    MODULE SUBROUTINE factorial()
    END SUBROUTINE factorial

    ELEMENTAL REAL(KIND=RP) MODULE FUNCTION norm_fac(e,n)
      REAL(KIND=RP), INTENT(IN) :: e
      INTEGER      , INTENT(IN) :: n
    END FUNCTION norm_fac

    ELEMENTAL REAL(KIND=RP) MODULE FUNCTION y1y2y3(l1,l2,l3,m1,m2,m3)
      INTEGER,INTENT(IN) :: l1,l2,l3,m1,m2,m3
    END FUNCTION y1y2y3

    MODULE SUBROUTINE ode_second_dw(km,lmax,rc,z,f,s,delta)
      INTEGER, INTENT(IN) :: lmax, z
      REAL(KIND=RP), INTENT(IN)  :: f(0:,0:),rc,km
      REAL(KIND=RP), INTENT(OUT) :: s(0:,0:),delta(0:lmax)
    END SUBROUTINE ode_second_dw

    ELEMENTAL REAL(KIND=RP) MODULE FUNCTION Uij(ni, ei, nj, ej, r)
      REAL(KIND=RP), INTENT(IN) :: ei, ej, r
      INTEGER      , INTENT(IN) :: ni, nj
    END FUNCTION Uij

    MODULE SUBROUTINE calculate_U(Atom, Orbit, r, U )
      CHARACTER(LEN=2), INTENT(IN) :: Atom, Orbit
      REAL(KIND=RP)   , INTENT(IN) :: r(:)
      REAL(KIND=RP)   , INTENT(OUT) :: U(:)
    END SUBROUTINE calculate_U

    MODULE SUBROUTINE INTRPL(X, Y, U, V )
      REAL(RP), INTENT(IN) :: X(:), Y(:), U(:)
      REAL(RP), INTENT(OUT) :: V(:)
    END SUBROUTINE INTRPL

  END INTERFACE

END MODULE utils
