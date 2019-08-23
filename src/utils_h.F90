MODULE utils
  USE constants, ONLY : RP
  IMPLICIT NONE

  INTERFACE

    ELEMENTAL REAL(RP) MODULE FUNCTION norm_fac(e,n)
      REAL(RP), INTENT(IN) :: e
      INTEGER , INTENT(IN) :: n
    END FUNCTION

    ELEMENTAL REAL(RP) MODULE FUNCTION y1y2y3(l1, l2, l3, m1, m2, m3 )
      INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
    END FUNCTION y1y2y3

    !> ode_second_dw
    !!
    !! This subroutine solve Equation of the form
    !! s_l''(r) +f_l(r)*s_l(r) = km**2*s_l(r)

    MODULE SUBROUTINE ode_second_dw(km, lmax, rc, z, f, s, delta )
      INTEGER, INTENT(IN) :: lmax, z
      REAL(RP), INTENT(IN)  :: f(0:,0:), rc, km
      REAL(RP), INTENT(OUT) :: s(0:,0:), delta(0:lmax)
    END SUBROUTINE ode_second_dw

    MODULE SUBROUTINE calculate_U(Atom, Orbit, r, U, state )
      CHARACTER(LEN=2), INTENT(IN) :: Atom, Orbit
      REAL(RP)   , INTENT(IN) :: r(:)
      REAL(RP)   , INTENT(OUT) :: U(:)
      INTEGER :: state
    END SUBROUTINE calculate_U

    PURE MODULE SUBROUTINE INTRPL(X, Y, U, V )
      REAL(RP), INTENT(IN) :: X(:), Y(:), U(:)
      REAL(RP), INTENT(OUT) :: V(:)
    END SUBROUTINE INTRPL

  END INTERFACE

END MODULE utils