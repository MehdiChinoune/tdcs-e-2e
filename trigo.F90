!> Trigonometric Transformations

!> This module contains some trigonometric transformations
SUBMODULE (trigo) trigo

CONTAINS

  PURE MODULE SUBROUTINE spher2cartez( km, theta, phi, k )
    REAL(KIND=RP), INTENT(IN)  :: km, phi, theta
    REAL(KIND=RP), INTENT(OUT) :: k(3)
    k = km *[ SIN(theta)*COS(phi), SIN(theta)*SIN(phi), COS(theta) ]
  END SUBROUTINE spher2cartez

  PURE MODULE SUBROUTINE cartez2spher( k, km, theta, phi )
    REAL(KIND=RP), INTENT(IN)  :: k(3)
    REAL(KIND=RP), INTENT(OUT) :: km, theta, phi
    km = NORM2(k)
    theta = ACOS( k(3)/km )
    phi = ATAN2( k(2), k(1) )
  END SUBROUTINE cartez2spher

END SUBMODULE trigo
