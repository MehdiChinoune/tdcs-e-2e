!> Trigonometric Transformations

!> This module contains some trigonometric transformations
MODULE trigo
  USE constants ,ONLY: RP
  IMPLICIT NONE

CONTAINS

  PURE REAL(RP) FUNCTION nrm2(a)
    REAL(RP), INTENT(IN) :: a(:)
#if !defined(__GNUC) || !defined(__INTEL)
    nrm2 = SQRT( SUM(ABS(a)**2) )
#else
    nrm2 = NORM2(a)
#endif
  END FUNCTION

  PURE SUBROUTINE spher2cartez( km, theta, phi, k )
    REAL(KIND=RP), INTENT(IN)  :: km, phi, theta
    REAL(KIND=RP), INTENT(OUT) :: k(3)
    k = km *[ SIN(theta)*COS(phi), SIN(theta)*SIN(phi), COS(theta) ]
  END SUBROUTINE spher2cartez

  PURE SUBROUTINE cartez2spher( k, km, theta, phi )
    REAL(KIND=RP), INTENT(IN)  :: k(3)
    REAL(KIND=RP), INTENT(OUT) :: km, theta, phi
    km = nrm2(k)
    theta = ACOS( k(3)/km )
    phi = ATAN2( k(2), k(1) )
  END SUBROUTINE cartez2spher

END MODULE trigo
