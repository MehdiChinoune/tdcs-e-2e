!> Trigonometric Transformations

!> This module contains some trigonometric transformations
MODULE trigo
  USE constants ,ONLY : RP
  IMPLICIT NONE

  INTERFACE

    PURE MODULE SUBROUTINE spher2cartez( km ,theta ,phi ,k )
      REAL(KIND=RP) ,INTENT(IN)  :: km ,phi ,theta
      REAL(KIND=RP) ,INTENT(OUT) :: k(3)
    END SUBROUTINE spher2cartez

    PURE MODULE SUBROUTINE cartez2spher( k ,km ,theta ,phi )
      REAL(KIND=RP) ,INTENT(IN)  :: k(3)
      REAL(KIND=RP) ,INTENT(OUT) :: km ,theta ,phi
    END SUBROUTINE cartez2spher

  END INTERFACE

END MODULE trigo
