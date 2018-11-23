MODULE integration
  USE constants ,ONLY: RP
  IMPLICIT NONE

  INTERFACE

    MODULE SUBROUTINE clenshaw_curtis( a, b, x, w, n )
      REAL(RP), INTENT(IN) :: a, b
      INTEGER, INTENT(IN) :: n
      REAL(RP), INTENT(OUT), ALLOCATABLE :: w(:), x(:)
    END SUBROUTINE clenshaw_curtis

    PURE MODULE SUBROUTINE gauleg(a,b,x,w,n)
      INTEGER, INTENT(IN) :: n
      REAL(RP), INTENT(IN) :: a,b
      REAL(RP), INTENT(OUT), ALLOCATABLE :: x(:),w(:)
    END SUBROUTINE gauleg

    PURE MODULE SUBROUTINE pd(sx,sw,n)
      INTEGER, INTENT(IN) :: n
      REAL(RP), INTENT(OUT), CONTIGUOUS :: sw(:)
      REAL(RP), INTENT(INOUT), CONTIGUOUS :: sx(:)
    END SUBROUTINE pd

  END INTERFACE

END MODULE integration
