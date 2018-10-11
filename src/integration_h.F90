MODULE integration
  USE constants ,only: RP
  IMPLICIT NONE

  INTERFACE

    MODULE SUBROUTINE clenshaw_curtis( a, b, x, w, n )
      REAL(KIND=RP), INTENT(IN) :: a, b
      INTEGER, INTENT(IN) :: n
      REAL(KIND=RP), INTENT(OUT), ALLOCATABLE :: w(:), x(:)
    END SUBROUTINE clenshaw_curtis

    PURE MODULE SUBROUTINE gauleg(a,b,x,w,n)
      INTEGER, INTENT(IN) :: n
      REAL(rp), INTENT(IN) :: a,b
      REAL(rp), INTENT(OUT), ALLOCATABLE :: x(:),w(:)
    END SUBROUTINE gauleg

  END INTERFACE

END MODULE integration