MODULE integration
  USE constants ,ONLY: RP
  IMPLICIT NONE

CONTAINS

  ! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
  !  Discussion:
  !    Our convention is that the abscissas are numbered from left to right.
  !    The rule is defined on [-1,1].
  !    The integral to approximate:
  !      Integral ( -1 <= X <= 1 ) F(X) dX
  !    The quadrature rule:
  !      Sum ( 1 <= I <= n ) W(I)*F ( X(I) )
  !
  !  Licensing:
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !    15 February 2009
  !
  !  Author:
  !    John Burkardt

  !  Parameters:
  !    Input, integer ( kind = 4 ) n, the order of the rule.
  !    1 <= ORDER.
  !    Output, real ( kind = 8 ) X(n), the abscissas.
  !    Output, real ( kind = 8 ) W(n), the weights.

  SUBROUTINE clenshaw_curtis( a, b, x, w, n )
    USE constants ,ONLY: pi
    REAL(RP), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL(RP), INTENT(OUT), ALLOCATABLE :: w(:), x(:)

    REAL(RP) :: theta !,bj
    INTEGER :: i, j

!    IF ( n<3 ) THEN
!      ERROR STOP 'clenashaw_curtis  error : n<3 '
!    END IF

    ALLOCATE(x(n),w(n))

    x(1) = -1._RP
    x(2:n-1) = [ ( COS( (n-i)*pi/(n-1) ), i = 2,n-1 ) ]
    x(n) = 1._RP

    IF ( MOD(n,2)==1 ) THEN
      x((n+1)/2) = 0._RP
    END IF

    w = 1._RP
    DO i = 1,n
      theta = (i-1)*pi/(n-1)

      DO j = 1, (n-1) /2 -1
        w(i) = w(i) -2.*COS( 2.*j*theta ) / (4*j*j -1.)
      END DO

      IF( MOD(n,2)==0 ) THEN
        w(i) = w(i) -2.*COS( 2.*( (n-1)/2 )*theta ) / (4*( (n-1)/2 )**2 -1.)
      ELSE
        w(i) = w(i) -COS( (n-1)*theta ) /( (n-1)**2 -1.)
      END IF

    END DO

    w(2:n-1) = 2.*w(2:n-1)

    w = w /(n-1.)

    x = ( (a+b) +(b-a)*x ) /2.
    w = (b-a)*w/2.

  END SUBROUTINE clenshaw_curtis

  PURE SUBROUTINE gauleg(a,b,x,w,n)
    use constants ,only: pi
    INTEGER, INTENT(IN) :: n
    REAL(RP), INTENT(IN) :: a,b
    REAL(RP), INTENT(OUT), ALLOCATABLE :: x(:),w(:)
    INTEGER :: i

    ALLOCATE(x(n),w(n))

    x(1:(n+1)/2) = [(COS(pi*(4.*i-1.)/(4.*n+2.)), i = 1,(n+1)/2 )]

    CALL pd(x(1:(n+1)/2),w(1:(n+1)/2),n)

    DO i = 1,(n+1)/2
      x(n-i+1) = -x(i)
      w(n-i+1) = w(i)
    END DO

    x = ( (a-b)*x +b+a )/2.
    w = (b-a)*w/2.

  END SUBROUTINE gauleg

  PURE SUBROUTINE pd(sx,sw,n)
    INTEGER, INTENT(IN) :: n
    REAL(RP), INTENT(OUT), CONTIGUOUS :: sw(:)
    REAL(RP), INTENT(INOUT), CONTIGUOUS :: sx(:)

    REAL(RP),PARAMETER :: eps = EPSILON(eps)
    REAL(RP), DIMENSION((n+1)/2) :: dp0,dp1,dp2,dp
    INTEGER :: i

    dp2 = 0._RP
    DO
      dp0 = 1._RP
      dp1 = sx
      DO i = 1,n-1
        dp2 = ((2.*i+1._RP)*sx*dp1-i*dp0)/(i+1._RP)
        dp0 = dp1
        dp1 = dp2
      ENDDO
      dp = n*(dp0-sx*dp1)/(1._RP-sx**2)
      sx = sx-dp2/dp

      IF( ALL(ABS(dp2/dp)<=eps) ) EXIT
    ENDDO
    sw = 2._RP/((1._RP-sx**2)*dp**2)
  END SUBROUTINE pd

END MODULE integration
