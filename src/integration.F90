submodule(integration) integration
  implicit none

contains

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

  module subroutine clenshaw_curtis( a, b, x, w, n )
    use constants ,only: pi
    real(RP), intent(in) :: a, b
    integer, intent(in) :: n
    real(RP), intent(out), allocatable :: w(:), x(:)

    real(RP) :: theta !,bj
    integer :: i, j

!    IF ( n<3 ) THEN
!      ERROR STOP 'clenashaw_curtis  error : n<3 '
!    END IF

    allocate(x(n),w(n))

    x(1) = -1._RP
    x(2:n-1) = [ ( cos( (n-i)*pi/(n-1) ), i = 2,n-1 ) ]
    x(n) = 1._RP

    if ( mod(n,2)==1 ) then
      x((n+1)/2) = 0._RP
    end if

    w = 1._RP
    do i = 1,n
      theta = (i-1)*pi/(n-1)

      do j = 1, (n-1) /2 -1
        w(i) = w(i) -2.*cos( 2.*j*theta ) / (4*j*j -1.)
      end do

      if( mod(n,2)==0 ) then
        w(i) = w(i) -2.*cos( 2.*( (n-1)/2 )*theta ) / (4*( (n-1)/2 )**2 -1.)
      else
        w(i) = w(i) -cos( (n-1)*theta ) /( (n-1)**2 -1.)
      end if

    end do

    w(2:n-1) = 2.*w(2:n-1)

    w = w /(n-1.)

    x = ( (a+b) +(b-a)*x ) /2.
    w = (b-a)*w/2.

  end subroutine clenshaw_curtis

  pure module subroutine gauleg(a,b,x,w,n)
    use constants ,only: pi
    integer, intent(in) :: n
    real(RP), intent(in) :: a,b
    real(RP), intent(out), allocatable :: x(:),w(:)
    integer :: i

    allocate(x(n),w(n))

    x(1:(n+1)/2) = [(cos(pi*(4.*i-1.)/(4.*n+2.)), i = 1,(n+1)/2 )]

    call pd(x(1:(n+1)/2),w(1:(n+1)/2),n)

    do i = 1,(n+1)/2
      x(n-i+1) = -x(i)
      w(n-i+1) = w(i)
    end do

    x = ( (a-b)*x +b+a )/2.
    w = (b-a)*w/2.

  end subroutine gauleg

  pure module subroutine pd(sx,sw,n)
    integer, intent(in) :: n
    real(RP), intent(out), contiguous :: sw(:)
    real(RP), intent(inout), contiguous :: sx(:)

    real(RP),parameter :: eps = epsilon(eps)
    real(RP), dimension((n+1)/2) :: dp0,dp1,dp2,dp
    integer :: i

    dp2 = 0._RP
    do
      dp0 = 1._RP
      dp1 = sx
      do i = 1,n-1
        dp2 = ((2.*i+1._RP)*sx*dp1-i*dp0)/(i+1._RP)
        dp0 = dp1
        dp1 = dp2
      enddo
      dp = n*(dp0-sx*dp1)/(1._RP-sx**2)
      sx = sx-dp2/dp

      if( all(abs(dp2/dp)<=eps) ) exit
    enddo
    sw = 2._RP/((1._RP-sx**2)*dp**2)
  end subroutine pd

end submodule integration
