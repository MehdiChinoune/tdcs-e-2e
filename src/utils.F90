SUBMODULE (utils) utils
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE factorial()
    INTEGER :: i

    IF(fac_called) RETURN

    fac(0) = 1._RP
    DO i=1,34
      fac(i) = i*fac(i-1)
    END DO

    lnfac(0:34) = LOG( fac )
    DO i = 35, 400
      lnfac(i) = lnfac(i-1) +LOG( REAL(i,KIND=RP) )
    END DO

    fac_called = .true.

  END SUBROUTINE

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION norm_fac(e,n)
    REAL(KIND=RP), INTENT(IN) :: e
    INTEGER      , INTENT(IN) :: n

    norm_fac = SQRT( (2.*e)**(2*n+1) / fac(2*n) )

  END FUNCTION

  MODULE SUBROUTINE read_input(in_unit, Ei, Es, Ee, thetas, step, Atom, Orbit)
    INTEGER, INTENT(IN) :: in_unit
    REAL(KIND=RP)   , INTENT(OUT) :: Ei, Es, Ee, thetas
    INTEGER         , INTENT(OUT) :: step(3)
    CHARACTER(LEN=2), INTENT(OUT) :: Atom, Orbit

    READ( in_unit, * ) Atom
    READ( in_unit, * ) Orbit
    READ( in_unit, * ) Ei, Es, Ee
    READ( in_unit, * ) thetas
    READ( in_unit, * ) step

  END SUBROUTINE read_input

  MODULE SUBROUTINE read_orbit(orbit_file, lo, no, n, a, e )
    CHARACTER(LEN=5), INTENT(IN)  :: orbit_file
    INTEGER         , INTENT(OUT) :: lo, no
    INTEGER, ALLOCATABLE, INTENT(OUT) :: n(:)
    REAL(KIND=RP), ALLOCATABLE, INTENT(OUT) :: a(:), e(:)

    INTEGER :: IN

    OPEN( newunit=IN, FILE='Data/'//orbit_file//'.dat', STATUS='old', ACTION='read')

    READ( IN, * ) lo
    READ( IN, * ) no
    ALLOCATE ( a(no), e(no), n(no) )
    READ( IN, * ) n
    READ( IN, * ) a
    READ( IN, * ) e
    CLOSE(IN)

  END SUBROUTINE read_orbit

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION y1y2y3(l1, l2, l3, m1, m2, m3 )
    USE constants ,ONLY: pi
    USE ieee_arithmetic ,ONLY: ieee_is_nan, ieee_is_finite
    INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
    INTEGER :: t
    REAL(KIND=RP) :: s0, s1, cst_1, cst_0

    IF(.NOT. fac_called ) ERROR STOP 'you should call factorial before using y1y2y3'

    y1y2y3=0._RP
    IF( MOD(l1+l2+l3,2)/=0 .or. m1+m2+m3/=0 .or. l3<ABS(l1-l2) .or. l3>l1+l2 .or. ABS(m1)>l1 &
      .or. ABS(m2)>l2 .or. ABS(m3)>l3  ) RETURN

    !  / l_1 l_2 l_3 \
    !  |             |
    !  \ m_1 m_2 m_3 /
    s1 = 0._RP
    cst_1 = 0.5_rp*(lnfac(l1+m1) +lnfac(l1-m1) +lnfac(l2+m2) +lnfac(l2-m2) +lnfac(l3+m3) &
      +lnfac(l3-m3) +lnfac(l2+l3-l1) +lnfac(l3+l1-l2) &
      -0.7*lnfac(l1+l2+l3+1) )
    DO t = MAX(0, l2-l3-m1, l1-l3+m2 ), MIN(l1+l2-l3, l1-m1, l2+m2 )
      s1 = s1 +(-1)**t* EXP( cst_1 -( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+m1+t) &
        +lnfac(l3-l1-m2+t) +lnfac(l1-m1-t) +lnfac(l2+m2-t) ) )
    END DO

    ! / l_1 l_2 l_3 \
    ! |             |
    ! \ 0   0   0   /
    s0 = 0._RP
    cst_0 = lnfac(l1) +lnfac(l2) +lnfac(l3) +0.5_rp*( lnfac(l2+l3-l1) &
      +lnfac(l3+l1-l2) -0.7*lnfac(l1+l2+l3+1) )
    DO t = MAX(0, l2-l3, l1-l3 ), MIN(l1+l2-l3, l1, l2 )
      s0 = s0 +(-1)**t *EXP( cst_0 -( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+t) +lnfac(l3-l1+t) &
        +lnfac(l1-t) +lnfac(l2-t) ) )
    END DO

    y1y2y3 = (-1)**m3 *SQRT( (2*l1+1)*(2*l2+1)*(2*l3+1)/(4.*pi) ) *s1 &
      *( EXP(lnfac(l1+l2-l3)-0.3*lnfac(l1+l2+l3+1) ) *s0 )
    IF( ieee_is_nan(y1y2y3) .or. (.NOT. ieee_is_finite(y1y2y3)) ) ERROR STOP 'y1y2y3 overflow'

  END FUNCTION y1y2y3


  ! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
  !  Discussion:
  !    Our convention is that the abscissas are numbered from left to right.
  !    The rule is defined on [-1,1].
  !    The integral to approximate:
  !      Integral ( -1 <= X <= 1 ) F(X) dX
  !    The quadrature rule:
  !      Sum ( 1 <= I <= n ) W(I) * F ( X(I) )
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

  MODULE SUBROUTINE clenshaw_curtis( a, b, x, w, n )
    USE constants ,ONLY: pi
    REAL(KIND=RP), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL(KIND=RP), INTENT(OUT), ALLOCATABLE :: w(:), x(:)

    REAL(KIND=RP) :: theta !,bj
    INTEGER :: i, j

    IF ( n < 3 ) THEN
      ERROR STOP 'clenashaw_curtis  error : n < 3 '
    END IF

    ALLOCATE(x(n),w(n))

    !    IF ( n == 1 ) THEN
    !      x(1) = 0._RP
    !      w(1) = 2._RP
    !      RETURN
    !    END IF

    !    DO CONCURRENT( i = 1:n )
    !      x(i) = COS( (n-i) *pi /(n-1) )
    !    END DO
    x(1) = -1._RP
    x(2:n-1) = [ ( COS( (n-i)*pi/(n-1) ), i=2,n-1 ) ]
    x(n) = 1._RP

    IF ( MOD(n,2) == 1 ) THEN
      x((n+1)/2) = 0._RP
    END IF

    w = 1._RP
    DO CONCURRENT( i = 1:n )
      theta = (i-1) *pi/(n-1)

      DO j = 1, (n-1) /2 -1
        !        IF ( 2*j == (n-1) ) THEN
        !          bj = 1._RP
        !        ELSE
        !          bj = 2._RP
        !        END IF
        !        w(i) = w(i) - bj *COS ( 2.*j *theta ) / (4*j*j - 1.)
        w(i) = w(i) - 2.*COS ( 2.*j *theta ) / (4*j*j - 1.)
      END DO

      IF( MOD(n,2) == 0 ) THEN
        w(i) = w(i) -2.*COS ( 2.*( (n-1)/2 ) *theta ) / (4*( (n-1)/2 )**2 - 1.)
      ELSE
        w(i) = w(i) -COS ( (n-1)*theta ) /( (n-1)**2 - 1.)
      END IF

    END DO

    w(2:n-1) = 2.*w(2:n-1)

    w = w /(n-1.)

    x = ( (a+b) + (b-a)*x ) /2.
    w = (b-a)*w/2.

  END SUBROUTINE clenshaw_curtis

  PURE MODULE SUBROUTINE gauleg(a,b,x,w,n)
    use constants ,only: pi
    INTEGER, INTENT(IN) :: n
    REAL(rp), INTENT(IN) :: a,b
    REAL(rp), INTENT(OUT), ALLOCATABLE :: x(:),w(:)
    INTEGER :: i

    ALLOCATE(x(n),w(n))

    x(1:(n+1)/2) = [(COS(pi*(4.*i-1.)/(4.*n+2.)), i=1,(n+1)/2 )]

    CALL pd(x(1:(n+1)/2),w(1:(n+1)/2),n)

    DO CONCURRENT (i=1:(n+1)/2)
      x(n-i+1)=-x(i)
      w(n-i+1)=w(i)
    END DO

    x = ( (a-b)*x +b+a )/2.
    w = (b-a)*w/2.

  END SUBROUTINE gauleg

  PURE SUBROUTINE pd(sx,sw,n)
    INTEGER, INTENT(IN) :: n
    REAL(rp), INTENT(OUT),CONTIGUOUS :: sw(:)
    REAL(rp), INTENT(INOUT),CONTIGUOUS :: sx(:)

    REAL(RP),PARAMETER :: eps=EPSILON(eps)
    REAL(rp), DIMENSION((n+1)/2) :: dp0,dp1,dp2,dp
    INTEGER :: i

    dp2 = 0._RP
    DO
      dp0 = 1._RP
      dp1 = sx
      DO i=1,n-1
        dp2 = ((2.*i+1._RP)*sx*dp1-i*dp0)/(i+1._RP)
        dp0 = dp1
        dp1 = dp2
      ENDDO
      dp = n*(dp0-sx*dp1)/(1._RP-sx**2)
      sx = sx-dp2/dp

      IF( ALL(ABS(dp2/dp)<=eps) ) EXIT
    ENDDO
    sw = 2._rp/((1._rp-sx**2)*dp**2)
  END SUBROUTINE pd

  !> ode_second_dw
  !!
  !! This subroutine solve Equation of the form
  !! s_l''(r) + f_l(r) *s_l(r) = km**2 *s_l(r)

  MODULE SUBROUTINE ode_second_dw(km,lmax,rc,f,s,delta)
    USE special_functions ,ONLY: coul90
    INTEGER, INTENT(IN) :: lmax
    REAL(KIND=RP), INTENT(IN)  :: f(0:,0:),rc,km
    REAL(KIND=RP), INTENT(OUT) :: s(0:,0:),delta(0:lmax)
    REAL(KIND=RP) :: h,rho,eta
    REAL(KIND=RP), DIMENSION(0:lmax) :: jl,gl,jpl,gpl,fn,betap,beta
    INTEGER :: i,ns,ifail,l

    ns=SIZE(s,1)-1
    h=rc/ns

    rho=km*h*(ns-2)

    eta=-1./km
    CALL coul90(rho,eta,0,lmax,jl,gl,jpl,gpl,0,ifail)

    DO l=0,lmax
      s(0,l) = 0.
      s(1,l) = h**(l+1)
      DO i=1,ns-1
        s(i+1,l) = ( (2.+f(i,l)*5.*h**2/6.)*s(i,l) - (1.-f(i-1,l)*h**2/12.)*s(i-1,l) ) &
          /( 1.-f(i+1,l)*h**2/12. )
      END DO
    END DO

    betap = ( s(ns-4,:) -8.*s(ns-3,:) +8.*s(ns-1,:) -s(ns,:) ) / ( 12.*h*km )
    beta = s(ns-2,:)
    delta = ATAN( (-jl*betap+jpl*beta) /(gl*betap-gpl*beta ) )
    fn = (COS(delta)*jl+SIN(delta)*gl) /beta

    DO l=0,lmax
      s(:,l)=s(:,l)*fn(l)
    END DO

  END SUBROUTINE ode_second_dw

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION Uij(ni, ei, nj, ej, r)
    REAL(KIND=RP), INTENT(IN) :: ei, ej, r
    INTEGER      , INTENT(IN) :: ni, nj
    INTEGER :: k
    REAL(KIND=RP) :: a

    a = 1._RP
    DO CONCURRENT( k = 1:ni+nj-1 )
      a = a + ((ei+ej)*r)**k *(ni+nj-k) /( (ni+nj)*fac(k) )
    END DO
    Uij = a *EXP(-(ei+ej)*r) *fac(ni+nj) *ei**ni *ej**nj *SQRT(ei*ej) &
      /( SQRT( fac(2*ni) *fac(2*nj) ) *( (ei+ej)/2. )**(ni+nj+1) )

  END FUNCTION Uij

  MODULE SUBROUTINE calculate_U(Atom, Orbit, r, U )
    CHARACTER(LEN=2), INTENT(IN) :: Atom, Orbit
    REAL(KIND=RP)   , INTENT(IN) :: r(:)
    REAL(KIND=RP)   , INTENT(OUT) :: U(:)

    INTEGER :: IN

    OPEN( newunit=IN, FILE='Data/'//Atom//'.dat', STATUS='old', ACTION='read')

    U = 0._RP
    DO
      BLOCK
        REAL(KIND=rp), ALLOCATABLE :: a(:),e(:)
        INTEGER, ALLOCATABLE :: n(:)
        INTEGER :: lo, no, i1, i2, nocup
        CHARACTER(LEN=2) :: orbit_i

        READ(IN, FMT=*, IOSTAT=lo ) orbit_i
        IF(lo<0) EXIT

        CALL read_orbit(Atom//'_'//orbit_i, lo, no, n, a, e)
        nocup = 2*(2*lo+1)
        IF(orbit_i==Orbit) nocup = nocup -1

        DO CONCURRENT( i1 = 1:no, i2 = 1:no )
          U = U + nocup *a(i1) *a(i2) *Uij( n(i1), e(i1), n(i2), e(i2), r )
        END DO

      END BLOCK
    END DO

    CLOSE(IN)

  END SUBROUTINE calculate_U

END SUBMODULE utils
