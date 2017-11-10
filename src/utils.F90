MODULE utils
  USE constants ,ONLY: RP
  IMPLICIT NONE

CONTAINS

  ELEMENTAL REAL(KIND=RP) FUNCTION norm_fac(e,n)
    USE special_functions ,ONLY: fac
    REAL(KIND=RP), INTENT(IN) :: e
    INTEGER      , INTENT(IN) :: n

    norm_fac = SQRT( (2.*e)**(2*n+1) / fac(2*n) )

  END FUNCTION

  ELEMENTAL REAL(KIND=RP) FUNCTION y1y2y3(l1, l2, l3, m1, m2, m3 )
    USE constants ,ONLY: pi
    !USE ieee_arithmetic ,ONLY: ieee_is_nan, ieee_is_finite
    USE special_functions ,ONLY: lnfac!, fac_called
    INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
    INTEGER :: t
    REAL(KIND=RP) :: s0, s1, cst_1, cst_0

    !IF(.NOT. fac_called ) ERROR STOP 'you should call factorial before using y1y2y3'

    y1y2y3=0._RP
    IF( MOD(l1+l2+l3,2)/=0 .or. m1+m2+m3/=0 .or. l3<ABS(l1-l2) .or. l3>l1+l2 .or. ABS(m1)>l1 &
      .or. ABS(m2)>l2 .or. ABS(m3)>l3  ) RETURN

    !  / l_1 l_2 l_3 \
    !  |             |
    !  \ m_1 m_2 m_3 /
    s1 = 0._RP
    cst_1 = 0.5_RP*(lnfac(l1+m1) +lnfac(l1-m1) +lnfac(l2+m2) +lnfac(l2-m2) +lnfac(l3+m3) &
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
    cst_0 = lnfac(l1) +lnfac(l2) +lnfac(l3) +0.5_RP*( lnfac(l2+l3-l1) &
      +lnfac(l3+l1-l2) -0.7*lnfac(l1+l2+l3+1) )
    DO t = MAX(0, l2-l3, l1-l3 ), MIN(l1+l2-l3, l1, l2 )
      s0 = s0 +(-1)**t *EXP( cst_0 -( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+t) +lnfac(l3-l1+t) &
        +lnfac(l1-t) +lnfac(l2-t) ) )
    END DO

    y1y2y3 = (-1)**m3 *SQRT( (2*l1+1)*(2*l2+1)*(2*l3+1)/(4.*pi) ) *s1 &
      *( EXP(lnfac(l1+l2-l3)-0.3*lnfac(l1+l2+l3+1) ) *s0 )
    !IF( ieee_is_nan(y1y2y3) .or. (.NOT. ieee_is_finite(y1y2y3)) ) ERROR STOP 'y1y2y3 overflow'

  END FUNCTION y1y2y3

  !> ode_second_dw
  !!
  !! This subroutine solve Equation of the form
  !! s_l''(r) + f_l(r) *s_l(r) = km**2 *s_l(r)

  SUBROUTINE ode_second_dw(km, lmax, rc, z, f, s, delta )
    USE special_functions ,ONLY: coul90, ricbes
    INTEGER, INTENT(IN) :: lmax, z
    REAL(KIND=RP), INTENT(IN)  :: f(0:,0:), rc, km
    REAL(KIND=RP), INTENT(OUT) :: s(0:,0:), delta(0:lmax)
    REAL(KIND=RP) :: Big = SQRT(HUGE(1._RP)), eps = EPSILON(1._RP)
    REAL(KIND=RP) :: h, rho, eta
    REAL(KIND=RP), DIMENSION(0:lmax) :: jl, gl, jpl, gpl, fn, betap, beta
    INTEGER :: i, ns, l, is

    ns = SIZE(s,1)-1
    h = rc/ns

    rho = km*h*(ns-2)

    IF(z/=0) THEN
      eta = -z/km
      CALL coul90(rho, eta, 0, lmax, jl, gl, jpl, gpl, 0 )
    ELSE
      CALL ricbes(rho, lmax, jl, gl, jpl, gpl )
    END IF


    DO l=0,lmax
      s(0,l) = 0._RP
      DO is=1,ns-1
        s(is,l) = (is*h)**(l+1)
        IF(s(is,l)>0._RP) EXIT
      END DO
      DO i=is,ns-1
        s(i+1,l) = ( (2._RP+f(i,l)*5._RP*h**2/6.)*s(i,l) - (1._RP-f(i-1,l)*h**2/12.)*s(i-1,l) ) &
          /( 1._RP-f(i+1,l)*h**2/12. )
        IF(ABS(s(i+1,l))>=Big) s(1:i+1,l) = s(1:i+1,l)*eps
      END DO
    END DO

    betap = ( s(ns-4,:) -8.*s(ns-3,:) +8.*s(ns-1,:) -s(ns,:) ) / ( 12.*h*km )
    beta = s(ns-2,:)
    delta = ATAN( (-jl*betap+jpl*beta) /(gl*betap-gpl*beta ) )
    fn = (COS(delta)*jl+SIN(delta)*gl) /beta

    DO l=0,lmax
      s(:,l) = s(:,l)*fn(l)
    END DO

  END SUBROUTINE ode_second_dw

  ELEMENTAL REAL(KIND=RP) FUNCTION Uij(ni, ei, nj, ej, r)
    USE special_functions ,ONLY: fac
    REAL(KIND=RP), INTENT(IN) :: ei, ej, r
    INTEGER      , INTENT(IN) :: ni, nj
    INTEGER :: k
    REAL(KIND=RP) :: a

    a = 1._RP
    DO k = 1,ni+nj-1
      a = a + ((ei+ej)*r)**k *(ni+nj-k) /( (ni+nj)*fac(k) )
    END DO
    Uij = a *EXP(-(ei+ej)*r) *fac(ni+nj) *ei**ni *ej**nj *SQRT(ei*ej) &
      /( SQRT( fac(2*ni) *fac(2*nj) ) *( (ei+ej)/2. )**(ni+nj+1) )

  END FUNCTION Uij

  SUBROUTINE calculate_U(Atom, Orbit, r, U, state )
    USE input ,ONLY: read_orbit
    CHARACTER(LEN=2), INTENT(IN) :: Atom, Orbit
    REAL(KIND=RP)   , INTENT(IN) :: r(:)
    REAL(KIND=RP)   , INTENT(OUT) :: U(:)
    INTEGER :: state

    INTEGER :: IN

    REAL(KIND=RP), ALLOCATABLE :: a(:),e(:)
    INTEGER, ALLOCATABLE :: n(:)
    INTEGER :: nelec, lo, no, i1, i2, nocup
    CHARACTER(LEN=2) :: orbit_i

    OPEN( newunit=IN, FILE='Data/'//Atom//'.dat', STATUS='old', ACTION='read')

    U = 0._RP
    DO
      READ(IN, FMT=*, IOSTAT=lo ) orbit_i
      IF(lo<0) EXIT

      CALL read_orbit(Atom//'_'//orbit_i, nelec, lo, no, n, a, e)
      nocup = nelec*(2*lo+1)
      IF(orbit_i==Orbit ) nocup = nocup -state

      IF(nocup==0) CYCLE

      DO i1 = 1,no
        U = U + nocup *a(i1)**2 *Uij( n(i1), e(i1), n(i1), e(i1), r )
        if(i1==1) cycle
        DO i2 = 1,i1-1
          U = U + 2*nocup *a(i1) *a(i2) *Uij( n(i1), e(i1), n(i2), e(i2), r )
        END DO
      END DO

    END DO

    CLOSE(IN)

  END SUBROUTINE calculate_U

  PURE SUBROUTINE INTRPL(X, Y, U, V )
    USE CONSTANTS ,ONLY: RP
    IMPLICIT NONE
    REAL(RP), INTENT(IN) :: X(:), Y(:), U(:)
    REAL(RP), INTENT(OUT) :: V(:)
    !
    !  REAL(WP) INTERPOLATION OF A SINGLE VALUED FUNCTION
    !  THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
    !  GIVEN  AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
    !  AND FOR A GIVEN SET OF X VALUES(ABSCISSAE),THE VALUES OF
    !  A SINGLE VALUED FUNCTION Y=Y(X).
    !
    !  THE INPUT PARAMETERS ARE:
    !
    !  L = NUMBER OF DATA POINTS (MUST BE TWO OR GREATER)
    !  X = ARRAY OF DIMENSION L STORING THE X VALUES OF INPUT DATA POINTS (IN ASCENDING ORDER)
    !  Y = ARRAY OF DIMENSION L STORING THE Y VALUES OF INPUT DATA POINTS
    !  N = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE Y-VALUES IS REQUIRED (MUST BE 1 OR GREATER)
    !  U = ARRAY OF DIMENSION N STORING THE X VALUES OF THE DESIRED POINTS
    !
    !  THE OUTPUT PARAMETER IS:
    !
    !  V = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y VALUES ARE TO BE DISPLAYED
    !
    !  DECLARATION STATEMENTS
    REAL(RP) :: A2, A3, A4, T4, TM2, TM3, TM4, X4, Y4
    INTEGER :: I, IMN, IMX, IPV, J, K, L, N
    REAL(RP), TARGET :: P0, Q0, Q1, UK, X2, X5, SW, Y2, Y5
    REAL(RP), POINTER :: DX, X3, Y3, T3, A1, TM1, A5, TM5, SA, W2, W4, Q2, W3, Q3

    X3 => P0
    P0 = 0._RP
    Y3 => Q0
    Q0 = 0._RP
    T3 => Q1
    Q1 = 0._RP
    DX => UK
    A1 => X2; TM1 => X2
    X2 = 0._RP
    A5 => X5; TM5 => X5
    X5 = 0._RP
    SA => SW
    W2 => Y2; W4 => Y2; Q2 => Y2
    Y2 = 0._RP
    W3 => Y5; Q3 => Y5
    Y5 = 0._RP


    !  PRELIMINARY PROCESSING
    L = SIZE(X)
    !IF(SIZE(Y)/=L) ERROR STOP 'size(Y)/=size(X)'
    N = SIZE(U)
    !IF(SIZE(V)/=N) ERROR STOP 'INTRPL : size(V)/=size(U)'

    TM4 = 0._RP
    A4 = 0._RP
    A2 = 0._RP

    IPV=0
    !  MAIN LOOP
    DO K = 1,N
      UK = U(K)
      !  ROUTINE TO LOCATE THE DESIRED POINT
      IF(UK<X(1)) THEN
        I = 1
      ELSEIF(UK>=X(L)) THEN
        I = L+1
      ELSE
        IMN = 2
        IMX = L
        DO
          I = (IMN+IMX)/2
          IF(UK<X(I)) THEN
            IMX = I
          ELSE
            IMN = I+1
          END IF
          IF(IMX<=IMN) EXIT
        END DO
        I = IMX
      END IF
      !  CHECK IF I=IPV
      IF(I/=IPV) THEN
        IPV = I
        !  ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND TO ESTIMATE THEM IF NECESSARY
        J = I
        IF(J==1) J = 2
        IF(J==L+1) J = L
        X3 = X(J-1)
        Y3 = Y(J-1)
        X4 = X(J)
        Y4 = Y(J)
        A3 = X4-X3
        TM3 = (Y4-Y3)/A3
        TM2 = 0._RP
        IF( L/=2 ) THEN
          IF( J/=2 ) THEN
            X2 = X(J-2)
            Y2 = Y(J-2)
            A2 = X3-X2
            TM2 = (Y3-Y2)/A2
          END IF
          IF(J/=L) THEN
            X5 = X(J+1)
            Y5 = Y(J+1)
            A4 = X5-X4
            TM4 = (Y5-Y4)/A4
            IF(J==2) TM2 = TM3+TM3-TM4
          ELSE
            TM4 = TM3+TM3-TM2
          END IF
        ELSE
          TM2 = TM3
        END IF
        IF(J>3) THEN
          A1 = X2-X(J-3)
          TM1 = (Y2-Y(J-3))/A1
        ELSE
          TM1 = TM2+TM2-TM3
        END IF
        IF(J<L-1) THEN
          A5 = X(J+2)-X5
          TM5 = (Y(J+2)-Y5)/A5
        ELSE
          TM5 = TM4+TM4-TM3
        END IF
        !  NUMERICAL DIFFERENTIATION
        IF(I/=L+1) THEN
          W2 = ABS(TM4-TM3)
          W3 = ABS(TM2-TM1)
          SW = W2+W3
          IF(SW==0._RP) THEN
            W2 = 0.5_RP
            W3 = 0.5_RP
            SW = 1._RP
          END IF
          T3 = (W2*TM2+W3*TM3)/SW
        END IF
        IF(I==1) THEN
          T4 = T3
          SA = A3+A4
          T3 = 0.5_RP*(TM1+TM2-A4*(A3-A4)*(TM3-TM4)/(SA*SA))
          X3 = X3-A4
          Y3 = Y3-TM2*A4
          A3 = A4
          TM3 = TM2
        ELSE
          W3 = ABS(TM5-TM4)
          W4 = ABS(TM3-TM2)
          SW = W3+W4
          IF(SW==0._RP) THEN
            W3 = 0.5_RP
            W4 = 0.5_RP
            SW = 1._RP
          END IF
          T4 = (W3*TM3+W4*TM4)/SW
          IF(I==L+1) THEN
            T3 = T4
            SA = A2+A3
            T4 = 0.5_RP*(TM4+TM5-A2*(A2-A3)*(TM2-TM3)/(SA*SA))
            X3 = X4
            Y3 = Y4
            A3 = A2
            TM3 = TM4
          END IF
        END IF
        !  COMPUTATION OF THE POLYNOMIAL
        Q2 = (2._RP*(TM3-T3)+TM3-T4)/A3
        Q3 = (-TM3-TM3+T3+T4)/(A3*A3)
      END IF
      DX = UK-P0
      V(K) = Q0+DX*(Q1+DX*(Q2+DX*Q3))

    END DO

    RETURN
  END SUBROUTINE INTRPL

END MODULE utils
