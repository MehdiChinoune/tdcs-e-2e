SUBMODULE(special_functions) special_functions
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE factorial()
    INTEGER :: i

    IF(fac_called) RETURN

    fac(0) = 1._RP
    DO i = 1,34
      fac(i) = i*fac(i-1)
    END DO

    lnfac(0:34) = LOG( fac )
    DO i = 35, 400
      lnfac(i) = lnfac(i-1) +LOG( REAL(i,KIND=RP) )
    END DO

    fac_called = .TRUE.

  END SUBROUTINE

  !-----------------------------------------------------------------------
  !        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS
  !                        ---------------
  !     MO IS AN INTEGER, Z A COMPLEX ARGUMENT, AND W A COMPLEX VARIABLE.
  !                 W = GAMMA(Z)       IF MO = 0
  !                 W = LN(GAMMA(Z))   OTHERWISE
  !-----------------------------------------------------------------------
  !     WRITTEN BY ALFRED H. MORRIS, JR.
  !        NAVAL SURFACE WARFARE CENTER
  !        DAHLGREN, VIRGINIA
  !     This version, in a subset of Fortran 90, prepared by
  !     Alan.Miller @ vic.cmis.csiro.au
  !     http://www.ozemail.com.au/~milleraj
  !
  !     This version is accurate to within 5 in the 14th significant
  !     decimal digit.
  !-----------------------------------------------------------------------
  ELEMENTAL MODULE FUNCTION cgamma(z, mo) RESULT(w)
    COMPLEX(RP) :: w
    COMPLEX(RP), INTENT(IN)  :: z
    INTEGER, INTENT(IN) :: mo

    ! Local variables
    COMPLEX(RP) :: eta, eta2, sum_
    REAL(RP), PARAMETER :: c0(12) = [ .833333333333333E-01_RP, -.277777777777778E-02_RP &
      , .793650793650794E-03_RP, -.595238095238095E-03_RP,  .841750841750842E-03_RP &
      , -.191752691752692E-02_RP,  .641025641025641E-02_RP, -.295506535947712E-01_RP &
      , .179644372368831_RP    , -1.39243221690590_RP    , 13.4028640441684_RP      &
      , -156.848284626002_RP ], pi = 3.14159265358979_RP, pi2  = 6.28318530717959_RP &
      , alpi = 1.14472988584940_RP, hl2p = .918938533204673_RP, half = 0.5_RP
    REAL(RP) :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
      s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
    INTEGER :: j, k, l, m, n, nm1
    !---------------------------
    !     ALPI = LOG(PI)
    !     HL2P = 0.5*LOG(2*PI)
    !---------------------------

    !     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
    !            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
    !            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
    !            SUCH THAT 1.0 +EPS>1.0.

    !                      MAX = IPMPAR(3)
    !MAX = HUGE(1)
    eps = EPSILON(1._RP)
    !---------------------------
    x = REAL(z, KIND=RP)
    y = AIMAG(z)
    IF( y==0._RP ) THEN
      IF( mo==0 ) THEN
        w = CMPLX( gamma(x), 0._RP, RP )
      ELSE
        w = CMPLX( log_gamma(x), 0._RP, RP )
      END IF
      RETURN
    END IF
    !IF( ABS(x)>=MIN( REAL(MAX,RP), 1._RP/eps ) ) THEN! GO TO 70
    w = CMPLX(0.,0.,RP)

    s2 = 0._RP
    s1 = 0._RP
    s = 0._RP
    h2 = 0._RP
    h1 = 0._RP
    ! RETURN
    !END IF
    IF( x<0._RP ) THEN
      ! CASE WHEN THE REAL PART OF Z IS NEGATIVE
      y = ABS(y)
      t = -pi*y
      et = EXP(t)
      e2t = et*et

      ! SET A1 = (1 +E2T)/2  AND  A2 = (1 -E2T)/2
      a1 = half*(1._RP +e2t)
      t2 = t +t
      IF(t2>=-0.15_RP) THEN
        a2 = -half*rexp(t2)
      ELSE
        a2 = half*(half +(half -e2t))
      END IF

      !  COMPUTE SIN(PI*X) AND COS(PI*X)
      !IF(ABS(x)>=MIN(REAL(MAX,RP), 1._RP/eps)) GO TO 70
      k = INT( ABS(x) )
      u = x +k
      k = MOD(k,2)
      IF(u <= -half) THEN
        u = half +(half +u)
        k = k +1
      END IF
      u = pi*u
      sn = SIN(u)
      cn = COS(u)
      IF(k==1) THEN
        sn = -sn
        cn = -cn
      END IF
      !  SET  H1 +H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))
      a1 = sn*a1
      a2 = cn*a2
      a = a1*a1 +a2*a2
      IF( a==0._RP ) RETURN!GO TO 70
      IF( mo==0 ) THEN
        h1 = a1 / a
        h2 = -a2 / a
        c = pi*et
        h1 = c*h1
        h2 = c*h2
      ELSE
        h1 = (alpi+t) -half*LOG(a)
        h2 = -ATAN2(a2,a1)
      END IF
      IF( AIMAG(z)>=0._RP ) THEN
        x = 1._RP -x
        y = -y
      ELSE
        h2 = -h2
        x = 1._RP -x
      END IF
    END IF
    !-----------------------------------------------------------------------
    !           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
    !-----------------------------------------------------------------------
    w1 = 0._RP
    w2 = 0._RP
    n = 0
    t = x
    y2 = y*y
    a = t*t +y2
    cut = 36._RP
    IF( eps>1.E-8_RP ) cut = 16._RP
    IF( a<cut ) THEN
      IF( a==0._RP) RETURN
      DO WHILE( a<cut )
        n = n +1
        t = t +1._RP
        a = t*t +y2
      END DO
      !     LET S1 +S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)
      u1 = (x*t+y2) / a
      u2 = y / a
      s1 = u1
      s2 = n*u2
      IF( n>=2 ) THEN
        u = t / a
        nm1 = n -1
        DO j = 1, nm1
          v1 = u1 +j*u
          v2 = (n-j)*u2
          c = s1*v1 -s2*v2
          d = s1*v2 +s2*v1
          s1 = c
          s2 = d
        END DO
      END IF

      !     SET  W1 +W2*I = LOG(S1 +S2*I)  WHEN MO IS NONZERO

      s = s1*s1 +s2*s2
      IF( mo/=0 ) THEN
        w1 = half*LOG(s)
        w2 = ATAN2(s2,s1)
      END IF
    END IF

    !     SET  V1 +V2*I = (Z -0.5)*LOG(Z +N) -Z

    t1 = half*LOG(a) -1._RP
    t2 = ATAN2(y,t)
    u = x -half
    v1 = (u*t1-half) -y*t2
    v2 = u*t2 +y*t1

    !     LET A1 +A2*I BE THE ASYMPTOTIC SUM

    eta = CMPLX(t/a, -y/a, KIND=RP)
    eta2 = eta*eta
    m = 12
    IF( a>=289._RP ) m = 6
    IF( eps>1.e-8_RP ) m = m / 2
    sum_ = CMPLX(c0(m), 0._RP, KIND=RP)
    l = m
    DO j = 2, m
      l = l -1
      sum_ = CMPLX(c0(l), 0._RP, KIND=RP) +sum_*eta2
    END DO
    sum_ = sum_*eta
    a1 = REAL(sum_, KIND=RP)
    a2 = AIMAG(sum_)
    !-----------------------------------------------------------------------
    !                 GATHERING TOGETHER THE RESULTS
    !-----------------------------------------------------------------------
    w1 = (((a1 +hl2p) -w1) +v1) -n
    w2 = (a2 -w2) +v2
    IF( x>=0._RP .AND. mo==0 ) THEN ! CASE WHEN Re(Z)>0 AND MO==0
      a = EXP(w1)
      w1 = a*COS(w2)
      w2 = a*SIN(w2)
      IF( n/=0 ) THEN
        c = (s1*w1 +s2*w2) / s
        d = (s1*w2 -s2*w1) / s
        w1 = c
        w2 = d
      END IF
    ELSEIF(M/=0) THEN ! CASE WHEN MO/=0
      IF( x<0._RP ) THEN ! CASE WHEN Re(Z)<0 AND M/=0
        w1 = h1 -w1
        w2 = h2 -w2
      END IF
      ! CASE WHEN Re(Z)>0 AND MO/=0
      ! THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI<W2 <= PI.
      IF( w2<=pi ) THEN
        k = INT( half -w2 / pi2 )
        w2 = w2 +pi2*k
      ElSE
        k = INT ( w2 / pi2 -half )
        w2 = w2 -pi2*REAL(k+1,KIND=RP)
        IF( w2<=-pi ) w2 = pi
      END IF
    ELSE ! CASE WHEN Re(Z)<0 AND MO==0
      a = EXP(-w1)
      t1 = a*COS(-w2)
      t2 = a*SIN(-w2)
      w1 = h1*t1 -h2*t2
      w2 = h1*t2 +h2*t1
      IF( n/=0 ) THEN
        c = w1*s1 -w2*s2
        d = w1*s2 +w2*s1
        w1 = c
        w2 = d
      END IF
    END IF

    w = CMPLX(w1, w2, KIND=RP)

    RETURN

  CONTAINS
    !-----------------------------------------------------------------------
    !            EVALUATION OF THE FUNCTION EXP(X) -1
    !-----------------------------------------------------------------------
    ELEMENTAL REAL(RP) FUNCTION rexp(x)
      REAL(RP), INTENT(IN) :: x
      ! Local variables
      REAL(RP), PARAMETER  :: p1 = .914041914819518E-9_RP, p2 = .238082361044469E-1_RP &
        , q1 = -.499999999085958_RP, q2 = .107141568980644_RP &
        , q3 = -.119041179760821E-1_RP, q4 = 0.595130811860248E-3_RP
      REAL(RP) :: e
      !-----------------------
      IF( ABS(x)<=0.15_RP ) THEN
        rexp = x*(((p2*x +p1)*x +1._RP) / ( ( ( (q4*x +q3)*x +q2)*x +q1 )*x +1._RP) )
      ELSEIF( x>=0._RP ) THEN
        e = EXP(x)
        rexp = e*(half +(half -1._RP/e))
      ELSEIF( x>=-37._RP ) THEN
        rexp = (EXP(x) -half) -half
      ELSE
        rexp = -1._RP
      END IF

      RETURN
    END FUNCTION rexp

  END FUNCTION cgamma

  !> assoc_legendre
  !! Computes the associated Legendre polynomial P_l^m (x). Here m and l are integers
  !! satisfying 0 <= m <= l,  while x lies in the range −1 <= x <= 1.
  ELEMENTAL REAL(RP) MODULE FUNCTION assoc_legendre(l,m,x)
    USE ieee_arithmetic ,only: ieee_is_finite
    INTEGER, INTENT(IN) :: l, m
    REAL(RP), INTENT(IN) :: x

    INTEGER :: i, l1
    REAL(RP) :: fact, pmm, pmm1, pmm2

    !IF( m<0 .OR. m>l ) ERROR STOP 'M is out of range [0,L]'
    !IF( ABS(x)>1._RP ) ERROR STOP 'X is out of range [-1., 1.]'

    !Compute P_m^m
    pmm = 1._RP
    IF( m>0 ) THEN
      fact = 1._RP
      DO i = 1, m
        pmm = pmm *fact
        fact = fact +2._RP
      END DO
      pmm = pmm *( -SQRT(1._RP-x**2) )**m
    ENDIF

    IF( m==l ) THEN ! Compute P_l^l
      assoc_legendre = pmm
    ELSEIF( m==l-1 ) THEN  ! Compute P_l^{l-1}
      assoc_legendre = x*(2*m+1) *pmm
    ELSE ! Compute P_l^m, m<l-1
      pmm1 = 0._RP
      DO l1 = m+1, l
        pmm2 = pmm1
        pmm1 = pmm
        pmm = ( x*(2*l1-1)*pmm1 -(l1+m-1)*pmm2 ) /(l1-m)
      END DO
      assoc_legendre = pmm
    ENDIF
    IF( .NOT. ieee_is_finite(pmm) ) assoc_legendre = huge(1._RP)

    RETURN
  END FUNCTION assoc_legendre

  ELEMENTAL COMPLEX(RP) MODULE FUNCTION spherical_harmonic( l, m, theta, phi )
    use constants ,only: pi
    INTEGER, INTENT(IN) :: l, m
    REAL(RP), INTENT(IN) :: theta, phi
    REAL(RP), PARAMETER :: Tinye = log(tiny(1._RP))
    INTEGER :: ma

    !if(.not. fac_called ) error stop 'you should call factorial before using &
    !  &y_l^m(\theta,\phi)'

    ma = ABS(m)
    IF( (lnfac(l-ma)-lnfac(l+ma))<2*Tinye ) THEN
      spherical_harmonic = (0._RP, 0._RP)
      RETURN
    ELSEIF( l==0 .AND. m==0 ) THEN
      spherical_harmonic = 1._RP/SQRT(4._RP*pi)
      RETURN
    ELSEIF( m==0 ) THEN
      spherical_harmonic = SQRT( (2.*l+1._RP)/(4.*pi) ) *assoc_legendre(l, 0, COS(theta) )
      RETURN
    END IF

    spherical_harmonic = EXP( 0.5*( lnfac(l-ma) -lnfac(l+ma) ) ) &
      *SQRT( (2.*l+1._RP)/(4.*pi) ) *CMPLX( COS(m*phi), SIN(m*phi), KIND=RP ) &
      *assoc_legendre(l, ma, COS(theta) )
    IF( MOD(m,2)>0 ) spherical_harmonic = -spherical_harmonic
    IF( MODULO(theta,2*pi)>pi .AND. MOD(m,2)/=0 ) spherical_harmonic = -spherical_harmonic

  END FUNCTION spherical_harmonic

  !>  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
  !!
  !!  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
  !!   FOR REAL X>0., REAL ETA (INCLUDING 0.), AND REAL XLMIN >-1.
  !!   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
  !!   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
  !!   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
  !!   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
  !!   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
  !!
  !!   CALLING VARIABLES; ALL REALS ARE REAL (REAL*8)
  !!
  !! @param[in]  X     -REAL ARGUMENT FOR COULOMB FUNCTIONS>0.0 [ X>SQRT(ACCUR) :
  !!                     ACCUR IS TARGET ACCURACY 1.0D-14 ]
  !! @param[out] ETA   -REAL SOMMERFELD PARAMETER, UNRESTRICTED>=< 0.0
  !! @param[in] XLMIN  -REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER), GENERALLY IN RANGE
  !!                      0.0 -1.0 AND MOST  USUALLY 0.0
  !! @param[in] LRANGE -INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED FOR
  !!                      L-VALUES XLMIN TO XLMIN +LRANGE INCLUSIVE
  !! @param[out] FC    -REAL VECTORS F OF REGULAR COULOMB FUNCTIONS
  !! @param[out] GC    -REAL VECTORS G OF IRREGULAR COULOMB FUNCTIONS
  !! @param[out] FCP   -REAL VECTOR FOR THE X-DERIVATIVES OF  F THIS VECTOR TO BE OF
  !!                      LENGTH AT LEAST IDUM3 +LRANGE STARTING ELEMENT
  !!                      IDUM3 = MAX( INT(XLMIN+ACCUR),0 )
  !! @param[out] GCP   -REAL VECTOR FOR THE X-DERIVATIVES OF  G THIS VECTOR TO BE OF
  !!                      LENGTH AT LEAST IDUM3 +LRANGE STARTING ELEMENT
  !!                      IDUM3 = MAX( INT(XLMIN+ACCUR),0 )
  !! @param[in]  KFN   -INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
  !!                     = 0      REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
  !!                     = 1    SPHERICAL BESSEL      "      "     "     J & Y
  !!                     = 2  CYLINDRICAL BESSEL      "      "     "     J & Y
  !! @param[in]   IFAIL ON INPUT IS SET TO 0 (LIMIT = 20000)
  !! @param[out]  IFAIL IN OUTPUT
  !!                =  0 : CALCULATIONS SATISFACTORY
  !!                =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
  !!                =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
  !!                = -1 : X<1D-7 = SQRT(ACCUR)
  !!                = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES)
  !!
  !!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
  !!   IN OSCILLATING REGION X .GE. [ETA +SQRT{ETA**2 +XLM*(XLM+1)}]
  !!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100*ACC8 WHERE ACC8 IS
  !!   THE SMALLEST NUMBER WITH 1.+ACC8/=1. FOR OUR WORKING PRECISION.
  !!   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
  !!   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
  !!   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADLY WORSE :
  !!   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
  !!----------------------------------------------------------------------
  !!   ERROR RETURNS      THE USER SHOULD TEST IFAIL ON EXIT
  !!----------------------------------------------------------------------
  !!  MACHINE-DEPENDENT PARAMETERS:    ACCUR -SEE ABOVE
  !!      SMALL -OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
  !!      IE 1D-30 FOR IBM REAL*8,    1E-150 FOR REAL
  !!----------------------------------------------------------------------
  !!  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
  !!  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
  !!                 + RCWFF       IN    CPC 11 (1976) 141-142
  !!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
  !!  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
  !!  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188
  !!  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
  !!  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
  !!  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
  !!  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509
  !!----------------------------------------------------------------------
  !!  AUTHOR: A. R. BARNETT      MANCHESTER  MARCH   1981
  !!                             AUCKLAND    MARCH   1991
  !!----------------------------------------------------------------------
  PURE MODULE SUBROUTINE coul90(x, eta, lmin, lrange, fc, gc, fcp, gcp, kfn, ifail )
    INTEGER, INTENT(IN)  :: lmin, lrange, kfn
    INTEGER, INTENT(OUT), OPTIONAL :: ifail
    REAL(RP), INTENT(IN) :: x, eta
    REAL(RP), INTENT(OUT), DIMENSION(lmin:lmin+lrange) :: fc,  gc,  fcp, gcp

    INTEGER, PARAMETER :: limit = 20000
    REAL(RP), PARAMETER :: small = SQRT(TINY(1._RP)), ten2 = 100._RP, half = 0.5_RP &
      , rt2dpi = 0.797884560802865_RP
    !!---- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM IDUM3
    REAL(RP) :: accur, acch, xinv, pk, cf1, c, d, pk1, etak, rk2, tk, dcf1, den, xlm &
      , xll, el, xl, rl, sl, f, fcmaxl, fcminl, gcminl, omega, wi, a, b, ar, ai, br &
      , bi, dr, di, dp, dq, alpha, beta, e2mm1, fjwkb,  gjwkb, p, q, gamma_, gammai, ERR
    INTEGER :: l, maxl, idum2, idum3
    LOGICAL :: etane0, xlturn
    !----------------------------------------------------------------------
    !     COUL90 HAS CALLS TO: SQRT,ABS,MAX,INT,SIGN,DBLE,MIN
    !----------------------------------------------------------------------
    !Q    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
    !---- THIS CONSTANT IS  SQRT(2._RP / PI):
    !---- USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND REAL
    !---- CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED

    accur = EPSILON(1._RP)
    !    accur = 1.E-7_RP
    IF(present(ifail)) ifail = 0
    idum2 = 1
    !idum1 = 0
    gjwkb = 0._RP
    ERR = 1._RP
    ! IF(KFN/=0) ETA = 0._RP
    etane0 = ( eta/=0._RP .and. kfn==0 )
    acch = SQRT(accur)
    gammai = 0._RP
    !!---- TEST RANGE OF X, EXIT IF<=SQRT(ACCUR) OR IF NEGATIVE
    IF( x <= acch .and. present(ifail) ) THEN
      ifail = -1
      !WRITE(6,1000) X,ACCH
      !1000 FORMAT(' FOR X = ',1P,D12.3,'     TRY SMALL-X  SOLUTIONS' &
      ! ,' OR X IS NEGATIVE'/, ' SQUARE ROOT (ACCURACY) =  ',D12.3/)
      RETURN
    ENDIF

    IF( kfn==2 ) THEN
      xlm = lmin -half
    ELSE
      xlm = lmin
    ENDIF
    IF( xlm <= -1._RP .OR. lrange<0 .and. present(ifail) ) THEN
      ifail = -2
      !  WRITE (6,1001) LRANGE,LMIN,XLM
      !1001 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',2I10,1P,D15.6/)
      RETURN
    ENDIF
    e2mm1 = xlm*xlm +xlm
    xlturn = x*(x -2._RP*eta)<e2mm1
    e2mm1 = e2mm1 +eta*eta
    xll = xlm +REAL(lrange,KIND=RP)
    !---- LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
    !---- XLL    IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
    !---- DETERMINE STARTING ARRAY ELEMENT (IDUM3) FROM XLMIN
    idum3 = lmin!MAX( INT(XLMIN +ACCUR),0 )
    maxl = idum3 +lrange
    !---- EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
    xinv = 1._RP / x
    !---- UNNORMALISED F(MAXL,ETA,X)
    den = 1._RP
    pk = xll +1._RP
    cf1 = eta / pk +pk*xinv
    IF( ABS(cf1)<small ) cf1 = small
    rk2 = 1._RP
    d = 0._RP
    c = cf1
    !---- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA +1: LENTZ-THOMPSON
    DO l = 1,  limit
      pk1 = pk +1._RP
      IF( etane0 ) THEN
        etak = eta / pk
        rk2 = 1._RP +etak*etak
        tk = (pk +pk1)*(xinv +etak / pk1)
      ELSE
        tk = (pk +pk1)*xinv
      ENDIF
      !---- DIRECT RATIO OF B CONVERGENTS
      d = tk -rk2*d
      !---- INVERSE RATIO OF A CONVERGENTS
      c = tk -rk2 / c
      IF( ABS(c)<small ) c = small
      IF( ABS(d)<small ) d = small
      d = 1._RP / d
      dcf1 =  d*c
      cf1 = cf1*dcf1
      IF( d<0._RP ) den = -den
      pk = pk1
      !----  PROPER EXIT
      IF( ABS(dcf1-1._RP)<accur ) EXIT
    END DO
    !---- ERROR EXIT
    IF( ABS(dcf1-1._RP)>=accur .and. present(ifail) ) THEN
      ifail = 1
      !WRITE (6,1002) LIMIT, CF1,DCF1, PK,ACCUR
      !1002  FORMAT('CF1 HAS FAILED TO CONVERGE AFTER',I10,'ITERATIONS',/'CF1,DCF1,PK
      !,ACCUR = ',1P,4D12.3/)
      RETURN
    ENDIF

    !nfp = INT(pk -xll -1)
    f = cf1
    !---- DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
    IF( lrange>0 ) THEN
      fcmaxl = small*den
      fcp(maxl) = fcmaxl*cf1
      fc (maxl) = fcmaxl
      xl = xll
      rl = 1._RP
      DO l = maxl, idum3+1, -1
        IF( etane0 ) THEN
          el = eta / xl
          rl = SQRT( 1._RP +el*el )
          sl = xl*xinv +el
          gc (l) = rl
          gcp(l) = sl
        ELSE
          sl = xl*xinv
        ENDIF
        fc (l-1) = ( fc(l)*sl +fcp(l) ) / rl
        fcp(l-1) = fc(l-1)*sl -fc (l)*rl
        !---- END VALUE IS XLM
        xl = xl -1._RP
      END DO
      IF( ABS(fc(idum3))<accur*small ) fc(idum3) = accur*small
      f = fcp(idum3) / fc(idum3)
      den = fc (idum3)
    ENDIF
    !---------------------------------------------------------------------
    !---- NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
    !---- EVALUATE CF2 = P +I.Q  USING STEED'S ALGORITHM (NO ZEROS)
    !---------------------------------------------------------------------
    IF( xlturn ) CALL jwkb(x,eta,MAX(xlm,0._RP),fjwkb,gjwkb,idum2)
    IF( idum2>1 .OR. gjwkb>1._RP/(acch*ten2) ) THEN
      omega = fjwkb
      gamma_ = gjwkb*omega
      p = f
      q = 1._RP
    ELSE
      xlturn = .false.
      pk = 0._RP
      wi = eta +eta
      p = 0._RP
      q = 1._RP -eta*xinv
      ar = -e2mm1
      ai = eta
      br = 2._RP*(x -eta)
      bi = 2._RP
      dr = br / (br*br +bi*bi)
      di = -bi / (br*br +bi*bi)
      dp = -xinv*(ar*di +ai*dr)
      dq = xinv*(ar*dr -ai*di)
      DO l = 1, limit
        p = p +dp
        q = q +dq
        pk = pk +2._RP
        ar = ar +pk
        ai = ai +wi
        bi = bi +2._RP
        d = ar*dr -ai*di +br
        di = ai*dr +ar*di +bi
        c = 1._RP / (d*d +di*di)
        dr = c*d
        di = -c*di
        a = br*dr -bi*di -1._RP
        b = bi*dr +br*di
        c = dp*a -dq*b
        dq = dp*b +dq*a
        dp = c
        IF( ABS(dp)+ABS(dq)<(ABS(p)+ABS(q))*accur ) EXIT
      END DO
      !---- ERROR EXIT
      IF( ABS(dp)+ABS(dq)>=(ABS(p)+ABS(q))*accur .and. present(ifail) ) THEN
        ifail = 2
        !WRITE (6,1003) LIMIT,P,Q,DP,DQ,ACCUR
        !1003    FORMAT('CF2 HAS FAILED TO CONVERGE AFTER',I7,'ITERATIONS',&
          !               &/'P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3/)
        RETURN
      ENDIF

      !idum1 = INT(pk / 2._RP)
      ERR = half*accur / MIN( ABS(q),1._RP )
      IF( ABS(p)>ABS(q) ) ERR = ERR*ABS(p)
      !---------------------------------------------------------------------
      !    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
      !---------------------------------------------------------------------
      gamma_ = (f -p)/q
      gammai = 1._RP/gamma_
      IF( ABS(gamma_) <= 1._RP ) THEN
        omega = SQRT( 1._RP +gamma_*gamma_ )
      ELSE
        omega = SQRT( 1._RP +gammai* gammai)*ABS(gamma_)
      ENDIF
      omega = 1._RP / ( omega*SQRT(q) )
      ! wronsk = omega
    ENDIF
    !---------------------------------------------------------------------
    !---- RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
    !---------------------------------------------------------------------
    IF(kfn==1) THEN !---- SPHERICAL
      alpha = xinv
      beta = xinv
    ELSE IF(kfn==2) THEN !----   CYLINDRICAL
      alpha = half*xinv
      beta = SQRT( xinv )*rt2dpi
    ELSE !----   KFN = 0, COULOMB FUNCTIONS
      alpha = 0._RP
      beta = 1._RP
    ENDIF
    fcminl = SIGN( omega,den )*beta
    IF( xlturn ) THEN
      gcminl = gjwkb*beta
    ELSE
      gcminl = fcminl*gamma_
    ENDIF
    !---- BESSEL SIGN
    IF( kfn/=0 ) gcminl = -gcminl
    fc (idum3) = fcminl
    gc (idum3) = gcminl
    gcp(idum3) = gcminl*(p -q*gammai -alpha)
    fcp(idum3) = fcminl*(f -alpha)
    IF( lrange==0 ) RETURN
    !---------------------------------------------------------------------
    !---- UPWARD RECURRENCE FROM GC(IDUM3),GCP(IDUM3) STORED VALUES ARE RL,SL
    !---- RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
    !---- XL   = XLM HERE  AND RL = 1._RP,  EL = 0._RP FOR BESSELS
    !---------------------------------------------------------------------
    omega = beta*omega / ABS(den)
    xl = xlm
    rl = 1._RP
    DO l = idum3+1,  maxl
      xl = xl +1._RP
      IF( etane0 ) THEN
        rl = gc (l)
        sl = gcp(l)
      ELSE
        sl = xl*xinv
      ENDIF
      gc (l) = ( (sl -alpha)*gc(l-1) -gcp(l-1) ) / rl
      gcp(l) = rl*gc(l-1) -(sl +alpha)*gc(l)
      fcp(l) = omega*( fcp(l) -alpha*fc(l) )
      fc (l) = omega*fc (l)
    END DO

    if(present(ifail) ) ifail = 0

    RETURN
  END SUBROUTINE coul90
  !----------------------------------------------------------------------
  !> COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
  !! AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
  !! CALCULATED IN SINGLE, RETURNED IN REAL VARIABLES
  !! CALLS MAX, SQRT, LOG, EXP, ATAN2, REAL, INT
  !!     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
  !----------------------------------------------------------------------
  ELEMENTAL SUBROUTINE  jwkb( x, eta, xl, fjwkb, gjwkb, iexp )
    REAL(RP), INTENT(IN)  :: x, eta, xl
    REAL(RP), INTENT(OUT) :: fjwkb, gjwkb
    INTEGER, INTENT(OUT) :: iexp

    INTEGER, PARAMETER :: maxexp = 300
    REAL(RP), PARAMETER :: half = 0._RP, six = 6._RP, ten = 10._RP &
      , dzero = 0._RP, rl35 = 35._RP, aloge = 0.4342945_RP

    REAL(RP) :: phi, phi10, gh2, xll1, hll, hl, sl, rl2, gh
    !----------------------------------------------------------------------
    !---- CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR REAL
    !----------------------------------------------------------------------
    gh2   = x*(eta +eta -x)
    xll1  = MAX( xl*xl +xl, dzero )
    IF( gh2 +xll1 <= 0._RP ) RETURN
    hll = xll1 +six / rl35
    hl = SQRT(hll)
    sl = eta / hl +hl / x
    rl2 = 1._RP +eta*eta / hll
    gh = SQRT(gh2 +hll) / x
    phi = x*gh -half*( hl*LOG((gh +sl)**2 / rl2) -LOG(gh) )
    IF( eta/=0._RP ) phi = phi -eta*ATAN2(x*gh,x -eta)
    phi10 = -phi*aloge
    iexp = INT(phi10)
    IF( iexp>maxexp ) THEN
      gjwkb = ten**(phi10 -REAL(iexp,KIND=RP))
    ELSE
      gjwkb = EXP(-phi)
      iexp = 0
    ENDIF
    fjwkb = half / (gh*gjwkb)

    RETURN
  END SUBROUTINE  jwkb

  !>   REAL RICCATI-BESSEL FUNCTIONS AND X-DERIVATIVES :
  !!   PSI = X . J/L/(X),  CHI = X . Y/L/(X)    FROM L = 0 TO L = LMAX
  !!      FOR REAL X>SQRT(ACCUR) (E.G. 1D-7)  AND INTEGER LMAX
  !! PSI (L)  =      PSI/L/(X) STORES   REGULAR RICCATI-BESSEL FUNCTION:
  !! PSID(L)  = D/DX PSI/L/(X)          PSI(0) =  SIN(X)
  !! CHI (L)  =      CHI/L/(X) STORES IRREGULAR RICCATI-BESSEL FUNCTION:
  !! CHID(L)  = D/DX CHI/L/(X)          CHI(0) = -COS(X) [HANDBOOK DEFN.]
  !!
  !!    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
  !!          =  0 FOR ALL RESULTS SATISFACTORY
  !!
  !!   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
  !!   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
  !!   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
  !!   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
  !!   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
  !!         SO 1E-150 FOR REAL ON VAX, PCS ETC.
  !!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
  !!   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
  !!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100*ACC8 WHERE ACC8
  !!   IS THE SMALLEST NUMBER WITH 1+ACC8/=1 FOR OUR WORKING PRECISION
  !!   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
  !!   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
  !!   IF X IS SMALLER THAN [ ] ABOVE, THE ACCURACY BECOMES STEADLY WORSE:
  !!   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
  !!
  !!   NOTE: FOR X = 1 AND L = 100   PSI = 7.4 E-190     CHI = -6.7+E186
  !!---------------------------------------------------------------------
  !!   AUTHOR :   A.R.BARNETT      MANCHESTER    12 MARCH 1990.
  !!                               AUCKLAND      12 MARCH 1991.
  !!---------------------------------------------------------------------
  PURE MODULE SUBROUTINE ricbes( x, lmax, psi, chi, psid, chid, ifail )
    REAL(RP), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: lmax
    INTEGER, INTENT(INOUT), OPTIONAL :: ifail
    REAL(RP), INTENT(OUT) :: psi(0:lmax), chi(0:lmax), psid(0:lmax), chid(0:lmax)

    INTEGER, PARAMETER :: limit = 20000 !,maxl = 1001
    REAL(RP), PARAMETER :: small = SQRT(TINY(1._RP)) , three = 3._RP

    INTEGER :: l
    REAL(RP) :: accur, tk, sl!, ERR
    REAL(RP) :: xinv, cf1, dcf1, den, c, d, omega, twoxi
    !COMMON / STEED  / ERR,NFP,IDUM1,IDUM2,IDUM3
    !!----
    !!---- CALCULATE THE L = 0   RICCATI-BESSEL FUNCTIONS DIRECTLY
    psi (0) = SIN(x)
    chi (0) = -COS(x)
    psid(0) = -chi(0)
    chid(0) = psi(0)
    accur = EPSILON(1._RP)!1.E-17
    IF( present(ifail) ) ifail = -1
    IF( x<SQRT(accur) ) GOTO 50
    !!---- TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !!---- BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
    xinv = 1._RP / x
    IF( lmax>0 ) THEN
      twoxi = xinv +xinv
      sl = REAL(lmax,KIND=RP)* xinv
      tk = 2._RP*sl +xinv*three
      cf1 = sl +xinv
      den = 1._RP
      IF( ABS(cf1)<small ) cf1 = small
      !!---- INVERSE RATIO OF A CONVERGENTS
      c = cf1
      !!---- DIRECT  RATIO OF B CONVERGENTS
      d = 0._RP
      DO l = 1,limit
        c = tk -1._RP / c
        d = tk -d
        IF( ABS(c)<small ) c = small
        IF( ABS(d)<small ) d = small
        d = 1._RP / d
        dcf1 =  d*c
        cf1 = cf1*dcf1
        IF( d<0._RP ) den = -den
        IF( ABS(dcf1 -1._RP)<=accur ) EXIT
        tk = tk +twoxi
      END DO
      !!---- ERROR EXIT, NO CONVERGENCE
      IF( l>=limit ) GOTO 50
      ! 20 CONTINUE !   nfp = l
      !!---- ERROR ESTIMATE
      !ERR = accur*SQRT(REAL(nfp,KIND=RP))
      psi (lmax) = den
      psid(lmax) = cf1*den
      !!---- DOWNWARD RECURSION TO L = 0  AS RICCATI-BESSEL FUNCTIONS
      DO l = lmax,  2, -1
        psi (l-1) = sl*psi(l) +psid(l)
        psid(l-1) = sl*psi(l-1) -psi (l)
        sl = sl -xinv
      END DO
      den  = sl*psi(1) +psid(1)!PSI(0)
      !!---- END LOOP FOR LMAX>0
      !ENDIF
      !IF(LMAX>0) THEN
      omega  = psi(0) / den
      DO l = 1,  lmax
        psi (l) = omega*psi (l)
        psid(l) = omega*psid(l)
        sl = xinv*REAL(l,KIND=RP)
        chi (l) = sl*chi(l-1) -chid(l-1)
        chid(l) = chi(l-1) -sl*chi (l)
      END DO
    ENDIF
    !!---- CALCULATIONS SUCCESSFUL
    if(present(ifail)) ifail = 0
    RETURN
    !!---------------------------------------------------------------------
    !!---- ERROR TRAPS
    !!---------------------------------------------------------------------
    50 CONTINUE
    IF( x<0._RP ) THEN
      !    WRITE(6,1000) x
      !    1000   FORMAT(' X NEGATIVE !',1p,e15.5,' ... USE REFLECTION FORMULA'/)
    ELSE IF( x==0._RP ) THEN
      ifail = 0
      chi(0) = 1._RP
      DO l = 1, lmax
        chi(l) = 0._RP
      END DO
    ELSE
      ! WRITE(6,1001) x
      ! 1001   FORMAT(' WITH X = ',1p,e15.5,'    TRY SMALL-X SOLUTIONS',/ &
        !   , '  PSI_L(X)=>X**(L+1)/(2L+1) AND',/,'  CHI/L/(X)-> -(2L-1)/X**L'/)
    ENDIF

    RETURN
  END SUBROUTINE ricbes

  ELEMENTAL REAL(RP) MODULE FUNCTION symbol_3j(l1, l2, l3, m1, m2, m3)
    INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3

    REAL(RP) :: s, cst
    INTEGER :: t, ti, tf

    symbol_3j = 0._RP
    IF( ABS(M1)>L1 .OR. ABS(M2)>L2 .OR. ABS(M3)>L3 .OR. L3>(L1+L2) .OR. L3<ABS(L1-L2) &
      .OR. M1+M2+M3/=0 ) RETURN

    s = 0._RP
    cst = 0.5_RP*(lnfac(l1+m1) +lnfac(l1-m1) +lnfac(l2+m2) +lnfac(l2-m2) &
      +lnfac(l3+m3) +lnfac(l3-m3) +lnfac(l2+l3-l1) +lnfac(l3+l1-l2) +lnfac(l1+l2-l3) &
      -lnfac(l1+l2+l3+1) )

    ti = MAX(0, l2-l3-m1, l1-l3+m2 )
    tf = MIN(l1+l2-l3, l1-m1, l2+m2 )

    DO t = ti,tf
      s = -s +(-1)**tf *EXP( cst-( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+m1+t) &
        +lnfac(l3-l1-m2+t) +lnfac(l1-m1-t) +lnfac(l2+m2-t) ) )
    END DO

    symbol_3j = s *(-1)**(MOD((l3-m3),2))

  END FUNCTION symbol_3j

  !> Conhyp_opt
  !!
  !> Optimized Confulent Hypergeometric Function when b = 1.
  !! and (a and z) are pure imaginary
  !! @param[in] ai the imaginary part of the original A
  !! @param[in] zi the imaginary part of the original Z
  ELEMENTAL COMPLEX(RP) MODULE FUNCTION conhyp_opt(ai,zi)
    REAL(RP), INTENT(IN) :: ai, zi
    INTEGER :: i
    COMPLEX(RP) :: u

    conhyp_opt = (1._RP,0._RP)
    IF( zi==0._RP ) RETURN
    !
    u = (1._RP, 0._RP)
    DO i = 1,512
      u = (zi/i**2)*CMPLX( -ai, (i-1), RP )*u
      conhyp_opt = conhyp_opt +u
      IF( ABS(u/conhyp_opt)<=EPSILON(1._RP) ) EXIT
    END DO

    RETURN
  END FUNCTION conhyp_opt

END SUBMODULE special_functions