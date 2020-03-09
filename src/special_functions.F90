submodule(special_functions) special_functions
  implicit none

contains

  module subroutine factorial()
    integer :: i

    if(fac_called) return

    fac(0) = 1._wp
    do i = 1,34
      fac(i) = i*fac(i-1)
    end do

    lnfac(0:34) = log( fac )
    do i = 35, 400
      lnfac(i) = lnfac(i-1) +log( real(i,wp) )
    end do

    fac_called = .TRUE.

  end subroutine

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
  elemental module function cgamma(z, mo) result(w)
    complex(wp) :: w
    complex(wp), intent(in)  :: z
    integer, intent(in) :: mo

    ! Local variables
    complex(wp) :: eta, eta2, sum_
    real(wp), parameter :: c0(12) = [ .833333333333333E-01_wp, -.277777777777778E-02_wp &
      , .793650793650794E-03_wp, -.595238095238095E-03_wp,  .841750841750842E-03_wp &
      , -.191752691752692E-02_wp,  .641025641025641E-02_wp, -.295506535947712E-01_wp &
      , .179644372368831_wp    , -1.39243221690590_wp    , 13.4028640441684_wp      &
      , -156.848284626002_wp ], pi = 3.14159265358979_wp, pi2  = 6.28318530717959_wp &
      , alpi = 1.14472988584940_wp, hl2p = .918938533204673_wp, half = 0.5_wp
    real(wp) :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
      s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
    integer :: j, k, l, m, n, nm1
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
    eps = epsilon(1._wp)
    !---------------------------
    x = real(z, wp)
    y = aimag(z)
    if( y==0._wp ) then
      if( mo==0 ) then
        w = cmplx( gamma(x), 0._wp, wp )
      else
        w = cmplx( log_gamma(x), 0._wp, wp )
      end if
      return
    end if
    !IF( ABS(x)>=MIN( REAL(MAX,wp), 1._wp/eps ) ) THEN! GO TO 70
    w = (0._wp,0._wp)

    s2 = 0._wp
    s1 = 0._wp
    s = 0._wp
    h2 = 0._wp
    h1 = 0._wp
    ! RETURN
    !END IF
    if( x<0._wp ) then
      ! CASE WHEN THE REAL PART OF Z IS NEGATIVE
      y = abs(y)
      t = -pi*y
      et = exp(t)
      e2t = et*et

      ! SET A1 = (1 +E2T)/2  AND  A2 = (1 -E2T)/2
      a1 = half*(1._wp +e2t)
      t2 = t +t
      if(t2>=-0.15_wp) then
        a2 = -half*rexp(t2)
      else
        a2 = half*(half +(half -e2t))
      end if

      !  COMPUTE SIN(PI*X) AND COS(PI*X)
      !IF(ABS(x)>=MIN(REAL(MAX,wp), 1._wp/eps)) GO TO 70
      k = int( abs(x) )
      u = x +k
      k = mod(k,2)
      if(u <= -half) then
        u = half +(half +u)
        k = k +1
      end if
      u = pi*u
      sn = sin(u)
      cn = cos(u)
      if(k==1) then
        sn = -sn
        cn = -cn
      end if
      !  SET  H1 +H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))
      a1 = sn*a1
      a2 = cn*a2
      a = a1*a1 +a2*a2
      if( a==0._wp ) return!GO TO 70
      if( mo==0 ) then
        h1 = a1 / a
        h2 = -a2 / a
        c = pi*et
        h1 = c*h1
        h2 = c*h2
      else
        h1 = (alpi+t) -half*log(a)
        h2 = -atan2(a2,a1)
      end if
      if( aimag(z)>=0._wp ) then
        x = 1._wp -x
        y = -y
      else
        h2 = -h2
        x = 1._wp -x
      end if
    end if
    !-----------------------------------------------------------------------
    !           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
    !-----------------------------------------------------------------------
    w1 = 0._wp
    w2 = 0._wp
    n = 0
    t = x
    y2 = y*y
    a = t*t +y2
    cut = 36._wp
    if( eps>1.E-8_wp ) cut = 16._wp
    if( a<cut ) then
      if( a==0._wp) return
      do while( a<cut )
        n = n +1
        t = t +1._wp
        a = t*t +y2
      end do
      !     LET S1 +S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)
      u1 = (x*t+y2) / a
      u2 = y / a
      s1 = u1
      s2 = n*u2
      if( n>=2 ) then
        u = t / a
        nm1 = n -1
        do j = 1, nm1
          v1 = u1 +j*u
          v2 = (n-j)*u2
          c = s1*v1 -s2*v2
          d = s1*v2 +s2*v1
          s1 = c
          s2 = d
        end do
      end if

      !     SET  W1 +W2*I = LOG(S1 +S2*I)  WHEN MO IS NONZERO

      s = s1*s1 +s2*s2
      if( mo/=0 ) then
        w1 = half*log(s)
        w2 = atan2(s2,s1)
      end if
    end if

    !     SET  V1 +V2*I = (Z -0.5)*LOG(Z +N) -Z

    t1 = half*log(a) -1._wp
    t2 = atan2(y,t)
    u = x -half
    v1 = (u*t1-half) -y*t2
    v2 = u*t2 +y*t1

    !     LET A1 +A2*I BE THE ASYMPTOTIC SUM

    eta = cmplx(t/a, -y/a, wp)
    eta2 = eta*eta
    m = 12
    if( a>=289._wp ) m = 6
    if( eps>1.e-8_wp ) m = m / 2
    sum_ = cmplx(c0(m), 0._wp, wp)
    l = m
    do j = 2, m
      l = l -1
      sum_ = cmplx(c0(l), 0._wp, wp) +sum_*eta2
    end do
    sum_ = sum_*eta
    a1 = real(sum_, wp)
    a2 = aimag(sum_)
    !-----------------------------------------------------------------------
    !                 GATHERING TOGETHER THE RESULTS
    !-----------------------------------------------------------------------
    w1 = (((a1 +hl2p) -w1) +v1) -n
    w2 = (a2 -w2) +v2
    if( x>=0._wp .AND. mo==0 ) then ! CASE WHEN Re(Z)>0 AND MO==0
      a = exp(w1)
      w1 = a*cos(w2)
      w2 = a*sin(w2)
      if( n/=0 ) then
        c = (s1*w1 +s2*w2) / s
        d = (s1*w2 -s2*w1) / s
        w1 = c
        w2 = d
      end if
    elseif(M/=0) then ! CASE WHEN MO/=0
      if( x<0._wp ) then ! CASE WHEN Re(Z)<0 AND M/=0
        w1 = h1 -w1
        w2 = h2 -w2
      end if
      ! CASE WHEN Re(Z)>0 AND MO/=0
      ! THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI<W2 <= PI.
      if( w2<=pi ) then
        k = int( half -w2 / pi2 )
        w2 = w2 +pi2*k
      else
        k = int ( w2 / pi2 -half )
        w2 = w2 -pi2*real(k+1,wp)
        if( w2<=-pi ) w2 = pi
      end if
    else ! CASE WHEN Re(Z)<0 AND MO==0
      a = exp(-w1)
      t1 = a*cos(-w2)
      t2 = a*sin(-w2)
      w1 = h1*t1 -h2*t2
      w2 = h1*t2 +h2*t1
      if( n/=0 ) then
        c = w1*s1 -w2*s2
        d = w1*s2 +w2*s1
        w1 = c
        w2 = d
      end if
    end if

    w = cmplx(w1, w2, wp)

    return

  contains
    !-----------------------------------------------------------------------
    !            EVALUATION OF THE FUNCTION EXP(X) -1
    !-----------------------------------------------------------------------
    elemental real(wp) function rexp(x)
      real(wp), intent(in) :: x
      ! Local variables
      real(wp), parameter  :: p1 = .914041914819518E-9_wp, p2 = .238082361044469E-1_wp &
        , q1 = -.499999999085958_wp, q2 = .107141568980644_wp &
        , q3 = -.119041179760821E-1_wp, q4 = 0.595130811860248E-3_wp
      real(wp) :: e
      !-----------------------
      if( abs(x)<=0.15_wp ) then
        rexp = x*(((p2*x +p1)*x +1._wp) / ( ( ( (q4*x +q3)*x +q2)*x +q1 )*x +1._wp) )
      elseif( x>=0._wp ) then
        e = exp(x)
        rexp = e*(half +(half -1._wp/e))
      elseif( x>=-37._wp ) then
        rexp = (exp(x) -half) -half
      else
        rexp = -1._wp
      end if

      return
    end function rexp

  end function cgamma

  !> assoc_legendre
  !! Computes the associated Legendre polynomial P_l^m (x). Here m and l are integers
  !! satisfying 0 <= m <= l,  while x lies in the range −1 <= x <= 1.
  elemental real(wp) module function assoc_legendre(l,m,x)
    use ieee_arithmetic, only : ieee_is_finite
    integer, intent(in) :: l, m
    real(wp), intent(in) :: x

    integer :: i, l1
    real(wp) :: fact, pmm, pmm1, pmm2

    !IF( m<0 .OR. m>l ) ERROR STOP 'M is out of range [0,L]'
    !IF( ABS(x)>1._wp ) ERROR STOP 'X is out of range [-1., 1.]'

    !Compute P_m^m
    pmm = 1._wp
    if( m>0 ) then
      fact = 1._wp
      do i = 1, m
        pmm = pmm *fact
        fact = fact +2._wp
      end do
      pmm = pmm *( -sqrt(1._wp-x**2) )**m
    endif

    if( m==l ) then ! Compute P_l^l
      assoc_legendre = pmm
    elseif( m==l-1 ) then  ! Compute P_l^{l-1}
      assoc_legendre = x*(2*m+1) *pmm
    else ! Compute P_l^m, m<l-1
      pmm1 = 0._wp
      do l1 = m+1, l
        pmm2 = pmm1
        pmm1 = pmm
        pmm = ( x*(2*l1-1)*pmm1 -(l1+m-1)*pmm2 ) /(l1-m)
      end do
      assoc_legendre = pmm
    endif
    if( .not. ieee_is_finite(pmm) ) assoc_legendre = huge(1._wp)

    return
  end function assoc_legendre

  elemental complex(wp) module function spherical_harmonic( l, m, theta, phi )
    use constants, only : pi
    integer, intent(in) :: l, m
    real(wp), intent(in) :: theta, phi
    real(wp), parameter :: Tinye = log(tiny(1._wp))
    integer :: ma

    !if(.not. fac_called ) error stop 'you should call factorial before using &
    !  &y_l^m(\theta,\phi)'

    ma = abs(m)
    if( (lnfac(l-ma)-lnfac(l+ma))<2*Tinye ) then
      spherical_harmonic = (0._wp, 0._wp)
      return
    elseif( l==0 .AND. m==0 ) then
      spherical_harmonic = 1._wp/sqrt(4._wp*pi)
      return
    elseif( m==0 ) then
      spherical_harmonic = sqrt( (2.*l+1._wp)/(4.*pi) ) *assoc_legendre(l, 0, cos(theta) )
      return
    end if

    spherical_harmonic = exp( 0.5*( lnfac(l-ma) -lnfac(l+ma) ) ) &
      *sqrt( (2.*l+1._wp)/(4.*pi) ) *cmplx( cos(m*phi), sin(m*phi), wp ) &
      *assoc_legendre(l, ma, cos(theta) )
    if( mod(m,2)>0 ) spherical_harmonic = -spherical_harmonic
    if( modulo(theta,2*pi)>pi .AND. mod(m,2)/=0 ) spherical_harmonic = -spherical_harmonic

  end function spherical_harmonic

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
  pure module subroutine coul90(x, eta, lmin, lrange, fc, gc, fcp, gcp, kfn, ifail )
    integer, intent(in)  :: lmin, lrange, kfn
    integer, intent(out), optional :: ifail
    real(wp), intent(in) :: x, eta
    real(wp), intent(out), dimension(lmin:lmin+lrange) :: fc,  gc,  fcp, gcp

    integer, parameter :: limit = 20000
    real(wp), parameter :: small = sqrt(tiny(1._wp)), ten2 = 100._wp, half = 0.5_wp &
      , rt2dpi = 0.797884560802865_wp
    !!---- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM IDUM3
    real(wp) :: accur, acch, xinv, pk, cf1, c, d, pk1, etak, rk2, tk, dcf1, den, xlm &
      , xll, el, xl, rl, sl, f, fcmaxl, fcminl, gcminl, omega, wi, a, b, ar, ai, br &
      , bi, dr, di, dp, dq, alpha, beta, e2mm1, fjwkb,  gjwkb, p, q, gamma_, gammai, err
    integer :: l, maxl, idum2, idum3
    logical :: etane0, xlturn
    !----------------------------------------------------------------------
    !     COUL90 HAS CALLS TO: SQRT,ABS,MAX,INT,SIGN,DBLE,MIN
    !----------------------------------------------------------------------
    !Q    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
    !---- THIS CONSTANT IS  SQRT(2._wp / PI):
    !---- USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND REAL
    !---- CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED

    accur = epsilon(1._wp)
    !    accur = 1.E-7_wp
    if(present(ifail)) ifail = 0
    idum2 = 1
    !idum1 = 0
    gjwkb = 0._wp
    err = 1._wp
    ! IF(KFN/=0) ETA = 0._wp
    etane0 = ( eta/=0._wp .and. kfn==0 )
    acch = sqrt(accur)
    gammai = 0._wp
    !!---- TEST RANGE OF X, EXIT IF<=SQRT(ACCUR) OR IF NEGATIVE
    if( x <= acch .and. present(ifail) ) then
      ifail = -1
      !WRITE(6,1000) X,ACCH
      !1000 FORMAT(' FOR X = ',1P,D12.3,'     TRY SMALL-X  SOLUTIONS' &
      ! ,' OR X IS NEGATIVE'/, ' SQUARE ROOT (ACCURACY) =  ',D12.3/)
      return
    endif

    if( kfn==2 ) then
      xlm = lmin -half
    else
      xlm = lmin
    endif
    if( xlm <= -1._wp .OR. lrange<0 .and. present(ifail) ) then
      ifail = -2
      !  WRITE (6,1001) LRANGE,LMIN,XLM
      !1001 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',2I10,1P,D15.6/)
      return
    endif
    e2mm1 = xlm*xlm +xlm
    xlturn = x*(x -2._wp*eta)<e2mm1
    e2mm1 = e2mm1 +eta*eta
    xll = xlm +real(lrange,wp)
    !---- LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
    !---- XLL    IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
    !---- DETERMINE STARTING ARRAY ELEMENT (IDUM3) FROM XLMIN
    idum3 = lmin!MAX( INT(XLMIN +ACCUR),0 )
    maxl = idum3 +lrange
    !---- EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
    xinv = 1._wp / x
    !---- UNNORMALISED F(MAXL,ETA,X)
    den = 1._wp
    pk = xll +1._wp
    cf1 = eta / pk +pk*xinv
    if( abs(cf1)<small ) cf1 = small
    rk2 = 1._wp
    d = 0._wp
    c = cf1
    !---- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA +1: LENTZ-THOMPSON
    do l = 1,  limit
      pk1 = pk +1._wp
      if( etane0 ) then
        etak = eta / pk
        rk2 = 1._wp +etak*etak
        tk = (pk +pk1)*(xinv +etak / pk1)
      else
        tk = (pk +pk1)*xinv
      endif
      !---- DIRECT RATIO OF B CONVERGENTS
      d = tk -rk2*d
      !---- INVERSE RATIO OF A CONVERGENTS
      c = tk -rk2 / c
      if( abs(c)<small ) c = small
      if( abs(d)<small ) d = small
      d = 1._wp / d
      dcf1 =  d*c
      cf1 = cf1*dcf1
      if( d<0._wp ) den = -den
      pk = pk1
      !----  PROPER EXIT
      if( abs(dcf1-1._wp)<accur ) exit
    end do
    !---- ERROR EXIT
    if( abs(dcf1-1._wp)>=accur .and. present(ifail) ) then
      ifail = 1
      !WRITE (6,1002) LIMIT, CF1,DCF1, PK,ACCUR
      !1002  FORMAT('CF1 HAS FAILED TO CONVERGE AFTER',I10,'ITERATIONS',/'CF1,DCF1,PK
      !,ACCUR = ',1P,4D12.3/)
      return
    endif

    !nfp = INT(pk -xll -1)
    f = cf1
    !---- DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
    if( lrange>0 ) then
      fcmaxl = small*den
      fcp(maxl) = fcmaxl*cf1
      fc (maxl) = fcmaxl
      xl = xll
      rl = 1._wp
      do l = maxl, idum3+1, -1
        if( etane0 ) then
          el = eta / xl
          rl = sqrt( 1._wp +el*el )
          sl = xl*xinv +el
          gc (l) = rl
          gcp(l) = sl
        else
          sl = xl*xinv
        endif
        fc (l-1) = ( fc(l)*sl +fcp(l) ) / rl
        fcp(l-1) = fc(l-1)*sl -fc (l)*rl
        !---- END VALUE IS XLM
        xl = xl -1._wp
      end do
      if( abs(fc(idum3))<accur*small ) fc(idum3) = accur*small
      f = fcp(idum3) / fc(idum3)
      den = fc (idum3)
    endif
    !---------------------------------------------------------------------
    !---- NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
    !---- EVALUATE CF2 = P +I.Q  USING STEED'S ALGORITHM (NO ZEROS)
    !---------------------------------------------------------------------
    if( xlturn ) call jwkb(x,eta,max(xlm,0._wp),fjwkb,gjwkb,idum2)
    if( idum2>1 .OR. gjwkb>1._wp/(acch*ten2) ) then
      omega = fjwkb
      gamma_ = gjwkb*omega
      p = f
      q = 1._wp
    else
      xlturn = .false.
      pk = 0._wp
      wi = eta +eta
      p = 0._wp
      q = 1._wp -eta*xinv
      ar = -e2mm1
      ai = eta
      br = 2._wp*(x -eta)
      bi = 2._wp
      dr = br / (br*br +bi*bi)
      di = -bi / (br*br +bi*bi)
      dp = -xinv*(ar*di +ai*dr)
      dq = xinv*(ar*dr -ai*di)
      do l = 1, limit
        p = p +dp
        q = q +dq
        pk = pk +2._wp
        ar = ar +pk
        ai = ai +wi
        bi = bi +2._wp
        d = ar*dr -ai*di +br
        di = ai*dr +ar*di +bi
        c = 1._wp / (d*d +di*di)
        dr = c*d
        di = -c*di
        a = br*dr -bi*di -1._wp
        b = bi*dr +br*di
        c = dp*a -dq*b
        dq = dp*b +dq*a
        dp = c
        if( abs(dp)+abs(dq)<(abs(p)+abs(q))*accur ) exit
      end do
      !---- ERROR EXIT
      if( abs(dp)+abs(dq)>=(abs(p)+abs(q))*accur .and. present(ifail) ) then
        ifail = 2
        !WRITE (6,1003) LIMIT,P,Q,DP,DQ,ACCUR
        !1003    FORMAT('CF2 HAS FAILED TO CONVERGE AFTER',I7,'ITERATIONS',&
          !               &/'P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3/)
        return
      endif

      !idum1 = INT(pk / 2._wp)
      err = half*accur / min( abs(q),1._wp )
      if( abs(p)>abs(q) ) err = err*abs(p)
      !---------------------------------------------------------------------
      !    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
      !---------------------------------------------------------------------
      gamma_ = (f -p)/q
      gammai = 1._wp/gamma_
      if( abs(gamma_) <= 1._wp ) then
        omega = sqrt( 1._wp +gamma_*gamma_ )
      else
        omega = sqrt( 1._wp +gammai* gammai)*abs(gamma_)
      endif
      omega = 1._wp / ( omega*sqrt(q) )
      ! wronsk = omega
    endif
    !---------------------------------------------------------------------
    !---- RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
    !---------------------------------------------------------------------
    if(kfn==1) then !---- SPHERICAL
      alpha = xinv
      beta = xinv
    else if(kfn==2) then !----   CYLINDRICAL
      alpha = half*xinv
      beta = sqrt( xinv )*rt2dpi
    else !----   KFN = 0, COULOMB FUNCTIONS
      alpha = 0._wp
      beta = 1._wp
    endif
    fcminl = sign( omega,den )*beta
    if( xlturn ) then
      gcminl = gjwkb*beta
    else
      gcminl = fcminl*gamma_
    endif
    !---- BESSEL SIGN
    if( kfn/=0 ) gcminl = -gcminl
    fc (idum3) = fcminl
    gc (idum3) = gcminl
    gcp(idum3) = gcminl*(p -q*gammai -alpha)
    fcp(idum3) = fcminl*(f -alpha)
    if( lrange==0 ) return
    !---------------------------------------------------------------------
    !---- UPWARD RECURRENCE FROM GC(IDUM3),GCP(IDUM3) STORED VALUES ARE RL,SL
    !---- RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
    !---- XL   = XLM HERE  AND RL = 1._wp,  EL = 0._wp FOR BESSELS
    !---------------------------------------------------------------------
    omega = beta*omega / abs(den)
    xl = xlm
    rl = 1._wp
    do l = idum3+1,  maxl
      xl = xl +1._wp
      if( etane0 ) then
        rl = gc (l)
        sl = gcp(l)
      else
        sl = xl*xinv
      endif
      gc (l) = ( (sl -alpha)*gc(l-1) -gcp(l-1) ) / rl
      gcp(l) = rl*gc(l-1) -(sl +alpha)*gc(l)
      fcp(l) = omega*( fcp(l) -alpha*fc(l) )
      fc (l) = omega*fc (l)
    end do

    if(present(ifail) ) ifail = 0

    return
  end subroutine coul90
  !----------------------------------------------------------------------
  !> COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
  !! AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
  !! CALCULATED IN SINGLE, RETURNED IN REAL VARIABLES
  !! CALLS MAX, SQRT, LOG, EXP, ATAN2, REAL, INT
  !!     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
  !----------------------------------------------------------------------
  elemental subroutine  jwkb( x, eta, xl, fjwkb, gjwkb, iexp )
    real(wp), intent(in)  :: x, eta, xl
    real(wp), intent(out) :: fjwkb, gjwkb
    integer, intent(out) :: iexp

    integer, parameter :: maxexp = 300
    real(wp), parameter :: half = 0._wp, six = 6._wp, ten = 10._wp &
      , dzero = 0._wp, rl35 = 35._wp, aloge = 0.4342945_wp

    real(wp) :: phi, phi10, gh2, xll1, hll, hl, sl, rl2, gh
    !----------------------------------------------------------------------
    !---- CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR REAL
    !----------------------------------------------------------------------
    gh2   = x*(eta +eta -x)
    xll1  = max( xl*xl +xl, dzero )
    if( gh2 +xll1 <= 0._wp ) return
    hll = xll1 +six / rl35
    hl = sqrt(hll)
    sl = eta / hl +hl / x
    rl2 = 1._wp +eta*eta / hll
    gh = sqrt(gh2 +hll) / x
    phi = x*gh -half*( hl*log((gh +sl)**2 / rl2) -log(gh) )
    if( eta/=0._wp ) phi = phi -eta*atan2(x*gh,x -eta)
    phi10 = -phi*aloge
    iexp = int(phi10)
    if( iexp>maxexp ) then
      gjwkb = ten**(phi10 -real(iexp,wp))
    else
      gjwkb = exp(-phi)
      iexp = 0
    endif
    fjwkb = half / (gh*gjwkb)

    return
  end subroutine  jwkb

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
  pure module subroutine ricbes( x, lmax, psi, chi, psid, chid, ifail )
    real(wp), intent(in) :: x
    integer, intent(in) :: lmax
    integer, intent(inout), optional :: ifail
    real(wp), intent(out) :: psi(0:lmax), chi(0:lmax), psid(0:lmax), chid(0:lmax)

    integer, parameter :: limit = 20000 !,maxl = 1001
    real(wp), parameter :: small = sqrt(tiny(1._wp)) , three = 3._wp

    integer :: l
    real(wp) :: accur, tk, sl!, ERR
    real(wp) :: xinv, cf1, dcf1, den, c, d, omega, twoxi
    !COMMON / STEED  / ERR,NFP,IDUM1,IDUM2,IDUM3
    !!----
    !!---- CALCULATE THE L = 0   RICCATI-BESSEL FUNCTIONS DIRECTLY
    psi (0) = sin(x)
    chi (0) = -cos(x)
    psid(0) = -chi(0)
    chid(0) = psi(0)
    accur = epsilon(1._wp)!1.E-17
    if( present(ifail) ) ifail = -1
    if( x<sqrt(accur) ) goto 50
    !!---- TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !!---- BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
    xinv = 1._wp / x
    if( lmax>0 ) then
      twoxi = xinv +xinv
      sl = real(lmax,wp)* xinv
      tk = 2._wp*sl +xinv*three
      cf1 = sl +xinv
      den = 1._wp
      if( abs(cf1)<small ) cf1 = small
      !!---- INVERSE RATIO OF A CONVERGENTS
      c = cf1
      !!---- DIRECT  RATIO OF B CONVERGENTS
      d = 0._wp
      do l = 1,limit
        c = tk -1._wp / c
        d = tk -d
        if( abs(c)<small ) c = small
        if( abs(d)<small ) d = small
        d = 1._wp / d
        dcf1 =  d*c
        cf1 = cf1*dcf1
        if( d<0._wp ) den = -den
        if( abs(dcf1 -1._wp)<=accur ) exit
        tk = tk +twoxi
      end do
      !!---- ERROR EXIT, NO CONVERGENCE
      if( l>=limit ) goto 50
      ! 20 CONTINUE !   nfp = l
      !!---- ERROR ESTIMATE
      !ERR = accur*SQRT(REAL(nfp,KIND=RP))
      psi (lmax) = den
      psid(lmax) = cf1*den
      !!---- DOWNWARD RECURSION TO L = 0  AS RICCATI-BESSEL FUNCTIONS
      do l = lmax,  2, -1
        psi (l-1) = sl*psi(l) +psid(l)
        psid(l-1) = sl*psi(l-1) -psi (l)
        sl = sl -xinv
      end do
      den  = sl*psi(1) +psid(1)!PSI(0)
      !!---- END LOOP FOR LMAX>0
      !ENDIF
      !IF(LMAX>0) THEN
      omega  = psi(0) / den
      do l = 1,  lmax
        psi (l) = omega*psi (l)
        psid(l) = omega*psid(l)
        sl = xinv*real(l,wp)
        chi (l) = sl*chi(l-1) -chid(l-1)
        chid(l) = chi(l-1) -sl*chi (l)
      end do
    endif
    !!---- CALCULATIONS SUCCESSFUL
    if(present(ifail)) ifail = 0
    return
    !!---------------------------------------------------------------------
    !!---- ERROR TRAPS
    !!---------------------------------------------------------------------
    50 continue
    if( x<0._wp ) then
      !    WRITE(6,1000) x
      !    1000   FORMAT(' X NEGATIVE !',1p,e15.5,' ... USE REFLECTION FORMULA'/)
    else if( x==0._wp ) then
      ifail = 0
      chi(0) = 1._wp
      do l = 1, lmax
        chi(l) = 0._wp
      end do
    else
      ! WRITE(6,1001) x
      ! 1001   FORMAT(' WITH X = ',1p,e15.5,'    TRY SMALL-X SOLUTIONS',/ &
        !   , '  PSI_L(X)=>X**(L+1)/(2L+1) AND',/,'  CHI/L/(X)-> -(2L-1)/X**L'/)
    endif

    return
  end subroutine ricbes

  elemental real(wp) module function symbol_3j(l1, l2, l3, m1, m2, m3)
    integer, intent(in) :: l1, l2, l3, m1, m2, m3

    real(wp) :: s, cst
    integer :: t, ti, tf

    symbol_3j = 0._wp
    if( abs(M1)>L1 .OR. abs(M2)>L2 .OR. abs(M3)>L3 .OR. L3>(L1+L2) .OR. L3<abs(L1-L2) &
      .OR. M1+M2+M3/=0 ) return

    s = 0._wp
    cst = 0.5_wp*(lnfac(l1+m1) +lnfac(l1-m1) +lnfac(l2+m2) +lnfac(l2-m2) &
      +lnfac(l3+m3) +lnfac(l3-m3) +lnfac(l2+l3-l1) +lnfac(l3+l1-l2) +lnfac(l1+l2-l3) &
      -lnfac(l1+l2+l3+1) )

    ti = max(0, l2-l3-m1, l1-l3+m2 )
    tf = min(l1+l2-l3, l1-m1, l2+m2 )

    do t = ti,tf
      s = -s +(-1)**tf *exp( cst-( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+m1+t) &
        +lnfac(l3-l1-m2+t) +lnfac(l1-m1-t) +lnfac(l2+m2-t) ) )
    end do

    symbol_3j = s *(-1)**(mod((l3-m3),2))

  end function symbol_3j

  !> Conhyp_opt
  !!
  !> Optimized Confulent Hypergeometric Function when b = 1.
  !! and (a and z) are pure imaginary
  !! @param[in] ai the imaginary part of the original A
  !! @param[in] zi the imaginary part of the original Z
  elemental complex(wp) module function conhyp_opt(ai,zi)
    real(wp), intent(in) :: ai, zi
    integer :: i
    complex(wp) :: u
    real(wp), parameter :: eps = sqrt( epsilon(1._wp) )

    conhyp_opt = (1._wp,0._wp)
    if( zi==0._wp ) return
    !
    u = (1._wp, 0._wp)
    do i = 1, 10240
      u = (zi/i**2)*cmplx( -ai, (i-1), wp )*u
      conhyp_opt = conhyp_opt +u
      if( abs(u/conhyp_opt)<=eps ) exit
    end do

    return
  end function conhyp_opt

end submodule special_functions
