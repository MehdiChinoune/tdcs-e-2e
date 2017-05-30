SUBMODULE (fdcs_e2e) fdcs_e2e
  USE constants         ,ONLY : RP ,pi ,deg ,eV
  USE trigo             ,ONLY : cartez2spher ,spher2cartez
  USE utils             ,ONLY : norm_fac ,read_input ,read_orbit ,fac ,factorial ,y1y2y3 &
    ,ode_second_dw ,calculate_U
  USE special_functions ,ONLY : spherical_harmonic ,cgamma ,coul90
  USE iso_fortran_env   ,ONLY : INPUT_UNIT ,OUTPUT_UNIT ,ERROR_UNIT
  IMPLICIT NONE
CONTAINS

  MODULE SUBROUTINE fdcs_fba_pw(in_unit,out_unit)
    INTEGER ,INTENT(IN) :: in_unit
    INTEGER ,INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom ,Orbit
    REAL(KIND=RP) :: Ei ,Es ,Ee
    REAL(KIND=RP) :: thetas
    INTEGER       :: lo ,no
    REAL(KIND=RP) ,ALLOCATABLE :: a(:) ,e(:)
    INTEGER ,ALLOCATABLE :: n(:)
    INTEGER       :: step(3)

    REAL(KIND=RP) ,PARAMETER   :: phie = 0._RP
    REAL(KIND=RP) :: kim ,ksm ,kem ,km
    REAL(KIND=RP) :: ki(3) ,ks(3) ,ke(3) ,k(3)

    REAL(KIND=RP) :: factor ,sigma
    COMPLEX(KIND=RP) :: term
    INTEGER :: i,io,mo

    CALL factorial()

    CALL read_input(in_unit,Ei ,Es ,Ee ,thetas ,step ,Atom ,orbit)

    CALL read_orbit(Atom//'_'//Orbit ,lo ,no ,n ,a ,e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    CALL spher2cartez( kim ,0._RP ,0._RP ,ki )
    CALL spher2cartez( ksm ,thetas*deg ,pi ,ks )
    k = ki -ks
    km = NORM2(k)

    factor=8.*ksm*kem/(kim*km**4)

    DO i = step(1) ,step(2) ,step(3)

      CALL spher2cartez( kem ,i*deg ,phie ,ke )

      sigma = 0._RP
      DO mo = 0 ,lo
        term=(0._RP,0._RP)
        DO io = 1 ,no
          term = term +a(io) *( tpw( n(io) ,lo ,mo ,e(io) ,ke ,k ) -tpw( n(io) ,lo ,mo ,e(io) ,ke ) )
        END DO
        sigma = sigma + (1+mo)*ABS(term)**2
      END DO

      sigma=factor*sigma

      PRINT'(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit ,'(1x,i4,1x,es15.8)' ) i ,sigma
    END DO

  END SUBROUTINE fdcs_fba_pw

  MODULE SUBROUTINE fdcs_fba_cw(in_unit,out_unit)
    INTEGER ,INTENT(IN) :: in_unit
    INTEGER ,INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom ,Orbit
    REAL(KIND=RP) :: Ei ,Es ,Ee
    REAL(KIND=RP) :: thetas
    INTEGER       :: lo ,no
    REAL(KIND=RP) ,ALLOCATABLE :: a(:) ,e(:)
    INTEGER ,ALLOCATABLE :: n(:)
    INTEGER       :: step(3)

    REAL(KIND=RP) ,PARAMETER   :: phie = 0._RP
    REAL(KIND=RP) :: kim ,ksm ,kem ,km
    REAL(KIND=RP) :: ki(3) ,ks(3) ,ke(3) ,k(3)

    REAL(KIND=RP) :: factor ,sigma
    COMPLEX(KIND=RP) :: term
    INTEGER :: i ,io ,mo

    CALL factorial()

    CALL read_input(in_unit,Ei ,Es ,Ee ,thetas ,step ,Atom ,orbit)

    CALL read_orbit(Atom//'_'//Orbit ,lo ,no ,n ,a ,e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    CALL spher2cartez( kim ,0._RP ,0._RP ,ki )
    CALL spher2cartez( ksm ,thetas*deg ,pi ,ks )
    k = ki -ks
    km = NORM2(k)

    factor = 8. *ksm *kem / (kim*km**4)
    factor = factor *2.*pi/(kem*(1.-EXP(-2.*pi/kem))) ! \abs{ \exp{\pi\alpha/2} *\Gamma(1-i\alpha) }^2

    DO CONCURRENT( i = step(1) :step(2) :step(3) )

      CALL spher2cartez( kem ,i*deg ,phie ,ke )

      sigma = 0._RP
      DO mo = 0 ,lo
        term=(0._RP,0._RP)
        DO io = 1 ,no
          term = term +a(io) *( tcw( n(io) ,lo ,mo ,e(io) ,ke ,k ) -tcw0( n(io) ,lo ,mo ,e(io) ,ke ) )
        END DO
        sigma = sigma +(mo+1)*ABS(term)**2
      END DO

      sigma=factor*sigma

      PRINT'(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit ,'(1x,i4,1x,es15.8)' ) i ,sigma
    END DO

  END SUBROUTINE fdcs_fba_cw

  MODULE SUBROUTINE fdcs_fba_dw(in_unit,out_unit)
    INTEGER ,INTENT(IN) :: in_unit
    INTEGER ,INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom ,Orbit
    REAL(KIND=RP) :: Ei ,Es ,Ee
    REAL(KIND=RP) :: thetas
    INTEGER       :: lo ,no
    REAL(KIND=RP) ,ALLOCATABLE :: a(:) ,e(:)
    INTEGER ,ALLOCATABLE :: n(:)
    INTEGER       :: step(3)

    REAL(KIND=RP) ,PARAMETER   :: phie = 0._RP
    REAL(KIND=RP) :: kim ,ksm ,kem ,km
    REAL(KIND=RP) :: theta ,phi
    REAL(KIND=RP) :: ki(3) ,ks(3) ,ke(3) ,k(3)

    INTEGER ,PARAMETER :: lmax = 15 ,lemax = 7 ,np = 800
    REAL(KIND=RP) :: x(0:np) ,wf(0:np) ,chi_b(0:np,0:lemax) ,chi_0a(0:np,0:lmax) ,sigma_l(0:lemax)
    REAL(KIND=RP) :: etae ,rc ,h ,integral(0:lemax,0:lmax)
    COMPLEX(KIND=RP) :: tmp_z
    COMPLEX(KIND=RP) ,PARAMETER :: zi = (0._RP ,1._RP)
    INTEGER :: le ,l ,me ,z=1

    REAL(KIND=RP) :: factor ,sigma
    COMPLEX(KIND=RP) :: term ,term0
    INTEGER :: i ,io ,mo

    CALL factorial()

    CALL read_input(in_unit,Ei ,Es ,Ee ,thetas ,step ,Atom ,orbit)

    CALL read_orbit(Atom//'_'//Orbit ,lo ,no ,n ,a ,e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    CALL spher2cartez( kim ,0._RP ,0._RP ,ki )
    CALL spher2cartez( ksm ,thetas*deg ,pi ,ks )
    k = ki -ks
    CALL cartez2spher( k ,km ,theta ,phi)

    factor = 8. *ksm *kem / (kim*km**4)
    factor = factor *2./(pi*kem**2)

    rc = 10._RP
    h = rc/np

    etae = -z/kem

    x = [(i*h ,i = 0 ,np) ]
    chi_b(0,:) = 0._RP
    chi_0a(0,:) = 0._RP
    DO CONCURRENT(i = 1:np)
      BLOCK
        !REAL(KIND=RP) :: gc_a(0:lemax) ,fdc_a(0:lemax) ,gdc_a(0:lemax) ,fc_a(0:lemax)
        REAL(KIND=RP) :: gc_0a(0:lmax) ,fdc_0a(0:lmax) ,gdc_0a(0:lmax) ,fc_0a(0:lmax)
        !CALL coul90(kem*x(i) ,etae ,0 ,lemax ,fc_a ,gc_a ,fdc_a ,gdc_a ,1-z ,le)
        !chi_b(i,:) = fc_a
        CALL coul90(km*x(i) ,0._RP ,0 ,lmax ,fc_0a ,gc_0a ,fdc_0a ,gdc_0a ,1 ,le)
        chi_0a(i,:) = fc_0a
      END BLOCK
    END DO

    BLOCK
      REAL(KIND=RP) :: U(0:np,0:lemax) ,delta(0:lemax) ,U_tmp(0:np)
      INTEGER :: l

      CALL calculate_U(Atom ,Orbit ,x ,U_tmp)

      U(0,:)=HUGE(1._RP)
      U(1:np,lemax) = -kem**2 -2.*( 1._RP +U_tmp(1:np) ) /x(1:np)
      DO CONCURRENT( l = 0:lemax )
        U(1:np,l) = l*(l+1)/x(1:np)**2 +U(1:np,lemax)
      END DO

      CALL ode_second_dw(kem ,lemax ,rc ,U ,chi_b ,delta)

    END BLOCK

    wf = 0._RP
    DO io = 1 ,no
      wf = wf +a(io)*norm_fac(e(io) ,n(io) ) *x**n(io) *EXP(-e(io)*x)
    END DO

    DO CONCURRENT(le = 0:lemax)
      tmp_z = cgamma( CMPLX(le+1._RP ,etae ,KIND=RP) ,0 )
      sigma_l(le) = ATAN2( AIMAG(tmp_z) ,REAL(tmp_z,KIND=RP) )
    END DO

    DO CONCURRENT(le=0:lemax ,l=0:lmax ,MOD(le+l+lo,2)==0)
      integral(le ,l) = 0.5*h*SUM( wf(1:np)*chi_b(1:np,le)*chi_0a(1:np,l) &
        + wf(0:np-1)*chi_b(0:np-1,le)*chi_0a(0:np-1,l) )
    END DO

    term0 = zi**(-lo) *EXP( CMPLX(0._RP ,-sigma_l(lo) ,KIND=RP ) ) &
      *0.5*h*SUM( wf(1:np)*chi_b(1:np,lo) + wf(0:np-1)*chi_b(0:np-1,lo) )

    DO i = step(1) ,step(2) ,step(3)

      CALL spher2cartez( kem ,i*deg ,phie ,ke )

      sigma = 0._RP
      DO mo = 0 ,lo

        term = (0._RP ,0._RP)
        DO CONCURRENT(le=0:lemax ,l=0:lmax ,MOD(le+l+lo,2)==0 )
          tmp_z = (0._RP ,0._RP)
          DO CONCURRENT(me=-le:le ,ABS(me-mo)<=l )
            tmp_z = tmp_z + y1y2y3(le ,l ,lo ,-me ,me-mo ,mo) *spherical_harmonic(le ,me ,i*deg ,phie) &
              *spherical_harmonic(l ,-me+mo ,theta ,phi) *(-1)**(mo)
          END DO
          term = term +zi**(l-le) *EXP( CMPLX(0._RP ,-sigma_l(le) ,KIND=RP ) ) *integral(le ,l) *tmp_z
        END DO
        term = 4.*pi *term

        sigma = sigma +(mo+1)*ABS( term -spherical_harmonic(lo ,mo ,i*deg ,phie) *term0 )**2

      END DO

      sigma = factor*sigma

      PRINT'(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit ,'(1x,i4,1x,es15.8)' ) i ,sigma
    END DO

  END SUBROUTINE fdcs_fba_dw

  PURE COMPLEX(KIND=RP) FUNCTION tpw( n ,l ,m ,e ,ke ,k)
    INTEGER       ,INTENT(IN) :: n,l,m
    REAL(KIND=RP) ,INTENT(IN) :: e,ke(3)
    REAL(KIND=RP) ,INTENT(IN),OPTIONAL :: k(3)

    REAL(KIND=RP)    :: q(3)
    INTEGER          :: j
    REAL(KIND=RP)    :: kem,thetae,phie,A
    COMPLEX(KIND=RP) :: kec

    IF(PRESENT(k)) THEN
      q=k-ke
    ELSE
      q=-ke
    ENDIF

    CALL cartez2spher( q ,kem ,thetae ,phie)
    a = kem**2 +e**2
    kec = CMPLX( 0._RP ,kem ,KIND=RP )

    tpw = 0.
    DO j = 0 ,(n-l)/2
      tpw = tpw +(-0.25_RP*a/e**2)**j *fac(n-j)/( fac(j) *fac(n-l-2*j) )
    END DO


    tpw = tpw *SQRT(2./pi) *norm_fac(e,n) *fac(n-l) *(kec/e)**l *(2.*e/a)**n /a &
      *spherical_harmonic(l,m,thetae,phie)

  END FUNCTION tpw

  PURE COMPLEX(KIND=RP) FUNCTION tcw( n ,l ,m ,e ,ke ,k)
    INTEGER       ,INTENT(IN) :: n ,l ,m
    REAL(KIND=RP) ,INTENT(IN) :: e ,ke(3) ,k(3)

    REAL(KIND=RP)    :: kem ,km ,alpha ,a ,aj1 ,ke_t(3)
    COMPLEX(KIND=RP) :: w ,kec ,ekec ,gam(0:n) ,f21(0:n,0:n) ,alphac ,w1m ,kep ,kp
    COMPLEX(KIND=RP) :: tmp_j ,tmp_j1 ,tmp_m1 ,cst_j ,cst_j1 ,cst_m1
    REAL(KIND=RP)    :: tmp_s ,tmp_s3 ,tmp_s1 ,cst_s ,cst_s3 ,cst_s1 ,cst_s2
    COMPLEX(KIND=RP) ,PARAMETER :: zi = (0._RP ,1._RP)
    INTEGER          :: j ,j1 ,ma ,m1 ,s ,s1 ,s2 ,s3

    ma = ABS(m)
    kem = NORM2(ke)
    ke_t=-ke
    km = NORM2(k)
    IF( m>=0 ) THEN
      kp = CMPLX( k(1) ,k(2) ,KIND=RP )
      kep = CMPLX( ke_t(1) ,ke_t(2) ,KIND=RP )
    ELSE
      kp = CMPLX( k(1) ,-k(2) ,KIND=RP )
      kep = CMPLX( ke_t(1) ,-ke_t(2) ,KIND=RP )
    END IF

    alpha  = 1._RP / kem
    alphac = CMPLX( 0._RP ,alpha ,KIND=RP) ! i*\alpha
    kec    = CMPLX( 0._RP ,kem   ,KIND=RP) ! i*ke
    ekec   = CMPLX( e     ,-kem  ,KIND=RP) ! (\epsilon-ike)
    a      = NORM2(k+ke_t)**2 +e**2
    w      = CMPLX( NORM2(k+ke_t)**2 -km**2 +kem**2    ,2.*e*kem  ,KIND=RP) /a ! w = b/a
    w1m    = CMPLX( e**2 +km**2 -kem**2 ,-2.*e*kem ,KIND=RP) /a ! (1-w)

    gam(0) = 1.
    DO j1 = 1 ,n
      gam(j1) = gam(j1-1) *CMPLX( j1 ,-alpha ,KIND=RP)
    END DO

    DO CONCURRENT( j1= 0:n )
      aj1 = j1 +1._RP
      f21(j1 ,0) = 1._RP
      f21(j1 ,1) = 1. + alphac*w / ( aj1*w1m )
      IF( n<2 .or. j1>n-2 ) CYCLE
      DO j = 1 ,n -j1 -1
        f21( j1 , j +1 ) = (1. + ( j +alphac*w) / ( (aj1 +j)*w1m ) ) *f21(j1 ,j) &
          - j/( (aj1 +j)*w1m ) *f21(j1 ,j-1)
      END DO
    END DO

    tcw = 0.
    cst_j = -0.25*a/ekec**2
    cst_j1 = kec/ekec
    cst_m1 = kep/kp
    cst_s = -(km/k(3))**2
    cst_s3 = ke_t(3)/k(3)
    cst_s1 = kem/km
    cst_s2 = 2.*DOT_PRODUCT(k,ke_t)
    DO CONCURRENT ( j = 0:(n-l)/2 )
      tmp_j = fac(n-j) /fac(j) *cst_j**j!tmp_j = fac(n-j)/fac(j)*(-0.25*a/ekec**2)**j
      DO CONCURRENT( j1 = 0:(n-l-2*j) )
        !tmp_j1 = (kec/ekec)**j1 /( fac(j1) *fac(n-l-2*j-j1) ) *tmp_j
        tmp_j1 = cst_j1**j1 /( fac(j1) *fac(n-l-2*j-j1) ) *tmp_j
        DO CONCURRENT(m1 = 0:ma)
          tmp_m1 = cst_m1**m1 /( fac(ma-m1) *fac(m1) ) *tmp_j1
          !tmp_m1 = (kep/kp)**m1 /( fac(ma-m1) *fac(m1) ) *tmp_j1
          DO s = 0,(l-ma)/2
            tmp_s = cst_s**s /( fac(l-s) )!tmp_s = (-1)**s *(km/k(3))**(2*s) /( fac(l-s) )
            DO s3 = 0,l-ma-2*s
              tmp_s3 = cst_s3**s3 /( fac(s3) *fac(l-ma-2*s-s3) ) *tmp_s
              !tmp_s3 = (ke_t(3)/k(3))**s3 /( fac(s3) *fac(l-ma-2*s-s3) ) *tmp_s
              DO s1 = 0,s
                tmp_s1 = cst_s1**s1 /fac(s-s1) *tmp_s3!tmp_s1 = (kem/km)**s1 /fac(s-s1) *tmp_s3
                DO CONCURRENT( s2=0:s1)
                  tcw = tcw + gam(m1+s3+2*s1-s2+j1) *f21(m1+s3+2*s1-s2+j1 ,n -j -(m1+s3+2*s1-s2+j1) ) &
                    /( fac(s2) *fac(s1-s2) *fac(m1+s3+2*s1-s2+j1) ) &
                    *tmp_m1 *tmp_s1 *cst_s2**s2
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    tcw = tcw *norm_fac(e ,n ) *SQRT(l+0.5_RP) *fac(n-l) *SQRT( fac(l-ma)/fac(l+ma) ) *(-1)**ma &
      *fac(ma) *(zi*k(3)/ekec)**l *(2.*ekec/a)**n/a *(kp/k(3))**ma *w1m**(-alphac) /pi
    !*powcc(w1m,-alpha)/pi
    IF( MOD(m,2)<0 ) tcw = -tcw

  END FUNCTION tcw

  PURE COMPLEX(KIND=RP) FUNCTION tcw0( n ,l ,m ,e ,ke)
    INTEGER       ,INTENT(IN) :: n ,l ,m
    REAL(KIND=RP) ,INTENT(IN) :: e ,ke(3)

    REAL(KIND=RP)    :: kem ,alpha ,thetae ,phie ,a ,aj1
    COMPLEX(KIND=RP) :: w ,kec ,ekec ,gam(0:n) ,f21(0:n,0:(n-l)) ,alphac ,w1m ,tmp
    INTEGER          :: j ,j1

    CALL cartez2spher( -ke ,kem ,thetae ,phie)

    alpha  = 1._RP / kem
    alphac = CMPLX( 0._RP ,alpha ,KIND=RP) ! i*\alpha
    kec    = CMPLX( 0._RP ,kem   ,KIND=RP) ! i*ke
    ekec   = CMPLX( e     ,-kem  ,KIND=RP) ! (\epsilon-ike)
    a      = kem**2 +e**2
    w      = CMPLX( 2.*kem**2    ,2.*e*kem  ,KIND=RP) /a ! w = b/a
    w1m    = CMPLX( e**2 -kem**2 ,-2.*e*kem ,KIND=RP) /a ! (1-w)

    gam(0) = 1.
    DO j1 = 1 ,n
      gam(j1) = gam(j1-1) *CMPLX( j1 ,-alpha ,KIND=RP)
    END DO

    DO CONCURRENT( j1= 0:n-l )
      aj1 = l +j1 +1._RP
      f21(j1 ,0) = 1._RP
      f21(j1 ,1) = 1. + alphac*w / ( aj1*w1m )
      IF(n-l<2) CYCLE
      DO j = 1 ,n -l -1
        f21( j1 , j +1 ) = (1. + ( j +alphac*w) / ( (aj1 +j)*w1m ) ) *f21(j1 ,j) &
          - j/( (aj1 +j)*w1m ) *f21(j1 ,j-1)
      END DO
    END DO

    tcw0 = 0.
    DO CONCURRENT ( j = 0:(n-l)/2 )
      tmp = fac(n-j) /fac(j) *(-0.25*a/ekec**2)**j
      DO CONCURRENT( j1 = 0:(n-l-2*j) )!, j1+2*j <= n-l )
        tcw0 = tcw0 + 1./( fac(j1) *fac(n-l-2*j-j1) *fac(l+j1) ) *(kec/ekec)**j1 *gam(l+j1) &
          *f21(j1 ,n -l -j -j1 ) *tmp !*fac(n-j)/ fac(j) *(-0.25*a/ekec**2)**j
      END DO
    END DO

    tcw0 = tcw0 *norm_fac(e ,n ) *SQRT(2./pi) *fac(n-l) *(kec/ekec)**l *(2.*ekec/a)**n/a &
      *spherical_harmonic(l ,m ,thetae ,phie) *w1m**(-alphac)!*powcc(w1m,-alpha)

  END FUNCTION tcw0

  ELEMENTAL COMPLEX(KIND=rp) FUNCTION powcc(z1,y2)
    COMPLEX(KIND=rp),INTENT(IN) :: z1
    REAL(KIND=rp),INTENT(IN) :: y2
    REAL(KIND=rp) :: theta,zm

    theta = ATAN2( AIMAG(z1) ,REAL(z1 ,KIND=rp) )
    zm = LOG( ABS(z1) )
    powcc = EXP(-y2*theta) *CMPLX( COS(y2*zm) ,SIN(y2*zm) ,KIND=RP )

  END FUNCTION powcc

END SUBMODULE fdcs_e2e
