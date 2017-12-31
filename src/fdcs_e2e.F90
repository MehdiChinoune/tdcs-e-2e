MODULE fdcs_e2e
  USE constants ,ONLY: RP
  IMPLICIT NONE
CONTAINS

  SUBROUTINE fdcs_fba_pw(in_unit,out_unit)
    USE constants ,ONLY: pi, ev, deg
    USE special_functions ,ONLY: factorial
    USE input ,ONLY: read_input, read_orbit
    USE trigo ,ONLY: spher2cartez, nrm2

    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom, Orbit
    REAL(RP) :: Ei, Es, Ee
    REAL(RP) :: thetas
    INTEGER :: step(3), exchange
    INTEGER :: nelec
    INTEGER :: lo, no
    REAL(RP), ALLOCATABLE :: a(:), e(:)
    INTEGER, ALLOCATABLE :: n(:)

    REAL(RP), PARAMETER   :: phie = 0._RP
    REAL(RP) :: kim, ksm, kem, km
    REAL(RP) :: ki(3), ks(3), ke(3), k(3)

    REAL(RP) :: factor, sigma
    COMPLEX(RP) :: D_term
    INTEGER :: i, io, mo

    REAL(rp) :: k2(3), k2m = 0._RP
    COMPLEX(rp) :: E_term = (0._RP, 0._RP)

    CALL factorial()

    CALL read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange)

    CALL read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    CALL spher2cartez( kim, 0._RP, 0._RP, ki )
    CALL spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = nrm2(k)

    factor = nelec*4._RP*ksm*kem/(kim)

    DO i = step(1), step(2), step(3)

      CALL spher2cartez( kem, i*deg, phie, ke )

      IF(exchange==1) THEN
        k2 = ki -ke
        k2m = nrm2(k2)
      END IF

      sigma = 0._RP
      DO mo = 0, lo

        D_term = (0._RP,0._RP)
        DO io = 1, no
          D_term = D_term +a(io)*( tpw(n(io), lo, mo, e(io), ke, k ) &
            -tpw(n(io), lo, mo, e(io), ke ) )
        END DO

        IF(exchange==1) THEN
          E_term = (0._RP,0._RP)
          DO io = 1, no
            E_term = E_term &
              +a(io)*( tpw(n(io), lo, mo, e(io), ks, k2 ) &
                -tpw(n(io), lo, mo, e(io), ks ) )
          END DO
          sigma = sigma +(1+mo)*( ABS(D_term/km**2)**2 +ABS(E_term/k2m**2)**2 &
            -REAL( D_term*CONJG(E_term)/(km**2*k2m**2 ), rp ) )
        ELSE
          sigma = sigma +(1+mo)*ABS(D_term/km**2)**2
        END IF

      END DO

      sigma = factor*sigma

      PRINT '(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    END DO

  END SUBROUTINE fdcs_fba_pw

  SUBROUTINE fdcs_fba_cw(in_unit,out_unit)
    USE constants ,ONLY: ev, deg, pi
    USE trigo ,ONLY: spher2cartez, nrm2
    USE special_functions ,ONLY: factorial
    USE input ,ONLY: read_input, read_orbit
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom, Orbit
    REAL(RP) :: Ei, Es, Ee
    REAL(RP) :: thetas
    INTEGER :: step(3), exchange
    INTEGER :: nelec, ze, zs
    INTEGER :: lo, no
    REAL(RP), ALLOCATABLE :: a(:), e(:)
    INTEGER, ALLOCATABLE :: n(:)

    REAL(RP), PARAMETER   :: phie = 0._RP
    REAL(RP) :: kim, ksm, kem, km, alpha
    REAL(RP) :: ki(3), ks(3), ke(3), k(3)

    REAL(RP) :: factor, sigma
    COMPLEX(RP) :: term
    INTEGER :: i, io, mo

    CALL factorial()

    CALL read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange)

    CALL read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    zs = -1 ! Charge of the projectile
    kem = SQRT(2.*Ee*eV)
    ze = -1 ! Charge of ejected particle
    alpha = -ze/kem

    CALL spher2cartez( kim, 0._RP, 0._RP, ki )
    CALL spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = nrm2(k)

    factor = nelec*4._RP*ksm*kem / (kim*km**4)
    !\abs{ \exp{\pi\alpha/2}*\Gamma(1-i\alpha) }^2
    if(ze/=0) factor = factor*2.*pi*alpha/(1._RP-EXP(-2.*pi*alpha))

    DO i = step(1),step(2),step(3)

      CALL spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      DO mo = 0, lo
        term = (0._RP,0._RP)
        DO io = 1, no
          term = term +a(io)*( tcw(n(io), lo, mo, e(io), alpha, ke, k ) &
            +zs*tcw0(n(io), lo, mo, e(io), alpha, ke ) )
        END DO
        sigma = sigma +(mo+1)*ABS(term)**2
      END DO

      sigma = factor*sigma/(2*lo+1)

      PRINT '(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    END DO

  END SUBROUTINE fdcs_fba_cw

  SUBROUTINE fdcs_fba_dw(in_unit,out_unit)
    USE ieee_arithmetic ,only: ieee_is_nan
    USE constants ,ONLY: pi, deg, ev
    USE special_functions ,ONLY: cgamma, spherical_harmonic, ricbes, factorial!, coul90, symbol_3j
    USE utils ,ONLY: norm_fac, y1y2y3, calculate_U, ode_second_dw
    USE input ,ONLY: read_input, read_orbit
    USE trigo ,ONLY: spher2cartez, cartez2spher, nrm2

    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom, Orbit
    REAL(RP) :: Ei, Es, Ee
    REAL(RP) :: thetas
    INTEGER :: step(3), exchange
    INTEGER :: nelec
    INTEGER :: lo, no
    REAL(RP), ALLOCATABLE :: a(:), e(:)
    INTEGER, ALLOCATABLE :: n(:)

    REAL(RP), PARAMETER   :: phie = 0._RP
    REAL(RP) :: kim, ksm, kem, km
    REAL(RP) :: theta, phi
    REAL(RP) :: ki(3), ks(3), ke(3), k(3)

    INTEGER, PARAMETER :: lmax = 47, lemax = 15, np = 2400
    REAL(RP) :: x(0:np), wf(0:np)
    REAL(RP), ALLOCATABLE :: chi_b(:,:), chi_0a(:,:)
    REAL(RP) :: sigma_le(0:lemax), delta(0:lemax), integral(0:lemax,0:lmax)
    REAL(RP) :: etae, rc, h
    COMPLEX(RP) :: tmp_z
    COMPLEX(RP), PARAMETER :: zi = (0._RP, 1._RP)
    INTEGER :: le, l, me, ze

    REAL(RP) :: factor, sigma
    COMPLEX(RP) :: term, term0
    INTEGER :: i, io, mo

    REAL(RP) :: gc_0a(0:lmax), fdc_0a(0:lmax), gdc_0a(0:lmax), fc_0a(0:lmax)
    REAL(RP), ALLOCATABLE :: U(:,:), U_tmp(:)

    CALL factorial()

    CALL read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange )

    CALL read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)
    ze = 1

    CALL spher2cartez( kim, 0._RP, 0._RP, ki )
    CALL spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    CALL cartez2spher( k, km, theta, phi)

    factor = nelec*4._RP*ksm*kem / (kim*km**4)
    factor = factor*2./(pi*kem**2)

    rc = 10._RP
    h = rc/np

    etae = -ze/kem

    x = [(i*h, i = 0, np) ]
    ALLOCATE( chi_b(0:np,0:lemax), chi_0a(0:np,0:lmax) )
    chi_b(0,:) = 0._RP
    chi_0a(0,:) = 0._RP
    DO i = 1,np
      CALL ricbes(km*x(i), lmax, fc_0a, gc_0a, fdc_0a, gdc_0a, le)
      chi_0a(i,:) = fc_0a/(km*x(i))
    END DO

    WHERE(ieee_is_nan(chi_0a) )
      chi_0a = 0._RP
    END WHERE

    ALLOCATE( U(0:np,0:lemax), U_tmp(0:np) )
    CALL calculate_U(Atom, Orbit, x, U_tmp, 1 )

    IF(ze/=0) THEN
      U(0,:) = -HUGE(1._RP)
    ELSE
      U(0,:) = -kem**2
    END IF

    U(1:np,lemax) = -kem**2 -2.*( ze +U_tmp(1:np) ) /x(1:np)
    DO l = 0,lemax
      U(1:np,l) = l*(l+1)/x(1:np)**2 +U(1:np,lemax)
    END DO

    CALL ode_second_dw(kem, lemax, rc, ze, U, chi_b, delta)

    DEALLOCATE(U,U_tmp)

    wf = 0._RP
    DO io = 1, no
      wf = wf +a(io)*norm_fac(e(io), n(io) )*x**n(io)*EXP(-e(io)*x)
    END DO

    DO le = 0,lemax
      tmp_z = cgamma( CMPLX(le+1._RP, etae, KIND=RP), 1 )
      tmp_z = exp( zi*aimag(tmp_z) )
      sigma_le(le) = ATAN2( AIMAG(tmp_z), REAL(tmp_z,KIND=RP) )
    END DO

    integral = 0._RP
    DO l = 0,lmax
      DO le = 0,lemax
        IF( MOD(le+l+lo,2)/=0 ) CYCLE
        integral(le, l) = 0.5*h*SUM( wf(1:np)*chi_b(1:np,le)*chi_0a(1:np,l) &
          +wf(0:np-1)*chi_b(0:np-1,le)*chi_0a(0:np-1,l) )
      END DO
    END DO

    term0 = zi**(-lo)*EXP( CMPLX(0._RP, -(sigma_le(lo)+delta(lo)), KIND=RP ) ) &
      *0.5*h*SUM( wf(1:np)*chi_b(1:np,lo) +wf(0:np-1)*chi_b(0:np-1,lo) )

    DO i = step(1), step(2), step(3)

      CALL spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      DO mo = 0, lo

        term = (0._RP, 0._RP)
        DO l = 0,lmax
          DO le = 0,lemax
            IF( MOD(le+l+lo,2)/=0 ) CYCLE
            tmp_z = (0._RP, 0._RP)
            DO me = -le,le
              IF( ABS(me-mo)>l ) CYCLE
              tmp_z = tmp_z +y1y2y3(le, l, lo, -me, me-mo, mo) &
                *spherical_harmonic(le, me, i*deg, phie) &
                *spherical_harmonic(l, mo-me, theta, phi)
            END DO
            term = term +zi**(l-le)*integral(le, l)*tmp_z &
              *EXP( CMPLX(0._RP, -(sigma_le(le)+delta(le)), KIND=RP ) )

          END DO
        END DO
        term = 4.*pi*term

        sigma = sigma +(mo+1)*ABS( term &
          -spherical_harmonic(lo, mo, i*deg, phie)*term0 )**2

      END DO

      sigma = factor*sigma/(2*lo+1)

      PRINT '(1x,i4,1x,es15.8)',i,sigma
      WRITE(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    END DO

  END SUBROUTINE fdcs_fba_dw

  SUBROUTINE fdcs_dwb(in_unit,out_unit)
    USE constants ,ONLY: ev, deg, pi
    USE special_functions ,ONLY: spherical_harmonic, cgamma, factorial
    USE utils ,ONLY: norm_fac, calculate_U
    USE input ,ONLY: read_input, read_orbit
    USE trigo ,ONLY: spher2cartez, nrm2
    USE integration ,ONLY: gauleg

    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom, Orbit
    REAL(RP) :: Ei, Es, Ee
    REAL(RP) :: thetas
    INTEGER :: step(3), exchange, PCI = 2
    INTEGER :: nelec
    INTEGER :: lo, no
    REAL(RP), ALLOCATABLE :: a(:), e(:)
    INTEGER, ALLOCATABLE :: n(:)

    REAL(RP), PARAMETER   :: phie = 0._RP
    REAL(RP) :: kim, ksm, kem

    INTEGER, PARAMETER :: nx = 31*64, nr = 48000

    INTEGER, PARAMETER :: limax = 95
    REAL(RP), ALLOCATABLE, TARGET :: chi_0(:,:)
    REAL(RP), TARGET :: delta_li(0:limax)

    INTEGER, PARAMETER :: lsmax = 95
    REAL(RP), ALLOCATABLE, TARGET :: chi_a(:,:)
    REAL(RP) :: sigma_ls(0:lsmax), etas
    REAL(RP), TARGET :: delta_ls(0:lsmax)
    COMPLEX(RP) :: ylms(0:lsmax,-lsmax:lsmax)
    INTEGER :: ls, zs = 1

    INTEGER, PARAMETER :: lemax = 23
    REAL(RP), ALLOCATABLE, TARGET :: chi_b(:,:)
    REAL(RP) :: sigma_le(0:lemax), etae
    REAL(RP), TARGET :: delta_le(0:lemax)
    COMPLEX(RP), ALLOCATABLE :: ylme(:,:,:)
    INTEGER :: le, me, ze = 1

    REAL(RP), ALLOCATABLE :: r(:), wf(:)
    REAL(RP) :: rc, h
    COMPLEX(RP), ALLOCATABLE :: integral(:,:,:,:), integralx(:,:,:,:)
    COMPLEX(RP) :: tmp_z

    REAL(RP) :: factor, sigma
    COMPLEX(RP) :: termd, termx
    INTEGER :: i, io, mo = 0

    REAL(RP), ALLOCATABLE :: x(:), w(:)

    REAL(RP), ALLOCATABLE :: U_tmp(:)
    REAL(RP), ALLOCATABLE :: xt(:), wt(:)
    INTEGER, PARAMETER :: nc = nx/64

    CALL factorial()

    CALL read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange )

    CALL read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    etas = -zs/ksm
    etae = -ze/kem

    factor = nelec* (2.*pi)**4*ksm*kem / kim
    factor = factor*2._RP/pi**4 /(kem*ksm*kim)**2

    rc = 250._RP

    ALLOCATE(x(nx), w(nx) )

      DO i = 1,nc
        CALL gauleg( (i-1)*rc/nc, i*rc/nc, xt, wt, 64 )
        x( (i-1)*64+1: i*64 ) = xt
        w( (i-1)*64+1: i*64 ) = wt
      END DO
    DEALLOCATE( xt, wt )

    ALLOCATE(chi_0(nx,0:limax), chi_a(nx,0:lsmax), chi_b(nx,0:lemax) )

    h = rc/nr
    ALLOCATE( r(0:nr) )
    r = [(i*h,i = 0,nr)]

      ALLOCATE( U_tmp(0:nr) )

      CALL calculate_U(Atom, Orbit, r, U_tmp, 0)
      CALL calculate_chi( kim, r, U_tmp, 0, x, chi_0, delta_li )

      CALL calculate_U(Atom, Orbit, r, U_tmp, 1)
      CALL calculate_chi( ksm, r, U_tmp, zs, x, chi_a, delta_ls )
      CALL calculate_chi( kem, r, U_tmp, ze, x, chi_b, delta_le )

    DEALLOCATE( U_tmp )

    ALLOCATE(wf(nx))
    wf = 0._RP
    DO io = 1, no
      wf = wf +a(io)*norm_fac(e(io), n(io) )*x**n(io)*EXP(-e(io)*x)
    END DO

    chi_b(:,lo) = chi_b(:,lo) -wf*sum( w*chi_b(:,lo)*wf )
    chi_a(:,lo) = chi_a(:,lo) -wf*sum( w*chi_a(:,lo)*wf )

    tmp_z = cgamma( CMPLX( 1._RP, etae, RP), 0 )
    sigma_le(0) = ATAN2( AIMAG(tmp_z), REAL(tmp_z, RP) )
    DO le = 0,lemax-1
      sigma_le(le+1) = sigma_le(le) +atan2( etae, le+1._RP )
    END DO
    sigma_le = sigma_le +delta_le

    tmp_z = cgamma( CMPLX( 1._RP, etas, RP), 0 )
    sigma_ls(0) = ATAN2( AIMAG(tmp_z), REAL(tmp_z, RP) )
    DO ls = 0,lsmax-1
      sigma_ls(ls+1) = sigma_ls(ls) +atan2( etas, ls+1._RP )
    END DO
    sigma_ls = sigma_ls +delta_ls

    CALL dwb_integrals(chi_0, chi_a, chi_b, delta_li, sigma_ls, sigma_le &
      , wf, x, w, lo, integral)
    IF(exchange==1) THEN
      CALL dwb_integrals(chi_0, chi_b, chi_a, delta_li, sigma_le, sigma_ls &
        , wf, x, w, lo, integralx)
    END IF

    ylms = (0._RP, 0._RP)
    DO ls = 0,lsmax
      ylms(ls,0) = spherical_harmonic(ls,0,thetas*deg,pi)
      DO me = 1,ls
        ylms(ls,me) = spherical_harmonic(ls,me,thetas*deg,pi)
      END DO
      ylms(ls,-2:-ls:-2) = ylms(ls,2:ls:2)
      ylms(ls,-1:-ls:-2) = -ylms(ls,1:ls:2)
    END DO

    ALLOCATE(ylme(0:lemax,-lemax:lemax,step(1)/step(3):step(2)/step(3)))
    ylme = (0._RP, 0._RP)
    DO i = lbound(ylme,3),ubound(ylme,3)
      DO le = 0,lemax
        ylme(le,0,i) = spherical_harmonic(le,0,i*step(3)*deg,phie)
        DO me = 1,le
          ylme(le,me,i) = spherical_harmonic(le,me,i*step(3)*deg,phie)
        END DO
        ylme(le,-2:-le:-2,i) = ylme(le,2:le:2,i)
        ylme(le,-1:-le:-2,i) = -ylme(le,1:le:2,i)
      END DO
    END DO
!!$OMP PARALLEL DO PRIVATE(sigma, mo, termd, termx, le, ls, me)
    DO i = step(1), step(2), step(3)

      sigma = 0._RP
      DO mo = 0, lo

        termd = (0._RP, 0._RP)
        DO le = 0,lemax
          DO ls = 0,lsmax
            DO me = -le,le
              IF( ABS(me+mo)>ls ) CYCLE
              termd = termd +ylme(le,me,i/step(3))*ylms(ls,mo+me) &
                *integral(ls, le, me, mo )
            END DO
          END DO
        END DO

        termx = (0._RP, 0._RP)
        IF(exchange==1) THEN
          DO ls = 0,lemax
            DO le = 0,lsmax
              DO me = -le,le
                IF( ABS(me+mo)>ls ) CYCLE
                termx = termx +ylms(le,me)*ylme(ls,mo+me,i/step(3)) &
                  *integralx(ls, le, me, mo )
              END DO
            END DO
          END DO
        END IF

        sigma = sigma +(mo+1)*( ABS(termd)**2 +ABS(termx)**2 &
          -REAL( CONJG(termd)*termx, RP ) )
      END DO

      sigma = factor*sigma/(2*lo+1)
      IF(PCI>=1) THEN
        CALL PCI_EFFECTS(i,sigma)
      END IF

      PRINT '(1x,i4,1x,es15.8)', i, sigma
      WRITE(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    END DO

  CONTAINS

    ELEMENTAL SUBROUTINE PCI_EFFECTS(i,sigma)
      USE special_functions ,ONLY: conhyp_opt
      INTEGER, INTENT(IN) :: i
      REAL(RP), INTENT(INOUT) :: sigma
      REAL(RP) :: kesm
      REAL(RP) :: r12_av, Et
      kesm = SQRT( kem**2 +ksm**2 -2*kem*ksm*COS( (thetas+i)*deg ) )/2._RP
      sigma = pi/(kesm*(EXP(pi/kesm) -1._RP ) )*sigma
      IF(PCI==2) THEN
        Et = (Es +Ee)*eV
        r12_av = pi**2/(16.*Et)*( 1._RP +(0.627/pi)*SQRT(Et)*LOG(Et) )**2
        sigma = ABS( conhyp_opt( 0.5/kesm, -2.*kesm*r12_av  ) )**2*sigma
      END IF
    END SUBROUTINE

  END SUBROUTINE fdcs_dwb

  SUBROUTINE fdcs_bbk(in_unit,out_unit)
    USE constants ,ONLY: pi, ev, deg
    USE special_functions ,ONLY: factorial
    USE input ,ONLY: read_input, read_orbit
    USE trigo ,ONLY: spher2cartez, nrm2

    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit

    CHARACTER(LEN=2) :: Atom, Orbit
    REAL(RP) :: Ei, Es, Ee
    REAL(RP) :: thetas
    INTEGER :: step(3)
    INTEGER :: nelec
    INTEGER :: lo, no
    REAL(RP), ALLOCATABLE :: a(:), e(:)
    INTEGER, ALLOCATABLE :: n(:)

    REAL(RP), PARAMETER   :: phie = 0._RP
    REAL(RP) :: kim, ksm, kem, km
    REAL(RP) :: ki(3), ks(3), ke(3), k(3)

    REAL(RP) :: factor, sigma
    COMPLEX(RP) :: D_term
    INTEGER :: i, io, mo

    CALL factorial()

    CALL read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit)

    CALL read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = SQRT(2.*Ei*eV)
    ksm = SQRT(2.*Es*eV)
    kem = SQRT(2.*Ee*eV)

    CALL spher2cartez( kim, 0._RP, 0._RP, ki )
    CALL spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = nrm2(k)

    factor = nelec*4._RP*ksm*kem/(kim)

    DO i = step(1), step(2), step(3)

      CALL spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      DO mo = 0, lo

        D_term = (0._RP,0._RP)
        DO io = 1, no
          D_term = D_term +a(io)*( tpw(n(io), lo, mo, e(io), ke, k ) &
            -tpw(n(io), lo, mo, e(io), ke ) )
        END DO
        sigma = sigma +(1+mo)*ABS(D_term/km**2)**2
      END DO

      sigma = factor*sigma

      PRINT '(1x,i4,1x,es15.8)',i,sigma
      WRITE( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    END DO

  END SUBROUTINE fdcs_bbk

  COMPLEX(RP) FUNCTION U_bbk(alpha1, alpha2, alpha3, k1, k2, k3, lam1, lam2, lam3, p1, p2)
    USE constants ,ONLY: pi
    USE integration ,ONLY: gauleg
    USE trigo ,ONLY: nrm2

    REAL(RP), INTENT(IN) :: alpha1, alpha2, alpha3
    REAL(RP), INTENT(IN) :: k1(3), k2(3), k3(3)
    REAL(RP), INTENT(IN) :: lam1, lam2, lam3
    REAL(RP), INTENT(IN) :: p1(3), p2(3)

    COMPLEX(RP), PARAMETER :: zi = (0._RP, 1._RP)

    REAL(RP),ALLOCATABLE :: t3(:), y(:), wy(:), s(:), ws(:)
    REAL(RP) :: p11(3), p22(3), q1(3), q2(3)
    REAL(RP) :: q1m, q2m
    REAL(RP) :: k1m, k2m, k3m
    COMPLEX(RP) :: lam33
    COMPLEX(RP) :: sig0, sig1, sig2, sig12
    COMPLEX(RP) :: sig0_a(3), sig1_a(3), sig2_a(3), sig12_a(3)

    INTEGER :: i, j
    COMPLEX(RP) :: N_t3, Integ

    U_bbk = ( 0._RP, 0._RP)

    k1m = nrm2(k1)
    k2m = nrm2(k2)
    k3m = nrm2(k3)

    CALL gauleg(-10._RP, 10._RP, y, wy, 64)
    t3 = 1._RP /(1._RP+ EXP(y) )

    CALL gauleg(0._RP, 10._RP, s, ws, 64)

    Integ = (0._RP, 0._RP)
    DO i = 1,64
      p11 = p1 -t3(i)*k3
      q1 = k1 +p11
      q1m = nrm2(q1)
      p22 = p2 -t3(i)*k3
      q2 = k2 +p22
      q2m = nrm2(q2)
      lam33 = CMPLX(lam3, -t3(i)*k3m, RP )

      sig0_a = [ CMPLX( (lam1+lam2)**2+nrm2(q1-q2)**2, 0._RP, RP) &
        , 2*(lam2*(lam1**2+lam33**2+q1m**2) +lam1*(lam2**2+lam33**2+q2m**2) ) &
        , ( (lam1+lam33)**2+q1m**2)*((lam2+lam33)**2+q2m**2) ]
      sig1_a = -2*[ DOT_PRODUCT(q1-q2,k1) +zi*(lam1+lam2)*k1m &
        , 2*lam2*DOT_PRODUCT(q1,k1) +zi*( (lam2**2+lam33**2+q2m**2)*k1m+2*lam1*lam2*k1m ) &
        , ((lam2+lam33)**2 +q2m**2)*( DOT_PRODUCT(q1,k1) +zi*(lam1+lam33)*k1m ) ]
      sig2_a = -2*[ DOT_PRODUCT(q2-q1,k2) +zi*(lam1+lam2)*k2m &
        , 2*lam1*DOT_PRODUCT(q2,k2) +zi*( (lam1**2+lam33**2+q1m**2)*k2m+2*lam1*lam2*k2m ) &
        , ((lam1+lam33)**2 +q1m**2)*( DOT_PRODUCT(q2,k2) +zi*(lam2+lam33)*k2m ) ]
      sig12_a = -2*[ CMPLX( k1m*k2m +DOT_PRODUCT(k1,k2), 0._RP, RP) &
        , 2*( (lam1*k1m-zi*DOT_PRODUCT(q1,k1))*k2m +(lam2*k2m-zi*DOT_PRODUCT(q2,k2))*k1m) &
        , -2*( DOT_PRODUCT(q1,k1)+zi*(lam1+lam33)*k1m) &
          *( DOT_PRODUCT(q2,k2)+zi*(lam2+lam33)*k2m) ]

      N_t3 = (0._RP, 0._RP)
      DO j = 1,64
        sig0 = sig0_a(1)*(s(j)+2*lam33)*s(j) +sig0_a(2)*s(j) +sig0_a(3)
        sig1 = sig1_a(1)*(s(j)+2*lam33)*s(j) +sig1_a(2)*s(j) +sig1_a(3)
        sig2 = sig2_a(1)*(s(j)+2*lam33)*s(j) +sig2_a(2)*s(j) +sig2_a(3)
        sig12 = sig12_a(1)*(s(j)+2*lam33)*s(j) +sig12_a(2)*s(j) +sig12_a(3)
        N_t3 = (sig0/(sig0+sig1) )**(-zi*alpha1)*(sig0/(sig0+sig2) )**(-zi*alpha2)/sig0 &
          +N_t3
      END DO
      Integ = Integ +SINH(pi*alpha3)/(2*pi*zi)*(4*pi)**2*N_t3
    END DO


  END FUNCTION

  SUBROUTINE calculate_chi( km, r, U_tmp, z, x, chi, delta )
    USE utils ,ONLY: ode_second_dw, intrpl
    REAL(RP), INTENT(IN) :: U_tmp(0:), r(0:), x(:), km
    INTEGER, INTENT(IN) :: z
    REAL(RP), INTENT(OUT) :: chi(:,0:), delta(0:)

    REAL(RP), ALLOCATABLE :: U(:,:), chi_tmp(:,:)
    REAL(RP) :: rc
    INTEGER :: l, lmax, nr

    lmax = size(delta) -1
    nr = size(r) -1
    rc = r(nr)

    ALLOCATE( U(0:nr,0:lmax), chi_tmp(0:nr,0:lmax) )

    IF(z==0) THEN
      U(0,0:lmax) = -km**2
    ELSE
      U(0,0:lmax) = -HUGE(1._RP)
    END IF

    U(1:nr,lmax) = -km**2 -2.*( z*1._RP +U_tmp(1:nr) )/r(1:nr)
    DO l = 0,lmax
      U(1:nr,l) = l*(l+1)/r(1:nr)**2 +U(1:nr,lmax)
    END DO

    CALL ode_second_dw(km, lmax, rc, z, U, chi_tmp, delta )

    DO l = 0,lmax
      CALL INTRPL(r, chi_tmp(:,l), x, chi(:,l))
    END DO

  END SUBROUTINE calculate_chi

  SUBROUTINE dwb_integrals( chi_0, chi_a, chi_b, sig_0, sig_a, sig_b, wf, x, w, lo &
    , integral )
    USE special_functions ,ONLY: symbol_3j
    REAL(RP), INTENT(IN) :: chi_0(:,0:), chi_a(:,0:), chi_b(:,0:), wf(:)
    REAL(RP), INTENT(IN) :: x(:), w(:)
    REAL(RP), INTENT(IN) :: sig_0(0:), sig_a(0:), sig_b(0:)
    INTEGER, INTENT(IN) :: lo
    COMPLEX(RP), ALLOCATABLE, INTENT(OUT) :: integral(:,:,:,:)

    REAL(RP), PARAMETER :: VSmall = TINY(1._RP), VBig = HUGE(1._RP)
    COMPLEX(RP), PARAMETER :: zi = (0._RP, 1._RP)
    INTEGER :: limax, lsmax, lemax, nx, nmin
    COMPLEX(RP), ALLOCATABLE :: integ(:,:)
    INTEGER :: li, ls, le, l, me, mo
    ! Block variables
    REAL(RP) :: Ti, Si
    REAL(RP) :: integ0
    REAL(RP) :: xil, xil1
    INTEGER :: i, is, nmax

    nx = size(x)
    limax = size(sig_0) -1
    lsmax = size(sig_a) -1
    lemax = size(sig_b) -1

    nmin = count(x<=25._RP)

    ALLOCATE( integral(0:lsmax, 0:lemax, -lemax:lemax, 0:lo) )
    ALLOCATE( integ( 0:MIN(limax+lsmax,lemax+lo),  0:limax ) )
    integral = 0._RP
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(integ,li,l,me,mo ,ti,si,integ0,xil,xil1,i,is,nmax)
    DO ls = 0,lsmax; DO le = 0,lemax
      integ = 0._RP
      DO li = 0,limax
        DO l = MAX(ABS(ls-li), ABS(le-lo)), MIN(ls+li, le+lo)
          IF( MOD(ls+l+li,2)/=0 .OR. MOD(le+l+lo,2)/=0 ) CYCLE

          IF(l<=5) nmax = nmin +(5-l)*(nx-nmin)/5

          Ti = 0._RP
          integ0 = 0._RP
          Si = 0._RP

          is = 1
          xil = x(is)**l
          xil1 = xil*x(is)
          IF( xil1>VSmall ) THEN
            Ti = Ti +w(is)*chi_a(is,ls)*chi_0(is,li)*xil
            integ0 = integ0 +w(is)*Ti*wf(is)*chi_b(is,le)/xil1
          END IF

          DO i = is+1,nmax
            Si = Si +w(i-1)*chi_b(i-1,le)*wf(i-1)*xil
            xil = x(i)**l
            IF( xil>=VBig ) EXIT
            Ti = Ti +w(i)*chi_a(i,ls)*chi_0(i,li)*xil
            xil1 = xil*x(i)
            IF( xil1<VSmall ) CYCLE
            integ0 = integ0 +( Ti*w(i)*wf(i)*chi_b(i,le) &
              +Si*w(i)*chi_a(i,ls)*chi_0(i,li) )/xil1
          END DO

          integ(l, li ) = integ(l, li ) +integ0*symbol_3j(ls, l, li, 0, 0, 0 ) &
            *symbol_3j(lo, l, le, 0, 0, 0 )*(2*li+1)*zi**li*exp( -zi*sig_0(li) )
        END DO
      END DO

      DO mo = 0,lo; DO me = -le,le
        IF( ABS(mo+me)>ls ) CYCLE
        DO li = 0,limax
          DO l = MAX(ABS(ls-li), ABS(le-lo)), MIN(ls+li, le+lo)
            IF( MOD(ls+l+li,2)/=0 .OR. MOD(le+l+lo,2)/=0 ) CYCLE
            integral(ls,le,me,mo) = integral(ls,le,me,mo) +integ(l, li ) &
              *symbol_3j(ls, l, li, mo+me, -me-mo, 0 ) &
              *symbol_3j(lo, l, le, mo, -me-mo, me)
          END DO
        END DO
      END DO; END DO

      integral(ls,le,:,:) = integral(ls,le,:,:)*zi**(-ls-le)*SQRT((2*lo+1)*(2*ls+1) &
        *(2*le+1._RP) )*EXP( zi*( sig_b(le) +sig_a(ls) ) )

    END DO; END DO

  END SUBROUTINE dwb_integrals

  PURE COMPLEX(RP) FUNCTION tpw( n, l, m, e, ke, k)
    USE constants ,ONLY: pi
    USE trigo ,ONLY: cartez2spher
    USE utils ,ONLY: norm_fac
    USE special_functions ,ONLY: spherical_harmonic, fac
    INTEGER, INTENT(IN) :: n,l,m
    REAL(RP), INTENT(IN) :: e,ke(3)
    REAL(RP), INTENT(IN),OPTIONAL :: k(3)

    REAL(RP) :: q(3)
    INTEGER :: j
    REAL(RP) :: kem,thetae,phie,A
    COMPLEX(RP) :: kec

    IF(PRESENT(k)) THEN
      q = k-ke
    ELSE
      q = -ke
    ENDIF

    CALL cartez2spher( q, kem, thetae, phie)
    a = kem**2 +e**2
    kec = CMPLX( 0._RP, kem, KIND=RP )

    tpw = 0._RP
    DO j = 0, (n-l)/2
      tpw = tpw +(-0.25_RP*a/e**2)**j*fac(n-j)/( fac(j)*fac(n-l-2*j) )
    END DO

    tpw = tpw*SQRT(2./pi)*norm_fac(e,n)*fac(n-l)*(kec/e)**l*(2.*e/a)**n /a &
      *spherical_harmonic(l, m, thetae, phie )

  END FUNCTION tpw

  PURE COMPLEX(RP) FUNCTION tcw( n, l, m, e, alpha, ke, k)
    USE constants ,ONLY: pi
    USE special_functions ,ONLY: fac
    USE utils ,ONLY: norm_fac
    USE trigo ,ONLY: nrm2

    INTEGER, INTENT(IN) :: n, l, m
    REAL(RP), INTENT(IN) :: e, alpha, ke(3), k(3)

    REAL(RP)    :: kem, km, a, aj1, ke_t(3)
    COMPLEX(RP) :: w, kec, ekec, gam(0:n), f21(0:n,0:n), alphac, w1m, kep, kp
    COMPLEX(RP) :: tmp_j, tmp_j1, tmp_m1, cst_j, cst_j1, cst_m1
    REAL(RP)    :: tmp_s, tmp_s3, tmp_s1, cst_s, cst_s3, cst_s1, cst_s2
    COMPLEX(RP), PARAMETER :: zi = (0._RP, 1._RP)
    INTEGER :: j, j1, ma, m1, s, s1, s2, s3

    ma = ABS(m)
    kem = nrm2(ke)
    ke_t = -ke
    km = nrm2(k)
    IF( m>=0 ) THEN
      kp = CMPLX( k(1), k(2), KIND=RP )
      kep = CMPLX( ke_t(1), ke_t(2), KIND=RP )
    ELSE
      kp = CMPLX( k(1), -k(2), KIND=RP )
      kep = CMPLX( ke_t(1), -ke_t(2), KIND=RP )
    END IF

    alphac = CMPLX( 0._RP, alpha, KIND=RP) ! i*\alpha
    kec    = CMPLX( 0._RP, kem  , KIND=RP) ! i*ke
    ekec   = CMPLX( e    , -kem , KIND=RP) ! (\epsilon-ike)
    a      = nrm2(k+ke_t)**2 +e**2
    w      = CMPLX( nrm2(k+ke_t)**2 -km**2 +kem**2   , 2.*e*kem , KIND=RP) /a ! w = b/a
    w1m    = CMPLX( e**2 +km**2 -kem**2, -2.*e*kem, KIND=RP) /a ! (1-w)

    gam(0) = 1._RP
    DO j1 = 1, n
      gam(j1) = gam(j1-1)*CMPLX( j1, -alpha, KIND=RP)
    END DO

    DO j1 =  0,n
      aj1 = j1 +1._RP
      f21(j1, 0) = 1._RP
      f21(j1, 1) = 1. +alphac*w / ( aj1*w1m )
      IF( n<2 .OR. j1>n-2 ) CYCLE
      DO j = 1, n -j1 -1
        f21( j1,  j +1 ) = (1. +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      END DO
    END DO

    tcw = 0._RP
    cst_j = -0.25*a/ekec**2
    cst_j1 = kec/ekec
    cst_m1 = kep/kp
    cst_s = -(km/k(3))**2
    cst_s3 = ke_t(3)/k(3)
    cst_s1 = kem/km
    cst_s2 = 2.*DOT_PRODUCT(k,ke_t)
    DO j = 0,(n-l)/2
      tmp_j = fac(n-j) /fac(j)*cst_j**j
      !tmp_j = fac(n-j)/fac(j)*(-0.25*a/ekec**2)**j
      DO j1 = 0,(n-l-2*j)
        !tmp_j1 = (kec/ekec)**j1 /( fac(j1)*fac(n-l-2*j-j1) )*tmp_j
        tmp_j1 = cst_j1**j1 /( fac(j1)*fac(n-l-2*j-j1) )*tmp_j
        DO m1 = 0,ma
          tmp_m1 = cst_m1**m1 /( fac(ma-m1)*fac(m1) )*tmp_j1
          !tmp_m1 = (kep/kp)**m1 /( fac(ma-m1)*fac(m1) )*tmp_j1
          DO s = 0,(l-ma)/2
            tmp_s = cst_s**s /( fac(l-s) )
            !tmp_s = (-1)**s*(km/k(3))**(2*s) /( fac(l-s) )
            DO s3 = 0,l-ma-2*s
              tmp_s3 = cst_s3**s3 /( fac(s3)*fac(l-ma-2*s-s3) )*tmp_s
              !tmp_s3 = (ke_t(3)/k(3))**s3 /( fac(s3)*fac(l-ma-2*s-s3) )*tmp_s
              DO s1 = 0,s
                tmp_s1 = cst_s1**s1 /fac(s-s1)*tmp_s3
                !tmp_s1 = (kem/km)**s1 /fac(s-s1)*tmp_s3
                DO s2 = 0,s1
                  tcw = tcw +gam(m1+s3+2*s1-s2+j1)*f21(m1+s3+2*s1-s2+j1, n -j &
                    -(m1+s3+2*s1-s2+j1) )/( fac(s2)*fac(s1-s2)*fac(m1+s3+2*s1-s2+j1) ) &
                    *tmp_m1*tmp_s1*cst_s2**s2
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    tcw = tcw*norm_fac(e, n )*SQRT(l+0.5_RP)*fac(n-l)*SQRT( fac(l-ma)/fac(l+ma) ) &
      *(-1)**ma*fac(ma)*(zi*k(3)/ekec)**l*(2.*ekec/a)**n/a*(kp/k(3))**ma &
      *w1m**(-alphac) /pi
      !*powcc(w1m,-alpha)/pi
    IF( MOD(m,2)<0 ) tcw = -tcw

  END FUNCTION tcw

  PURE COMPLEX(RP) FUNCTION tcw0( n, l, m, e, alpha, ke)
    USE constants ,ONLY: pi
    USE trigo ,ONLY: cartez2spher
    USE utils ,ONLY: norm_fac
    USE special_functions ,ONLY: spherical_harmonic, fac
    INTEGER, INTENT(IN) :: n, l, m
    REAL(RP), INTENT(IN) :: e, alpha, ke(3)

    REAL(RP) :: kem, thetae, phie, a, aj1
    COMPLEX(RP) :: w, kec, ekec, gam(0:n), f21(0:n,0:(n-l)), alphac, w1m, tmp
    INTEGER :: j, j1

    CALL cartez2spher( -ke, kem, thetae, phie)

    alphac = CMPLX( 0._RP, alpha, KIND=RP) ! i*\alpha
    kec    = CMPLX( 0._RP, kem  , KIND=RP) ! i*ke
    ekec   = CMPLX( e    , -kem , KIND=RP) ! (\epsilon-ike)
    a      = kem**2 +e**2
    w      = CMPLX( 2.*kem**2   , 2.*e*kem , KIND=RP) /a ! w = b/a
    w1m    = CMPLX( e**2 -kem**2, -2.*e*kem, KIND=RP) /a ! (1-w)

    gam(0) = 1._RP
    DO j1 = 1, n
      gam(j1) = gam(j1-1)*CMPLX( j1, -alpha, KIND=RP)
    END DO

    DO j1 =  0,n-l
      aj1 = l +j1 +1._RP
      f21(j1, 0) = 1._RP
      f21(j1, 1) = 1._RP +alphac*w / ( aj1*w1m )
      IF(n-l<2) CYCLE
      DO j = 1, n -l -1
        f21( j1,  j +1 ) = (1._RP +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      END DO
    END DO

    tcw0 = 0._RP
    DO j = 0,(n-l)/2
      tmp = fac(n-j) /fac(j)*(-0.25*a/ekec**2)**j
      DO j1 = 0,(n-l-2*j) !, j1+2*j <= n-l )
        tcw0 = tcw0 +1./( fac(j1)*fac(n-l-2*j-j1)*fac(l+j1) )*(kec/ekec)**j1*gam(l+j1) &
          *f21(j1, n -l -j -j1 )*tmp !*fac(n-j)/ fac(j)*(-0.25*a/ekec**2)**j
      END DO
    END DO

    tcw0 = tcw0*norm_fac(e, n )*SQRT(2./pi)*fac(n-l)*(kec/ekec)**l*(2.*ekec/a)**n/a &
      *spherical_harmonic(l, m, thetae, phie)*w1m**(-alphac)!*powcc(w1m,-alpha)

  END FUNCTION tcw0

  ELEMENTAL COMPLEX(RP) FUNCTION powcc(z1, y2)
    COMPLEX(RP), INTENT(IN) :: z1
    REAL(RP), INTENT(IN) :: y2
    REAL(RP) :: theta,zm

    theta = ATAN2( AIMAG(z1), REAL(z1, KIND=RP) )
    zm = LOG( ABS(z1) )
    powcc = EXP(-y2*theta)*CMPLX( COS(y2*zm), SIN(y2*zm), KIND=RP )

  END FUNCTION powcc

END MODULE fdcs_e2e
