submodule(fdcs_e2e) fdcs_e2e
  implicit none
contains

  module subroutine fdcs_fba_pw(in_unit,out_unit)
    use constants ,only: pi, ev, deg
    use special_functions ,only: factorial
    use input ,only: read_input, read_orbit
    use trigo ,only: spher2cartez

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom, Orbit
    real(RP) :: Ei, Es, Ee
    real(RP) :: thetas
    integer :: step(3), exchange
    integer :: nelec
    integer :: lo, no
    real(RP), allocatable :: a(:), e(:)
    integer, allocatable :: n(:)

    real(RP), parameter   :: phie = 0._RP
    real(RP) :: kim, ksm, kem, km
    real(RP) :: ki(3), ks(3), ke(3), k(3)

    real(RP) :: factor, sigma
    complex(RP) :: D_term
    integer :: i, io, mo

    real(RP) :: k2(3), k2m = 0._RP
    complex(RP) :: E_term = (0._RP, 0._RP)

    call factorial()

    call read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange)

    call read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    call spher2cartez( kim, 0._RP, 0._RP, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = nelec*4._RP*ksm*kem/kim

    write( out_unit, * ) "Theta TDCS_PW"
    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      if(exchange==1) then
        k2 = ki -ke
        k2m = norm2(k2)
      end if

      sigma = 0._RP
      do mo = 0, lo

        D_term = (0._RP,0._RP)
        do io = 1, no
          D_term = D_term +a(io)*( tpw(n(io), lo, mo, e(io), ke, k ) &
            -tpw(n(io), lo, mo, e(io), ke ) )
        end do

        if(exchange==1) then
          E_term = (0._RP,0._RP)
          do io = 1, no
            E_term = E_term &
              +a(io)*( tpw(n(io), lo, mo, e(io), ks, k2 ) &
              -tpw(n(io), lo, mo, e(io), ks ) )
          end do
          sigma = sigma +(1+mo)*( abs(D_term/km**2)**2 +abs(E_term/k2m**2)**2 &
            -real( D_term*conjg(E_term)/(km**2*k2m**2 ), RP ) )
        else
          sigma = sigma +(1+mo)*abs(D_term/km**2)**2
        end if

      end do

      sigma = factor*sigma

      if( show_output ) print '(1x,i4,1x,es15.8)',i,sigma
      write( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_fba_pw


  module subroutine fdcs_fba_cw(in_unit,out_unit)
    use constants ,only: ev, deg, pi
    use trigo ,only: spher2cartez
    use special_functions ,only: factorial
    use input ,only: read_input, read_orbit
    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom, Orbit
    real(RP) :: Ei, Es, Ee
    real(RP) :: thetas
    integer :: step(3), exchange
    integer :: nelec, ze, zs
    integer :: lo, no
    real(RP), allocatable :: a(:), e(:)
    integer, allocatable :: n(:)

    real(RP), parameter   :: phie = 0._RP
    real(RP) :: kim, ksm, kem, km, alpha
    real(RP) :: ki(3), ks(3), ke(3), k(3)

    real(RP) :: factor, sigma
    complex(RP) :: term
    integer :: i, io, mo

    call factorial()

    call read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange)

    call read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    zs = -1 ! Charge of the projectile
    kem = sqrt(2.*Ee*eV)
    ze = -1 ! Charge of ejected particle
    alpha = -ze/kem

    call spher2cartez( kim, 0._RP, 0._RP, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = nelec*4._RP*ksm*kem / (kim*km**4)
    !\abs{ \exp{\pi\alpha/2}*\Gamma(1-i\alpha) }^2
    if(ze/=0) factor = factor*2.*pi*alpha/(1._RP-exp(-2.*pi*alpha))


    write( out_unit, * ) "Theta TDCS_CW"
    do i = step(1),step(2),step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      do mo = 0, lo
        term = (0._RP,0._RP)
        do io = 1, no
          term = term +a(io)*( tcw(n(io), lo, mo, e(io), alpha, ke, k ) &
            +zs*tcw0(n(io), lo, mo, e(io), alpha, ke ) )
        end do
        sigma = sigma +(mo+1)*abs(term)**2
      end do

      sigma = factor*sigma/(2*lo+1)

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_fba_cw

  module subroutine fdcs_fba_dw(in_unit,out_unit)
    use ieee_arithmetic ,only: ieee_is_nan
    use constants ,only: pi, deg, ev
    use special_functions ,only: cgamma, spherical_harmonic, ricbes, factorial !&
      !, coul90, symbol_3j
    use utils ,only: norm_fac, y1y2y3, calculate_U, ode_second_dw
    use input ,only: read_input, read_orbit
    use trigo ,only: spher2cartez, cartez2spher

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom, Orbit
    real(RP) :: Ei, Es, Ee
    real(RP) :: thetas
    integer :: step(3), exchange
    integer :: nelec
    integer :: lo, no
    real(RP), allocatable :: a(:), e(:)
    integer, allocatable :: n(:)

    real(RP), parameter   :: phie = 0._RP
    real(RP) :: kim, ksm, kem, km
    real(RP) :: theta, phi
    real(RP) :: ki(3), ks(3), ke(3), k(3)

    integer, parameter :: lmax = 47, lemax = 15, np = 2400
    real(RP) :: x(0:np), wf(0:np)
    real(RP), allocatable :: chi_b(:,:), chi_0a(:,:)
    real(RP) :: sigma_le(0:lemax), delta(0:lemax), integral(0:lemax,0:lmax)
    real(RP) :: etae, rc, h
    complex(RP) :: tmp_z
    complex(RP), parameter :: zi = (0._RP, 1._RP)
    integer :: le, l, me, ze

    real(RP) :: factor, sigma
    complex(RP) :: term, term0
    integer :: i, io, mo

    real(RP) :: gc_0a(0:lmax), fdc_0a(0:lmax), gdc_0a(0:lmax), fc_0a(0:lmax)
    real(RP), allocatable :: U(:,:), U_tmp(:)

    call factorial()

    call read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange )

    call read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)
    ze = 1

    call spher2cartez( kim, 0._RP, 0._RP, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    call cartez2spher( k, km, theta, phi)

    factor = nelec*4._RP*ksm*kem / (kim*km**4)
    factor = factor*2./(pi*kem**2)

    rc = 10._RP
    h = rc/np

    etae = -ze/kem

    x = [(i*h, i = 0, np) ]
    allocate( chi_b(0:np,0:lemax), chi_0a(0:np,0:lmax) )
    chi_b(0,:) = 0._RP
    chi_0a(0,:) = 0._RP
    do i = 1,np
      call ricbes(km*x(i), lmax, fc_0a, gc_0a, fdc_0a, gdc_0a, le)
      chi_0a(i,:) = fc_0a/(km*x(i))
    end do

    where(ieee_is_nan(chi_0a) )
      chi_0a = 0._RP
    end where

    allocate( U(0:np,0:lemax), U_tmp(0:np) )
    call calculate_U(Atom, Orbit, x, U_tmp, 1 )

    if(ze/=0) then
      U(0,:) = -huge(1._RP)
    else
      U(0,:) = -kem**2
    end if

    U(1:np,lemax) = -kem**2 -2.*( ze +U_tmp(1:np) ) /x(1:np)
    do l = 0,lemax
      U(1:np,l) = l*(l+1)/x(1:np)**2 +U(1:np,lemax)
    end do

    call ode_second_dw(kem, lemax, rc, ze, U, chi_b, delta)

    deallocate(U,U_tmp)

    wf = 0._RP
    do io = 1, no
      wf = wf +a(io)*norm_fac(e(io), n(io) )*x**n(io)*exp(-e(io)*x)
    end do

    do le = 0,lemax
      tmp_z = cgamma( cmplx(le+1._RP, etae, kind=RP), 1 )
      tmp_z = exp( zi*aimag(tmp_z) )
      sigma_le(le) = atan2( aimag(tmp_z), real(tmp_z,kind=RP) )
    end do

    integral = 0._RP
    do l = 0,lmax
      do le = 0,lemax
        if( mod(le+l+lo,2)/=0 ) cycle
        integral(le, l) = 0.5*h*sum( wf(1:np)*chi_b(1:np,le)*chi_0a(1:np,l) &
          +wf(0:np-1)*chi_b(0:np-1,le)*chi_0a(0:np-1,l) )
      end do
    end do

    term0 = zi**(-lo)*exp( cmplx(0._RP, -(sigma_le(lo)+delta(lo)), kind=RP ) ) &
      *0.5*h*sum( wf(1:np)*chi_b(1:np,lo) +wf(0:np-1)*chi_b(0:np-1,lo) )


    write( out_unit, * ) "Theta TDCS_DW"
    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      do mo = 0, lo

        term = (0._RP, 0._RP)
        do l = 0,lmax
          do le = 0,lemax
            if( mod(le+l+lo,2)/=0 ) cycle
            tmp_z = (0._RP, 0._RP)
            do me = -le,le
              if( abs(me-mo)>l ) cycle
              tmp_z = tmp_z +y1y2y3(le, l, lo, -me, me-mo, mo) &
                *spherical_harmonic(le, me, i*deg, phie) &
                *spherical_harmonic(l, mo-me, theta, phi)
            end do
            term = term +zi**(l-le)*integral(le, l)*tmp_z &
              *exp( cmplx(0._RP, -(sigma_le(le)+delta(le)), kind=RP ) )

          end do
        end do
        term = 4.*pi*term

        sigma = sigma +(mo+1)*abs( term &
          -spherical_harmonic(lo, mo, i*deg, phie)*term0 )**2

      end do

      sigma = factor*sigma/(2*lo+1)

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_fba_dw

  module subroutine fdcs_dwb(in_unit,out_unit)
    use constants ,only: ev, deg, pi
    use special_functions ,only: spherical_harmonic, cgamma, factorial
    use utils ,only: norm_fac, calculate_U
    use input ,only: read_input, read_orbit
    use trigo ,only: spher2cartez
    use integration ,only: gauleg

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom, Orbit
    real(RP) :: Ei, Es, Ee
    real(RP) :: thetas
    integer :: step(3), exchange, PCI = 2
    integer :: nelec
    integer :: lo, no
    real(RP), allocatable :: a(:), e(:)
    integer, allocatable :: n(:)

    real(RP), parameter   :: phie = 0._RP
    real(RP) :: kim, ksm, kem

    integer, parameter :: nx = 31*64, nr = 48000

    integer, parameter :: limax = 95
    real(RP), allocatable, target :: chi_0(:,:)
    real(RP), target :: delta_li(0:limax)

    integer, parameter :: lsmax = 95
    real(RP), allocatable, target :: chi_a(:,:)
    real(RP) :: sigma_ls(0:lsmax), etas
    real(RP), target :: delta_ls(0:lsmax)
    complex(RP) :: ylms(0:lsmax,-lsmax:lsmax)
    integer :: ls, zs = 1

    integer, parameter :: lemax = 23
    real(RP), allocatable, target :: chi_b(:,:)
    real(RP) :: sigma_le(0:lemax), etae
    real(RP), target :: delta_le(0:lemax)
    complex(RP), allocatable :: ylme(:,:,:)
    integer :: le, me, ze = 1

    real(RP), allocatable :: r(:), wf(:)
    real(RP) :: rc, h
    complex(RP), allocatable :: integral(:,:,:,:), integralx(:,:,:,:)
    complex(RP) :: tmp_z

    real(RP) :: factor, sigma
    complex(RP) :: termd, termx
    integer :: i, io, mo = 0

    real(RP), allocatable :: x(:), w(:)

    real(RP), allocatable :: U_tmp(:)
    real(RP), allocatable :: xt(:), wt(:)
    integer, parameter :: nc = nx/64

    call factorial()

    call read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit, exchange )

    call read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    etas = -zs/ksm
    etae = -ze/kem

    factor = nelec* (2.*pi)**4*ksm*kem / kim
    factor = factor*2._RP/pi**4 /(kem*ksm*kim)**2

    rc = 250._RP

    allocate(x(nx), w(nx) )

    do i = 1,nc
      call gauleg( (i-1)*rc/nc, i*rc/nc, xt, wt, 64 )
      x( (i-1)*64+1: i*64 ) = xt
      w( (i-1)*64+1: i*64 ) = wt
    end do
    deallocate( xt, wt )

    allocate(chi_0(nx,0:limax), chi_a(nx,0:lsmax), chi_b(nx,0:lemax) )

    h = rc/nr
    allocate( r(0:nr) )
    r = [(i*h,i = 0,nr)]

    allocate( U_tmp(0:nr) )

    call calculate_U(Atom, Orbit, r, U_tmp, 0)
    call calculate_chi( kim, r, U_tmp, 0, x, chi_0, delta_li )

    call calculate_U(Atom, Orbit, r, U_tmp, 1)
    call calculate_chi( ksm, r, U_tmp, zs, x, chi_a, delta_ls )
    call calculate_chi( kem, r, U_tmp, ze, x, chi_b, delta_le )

    deallocate( U_tmp )

    allocate(wf(nx))
    wf = 0._RP
    do io = 1, no
      wf = wf +a(io)*norm_fac(e(io), n(io) )*x**n(io)*exp(-e(io)*x)
    end do

    chi_b(:,lo) = chi_b(:,lo) -wf*sum( w*chi_b(:,lo)*wf )
    chi_a(:,lo) = chi_a(:,lo) -wf*sum( w*chi_a(:,lo)*wf )

    tmp_z = cgamma( cmplx( 1._RP, etae, RP), 0 )
    sigma_le(0) = atan2( aimag(tmp_z), real(tmp_z, RP) )
    do le = 0,lemax-1
      sigma_le(le+1) = sigma_le(le) +atan2( etae, le+1._RP )
    end do
    sigma_le = sigma_le +delta_le

    tmp_z = cgamma( cmplx( 1._RP, etas, RP), 0 )
    sigma_ls(0) = atan2( aimag(tmp_z), real(tmp_z, RP) )
    do ls = 0,lsmax-1
      sigma_ls(ls+1) = sigma_ls(ls) +atan2( etas, ls+1._RP )
    end do
    sigma_ls = sigma_ls +delta_ls

    call dwb_integrals(chi_0, chi_a, chi_b, delta_li, sigma_ls, sigma_le &
      , wf, x, w, lo, integral)
    if(exchange==1) then
      call dwb_integrals(chi_0, chi_b, chi_a, delta_li, sigma_le, sigma_ls &
        , wf, x, w, lo, integralx)
    end if

    ylms = (0._RP, 0._RP)
    do ls = 0,lsmax
      ylms(ls,0) = spherical_harmonic(ls,0,thetas*deg,pi)
      do me = 1,ls
        ylms(ls,me) = spherical_harmonic(ls,me,thetas*deg,pi)
      end do
      ylms(ls,-2:-ls:-2) = ylms(ls,2:ls:2)
      ylms(ls,-1:-ls:-2) = -ylms(ls,1:ls:2)
    end do

    allocate(ylme(0:lemax,-lemax:lemax,step(1)/step(3):step(2)/step(3)))
    ylme = (0._RP, 0._RP)
    do i = lbound(ylme,3),ubound(ylme,3)
      do le = 0,lemax
        ylme(le,0,i) = spherical_harmonic(le,0,i*step(3)*deg,phie)
        do me = 1,le
          ylme(le,me,i) = spherical_harmonic(le,me,i*step(3)*deg,phie)
        end do
        ylme(le,-2:-le:-2,i) = ylme(le,2:le:2,i)
        ylme(le,-1:-le:-2,i) = -ylme(le,1:le:2,i)
      end do
    end do
!!$OMP PARALLEL DO PRIVATE(sigma, mo, termd, termx, le, ls, me)

    write( out_unit, * ) "Theta TDCS_DWBA"
    do i = step(1), step(2), step(3)

      sigma = 0._RP
      do mo = 0, lo

        termd = (0._RP, 0._RP)
        do le = 0,lemax
          do ls = 0,lsmax
            do me = -le,le
              if( abs(me+mo)>ls ) cycle
              termd = termd +ylme(le,me,i/step(3))*ylms(ls,mo+me) &
                *integral(ls, le, me, mo )
            end do
          end do
        end do

        termx = (0._RP, 0._RP)
        if(exchange==1) then
          do ls = 0,lemax
            do le = 0,lsmax
              do me = -le,le
                if( abs(me+mo)>ls ) cycle
                termx = termx +ylms(le,me)*ylme(ls,mo+me,i/step(3)) &
                  *integralx(ls, le, me, mo )
              end do
            end do
          end do
        end if

        sigma = sigma +(mo+1)*( abs(termd)**2 +abs(termx)**2 &
          -real( conjg(termd)*termx, RP ) )
      end do

      sigma = factor*sigma/(2*lo+1)
      if(PCI>=1) then
        call PCI_EFFECTS(i,sigma)
      end if

      if( show_output ) print'(1x,i4,1x,es15.8)', i, sigma
      write(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  contains

    elemental subroutine PCI_EFFECTS(i,sigma)
      use special_functions ,only: conhyp_opt
      integer, intent(in) :: i
      real(RP), intent(inout) :: sigma
      real(RP) :: kesm
      real(RP) :: r12_av, Et
      kesm = sqrt( kem**2 +ksm**2 -2*kem*ksm*cos( (thetas+i)*deg ) )/2._RP
      sigma = pi/(kesm*(exp(pi/kesm) -1._RP ) )*sigma
      if(PCI==2) then
        Et = (Es +Ee)*eV
        r12_av = pi**2/(16.*Et)*( 1._RP +(0.627/pi)*sqrt(Et)*log(Et) )**2
        sigma = abs( conhyp_opt( 0.5/kesm, -2.*kesm*r12_av  ) )**2*sigma
      end if
    end subroutine

  end subroutine fdcs_dwb

  module subroutine fdcs_bbk(in_unit,out_unit)
    use constants ,only: pi, ev, deg
    use special_functions ,only: factorial
    use input ,only: read_input, read_orbit
    use trigo ,only: spher2cartez

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom, Orbit
    real(RP) :: Ei, Es, Ee
    real(RP) :: thetas
    integer :: step(3)
    integer :: nelec
    integer :: lo, no
    real(RP), allocatable :: a(:), e(:)
    integer, allocatable :: n(:)

    real(RP), parameter   :: phie = 0._RP
    real(RP) :: kim, ksm, kem, km
    real(RP) :: ki(3), ks(3), ke(3), k(3)

    real(RP) :: factor, sigma
    complex(RP) :: D_term
    integer :: i, io, mo

    call factorial()

    call read_input(in_unit,Ei, Es, Ee, thetas, step, Atom, orbit)

    call read_orbit(Atom//'_'//Orbit, nelec, lo, no, n, a, e )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    call spher2cartez( kim, 0._RP, 0._RP, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = nelec*4._RP*ksm*kem/(kim)

    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._RP
      do mo = 0, lo

        D_term = (0._RP,0._RP)
        do io = 1, no
          D_term = D_term +a(io)*( tpw(n(io), lo, mo, e(io), ke, k ) &
            -tpw(n(io), lo, mo, e(io), ke ) )
        end do
        sigma = sigma +(1+mo)*abs(D_term/km**2)**2
      end do

      sigma = factor*sigma

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_bbk

  complex(RP) function U_bbk(alpha1, alpha2, alpha3, k1, k2, k3, lam1, lam2, lam3 &
    , p1, p2)
    use constants ,only: pi
    use integration ,only: gauleg

    real(RP), intent(in) :: alpha1, alpha2, alpha3
    real(RP), intent(in) :: k1(3), k2(3), k3(3)
    real(RP), intent(in) :: lam1, lam2, lam3
    real(RP), intent(in) :: p1(3), p2(3)

    complex(RP), parameter :: zi = (0._RP, 1._RP)

    real(RP),allocatable :: t3(:), y(:), wy(:), s(:), ws(:)
    real(RP) :: p11(3), p22(3), q1(3), q2(3)
    real(RP) :: q1m, q2m
    real(RP) :: k1m, k2m, k3m
    complex(RP) :: lam33
    complex(RP) :: sig0, sig1, sig2, sig12
    complex(RP) :: sig0_a(3), sig1_a(3), sig2_a(3), sig12_a(3)

    integer :: i, j
    complex(RP) :: N_t3, Integ

    U_bbk = ( 0._RP, 0._RP)

    k1m = norm2(k1)
    k2m = norm2(k2)
    k3m = norm2(k3)

    call gauleg(-10._RP, 10._RP, y, wy, 64)
    allocate( t3(size(y)) )
    t3 = 1._RP /(1._RP+ exp(y) )

    call gauleg(0._RP, 10._RP, s, ws, 64)

    Integ = (0._RP, 0._RP)
    do i = 1,64
      p11 = p1 -t3(i)*k3
      q1 = k1 +p11
      q1m = norm2(q1)
      p22 = p2 -t3(i)*k3
      q2 = k2 +p22
      q2m = norm2(q2)
      lam33 = cmplx(lam3, -t3(i)*k3m, RP )

      sig0_a = [ cmplx( (lam1+lam2)**2+norm2(q1-q2)**2, 0._RP, RP) &
        , 2*(lam2*(lam1**2+lam33**2+q1m**2) +lam1*(lam2**2+lam33**2+q2m**2) ) &
        , ( (lam1+lam33)**2+q1m**2)*((lam2+lam33)**2+q2m**2) ]
      sig1_a = -2*[ dot_product(q1-q2,k1) +zi*(lam1+lam2)*k1m &
        , 2*lam2*dot_product(q1,k1) +zi*( (lam2**2+lam33**2+q2m**2)*k1m+2*lam1*lam2*k1m ) &
        , ((lam2+lam33)**2 +q2m**2)*( dot_product(q1,k1) +zi*(lam1+lam33)*k1m ) ]
      sig2_a = -2*[ dot_product(q2-q1,k2) +zi*(lam1+lam2)*k2m &
        , 2*lam1*dot_product(q2,k2) +zi*( (lam1**2+lam33**2+q1m**2)*k2m+2*lam1*lam2*k2m ) &
        , ((lam1+lam33)**2 +q1m**2)*( dot_product(q2,k2) +zi*(lam2+lam33)*k2m ) ]
      sig12_a = -2*[ cmplx( k1m*k2m +dot_product(k1,k2), 0._RP, RP) &
        , 2*( (lam1*k1m-zi*dot_product(q1,k1))*k2m +(lam2*k2m-zi*dot_product(q2,k2))*k1m) &
        , -2*( dot_product(q1,k1)+zi*(lam1+lam33)*k1m) &
        *( dot_product(q2,k2)+zi*(lam2+lam33)*k2m) ]

      N_t3 = (0._RP, 0._RP)
      do j = 1,64
        sig0 = sig0_a(1)*(s(j)+2*lam33)*s(j) +sig0_a(2)*s(j) +sig0_a(3)
        sig1 = sig1_a(1)*(s(j)+2*lam33)*s(j) +sig1_a(2)*s(j) +sig1_a(3)
        sig2 = sig2_a(1)*(s(j)+2*lam33)*s(j) +sig2_a(2)*s(j) +sig2_a(3)
        sig12 = sig12_a(1)*(s(j)+2*lam33)*s(j) +sig12_a(2)*s(j) +sig12_a(3)
        N_t3 = (sig0/(sig0+sig1) )**(-zi*alpha1)*(sig0/(sig0+sig2) )**(-zi*alpha2)/sig0 &
          +N_t3
      end do
      Integ = Integ +sinh(pi*alpha3)/(2*pi*zi)*(4*pi)**2*N_t3
    end do

  end function

  subroutine calculate_chi( km, r, U_tmp, z, x, chi, delta )
    use utils ,only: ode_second_dw, intrpl
    real(RP), intent(in) :: U_tmp(0:), r(0:), x(:), km
    integer, intent(in) :: z
    real(RP), intent(out) :: chi(:,0:), delta(0:)

    real(RP), allocatable :: U(:,:), chi_tmp(:,:)
    real(RP) :: rc
    integer :: l, lmax, nr

    lmax = size(delta) -1
    nr = size(r) -1
    rc = r(nr)

    allocate( U(0:nr,0:lmax), chi_tmp(0:nr,0:lmax) )

    if(z==0) then
      U(0,0:lmax) = -km**2
    else
      U(0,0:lmax) = -huge(1._RP)
    end if

    U(1:nr,lmax) = -km**2 -2.*( z*1._RP +U_tmp(1:nr) )/r(1:nr)
    do l = 0,lmax
      U(1:nr,l) = l*(l+1)/r(1:nr)**2 +U(1:nr,lmax)
    end do

    call ode_second_dw(km, lmax, rc, z, U, chi_tmp, delta )

    do l = 0,lmax
      call INTRPL(r, chi_tmp(:,l), x, chi(:,l))
    end do

  end subroutine calculate_chi

  subroutine dwb_integrals( chi_0, chi_a, chi_b, sig_0, sig_a, sig_b, wf, x, w, lo &
    , integral )
    use special_functions ,only: symbol_3j
    real(RP), intent(in) :: chi_0(:,0:), chi_a(:,0:), chi_b(:,0:), wf(:)
    real(RP), intent(in) :: x(:), w(:)
    real(RP), intent(in) :: sig_0(0:), sig_a(0:), sig_b(0:)
    integer, intent(in) :: lo
    complex(RP), allocatable, intent(out) :: integral(:,:,:,:)

    real(RP), parameter :: VSmall = tiny(1._RP), VBig = huge(1._RP)
    complex(RP), parameter :: zi = (0._RP, 1._RP)
    integer :: limax, lsmax, lemax, nx, nmin
    complex(RP), allocatable :: integ(:,:)
    integer :: li, ls, le, l, me, mo
    ! Block variables
    real(RP) :: Ti, Si
    real(RP) :: integ0
    real(RP) :: xil, xil1
    integer :: i, is, nmax

    nx = size(x)
    limax = size(sig_0) -1
    lsmax = size(sig_a) -1
    lemax = size(sig_b) -1

    nmin = count(x<=25._RP)

    allocate( integral(0:lsmax, 0:lemax, -lemax:lemax, 0:lo) )
    allocate( integ( 0:min(limax+lsmax,lemax+lo),  0:limax ) )
    integral = 0._RP
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(integ,li,l,me,mo ,ti,si,integ0,xil,xil1,i,is,nmax)
    do ls = 0,lsmax; do le = 0,lemax
      integ = 0._RP
      do li = 0,limax
        do l = max(abs(ls-li), abs(le-lo)), min(ls+li, le+lo)
          if( mod(ls+l+li,2)/=0 .OR. mod(le+l+lo,2)/=0 ) cycle

          if(l<=5) nmax = nmin +(5-l)*(nx-nmin)/5

          Ti = 0._RP
          integ0 = 0._RP
          Si = 0._RP

          is = 1
          xil = x(is)**l
          xil1 = xil*x(is)
          if( xil1>VSmall ) then
            Ti = Ti +w(is)*chi_a(is,ls)*chi_0(is,li)*xil
            integ0 = integ0 +w(is)*Ti*wf(is)*chi_b(is,le)/xil1
          end if

          do i = is+1,nmax
            Si = Si +w(i-1)*chi_b(i-1,le)*wf(i-1)*xil
            xil = x(i)**l
            if( xil>=VBig ) exit
            Ti = Ti +w(i)*chi_a(i,ls)*chi_0(i,li)*xil
            xil1 = xil*x(i)
            if( xil1<VSmall ) cycle
            integ0 = integ0 +( Ti*w(i)*wf(i)*chi_b(i,le) &
              +Si*w(i)*chi_a(i,ls)*chi_0(i,li) )/xil1
          end do

          integ(l, li ) = integ(l, li ) +integ0*symbol_3j(ls, l, li, 0, 0, 0 ) &
            *symbol_3j(lo, l, le, 0, 0, 0 )*(2*li+1)*zi**li*exp( -zi*sig_0(li) )
        end do
      end do

      do mo = 0,lo; do me = -le,le
        if( abs(mo+me)>ls ) cycle
        do li = 0,limax
          do l = max(abs(ls-li), abs(le-lo)), min(ls+li, le+lo)
            if( mod(ls+l+li,2)/=0 .OR. mod(le+l+lo,2)/=0 ) cycle
            integral(ls,le,me,mo) = integral(ls,le,me,mo) +integ(l, li ) &
              *symbol_3j(ls, l, li, mo+me, -me-mo, 0 ) &
              *symbol_3j(lo, l, le, mo, -me-mo, me)
          end do
        end do
      end do; end do

      integral(ls,le,:,:) = integral(ls,le,:,:)*zi**(-ls-le)*sqrt((2*lo+1)*(2*ls+1) &
        *(2*le+1._RP) )*exp( zi*( sig_b(le) +sig_a(ls) ) )

    end do; end do

  end subroutine dwb_integrals

  pure complex(RP) function tpw( n, l, m, e, ke, k)
    use constants ,only: pi
    use trigo ,only: cartez2spher
    use utils ,only: norm_fac
    use special_functions ,only: spherical_harmonic, fac
    integer, intent(in) :: n,l,m
    real(RP), intent(in) :: e,ke(3)
    real(RP), intent(in),optional :: k(3)

    real(RP) :: q(3)
    integer :: j
    real(RP) :: kem,thetae,phie,A
    complex(RP) :: kec

    if(present(k)) then
      q = k-ke
    else
      q = -ke
    endif

    call cartez2spher( q, kem, thetae, phie)
    a = kem**2 +e**2
    kec = cmplx( 0._RP, kem, kind=RP )

    tpw = 0._RP
    do j = 0, (n-l)/2
      tpw = tpw +(-0.25_RP*a/e**2)**j*fac(n-j)/( fac(j)*fac(n-l-2*j) )
    end do

    tpw = tpw*sqrt(2./pi)*norm_fac(e,n)*fac(n-l)*(kec/e)**l*(2.*e/a)**n /a &
      *spherical_harmonic(l, m, thetae, phie )

  end function tpw

  pure complex(RP) function tcw( n, l, m, e, alpha, ke, k)
    use constants ,only: pi
    use special_functions ,only: fac
    use utils ,only: norm_fac

    integer, intent(in) :: n, l, m
    real(RP), intent(in) :: e, alpha, ke(3), k(3)

    real(RP)    :: kem, km, a, aj1, ke_t(3)
    complex(RP) :: w, kec, ekec, gam(0:n), f21(0:n,0:n), alphac, w1m, kep, kp
    complex(RP) :: tmp_j, tmp_j1, tmp_m1, cst_j, cst_j1, cst_m1
    real(RP)    :: tmp_s, tmp_s3, tmp_s1, cst_s, cst_s3, cst_s1, cst_s2
    complex(RP), parameter :: zi = (0._RP, 1._RP)
    integer :: j, j1, ma, m1, s, s1, s2, s3

    ma = abs(m)
    kem = norm2(ke)
    ke_t = -ke
    km = norm2(k)
    if( m>=0 ) then
      kp = cmplx( k(1), k(2), kind=RP )
      kep = cmplx( ke_t(1), ke_t(2), kind=RP )
    else
      kp = cmplx( k(1), -k(2), kind=RP )
      kep = cmplx( ke_t(1), -ke_t(2), kind=RP )
    end if

    alphac = cmplx( 0._RP, alpha, kind=RP) ! i*\alpha
    kec    = cmplx( 0._RP, kem  , kind=RP) ! i*ke
    ekec   = cmplx( e    , -kem , kind=RP) ! (\epsilon-ike)
    a      = norm2(k+ke_t)**2 +e**2
    w      = cmplx( norm2(k+ke_t)**2 -km**2 +kem**2   , 2.*e*kem , kind=RP) /a ! w = b/a
    w1m    = cmplx( e**2 +km**2 -kem**2, -2.*e*kem, kind=RP) /a ! (1-w)

    gam(0) = 1._RP
    do j1 = 1, n
      gam(j1) = gam(j1-1)*cmplx( j1, -alpha, kind=RP)
    end do

    do j1 =  0,n
      aj1 = j1 +1._RP
      f21(j1, 0) = 1._RP
      f21(j1, 1) = 1. +alphac*w / ( aj1*w1m )
      if( n<2 .OR. j1>n-2 ) cycle
      do j = 1, n -j1 -1
        f21( j1,  j +1 ) = (1. +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      end do
    end do

    tcw = 0._RP
    cst_j = -0.25*a/ekec**2
    cst_j1 = kec/ekec
    cst_m1 = kep/kp
    cst_s = -(km/k(3))**2
    cst_s3 = ke_t(3)/k(3)
    cst_s1 = kem/km
    cst_s2 = 2.*dot_product(k,ke_t)
    do j = 0,(n-l)/2
      tmp_j = fac(n-j) /fac(j)*cst_j**j
      !tmp_j = fac(n-j)/fac(j)*(-0.25*a/ekec**2)**j
      do j1 = 0,(n-l-2*j)
        !tmp_j1 = (kec/ekec)**j1 /( fac(j1)*fac(n-l-2*j-j1) )*tmp_j
        tmp_j1 = cst_j1**j1 /( fac(j1)*fac(n-l-2*j-j1) )*tmp_j
        do m1 = 0,ma
          tmp_m1 = cst_m1**m1 /( fac(ma-m1)*fac(m1) )*tmp_j1
          !tmp_m1 = (kep/kp)**m1 /( fac(ma-m1)*fac(m1) )*tmp_j1
          do s = 0,(l-ma)/2
            tmp_s = cst_s**s /( fac(l-s) )
            !tmp_s = (-1)**s*(km/k(3))**(2*s) /( fac(l-s) )
            do s3 = 0,l-ma-2*s
              tmp_s3 = cst_s3**s3 /( fac(s3)*fac(l-ma-2*s-s3) )*tmp_s
              !tmp_s3 = (ke_t(3)/k(3))**s3 /( fac(s3)*fac(l-ma-2*s-s3) )*tmp_s
              do s1 = 0,s
                tmp_s1 = cst_s1**s1 /fac(s-s1)*tmp_s3
                !tmp_s1 = (kem/km)**s1 /fac(s-s1)*tmp_s3
                do s2 = 0,s1
                  tcw = tcw +gam(m1+s3+2*s1-s2+j1)*f21(m1+s3+2*s1-s2+j1, n -j &
                    -(m1+s3+2*s1-s2+j1) )/( fac(s2)*fac(s1-s2)*fac(m1+s3+2*s1-s2+j1) ) &
                    *tmp_m1*tmp_s1*cst_s2**s2
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    tcw = tcw*norm_fac(e, n )*sqrt(l+0.5_RP)*fac(n-l)*sqrt( fac(l-ma)/fac(l+ma) ) &
      *(-1)**ma*fac(ma)*(zi*k(3)/ekec)**l*(2.*ekec/a)**n/a*(kp/k(3))**ma &
      *w1m**(-alphac) /pi
      !*powcc(w1m,-alpha)/pi
    if( mod(m,2)<0 ) tcw = -tcw

  end function tcw

  pure complex(RP) function tcw0( n, l, m, e, alpha, ke)
    use constants ,only: pi
    use trigo ,only: cartez2spher
    use utils ,only: norm_fac
    use special_functions ,only: spherical_harmonic, fac
    integer, intent(in) :: n, l, m
    real(RP), intent(in) :: e, alpha, ke(3)

    real(RP) :: kem, thetae, phie, a, aj1
    complex(RP) :: w, kec, ekec, gam(0:n), f21(0:n,0:(n-l)), alphac, w1m, tmp
    integer :: j, j1

    call cartez2spher( -ke, kem, thetae, phie)

    alphac = cmplx( 0._RP, alpha, kind=RP) ! i*\alpha
    kec    = cmplx( 0._RP, kem  , kind=RP) ! i*ke
    ekec   = cmplx( e    , -kem , kind=RP) ! (\epsilon-ike)
    a      = kem**2 +e**2
    w      = cmplx( 2.*kem**2   , 2.*e*kem , kind=RP) /a ! w = b/a
    w1m    = cmplx( e**2 -kem**2, -2.*e*kem, kind=RP) /a ! (1-w)

    gam(0) = 1._RP
    do j1 = 1, n
      gam(j1) = gam(j1-1)*cmplx( j1, -alpha, kind=RP)
    end do

    do j1 =  0,n-l
      aj1 = l +j1 +1._RP
      f21(j1, 0) = 1._RP
      f21(j1, 1) = 1._RP +alphac*w / ( aj1*w1m )
      if(n-l<2) cycle
      do j = 1, n -l -1
        f21( j1,  j +1 ) = (1._RP +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      end do
    end do

    tcw0 = 0._RP
    do j = 0,(n-l)/2
      tmp = fac(n-j) /fac(j)*(-0.25*a/ekec**2)**j
      do j1 = 0,(n-l-2*j) !, j1+2*j <= n-l )
        tcw0 = tcw0 +1./( fac(j1)*fac(n-l-2*j-j1)*fac(l+j1) )*(kec/ekec)**j1*gam(l+j1) &
          *f21(j1, n -l -j -j1 )*tmp !*fac(n-j)/ fac(j)*(-0.25*a/ekec**2)**j
      end do
    end do

    tcw0 = tcw0*norm_fac(e, n )*sqrt(2./pi)*fac(n-l)*(kec/ekec)**l*(2.*ekec/a)**n/a &
      *spherical_harmonic(l, m, thetae, phie)*w1m**(-alphac)!*powcc(w1m,-alpha)

  end function tcw0

  elemental complex(RP) function powcc(z1, y2)
    complex(RP), intent(in) :: z1
    real(RP), intent(in) :: y2
    real(RP) :: theta,zm

    theta = atan2( aimag(z1), real(z1, kind=RP) )
    zm = log( abs(z1) )
    powcc = exp(-y2*theta)*cmplx( cos(y2*zm), sin(y2*zm), kind=RP )

  end function powcc

end submodule fdcs_e2e
