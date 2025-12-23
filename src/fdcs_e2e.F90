submodule(fdcs_e2e) fdcs_e2e
  implicit none
contains

  module subroutine fdcs_fba_pw(in_unit,out_unit)
    use constants, only : pi, ev, deg
    use special_functions, only : factorial
    use input, only : read_fdcs_input, read_orbit
    use trigo, only : spher2cartez
    use types, only: orbit

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom_name, Orbit_name
    real(wp) :: Ei, Es, Ee
    real(wp) :: thetas
    integer :: step(3), exchange
    type(orbit) :: orbit_target

    real(wp), parameter   :: phie = 0._wp
    real(wp) :: kim, ksm, kem, km
    real(wp) :: ki(3), ks(3), ke(3), k(3)

    real(wp) :: factor, sigma
    complex(wp) :: D_term
    integer :: i, io, mo

    real(wp) :: k2(3), k2m = 0._wp
    complex(wp) :: E_term = (0._wp, 0._wp)

    call factorial()

    call read_fdcs_input(in_unit,Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange)

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    call spher2cartez( kim, 0._wp, 0._wp, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) *4._wp*ksm*kem/kim

    write( out_unit, * ) "Theta TDCS_PW"
    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      if(exchange==1) then
        k2 = ki -ke
        k2m = norm2(k2)
      end if

      sigma = 0._wp
      do mo = 0, orbit_target%l

        D_term = (0._wp,0._wp)
        do io = 1, orbit_target%nf
          D_term = D_term +orbit_target%a(io)*( tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ke, k ) &
            -tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ke ) )
        end do

        if(exchange==1) then
          E_term = (0._wp,0._wp)
          do io = 1, orbit_target%nf
            E_term = E_term &
              +orbit_target%a(io)*( tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ks, k2 ) &
              -tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ks ) )
          end do
          sigma = sigma +(1+mo)*( abs(D_term/km**2)**2 +abs(E_term/k2m**2)**2 &
            -real( D_term*conjg(E_term)/(km**2*k2m**2 ), wp ) )
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
    use constants, only : ev, deg, pi
    use trigo, only : spher2cartez
    use special_functions, only : factorial
    use input, only : read_fdcs_input, read_orbit
    use types, only: orbit
    !
    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom_name, Orbit_name
    real(wp) :: Ei, Es, Ee
    real(wp) :: thetas
    integer :: step(3), exchange
    integer :: ze, zs
    type(orbit) :: orbit_target

    real(wp), parameter   :: phie = 0._wp
    real(wp) :: kim, ksm, kem, km, alpha
    real(wp) :: ki(3), ks(3), ke(3), k(3)

    real(wp) :: factor, sigma
    complex(wp) :: term
    integer :: i, io, mo

    call factorial()

    call read_fdcs_input(in_unit,Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange)

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    zs = -1 ! Charge of the projectile
    kem = sqrt(2.*Ee*eV)
    ze = -1 ! Charge of ejected particle
    alpha = -ze/kem

    call spher2cartez( kim, 0._wp, 0._wp, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) *4._wp*ksm*kem / (kim*km**4)
    !\abs{ \exp{\pi\alpha/2}*\Gamma(1-i\alpha) }^2
    if(ze/=0) factor = factor*2.*pi*alpha/(1._wp-exp(-2.*pi*alpha))


    write( out_unit, * ) "Theta TDCS_CW"
    do i = step(1),step(2),step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._wp
      do mo = 0, orbit_target%l
        term = (0._wp,0._wp)
        do io = 1, orbit_target%nf
          term = term +orbit_target%a(io)*( tcw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), alpha, ke, k ) &
            +zs*tcw0(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), alpha, ke ) )
        end do
        sigma = sigma +(mo+1)*abs(term)**2
      end do

      sigma = factor*sigma

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_fba_cw

  module subroutine fdcs_fba_dw(in_unit,out_unit)
    use ieee_arithmetic, only : ieee_is_nan
    use constants, only : pi, deg, ev
    use special_functions, only : cgamma, spherical_harmonic, ricbes, factorial !&
      !, coul90, symbol_3j
    use utils, only : norm_fac, y1y2y3, calculate_U, ode_second_dw
    use input, only : read_fdcs_input, read_orbit
    use trigo, only : spher2cartez, cartez2spher
    use types, only: orbit

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom_name, Orbit_name
    real(wp) :: Ei, Es, Ee
    real(wp) :: thetas
    integer :: step(3), exchange
    type(orbit) :: orbit_target

    real(wp), parameter   :: phie = 0._wp
    real(wp) :: kim, ksm, kem, km
    real(wp) :: theta, phi
    real(wp) :: ki(3), ks(3), ke(3), k(3)

    integer, parameter :: lmax = 47, lemax = 15, np = 2400
    real(wp) :: x(0:np), wf(0:np)
    real(wp), allocatable :: chi_b(:,:), chi_0a(:,:)
    real(wp) :: sigma_le(0:lemax), delta(0:lemax), integral(0:lemax,0:lmax)
    real(wp) :: etae, rc, h
    complex(wp) :: tmp_z
    complex(wp), parameter :: zi = (0._wp, 1._wp)
    integer :: le, l, me, ze

    real(wp) :: factor, sigma
    complex(wp) :: term, term0
    integer :: i, io, mo

    real(wp) :: gc_0a(0:lmax), fdc_0a(0:lmax), gdc_0a(0:lmax), fc_0a(0:lmax)
    real(wp), allocatable :: U(:,:), U_tmp(:)

    call factorial()

    call read_fdcs_input(in_unit,Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange )

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)
    ze = 1

    call spher2cartez( kim, 0._wp, 0._wp, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    call cartez2spher( k, km, theta, phi)

    factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) *4._wp*ksm*kem / (kim*km**4)
    factor = factor*2./(pi*kem**2)

    rc = 10._wp
    h = rc/np

    etae = -ze/kem

    x = [(i*h, i = 0, np) ]
    allocate( chi_b(0:np,0:lemax), chi_0a(0:np,0:lmax) )
    chi_b(0,:) = 0._wp
    chi_0a(0,:) = 0._wp
    do i = 1,np
      call ricbes(km*x(i), lmax, fc_0a, gc_0a, fdc_0a, gdc_0a, le)
      chi_0a(i,:) = fc_0a/(km*x(i))
    end do

    where(ieee_is_nan(chi_0a) )
      chi_0a = 0._wp
    end where

    allocate( U(0:np,0:lemax), U_tmp(0:np) )
    call calculate_U(Atom_name, Orbit_name, x, U_tmp, 1 )

    if(ze/=0) then
      U(0,:) = -huge(1._wp)
    else
      U(0,:) = -kem**2
    end if

    U(1:np,lemax) = -kem**2 -2.*( ze +U_tmp(1:np) ) /x(1:np)
    do l = 0,lemax
      U(1:np,l) = l*(l+1)/x(1:np)**2 +U(1:np,lemax)
    end do

    call ode_second_dw(kem, lemax, rc, ze, U, chi_b, delta)

    deallocate(U,U_tmp)

    wf = 0._wp
    do io = 1, orbit_target%nf
      wf = wf +orbit_target%a(io)*norm_fac(orbit_target%e(io), orbit_target%n(io) )*x**orbit_target%n(io)*exp(-orbit_target%e(io)*x)
    end do

    do le = 0,lemax
      tmp_z = cgamma( cmplx(le+1._wp, etae, wp), 1 )
      tmp_z = exp( zi*aimag(tmp_z) )
      sigma_le(le) = atan2( aimag(tmp_z), real(tmp_z,wp) )
    end do

    integral = 0._wp
    do l = 0,lmax
      do le = 0,lemax
        if( mod(le+l+orbit_target%l,2)/=0 ) cycle
        integral(le, l) = 0.5*h*sum( wf(1:np)*chi_b(1:np,le)*chi_0a(1:np,l) &
          +wf(0:np-1)*chi_b(0:np-1,le)*chi_0a(0:np-1,l) )
      end do
    end do

    term0 = zi**(-orbit_target%l)*exp( cmplx(0._wp, -(sigma_le(orbit_target%l)+delta(orbit_target%l)), wp ) ) &
      *0.5*h*sum( wf(1:np)*chi_b(1:np,orbit_target%l) +wf(0:np-1)*chi_b(0:np-1,orbit_target%l) )


    write( out_unit, * ) "Theta TDCS_DW"
    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._wp
      do mo = 0, orbit_target%l

        term = (0._wp, 0._wp)
        do l = 0,lmax
          do le = 0,lemax
            if( mod(le+l+orbit_target%l,2)/=0 ) cycle
            tmp_z = (0._wp, 0._wp)
            do me = -le,le
              if( abs(me-mo)>l ) cycle
              tmp_z = tmp_z +y1y2y3(le, l, orbit_target%l, -me, me-mo, mo) &
                *spherical_harmonic(le, me, i*deg, phie) &
                *spherical_harmonic(l, mo-me, theta, phi)
            end do
            term = term +zi**(l-le)*integral(le, l)*tmp_z &
              *exp( cmplx(0._wp, -(sigma_le(le)+delta(le)), wp ) )

          end do
        end do
        term = 4.*pi*term

        sigma = sigma +(mo+1)*abs( term &
          -spherical_harmonic(orbit_target%l, mo, i*deg, phie)*term0 )**2

      end do

      sigma = factor*sigma

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_fba_dw

  module subroutine fdcs_dwb(in_unit,out_unit)
    use constants, only : ev, deg, pi
    use special_functions, only : spherical_harmonic, cgamma, factorial
    use utils, only : norm_fac, calculate_U
    use input, only : read_fdcs_input, read_orbit
    use trigo, only : spher2cartez
    use integration, only : gauleg
    use types, only: orbit

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom_name, Orbit_name
    real(wp) :: Ei, Es, Ee
    real(wp) :: thetas
    integer :: step(3), exchange, PCI = 2
    type(orbit) :: orbit_target

    real(wp), parameter   :: phie = 0._wp
    real(wp) :: kim, ksm, kem

    integer, parameter :: nx = 31*64, nr = 48000

    integer, parameter :: limax = 95
    real(wp), allocatable, target :: chi_0(:,:)
    real(wp), target :: delta_li(0:limax)

    integer, parameter :: lsmax = 95
    real(wp), allocatable, target :: chi_a(:,:)
    real(wp) :: sigma_ls(0:lsmax), etas
    real(wp), target :: delta_ls(0:lsmax)
    complex(wp) :: ylms(0:lsmax,-lsmax:lsmax)
    integer :: ls, zs = 1

    integer, parameter :: lemax = 23
    real(wp), allocatable, target :: chi_b(:,:)
    real(wp) :: sigma_le(0:lemax), etae
    real(wp), target :: delta_le(0:lemax)
    complex(wp), allocatable :: ylme(:,:,:)
    integer :: le, me, ze = 1

    real(wp), allocatable :: r(:), wf(:)
    real(wp) :: rc, h
    complex(wp), allocatable :: integral(:,:,:,:), integralx(:,:,:,:)
    complex(wp) :: tmp_z

    real(wp) :: factor, sigma
    complex(wp) :: termd, termx
    integer :: i, io, mo = 0

    real(wp), allocatable :: x(:), w(:)

    real(wp), allocatable :: U_tmp(:)
    real(wp), allocatable :: xt(:), wt(:)
    integer, parameter :: nc = nx/64

    call factorial()

    call read_fdcs_input(in_unit,Ei, Es, Ee, thetas, step, Atom_name, Orbit_name, exchange )

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    etas = -zs/ksm
    etae = -ze/kem

    factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) *(2.*pi)**4*ksm*kem / kim
    factor = factor*2._wp/pi**4 /(kem*ksm*kim)**2

    rc = 250._wp

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

    call calculate_U(Atom_name, Orbit_name, r, U_tmp, 0)
    call calculate_chi( kim, r, U_tmp, 0, x, chi_0, delta_li )

    call calculate_U(Atom_name, Orbit_name, r, U_tmp, 1)
    call calculate_chi( ksm, r, U_tmp, zs, x, chi_a, delta_ls )
    call calculate_chi( kem, r, U_tmp, ze, x, chi_b, delta_le )

    deallocate( U_tmp )

    allocate(wf(nx))
    wf = 0._wp
    do io = 1, orbit_target%nf
      wf = wf +orbit_target%a(io)*norm_fac(orbit_target%e(io), orbit_target%n(io) )*x**orbit_target%n(io)*exp(-orbit_target%e(io)*x)
    end do

    chi_b(:,orbit_target%l) = chi_b(:,orbit_target%l) -wf*sum( w*chi_b(:,orbit_target%l)*wf )
    chi_a(:,orbit_target%l) = chi_a(:,orbit_target%l) -wf*sum( w*chi_a(:,orbit_target%l)*wf )

    tmp_z = cgamma( cmplx( 1._wp, etae, wp), 0 )
    sigma_le(0) = atan2( aimag(tmp_z), real(tmp_z, wp) )
    do le = 0,lemax-1
      sigma_le(le+1) = sigma_le(le) +atan2( etae, le+1._wp )
    end do
    sigma_le = sigma_le +delta_le

    tmp_z = cgamma( cmplx( 1._wp, etas, wp), 0 )
    sigma_ls(0) = atan2( aimag(tmp_z), real(tmp_z, wp) )
    do ls = 0,lsmax-1
      sigma_ls(ls+1) = sigma_ls(ls) +atan2( etas, ls+1._wp )
    end do
    sigma_ls = sigma_ls +delta_ls

    call dwb_integrals(chi_0, chi_a, chi_b, delta_li, sigma_ls, sigma_le &
      , wf, x, w, orbit_target%l, integral)
    if(exchange==1) then
      call dwb_integrals(chi_0, chi_b, chi_a, delta_li, sigma_le, sigma_ls &
        , wf, x, w, orbit_target%l, integralx)
    end if

    ylms = (0._wp, 0._wp)
    do ls = 0,lsmax
      ylms(ls,0) = spherical_harmonic(ls,0,thetas*deg,pi)
      do me = 1,ls
        ylms(ls,me) = spherical_harmonic(ls,me,thetas*deg,pi)
      end do
      ylms(ls,-2:-ls:-2) = ylms(ls,2:ls:2)
      ylms(ls,-1:-ls:-2) = -ylms(ls,1:ls:2)
    end do

    allocate(ylme(0:lemax,-lemax:lemax,step(1)/step(3):step(2)/step(3)))
    ylme = (0._wp, 0._wp)
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
    write( out_unit, * ) "Theta TDCS_DWBA"

    !$OMP PARALLEL DO PRIVATE(sigma, mo, termd, termx, le, ls, me)
    do i = step(1), step(2), step(3)

      sigma = 0._wp
      do mo = 0, orbit_target%l

        termd = (0._wp, 0._wp)
        do le = 0,lemax
          do ls = 0,lsmax
            do me = -le,le
              if( abs(me+mo)>ls ) cycle
              termd = termd +ylme(le,me,i/step(3))*ylms(ls,mo+me) &
                *integral(ls, le, me, mo )
            end do
          end do
        end do

        termx = (0._wp, 0._wp)
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
          -real( conjg(termd)*termx, wp ) )
      end do

      sigma = factor*sigma
      if(PCI>=1) then
        call PCI_EFFECTS(i,sigma)
      end if

      if( show_output ) print'(1x,i4,1x,es15.8)', i, sigma
      write(out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  contains

    elemental subroutine PCI_EFFECTS(i,sigma)
      use special_functions, only : conhyp_opt
      integer, intent(in) :: i
      real(wp), intent(inout) :: sigma
      real(wp) :: kesm
      real(wp) :: r12_av, Et
      kesm = sqrt( kem**2 +ksm**2 -2*kem*ksm*cos( (thetas+i)*deg ) )/2._wp
      sigma = pi/(kesm*(exp(pi/kesm) -1._wp ) )*sigma
      if(PCI==2) then
        Et = (Es +Ee)*eV
        r12_av = pi**2/(16.*Et)*( 1._wp +(0.627/pi)*sqrt(Et)*log(Et) )**2
        sigma = abs( conhyp_opt( 0.5/kesm, -2.*kesm*r12_av  ) )**2*sigma
      end if
    end subroutine

  end subroutine fdcs_dwb

  module subroutine fdcs_bbk(in_unit,out_unit)
    use constants, only : pi, ev, deg
    use special_functions, only : factorial
    use input, only : read_fdcs_input, read_orbit
    use trigo, only : spher2cartez
    use types, only: orbit

    integer, intent(in) :: in_unit
    integer, intent(in) :: out_unit

    character(len=2) :: Atom_name, Orbit_name
    real(wp) :: Ei, Es, Ee
    real(wp) :: thetas
    integer :: step(3)
    type(orbit) :: orbit_target

    real(wp), parameter   :: phie = 0._wp
    real(wp) :: kim, ksm, kem, km
    real(wp) :: ki(3), ks(3), ke(3), k(3)

    real(wp) :: factor, sigma
    complex(wp) :: D_term
    integer :: i, io, mo

    call factorial()

    call read_fdcs_input(in_unit,Ei, Es, Ee, thetas, step, Atom_name, Orbit_name)

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    kim = sqrt(2.*Ei*eV)
    ksm = sqrt(2.*Es*eV)
    kem = sqrt(2.*Ee*eV)

    call spher2cartez( kim, 0._wp, 0._wp, ki )
    call spher2cartez( ksm, thetas*deg, pi, ks )
    k = ki -ks
    km = norm2(k)

    factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) *4._wp*ksm*kem/(kim)

    do i = step(1), step(2), step(3)

      call spher2cartez( kem, i*deg, phie, ke )

      sigma = 0._wp
      do mo = 0, orbit_target%l

        D_term = (0._wp,0._wp)
        do io = 1, orbit_target%nf
          D_term = D_term +orbit_target%a(io)*( tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ke, k ) &
            -tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ke ) )
        end do
        sigma = sigma +(1+mo)*abs(D_term/km**2)**2
      end do

      sigma = factor*sigma

      if( show_output ) print'(1x,i4,1x,es15.8)',i,sigma
      write( out_unit, '(1x,i4,1x,es15.8)' ) i, sigma
    end do

  end subroutine fdcs_bbk

  complex(wp) function U_bbk(alpha1, alpha2, alpha3, k1, k2, k3, lam1, lam2, lam3 &
    , p1, p2)
    use constants, only : pi
    use integration, only : gauleg

    real(wp), intent(in) :: alpha1, alpha2, alpha3
    real(wp), intent(in) :: k1(3), k2(3), k3(3)
    real(wp), intent(in) :: lam1, lam2, lam3
    real(wp), intent(in) :: p1(3), p2(3)

    complex(wp), parameter :: zi = (0._wp, 1._wp)

    real(wp),allocatable :: t3(:), y(:), wy(:), s(:), ws(:)
    real(wp) :: p11(3), p22(3), q1(3), q2(3)
    real(wp) :: q1m, q2m
    real(wp) :: k1m, k2m, k3m
    complex(wp) :: lam33
    complex(wp) :: sig0, sig1, sig2, sig12
    complex(wp) :: sig0_a(3), sig1_a(3), sig2_a(3), sig12_a(3)

    integer :: i, j
    complex(wp) :: N_t3, Integ

    U_bbk = ( 0._wp, 0._wp)

    k1m = norm2(k1)
    k2m = norm2(k2)
    k3m = norm2(k3)

    call gauleg(-10._wp, 10._wp, y, wy, 64)
    allocate( t3(size(y)) )
    t3 = 1._wp /(1._wp+ exp(y) )

    call gauleg(0._wp, 10._wp, s, ws, 64)

    Integ = (0._wp, 0._wp)
    do i = 1,64
      p11 = p1 -t3(i)*k3
      q1 = k1 +p11
      q1m = norm2(q1)
      p22 = p2 -t3(i)*k3
      q2 = k2 +p22
      q2m = norm2(q2)
      lam33 = cmplx(lam3, -t3(i)*k3m, wp )

      sig0_a = [ cmplx( (lam1+lam2)**2+norm2(q1-q2)**2, 0._wp, wp) &
        , 2*(lam2*(lam1**2+lam33**2+q1m**2) +lam1*(lam2**2+lam33**2+q2m**2) ) &
        , ( (lam1+lam33)**2+q1m**2)*((lam2+lam33)**2+q2m**2) ]
      sig1_a = -2*[ dot_product(q1-q2,k1) +zi*(lam1+lam2)*k1m &
        , 2*lam2*dot_product(q1,k1) +zi*( (lam2**2+lam33**2+q2m**2)*k1m+2*lam1*lam2*k1m ) &
        , ((lam2+lam33)**2 +q2m**2)*( dot_product(q1,k1) +zi*(lam1+lam33)*k1m ) ]
      sig2_a = -2*[ dot_product(q2-q1,k2) +zi*(lam1+lam2)*k2m &
        , 2*lam1*dot_product(q2,k2) +zi*( (lam1**2+lam33**2+q1m**2)*k2m+2*lam1*lam2*k2m ) &
        , ((lam1+lam33)**2 +q1m**2)*( dot_product(q2,k2) +zi*(lam2+lam33)*k2m ) ]
      sig12_a = -2*[ cmplx( k1m*k2m +dot_product(k1,k2), 0._wp, wp) &
        , 2*( (lam1*k1m-zi*dot_product(q1,k1))*k2m +(lam2*k2m-zi*dot_product(q2,k2))*k1m) &
        , -2*( dot_product(q1,k1)+zi*(lam1+lam33)*k1m) &
        *( dot_product(q2,k2)+zi*(lam2+lam33)*k2m) ]

      N_t3 = (0._wp, 0._wp)
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
    use utils, only : ode_second_dw, intrpl
    real(wp), intent(in) :: U_tmp(0:), r(0:), x(:), km
    integer, intent(in) :: z
    real(wp), intent(out) :: chi(:,0:), delta(0:)

    real(wp), allocatable :: U(:,:), chi_tmp(:,:)
    real(wp) :: rc
    integer :: l, lmax, nr

    lmax = size(delta) -1
    nr = size(r) -1
    rc = r(nr)

    allocate( U(0:nr,0:lmax), chi_tmp(0:nr,0:lmax) )

    if(z==0) then
      U(0,0:lmax) = -km**2
    else
      U(0,0:lmax) = -huge(1._wp)
    end if

    U(1:nr,lmax) = -km**2 -2.*( z*1._wp +U_tmp(1:nr) )/r(1:nr)
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
    use special_functions, only : symbol_3j
    real(wp), intent(in) :: chi_0(:,0:), chi_a(:,0:), chi_b(:,0:), wf(:)
    real(wp), intent(in) :: x(:), w(:)
    real(wp), intent(in) :: sig_0(0:), sig_a(0:), sig_b(0:)
    integer, intent(in) :: lo
    complex(wp), allocatable, intent(out) :: integral(:,:,:,:)

    real(wp), parameter :: VSmall = tiny(1._wp), VBig = huge(1._wp)
    complex(wp), parameter :: zi = (0._wp, 1._wp)
    integer :: limax, lsmax, lemax, nx, nmin
    complex(wp), allocatable :: integ(:,:)
    integer :: li, ls, le, l, me, mo
    ! Block variables
    real(wp) :: Ti, Si
    real(wp) :: integ0
    real(wp) :: xil, xil1
    integer :: i, is, nmax

    nx = size(x)
    limax = size(sig_0) -1
    lsmax = size(sig_a) -1
    lemax = size(sig_b) -1

    nmin = count(x<=25._wp)

    allocate( integral(0:lsmax, 0:lemax, -lemax:lemax, 0:lo) )
    allocate( integ( 0:min(limax+lsmax,lemax+lo),  0:limax ) )
    integral = 0._wp
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(integ,li,l,me,mo ,ti,si,integ0,xil,xil1,i,is,nmax)
    do ls = 0,lsmax; do le = 0,lemax
      integ = 0._wp
      do li = 0,limax
        do l = max(abs(ls-li), abs(le-lo)), min(ls+li, le+lo)
          if( mod(ls+l+li,2)/=0 .OR. mod(le+l+lo,2)/=0 ) cycle

          if(l<=5) nmax = nmin +(5-l)*(nx-nmin)/5

          Ti = 0._wp
          integ0 = 0._wp
          Si = 0._wp

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
        *(2*le+1._wp) )*exp( zi*( sig_b(le) +sig_a(ls) ) )

    end do; end do

  end subroutine dwb_integrals

  module pure complex(wp) function tpw( n, l, m, e, ke, k)
    use constants, only : pi
    use trigo, only : cartez2spher
    use utils, only : norm_fac
    use special_functions, only : spherical_harmonic, fac
    integer, intent(in) :: n,l,m
    real(wp), intent(in) :: e,ke(3)
    real(wp), intent(in),optional :: k(3)

    real(wp) :: q(3)
    integer :: j
    real(wp) :: kem,thetae,phie,A
    complex(wp) :: kec

    if(present(k)) then
      q = k-ke
    else
      q = -ke
    endif

    call cartez2spher( q, kem, thetae, phie)
    a = kem**2 +e**2
    kec = cmplx( 0._wp, kem, wp )

    tpw = 0._wp
    do j = 0, (n-l)/2
      tpw = tpw +(-0.25_wp*a/e**2)**j*fac(n-j)/( fac(j)*fac(n-l-2*j) )
    end do

    tpw = tpw*sqrt(2./pi)*norm_fac(e,n)*fac(n-l)*(kec/e)**l*(2.*e/a)**n /a &
      *spherical_harmonic(l, m, thetae, phie )

  end function tpw

  module pure complex(wp) function tcw( n, l, m, e, alpha, ke, k)
    use constants, only : pi
    use special_functions, only : fac
    use utils, only : norm_fac

    integer, intent(in) :: n, l, m
    real(wp), intent(in) :: e, alpha, ke(3), k(3)

    real(wp)    :: kem, km, a, aj1, ke_t(3)
    complex(wp) :: w, kec, ekec, gam(0:n), f21(0:n,0:n), alphac, w1m, kep, kp
    complex(wp) :: tmp_j, tmp_j1, tmp_m1, cst_j, cst_j1, cst_m1
    real(wp)    :: tmp_s, tmp_s3, tmp_s1, cst_s, cst_s3, cst_s1, cst_s2
    complex(wp), parameter :: zi = (0._wp, 1._wp)
    integer :: j, j1, ma, m1, s, s1, s2, s3

    ma = abs(m)
    kem = norm2(ke)
    ke_t = -ke
    km = norm2(k)
    if( m>=0 ) then
      kp = cmplx( k(1), k(2), wp )
      kep = cmplx( ke_t(1), ke_t(2), wp )
    else
      kp = cmplx( k(1), -k(2), wp )
      kep = cmplx( ke_t(1), -ke_t(2), wp )
    end if

    alphac = cmplx( 0._wp, alpha, wp) ! i*\alpha
    kec    = cmplx( 0._wp, kem  , wp) ! i*ke
    ekec   = cmplx( e    , -kem , wp) ! (\epsilon-ike)
    a      = norm2(k+ke_t)**2 +e**2
    w      = cmplx( norm2(k+ke_t)**2 -km**2 +kem**2   , 2.*e*kem , wp) /a ! w = b/a
    w1m    = cmplx( e**2 +km**2 -kem**2, -2.*e*kem, wp) /a ! (1-w)

    gam(0) = 1._wp
    do j1 = 1, n
      gam(j1) = gam(j1-1)*cmplx( j1, -alpha, wp)
    end do

    do j1 =  0,n
      aj1 = j1 +1._wp
      f21(j1, 0) = 1._wp
      f21(j1, 1) = 1. +alphac*w / ( aj1*w1m )
      if( n<2 .OR. j1>n-2 ) cycle
      do j = 1, n -j1 -1
        f21( j1,  j +1 ) = (1. +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      end do
    end do

    tcw = 0._wp
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

    tcw = tcw*norm_fac(e, n )*sqrt(l+0.5_wp)*fac(n-l)*sqrt( fac(l-ma)/fac(l+ma) ) &
      *(-1)**ma*fac(ma)*(zi*k(3)/ekec)**l*(2.*ekec/a)**n/a*(kp/k(3))**ma &
      *w1m**(-alphac) /pi
      !*powcc(w1m,-alpha)/pi
    if( mod(m,2)<0 ) tcw = -tcw

  end function tcw

  module pure complex(wp) function tcw0( n, l, m, e, alpha, ke)
    use constants, only : pi
    use trigo, only : cartez2spher
    use utils, only : norm_fac
    use special_functions, only : spherical_harmonic, fac
    integer, intent(in) :: n, l, m
    real(wp), intent(in) :: e, alpha, ke(3)

    real(wp) :: kem, thetae, phie, a, aj1
    complex(wp) :: w, kec, ekec, gam(0:n), f21(0:n,0:(n-l)), alphac, w1m, tmp
    integer :: j, j1

    call cartez2spher( -ke, kem, thetae, phie)

    alphac = cmplx( 0._wp, alpha, wp) ! i*\alpha
    kec    = cmplx( 0._wp, kem  , wp) ! i*ke
    ekec   = cmplx( e    , -kem , wp) ! (\epsilon-ike)
    a      = kem**2 +e**2
    w      = cmplx( 2.*kem**2   , 2.*e*kem , wp) /a ! w = b/a
    w1m    = cmplx( e**2 -kem**2, -2.*e*kem, wp) /a ! (1-w)

    gam(0) = 1._wp
    do j1 = 1, n
      gam(j1) = gam(j1-1)*cmplx( j1, -alpha, wp)
    end do

    do j1 =  0,n-l
      aj1 = l +j1 +1._wp
      f21(j1, 0) = 1._wp
      f21(j1, 1) = 1._wp +alphac*w / ( aj1*w1m )
      if(n-l<2) cycle
      do j = 1, n -l -1
        f21( j1,  j +1 ) = (1._wp +( j +alphac*w) / ( (aj1 +j)*w1m ) )*f21(j1, j) &
          -j/( (aj1 +j)*w1m )*f21(j1, j-1)
      end do
    end do

    tcw0 = 0._wp
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

  elemental complex(wp) function powcc(z1, y2)
    complex(wp), intent(in) :: z1
    real(wp), intent(in) :: y2
    real(wp) :: theta,zm

    theta = atan2( aimag(z1), real(z1, wp) )
    zm = log( abs(z1) )
    powcc = exp(-y2*theta)*cmplx( cos(y2*zm), sin(y2*zm), wp )

  end function powcc

end submodule fdcs_e2e
