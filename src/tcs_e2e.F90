submodule(tcs_e2e) tcs_e2e
  implicit none
  real(wp), parameter :: prec = 0.001_wp
contains

  module subroutine tcs_fba_pw(in_unit)
    use constants, only : pi, ev
    use special_functions, only : factorial
    use input, only : read_tcs_input, read_orbit
    use trigo, only : spher2cartez
    use fdcs_e2e, only: tpw
    use integration, only: gauleg
    use types, only: orbit
    !
    integer, intent(in) :: in_unit
    !
    character(len=2) :: Atom_name, Orbit_name
    integer :: out_unit, exchange
    type(orbit) :: orbit_target

    real(wp) :: Ei, Es, Ee
    real(wp) :: kim, ksm, kem, km
    real(wp) :: ki(3), ks(3), ke(3), k(3)
    real(wp) :: Ei_values(46) = &
      [14, 15, 16, 18, 20, 22, 25, 28, 30, 32, 35, 40, 45, 50, 55, 60, 65, 70, &
      75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200, 250, 300, 400, &
      500, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000]
    real(wp) :: Ee_max
    integer :: ne = 32, nts = 16, nte = 8, nps = 16, npe = 8
    integer :: iei, iee, its, ips, ite, ipe
    real(wp), allocatable :: xe(:), we(:), xts(:), wts(:), xte(:), wte(:), &
      xps(:), wps(:), xpe(:), wpe(:)

    real(wp) :: factor, sigma, tcs
    complex(wp) :: D_term
    integer :: io, mo

    real(wp) :: k2(3), k2m = 0._wp
    complex(wp) :: E_term = (0._wp, 0._wp)

    call factorial()

    call read_tcs_input(in_unit, Atom_name, Orbit_name, exchange)

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    open( newunit = out_unit, file = 'tcs_pw_'//trim(Atom_name)//'_'//Orbit_name//'.dat', &
      status = 'replace', action = 'write' )
    write( out_unit, * ) "Ei TCS_PW"

    call gauleg(0._wp, pi, xts, wts, nts)
    call gauleg(0._wp, pi, xte, wte, nte)
    call gauleg(0._wp, 2._wp*pi, xps, wps, nps)
    call gauleg(0._wp, 2._wp*pi, xpe, wpe, npe)

    do iei = 1, size(Ei_values)

      Ei = Ei_values(iei)*Ev
      kim = sqrt(2.*Ei)
      call spher2cartez( kim, 0._wp, 0._wp, ki )

      Ee_max = (Ei - orbit_target%Ie)/2
      if(Ee_max<0._wp) cycle

      call gauleg(0._wp, Ee_max, xe, we, ne)
      tcs = 0._wp

      do iee = 1, ne

        Ee = xe(iee)
        kem = sqrt(2.*Ee)

        Es = Ei - orbit_target%Ie - Ee
        ksm = sqrt(2.*Es)

        factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp)*4._wp*ksm*kem/kim

        !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:tcs) &
        !$OMP PRIVATE(ks, k, km, ite, ipe, ke, k2, k2m, sigma, mo, D_term, io, E_term)
        do its = 1, nts; do ips = 1, nps
          call spher2cartez( ksm, xts(its), xps(ips), ks )
          k = ki -ks
          km = norm2(k)

          do ite = 1, nte; do ipe = 1, npe
            call spher2cartez( kem, xte(ite), xpe(ipe), ke )

            if(exchange==1) then
              k2 = ki -ke
              k2m = norm2(k2)
            end if

            sigma = 0._wp
            do mo = 0, orbit_target%l

              D_term = (0._wp,0._wp)
              do io = 1, orbit_target%nf
                D_term = D_term +orbit_target%a(io)*( tpw(orbit_target%n(io), orbit_target%l, &
                  mo, orbit_target%e(io), ke, k ) &
                  -tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ke ) )
              end do

              if(exchange==1) then
                E_term = (0._wp,0._wp)
                do io = 1, orbit_target%nf
                  E_term = E_term &
                    +orbit_target%a(io)*( tpw(orbit_target%n(io), orbit_target%l, mo, &
                    orbit_target%e(io), ks, k2 ) &
                    -tpw(orbit_target%n(io), orbit_target%l, mo, orbit_target%e(io), ks ) )
                end do
                sigma = sigma +(1+mo)*( abs(D_term/km**2)**2 +abs(E_term/k2m**2)**2 &
                  -real( D_term*conjg(E_term)/(km**2*k2m**2 ), wp ) )
              else
                sigma = sigma +(1+mo)*abs(D_term/km**2)**2
              end if

            end do

            tcs = tcs + factor*sigma * we(iee) * sin(xts(its))*wts(its)*wps(ips) * &
              sin(xte(ite))*wte(ite)*wpe(ipe)

          enddo; enddo
        enddo; enddo
      enddo

      if( show_output ) print '(1x,f6.0,1x,es15.8)', Ei_values(iei), tcs
      write( out_unit, '(1x,f6.0,1x,es15.8)' ) Ei_values(iei), tcs
    end do
    !
    close(out_unit)
    !
  end subroutine tcs_fba_pw

  module subroutine tcs_fba_cw(in_unit)
    use constants, only : pi, ev
    use types, only: atom, orbit
    use special_functions, only : factorial
    use input, only : read_tcs_input, read_orbit
    use trigo, only : spher2cartez
    use fdcs_e2e, only: tcw, tcw0
    use integration, only: gauleg, clenshaw_curtis
    !
    integer, intent(in) :: in_unit

    character(len=2) :: Atom_name, Orbit_name
    type(orbit) :: orbit_target
    integer :: out_unit, exchange, iei
    real(wp) :: Ei, kim, ki(3), zs, ze, Ee_max, tcs
    real(wp) :: Ei_values(55) = &
      [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 30, 32 &
      , 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140 &
      , 150, 160, 180, 200, 250, 300, 400, 500, 600, 800, 1000, 1500, 2000 &
      , 2500, 3000, 4000, 5000, 6000, 8000, 10000]
    !
    call factorial()

    call read_tcs_input(in_unit, Atom_name, Orbit_name, exchange)

    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )

    zs = -1 ! Charge of the projectile
    ze = -1 ! Charge of ejected particle

    open( newunit = out_unit, file = 'tcs_cw_'//trim(Atom_name)//'_'//Orbit_name//'.dat', &
      status = 'replace', action = 'write' )
    write( out_unit, * ) "Ei TCS_CW"

    do iei = 1, size(Ei_values)
      !
      Ei = Ei_values(iei)*eV
      kim = sqrt(2.*Ei)
      call spher2cartez( kim, 0._wp, 0._wp, ki )
      !
      Ee_max = (Ei - orbit_target%Ie)/2
      if(Ee_max<0._wp) cycle
      !
      call int_adaptive_cw([0._wp,0._wp,0._wp,0._wp],[Ee_max,pi,pi,2._wp*pi], ze, zs, &
        orbit_target, kim, ki, Ei, tcs, 1.e-3_wp)
      !
      if( show_output ) print '(1x,f6.0,1x,es15.8)', Ei_values(iei), tcs
      write( out_unit, '(1x,f6.0,1x,es15.8)' ) Ei_values(iei), tcs
      !
    end do
    !
    close(out_unit)
    !
  end subroutine tcs_fba_cw

  recursive subroutine int_adaptive_cw(a, b, ze, zs, orbit_target, &
      kim, ki, Ei, integ, tol)
    use constants, only: wg, wk, xk, pi
    use fdcs_e2e, only: tcw, tcw0
    use trigo, only: spher2cartez
    use types, only: orbit
    !
    real(wp), intent(in) :: a(4), b(4), tol
    real(wp), intent(out) :: integ
    type(orbit), intent(in) :: orbit_target
    real(wp), intent(in) :: ze, zs, kim, ki(3), Ei
    !
    real(wp) :: Ee, ke(3), kem, alpha, Es, ksm, sigma, factor, ks(3), k(3), km
    complex(wp) :: term
    real(wp) :: integ_g, integ_k(4), integ_1, integ_2, b1(4), a2(4), atol
    integer :: i1, i2, i3, i4, iee, its, ite, ips, io, mq, i_d
    real(wp) :: wgee(7), wgts(7), wgte(7), wgps(7), wkee(15), wkts(15), wkte(15), &
      wkps(15), xkee(15), xkts(15), xkte(15), xkps(15), funct(15,15,15,15)
    !
    wgee = (b(1)-a(1))*wg/2._wp
    wgts = (b(2)-a(2))*wg/2._wp
    wgte = (b(3)-a(3))*wg/2._wp
    wgps = (b(4)-a(4))*wg/2._wp
    !
    wkee = (b(1)-a(1))*wk/2._wp
    wkts = (b(2)-a(2))*wk/2._wp
    wkte = (b(3)-a(3))*wk/2._wp
    wkps = (b(4)-a(4))*wk/2._wp
    !
    xkee = ( (a(1)-b(1))*xk +b(1)+a(1) )/2._wp
    xkts = ( (a(2)-b(2))*xk +b(2)+a(2) )/2._wp
    xkte = ( (a(3)-b(3))*xk +b(3)+a(3) )/2._wp
    xkps = ( (a(4)-b(4))*xk +b(4)+a(4) )/2._wp
    !
    funct = 0._wp
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE(Ee, Es, kem, ke, alpha, ksm, ks, factor, km, k, sigma, mq, term, io)
    do iee = 1, 15; do its = 1, 15; do ips = 1, 15; do ite = 1, 15
      Ee = xkee(iee)
      kem = sqrt(2.*Ee)
      call spher2cartez( kem, xkte(ite), 0._wp, ke )
      alpha = -ze/kem
      !
      Es = Ei - orbit_target%Ie - Ee
      ksm = sqrt(2.*Es)
      call spher2cartez( ksm, xkts(its), xkps(ips), ks )
      !
      factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) &
        *4._wp*ksm*kem/kim &
        *2.*pi*alpha/(1._wp-exp(-2.*pi*alpha))
      if(factor /= factor) cycle
      !
      k = ki - ks
      km = norm2(k)
      !
      sigma = 0._wp
      do mq = 0, orbit_target%l
        term = (0._wp,0._wp)
        do io = 1, orbit_target%nf
          term = term +orbit_target%a(io)*( tcw(orbit_target%n(io), orbit_target%l, &
            mq, orbit_target%e(io), alpha, ke, k ) &
            +zs*tcw0(orbit_target%n(io), orbit_target%l, mq, orbit_target%e(io), alpha, ke ) )
        end do
        sigma = sigma +(mq+1)*abs(term)**2
      end do
      !
      funct(iee, its, ite, ips) = 2._wp*pi*factor*sigma/km**4*sin(xkts(its))*sin(xkte(ite))
    enddo; enddo; enddo ;enddo
    integ_g = 0._wp
    integ_k = 0._wp
    !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:integ_g,integ_k) &
    !$OMP PRIVATE(i1)
    do i2 = 2,15,2; do i3 = 2,15,2; do i4 = 2,15,2
      do i1=2,15,2
        integ_g = integ_g + funct(i1, i2, i3, i4)*wgee(i1/2)*wgts(i2/2)*wgte(i3/2)*wgps(i4/2)
      end do
      do i1= 1, 15
        integ_k(1) = integ_k(1) + funct(i1, i2, i3, i4)*wkee(i1)*wgts(i2/2)*wgte(i3/2)*wgps(i4/2)
        integ_k(2) = integ_k(2) + funct(i2, i1, i3, i4)*wgee(i2/2)*wkts(i1)*wgte(i3/2)*wgps(i4/2)
        integ_k(3) = integ_k(3) + funct(i3, i2, i1, i4)*wgee(i3/2)*wgts(i2/2)*wkte(i1)*wgps(i4/2)
        integ_k(4) = integ_k(4) + funct(i4, i2, i3, i1)*wgee(i4/2)*wgts(i2/2)*wgte(i3/2)*wkps(i1)
      end do
    end do; end do; end do
    !
    i_d = maxloc(abs(integ_g-integ_k),1)
    if( abs(integ_g)<=epsilon(1._wp) .AND. abs(integ_k(i_d))<=epsilon(1._wp) ) then
      atol = 0._wp
    elseif( abs(integ_g) <= epsilon(1._wp) ) then
      atol = abs(integ_g-integ_k(i_d))
    else
      atol = abs(integ_g-integ_k(i_d))/abs(integ_g)
    end if
    !
    if(atol<tol) then
      integ = integ_k(i_d)
    else
      b1 = b
      b1(i_d) = (a(i_d)+b(i_d))/2
      call int_adaptive_cw(a, b1, ze, zs, orbit_target, &
        kim, ki, Ei, integ_1, tol)
      a2 = a
      a2(i_d) = (a(i_d)+b(i_d))/2
      call int_adaptive_cw(a2, b, ze, zs, orbit_target, &
        kim, ki, Ei, integ_2, tol)
      integ = integ_1 + integ_2
    endif

  end subroutine int_adaptive_cw

  module subroutine tcs_fba_ocw(in_unit)
    use constants, only : pi, ev
    use types, only: atom, orbit
    use special_functions, only : factorial
    use input, only : read_tcs_input, read_orbit
    use trigo, only : spher2cartez
    use fdcs_e2e, only: tcw, tcw0
    use integration, only: gauleg, clenshaw_curtis
    !
    integer, intent(in) :: in_unit
    !
#define OCW
    character(len=2) :: Atom_name, Orbit_name
    type(orbit) :: orbit_target
    integer :: out_unit, exchange, iei
    real(wp) :: Ei, kim, ki(3), zs, ze, Ee_max, tcs
#ifdef OCW
    type(orbit) :: tmp_orbit
    type(atom) :: atom_target
    integer :: i, inp
    character(len=2) :: orbit_i
#endif
    real(wp) :: Ei_values(55) = &
      [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 30, 32 &
      , 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140 &
      , 150, 160, 180, 200, 250, 300, 400, 500, 600, 800, 1000, 1500, 2000 &
      , 2500, 3000, 4000, 5000, 6000, 8000, 10000]
    !
    call factorial()
    !
    call read_tcs_input(in_unit, Atom_name, Orbit_name, exchange)
    !
#ifdef OCW
    open( newunit = inp, file = 'data/'//trim(Atom_name)//'.dat', status = 'old', action = 'read')
    i = 0
    allocate(atom_target%orbits(0))
    do
      read(inp, fmt=*, iostat = iei ) orbit_i
      if(iei<0) exit
      i = i + 1
      call read_orbit(trim(Atom_name)//'_'//orbit_i, tmp_orbit)
      atom_target%orbits = [ atom_target%orbits, tmp_orbit]
    end do
#endif
    !
    call read_orbit(trim(Atom_name)//'_'//Orbit_name, orbit_target )
    !
    zs = -1 ! Charge of the projectile
    ze = -1 ! Charge of ejected particle
    !
#ifdef SOCW
    open( newunit = out_unit, file = 'tcs_socw_'//trim(Atom_name)//'_'//Orbit_name//'.dat', &
      status = 'replace', action = 'write' )
    write( out_unit, * ) "Ei TCS_SOCW"
#else
    open( newunit = out_unit, file = 'tcs_ocw_'//trim(Atom_name)//'_'//Orbit_name//'.dat', &
      status = 'replace', action = 'write' )
    write( out_unit, * ) "Ei TCS_OCW"
#endif

    do iei = 1, size(Ei_values)
      !
      Ei = Ei_values(iei)*eV
      kim = sqrt(2.*Ei)
      call spher2cartez( kim, 0._wp, 0._wp, ki )
      !
      Ee_max = (Ei - orbit_target%Ie)/2
      if(Ee_max<0._wp) cycle
#ifdef SOCW
      call int_adaptive_ocw([0._wp,0._wp,0._wp,0._wp],[Ee_max,pi,pi,2._wp*pi], ze, zs, &
        orbit_target, kim, ki, Ei, tcs, 1.e-3_wp)
#else
      call int_adaptive_ocw([0._wp,0._wp,0._wp,0._wp],[Ee_max,pi,pi,2._wp*pi], ze, zs, &
        orbit_target, atom_target, kim, ki, Ei, tcs, 1.e-3_wp)
#endif
      !
      if( show_output ) print '(1x,f6.0,1x,es15.8)', Ei_values(iei), tcs
      write( out_unit, '(1x,f6.0,1x,es15.8)' ) Ei_values(iei), tcs
      !
    end do
    !
    close(out_unit)
    !
  end subroutine tcs_fba_ocw

#ifdef SOCW
  recursive subroutine int_adaptive_ocw(a, b, ze, zs, orbit_target, &
    kim, ki, Ei, integ, tol)
#else
  recursive subroutine int_adaptive_ocw(a, b, ze, zs, orbit_target, atom_target, &
    kim, ki, Ei, integ, tol)
#endif
    use constants, only: wg, wk, xk, pi
    use types, only: orbit, atom
    use fdcs_e2e, only: tcw, tcw0
    use trigo, only: spher2cartez
    !
    real(wp), intent(in) :: a(4), b(4), tol
    real(wp), intent(out) :: integ
    type(orbit), intent(in) :: orbit_target
#ifdef OCW
    type(atom), intent(in) :: atom_target
#endif
    real(wp), intent(in) :: ze, zs, kim, ki(3), Ei
    !
    real(wp) :: Ee, ke(3), kem, alpha, Es, ksm, sigma, factor, ks(3), k(3), km
    complex(wp) :: term, term1, term2
#ifdef SOCW
    complex(wp), allocatable :: term_orth(:)
#else
    complex(wp) :: term_orth
#endif
    real(wp) :: integ_g, integ_k(4), integ_1, integ_2, b1(4), a2(4), atol
    integer :: i1, i2, i3, i4, iee, its, ite, ips, io, io2, mq, mq2, i_d, iorbits
    real(wp) :: wgee(7), wgts(7), wgte(7), wgps(7), wkee(15), wkts(15), wkte(15), &
      wkps(15), xkee(15), xkts(15), xkte(15), xkps(15), funct(15,15,15,15)
    !
    wgee = (b(1)-a(1))*wg/2._wp
    wgts = (b(2)-a(2))*wg/2._wp
    wgte = (b(3)-a(3))*wg/2._wp
    wgps = (b(4)-a(4))*wg/2._wp
    !
    wkee = (b(1)-a(1))*wk/2._wp
    wkts = (b(2)-a(2))*wk/2._wp
    wkte = (b(3)-a(3))*wk/2._wp
    wkps = (b(4)-a(4))*wk/2._wp
    !
    xkee = ( (a(1)-b(1))*xk +b(1)+a(1) )/2._wp
    xkts = ( (a(2)-b(2))*xk +b(2)+a(2) )/2._wp
    xkte = ( (a(3)-b(3))*xk +b(3)+a(3) )/2._wp
    xkps = ( (a(4)-b(4))*xk +b(4)+a(4) )/2._wp
    !
    funct = 0._wp
#ifdef SOCW
    allocate(term_orth(0:orbit_target%l))
#endif
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE(Ee, Es, kem, ke, alpha, ksm, ks, factor, km, k, sigma, term_orth, &
    !$OMP mq, mq2, term, term1, term2, io, io2, iorbits)
    do iee = 1, 15; do its = 1, 15; do ips = 1, 15; do ite = 1, 15
      Ee = xkee(iee)
      kem = sqrt(2.*Ee)
      call spher2cartez( kem, xkte(ite), 0._wp, ke )
      alpha = -ze/kem
      !
      Es = Ei - orbit_target%Ie - Ee
      ksm = sqrt(2.*Es)
      call spher2cartez( ksm, xkts(its), xkps(ips), ks )
      !
      factor = orbit_target%nelec/(2._wp*orbit_target%l+1._wp) &
        *4._wp*ksm*kem/kim &
        *2.*pi*alpha/(1._wp-exp(-2.*pi*alpha))
      if(factor /= factor) cycle
      !
      k = ki - ks
      km = norm2(k)
#ifdef SOCW
      term_orth = (0._wp, 0._wp)
      do mq = 0, orbit_target%l; do io = 1, orbit_target%nf; do io2 = 1, orbit_target%nf
        term_orth(mq) = term_orth(mq) + orbit_target%a(io)*orbit_target%a(io2) * &
          eikr_element(orbit_target%n(io), orbit_target%l, mq, orbit_target%e(io), &
          orbit_target%n(io2), orbit_target%l, mq, orbit_target%e(io2), k)
      end do; end do; end do
#endif
      !
      sigma = 0._wp
      do mq = 0, orbit_target%l
        term = (0._wp,0._wp)
        do io = 1, orbit_target%nf
#ifdef OCW
          term_orth = (0._wp,0._wp)
          do iorbits = 1, size(atom_target%orbits)
            associate (at => atom_target%orbits(iorbits))
              do mq2 = -at%l, at%l
                term1 = (0._wp, 0._wp)
                term2 = (0._wp, 0._wp)
                do io2 = 1, at%nf
                  term1 = term1 + at%a(io2) * &
                    eikr_element(at%n(io2), at%l, mq2, at%e(io2), orbit_target%n(io), &
                    orbit_target%l, mq, orbit_target%e(io), k)
                  term2 = term2 + at%a(io2)*tcw0(at%n(io2), at%l, mq2, at%e(io2), alpha, ke)
                end do
                term_orth = term_orth + term1*term2
              end do
            end associate
          end do
#endif

          term = term +orbit_target%a(io)*( tcw(orbit_target%n(io), orbit_target%l, &
            mq, orbit_target%e(io), alpha, ke, k ) &
#ifdef SOCW
            - term_orth(mq)*tcw0(orbit_target%n(io), orbit_target%l, mq, orbit_target%e(io), alpha, ke))
#else
            - term_orth)
#endif
        end do
        sigma = sigma +(mq+1)*abs(term)**2
      end do
      !
      funct(iee, its, ite, ips) = 2._wp*pi*factor*sigma/km**4*sin(xkts(its))*sin(xkte(ite))
    enddo; enddo; enddo ;enddo
    integ_g = 0._wp
    integ_k = 0._wp
    !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:integ_g,integ_k) &
    !$OMP PRIVATE(i1)
    do i2 = 2,15,2; do i3 = 2,15,2; do i4 = 2,15,2
      do i1=2,15,2
        integ_g = integ_g + funct(i1, i2, i3, i4)*wgee(i1/2)*wgts(i2/2)*wgte(i3/2)*wgps(i4/2)
      end do
      do i1= 1, 15
        integ_k(1) = integ_k(1) + funct(i1, i2, i3, i4)*wkee(i1)*wgts(i2/2)*wgte(i3/2)*wgps(i4/2)
        integ_k(2) = integ_k(2) + funct(i2, i1, i3, i4)*wgee(i2/2)*wkts(i1)*wgte(i3/2)*wgps(i4/2)
        integ_k(3) = integ_k(3) + funct(i3, i2, i1, i4)*wgee(i3/2)*wgts(i2/2)*wkte(i1)*wgps(i4/2)
        integ_k(4) = integ_k(4) + funct(i4, i2, i3, i1)*wgee(i4/2)*wgts(i2/2)*wgte(i3/2)*wkps(i1)
      end do
    end do; end do; end do
    !
    i_d = maxloc(abs(integ_g-integ_k),1)
    if( abs(integ_g)<=epsilon(1._wp) .AND. abs(integ_k(i_d))<=epsilon(1._wp) ) then
      atol = 0._wp
    elseif( abs(integ_g) <= epsilon(1._wp) ) then
      atol = abs(integ_g-integ_k(i_d))
    else
      atol = abs(integ_g-integ_k(i_d))/abs(integ_g)
    end if
    !
    if(atol<tol) then
      integ = integ_k(i_d)
    else
      b1 = b
      b1(i_d) = (a(i_d)+b(i_d))/2
      a2 = a
      a2(i_d) = (a(i_d)+b(i_d))/2
#ifdef SOCW
      call int_adaptive_ocw(a, b1, ze, zs, orbit_target, &
        kim, ki, Ei, integ_1, tol)
      call int_adaptive_ocw(a2, b, ze, zs, orbit_target, &
        kim, ki, Ei, integ_2, tol)
#else
      call int_adaptive_ocw(a, b1, ze, zs, orbit_target, atom_target, &
        kim, ki, Ei, integ_1, tol)
      call int_adaptive_ocw(a2, b, ze, zs, orbit_target, atom_target, &
        kim, ki, Ei, integ_2, tol)
#endif
      integ = integ_1 + integ_2
    endif

  end subroutine int_adaptive_ocw

  module complex(wp) function eikr_element(n1, l1, m1, e1, n2, l2, m2, e2, k)
    use constants, only: pi, fac
    use utils, only: norm_fac, y1y2y3
    use trigo, only: cartez2spher
    use special_functions, only: factorial, spherical_harmonic
    !
    integer, intent(in) :: n1, l1, m1, n2, l2, m2
    real(wp), intent(in) :: e1, e2, k(3)
    !
    integer :: n, l, j
    real(wp) :: e, km, theta, phi
    complex(wp), parameter :: xi = (0._wp, 1._wp)
    complex(wp) :: term_l, term_j
    !
    call factorial()
    !
    n = n1 + n2
    e = e1 + e2
    call cartez2spher(k, km, theta, phi)
    !
    term_l = (0._wp, 0._wp)
    do l = abs(l1-l2), l1 + l2
      if( mod(l1+l2+l,2)/=0 .OR. abs(m1-m2)>l ) cycle
      term_j = 0._wp
      do j = 0, (n-l-1)/2
        term_j = term_j + (-1)**j * (real(fac(n-l-1)*fac(n-j-1),wp)*(2._wp*e)**(n-l-2*j-1)) &
          / (real(fac(j)*fac(n-l-2*j-1),wp)*(e**2+km**2)**(n-j))
      enddo
      term_l = term_l + (2._wp*km*xi)**l*conjg(spherical_harmonic( l, m1-m2, theta, phi))* &
        (-1)**m1*y1y2y3(l1, l2, l, -m1, m2, m1-m2) * term_j
    enddo
    !
    eikr_element = 4._wp*pi*norm_fac(e1,n1)*norm_fac(e2,n2)*term_l
  end function eikr_element

end submodule
