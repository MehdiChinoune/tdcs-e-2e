program test_eikr_element
  use constants, only: wp, pi, deg
  use tcs_e2e, only: eikr_element
  use special_functions, only: spherical_harmonic
  use trigo, only: spher2cartez
  use integration, only: gauleg
  use utils, only: norm_fac
  implicit none
  !
  integer :: n1, n2, l1, l2, m1, m2
  real(wp) :: e(6) = [0.5_wp, 1._wp, 2._wp, 5._wp, 10._wp, 20._wp]
  real(wp) :: km(6) = [0.1_wp, 0.2_wp, 0.5_wp, 1._wp, 2._wp, 5._wp]
  real(wp) :: k(3)
  !
  integer :: ie1, ie2, ik, ikt
  integer, parameter :: nr = 32, nt = 64, nph = 96
  real(wp), allocatable :: xr(:), wr(:), xt(:), wt(:), xph(:), wph(:)
  real(wp) :: r(3), kr
  integer :: ir, it, iph
  complex(wp) :: integ
  !
  call gauleg(0._wp, pi, xt, wt, nt)
  call gauleg(0._wp, 2._wp*pi, xph, wph, nph)
  !
  do ik = 1, size(km); do ikt = 0, 4
    call spher2cartez(km(ik), ikt*pi/4, 0._wp, k)
    do ie1 = 1, size(e); do ie2 = 1, size(e)
      do n1 = 1, 2; do n2 = 1, 2
        call gauleg(0._wp, ((n1+n2)/2._wp+2)*5._wp/(e(ie1)+e(ie2)), xr, wr, nr)
        do l1 = 0, n1-1; do l2 = 0, n2-1
          do m1 = -l1, l1; do m2 = -l2, l2
            integ = (0._wp, 0._wp)
            !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(r, kr) REDUCTION(+:integ)
            do ir = 1, nr; do it = 1, nt; do iph = 1, nph
              call spher2cartez(xr(ir), xt(it), xph(iph), r)
              kr = dot_product(k, r)
              integ = integ + norm_fac(e(ie1), n1)* norm_fac(e(ie2), n2) * &
                cmplx(cos(kr), sin(kr), wp) * exp(-(e(ie1)+e(ie2))*xr(ir)) * &
                conjg( spherical_harmonic(l1, m1, xt(it), xph(iph)) ) * &
                spherical_harmonic(l2, m2, xt(it), xph(iph)) * &
                xr(ir)**(n1+n2) * sin(xt(it)) * wr(ir) * wt(it) * wph(iph)
            end do; end do; end do
            if( abs((eikr_element(n1, l1, m1, e(ie1), n2, l2, m2, e(ie2), k)-integ)) &
                > 1.e-3_wp ) then
                print*, nr, km(ik), (ikt*pi/4)/deg
                print*, n1, l1, m1, e(ie1)
                print*, n2, l2, m2, e(ie2)
                print*, integ
                print*, eikr_element(n1, l1, m1, e(ie1), n2, l2, m2, e(ie2), k)
            end if
          end do; end do
        end do; end do
      end do; end do
    end do; end do
  end do; end do
end program
