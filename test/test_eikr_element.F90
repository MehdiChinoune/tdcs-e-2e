program test_eikr_element
  use constants, only: wp, pi, deg
  use tcs_e2e, only: eikr_element
  use trigo, only: spher2cartez
  implicit none
  !
  integer :: n1, n2, l1, l2, m1, m2
  real(wp) :: e(6) = [0.5_wp, 1._wp, 2._wp, 5._wp, 10._wp, 20._wp]
  real(wp) :: km(6) = [0.1_wp, 0.2_wp, 0.5_wp, 1._wp, 2._wp, 5._wp]
  real(wp) :: k(3)
  !
  integer :: ie1, ie2, ik, ikt
  real(wp) :: atol
  complex(wp) :: integ, eikr
  !
  !$OMP PARALLEL DO COLLAPSE(6) PRIVATE(k,l1,l2,m1,m2,integ,eikr,atol)
  do ik = 1, size(km); do ikt = 0, 4
    do ie1 = 1, size(e); do ie2 = 1, size(e)
      do n1 = 1, 2; do n2 = 1, 2
        call spher2cartez(km(ik), ikt*pi/4, 0._wp, k)
        do l1 = 0, n1-1; do l2 = 0, n2-1
          do m1 = -l1, l1; do m2 = -l2, l2
            integ = integrate([0._wp,0._wp,0._wp],[((n1+n2)/2._wp+2)*5._wp/(e(ie1)+e(ie2)),pi,2._wp*pi],&
              n1, l1, m1, e(ie1), n2, l2, m2, e(ie2), k, 1.e-4_wp)
            eikr=eikr_element(n1, l1, m1, e(ie1), n2, l2, m2, e(ie2), k)
            !
            if( abs(integ)<=epsilon(1._wp) .AND. abs(eikr)<=epsilon(1._wp) ) then
              atol = 0._wp
            elseif( abs(eikr) <= epsilon(1._wp) ) then
              atol = abs(eikr-integ)
            else
              atol = abs(eikr-integ)/abs(eikr)
            end if
            !
            if(atol>1.e-2_wp) then
              print*, n1, l1, m1, e(ie1)
              print*, n2, l2, m2, e(ie2)
              print*, km(ik), (ikt*pi/4)/deg
              print*, integ
              print*, eikr
              error stop "Failed: numerical value of <n_1l_1m_1|e^{ikr}|n_2l_2m_2> "&
                "is not equal to its anlaytical value"
            end if
          end do; end do
        end do; end do
      end do; end do
    end do; end do
  end do; end do
  !
contains
  recursive complex(wp) function integrate(a, b, n1, l1, m1, e1, n2, l2, m2, e2, k, tol) result(res)
    use constants, only: wp, wg, xk, wk
    use special_functions, only: spherical_harmonic
    use utils, only: norm_fac
    !
    real(wp), intent(in) :: a(3), b(3), e1, e2, k(3), tol
    integer, intent(in) :: n1, l1, m1, n2, l2, m2
    !
    real(wp) :: wig(3,7), xik(3,15), wik(3,15), a2(3), b1(3)
    complex(wp) :: func(15,15,15)
    complex(wp) :: integ_g, integ, integ_k(3), res_k(2)
    real(wp) :: kr, atol, r(3)
    integer :: i, i1, i2, i3, i_d
    !
    do i = 1, 3
      wig(i,:) = (b(i)-a(i))*wg/2._wp
      xik(i,:) = ( (a(i)-b(i))*xk +b(i)+a(i) )/2._wp
      wik(i,:) = (b(i)-a(i))*wk/2._wp
    end do
    func = 0._wp
    do i1=1,15; do i2=1,15; do i3=1,15
      call spher2cartez(xik(1,i1), xik(2,i2), xik(3,i3), r)
      kr = dot_product(k, r)
      func(i1,i2,i3) = norm_fac(e1, n1)* norm_fac(e2, n2) * &
        cmplx(cos(kr), sin(kr), wp) * exp(-(e1+e2)*xik(1,i1)) * &
        conjg( spherical_harmonic(l1, m1, xik(2,i2), xik(3,i3)) ) * &
        spherical_harmonic(l2, m2, xik(2,i2), xik(3,i3)) * &
        xik(1,i1)**(n1+n2) * sin(xik(2,i2))
    end do; end do; end do
    integ_g = 0._wp
    integ_k = 0._wp
    do i2=2,15,2; do i3=2,15,2
      do i1=2,15,2
        integ_g = integ_g + func(i1,i2,i3)*wig(1,i1/2)*wig(2,i2/2)*wig(3,i3/2)
      end do
      do i1=1,15
        integ_k(1) = integ_k(1) + func(i1,i2,i3)*wik(1,i1)*wig(2,i2/2)*wig(3,i3/2)
        integ_k(2) = integ_k(2) + func(i2,i1,i3)*wig(1,i2/2)*wik(2,i1)*wig(3,i3/2)
        integ_k(3) = integ_k(3) + func(i3,i2,i1)*wig(1,i3/2)*wig(2,i2/2)*wik(3,i1)
      end do
    end do; end do
    integ = 0._wp
    do i1=1,15; do i2=1,15; do i3=1,15
      integ = integ + func(i1,i2,i3)*wik(1,i1)*wik(2,i2)*wik(3,i3)
    end do; end do; end do
    i_d = maxloc(abs(integ_g-integ_k),1)
    if( abs(integ_g)<=epsilon(1._wp) .AND. abs(integ_k(i_d))<=epsilon(1._wp) ) then
      atol = 0._wp
    elseif( abs(integ_g) <= epsilon(1._wp) ) then
      atol = abs(integ_g-integ_k(i_d))
    else
      atol = abs(integ_g-integ_k(i_d))/abs(integ_g)
    end if
    if(atol<tol) then
      res = integ_k(i_d)
    else
      b1 = b
      b1(i_d) = (a(i_d)+b(i_d))/2
      res_k(1) = integrate(a, b1, n1, l1, m1, e1, n2, l2, m2, e2, k, tol)
      a2 = a
      a2(i_d) = (a(i_d)+b(i_d))/2
      res_k(2) = integrate(a2, b, n1, l1, m1, e1, n2, l2, m2, e2, k, tol)
      res = res_k(1) + res_k(2)
    endif
  end function integrate
end program
