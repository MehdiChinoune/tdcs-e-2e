module functions
  use constants, only: wp
  implicit none
  real(wp), parameter :: alphax = 5._wp, alphay = 5._wp
contains
  real(wp) function gauss2d(x,y)
    real(wp), intent(in) :: x, y
    gauss2d = exp(-alphax*x**2)*exp(-alphay*y**2)
  end function
end module

program test_adaptive_integration
  use functions, only: gauss2d
  use integration, only: gauleg
  use constants, only: wp, pi, xg, wg, xk, wk
  implicit none
  interface
    real(wp) function funct2d(x,y)
      import wp
      real(wp), intent(in) :: x, y
    end function
  end interface
  real(wp) :: integ, alpha(2), a(2), b(2)
  !
  alpha = [5._wp, 5._wp]
  a = [-10._wp, -10._wp]
  b = [10._wp, 10._wp]
  integ=integrate(gauss2d, a, b, 0.1_wp)
  print*, integ
  print*, (integ-pi/sqrt(alpha(1)*alpha(2)))
  !
contains
  recursive real(wp) function integrate(funct, a, b, tol) result(res)
    real(wp), intent(in) :: a(:), b(:), tol
    procedure(funct2d) :: funct
    !
    real(wp) :: wig(size(a),7), xik(size(a),15), wik(size(a),15), func(15,15)
    real(wp) :: integ_g, integ_k(size(a)), res_k(size(a)), a2(size(a)), b1(size(a))
    real(wp) :: atol
    integer :: i, i1, i2, i_d
    !
    do i = 1, size(a)
!      xig(i,:) = ( (a(i)-b(i))*xg +b(i)+a(i) )/2._wp
      wig(i,:) = (b(i)-a(i))*wg/2._wp
      xik(i,:) = ( (a(i)-b(i))*xk +b(i)+a(i) )/2._wp
      wik(i,:) = (b(i)-a(i))*wk/2._wp
    end do
    do i1=1,15; do i2=1,15
      func(i1,i2)=funct(xik(1,i1),xik(2,i2))
    end do; end do
    integ_g = 0._wp
    do i1=2,15,2; do i2=2,15,2
      integ_g = integ_g + func(i1,i2)*wig(1,i1/2)*wig(2,i2/2)
    end do; end do
    integ_k = 0._wp
    do i1=1,15; do i2=2,15,2
      integ_k(1) = integ_k(1) + func(i1,i2)*wik(1,i1)*wig(2,i2/2)
    end do; end do
    do i1=2,15,2; do i2=1,15
      integ_k(2) = integ_k(2) + func(i1,i2)*wig(1,i1/2)*wik(2,i2)
    end do; end do
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
      res_k(1) = integrate(funct, a, b1, tol)
      a2 = a
      a2(i_d) = (a(i_d)+b(i_d))/2
      res_k(2) = integrate(funct, a2, b, tol)
      res = res_k(1) + res_k(2)
    endif
  end function integrate
end program test_adaptive_integration
