module utils
  use constants, only : wp
  implicit none

  interface

    elemental real(wp) module function norm_fac(e,n)
      real(wp), intent(in) :: e
      integer , intent(in) :: n
    end function

    elemental real(wp) module function y1y2y3(l1, l2, l3, m1, m2, m3 )
      integer, intent(in) :: l1, l2, l3, m1, m2, m3
    end function y1y2y3

    !> ode_second_dw
    !!
    !! This subroutine solve Equation of the form
    !! s_l''(r) +f_l(r)*s_l(r) = km**2*s_l(r)

    module subroutine ode_second_dw(km, lmax, rc, z, f, s, delta )
      integer, intent(in) :: lmax, z
      real(wp), intent(in)  :: f(0:,0:), rc, km
      real(wp), intent(out) :: s(0:,0:), delta(0:lmax)
    end subroutine ode_second_dw

    module subroutine calculate_U(Atom_name, Orbit_name, r, U, state )
      character(len=2), intent(in) :: Atom_name, Orbit_name
      real(wp)   , intent(in) :: r(:)
      real(wp)   , intent(out) :: U(:)
      integer, intent(in) :: state
    end subroutine calculate_U

    pure module subroutine INTRPL(X, Y, U, V )
      real(wp), intent(in) :: X(:), Y(:), U(:)
      real(wp), intent(out) :: V(:)
    end subroutine INTRPL

  end interface

end module utils
