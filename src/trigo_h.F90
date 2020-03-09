!> Trigonometric Transformations

!> This module contains some trigonometric transformations
module trigo
  use constants, only : wp
  implicit none

  interface

    pure module subroutine spher2cartez( km, theta, phi, k )
      real(wp), intent(in)  :: km, phi, theta
      real(wp), intent(out) :: k(3)
    end subroutine spher2cartez

    pure module subroutine cartez2spher( k, km, theta, phi )
      real(wp), intent(in)  :: k(3)
      real(wp), intent(out) :: km, theta, phi
    end subroutine cartez2spher

  end interface

end module trigo
