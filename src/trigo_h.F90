!> Trigonometric Transformations

!> This module contains some trigonometric transformations
module trigo
  use constants, only : RP
  implicit none

  interface

    pure module subroutine spher2cartez( km, theta, phi, k )
      real(kind=RP), intent(in)  :: km, phi, theta
      real(kind=RP), intent(out) :: k(3)
    end subroutine spher2cartez

    pure module subroutine cartez2spher( k, km, theta, phi )
      real(kind=RP), intent(in)  :: k(3)
      real(kind=RP), intent(out) :: km, theta, phi
    end subroutine cartez2spher

  end interface

end module trigo
