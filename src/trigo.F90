!> Trigonometric Transformations

!> This module contains some trigonometric transformations
submodule(trigo) trigo
  implicit none

contains

  pure module subroutine spher2cartez( km, theta, phi, k )
    real(RP), intent(in)  :: km, phi, theta
    real(RP), intent(out) :: k(3)
    k = km*[ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) ]
  end subroutine spher2cartez

  pure module subroutine cartez2spher( k, km, theta, phi )
    real(RP), intent(in)  :: k(3)
    real(RP), intent(out) :: km, theta, phi
    km = norm2(k)
    theta = acos( k(3)/km )
    phi = atan2( k(2), k(1) )
  end subroutine cartez2spher

end submodule trigo
