module special_functions
  use constants, only : wp
  implicit none
  real(wp), protected :: fac(0:34), lnfac(0:400)
  logical, protected :: fac_called = .FALSE.

  interface

    module subroutine factorial()
    end subroutine

    elemental module function cgamma(z, mo) result(w)
      complex(wp) :: w
      complex(wp), intent(in)  :: z
      integer, intent(in) :: mo
    end function cgamma

    elemental real(wp) module function assoc_legendre(l,m,x)
      integer, intent(in) :: l, m
      real(wp), intent(in) :: x
    end function assoc_legendre

    elemental complex(wp) module function spherical_harmonic( l, m, theta, phi )
      integer, intent(in) :: l, m
      real(wp), intent(in) :: theta, phi
    end function spherical_harmonic

    pure module subroutine coul90(x, eta, lmin, lrange, fc, gc, fcp, gcp, kfn, ifail )
      integer, intent(in)  :: lmin, lrange, kfn
      integer, intent(out), optional :: ifail
      real(wp), intent(in) :: x, eta
      real(wp), intent(out), dimension(lmin:lmin+lrange) :: fc,  gc,  fcp, gcp
    end subroutine coul90

    pure module subroutine ricbes( x, lmax, psi, chi, psid, chid, ifail )
      real(wp), intent(in) :: x
      integer, intent(in) :: lmax
      integer, intent(inout), optional :: ifail
      real(wp), intent(out) :: psi(0:lmax), chi(0:lmax), psid(0:lmax), chid(0:lmax)
    end subroutine ricbes

    elemental real(wp) module function symbol_3j(l1, l2, l3, m1, m2, m3)
      integer, intent(in) :: l1, l2, l3, m1, m2, m3
    end function symbol_3j

    elemental complex(wp) module function conhyp_opt(ai,zi)
      real(wp), intent(in) :: ai, zi
    end function conhyp_opt

  end interface

end module special_functions
