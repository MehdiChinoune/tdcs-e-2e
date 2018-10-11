MODULE special_functions
  USE constants ,ONLY: RP
  IMPLICIT NONE
!  REAL(RP), PROTECTED :: fac(0:34), lnfac(0:400) !bug with Flang and Intel Compilers
  REAL(RP) :: fac(0:34), lnfac(0:400)
  LOGICAL, PROTECTED :: fac_called

  INTERFACE

    MODULE SUBROUTINE factorial()
    END SUBROUTINE

    ELEMENTAL MODULE FUNCTION cgamma(z, mo) RESULT(w)
      COMPLEX(RP) :: w
      COMPLEX(RP), INTENT(IN)  :: z
      INTEGER, INTENT(IN),OPTIONAL  :: mo
    END FUNCTION cgamma

    ELEMENTAL REAL(RP) MODULE FUNCTION assoc_legendre(l,m,x)
      INTEGER, INTENT(IN) :: l, m
      REAL(RP), INTENT(IN) :: x
    END FUNCTION assoc_legendre

    ELEMENTAL COMPLEX(RP) MODULE FUNCTION spherical_harmonic( l, m, theta, phi )
      INTEGER, INTENT(IN) :: l, m
      REAL(RP), INTENT(IN) :: theta, phi
    END FUNCTION spherical_harmonic

    PURE MODULE SUBROUTINE coul90(x, eta, lmin, lrange, fc, gc, fcp, gcp, kfn, ifail )
      INTEGER, INTENT(IN)  :: lmin, lrange, kfn
      INTEGER, INTENT(OUT), OPTIONAL :: ifail
      REAL(RP), INTENT(IN) :: x, eta
      REAL(RP), INTENT(OUT), DIMENSION(lmin:lmin+lrange) :: fc,  gc,  fcp, gcp
    END SUBROUTINE coul90

    ELEMENTAL MODULE SUBROUTINE  jwkb( x, eta, xl, fjwkb, gjwkb, iexp )
      REAL(RP), INTENT(IN)  :: x, eta, xl
      REAL(RP), INTENT(OUT) :: fjwkb, gjwkb
      INTEGER, INTENT(OUT) :: iexp
    END SUBROUTINE  jwkb

    PURE MODULE SUBROUTINE ricbes( x, lmax, psi, chi, psid, chid, ifail )
      REAL(RP), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: lmax
      INTEGER, INTENT(INOUT), OPTIONAL :: ifail
      REAL(RP), INTENT(OUT) :: psi(0:lmax), chi(0:lmax), psid(0:lmax), chid(0:lmax)
    END SUBROUTINE ricbes

    ELEMENTAL REAL(RP) MODULE FUNCTION symbol_3j(l1, l2, l3, m1, m2, m3)
      INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
    END FUNCTION symbol_3j

    ELEMENTAL COMPLEX(RP) MODULE FUNCTION conhyp_opt(a,z)
      REAL(RP), INTENT(IN) :: a, z
    END FUNCTION conhyp_opt

  END INTERFACE

END MODULE special_functions
