MODULE special_functions
  USE constants, ONLY : RP
  IMPLICIT NONE

  INTERFACE

    ELEMENTAL MODULE FUNCTION cgamma( z, mo ) RESULT(w)
      COMPLEX(KIND=RP) :: w
      COMPLEX(KIND=RP), INTENT(IN)  :: z
      INTEGER,OPTIONAL, INTENT(IN)  :: MO
    END FUNCTION cgamma

    ELEMENTAL REAL(KIND=RP) MODULE FUNCTION assoc_legendre( l, m, x )
      INTEGER      , INTENT(IN) :: l, m
      REAL(KIND=RP), INTENT(IN) :: x
    END FUNCTION assoc_legendre

    ELEMENTAL COMPLEX(KIND=RP) MODULE FUNCTION spherical_harmonic( l, m, theta, phi )
      INTEGER      , INTENT(IN) :: l, m
      REAL(KIND=RP), INTENT(IN) :: theta, phi
    END FUNCTION spherical_harmonic

    PURE MODULE SUBROUTINE coul90(x, eta, lmin, lrange, fc, gc, fcp, gcp, kfn, ifail )
      INTEGER      , INTENT(IN)  :: lmin,lrange, kfn
      INTEGER      , INTENT(OUT), OPTIONAL :: ifail
      REAL(KIND=RP), INTENT(IN)  :: x,eta
      REAL(KIND=RP), INTENT(OUT),DIMENSION(lmin:lmin+lrange) :: fc,  gc,  fcp, gcp
    END SUBROUTINE coul90

    PURE MODULE SUBROUTINE ricbes( x, lmax, psi, chi, psid, chid, ifail )
      REAL(KIND=RP), INTENT(IN)    :: x
      INTEGER      , INTENT(IN)    :: lmax
      INTEGER      , INTENT(INOUT) :: ifail
      REAL(KIND=RP), INTENT(OUT), DIMENSION(0:lmax)   :: psi, chi, psid, chid
    END SUBROUTINE ricbes

    ELEMENTAL REAL(RP) MODULE FUNCTION symbol_3j(l1, l2, l3, m1, m2, m3)
      INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
    END FUNCTION symbol_3j

  END INTERFACE

END MODULE special_functions
