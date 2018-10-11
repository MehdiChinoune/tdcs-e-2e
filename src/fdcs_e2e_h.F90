MODULE fdcs_e2e
  USE constants ,ONLY: RP
  IMPLICIT NONE

  INTERFACE

  MODULE SUBROUTINE fdcs_fba_pw(in_unit,out_unit)
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit
  END SUBROUTINE fdcs_fba_pw

  MODULE SUBROUTINE fdcs_fba_cw(in_unit,out_unit)
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit
  END SUBROUTINE fdcs_fba_cw

  MODULE SUBROUTINE fdcs_fba_dw(in_unit,out_unit)
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit
  END SUBROUTINE fdcs_fba_dw

  MODULE SUBROUTINE fdcs_dwb(in_unit,out_unit)
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit
  END SUBROUTINE fdcs_dwb

  MODULE SUBROUTINE fdcs_bbk(in_unit,out_unit)
    INTEGER, INTENT(IN) :: in_unit
    INTEGER, INTENT(IN) :: out_unit
  END SUBROUTINE fdcs_bbk

  COMPLEX(RP) MODULE FUNCTION U_bbk(alpha1, alpha2, alpha3, k1, k2, k3, lam1, lam2, lam3 &
    , p1, p2)
    REAL(RP), INTENT(IN) :: alpha1, alpha2, alpha3
    REAL(RP), INTENT(IN) :: k1(3), k2(3), k3(3)
    REAL(RP), INTENT(IN) :: lam1, lam2, lam3
    REAL(RP), INTENT(IN) :: p1(3), p2(3)
  END FUNCTION

  MODULE SUBROUTINE calculate_chi( km, r, U_tmp, z, x, chi, delta )
    REAL(RP), INTENT(IN) :: U_tmp(0:), r(0:), x(:), km
    INTEGER, INTENT(IN) :: z
    REAL(RP), INTENT(OUT) :: chi(:,0:), delta(0:)
  END SUBROUTINE calculate_chi

  MODULE SUBROUTINE dwb_integrals( chi_0, chi_a, chi_b, sig_0, sig_a, sig_b, wf, x, w, lo &
    , integral )
    REAL(RP), INTENT(IN) :: chi_0(:,0:), chi_a(:,0:), chi_b(:,0:), wf(:)
    REAL(RP), INTENT(IN) :: x(:), w(:)
    REAL(RP), INTENT(IN) :: sig_0(0:), sig_a(0:), sig_b(0:)
    INTEGER, INTENT(IN) :: lo
    COMPLEX(RP), ALLOCATABLE, INTENT(OUT) :: integral(:,:,:,:)
  END SUBROUTINE dwb_integrals

  PURE COMPLEX(RP) MODULE FUNCTION tpw( n, l, m, e, ke, k)
    INTEGER, INTENT(IN) :: n,l,m
    REAL(RP), INTENT(IN) :: e,ke(3)
    REAL(RP), INTENT(IN),OPTIONAL :: k(3)
  END FUNCTION tpw

  PURE COMPLEX(RP) MODULE FUNCTION tcw( n, l, m, e, alpha, ke, k)
    INTEGER, INTENT(IN) :: n, l, m
    REAL(RP), INTENT(IN) :: e, alpha, ke(3), k(3)
  END FUNCTION tcw

  PURE COMPLEX(RP) MODULE FUNCTION tcw0( n, l, m, e, alpha, ke)
    INTEGER, INTENT(IN) :: n, l, m
    REAL(RP), INTENT(IN) :: e, alpha, ke(3)
  END FUNCTION tcw0

  ELEMENTAL COMPLEX(RP) MODULE FUNCTION powcc(z1, y2)
    COMPLEX(RP), INTENT(IN) :: z1
    REAL(RP), INTENT(IN) :: y2
  END FUNCTION powcc

  END INTERFACE

END MODULE fdcs_e2e
