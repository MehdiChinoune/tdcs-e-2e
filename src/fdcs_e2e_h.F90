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

  END INTERFACE

END MODULE fdcs_e2e