MODULE fdcs_e2e
  USE constants ,ONLY : RP

  INTERFACE

    MODULE SUBROUTINE fdcs_fba_pw(in_unit,out_unit)
      INTEGER ,INTENT(IN) :: in_unit
      INTEGER ,INTENT(IN) :: out_unit
    END SUBROUTINE fdcs_fba_pw

    MODULE SUBROUTINE fdcs_fba_cw(in_unit,out_unit)
      INTEGER ,INTENT(IN) :: in_unit
      INTEGER ,INTENT(IN) :: out_unit
    END SUBROUTINE fdcs_fba_cw

    MODULE SUBROUTINE fdcs_fba_dw(in_unit,out_unit)
      INTEGER ,INTENT(IN) :: in_unit
      INTEGER ,INTENT(IN) :: out_unit
    END SUBROUTINE fdcs_fba_dw

  END INTERFACE

END MODULE fdcs_e2e
