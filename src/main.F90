PROGRAM main
  USE constants        ,ONLY: RP
  USE fdcs_e2e         ,ONLY: fdcs_fba_dw, fdcs_fba_cw, fdcs_fba_pw, fdcs_dwb

  IMPLICIT NONE
  INTEGER(KIND = SELECTED_INT_KIND(6)) :: start, finish
  REAL(RP) :: rate
  INTEGER :: in_unit, out_unit, narg
  CHARACTER(LEN=4) :: arg1

  OPEN( newunit = in_unit, FILE = 'input.dat', STATUS = 'old', ACTION = 'read')
  OPEN( newunit = out_unit, FILE = 'output.dat', STATUS = 'replace', ACTION = 'write')

  CALL SYSTEM_CLOCK(start,rate)

  narg = COMMAND_ARGUMENT_COUNT()
  IF(narg>0) THEN

    CALL GET_COMMAND_ARGUMENT(1, arg1 )

    IF(TRIM(arg1)=='dw') THEN
      CALL fdcs_fba_dw(in_unit,out_unit)
    ELSEIF(TRIM(arg1)=='cw') THEN
      CALL fdcs_fba_cw(in_unit,out_unit)
    ELSEIF(TRIM(arg1)=='pw') THEN
      CALL fdcs_fba_pw(in_unit,out_unit)
    END IF
  ELSE
    CALL fdcs_fba_pw(in_unit,out_unit)
  ENDIF

  CALL SYSTEM_CLOCK(finish)

  CLOSE(in_unit)
  CLOSE(out_unit)

  PRINT*, (finish-start)/rate

END PROGRAM main
