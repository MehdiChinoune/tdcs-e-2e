PROGRAM main
  USE constants        , ONLY : RP
  USE fdcs_e2e         , ONLY : fdcs_fba_dw, fdcs_fba_cw, fdcs_fba_pw, fdcs_dwb

  IMPLICIT NONE
  integer(KIND=selected_int_kind(6)) :: start, finish
  REAL(KIND=RP) :: rate
  integer :: in_unit, out_unit, narg
  character(len=4) :: arg1

  open( newunit=in_unit, file='input.dat', status='old', action='read')
  open( newunit=out_unit, file='output.dat', status='replace', action='write')

  call SYSTEM_CLOCK(start,rate)

  narg = COMMAND_ARGUMENT_COUNT()
  if(narg>0) then

    call get_command_argument(1, arg1 )

    if(trim(arg1)=='dw') then
      call fdcs_fba_dw(in_unit,out_unit)
    elseif(trim(arg1)=='cw') then
      call fdcs_fba_cw(in_unit,out_unit)
    elseif(trim(arg1)=='pw') then
      call fdcs_fba_pw(in_unit,out_unit)
    end if
  else
    block

      !call fdcs_fba_cw(in_unit,out_unit)
      !rewind in_unit
      !call fdcs_fba_dw(in_unit,out_unit)
      !rewind in_unit
      call fdcs_dwb(in_unit,out_unit)

    end block
  endif

  call SYSTEM_CLOCK(finish)

  close(in_unit)
  close(out_unit)

  print*,(finish-start)/rate

END PROGRAM main
