program main
  use constants        ,only: RP
  use fdcs_e2e         ,only: fdcs_fba_dw, fdcs_fba_cw, fdcs_fba_pw, fdcs_dwb, &
                                fdcs_bbk
  !
  implicit none
  integer(kind = selected_int_kind(6)) :: start, finish
  real(RP) :: rate
  integer :: in_unit, out_unit, narg
  character(len=4) :: arg1
  logical :: input_found
  !
  !
  inquire( file = "input.dat", exist = input_found )
  if( .not. input_found ) error stop "The input file 'input.dat' wasn't found, write your own"
  open( newunit = in_unit, file = 'input.dat', status = 'old', action = 'read' )
  open( newunit = out_unit, file = 'output.dat', status = 'replace', action = 'write' )
  !
  call system_clock(start,rate)
  !
  narg = command_argument_count()
  if(narg>0) then
    !
    call get_command_argument(1, arg1 )
    !
    if(trim(arg1)=='dw') then
      call fdcs_fba_dw(in_unit,out_unit)
    elseif(trim(arg1)=='cw') then
      call fdcs_fba_cw(in_unit,out_unit)
    elseif(trim(arg1)=='pw') then
      call fdcs_fba_pw(in_unit,out_unit)
    elseif(trim(arg1)=='dwb') then
      call fdcs_dwb(in_unit,out_unit)
    elseif(trim(arg1)=='bbk') then
      call fdcs_bbk(in_unit,out_unit)
    end if
    !
  else
    call fdcs_fba_pw(in_unit,out_unit)
  endif

  call system_clock(finish)

  close(in_unit)
  close(out_unit)

  print*, (finish-start)/rate

end program main
