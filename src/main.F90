program main
  use fdcs_e2e        ,only: fdcs_fba_dw, fdcs_fba_cw, fdcs_fba_pw, fdcs_dwb, &
                             fdcs_bbk
  use tcs_e2e         ,only: tcs_fba_pw, tcs_fba_cw, tcs_fba_ocw
  use constants       ,only: wp
  !
  implicit none
  integer(kind = selected_int_kind(6)) :: start, finish
  real(wp) :: rate
  integer :: in_unit, narg
  character(len=4) :: arg1, arg2
  logical :: input_found
  !
  inquire( file = "input.dat", exist = input_found )
  if( .not. input_found ) error stop "The input file 'input.dat' wasn't found, write your own"
  open( newunit = in_unit, file = 'input.dat', status = 'old', action = 'read' )
  !
  narg = command_argument_count()
  if(narg>0) then
    call get_command_argument(1, arg1 )
    if(narg>1) then
      call get_command_argument(2, arg2 )
    endif
  else
    print*, "Which Cross section do you want to calculate?"
    read*, arg1
    print*, "Which Method do you want to use?"
    read*, arg2
  endif
  !
  call system_clock(start,rate)
  !
  if (trim(arg1)=='fdcs') then
    if(trim(arg2)=='pw') then
      call fdcs_fba_pw(in_unit)
    elseif(trim(arg2)=='cw') then
      call fdcs_fba_cw(in_unit)
    elseif(trim(arg2)=='dw') then
      call fdcs_fba_dw(in_unit)
    elseif(trim(arg2)=='dwb') then
      call fdcs_dwb(in_unit)
    elseif(trim(arg2)=='bbk') then
      call fdcs_bbk(in_unit)
    end if
  elseif (trim(arg1)=='tcs') then
    if(trim(arg2)=='pw') then
      call tcs_fba_pw(in_unit)
    elseif(trim(arg2)=='cw') then
      call tcs_fba_cw(in_unit)
    elseif(trim(arg2)=='ocw') then
      call tcs_fba_ocw(in_unit)
    endif
  end if
  !
  close(in_unit)
  !
  call system_clock(finish)
  print*, (finish-start)/rate
  !
end program main
