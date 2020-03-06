module fdcs_e2e
  use constants ,only: RP
  implicit none

  logical :: show_output = .true.

  interface

    module subroutine fdcs_fba_pw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine fdcs_fba_pw

    module subroutine fdcs_fba_cw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine fdcs_fba_cw

    module subroutine fdcs_fba_dw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine fdcs_fba_dw

    module subroutine fdcs_dwb(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine fdcs_dwb

    module subroutine fdcs_bbk(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine fdcs_bbk


  end interface

end module fdcs_e2e
