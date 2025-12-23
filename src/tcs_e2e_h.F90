module tcs_e2e
  use constants, only : wp
  implicit none

  logical :: show_output = .true.

  interface

    module subroutine tcs_fba_pw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine tcs_fba_pw

    module subroutine tcs_fba_cw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine tcs_fba_cw

    module subroutine tcs_fba_ocw(in_unit,out_unit)
      integer, intent(in) :: in_unit
      integer, intent(in) :: out_unit
    end subroutine tcs_fba_ocw

    module complex(wp) function eikr_element(n1, l1, m1, e1, n2, l2, m2, e2, k)
      integer, intent(in) :: n1, l1, m1, n2, l2, m2
      real(wp), intent(in) :: e1, e2, k(3)
    end function eikr_element

  end interface

end module tcs_e2e
