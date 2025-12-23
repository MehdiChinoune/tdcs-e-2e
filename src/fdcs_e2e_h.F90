module fdcs_e2e
  use constants, only : wp
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

    module pure complex(wp) function tpw( n, l, m, e, ke, k)
      integer, intent(in) :: n,l,m
      real(wp), intent(in) :: e,ke(3)
      real(wp), intent(in),optional :: k(3)
    end function tpw

    module pure complex(wp) function tcw( n, l, m, e, alpha, ke, k)
      integer, intent(in) :: n, l, m
      real(wp), intent(in) :: e, alpha, ke(3), k(3)
    end function tcw

    module pure complex(wp) function tcw0( n, l, m, e, alpha, ke)
      integer, intent(in) :: n, l, m
      real(wp), intent(in) :: e, alpha, ke(3)
    end function tcw0


  end interface

end module fdcs_e2e
