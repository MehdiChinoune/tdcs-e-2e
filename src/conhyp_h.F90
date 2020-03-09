module conhyp_m
  use constants, only : wp
  implicit none

  interface

    elemental complex(wp) module function conhyp(A,B,Z,LNCHF,IP)
      complex(wp), intent(in) :: A, B, Z
      integer, intent(in) :: Lnchf, IP
    end function conhyp

  end interface

end module conhyp_m
