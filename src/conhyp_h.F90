module conhyp_m
  use constants, only: RP
  implicit none

  interface

    elemental complex(RP) module function conhyp(A,B,Z,LNCHF,IP)
      complex(RP), intent(in) :: A, B, Z
      integer, intent(in) :: Lnchf, IP
    end function conhyp

  end interface

end module conhyp_m
