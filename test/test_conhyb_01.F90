program test_conhyp_01
  use constants, only : RP
  use conhyp_m, only : conhyp
  implicit none
  !
  integer :: i
  complex(RP) :: a, b, z
  !
  z = (0._RP,0._RP)
  do i = 1, 200
    a = cmplx( i*1._RP, 0._RP, RP )
    b = cmplx( 0._RP, i*1._RP, RP )
    if( conhyp( a, b, z, 0, 10 )/=(1._RP,0._RP) ) error stop " Failed : 1_F_1(a,b,z=0)/=1"
  end do
  !
end program test_conhyp_01
