program test_conhyp_01
  use constants, only : wp
  use conhyp_m, only : conhyp
  implicit none
  !
  integer :: i
  complex(wp) :: a, b, z
  !
  z = (0._wp,0._wp)
  do i = 1, 200
    a = cmplx( i*1._wp, 0._wp, wp )
    b = cmplx( 0._wp, i*1._wp, wp )
    if( conhyp( a, b, z, 0, 10 )/=(1._wp,0._wp) ) error stop " Failed : 1_F_1(a,b,z=0)/=1"
  end do
  !
end program test_conhyp_01
