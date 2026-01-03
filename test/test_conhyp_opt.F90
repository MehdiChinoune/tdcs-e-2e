program test_conhyp_opt
  use constants, only : wp, dp
  use arb_functions, only : zhypgeom_1f1
  use special_functions, only : conhyp_opt
  implicit none
  !
  integer :: ia, iz, i_tst(8) = [ 1, 2, 5, 10, 20, 50, 100, 200 ]
  complex(wp) :: a, b, z, cho
  complex(dp) :: a_dp, b_dp, z_dp, zh
  ! Test My simplified function
  b = cmplx(1._wp, 0._wp, wp)
  b_dp = cmplx(1._dp, 0._dp, dp)
  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(a, z, a_dp, z_dp, zh, cho)
  do ia = 1, 4
    do iz = 1, 4
      a = cmplx( 0._wp, i_tst(ia)*1._wp, wp)
      z = cmplx( 0._wp, i_tst(iz)*1._wp, wp)
      cho = conhyp_opt(a%im, z%im)
      a_dp = cmplx( 0._dp, i_tst(ia)*1._dp, dp)
      z_dp = cmplx( 0._dp, i_tst(iz)*1._dp, dp)
      zh = zhypgeom_1f1(a_dp, b_dp, z_dp)
      if( abs(cho-zh)/abs(zh)>1.e-7_wp ) then
        print*, "a= i*", i_tst(ia)
        print*, "z= i*", i_tst(iz)
        print*, "zhypgeom_1f1=", zh
        print*, "conhyp_opt=", cho
        error stop " Failed : conhyp_opt!=zhypgeom_1f1"
      endif
    end do
  end do
end program test_conhyp_opt

