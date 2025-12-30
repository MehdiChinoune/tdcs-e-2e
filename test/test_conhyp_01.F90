program test_conhyp_01
  use constants, only : wp, dp
  use arb_functions, only : zhypgeom_1f1
  use conhyp_m, only : conhyp
  implicit none
  !
  complex(wp) :: a, b, z, zh, ch
  integer :: i, ia, ib, iz, ja, jb, jz
  integer, dimension(9) :: i_tst = [ 0, 1, 2, 5, 10, 20, 50, 100, 200 ]
  ! Test 1F1(a,b,0)=1 when a is real and b is pure complex
  z = (0._wp,0._wp)
  do i = 1, 200
    a = cmplx( i*1._wp, 0._wp, wp )
    b = cmplx( 0._wp, i*1._wp, wp )
    if( conhyp( a, b, z, 0, 10 )/=(1._wp,0._wp) ) error stop " Failed : 1_F_1(a,b,z=0)/=1"
  end do
  ! Test against flint/arb hypgeom_1f1
  !$OMP PARALLEL DO COLLAPSE(6) PRIVATE(a, b, z, zh ,ch)
  do ia = 1, 7
    do ja = 1, 7
      do ib = 2, 9
        do jb = 2, 9
          do iz = 1, 8
            do jz = 1, 8
              !if(ib==1 .AND. jb==1) cycle ! b/=0
              a = cmplx( i_tst(ia)*1._dp, i_tst(ja)*1._dp, dp)
              b = cmplx( i_tst(ib)*1._dp, i_tst(jb)*1._dp, dp)
              z = cmplx( i_tst(iz)*1._dp, i_tst(jz)*1._dp, dp)
              zh = zhypgeom_1f1(a, b, z)
              ch = conhyp(a, b, z, 0, 10)
              if( abs(ch-zh)/abs(zh)>1.e-7_dp ) then
                print*, "a= ", a
                print*, "b= ", b
                print*, "z= ", z
                print*, "zhypgeom_1f1=", zh
                print*, "conhyp=", ch
                error stop " Failed : conhyp!=zhypgeom_1f1"
              endif
            end do
          end do
        end do
      end do
    end do
  end do
end program test_conhyp_01
