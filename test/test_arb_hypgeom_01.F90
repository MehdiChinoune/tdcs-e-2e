
program test_arb_hypgeom_01
  use constants, only : dp
  use arb_functions, only : zhypgeom_1f1
  implicit none
  !
  integer :: i, ia, ib, iz, ja, jb, jz
  complex(dp) :: a, b, z, zh1, zh2
  integer, dimension(8) :: i_tst = [ 1, 2, 5, 10, 20, 50, 100, 200 ]
  ! Test 1F1(a,b,0)=1
  z = cmplx(0._dp,0._dp, dp)
  do i = 1, 200
    a = cmplx( i*1._dp, i*1._dp, dp )
    b = cmplx( i*1._dp, i*1._dp, dp )
    if( zhypgeom_1f1(a, b, z) /= (1._dp,0._dp) ) then
      print*, i
      print*, zhypgeom_1f1(a, b, z)
      error stop " Failed : 1_F_1(a,b,z=0)/=1"
    endif
  end do
  ! Test transformation 1F1(a,b,z)=exp(z)*1F1(b-a,b,-z)
  do ia = 1, 8
    do ja = 1, 8
      do ib = 1, 8
        do jb = 1, 8
          do iz = 1, 8
            do jz = 1, 8
              a = cmplx( i_tst(ia)*1._dp, i_tst(ja)*1._dp, dp)
              b = cmplx( i_tst(ib)*1._dp, i_tst(jb)*1._dp, dp)
              z = cmplx( i_tst(iz)*1._dp, i_tst(jz)*1._dp, dp)
              zh1 = zhypgeom_1f1(a, b, z)
              zh2 = exp(z)*zhypgeom_1f1(b-a, b, -z)
              if( abs(zh1-zh2)/abs(zh1)>1.e-7 ) then
                print*, "a=", a
                print*, "b=", b
                print*, "z=", z
                print*, "1F1(a,b,z)=", zh1
                print*, "exp(z)*1F1(b-a,b,z)=", zh2
                error stop " Failed : 1F1(a,b,z)/=exp(z)*1F1(b-a,b,z)="
              endif
            end do
          end do
        end do
      end do
    end do
  end do
end program test_arb_hypgeom_01
