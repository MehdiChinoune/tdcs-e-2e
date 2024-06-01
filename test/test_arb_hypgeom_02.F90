
program test_arb_hypgeom_02
  use constants, only : dp
  use arb_functions, only : zhypgeom_1f1
  implicit none
  !
  integer :: ia, ib, iz, ja, jb, jz
  complex(dp) :: a, b, z, zh1, zh2, zh3
  integer, dimension(8) :: i_tst = [ 1, 2, 5, 10, 20, 50, 100, 200 ]
  ! Test recurrence relation (b-a)*1F1(a-1,b,z)+(2a-b+z)*1F1(a,b,z)-a*1F1(a+1,b,z)=0
  do ia = 1, 8
    do ja = 1, 8
      do ib = 1, 8
        do jb = 1, 8
          do iz = 1, 8
            do jz = 1, 8
              a = cmplx( i_tst(ia)*1._dp, i_tst(ja)*1._dp, dp)
              b = cmplx( i_tst(ib)*1._dp, i_tst(jb)*1._dp, dp)
              z = cmplx( i_tst(iz)*1._dp, i_tst(jz)*1._dp, dp)
              zh1 = (b-a)*zhypgeom_1f1(a-1, b, z)
              zh2 = (2*a-b+z)*zhypgeom_1f1(a, b, z)
              zh3 = -a*zhypgeom_1f1(a+1, b, z)
              if( abs(zh1+zh2+zh3)/max(abs(zh1),abs(zh2),abs(zh3))>1.e-7 ) then
                print*, "a=", a
                print*, "b=", b
                print*, "z=", z
                print*, "(b-a)*1F1(a-1,b,z)=", zh1
                print*, "(2a-b+z)*1F1(a,b,z)=", zh2
                print*, "-a*1F1(a+1,b,z)=", zh3
                error stop " Failed : (b-a)*1F1(a-1,b,z)+(2a-b+z)*1F1(a,b,z)-a*1F1(a+1,b,z)/=0"
              endif
            end do
          end do
        end do
      end do
    end do
  end do
end program test_arb_hypgeom_02
