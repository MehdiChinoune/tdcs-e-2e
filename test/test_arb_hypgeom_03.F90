
program test_arb_hypgeom_03
  use constants, only : dp
  use arb_functions, only : zhypgeom_1f1
  implicit none
  !
  integer :: ia, ib, iz, ja, jb, jz
  complex(dp) :: a, b, z, zh1, zh2, zh3
  integer, dimension(8) :: i_tst = [ 1, 2, 5, 10, 20, 50, 100, 200 ]
  ! Test recurrence relation b(b-1)*1F1(a,b-1,z)+b(1-b-z)*1F1(a,b,z)+z(b-a)*1F1(a,b+1,z)=0
  do ia = 1, 8
    do ja = 1, 8
      do ib = 1, 8
        do jb = 1, 8
          do iz = 1, 8
            do jz = 1, 8
              a = cmplx( i_tst(ia)*1._dp, i_tst(ja)*1._dp, dp)
              b = cmplx( i_tst(ib)*1._dp, i_tst(jb)*1._dp, dp)
              z = cmplx( i_tst(iz)*1._dp, i_tst(jz)*1._dp, dp)
              zh1 = b*(b-1)*zhypgeom_1f1(a, b-1, z)
              zh2 = b*(1-b-z)*zhypgeom_1f1(a, b, z)
              zh3 = z*(b-a)*zhypgeom_1f1(a, b+1, z)
              if( abs(zh1+zh2+zh3)/max(abs(zh1),abs(zh2),abs(zh3))>1.e-7 ) then
                print*, "a=", a
                print*, "b=", b
                print*, "z=", z
                print*, "b(b-1)*1F1(a,b-1,z)=", zh1
                print*, "b(1-b-z)*1F1(a,b,z)=", zh2
                print*, "z(b-a)*1F1(a,b+1,z)=", zh3
                error stop " Failed : b(b-1)*1F1(a,b-1,z)+b(1-b-z)*1F1(a,b,z)+z(b-a)*1F1(a,b+1,z)/=0"
              endif
            end do
          end do
        end do
      end do
    end do
  end do
end program test_arb_hypgeom_03
