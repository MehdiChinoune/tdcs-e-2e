program test_orbital_wave_functions
  use constants, only: wp, fac
  use types, only: orbit, atom
  use input, only: read_orbit
  use utils, only: norm_fac
  implicit none
  !
  character(len=2), dimension(7) :: atom_names = ["H ", "He", "Li", "Be", "B ", "C ", "Ne"]
  character(len=2) :: orbit_name
  type(atom) :: atom_target
  type(orbit) :: tmp_orbit
  integer :: i, i1, i2, j1, j2, ia, io, inp
  real(wp) :: integ
  !
  do i = 1, 7
    open( newunit = inp, file = 'data/'//trim(atom_names(i))//'.dat', status = 'old', action = 'read')
    ia = 0
    allocate(atom_target%orbits(0))
    do
      read(inp, fmt=*, iostat = io ) orbit_name
      if(io<0) exit
      ia = ia + 1
      call read_orbit(trim(atom_names(i))//'_'//orbit_name, tmp_orbit)
      atom_target%orbits = [ atom_target%orbits, tmp_orbit]
    end do
    !
    do i1=1, size(atom_target%orbits)
      do i2 = 1, i1
        associate (at1 => atom_target%orbits(i1), at2 => atom_target%orbits(i2))
          integ = 0._wp
          if (at1%l == at2%l) then
            integ = 0._wp
            do j1 = 1, at1%nf; do j2 = 1, at2%nf
              integ = integ + at1%a(j1)*at2%a(j2)*norm_fac(at1%e(j1),at1%n(j1))*&
                norm_fac(at2%e(j2),at2%n(j2))*&
                real(fac(at1%n(j1)+at2%n(j2)),wp)/(at1%e(j1)+at2%e(j2))**(at1%n(j1)+at2%n(j2)+1)
            end do; end do
            if ( i1==i2 .AND. abs(integ-1._wp) >= 1.e-5_wp) then
              print*, "orbit number: ", i1, "from ", atom_names(i), " is not normalized"
              print*, at1%l, integ
            elseif( i1/=i2 .AND. abs(integ) >= 1.e-5_wp) then
              print*, "orbitals number: ", i1, " and ", i2, "from ", atom_names(i), " are not orthogonal"
              print*, at1%l, integ
            endif
          endif
        end associate
      end do
    end do
    !
    deallocate(atom_target%orbits)
    !
  end do
  !
end program test_orbital_wave_functions
