module types
  use constants, only: wp
  implicit none
  type :: orbit
    integer :: nf ! number of wave-functions
    integer :: nelec ! number of electrons in the orbit
    integer :: l
    integer, allocatable :: n(:)
    real(wp), allocatable :: a(:), e(:)
    real(wp) :: Ie ! Ionization Energy
  end type
  type :: atom
    type(orbit), allocatable :: orbits(:)
  end type
end module
