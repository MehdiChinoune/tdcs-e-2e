module integration
  use constants, only : wp
  implicit none

  interface

    module subroutine clenshaw_curtis( a, b, x, w, n )
      real(wp), intent(in) :: a, b
      integer, intent(in) :: n
      real(wp), intent(out), allocatable :: w(:), x(:)
    end subroutine clenshaw_curtis

    pure module subroutine gauleg(a,b,x,w,n)
      integer, intent(in) :: n
      real(wp), intent(in) :: a,b
      real(wp), intent(out), allocatable :: x(:),w(:)
    end subroutine gauleg

    pure module subroutine pd(sx,sw,n)
      integer, intent(in) :: n
      real(wp), intent(out), contiguous :: sw(:)
      real(wp), intent(inout), contiguous :: sx(:)
    end subroutine pd

  end interface

end module integration
