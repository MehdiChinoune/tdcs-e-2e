module integration
  use constants ,only: RP
  implicit none

  interface

    module subroutine clenshaw_curtis( a, b, x, w, n )
      real(RP), intent(in) :: a, b
      integer, intent(in) :: n
      real(RP), intent(out), allocatable :: w(:), x(:)
    end subroutine clenshaw_curtis

    pure module subroutine gauleg(a,b,x,w,n)
      integer, intent(in) :: n
      real(RP), intent(in) :: a,b
      real(RP), intent(out), allocatable :: x(:),w(:)
    end subroutine gauleg

    pure module subroutine pd(sx,sw,n)
      integer, intent(in) :: n
      real(RP), intent(out), contiguous :: sw(:)
      real(RP), intent(inout), contiguous :: sx(:)
    end subroutine pd

  end interface

end module integration
