module arb_types
  use ISO_C_BINDING, only: c_double
  implicit none
  type, bind(c) :: complex_double
    real(c_double) :: real, imag
  end type
end module

module arb_functions
  use constants, only: dp
  use ISO_C_BINDING, only: c_int, c_double
  use arb_types, only: complex_double
  implicit none
  !int arb_fpwrap_double_hypgeom_pfq(double *res, const double *a, slong p, const double *b, &
  ! slong q, double z, int regularized, int flags)
  !int arb_fpwrap_cdouble_hypgeom_pfq(complex_double *res, const complex_double *a, slong p, &
  ! const complex_double *b, slong q, complex_double z, int regularized, int flags)

  interface
    integer(c_int) function arb_fpwrap_double_hypgeom_1f1(res, a, b, x, regularized, flags) bind(c)
      use ISO_C_BINDING
      integer(c_int), intent(in), value :: regularized, flags
      real(c_double), intent(in), value :: a, b, x
      real(c_double), intent(out) :: res
    end function
    integer(c_int) function arb_fpwrap_cdouble_hypgeom_1f1(res, a, b, x, regularized, flags) bind(c)
      use ISO_C_BINDING
      use arb_types
      integer(c_int), intent(in), value :: regularized, flags
      type(complex_double), intent(in), value :: a, b, x
      type(complex_double), intent(out) :: res
    end function
    integer(c_int) function arb_fpwrap_double_hypgeom_2f1(res, a, b, c, x, regularized, flags) bind(c)
      use ISO_C_BINDING
      integer(c_int), intent(in), value :: regularized, flags
      real(c_double), intent(in), value :: a, b, c, x
      real(c_double), intent(out) :: res
    end function
    integer(c_int) function arb_fpwrap_cdouble_hypgeom_2f1(res, a, b, c, x, regularized, flags) bind(c)
      use ISO_C_BINDING
      use arb_types
      integer(c_int), intent(in), value :: regularized, flags
      type(complex_double), intent(in), value :: a, b, c, x
      type(complex_double), intent(out) :: res
    end function
  end interface
contains
  real(dp) elemental impure function dhypgeom_1f1(a, b, x)
    real(dp), intent(in) :: a, b , x
    integer(c_int) :: arb_status
    arb_status = arb_fpwrap_double_hypgeom_1f1(dhypgeom_1f1, a, b, x, 0, 0)
  end function
  complex(dp) elemental impure function zhypgeom_1f1(a, b, z)
    complex(dp), intent(in) :: a, b , z
    type(complex_double) :: a_arb, b_arb, z_arb, zhypgeom_1f1_arb
    integer(c_int) :: arb_status
    a_arb = complex_double(a%re, a%im)
    b_arb = complex_double(b%re, b%im)
    z_arb = complex_double(z%re, z%im)
    arb_status = arb_fpwrap_cdouble_hypgeom_1f1(zhypgeom_1f1_arb, a_arb, b_arb, z_arb, 0, 0)
    zhypgeom_1f1 = cmplx(zhypgeom_1f1_arb%real, zhypgeom_1f1_arb%imag, dp)
  end function
end module
