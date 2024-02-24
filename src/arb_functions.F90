module arb_types
  use ISO_C_BINDING, only: c_double
  implicit none
  type, bind(c) :: complex_double
    real(c_double) :: real, imag
  end type
end module

module arb_functions
  implicit none
  !int arb_fpwrap_double_hypgeom_pfq(double *res, const double *a, slong p, const double *b, slong q, double z, int regularized, int flags)
  !int arb_fpwrap_cdouble_hypgeom_pfq(complex_double *res, const complex_double *a, slong p, const complex_double *b, slong q, complex_double z, int regularized, int flags)

  interface
    integer(c_int) function dhypgeom_1f1(res, a, b, x, regularized, flags) bind(c, name="arb_fpwrap_double_hypgeom_1f1")
      use ISO_C_BINDING
      integer(c_int), intent(in) :: regularized, flags
      real(c_double), intent(in), value :: a, b, x
      real(c_double), intent(out) :: res
    end function
    integer(c_int) function zhypgeom_1f1(res, a, b, x, regularized, flags) bind(c, name="arb_fpwrap_cdouble_hypgeom_1f1")
      use ISO_C_BINDING
      use arb_types
      integer(c_int), intent(in) :: regularized, flags
      type(complex_double), intent(in), value :: a, b, x
      type(complex_double), intent(out) :: res
    end function
    integer(c_int) function dhypgeom_2f1(res, a, b, c, x, regularized, flags) bind(c, name="arb_fpwrap_double_hypgeom_2f1")
      use ISO_C_BINDING
      integer(c_int), intent(in) :: regularized, flags
      real(c_double), intent(in), value :: a, b, c, x
      real(c_double), intent(out) :: res
    end function
    integer(c_int) function zhypgeom_2f1(res, a, b, c, x, regularized, flags) bind(c, name="arb_fpwrap_cdouble_hypgeom_2f1")
      use ISO_C_BINDING
      use arb_types
      integer(c_int), intent(in) :: regularized, flags
      type(complex_double), intent(in), value :: a, b, c, x
      type(complex_double), intent(out) :: res
    end function
  end interface
end module
