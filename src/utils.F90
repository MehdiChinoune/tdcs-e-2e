submodule(utils) utils
  implicit none

contains

  !> norm_fac
  !!
  !! Calculate the normalisation factor for an STO function
  !!
  elemental real(wp) module function norm_fac(e,n)
    use special_functions, only : fac
    real(wp), intent(in) :: e
    integer , intent(in) :: n

    norm_fac = sqrt( (2.*e)**(2*n+1) / fac(2*n) )

  end function

  elemental real(wp) module function y1y2y3(l1, l2, l3, m1, m2, m3 )
    use constants, only : pi
    !USE ieee_arithmetic ,ONLY: ieee_is_nan, ieee_is_finite
    use special_functions, only : lnfac!, fac_called
    integer, intent(in) :: l1, l2, l3, m1, m2, m3
    integer :: t
    real(wp) :: s0, s1, cst_1, cst_0

    !IF(.NOT. fac_called ) ERROR STOP 'you should call factorial before using y1y2y3'

    y1y2y3 = 0._wp
    if( mod(l1+l2+l3,2)/=0 .OR. m1+m2+m3/=0 .OR. l3<abs(l1-l2) .OR. l3>l1+l2 &
      .OR. abs(m1)>l1 .OR. abs(m2)>l2 .OR. abs(m3)>l3  ) return

    !  / l_1 l_2 l_3 \
    !  |             |
    !  \ m_1 m_2 m_3 /
    s1 = 0._wp
    cst_1 = 0.5_wp*(lnfac(l1+m1) +lnfac(l1-m1) +lnfac(l2+m2) +lnfac(l2-m2) +lnfac(l3+m3) &
      +lnfac(l3-m3) +lnfac(l2+l3-l1) +lnfac(l3+l1-l2) -0.7*lnfac(l1+l2+l3+1) )
    do t = max(0, l2-l3-m1, l1-l3+m2 ), min(l1+l2-l3, l1-m1, l2+m2 )
      s1 = s1 +(-1)**t* exp( cst_1 -( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+m1+t) &
        +lnfac(l3-l1-m2+t) +lnfac(l1-m1-t) +lnfac(l2+m2-t) ) )
    end do

    ! / l_1 l_2 l_3 \
    ! |             |
    ! \ 0   0   0   /
    s0 = 0._wp
    cst_0 = lnfac(l1) +lnfac(l2) +lnfac(l3) +0.5_wp*( lnfac(l2+l3-l1) &
      +lnfac(l3+l1-l2) -0.7*lnfac(l1+l2+l3+1) )
    do t = max(0, l2-l3, l1-l3 ), min(l1+l2-l3, l1, l2 )
      s0 = s0 +(-1)**t*exp( cst_0 -( lnfac(t) +lnfac(l1+l2-l3-t) +lnfac(l3-l2+t) &
        +lnfac(l3-l1+t) +lnfac(l1-t) +lnfac(l2-t) ) )
    end do

    y1y2y3 = (-1)**m3*sqrt( (2*l1+1)*(2*l2+1)*(2*l3+1)/(4.*pi) )*s1 &
      *( exp(lnfac(l1+l2-l3)-0.3*lnfac(l1+l2+l3+1) )*s0 )
    !IF( ieee_is_nan(y1y2y3) .OR. (.NOT. ieee_is_finite(y1y2y3)) ) ERROR STOP 'y1y2y3 overflow'

  end function y1y2y3

  !> ode_second_dw
  !!
  !! This subroutine solve Equation of the form
  !! s_l''(r) +f_l(r)*s_l(r) = km**2*s_l(r)
  module subroutine ode_second_dw(km, lmax, rc, z, f, s, delta )
    use special_functions, only : coul90, ricbes
    integer, intent(in) :: lmax, z
    real(wp), intent(in)  :: f(0:,0:), rc, km
    real(wp), intent(out) :: s(0:,0:), delta(0:lmax)
    real(wp) :: Big = sqrt(huge(1._wp)), eps = epsilon(1._wp)
    real(wp) :: h, rho, eta
    real(wp), dimension(0:lmax) :: jl, gl, jpl, gpl, fn, betap, beta
    integer :: i, ns, l, i1

    ns = size(s,1)-1
    h = rc/ns

    rho = km*h*(ns-2)

    if(z/=0) then
      eta = -z/km
      call coul90(rho, eta, 0, lmax, jl, gl, jpl, gpl, 0 )
    else
      call ricbes(rho, lmax, jl, gl, jpl, gpl )
    end if


    do l = 0,lmax
      s(0,l) = 0._wp
      do i1 = 1,ns-1
        s(i1,l) = (i1*h)**(l+1)
        if(s(i1,l)>0._wp) exit
      end do
      do i = i1,ns-1
        s(i+1,l) = ( (2._wp+f(i,l)*5._wp*h**2/6.)*s(i,l) &
          -(1._wp-f(i-1,l)*h**2/12.)*s(i-1,l) )/( 1._wp-f(i+1,l)*h**2/12. )
        if(abs(s(i+1,l))>=Big) s(1:i+1,l) = s(1:i+1,l)*eps
      end do
    end do

    betap = ( s(ns-4,:) -8.*s(ns-3,:) +8.*s(ns-1,:) -s(ns,:) ) / ( 12.*h*km )
    beta = s(ns-2,:)
    delta = atan( (-jl*betap+jpl*beta) /(gl*betap-gpl*beta ) )
    fn = (cos(delta)*jl+sin(delta)*gl) /beta

    do l = 0,lmax
      s(:,l) = s(:,l)*fn(l)
    end do

  end subroutine ode_second_dw

  !> Uij
  !!
  !! Calculate the element <\phi(n_i,e_i)| \frac{1}{\vec{r}-\vec{r}_i} |\phi(n_j,e_j>
  !!
  elemental real(wp) function Uij(ni, ei, nj, ej, r)
    use special_functions, only : fac
    real(wp), intent(in) :: ei, ej, r
    integer, intent(in) :: ni, nj
    integer :: k
    real(wp) :: a

    a = 1._wp
    do k = 1,ni+nj-1
      a = a +((ei+ej)*r)**k*(ni+nj-k) /( (ni+nj)*fac(k) )
    end do
    Uij = a*exp(-(ei+ej)*r)*fac(ni+nj)*ei**ni*ej**nj*sqrt(ei*ej) &
      /( sqrt( fac(2*ni)*fac(2*nj) )*( (ei+ej)/2. )**(ni+nj+1) )

  end function Uij

  !> Calculate_U
  !!
  !! This subroutine calculate the short range potential around an atom.
  !! Input:
  !!   Atom: the abbreviation of the atom in two letters (ex: Hy for Hydrogen)
  !!   Orbit: The orbit we calculate the potential for (ex: 2s, 3p...etc)
  !!   r: an array of radius values
  !!   state: charge of the atom (0, 1, 2)
  !! Output:
  !!   U: an array that holds the values of the calculated potential.
  module subroutine calculate_U(Atom, Orbit, r, U, state )
    use input, only : read_orbit
    character(len=2), intent(in) :: Atom, Orbit
    real(wp)   , intent(in) :: r(:)
    real(wp)   , intent(out) :: U(:)
    integer, intent(in) :: state

    integer :: in

    real(wp), allocatable :: a(:),e(:)
    integer, allocatable :: n(:)
    integer :: nelec, lo, no, i1, i2, nocup
    character(len=2) :: orbit_i

    open( newunit = in, file = 'data/'//Atom//'.dat', status = 'old', action = 'read')

    U = 0._wp
    do
      read(in, fmt=*, iostat = lo ) orbit_i
      if(lo<0) exit

      call read_orbit(Atom//'_'//orbit_i, nelec, lo, no, n, a, e)
      nocup = nelec*(2*lo+1)
      if(orbit_i==Orbit ) nocup = nocup -state

      if(nocup==0) cycle

      do i1 = 1,no
        U = U +nocup*a(i1)**2*Uij( n(i1), e(i1), n(i1), e(i1), r )
        if(i1==1) cycle
        do i2 = 1,i1-1
          U = U +2*nocup*a(i1)*a(i2)*Uij( n(i1), e(i1), n(i2), e(i2), r )
        end do
      end do

    end do

    close(in)

  end subroutine calculate_U

  !> INTRPL
  !
  !  REAL(WP) INTERPOLATION OF A SINGLE VALUED FUNCTION
  !  THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
  !  GIVEN  AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
  !  AND FOR A GIVEN SET OF X VALUES(ABSCISSAE),THE VALUES OF
  !  A SINGLE VALUED FUNCTION Y = Y(X).
  !
  !  THE INPUT PARAMETERS ARE:
  !
  !  X = ARRAY OF DIMENSION L STORING THE X VALUES OF INPUT DATA POINTS (IN ASCENDING ORDER)
  !  Y = ARRAY OF DIMENSION L STORING THE Y VALUES OF INPUT DATA POINTS
  !  U = ARRAY OF DIMENSION N STORING THE X VALUES OF THE DESIRED POINTS
  !
  !  THE OUTPUT PARAMETER IS:
  !
  !  V = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y VALUES ARE TO BE DISPLAYED
  pure module subroutine INTRPL(X, Y, U, V )
    use constants, only : wp
    implicit none
    real(wp), intent(in) :: X(:), Y(:), U(:)
    real(wp), intent(out) :: V(:)
    !  DECLARATION STATEMENTS
    real(wp) :: A2, A3, A4, T4, TM2, TM3, TM4, X4, Y4
    integer :: I, IMN, IMX, IPV, J, K, L, N
    real(wp), target :: P0, Q0, Q1, UK, X2, X5, SW, Y2, Y5
    real(wp), pointer :: DX, X3, Y3, T3, A1, TM1, A5, TM5, SA, W2, W4, Q2, W3, Q3

    X3 => P0
    P0 = 0._wp
    Y3 => Q0
    Q0 = 0._wp
    T3 => Q1
    Q1 = 0._wp
    DX => UK
    A1 => X2; TM1 => X2
    X2 = 0._wp
    A5 => X5; TM5 => X5
    X5 = 0._wp
    SA => SW
    W2 => Y2; W4 => Y2; Q2 => Y2
    Y2 = 0._wp
    W3 => Y5; Q3 => Y5
    Y5 = 0._wp


    !  PRELIMINARY PROCESSING
    L = size(X)
    !IF(SIZE(Y)/=L) ERROR STOP 'size(Y)/=size(X)'
    N = size(U)
    !IF(SIZE(V)/=N) ERROR STOP 'INTRPL : size(V)/=size(U)'

    TM4 = 0._wp
    A4 = 0._wp
    A2 = 0._wp

    IPV = 0
    !  MAIN LOOP
    do K = 1,N
      UK = U(K)
      !  ROUTINE TO LOCATE THE DESIRED POINT
      if(UK<X(1)) then
        I = 1
      elseif(UK>=X(L)) then
        I = L+1
      else
        IMN = 2
        IMX = L
        do
          I = (IMN+IMX)/2
          if(UK<X(I)) then
            IMX = I
          else
            IMN = I+1
          end if
          if(IMX<=IMN) exit
        end do
        I = IMX
      end if
      !  CHECK IF I = IPV
      if(I/=IPV) then
        IPV = I
        !  ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND TO ESTIMATE THEM IF NECESSARY
        J = I
        if(J==1) J = 2
        if(J==L+1) J = L
        X3 = X(J-1)
        Y3 = Y(J-1)
        X4 = X(J)
        Y4 = Y(J)
        A3 = X4-X3
        TM3 = (Y4-Y3)/A3
        TM2 = 0._wp
        if( L/=2 ) then
          if( J/=2 ) then
            X2 = X(J-2)
            Y2 = Y(J-2)
            A2 = X3-X2
            TM2 = (Y3-Y2)/A2
          end if
          if(J/=L) then
            X5 = X(J+1)
            Y5 = Y(J+1)
            A4 = X5-X4
            TM4 = (Y5-Y4)/A4
            if(J==2) TM2 = TM3+TM3-TM4
          else
            TM4 = TM3+TM3-TM2
          end if
        else
          TM2 = TM3
        end if
        if(J>3) then
          A1 = X2-X(J-3)
          TM1 = (Y2-Y(J-3))/A1
        else
          TM1 = TM2+TM2-TM3
        end if
        if(J<L-1) then
          A5 = X(J+2)-X5
          TM5 = (Y(J+2)-Y5)/A5
        else
          TM5 = TM4+TM4-TM3
        end if
        !  NUMERICAL DIFFERENTIATION
        if(I/=L+1) then
          W2 = abs(TM4-TM3)
          W3 = abs(TM2-TM1)
          SW = W2+W3
          if(SW==0._wp) then
            W2 = 0.5_wp
            W3 = 0.5_wp
            SW = 1._wp
          end if
          T3 = (W2*TM2+W3*TM3)/SW
        end if
        if(I==1) then
          T4 = T3
          SA = A3+A4
          T3 = 0.5_wp*(TM1+TM2-A4*(A3-A4)*(TM3-TM4)/(SA*SA))
          X3 = X3-A4
          Y3 = Y3-TM2*A4
          A3 = A4
          TM3 = TM2
        else
          W3 = abs(TM5-TM4)
          W4 = abs(TM3-TM2)
          SW = W3+W4
          if(SW==0._wp) then
            W3 = 0.5_wp
            W4 = 0.5_wp
            SW = 1._wp
          end if
          T4 = (W3*TM3+W4*TM4)/SW
          if(I==L+1) then
            T3 = T4
            SA = A2+A3
            T4 = 0.5_wp*(TM4+TM5-A2*(A2-A3)*(TM2-TM3)/(SA*SA))
            X3 = X4
            Y3 = Y4
            A3 = A2
            TM3 = TM4
          end if
        end if
        !  COMPUTATION OF THE POLYNOMIAL
        Q2 = (2._wp*(TM3-T3)+TM3-T4)/A3
        Q3 = (-TM3-TM3+T3+T4)/(A3*A3)
      end if
      DX = UK-P0
      V(K) = Q0+DX*(Q1+DX*(Q2+DX*Q3))

    end do
  end subroutine INTRPL

end submodule utils
