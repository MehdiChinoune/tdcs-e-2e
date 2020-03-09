  !      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM.
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349.
submodule(conhyp_m) conhyp_m
  implicit none
  integer, parameter :: length = 777, bits = digits(1._wp) + 1
  integer, parameter :: min_exp = max( minexponent(1._wp), -74 )
  integer, parameter :: max_exp = min( maxexponent(1._wp), 74 )
  real(wp), parameter :: pi = acos(-1._wp)
  !
contains

  !     ****************************************************************
  !     *                                                              *
  !     *      SOLUTION TO THE CONFLUENT HYPERGEOMETRIC FUNCTION       *
  !     *                                                              *
  !     *                           by                                 *
  !     *                                                              *
  !     *                      MARK NARDIN,                            *
  !     *                                                              *
  !     *              W. F. PERGER and ATUL BHALLA                    *
  !     *                                                              *
  !     *                                                              *
  !     *  Michigan Technological University, Copyright 1989           *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : A numerical evaluator for the confluent       *
  !     *    hypergeometric function for complex arguments with large  *
  !     *    magnitudes using a direct summation of the Kummer series. *
  !     *    The method used allows an accuracy of up to thirteen      *
  !     *    decimal places through the use of large real arrays       *
  !     *    and a single final division.
  !     *    LNCHF is a variable which selects how the result should   *
  !     *    be represented. A '0' will return the value in standard   *
  !     *    exponential form. A '1' will return the LOG of the result.*
  !     *    IP is an integer variable that specifies how many array   *
  !     *    positions are desired (usually 10 is sufficient).         *
  !     *    Setting IP=0 causes the program to estimate the number of *
  !     *    array positions.
  !     *                                                              *
  !     *    The confluent hypergeometric function is the solution to  *
  !     *    the differential equation:                                *
  !     *                                                              *
  !     *             zf"(z) + (a-z)f'(z) - bf(z) = 0                  *
  !     *                                                              *
  !     *  Subprograms called: BITS, CHGF                              *
  !     *                                                              *
  !     ****************************************************************
  elemental complex(wp) module function conhyp(A,B,Z,Lnchf,Ip)
    !
    integer, intent(in) :: Lnchf, Ip
    complex(wp), intent(in) :: A, B, Z
    !
    integer :: i
    real(wp) :: nterm, fx, term1, max, term2, ang
    !
    if ( abs(Z)/=0._wp ) then
      ang = atan2(aimag(Z),real(Z,wp))
    else
      ang = 1._wp
    end if
    if ( abs(ang)<(pi*0.5) ) then
      ang = 1._wp
    else
      ang = sin(abs(ang)-(pi*0.5_wp)) + 1._wp
    end if
    max = 0
    nterm = 0
    fx = 0
    term1 = 0
    do
      nterm = nterm + 1
      term2 = abs((A+nterm-1)*Z/((B+nterm-1)*nterm))
      if ( term2==0._wp ) exit
      if ( term2<1._wp ) then
        if ( (real(A,wp)+nterm-1)>1._wp ) then
          if ( (real(B,wp)+nterm-1)>1._wp ) then
            if ( (term2-term1)<0._wp ) exit
          end if
        end if
      end if
      fx = fx + log(term2)
      if ( fx>max ) max = fx
      term1 = term2
    end do
    max = max*2/(bits*log(2._wp))
    i = int(max*ang) + 7
    if ( i<5 ) i = 5
    if ( Ip>i ) i = Ip
    conhyp = chgf(A,B,Z,i,Lnchf)
    !
  end function conhyp

  !     ****************************************************************
  !     *                                                              *
  !     *                   FUNCTION CHGF                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Function that sums the Kummer series and      *
  !     *    returns the solution of the confluent hypergeometric      *
  !     *    function.                                                 *
  !     *                                                              *
  !     *  Subprograms called: armult, arydiv, cmpadd, cmpmul          *
  !     *                                                              *
  !     ****************************************************************
  elemental complex(wp) function chgf(A,B,Z,L,Lnchf)
    !
    integer, intent(in) :: L, Lnchf
    complex(wp), intent(in) :: A, B, Z
    !
    real(wp) :: ar, ai, cr, ci, xr, xi, cnt, sigfig, mx1, mx2, rmax, ar2, &
      ai2, cr2, ci2, xr2, xi2
    complex(wp) :: final
    !
    real(wp) :: sumr(-1:length), sumi(-1:length), numr(-1:length), &
      numi(-1:length), denomr(-1:length), denomi(-1:length), &
      qr1(-1:length), qr2(-1:length), qi1(-1:length), qi2(-1:length), &
      dumr(-1:length), dumi(-1:length)
    !
    rmax = 2._wp**(bits/2)
    sigfig = 2._wp**(bits/4)
    !
    ar2 = real(A,wp)*sigfig
    ar = aint(ar2)
    ar2 = anint((ar2-ar)*rmax)
    ai2 = aimag(A)*sigfig
    ai = aint(ai2)
    ai2 = anint((ai2-ai)*rmax)
    !
    cr2 = real(B,wp)*sigfig
    cr = aint(cr2)
    cr2 = anint((cr2-cr)*rmax)
    ci2 = aimag(B)*sigfig
    ci = aint(ci2)
    ci2 = anint((ci2-ci)*rmax)
    !
    xr2 = real(Z,wp)*sigfig
    xr = aint(xr2)
    xr2 = anint((xr2-xr)*rmax)
    xi2 = aimag(Z)*sigfig
    xi = aint(xi2)
    xi2 = anint((xi2-xi)*rmax)
    !
    sumr(-1) = 1._wp
    sumi(-1) = 1._wp
    numr(-1) = 1._wp
    numi(-1) = 1._wp
    denomr(-1) = 1._wp
    denomi(-1) = 1._wp
    !
    sumr(0:L+1) = 0._wp
    sumi(0:L+1) = 0._wp
    numr(0:L+1) = 0._wp
    numi(0:L+1) = 0._wp
    denomr(0:L+1) = 0._wp
    denomi(0:L+1) = 0._wp
    !
    sumr(1) = 1._wp
    numr(1) = 1._wp
    denomr(1) = 1._wp
    !
    cnt = sigfig
    do
    if ( sumr(1)<0.5 ) then
      mx1 = sumi(L+1)
    else if ( sumi(1)<0.5 ) then
      mx1 = sumr(L+1)
    else
      mx1 = max(sumr(L+1),sumi(L+1))
    end if
    !
    if ( numr(1)<0.5 ) then
      mx2 = numi(L+1)
    else if ( numi(1)<0.5 ) then
      mx2 = numr(L+1)
    else
      mx2 = max(numr(L+1),numi(L+1))
    end if
    !
    if ( mx1-mx2>2.0 ) then
      if ( cr>0._wp ) then
        if ( abs(cmplx(ar,ai,wp)*cmplx(xr,xi,wp)/(cmplx(cr,ci,wp)*cnt)) &
            <=1._wp ) then
          call arydiv(sumr,sumi,denomr,denomi,final,L,Lnchf,rmax,bits)
          CHGF = final
          return
        end if
      end if
    end if
    call cmpmul(sumr,sumi,cr,ci,qr1,qi1,L,rmax)
    call cmpmul(sumr,sumi,cr2,ci2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    call cmpadd(qr1,qi1,qr2,qi2,sumr,sumi,L,rmax)
    !
    dumr(-1:L+1) = sumr(-1:L+1)
    dumi(-1:L+1) = sumi(-1:L+1)
    call armult(dumr,cnt,sumr,L,rmax)
    call armult(dumi,cnt,sumi,L,rmax)
    call cmpmul(denomr,denomi,cr,ci,qr1,qi1,L,rmax)
    call cmpmul(denomr,denomi,cr2,ci2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    call cmpadd(qr1,qi1,qr2,qi2,denomr,denomi,L,rmax)
    !
    dumr(-1:L+1) = denomr(-1:L+1)
    dumi(-1:L+1) = denomi(-1:L+1)
    call armult(dumr,cnt,denomr,L,rmax)
    call armult(dumi,cnt,denomi,L,rmax)
    call cmpmul(numr,numi,ar,ai,qr1,qi1,L,rmax)
    call cmpmul(numr,numi,ar2,ai2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    call cmpadd(qr1,qi1,qr2,qi2,numr,numi,L,rmax)

    call cmpmul(numr,numi,xr,xi,qr1,qi1,L,rmax)
    call cmpmul(numr,numi,xr2,xi2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    call cmpadd(qr1,qi1,qr2,qi2,numr,numi,L,rmax)
    !
    dumr(-1:L+1) = sumr(-1:L+1)
    dumi(-1:L+1) = sumi(-1:L+1)
    call cmpadd(dumr,dumi,numr,numi,sumr,sumi,L,rmax)
    cnt = cnt + sigfig
    ar = ar + sigfig
    cr = cr + sigfig
    end do
    !
    return
  end function chgf

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ARADD                             *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Accepts two arrays of numbers and returns     *
  !     *    the sum of the array.  Each array is holding the value    *
  !     *    of one number in the series.  The parameter L is the      *
  !     *    size of the array representing the number and RMAX is     *
  !     *    the actual number of digits needed to give the numbers    *
  !     *    the desired accuracy.                                     *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine aradd(A,B,C,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: A(-1:), B(-1:), Rmax
    real(wp), intent(out) :: C(-1:)
    !
    integer :: ediff, i, j
    real(wp) :: z(-1:length)
    !
    z(0:L+1) = 0._wp
    ediff = nint(A(L+1)-B(L+1))
    if ( abs(A(1))<0.5 .OR. ediff<=-L ) then
      C(-1:L+1) = B(-1:L+1)
      goto 311
    else
      if ( abs(B(1))<0.5 .OR. ediff>=L ) then
        C(-1:L+1) = A(-1:L+1)
        goto 311
      else
        z(-1) = A(-1)
        if ( abs(A(-1)-B(-1))>=0.5 ) then
          if ( ediff>0 ) then
            z(L+1) = A(L+1)
            goto 233
          end if
          if ( ediff<0 ) then
            z(L+1) = B(L+1)
            z(-1) = B(-1)
            goto 266
          end if
          do i = 1, L
            if ( A(i)>B(i) ) then
              z(L+1) = A(L+1)
              goto 233
            end if
            if ( A(i)<B(i) ) then
              z(L+1) = B(L+1)
              z(-1) = B(-1)
              goto 266
            end if
          end do

        else if ( ediff>0 ) then
          z(L+1) = A(L+1)
          do i = L, 1 + ediff, -1
            z(i) = A(i) + B(i-ediff) + z(i)
            if ( z(i)>=Rmax ) then
              z(i) = z(i) - Rmax
              z(i-1) = 1._wp
            end if
          end do
          do i = ediff, 1, -1
            z(i) = A(i) + z(i)
            if ( z(i)>=Rmax ) then
              z(i) = z(i) - Rmax
              z(i-1) = 1._wp
            end if
          end do
          if ( z(0)>0.5 ) then
            do i = L, 1, -1
              z(i) = z(i-1)
            end do
            z(L+1) = z(L+1) + 1
            z(0) = 0._wp
          end if
        else if ( ediff<0 ) then
          z(L+1) = B(L+1)
          do i = L, 1 - ediff, -1
            z(i) = A(i+ediff) + B(i) + z(i)
            if ( z(i)>=Rmax ) then
              z(i) = z(i) - Rmax
              z(i-1) = 1._wp
            end if
          end do
          do i = 0 - ediff, 1, -1
            z(i) = B(i) + z(i)
            if ( z(i)>=Rmax ) then
              z(i) = z(i) - Rmax
              z(i-1) = 1._wp
            end if
          end do
          if ( z(0)>0.5 ) then
            do i = L, 1, -1
              z(i) = z(i-1)
            end do
            z(L+1) = z(L+1) + 1._wp
            z(0) = 0._wp
          end if
        else
          z(L+1) = A(L+1)
          do i = L, 1, -1
            z(i) = A(i) + B(i) + z(i)
            if ( z(i)>=Rmax ) then
              z(i) = z(i) - Rmax
              z(i-1) = 1._wp
            end if
          end do
          if ( z(0)>0.5 ) then
            do i = L, 1, -1
              z(i) = z(i-1)
            end do
            z(L+1) = z(L+1) + 1._wp
            z(0) = 0._wp
          end if
        end if
        goto 300

  233   continue
        if ( ediff>0 ) then
          do i = L, 1 + ediff, -1
            z(i) = A(i) - B(i-ediff) + z(i)
            if ( z(i)<0._wp ) then
              z(i) = z(i) + Rmax
              z(i-1) = -1._wp
            end if
          end do
          do i = ediff, 1, -1
            z(i) = A(i) + z(i)
            if ( z(i)<0._wp ) then
              z(i) = z(i) + Rmax
              z(i-1) = -1._wp
            end if
          end do
        else
          do i = L, 1, -1
            z(i) = A(i) - B(i) + z(i)
            if ( z(i)<0._wp ) then
              z(i) = z(i) + Rmax
              z(i-1) = -1._wp
            end if
          end do
        end if
        goto 290
      end if

  266 continue
      if ( ediff<0 ) then
        do i = L, 1 - ediff, -1
          z(i) = B(i) - A(i+ediff) + z(i)
          if ( z(i)<0._wp ) then
            z(i) = z(i) + Rmax
            z(i-1) = -1._wp
          end if
        end do
        do i = 0 - ediff, 1, -1
          z(i) = B(i) + z(i)
          if ( z(i)<0._wp ) then
            z(i) = z(i) + Rmax
            z(i-1) = -1._wp
          end if
        end do
      else
        do i = L, 1, -1
          z(i) = B(i) - A(i) + z(i)
          if ( z(i)<0._wp ) then
            z(i) = z(i) + Rmax
            z(i-1) = -1._wp
          end if
        end do
      end if
    end if

290 continue
    if ( z(1)<=0.5 ) then
      i = 1
      do
        i = i + 1
        if ( z(i)>=0.5 .OR. i>=L+1 ) then
          if ( i==L+1 ) then
            z(-1) = 1._wp
            z(L+1) = 0._wp
            exit
          end if
          do j = 1, L + 1 - i
            z(j) = z(j+i-1)
          end do
          do j = L + 2 - i, L
            z(j) = 0._wp
          end do
          z(L+1) = z(L+1) - i + 1
          exit
        end if
      end do
    end if
300 continue
    C(-1:L+1) = z(-1:L+1)
    311 continue
    if ( C(1)<0.5 ) then
      C(-1) = 1._wp
      C(L+1) = 0._wp
    end if
    !
  end subroutine aradd

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ARSUB                             *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Accepts two arrays and subtracts each element *
  !     *    in the second array from the element in the first array   *
  !     *    and returns the solution.  The parameters L and RMAX are  *
  !     *    the size of the array and the number of digits needed for *
  !     *    the accuracy, respectively.                               *
  !     *                                                              *
  !     *  Subprograms called: aradd                                   *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine arsub(A,B,C,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: A(-1:), B(-1:), Rmax
    real(wp), intent(out) :: C(-1:)
    !
    real(wp) :: b2(-1:length)
    !
    b2(0:L+1) = B(0:L+1)
    b2(-1) = -B(-1)
    call aradd(A,b2,C,L,Rmax)
    !
  end subroutine arsub

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ARMULT                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Accepts two arrays and returns the product.   *
  !     *    L and RMAX are the size of the arrays and the number of   *
  !     *    digits needed to represent the numbers with the required  *
  !     *    accuracy.                                                 *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine armult(A,B,C,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: A(-1:), B, Rmax
    real(wp), intent(out) :: C(-1:)
    !
    integer :: i
    real(wp) :: b2, carry, rmax2
    real(wp) :: z(-1:length)
    !
    rmax2 = 1._wp/Rmax
    z(-1) = sign(1._wp,B)*A(-1)
    b2 = abs(B)
    z(L+1) = A(L+1)
    z(0:L) = 0._wp
    if ( b2<=1.0D-10 .OR. A(1)<=1.0D-10 ) then
      z(-1) = 1._wp
      z(L+1) = 0._wp
      goto 198
    end if
    do i = L, 1, -1
      z(i) = A(i)*b2 + z(i)
      if ( z(i)>=Rmax ) then
        carry = aint(z(i)/Rmax)
        z(i) = z(i) - carry*Rmax
        z(i-1) = carry
      end if
    end do
    if ( z(0)>=0.5 ) then
      do i = L, 1, -1
        z(i) = z(i-1)
      end do
      z(L+1) = z(L+1) + 1._wp
      z(0) = 0._wp
    end if
    !
198 continue
    C(-1:L+1) = z(-1:L+1)
    if ( C(1)<0.5 ) then
      C(-1) = 1._wp
      C(L+1) = 0._wp
    end if
    !
  end subroutine armult

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE CMPADD                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Takes two arrays representing one real and    *
  !     *    one imaginary part, and adds two arrays representing      *
  !     *    another complex number and returns two array holding the  *
  !     *    complex sum.                                              *
  !     *              (CR,CI) = (AR+BR, AI+BI)                        *
  !     *                                                              *
  !     *  Subprograms called: aradd                                   *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine cmpadd(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    real(wp), intent(out) :: Cr(-1:), Ci(-1:)
    !
    call aradd(Ar,Br,Cr,L,Rmax)
    call aradd(Ai,Bi,Ci,L,Rmax)
    !
  end subroutine cmpadd

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE CMPSUB                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Takes two arrays representing one real and    *
  !     *    one imaginary part, and subtracts two arrays representing *
  !     *    another complex number and returns two array holding the  *
  !     *    complex sum.                                              *
  !     *              (CR,CI) = (AR+BR, AI+BI)                        *
  !     *                                                              *
  !     *  Subprograms called: aradd                                   *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine cmpsub(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    real(wp), intent(out) :: Cr(-1:), Ci(-1:)
    !
    call arsub(Ar,Br,Cr,L,Rmax)
    call arsub(Ai,Bi,Ci,L,Rmax)
    !
  end subroutine cmpsub

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE CMPMUL                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Takes two arrays representing one real and    *
  !     *    one imaginary part, and multiplies it with two arrays     *
  !     *    representing another complex number and returns the       *
  !     *    complex product.                                          *
  !     *                                                              *
  !     *  Subprograms called: armult, arsub, aradd                    *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine cmpmul(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    integer, intent(in) :: L
    real(wp), intent(in) :: Ar(-1:), Ai(-1:), Br, Bi, Rmax
    real(wp), intent(out) :: Cr(-1:), Ci(-1:)
    !
    real(wp) :: d1(-1:length), d2(-1:length)
    !
    call armult(Ar,Br,d1,L,Rmax)
    call armult(Ai,Bi,d2,L,Rmax)
    call arsub(d1,d2,Cr,L,Rmax)
    call armult(Ar,Bi,d1,L,Rmax)
    call armult(Ai,Br,d2,L,Rmax)
    call aradd(d1,d2,Ci,L,Rmax)
    !
  end subroutine cmpmul

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ARYDIV                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Returns the double complex number   *
  !     *    resulting from the division of four arrays, representing  *
  !     *    two complex numbers.  The number returned will be in one  *
  !     *    two different forms.  Either standard scientific or as    *
  !     *    the log of the number.                                    *
  !     *                                                              *
  !     *  Subprograms called: conv21, conv12, eadd, ecpdiv, emult     *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine arydiv(Ar,Ai,Br,Bi,C,L,Lnchf,Rmax,Bit)
    !
    integer, intent(in) :: L, Bit, Lnchf
    real(wp), intent(in) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    complex(wp), intent(out) :: C
    !
    integer :: rexp, ir10, ii10
    real(wp) :: phi, n1, n2, n3, e1, e2, e3, rr10, ri10, x, x1, x2, dum1, dum2
    !
    real(wp) :: ae(2,2), be(2,2), ce(2,2)
    !
    rexp = Bit/2
    x = rexp*(Ar(L+1)-2)
    rr10 = x*log10(2._wp)/log10(10._wp)
    ir10 = int(rr10)
    rr10 = rr10 - ir10
    x = rexp*(Ai(L+1)-2)
    ri10 = x*log10(2._wp)/log10(10._wp)
    ii10 = int(ri10)
    ri10 = ri10 - ii10
    dum1 = sign(Ar(1)*Rmax*Rmax+Ar(2)*Rmax+Ar(3),Ar(-1))
    dum2 = sign(Ai(1)*Rmax*Rmax+Ai(2)*Rmax+Ai(3),Ai(-1))
    dum1 = dum1*10**rr10
    dum2 = dum2*10**ri10
    call conv12(cmplx(dum1,dum2,wp),ae)
    ae(1,2) = ae(1,2) + ir10
    ae(2,2) = ae(2,2) + ii10
    x = rexp*(Br(L+1)-2)
    rr10 = x*log10(2._wp)/log10(10._wp)
    ir10 = int(rr10)
    rr10 = rr10 - ir10
    x = rexp*(Bi(L+1)-2)
    ri10 = x*log10(2._wp)/log10(10._wp)
    ii10 = int(ri10)
    ri10 = ri10 - ii10
    dum1 = sign(Br(1)*Rmax*Rmax+Br(2)*Rmax+Br(3),Br(-1))
    dum2 = sign(Bi(1)*Rmax*Rmax+Bi(2)*Rmax+Bi(3),Bi(-1))
    dum1 = dum1*10**rr10
    dum2 = dum2*10**ri10
    call conv12(cmplx(dum1,dum2,wp),be)
    be(1,2) = be(1,2) + ir10
    be(2,2) = be(2,2) + ii10
    call ecpdiv(ae,be,ce)
    if ( Lnchf==0 ) then
      call conv21(ce,C)
    else
      call emult(ce(1,1),ce(1,2),ce(1,1),ce(1,2),n1,e1)
      call emult(ce(2,1),ce(2,2),ce(2,1),ce(2,2),n2,e2)
      call eadd(n1,e1,n2,e2,n3,e3)
      n1 = ce(1,1)
      e1 = ce(1,2) - ce(2,2)
      x2 = ce(2,1)
      if ( e1>max_exp ) then
        x1 = 1._wp * 10._wp**max_exp
      else if ( e1<min_exp ) then
        x1 = 0._wp
      else
        x1 = n1*(10**e1)
      end if
      phi = atan2(x2,x1)
      C = cmplx(0.50D0*(log(n3)+e3*log(10._wp)),phi,wp)
    end if
    !
  end subroutine arydiv

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE EMULT                             *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Takes one base and exponent and multiplies it *
  !     *    by another numbers base and exponent to give the product  *
  !     *    in the form of base and exponent.                         *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  elemental subroutine emult(N1,E1,N2,E2,Nf,Ef)
    !
    real(wp), intent(in) :: N1, E1, N2, E2
    real(wp), intent(out) :: Nf, Ef
    !
    Nf = N1*N2
    Ef = E1 + E2
    if ( abs(Nf)>=10._wp ) then
      Nf = Nf/10._wp
      Ef = Ef + 1._wp
    end if
    !
  end subroutine emult

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE EDIV                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : returns the solution in the form of base and  *
  !     *    exponent of the division of two exponential numbers.      *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  elemental subroutine ediv(N1,E1,N2,E2,Nf,Ef)
    !
    real(wp), intent(in) :: N1, E1, N2, E2
    real(wp), intent(out) :: Nf, Ef
    !
    Nf = N1/N2
    Ef = E1 - E2
    if ( (abs(Nf)<1._wp) .AND. (Nf/=0._wp) ) then
      Nf = Nf*10._wp
      Ef = Ef - 1._wp
    end if
    !
  end subroutine ediv

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE EADD                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Returns the sum of two numbers in the form    *
  !     *    of a base and an exponent.                                *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  elemental subroutine eadd(N1,E1,N2,E2,Nf,Ef)
    !
    real(wp), intent(in) :: N1, E1, N2, E2
    real(wp), intent(out) :: Nf, Ef
    !
    real(wp) :: ediff
    !
    ediff = E1 - E2
    if ( ediff>36._wp ) then
      Nf = N1
      Ef = E1
    else if ( ediff<-36._wp ) then
      Nf = N2
      Ef = E2
    else
      Nf = N1*(10._wp**ediff) + N2
      Ef = E2
      do while ( abs(Nf)>=10._wp )
        Nf = Nf/10._wp
        Ef = Ef + 1._wp
      end do
      do while ( (abs(Nf)<1._wp) .AND. (Nf/=0._wp) )
        Nf = Nf*10._wp
        Ef = Ef - 1._wp
      end do
    end if
    !
  end subroutine eadd

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ESUB                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Returns the solution to the subtraction of    *
  !     *    two numbers in the form of base and exponent.             *
  !     *                                                              *
  !     *  Subprograms called: eadd                                    *
  !     *                                                              *
  !     ****************************************************************
  elemental subroutine esub(N1,E1,N2,E2,Nf,Ef)
    !
    real(wp), intent(in) :: N1, E1, N2, E2
    real(wp), intent(out) :: Nf, Ef
    !
    call eadd(N1,E1,-N2,E2,Nf,Ef)
    !
  end subroutine esub

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE CONV12                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Converts a number from complex notation to a  *
  !     *    form of a 2x2 real array.                                 *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine conv12(Cn,Cae)
    !
    complex(wp), intent(in) :: Cn
    real(wp), intent(out) :: Cae(2,2)
    !
    Cae(1,1) = real(Cn,wp)
    Cae(1,2) = 0._wp
    do while ( abs(Cae(1,1))>=10._wp )
      Cae(1,1) = Cae(1,1)/10._wp
      Cae(1,2) = Cae(1,2) + 1._wp
    end do
    do while ( (abs(Cae(1,1))<1._wp) .AND. (Cae(1,1)/=0._wp) )
      Cae(1,1) = Cae(1,1)*10._wp
      Cae(1,2) = Cae(1,2) - 1._wp
    end do
    Cae(2,1) = aimag(Cn)
    Cae(2,2) = 0._wp
    do while ( abs(Cae(2,1))>=10._wp )
      Cae(2,1) = Cae(2,1)/10._wp
      Cae(2,2) = Cae(2,2) + 1._wp
    end do
    do while ( (abs(Cae(2,1))<1._wp) .AND. (Cae(2,1)/=0._wp) )
      Cae(2,1) = Cae(2,1)*10._wp
      Cae(2,2) = Cae(2,2) - 1._wp
    end do
    !
  end subroutine conv12

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE CONV21                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Converts a number represented in a 2x2 real   *
  !     *    array to the form of a complex number.                    *
  !     *                                                              *
  !     *  Subprograms called: none                                    *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine conv21(Cae,Cn)
    !
    real(wp), intent(in) :: Cae(2,2)
    complex(wp), intent(out) :: Cn
    !
    real(wp) :: a
    !
    if ( Cae(1,2)>max_exp .OR. Cae(2,2)>max_exp ) then
      a = 10._wp**max_exp
      Cn = cmplx( a, a, wp )
    else if ( Cae(2,2)<min_exp ) then
      Cn = cmplx(Cae(1,1)*(10**Cae(1,2)),0._wp,wp)
    else
      Cn = cmplx(Cae(1,1)*(10**Cae(1,2)),Cae(2,1)*(10**Cae(2,2)),wp)
    end if
    !
  end subroutine conv21

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ECPMUL                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Multiplies two numbers which are each         *
  !     *    represented in the form of a two by two array and returns *
  !     *    the solution in the same form.                            *
  !     *                                                              *
  !     *  Subprograms called: emult, esub, eadd                       *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine ecpmul(A,B,C)
    !
    real(wp), intent(in) :: A(2,2), B(2,2)
    real(wp), intent(out) :: C(2,2)
    !
    real(wp) :: n1, e1, n2, e2, c2(2,2)
    !
    call emult(A(1,1),A(1,2),B(1,1),B(1,2),n1,e1)
    call emult(A(2,1),A(2,2),B(2,1),B(2,2),n2,e2)
    call esub(n1,e1,n2,e2,c2(1,1),c2(1,2))
    call emult(A(1,1),A(1,2),B(2,1),B(2,2),n1,e1)
    call emult(A(2,1),A(2,2),B(1,1),B(1,2),n2,e2)
    call eadd(n1,e1,n2,e2,C(2,1),C(2,2))
    C(1,1) = c2(1,1)
    C(1,2) = c2(1,2)
    !
  end subroutine ecpmul

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ECPDIV                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Divides two numbers and returns the solution. *
  !     *    All numbers are represented by a 2x2 array.               *
  !     *                                                              *
  !     *  Subprograms called: eadd, ecpmul, ediv, emult               *
  !     *                                                              *
  !     ****************************************************************
  pure subroutine ecpdiv(A,B,C)
    !
    real(wp), intent(in) :: A(2,2), B(2,2)
    real(wp), intent(out) :: C(2,2)
    !
    real(wp) :: n1, e1, n2, e2, b2(2,2), n3, e3, c2(2,2)
    !
    b2(1,1) = B(1,1)
    b2(1,2) = B(1,2)
    b2(2,1) = -1._wp*B(2,1)
    b2(2,2) = B(2,2)
    call ecpmul(A,b2,c2)
    call emult(B(1,1),B(1,2),B(1,1),B(1,2),n1,e1)
    call emult(B(2,1),B(2,2),B(2,1),B(2,2),n2,e2)
    call eadd(n1,e1,n2,e2,n3,e3)
    call ediv(c2(1,1),c2(1,2),n3,e3,C(1,1),C(1,2))
    call ediv(c2(2,1),c2(2,2),n3,e3,C(2,1),C(2,2))
    !
  end subroutine ecpdiv
  !
end submodule conhyp_m
