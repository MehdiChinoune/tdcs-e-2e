  !      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM.
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349.
SUBMODULE(conhyp_m) conhyp_m
  IMPLICIT NONE
  INTEGER, PARAMETER :: RP = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: length = 777, bits = DIGITS(1._RP) + 1
  REAL(RP), PARAMETER :: pi = ACOS(-1._RP)
  !
CONTAINS

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
  !     *    and a single final division.  LNCHF is a variable which   *
  !     *    selects how the result should be represented.  A '0' will *
  !     *    return the value in standard exponential form.  A '1'     *
  !     *    will return the LOG of the result.  IP is an integer      *
  !     *    variable that specifies how many array positions are      *
  !     *    desired (usually 10 is sufficient).  Setting IP=0 causes  *
  !     *    the program to estimate the number of array positions.    *
  !     *                                                              *
  !     *    The confluent hypergeometric function is the solution to  *
  !     *    the differential equation:                                *
  !     *                                                              *
  !     *             zf"(z) + (a-z)f'(z) - bf(z) = 0                  *
  !     *                                                              *
  !     *  Subprograms called: BITS, CHGF                              *
  !     *                                                              *
  !     ****************************************************************
  ELEMENTAL COMPLEX(RP) MODULE FUNCTION conhyp(A,B,Z,Lnchf,Ip)
    !
    INTEGER, INTENT(IN) :: Lnchf, Ip
    COMPLEX(RP), INTENT(IN) :: A, B, Z
    !
    INTEGER :: i
    REAL(RP) :: nterm, fx, term1, MAX, term2, ang
    !
    IF ( ABS(Z)/=0._RP ) THEN
      ang = ATAN2(AIMAG(Z),REAL(Z,RP))
    ELSE
      ang = 1._RP
    END IF
    IF ( ABS(ang)<(pi*0.5) ) THEN
      ang = 1._RP
    ELSE
      ang = SIN(ABS(ang)-(pi*0.5_RP)) + 1._RP
    END IF
    MAX = 0
    nterm = 0
    fx = 0
    term1 = 0
    DO
      nterm = nterm + 1
      term2 = ABS((A+nterm-1)*Z/((B+nterm-1)*nterm))
      IF ( term2==0._RP ) EXIT
      IF ( term2<1._RP ) THEN
        IF ( (REAL(A,RP)+nterm-1)>1._RP ) THEN
          IF ( (REAL(B,RP)+nterm-1)>1._RP ) THEN
            IF ( (term2-term1)<0._RP ) EXIT
          END IF
        END IF
      END IF
      fx = fx + LOG(term2)
      IF ( fx>MAX ) MAX = fx
      term1 = term2
    END DO
    MAX = MAX*2/(bits*LOG(2._RP))
    i = INT(MAX*ang) + 7
    IF ( i<5 ) i = 5
    IF ( Ip>i ) i = Ip
    conhyp = chgf(A,B,Z,i,Lnchf)
    !
  END FUNCTION conhyp

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
  ELEMENTAL COMPLEX(RP) FUNCTION chgf(A,B,Z,L,Lnchf)
    !
    INTEGER, INTENT(IN) :: L, Lnchf
    COMPLEX(RP), INTENT(IN) :: A, B, Z
    !
    REAL(RP) :: ar, ai, cr, ci, xr, xi, cnt, sigfig, mx1, mx2, rmax, ar2, &
      ai2, cr2, ci2, xr2, xi2
    COMPLEX(RP) :: FINAL
    !
    REAL(RP) :: sumr(-1:length), sumi(-1:length), numr(-1:length), &
      numi(-1:length), denomr(-1:length), denomi(-1:length), &
      qr1(-1:length), qr2(-1:length), qi1(-1:length), qi2(-1:length), &
      dumr(-1:length), dumi(-1:length)
    !
    rmax = 2._RP**(bits/2)
    sigfig = 2._RP**(bits/4)
    !
    ar2 = REAL(A,RP)*sigfig
    ar = AINT(ar2)
    ar2 = ANINT((ar2-ar)*rmax)
    ai2 = AIMAG(A)*sigfig
    ai = AINT(ai2)
    ai2 = ANINT((ai2-ai)*rmax)
    !
    cr2 = REAL(B,RP)*sigfig
    cr = AINT(cr2)
    cr2 = ANINT((cr2-cr)*rmax)
    ci2 = AIMAG(B)*sigfig
    ci = AINT(ci2)
    ci2 = ANINT((ci2-ci)*rmax)
    !
    xr2 = REAL(Z,RP)*sigfig
    xr = AINT(xr2)
    xr2 = ANINT((xr2-xr)*rmax)
    xi2 = AIMAG(Z)*sigfig
    xi = AINT(xi2)
    xi2 = ANINT((xi2-xi)*rmax)
    !
    sumr(-1) = 1._RP
    sumi(-1) = 1._RP
    numr(-1) = 1._RP
    numi(-1) = 1._RP
    denomr(-1) = 1._RP
    denomi(-1) = 1._RP
    !
    sumr(0:L+1) = 0._RP
    sumi(0:L+1) = 0._RP
    numr(0:L+1) = 0._RP
    numi(0:L+1) = 0._RP
    denomr(0:L+1) = 0._RP
    denomi(0:L+1) = 0._RP
    !
    sumr(1) = 1._RP
    numr(1) = 1._RP
    denomr(1) = 1._RP
    !
    cnt = sigfig
    DO
    IF ( sumr(1)<0.5 ) THEN
      mx1 = sumi(L+1)
    ELSE IF ( sumi(1)<0.5 ) THEN
      mx1 = sumr(L+1)
    ELSE
      mx1 = MAX(sumr(L+1),sumi(L+1))
    END IF
    !
    IF ( numr(1)<0.5 ) THEN
      mx2 = numi(L+1)
    ELSE IF ( numi(1)<0.5 ) THEN
      mx2 = numr(L+1)
    ELSE
      mx2 = MAX(numr(L+1),numi(L+1))
    END IF
    !
    IF ( mx1-mx2>2.0 ) THEN
      IF ( cr>0._RP ) THEN
        IF ( ABS(CMPLX(ar,ai,RP)*CMPLX(xr,xi,RP)/(CMPLX(cr,ci,RP)*cnt)) &
            <=1._RP ) THEN
          CALL arydiv(sumr,sumi,denomr,denomi,FINAL,L,Lnchf,rmax,bits)
          CHGF = FINAL
          RETURN
        END IF
      END IF
    END IF
    CALL cmpmul(sumr,sumi,cr,ci,qr1,qi1,L,rmax)
    CALL cmpmul(sumr,sumi,cr2,ci2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    CALL cmpadd(qr1,qi1,qr2,qi2,sumr,sumi,L,rmax)
    !
    dumr(-1:L+1) = sumr(-1:L+1)
    dumi(-1:L+1) = sumi(-1:L+1)
    CALL armult(dumr,cnt,sumr,L,rmax)
    CALL armult(dumi,cnt,sumi,L,rmax)
    CALL cmpmul(denomr,denomi,cr,ci,qr1,qi1,L,rmax)
    CALL cmpmul(denomr,denomi,cr2,ci2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    CALL cmpadd(qr1,qi1,qr2,qi2,denomr,denomi,L,rmax)
    !
    dumr(-1:L+1) = denomr(-1:L+1)
    dumi(-1:L+1) = denomi(-1:L+1)
    CALL armult(dumr,cnt,denomr,L,rmax)
    CALL armult(dumi,cnt,denomi,L,rmax)
    CALL cmpmul(numr,numi,ar,ai,qr1,qi1,L,rmax)
    CALL cmpmul(numr,numi,ar2,ai2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    CALL cmpadd(qr1,qi1,qr2,qi2,numr,numi,L,rmax)

    CALL cmpmul(numr,numi,xr,xi,qr1,qi1,L,rmax)
    CALL cmpmul(numr,numi,xr2,xi2,qr2,qi2,L,rmax)
    qr2(L+1) = qr2(L+1) - 1
    qi2(L+1) = qi2(L+1) - 1
    CALL cmpadd(qr1,qi1,qr2,qi2,numr,numi,L,rmax)
    !
    dumr(-1:L+1) = sumr(-1:L+1)
    dumi(-1:L+1) = sumi(-1:L+1)
    CALL cmpadd(dumr,dumi,numr,numi,sumr,sumi,L,rmax)
    cnt = cnt + sigfig
    ar = ar + sigfig
    cr = cr + sigfig
    END DO
    !
    RETURN
  END FUNCTION chgf

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
  PURE SUBROUTINE aradd(A,B,C,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:), B(-1:), Rmax
    REAL(RP), INTENT(OUT) :: C(-1:)
    !
    INTEGER :: ediff, i, j
    REAL(RP) :: z(-1:length)
    !
    z(0:L+1) = 0._RP
    ediff = NINT(A(L+1)-B(L+1))
    IF ( ABS(A(1))<0.5 .OR. ediff<=-L ) THEN
      C(-1:L+1) = B(-1:L+1)
      GOTO 311
    ELSE
      IF ( ABS(B(1))<0.5 .OR. ediff>=L ) THEN
        C(-1:L+1) = A(-1:L+1)
        GOTO 311
      ELSE
        z(-1) = A(-1)
        IF ( ABS(A(-1)-B(-1))>=0.5 ) THEN
          IF ( ediff>0 ) THEN
            z(L+1) = A(L+1)
            GOTO 233
          END IF
          IF ( ediff<0 ) THEN
            z(L+1) = B(L+1)
            z(-1) = B(-1)
            GOTO 266
          END IF
          DO i = 1, L
            IF ( A(i)>B(i) ) THEN
              z(L+1) = A(L+1)
              GOTO 233
            END IF
            IF ( A(i)<B(i) ) THEN
              z(L+1) = B(L+1)
              z(-1) = B(-1)
              GOTO 266
            END IF
          END DO

        ELSE IF ( ediff>0 ) THEN
          z(L+1) = A(L+1)
          DO i = L, 1 + ediff, -1
            z(i) = A(i) + B(i-ediff) + z(i)
            IF ( z(i)>=Rmax ) THEN
              z(i) = z(i) - Rmax
              z(i-1) = 1._RP
            END IF
          END DO
          DO i = ediff, 1, -1
            z(i) = A(i) + z(i)
            IF ( z(i)>=Rmax ) THEN
              z(i) = z(i) - Rmax
              z(i-1) = 1._RP
            END IF
          END DO
          IF ( z(0)>0.5 ) THEN
            DO i = L, 1, -1
              z(i) = z(i-1)
            END DO
            z(L+1) = z(L+1) + 1
            z(0) = 0._RP
          END IF
        ELSE IF ( ediff<0 ) THEN
          z(L+1) = B(L+1)
          DO i = L, 1 - ediff, -1
            z(i) = A(i+ediff) + B(i) + z(i)
            IF ( z(i)>=Rmax ) THEN
              z(i) = z(i) - Rmax
              z(i-1) = 1._RP
            END IF
          END DO
          DO i = 0 - ediff, 1, -1
            z(i) = B(i) + z(i)
            IF ( z(i)>=Rmax ) THEN
              z(i) = z(i) - Rmax
              z(i-1) = 1._RP
            END IF
          END DO
          IF ( z(0)>0.5 ) THEN
            DO i = L, 1, -1
              z(i) = z(i-1)
            END DO
            z(L+1) = z(L+1) + 1._RP
            z(0) = 0._RP
          END IF
        ELSE
          z(L+1) = A(L+1)
          DO i = L, 1, -1
            z(i) = A(i) + B(i) + z(i)
            IF ( z(i)>=Rmax ) THEN
              z(i) = z(i) - Rmax
              z(i-1) = 1._RP
            END IF
          END DO
          IF ( z(0)>0.5 ) THEN
            DO i = L, 1, -1
              z(i) = z(i-1)
            END DO
            z(L+1) = z(L+1) + 1._RP
            z(0) = 0._RP
          END IF
        END IF
        GOTO 300

  233   CONTINUE
        IF ( ediff>0 ) THEN
          DO i = L, 1 + ediff, -1
            z(i) = A(i) - B(i-ediff) + z(i)
            IF ( z(i)<0._RP ) THEN
              z(i) = z(i) + Rmax
              z(i-1) = -1._RP
            END IF
          END DO
          DO i = ediff, 1, -1
            z(i) = A(i) + z(i)
            IF ( z(i)<0._RP ) THEN
              z(i) = z(i) + Rmax
              z(i-1) = -1._RP
            END IF
          END DO
        ELSE
          DO i = L, 1, -1
            z(i) = A(i) - B(i) + z(i)
            IF ( z(i)<0._RP ) THEN
              z(i) = z(i) + Rmax
              z(i-1) = -1._RP
            END IF
          END DO
        END IF
        GOTO 290
      END IF

  266 CONTINUE
      IF ( ediff<0 ) THEN
        DO i = L, 1 - ediff, -1
          z(i) = B(i) - A(i+ediff) + z(i)
          IF ( z(i)<0._RP ) THEN
            z(i) = z(i) + Rmax
            z(i-1) = -1._RP
          END IF
        END DO
        DO i = 0 - ediff, 1, -1
          z(i) = B(i) + z(i)
          IF ( z(i)<0._RP ) THEN
            z(i) = z(i) + Rmax
            z(i-1) = -1._RP
          END IF
        END DO
      ELSE
        DO i = L, 1, -1
          z(i) = B(i) - A(i) + z(i)
          IF ( z(i)<0._RP ) THEN
            z(i) = z(i) + Rmax
            z(i-1) = -1._RP
          END IF
        END DO
      END IF
    END IF

290 CONTINUE
    IF ( z(1)<=0.5 ) THEN
      i = 1
      DO
        i = i + 1
        IF ( z(i)>=0.5 .OR. i>=L+1 ) THEN
          IF ( i==L+1 ) THEN
            z(-1) = 1._RP
            z(L+1) = 0._RP
            EXIT
          END IF
          DO j = 1, L + 1 - i
            z(j) = z(j+i-1)
          END DO
          DO j = L + 2 - i, L
            z(j) = 0._RP
          END DO
          z(L+1) = z(L+1) - i + 1
          EXIT
        END IF
      END DO
    END IF
300 CONTINUE
    C(-1:L+1) = z(-1:L+1)
    311 CONTINUE
    IF ( C(1)<0.5 ) THEN
      C(-1) = 1._RP
      C(L+1) = 0._RP
    END IF
    !
  END SUBROUTINE aradd

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
  PURE SUBROUTINE arsub(A,B,C,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:), B(-1:), Rmax
    REAL(RP), INTENT(OUT) :: C(-1:)
    !
    REAL(RP) :: b2(-1:length)
    !
    b2(0:L+1) = B(0:L+1)
    b2(-1) = -B(-1)
    CALL aradd(A,b2,C,L,Rmax)
    !
  END SUBROUTINE arsub

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
  PURE SUBROUTINE armult(A,B,C,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:), B, Rmax
    REAL(RP), INTENT(OUT) :: C(-1:)
    !
    INTEGER :: i
    REAL(RP) :: b2, carry, rmax2
    REAL(RP) :: z(-1:length)
    !
    rmax2 = 1._RP/Rmax
    z(-1) = SIGN(1._RP,B)*A(-1)
    b2 = ABS(B)
    z(L+1) = A(L+1)
    z(0:L) = 0._RP
    IF ( b2<=1.0D-10 .OR. A(1)<=1.0D-10 ) THEN
      z(-1) = 1._RP
      z(L+1) = 0._RP
      GOTO 198
    END IF
    DO i = L, 1, -1
      z(i) = A(i)*b2 + z(i)
      IF ( z(i)>=Rmax ) THEN
        carry = AINT(z(i)/Rmax)
        z(i) = z(i) - carry*Rmax
        z(i-1) = carry
      END IF
    END DO
    IF ( z(0)>=0.5 ) THEN
      DO i = L, 1, -1
        z(i) = z(i-1)
      END DO
      z(L+1) = z(L+1) + 1._RP
      z(0) = 0._RP
    END IF
    !
198 CONTINUE
    C(-1:L+1) = z(-1:L+1)
    IF ( C(1)<0.5 ) THEN
      C(-1) = 1._RP
      C(L+1) = 0._RP
    END IF
    !
  END SUBROUTINE armult

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
  PURE SUBROUTINE cmpadd(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    REAL(RP), INTENT(OUT) :: Cr(-1:), Ci(-1:)
    !
    CALL aradd(Ar,Br,Cr,L,Rmax)
    CALL aradd(Ai,Bi,Ci,L,Rmax)
    !
  END SUBROUTINE cmpadd

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
  PURE SUBROUTINE cmpsub(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    REAL(RP), INTENT(OUT) :: Cr(-1:), Ci(-1:)
    !
    CALL arsub(Ar,Br,Cr,L,Rmax)
    CALL arsub(Ai,Bi,Ci,L,Rmax)
    !
  END SUBROUTINE cmpsub

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
  PURE SUBROUTINE cmpmul(Ar,Ai,Br,Bi,Cr,Ci,L,Rmax)
    !
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: Ar(-1:), Ai(-1:), Br, Bi, Rmax
    REAL(RP), INTENT(OUT) :: Cr(-1:), Ci(-1:)
    !
    REAL(RP) :: d1(-1:length), d2(-1:length)
    !
    CALL armult(Ar,Br,d1,L,Rmax)
    CALL armult(Ai,Bi,d2,L,Rmax)
    CALL arsub(d1,d2,Cr,L,Rmax)
    CALL armult(Ar,Bi,d1,L,Rmax)
    CALL armult(Ai,Br,d2,L,Rmax)
    CALL aradd(d1,d2,Ci,L,Rmax)
    !
  END SUBROUTINE cmpmul

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
  PURE SUBROUTINE arydiv(Ar,Ai,Br,Bi,C,L,Lnchf,Rmax,Bit)
    !
    INTEGER, INTENT(IN) :: L, Bit, Lnchf
    REAL(RP), INTENT(IN) :: Ar(-1:), Ai(-1:), Br(-1:), Bi(-1:), Rmax
    COMPLEX(RP), INTENT(OUT) :: C
    !
    INTEGER :: rexp, ir10, ii10
    REAL(RP) :: phi, n1, n2, n3, e1, e2, e3, rr10, ri10, x, x1, x2, dum1, dum2
    !
    REAL(RP) :: ae(2,2), be(2,2), ce(2,2)
    !
    rexp = Bit/2
    x = rexp*(Ar(L+1)-2)
    rr10 = x*LOG10(2._RP)/LOG10(10._RP)
    ir10 = INT(rr10)
    rr10 = rr10 - ir10
    x = rexp*(Ai(L+1)-2)
    ri10 = x*LOG10(2._RP)/LOG10(10._RP)
    ii10 = INT(ri10)
    ri10 = ri10 - ii10
    dum1 = SIGN(Ar(1)*Rmax*Rmax+Ar(2)*Rmax+Ar(3),Ar(-1))
    dum2 = SIGN(Ai(1)*Rmax*Rmax+Ai(2)*Rmax+Ai(3),Ai(-1))
    dum1 = dum1*10**rr10
    dum2 = dum2*10**ri10
    CALL conv12(CMPLX(dum1,dum2,RP),ae)
    ae(1,2) = ae(1,2) + ir10
    ae(2,2) = ae(2,2) + ii10
    x = rexp*(Br(L+1)-2)
    rr10 = x*LOG10(2._RP)/LOG10(10._RP)
    ir10 = INT(rr10)
    rr10 = rr10 - ir10
    x = rexp*(Bi(L+1)-2)
    ri10 = x*LOG10(2._RP)/LOG10(10._RP)
    ii10 = INT(ri10)
    ri10 = ri10 - ii10
    dum1 = SIGN(Br(1)*Rmax*Rmax+Br(2)*Rmax+Br(3),Br(-1))
    dum2 = SIGN(Bi(1)*Rmax*Rmax+Bi(2)*Rmax+Bi(3),Bi(-1))
    dum1 = dum1*10**rr10
    dum2 = dum2*10**ri10
    CALL conv12(CMPLX(dum1,dum2,RP),be)
    be(1,2) = be(1,2) + ir10
    be(2,2) = be(2,2) + ii10
    CALL ecpdiv(ae,be,ce)
    IF ( Lnchf==0 ) THEN
      CALL conv21(ce,C)
    ELSE
      CALL emult(ce(1,1),ce(1,2),ce(1,1),ce(1,2),n1,e1)
      CALL emult(ce(2,1),ce(2,2),ce(2,1),ce(2,2),n2,e2)
      CALL eadd(n1,e1,n2,e2,n3,e3)
      n1 = ce(1,1)
      e1 = ce(1,2) - ce(2,2)
      x2 = ce(2,1)
      IF ( e1>74._RP ) THEN
        x1 = 1.E75_RP
      ELSE IF ( e1<-74._RP ) THEN
        x1 = 0
      ELSE
        x1 = n1*(10**e1)
      END IF
      phi = ATAN2(x2,x1)
      C = CMPLX(0.50D0*(LOG(n3)+e3*LOG(10._RP)),phi,RP)
    END IF
    !
  END SUBROUTINE arydiv

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
  ELEMENTAL SUBROUTINE emult(N1,E1,N2,E2,Nf,Ef)
    !
    REAL(RP), INTENT(IN) :: N1, E1, N2, E2
    REAL(RP), INTENT(OUT) :: Nf, Ef
    !
    Nf = N1*N2
    Ef = E1 + E2
    IF ( ABS(Nf)>=10._RP ) THEN
      Nf = Nf/10._RP
      Ef = Ef + 1._RP
    END IF
    !
  END SUBROUTINE emult

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
  ELEMENTAL SUBROUTINE ediv(N1,E1,N2,E2,Nf,Ef)
    !
    REAL(RP), INTENT(IN) :: N1, E1, N2, E2
    REAL(RP), INTENT(OUT) :: Nf, Ef
    !
    Nf = N1/N2
    Ef = E1 - E2
    IF ( (ABS(Nf)<1._RP) .AND. (Nf/=0._RP) ) THEN
      Nf = Nf*10._RP
      Ef = Ef - 1._RP
    END IF
    !
  END SUBROUTINE ediv

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
  ELEMENTAL SUBROUTINE eadd(N1,E1,N2,E2,Nf,Ef)
    !
    REAL(RP), INTENT(IN) :: N1, E1, N2, E2
    REAL(RP), INTENT(OUT) :: Nf, Ef
    !
    REAL(RP) :: ediff
    !
    ediff = E1 - E2
    IF ( ediff>36._RP ) THEN
      Nf = N1
      Ef = E1
    ELSE IF ( ediff<-36._RP ) THEN
      Nf = N2
      Ef = E2
    ELSE
      Nf = N1*(10._RP**ediff) + N2
      Ef = E2
      DO WHILE ( ABS(Nf)>=10._RP )
        Nf = Nf/10._RP
        Ef = Ef + 1._RP
      END DO
      DO WHILE ( (ABS(Nf)<1._RP) .AND. (Nf/=0._RP) )
        Nf = Nf*10._RP
        Ef = Ef - 1._RP
      END DO
    END IF
    !
  END SUBROUTINE eadd

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
  ELEMENTAL SUBROUTINE esub(N1,E1,N2,E2,Nf,Ef)
    !
    REAL(RP), INTENT(IN) :: N1, E1, N2, E2
    REAL(RP), INTENT(OUT) :: Nf, Ef
    !
    CALL eadd(N1,E1,-N2,E2,Nf,Ef)
    !
  END SUBROUTINE esub

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
  PURE SUBROUTINE conv12(Cn,Cae)
    !
    COMPLEX(RP), INTENT(IN) :: Cn
    REAL(RP), INTENT(OUT) :: Cae(2,2)
    !
    Cae(1,1) = REAL(Cn,RP)
    Cae(1,2) = 0._RP
    DO WHILE ( ABS(Cae(1,1))>=10._RP )
      Cae(1,1) = Cae(1,1)/10._RP
      Cae(1,2) = Cae(1,2) + 1._RP
    END DO
    DO WHILE ( (ABS(Cae(1,1))<1._RP) .AND. (Cae(1,1)/=0._RP) )
      Cae(1,1) = Cae(1,1)*10._RP
      Cae(1,2) = Cae(1,2) - 1._RP
    END DO
    Cae(2,1) = AIMAG(Cn)
    Cae(2,2) = 0._RP
    DO WHILE ( ABS(Cae(2,1))>=10._RP )
      Cae(2,1) = Cae(2,1)/10._RP
      Cae(2,2) = Cae(2,2) + 1._RP
    END DO
    DO WHILE ( (ABS(Cae(2,1))<1._RP) .AND. (Cae(2,1)/=0._RP) )
      Cae(2,1) = Cae(2,1)*10._RP
      Cae(2,2) = Cae(2,2) - 1._RP
    END DO
    !
  END SUBROUTINE conv12

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
  PURE SUBROUTINE conv21(Cae,Cn)
    !
    REAL(RP), INTENT(IN) :: Cae(2,2)
    COMPLEX(RP), INTENT(OUT) :: Cn
    !
    IF ( Cae(1,2)>75 .OR. Cae(2,2)>75 ) THEN
      Cn = CMPLX(1.E75_RP,1.E75_RP,RP)
    ELSE IF ( Cae(2,2)<-75 ) THEN
      Cn = CMPLX(Cae(1,1)*(10**Cae(1,2)),0D0,RP)
    ELSE
      Cn = CMPLX(Cae(1,1)*(10**Cae(1,2)),Cae(2,1)*(10**Cae(2,2)),RP)
    END IF
    !
  END SUBROUTINE conv21

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
  PURE SUBROUTINE ecpmul(A,B,C)
    !
    REAL(RP), INTENT(IN) :: A(2,2), B(2,2)
    REAL(RP), INTENT(OUT) :: C(2,2)
    !
    REAL(RP) :: n1, e1, n2, e2, c2(2,2)
    !
    CALL emult(A(1,1),A(1,2),B(1,1),B(1,2),n1,e1)
    CALL emult(A(2,1),A(2,2),B(2,1),B(2,2),n2,e2)
    CALL esub(n1,e1,n2,e2,c2(1,1),c2(1,2))
    CALL emult(A(1,1),A(1,2),B(2,1),B(2,2),n1,e1)
    CALL emult(A(2,1),A(2,2),B(1,1),B(1,2),n2,e2)
    CALL eadd(n1,e1,n2,e2,C(2,1),C(2,2))
    C(1,1) = c2(1,1)
    C(1,2) = c2(1,2)
    !
  END SUBROUTINE ecpmul

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
  PURE SUBROUTINE ecpdiv(A,B,C)
    !
    REAL(RP), INTENT(IN) :: A(2,2), B(2,2)
    REAL(RP), INTENT(OUT) :: C(2,2)
    !
    REAL(RP) :: n1, e1, n2, e2, b2(2,2), n3, e3, c2(2,2)
    !
    b2(1,1) = B(1,1)
    b2(1,2) = B(1,2)
    b2(2,1) = -1._RP*B(2,1)
    b2(2,2) = B(2,2)
    CALL ecpmul(A,b2,c2)
    CALL emult(B(1,1),B(1,2),B(1,1),B(1,2),n1,e1)
    CALL emult(B(2,1),B(2,2),B(2,1),B(2,2),n2,e2)
    CALL eadd(n1,e1,n2,e2,n3,e3)
    CALL ediv(c2(1,1),c2(1,2),n3,e3,C(1,1),C(1,2))
    CALL ediv(c2(2,1),c2(2,2),n3,e3,C(2,1),C(2,2))
    !
  END SUBROUTINE ecpdiv
  !
END SUBMODULE conhyp_m