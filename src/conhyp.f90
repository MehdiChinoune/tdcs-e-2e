MODULE conhyp_m
  USE constants ,only: RP
  IMPLICIT NONE

CONTAINS
  !      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM.
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349.
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

  COMPLEX(RP) FUNCTION CONHYP (A,B,Z,LNCHF,IP)
    COMPLEX(RP), INTENT(IN) :: A, B, Z
    INTEGER, INTENT(IN), OPTIONAL :: LNCHF, IP
    INTEGER :: I
    REAL(RP) ::  NTERM,FX,TERM1,MAX,TERM2,ANG

    IF( ABS(Z) /= 0._RP) THEN
      ANG = ATAN2(AIMAG(Z),REAL(Z,RP))
    ELSE
      ANG = 1._RP
    ENDIF
    IF(ABS(ANG) < (3.14159_RP*0.5)) THEN
      ANG=1._RP
    ELSE
      ANG=SIN(ABS(ANG)-(3.14159265_RP*0.5_RP))+1._RP
    ENDIF
    MAX=0
    NTERM=0
    FX=0
    TERM1=0
    DO
      NTERM=NTERM+1
      TERM2=ABS((A+NTERM-1)*Z/((B+NTERM-1)*NTERM))
      IF(TERM2 == 0._RP) EXIT
      IF(TERM2 < 1._RP .AND. (REAL(A)+NTERM-1) > 1._RP .AND. (REAL(B)+NTERM-1) > 1._RP &
        .AND. (TERM2-TERM1) < 0._RP) EXIT
      FX=FX+LOG(TERM2)
      IF(FX > MAX) MAX=FX
      TERM1=TERM2
    END DO
    MAX = MAX*2/(BITS()*6.93147181E-1_RP)
    I = INT(MAX*ANG)+7
    IF(I < 5) I = 5

    IF( PRESENT(LNCHF) ) THEN
      IF( PRESENT(IP) .AND. IP>I ) THEN
        I = IP
        CONHYP = CHGF(A,B,Z,I,LNCHF)
      ELSE
        CONHYP = CHGF(A,B,Z,10,LNCHF)
      END IF
    ELSE
      CONHYP = CHGF(A,B,Z,10,0)
    END IF

    RETURN
  END FUNCTION CONHYP

  !     ****************************************************************
  !     *                                                              *
  !     *                   FUNCTION BITS                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Determines the number of significant figures  *
  !     *    of machine precision to arrive at the size of the array   *
  !     *    the numbers must must be stored in to get the accuracy    *
  !     *    of the solution.                                          *
  !     *                                                              *
  !     *  Subprogram called: STORE                                    *
  !     *                                                              *
  !     ****************************************************************

  INTEGER FUNCTION BITS()
    REAL(RP) ::  BIT,BIT2

    BIT = 1._RP
    BITS = 0
    DO
      BITS = BITS +1
      BIT2 = BIT*2._RP
      BIT = BIT2 + 1._RP
      IF( (BIT-BIT2)==0._RP ) EXIT
    END DO

    RETURN
  END FUNCTION BITS

  !***********************************************************
  !
  !
  !   This function forces its argument X to be stored in a
  ! memory location, thus providing a means of determining
  ! floating point number characteristics (such as the machine
  ! precision) when it is necessary to avoid computation in
  ! high precision registers.
  !
  ! On input:
  !
  !       X = Value to be stored.
  !
  ! X is not altered by this function.
  !
  ! On output:
  !
  !       STORE = Value of X after it has been stored and
  !               possibly truncated or rounded to the double
  !               precision word length.
  !
  ! Modules required by STORE:  None
  !
  !***********************************************************
  REAL(RP) FUNCTION STORE(X)
    REAL(RP), INTENT(IN) ::  X
    REAL(RP) ::  Y
    COMMON/STCOM/Y
    Y = X
    STORE = Y
    RETURN
  END FUNCTION STORE

  !     ****************************************************************
  !     *                                                              *
  !     *                   FUNCTION CHGF                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Function that sums the Kummer series and      *
  !     *    returns the solution of the confluent hypergeometric      *
  !     *    function.                                                 *
  !     *                                                              *
  !     *  Subprograms called: ARMULT, ARYDIV, BITS, CMPADD, CMPMUL    *
  !     *                                                              *
  !     ****************************************************************

  COMPLEX(RP) FUNCTION CHGF (A,B,Z,L,LNCHF)
    COMPLEX(RP), INTENT(IN) :: A, B, Z
    INTEGER, INTENT(IN) :: L, LNCHF

    INTEGER :: BIT
    COMPLEX(RP) :: FINAL
    REAL(RP) :: AR,AI,CR,CI,XR,XI,CNT,SIGFIG,MX1,MX2
    REAL(RP) :: RMAX
    REAL(RP) :: AR2,AI2,CR2,CI2,XR2,XI2
    REAL(RP) :: SUMR(-1:L+1), SUMI(-1:L+1), DENOMR(-1:L+1), DENOMI(-1:L+1)
    REAL(RP) :: SUM_tmp(-1:L+1), SUM_tmp2(-1:L+1)
    REAL(RP) :: NUMR(-1:L+1), NUMI(-1:L+1)
    REAL(RP) :: QR1(-1:L+1),QR2(-1:L+1),QI1(-1:L+1),QI2(-1:L+1)

    BIT = BITS()
    RMAX = 2._RP**(BIT/2)
    SIGFIG = 2._RP**(BIT/4)
    AR2 = REAL(A,RP)*SIGFIG
    AR = INT(AR2)
    AR2 = INT((AR2-AR)*RMAX)
    AI2 = AIMAG(A)*SIGFIG
    AI = INT(AI2)
    AI2 = INT((AI2-AI)*RMAX)
    CR2 = REAL(B,RP)*SIGFIG
    CR = INT(CR2)
    CR2 = INT((CR2-CR)*RMAX)
    CI2 = AIMAG(B)*SIGFIG
    CI = INT(CI2)
    CI2 = INT((CI2-CI)*RMAX)
    XR2 = REAL(Z,RP)*SIGFIG
    XR = INT(XR2)
    XR2 = INT((XR2-XR)*RMAX)
    XI2 = AIMAG(Z)*SIGFIG
    XI = INT(XI2)
    XI2 = INT((XI2-XI)*RMAX)
    SUMR(-1) = 1._RP
    SUMI(-1) = 1._RP
    NUMR(-1) = 1._RP
    NUMI(-1) = 1._RP
    DENOMR(-1) = 1._RP
    DENOMI(-1) = 1._RP
    SUMR(0:L+1) = 0._RP
    SUMI(0:L+1) = 0._RP
    NUMR(0:L+1) = 0._RP
    NUMI(0:L+1) = 0._RP
    DENOMR(0:L+1) = 0._RP
    DENOMI(0:L+1) = 0._RP
    SUMR(1) = 1._RP
    NUMR(1) = 1._RP
    DENOMR(1) = 1._RP
    CNT = SIGFIG
    DO
      IF(SUMR(1) < 0.5) THEN
        MX1=SUMI(L+1)
      ELSEIF(SUMI(1) < 0.5) THEN
        MX1=SUMR(L+1)
      ELSE
        MX1 = MAX(SUMR(L+1),SUMI(L+1))
      ENDIF
      IF(NUMR(1) < 0.5) THEN
        MX2=NUMI(L+1)
      ELSEIF(NUMI(1) < 0.5) THEN
        MX2=NUMR(L+1)
      ELSE
        MX2 = MAX(NUMR(L+1),NUMI(L+1))
      ENDIF
      IF(MX1-MX2 >  2._RP .AND. CR > 0._RP .AND. &
        ABS( CMPLX(AR,AI,RP)*CMPLX(XR,XI,RP)/(CMPLX(CR,CI,RP)*CNT)) <= 1._RP) EXIT
      CALL CMPMUL(SUMR,SUMI,CR,CI,QR1,QI1,L,RMAX)
      CALL CMPMUL(SUMR,SUMI,CR2,CI2,QR2,QI2,L,RMAX)
      QR2(L+1)=QR2(L+1)-1
      QI2(L+1)=QI2(L+1)-1
      CALL CMPADD(QR1,QI1,QR2,QI2,SUMR,SUMI,L,RMAX)

      CALL ARMULT(SUMR,CNT,SUM_tmp,L,RMAX)
      SUMR = SUM_tmp
      CALL ARMULT(SUMI,CNT,SUM_tmp,L,RMAX)
      SUMI = SUM_tmp
      CALL CMPMUL(DENOMR,DENOMI,CR,CI,QR1,QI1,L,RMAX)
      CALL CMPMUL(DENOMR,DENOMI,CR2,CI2,QR2,QI2,L,RMAX)
      QR2(L+1)=QR2(L+1)-1
      QI2(L+1)=QI2(L+1)-1
      CALL CMPADD(QR1,QI1,QR2,QI2,DENOMR,DENOMI,L,RMAX)

      CALL ARMULT(DENOMR,CNT,SUM_tmp,L,RMAX)
      DENOMR = SUM_tmp
      CALL ARMULT(DENOMI,CNT,SUM_tmp,L,RMAX)
      DENOMI = SUM_tmp
      CALL CMPMUL(NUMR,NUMI,AR,AI,QR1,QI1,L,RMAX)
      CALL CMPMUL(NUMR,NUMI,AR2,AI2,QR2,QI2,L,RMAX)
      QR2(L+1)=QR2(L+1)-1
      QI2(L+1)=QI2(L+1)-1
      CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,L,RMAX)

      CALL CMPMUL(NUMR,NUMI,XR,XI,QR1,QI1,L,RMAX)
      CALL CMPMUL(NUMR,NUMI,XR2,XI2,QR2,QI2,L,RMAX)
      QR2(L+1)=QR2(L+1)-1
      QI2(L+1)=QI2(L+1)-1
      CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,L,RMAX)

      CALL CMPADD(SUMR,SUMI,NUMR,NUMI,SUM_tmp,SUM_tmp2,L,RMAX)
      SUMR = SUM_tmp
      SUMI = SUM_tmp2
      CNT=CNT+SIGFIG
      AR=AR+SIGFIG
      CR=CR+SIGFIG
    END DO
    CALL ARYDIV(SUMR,SUMI,DENOMR,DENOMI,FINAL,L,LNCHF,RMAX,BIT)
    CHGF=FINAL
    RETURN
  END FUNCTION CHGF

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

  PURE SUBROUTINE ARADD(A,B,C,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:L+1), B(-1:L+1), RMAX
    REAL(RP), INTENT(OUT) :: C(-1:L+1)

    INTEGER EDIFF,I,J
    REAL(RP) :: Z(-1:L+1)

    Z(0:L+1) = 0._RP
    EDIFF = INT(A(L+1)-B(L+1))
    IF( ABS(A(1))<0.5_RP .OR. EDIFF<=-L ) THEN
      C = B
      GOTO 311
    ELSEIF( ABS(B(1)) < 0.5_RP .OR. EDIFF>=L ) THEN
      C = A
      GOTO 311
    ENDIF

    Z(-1) = A(-1)
    IF( ABS(A(-1)-B(-1))>=0.5 ) THEN
      IF( EDIFF>0 ) THEN
        Z(L+1) = A(L+1)
        GOTO 233
      ENDIF
      IF( EDIFF<0 ) THEN
        Z(L+1)=B(L+1)
        Z(-1)=B(-1)
        GOTO 266
      ENDIF
      DO I=1,L
        IF( A(I)>B(I) ) THEN
          Z(L+1) = A(L+1)
          GOTO 233
        ELSEIF( A(I)<B(I) ) THEN
          Z(L+1) = B(L+1)
          Z(-1) = B(-1)
          GOTO 266
        ENDIF
      END DO
      GOTO 300
    ELSEIF( EDIFF==0 ) THEN
      Z(L+1)=A(L+1)
      DO I=L,1,-1
        Z(I)=A(I)+B(I)+Z(I)
        IF(Z(I) >= RMAX) THEN
          Z(I)=Z(I)-RMAX
          Z(I-1)=1._RP
        ENDIF
      END DO
      IF(Z(0) > 0.5) THEN
        DO I=L,1,-1
          Z(I)=Z(I-1)
        END DO
        Z(L+1)=Z(L+1)+1._RP
        Z(0)=0._RP
      ENDIF
      GOTO 300
    ELSEIF( EDIFF>0 ) THEN
      Z(L+1)=A(L+1)
      DO I=L,1+EDIFF,-1
        Z(I)=A(I)+B(I-EDIFF)+Z(I)
        IF(Z(I) >= RMAX) THEN
          Z(I)=Z(I)-RMAX
          Z(I-1)=1._RP
        ENDIF
      END DO
      DO I=EDIFF,1,-1
        Z(I)=A(I)+Z(I)
        IF(Z(I) >= RMAX) THEN
          Z(I)=Z(I)-RMAX
          Z(I-1)=1._RP
        ENDIF
      END DO
      IF(Z(0) > 0.5) THEN
        DO I=L,1,-1
          Z(I)=Z(I-1)
        END DO
        Z(L+1)=Z(L+1)+1
        Z(0)=0._RP
      ENDIF
      GOTO 300
    ELSE
      Z(L+1)=B(L+1)
      DO I=L,1-EDIFF,-1
        Z(I)=A(I+EDIFF)+B(I)+Z(I)
        IF(Z(I) >= RMAX) THEN
          Z(I)=Z(I)-RMAX
          Z(I-1)=1._RP
        ENDIF
      END DO
      DO I=0-EDIFF,1,-1
        Z(I)=B(I)+Z(I)
        IF(Z(I) >= RMAX) THEN
          Z(I)=Z(I)-RMAX
          Z(I-1)=1._RP
        ENDIF
      END DO
      IF(Z(0) > 0.5) THEN
        DO I=L,1,-1
          Z(I)=Z(I-1)
        END DO
        Z(L+1)=Z(L+1)+1._RP
        Z(0)=0._RP
      ENDIF
      GOTO 300
    END IF
    233 CONTINUE
    IF( EDIFF<=0 ) THEN
      DO I=L,1,-1
        Z(I)=A(I)-B(I)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
      GOTO 290
    ELSE
      DO I=L,1+EDIFF,-1
        Z(I)=A(I)-B(I-EDIFF)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
      DO I=EDIFF,1,-1
        Z(I)=A(I)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
      GOTO 290
    END IF

    266 CONTINUE
    IF( EDIFF>=0 ) THEN
      DO I=L,1,-1
        Z(I)=B(I)-A(I)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
    ELSE
      DO I=L,1-EDIFF,-1
        Z(I)=B(I)-A(I+EDIFF)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
      DO I=0-EDIFF,1,-1
        Z(I)=B(I)+Z(I)
        IF(Z(I) < 0._RP) THEN
          Z(I)=Z(I)+RMAX
          Z(I-1)=-1._RP
        ENDIF
      END DO
    END IF
    290 IF(Z(1) > 0.5) GOTO 300
    I=1
    DO
      I=I+1
      IF( Z(I)>=0.5 .OR. I>=L+1 ) EXIT
    END DO
    IF( I==L+1 ) THEN
      Z(-1)=1._RP
      Z(L+1)=0._RP
      GOTO 300
    ENDIF

    DO J=1,L+1-I
      Z(J)=Z(J+I-1)
    END DO
    DO J=L+2-I,L
      Z(J) = 0._RP
    END DO
    Z(L+1)=Z(L+1)-I+1
    300 C = Z
    311 CONTINUE
    IF( C(1)<0.5 ) THEN
      C(-1)=1._RP
      C(L+1)=0._RP
    ENDIF

    RETURN
  END SUBROUTINE ARADD

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
  !     *  Subprograms called: ARADD                                   *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE ARSUB(A,B,C,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:L+1), B(-1:L+1), RMAX
    REAL(RP), INTENT(OUT) :: C(-1:L+1)

    REAL(RP) ::  B2(-1:L+1)

    B2 = B
    B2(-1) = (-1._RP)*B2(-1)
    CALL ARADD(A,B2,C,L,RMAX)

    RETURN
  END SUBROUTINE ARSUB

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

  PURE SUBROUTINE ARMULT(A,B,C,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: A(-1:L+1), B, RMAX
    REAL(RP), INTENT(OUT) :: C(-1:L+1)

    REAL(RP) ::  B2,CARRY, Z(-1:L+1)!,RMAX2
    INTEGER I

    !RMAX2 = 1._RP/RMAX
    Z(-1) = SIGN(1._RP,B)*A(-1)
    B2=ABS(B)
    Z(L+1) = A(L+1)
    Z(0:L) = 0._RP
    IF( B2<=1.0E-10_RP .OR. A(1)<=1.0E-10_RP ) THEN
      Z(-1) = 1._RP
      Z(L+1) = 0._RP
    ELSE
      DO I=L,1,-1
        Z(I) = A(I)*B2 + Z(I)
        IF( Z(I)>=RMAX ) THEN
          CARRY = INT(Z(I)/RMAX)
          Z(I) = Z(I)-CARRY*RMAX
          Z(I-1) = CARRY
        ENDIF
      END DO
      IF( Z(0)>=0.5 ) THEN
        Z(1:L) = Z(0:L-1)
        Z(L+1) = Z(L+1)+1._RP
        Z(0) = 0._RP
      END IF
    END IF
    C = Z
    IF( C(1)<0.5 ) THEN
      C(-1)=1._RP
      C(L+1)=0._RP
    ENDIF

    RETURN
  END SUBROUTINE ARMULT

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
  !     *  Subprograms called: ARADD                                   *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE CMPADD(AR,AI,BR,BI,CR,CI,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: AR(-1:L+1), AI(-1:L+1), BR(-1:L+1), BI(-1:L+1), RMAX
    REAL(RP), INTENT(OUT) :: CR(-1:L+1), CI(-1:L+1)

    CALL ARADD(AR,BR,CR,L,RMAX)
    CALL ARADD(AI,BI,CI,L,RMAX)
    RETURN
  END SUBROUTINE CMPADD

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
  !     *  Subprograms called: ARADD                                   *
  !     *                                                              *
  !     ****************************************************************

  SUBROUTINE CMPSUB(AR,AI,BR,BI,CR,CI,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: AR(-1:L+1), AI(-1:L+1), BR(-1:L+1), BI(-1:L+1), RMAX
    REAL(RP), INTENT(OUT) :: CR(-1:L+1), CI(-1:L+1)

    CALL ARSUB(AR,BR,CR,L,RMAX)
    CALL ARSUB(AI,BI,CI,L,RMAX)
    RETURN
  END SUBROUTINE CMPSUB

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
  !     *  Subprograms called: ARMULT, ARSUB, ARADD                    *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE CMPMUL(AR,AI,BR,BI,CR,CI,L,RMAX)
    INTEGER, INTENT(IN) :: L
    REAL(RP), INTENT(IN) :: AR(-1:L+1), AI(-1:L+1), BR, BI, RMAX
    REAL(RP), INTENT(OUT) :: CR(-1:L+1), CI(-1:L+1)
    REAL(RP) :: D1(-1:777),D2(-1:777)

    CALL ARMULT(AR,BR,D1,L,RMAX)
    CALL ARMULT(AI,BI,D2,L,RMAX)
    CALL ARSUB(D1,D2,CR,L,RMAX)
    CALL ARMULT(AR,BI,D1,L,RMAX)
    CALL ARMULT(AI,BR,D2,L,RMAX)
    CALL ARADD(D1,D2,CI,L,RMAX)

    RETURN
  END SUBROUTINE CMPMUL

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ARYDIV                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Returns the REAL(RP) ::  complex number   *
  !     *    resulting from the division of four arrays, representing  *
  !     *    two complex numbers.  The number returned will be in one  *
  !     *    two different forms.  Either standard scientific or as    *
  !     *    the log of the number.                                    *
  !     *                                                              *
  !     *  Subprograms called: CONV21, CONV12, EADD, ECPDIV, EMULT     *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE ARYDIV(AR,AI,BR,BI,C,L,LNCHF,RMAX,BIT)
    INTEGER, INTENT(IN) :: L, LNCHF, BIT
    REAL(RP), INTENT(IN) :: AR(-1:L+1), AI(-1:L+1), BR(-1:L+1), BI(-1:L+1), RMAX
    COMPLEX(RP), INTENT(OUT) :: C

    INTEGER :: REXP,IR10,II10
    REAL(RP) ::  PHI,N1,N2,N3,E1,E2,E3,RR10,RI10,X
    REAL(RP) ::  X1,X2,DUM1,DUM2
    REAL(RP) :: AE(2,2),BE(2,2),CE(2,2)

    REXP=BIT/2
    X=REXP*(AR(L+1)-2)
    RR10=X*LOG10(2._RP)/LOG10(10._RP)
    IR10=INT(RR10)
    RR10=RR10-IR10
    X=REXP*(AI(L+1)-2)
    RI10=X*LOG10(2._RP)/LOG10(10._RP)
    II10=INT(RI10)
    RI10=RI10-II10
    DUM1=SIGN(AR(1)*RMAX*RMAX+AR(2)*RMAX+AR(3),AR(-1))
    DUM2=SIGN(AI(1)*RMAX*RMAX+AI(2)*RMAX+AI(3),AI(-1))
    DUM1=DUM1*10**RR10
    DUM2=DUM2*10**RI10
    CALL CONV12(CMPLX(DUM1,DUM2,RP),AE)
    AE(1,2)=AE(1,2)+IR10
    AE(2,2)=AE(2,2)+II10
    X=REXP*(BR(L+1)-2)
    RR10=X*LOG10(2._RP)/LOG10(10._RP)
    IR10=INT(RR10)
    RR10=RR10-IR10
    X=REXP*(BI(L+1)-2)
    RI10=X*LOG10(2._RP)/LOG10(10._RP)
    II10=INT(RI10)
    RI10=RI10-II10
    DUM1 = SIGN(BR(1)*RMAX*RMAX+BR(2)*RMAX+BR(3),BR(-1))
    DUM2 = SIGN(BI(1)*RMAX*RMAX+BI(2)*RMAX+BI(3),BI(-1))
    DUM1=DUM1*10**RR10
    DUM2=DUM2*10**RI10
    CALL CONV12(CMPLX(DUM1,DUM2,RP),BE)
    BE(1,2)=BE(1,2)+IR10
    BE(2,2)=BE(2,2)+II10
    CALL ECPDIV(AE,BE,CE)
    IF(LNCHF == 0) THEN
      CALL CONV21(CE,C)
    ELSE
      CALL EMULT(CE(1,1),CE(1,2),CE(1,1),CE(1,2),N1,E1)
      CALL EMULT(CE(2,1),CE(2,2),CE(2,1),CE(2,2),N2,E2)
      CALL EADD(N1,E1,N2,E2,N3,E3)
      N1=CE(1,1)
      E1=CE(1,2)-CE(2,2)
      X2=CE(2,1)
      IF(E1 > 74._RP) THEN
        X1 = HUGE(1._RP) !1.E75_RP
      ELSEIF(E1 < -74._RP) THEN
        X1 = 0._RP
      ELSE
        X1 = N1*(10**E1)
      ENDIF
      PHI=ATAN2(X2,X1)
      C = CMPLX(0.5_RP*(LOG(N3)+E3*LOG(10._RP)),PHI,RP)
    ENDIF
    RETURN
  END SUBROUTINE ARYDIV

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

  PURE SUBROUTINE EMULT(N1,E1,N2,E2,NF,EF)

    REAL(RP), INTENT(IN) ::  N1,E1,N2,E2
    REAL(RP), INTENT(OUT) :: NF,EF

    NF=N1*N2
    EF=E1+E2
    IF(ABS(NF) >= 10._RP) THEN
      NF=NF/10._RP
      EF=EF+1._RP
    ENDIF
    RETURN
  END SUBROUTINE EMULT

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

  PURE SUBROUTINE EDIV(N1,E1,N2,E2,NF,EF)

    REAL(RP), INTENT(IN) ::  N1,E1,N2,E2
    REAL(RP), INTENT(OUT) :: NF,EF

    NF=N1/N2
    EF=E1-E2
    IF((ABS(NF) < 1._RP) .AND. (NF /= 0._RP)) THEN
      NF=NF*10._RP
      EF=EF-1._RP
    ENDIF
    RETURN
  END SUBROUTINE EDIV

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

  PURE SUBROUTINE EADD(N1,E1,N2,E2,NF,EF)

    REAL(RP),INTENT(IN) :: N1,E1,N2,E2
    REAL(RP),INTENT(OUT) :: NF,EF
    REAL(RP) :: EDIFF

    EDIFF=E1-E2
    IF(EDIFF > 36._RP) THEN
      NF=N1
      EF=E1
    ELSEIF( EDIFF<-36._RP ) THEN
      NF=N2
      EF=E2
    ELSE
      NF=N1*(10._RP**EDIFF)+N2
      EF=E2
      DO
        IF( ABS(NF)<10._RP ) EXIT
        NF=NF/10._RP
        EF=EF+1._RP
      END DO
      DO
        IF( ABS(NF)>=1._RP .OR. NF==0._RP ) EXIT
        NF=NF*10._RP
        EF=EF-1._RP
      END DO
    ENDIF

    RETURN
  END SUBROUTINE EADD

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ESUB                              *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Returns the solution to the subtraction of    *
  !     *    two numbers in the form of base and exponent.             *
  !     *                                                              *
  !     *  Subprograms called: EADD                                    *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE ESUB(N1,E1,N2,E2,NF,EF)

    REAL(RP),INTENT(IN) ::  N1,E1,N2,E2
    REAL(RP),INTENT(OUT) :: NF,EF

    CALL EADD(N1,E1,N2*(-1._RP),E2,NF,EF)
    RETURN
  END SUBROUTINE ESUB

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

  PURE SUBROUTINE CONV12(CN,CAE)

    COMPLEX(RP),INTENT(IN) ::  CN
    REAL(RP),INTENT(OUT) ::  CAE(2,2)

    CAE(1,1)=REAL(CN)
    CAE(1,2)=0._RP
    300 IF(ABS(CAE(1,1)) < 10._RP) GOTO 310
    CAE(1,1)=CAE(1,1)/10._RP
    CAE(1,2)=CAE(1,2)+1._RP
    GOTO 300
    310 IF((ABS(CAE(1,1)) >= 1._RP) .OR. (CAE(1,1) == 0._RP)) GOTO 320
    CAE(1,1)=CAE(1,1)*10._RP
    CAE(1,2)=CAE(1,2)-1._RP
    GOTO 310
    320 CAE(2,1)=AIMAG(CN)
    CAE(2,2)=0._RP
    330 IF(ABS(CAE(2,1)) < 10._RP) GOTO 340
    CAE(2,1)=CAE(2,1)/10._RP
    CAE(2,2)=CAE(2,2)+1._RP
    GOTO 330
    340 IF((ABS(CAE(2,1)) >= 1._RP) .OR. (CAE(2,1) == 0._RP)) GOTO 350
    CAE(2,1)=CAE(2,1)*10._RP
    CAE(2,2)=CAE(2,2)-1._RP
    GOTO 340
    350 RETURN
  END SUBROUTINE CONV12

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

  PURE SUBROUTINE CONV21(CAE,CN)
    !USE IEEE_ARITHMETIC ,ONLY: IEEE_IS_FINITE
    REAL(RP), INTENT(IN) ::  CAE(2,2)
    COMPLEX(RP), INTENT(OUT) :: CN

    IF( CAE(1,2) > 75 .OR. CAE(2,2) > 75) THEN
!      CN = CMPLX( CAE(1,1)*10._RP**75, CAE(2,1)*10._RP**75, RP)
      CN = CMPLX( CAE(1,1)*HUGE(1._RP), CAE(2,1)*HUGE(1._RP), RP)
    ELSEIF(CAE(2,2) < -75) THEN
      CN = CMPLX( CAE(1,1)*(10**CAE(1,2)), 0._RP, RP)
    ELSE
      CN = CMPLX( CAE(1,1)*(10**CAE(1,2)), CAE(2,1)*(10**CAE(2,2)), RP)
    ENDIF

    RETURN
  END SUBROUTINE CONV21

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ECPMUL                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Multiplies two numbers which are each         *
  !     *    represented in the form of a two by two array and returns *
  !     *    the solution in the same form.                            *
  !     *                                                              *
  !     *  Subprograms called: EMULT, ESUB, EADD                       *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE ECPMUL(A,B,C)
    REAL(RP),INTENT(IN) ::  A(2,2),B(2,2)
    REAL(RP),INTENT(OUT) ::  C(2,2)

    REAL(RP) ::  N1,E1,N2,E2
    REAL(RP) :: C2(2,2)

    CALL EMULT(A(1,1),A(1,2),B(1,1),B(1,2),N1,E1)
    CALL EMULT(A(2,1),A(2,2),B(2,1),B(2,2),N2,E2)
    CALL ESUB(N1,E1,N2,E2,C2(1,1),C2(1,2))
    CALL EMULT(A(1,1),A(1,2),B(2,1),B(2,2),N1,E1)
    CALL EMULT(A(2,1),A(2,2),B(1,1),B(1,2),N2,E2)
    CALL EADD(N1,E1,N2,E2,C(2,1),C(2,2))
    C(1,1) = C2(1,1)
    C(1,2) = C2(1,2)
    RETURN
  END SUBROUTINE ECPMUL

  !     ****************************************************************
  !     *                                                              *
  !     *                 SUBROUTINE ECPDIV                            *
  !     *                                                              *
  !     *                                                              *
  !     *  Description : Divides two numbers and returns the solution. *
  !     *    All numbers are represented by a 2x2 array.               *
  !     *                                                              *
  !     *  Subprograms called: EADD, ECPMUL, EDIV, EMULT               *
  !     *                                                              *
  !     ****************************************************************

  PURE SUBROUTINE ECPDIV(A,B,C)
    REAL(RP),INTENT(IN) ::  A(2,2),B(2,2)
    REAL(RP),INTENT(OUT) ::  C(2,2)

    REAL(RP) ::  N1,E1,N2,E2,N3,E3
    REAL(RP) ::  B2(2,2),C2(2,2)

    B2(1,1) = B(1,1)
    B2(1,2) = B(1,2)
    B2(2,1) = -1._RP*B(2,1)
    B2(2,2) = B(2,2)
    CALL ECPMUL(A,B2,C2)
    CALL EMULT(B(1,1),B(1,2),B(1,1),B(1,2),N1,E1)
    CALL EMULT(B(2,1),B(2,2),B(2,1),B(2,2),N2,E2)
    CALL EADD(N1,E1,N2,E2,N3,E3)
    CALL EDIV(C2(1,1),C2(1,2),N3,E3,C(1,1),C(1,2))
    CALL EDIV(C2(2,1),C2(2,2),N3,E3,C(2,1),C(2,2))
    RETURN
  END SUBROUTINE ECPDIV

END MODULE conhyp_m
