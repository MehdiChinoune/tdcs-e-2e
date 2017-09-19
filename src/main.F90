PROGRAM main
  use iso_fortran_env  ,only: real_kinds ,real64
  USE constants        , ONLY : RP ,DP ,pi
  USE trigo            , ONLY : spher2cartez
  USE special_functions, ONLY : assoc_legendre
  USE fdcs_e2e         , ONLY : fdcs_fba_dw, fdcs_fba_cw, fdcs_fba_pw, fdcs_dwb

  IMPLICIT NONE
  integer(KIND=selected_int_kind(6)) :: start, finish
  REAL(KIND=RP) :: rate
  integer :: in_unit, out_unit, narg
  character(len=4) :: arg1

  real(rp) :: fak(1:500)

  open( newunit=in_unit, file='input.dat', status='old', action='read')
  open( newunit=out_unit, file='output.dat', status='replace', action='write')

  call SYSTEM_CLOCK(start,rate)

  narg = COMMAND_ARGUMENT_COUNT()
  if(narg>0) then

    call get_command_argument(1, arg1 )

    if(trim(arg1)=='dw') then
      call fdcs_fba_dw(in_unit,out_unit)
    elseif(trim(arg1)=='cw') then
      call fdcs_fba_cw(in_unit,out_unit)
    elseif(trim(arg1)=='pw') then
      call fdcs_fba_pw(in_unit,out_unit)
    end if

  else

    block

       call fdcs_fba_cw(in_unit,out_unit)
       rewind in_unit
       call fdcs_dwb(in_unit,out_unit)

    end block
  endif

  call SYSTEM_CLOCK(finish)

  close(in_unit)
  close(out_unit)

  print*,(finish-start)/rate

CONTAINS

  SUBROUTINE FAKRED()
    REAL(rp)    :: ZI
    INTEGER :: I
!  TABLE OF LOGS OF N!/(20)**N
    FAK(1)=0._RP
    DO I=2,500
      ZI=REAL(I-1, rp)
      FAK(I)=FAK(I-1)+LOG(ZI)+LOG(0.05_RP)
    END DO
!
    RETURN
  END SUBROUTINE FAKRED

  PURE SUBROUTINE THREEJ(L1,L2,L3,M1,M2,M3,SU)
    INTEGER ,INTENT(IN) :: L1,L2,L3,M1,M2,M3
    REAL(RP) ,INTENT(OUT) :: SU
    REAL(RP) :: DEL,SYM
    INTEGER :: K,K1,K2,K3,K4,K5,K6,K7,J01,J02,J03,J04,MO,MU,M,J1,J2,J3,J4,J5,J6,J7,J8
!    REAL(RP) :: FAX(400)
!    FAK = LNFAC
    ! COMMON /REDFAK/ FAK(500)
!
!  CALCULATES THE WIGNER 3J COEFFICIENT
!    (J1, J2, J3)
!    (M1, M2, M3)
    SU=0._RP
    IF( ABS(M1)>L1 .OR. ABS(M2)>L2 .OR. ABS(M3)>L3  .OR. (L3-1)>(L1+L2) .OR. (L1-1)>(L2+L3) &
      .OR. (L2-1)>(L1+L3)  ) RETURN

    K1 = L3 -L2 +M1
    K2 = L3 -L1 -M2
    K3 = L1 +L2 -L3
    K4 = L1 -M1
    K5 = L2 +M2
    K6 = L1 +M1
    K7 = L3 -M3
    J01 = K3 +1
    J02 = K2 + K5 +1
    J03 = K1 + K4 +1
    J04 = J02 + K4 + K6 +1
    MO = MIN(K3,K4,K5) +1

    MU = MIN(K1,K2)
    M = 1
    IF( MU<0 ) M = 1 - MU
    IF(M>MO) RETURN

    DEL=LOG(0.05)+FAK(J01)+FAK(J02)+FAK(J03)-FAK(J04)
    J6 = K2 + K3 +1
    J7 = K2 + K4 +1
    J8 = K1 + K5 +1
    SYM=(DEL+FAK(K6+1)+FAK(K4+1)+FAK(K5+1)+FAK(J6)+FAK(J7) +FAK(J8))/2.
    DO K = M,MO
      J1 = K1 + K
      J2 = K2 + K
      J3 = K3 - K +2
      J4 = K4 - K +2
      J5 = K5 - K +2
      SU=EXP(SYM-FAK(K)-FAK(J1)-FAK(J2)-FAK(J3)-FAK(J4)-FAK(J5))-SU
    END DO
    SU = (-1)**(L3-M3+MO+1)*SU
    IF( ABS(SU)<epsilon(1._rp) ) SU = 0.E0
    RETURN
  END SUBROUTINE THREEJ

END PROGRAM main
