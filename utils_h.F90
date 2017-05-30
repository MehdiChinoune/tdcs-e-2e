MODULE utils
  USE constants ,ONLY : RP
  IMPLICIT NONE
! Common Variables
  !REAL(KIND=RP) ,PROTECTED :: fac(0:34) ,lnfac(0:400)
  REAL(KIND=RP) :: fac(0:34) ,lnfac(0:400)

INTERFACE

  MODULE SUBROUTINE factorial()
  END SUBROUTINE factorial

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION norm_fac(e,n)
    REAL(KIND=RP) ,INTENT(IN) :: e
    INTEGER       ,INTENT(IN) :: n
  END FUNCTION norm_fac

  MODULE SUBROUTINE read_input(in_unit ,Ei ,Es ,Ee ,thetas ,step ,Atom ,Orbit)
    INTEGER ,INTENT(IN) :: in_unit
    REAL(KIND=RP)    ,INTENT(OUT) :: Ei ,Es ,Ee ,thetas
    INTEGER          ,INTENT(OUT) :: step(3)
    CHARACTER(LEN=2) ,INTENT(OUT) :: Atom ,Orbit
  END SUBROUTINE read_input

  MODULE SUBROUTINE read_orbit(orbit_file ,lo ,no ,n ,a ,e )
    CHARACTER(LEN=5) ,INTENT(IN)  :: orbit_file
    INTEGER          ,INTENT(OUT) :: lo ,no
    INTEGER ,ALLOCATABLE ,INTENT(OUT) :: n(:)
    REAL(KIND=RP) ,ALLOCATABLE ,INTENT(OUT) :: a(:) ,e(:)
  END SUBROUTINE read_orbit

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION y1y2y3(l1,l2,l3,m1,m2,m3)
    INTEGER,INTENT(IN) :: l1,l2,l3,m1,m2,m3
  END FUNCTION y1y2y3

  MODULE SUBROUTINE clenshaw_curtis( a ,b ,x ,w ,n )
    REAL(KIND=RP) ,INTENT(IN) :: a ,b
    INTEGER ,INTENT(IN) :: n
    REAL(KIND=RP) ,INTENT(OUT) ,ALLOCATABLE :: w(:) ,x(:)
  END SUBROUTINE clenshaw_curtis

  MODULE SUBROUTINE ode_second_dw(km,lmax,rc,f,s,delta)
    INTEGER ,INTENT(IN) :: lmax
    REAL(KIND=RP) ,INTENT(IN)  :: f(0:,0:),rc,km
    REAL(KIND=RP) ,INTENT(OUT) :: s(0:,0:),delta(0:lmax)
  END SUBROUTINE ode_second_dw

  ELEMENTAL REAL(KIND=RP) MODULE FUNCTION Uij(ni ,ei ,nj ,ej ,r)
    REAL(KIND=RP) ,INTENT(IN) :: ei ,ej ,r
    INTEGER       ,INTENT(IN) :: ni ,nj
  END FUNCTION Uij

  MODULE SUBROUTINE calculate_U(Atom ,Orbit ,r ,U )
    CHARACTER(LEN=2) ,INTENT(IN) :: Atom ,Orbit
    REAL(KIND=RP)    ,INTENT(IN) :: r(:)
    REAL(KIND=RP)    ,INTENT(OUT) :: U(:)
  END SUBROUTINE calculate_U

END INTERFACE

END MODULE utils
