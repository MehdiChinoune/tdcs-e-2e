MODULE vectors
  USE constants ,ONLY: RP
  IMPLICIT NONE

  TYPE vec_c
    REAL(RP) :: x, y, z
  CONTAINS
    PROCEDURE :: norm
  END TYPE vec_c

  TYPE vec_s
    REAL(RP) :: r, theta, phi
  END TYPE vec_s

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE :: add_vec_c, add_vec_s
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE :: sub_vec_c, sub_vec_s
    MODULE PROCEDURE :: neg_vec_c, neg_vec_s
  END INTERFACE

CONTAINS

  ELEMENTAL REAL(RP) FUNCTION norm(cartez)
    CLASS(vec_c), INTENT(IN) :: cartez

    norm = SQRT( cartez%x**2 +cartez%y**2 +cartez%z**2 )

  END FUNCTION norm

  ELEMENTAL TYPE(vec_c) FUNCTION cartez(spher)
    TYPE(vec_s), INTENT(IN) :: spher

    cartez%x = spher%r*SIN(spher%theta)*COS(spher%phi)
    cartez%y = spher%r*SIN(spher%theta)*SIN(spher%phi)
    cartez%z = spher%r*COS(spher%theta)

  END FUNCTION cartez

  ELEMENTAL TYPE(vec_c) FUNCTION add_vec_c(vec1, vec2)
    TYPE(vec_c), INTENT(IN) :: vec1, vec2

    add_vec_c%x = vec1%x +vec2%x
    add_vec_c%y = vec1%y +vec2%y
    add_vec_c%z = vec1%z +vec2%z

  END FUNCTION add_vec_c

  ELEMENTAL TYPE(vec_c) FUNCTION sub_vec_c(vec1, vec2)
    TYPE(vec_c), INTENT(IN) :: vec1, vec2

    sub_vec_c%x = vec1%x -vec2%x
    sub_vec_c%y = vec1%y -vec2%y
    sub_vec_c%z = vec1%z -vec2%z

  END FUNCTION sub_vec_c

  ELEMENTAL TYPE(vec_c) FUNCTION neg_vec_c(vec)
    TYPE(vec_c), INTENT(IN) :: vec

    neg_vec_c%x = -vec%x
    neg_vec_c%y = -vec%y
    neg_vec_c%z = -vec%z

  END FUNCTION neg_vec_c

  ELEMENTAL TYPE(vec_s) FUNCTION spher(cartez)
    TYPE(vec_c), INTENT(IN) :: cartez

    spher%r = cartez%norm()
    spher%theta = ACOS( cartez%z/spher%r )
    spher%phi = ATAN2( cartez%y, cartez%x )

  END FUNCTION spher

  ELEMENTAL TYPE(vec_s) FUNCTION add_vec_s(vec1, vec2)
    TYPE(vec_s), INTENT(IN) :: vec1, vec2

    add_vec_s = spher( cartez(vec1) -cartez(vec2) )

  END FUNCTION add_vec_s

  ELEMENTAL TYPE(vec_s) FUNCTION sub_vec_s(vec1, vec2)
    TYPE(vec_s), INTENT(IN) :: vec1, vec2

    sub_vec_s = spher( cartez(vec1) -cartez(vec2) )

  END FUNCTION sub_vec_s

  ELEMENTAL TYPE(vec_s) FUNCTION neg_vec_s(vec)
    TYPE(vec_s), INTENT(IN) :: vec

    neg_vec_s = spher(-cartez(vec))

  END FUNCTION neg_vec_s

END MODULE vectors
