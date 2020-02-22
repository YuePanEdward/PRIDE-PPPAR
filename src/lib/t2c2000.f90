SUBROUTINE T2C2000(RPOM, THETA, RBPN, RMAT, RT2C)
!+
!  - - - - - - - -
!   T 2 C 2 0 0 0
!  - - - - - - - -
!
!  Form the TRS-to-CRS matrix, IAU 2000, from components.
!
!  Annex to IERS Conventions 2000, Chapter 5
!
!  Given:
!     RPOM     d(3,3)   polar motion matrix (W)
!
!  followed by either (for the CEO-based transformation):
!     THETA      d      Earth Rotation Angle (radians, giving matrix R)
!     RBPN     d(3,3)   intermediate-to-celestial matrix (Q)
!
!  or alternatively (for the classical, equinox-based, transformation):
!     THETA      d      Greenwich Sidereal Time (radians)
!     RBPN     d(3,3)   true-to-celestial matrix
!
!     RMAT     d(3,3)   rate matrix
!
!  Returned:
!     RT2C     d(3,3)   terrestrial-to-celestial matrix
!
!  Calls the SOFA routines iau_CR, iau_RZ, iau_RXR.
!
!  This revision:  2002 November 25
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INCLUDE '../header/const.h'

  REAL*8 RPOM(3, 3), THETA, RBPN(3, 3), RMAT(3, 3), RT2C(3, 3)
!
!! LOCAL
  REAL*8 SG, CG, R(3, 3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Rate matrix
  RMAT = 0.d0
  SG = DSIN(THETA)
  CG = DCOS(THETA)
  RMAT(1, 1) = -SG*OMEGA
  RMAT(2, 1) = CG*OMEGA
  RMAT(1, 2) = -CG*OMEGA
  RMAT(2, 2) = -SG*OMEGA

!  Polar motion.
  R = RPOM
  CALL MATMPY(RMAT, R, RMAT, 3, 3, 3)

!  Earth rotation.
  CALL ROT_Z(-THETA, R)

!  CIP motion.
  CALL MATMPY(RBPN, R, RT2C, 3, 3, 3)
  CALL MATMPY(RBPN, RMAT, RMAT, 3, 3, 3)

  RETURN
END
