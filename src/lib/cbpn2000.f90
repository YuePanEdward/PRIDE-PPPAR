SUBROUTINE CBPN2000(DATE1, DATE2, DPSI, DEPS, RBPNC)
!+
!  - - - - - - - - -
!   C B P N 2 0 0 0
!  - - - - - - - - -
!
!  Classical bias-precession-nutation matrix.
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Given:
!     DATE1,DATE2   d      TT date (JD = DATE1+DATE2)
!     DPSI,DEPS     d      nutation (luni-solar + planetary, radians)
!
!  Returned:
!     RBPNC       d(3,3)   true-to-celestial matrix
!
!  This revision:  2002 November 26
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 DATE1, DATE2, DPSI, DEPS, RBPNC(3, 3)

!  Arcseconds to radians
  REAL*8, PARAMETER :: DAS2R = 4.848136811095359935899141D-6

!  Reference epoch (J2000), JD
  REAL*8, PARAMETER :: DJ0 = 2451545D0

!  Days per Julian century
  REAL*8, PARAMETER :: DJC = 36525D0

!  J2000 obliquity (Lieske et al. 1977)
  REAL*8, PARAMETER :: EPS0 = 84381.448D0*DAS2R

!  The ICRS RA of the J2000 equinox (Chapront et al., 2002)
  REAL*8, PARAMETER :: DRA0 = -0.0146D0*DAS2R

  REAL*8 T

!  The precession and obliquity corrections (radians per century)
  REAL*8, PARAMETER :: PRECOR = -0.29965D0*DAS2R, OBLCOR = -0.02524D0*DAS2R

!  The frame bias corrections in longitude and obliquity
  REAL*8, PARAMETER :: DPSIBI = -0.041775D0*DAS2R, DEPSBI = -0.0068192D0*DAS2R

  REAL*8 DPSIPR, DEPSPR, EPSA80, PSIA77, OMA77, CHIA, PSIA, OMA, EPSA

  INTEGER*4 I
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and given date (JC).
  T = ((DATE1 - DJ0) + DATE2)/DJC

!  Precession rate contributions with respect to IAU 1976/80.
  DPSIPR = PRECOR*T
  DEPSPR = OBLCOR*T

!  IAU 1980 mean obliquity of date.
  EPSA80 = EPS0 + (-46.8150D0 + (-0.00059D0 + (0.001813D0)*T)*T)*T*DAS2R

!  Precession angles (Lieske et al. 1977)
  PSIA77 = (5038.7784D0 + (-1.07259D0 + (-0.001147D0)*T)*T)*T*DAS2R
  OMA77 = EPS0 + ((0.05127D0 + (-0.007726D0)*T)*T)*T*DAS2R
  CHIA = (10.5526D0 + (-2.38064D0 + (-0.001125D0)*T)*T)*T*DAS2R

!  Apply IAU 2000A precession corrections.
  PSIA = PSIA77 + DPSIPR
  OMA = OMA77 + DEPSPR
  EPSA = EPSA80 + DEPSPR

!  Initialize the true-to-celestial matrix.
  RBPNC = 0.D0
  DO I = 1, 3
    RBPNC(I, I) = 1.D0
  ENDDO

!  Remove IAU 2000A nutation (pure: luni-solar and planetary).
  CALL ROT_X(EPSA + DEPS, RBPNC)
  CALL ROT_Z(DPSI, RBPNC)
  CALL ROT_X(-EPSA, RBPNC)

!  Remove precession (Lieske et al. 1977 plus corrections).
  CALL ROT_Z(-CHIA, RBPNC)
  CALL ROT_X(OMA, RBPNC)
  CALL ROT_Z(PSIA, RBPNC)
  CALL ROT_X(-EPS0, RBPNC)

!  Remove frame bias.
  CALL ROT_X(DEPSBI, RBPNC)
  CALL ROT_Y(-DPSIBI*DSIN(EPS0), RBPNC)
  CALL ROT_Z(-DRA0, RBPNC)

  RETURN
END
