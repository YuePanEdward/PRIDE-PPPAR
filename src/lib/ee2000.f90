REAL*8 FUNCTION EE2000(DATE1, DATE2, DPSI)
!+
!  - - - - - - -
!   E E 2 0 0 0
!  - - - - - - -
!
!  The equation of the equinoxes, compatible with IAU 2000 resolutions,
!  given the nutation in longitude.
!
!  Annex to IERS Conventions 2000, Chapter 5
!
!  Given:
!     DATE1,DATE2    d      TT date (JD = DATE1+DATE2)
!     DPSI           d      nutation in longitude (radians)
!
!  Returned:
!     EE2000         d      equation of the equinoxes
!
!  Calls the IERS routine EECT2000.
!
!  This revision:  2002 November 26
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 DATE1, DATE2, DPSI

!  Arcseconds to radians
  REAL*8, PARAMETER :: DAS2R = 4.848136811095359935899141D-6

!  Reference epoch (J2000), JD
  REAL*8, PARAMETER :: DJ0 = 2451545D0

!  Days per Julian century
  REAL*8, PARAMETER :: DJC = 36525D0

!  J2000 obliquity (Lieske et al. 1977)
  REAL*8, PARAMETER :: EPS0 = 84381.448D0*DAS2R

  REAL*8 T, EPSA

  REAL*8 EECT2000

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and given date (JC).
  T = ((DATE1 - DJ0) + DATE2)/DJC

!  Mean obliquity from Chapter 5, expression (32).
  EPSA = EPS0 + (-46.8402D0 + (-0.00059D0 + (0.001813D0)*T)*T)*T*DAS2R

!  Equation of the equinoxes.
  EE2000 = DPSI*COS(EPSA) + EECT2000(DATE1, DATE2)

  RETURN
END
