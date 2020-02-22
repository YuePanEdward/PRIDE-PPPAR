!
REAL*8 FUNCTION GST2000(UTA, UTB, TTA, TTB, DPSI)
!+
!  - - - - - - - -
!   G S T 2 0 0 0
!  - - - - - - - -
!
!  Greenwich Sidereal Time (model consistent with IAU 2000 resolutions).
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Given:
!     UTA, UTB     d      UT1 date (JD = UTA+UTB)
!     TTA, TTB     d      TT date (JD = TTA+TTB)
!     DPSI         d      nutation in longitude (radians)
!
!  The result is the Greenwich Mean (apparent) Sidereal Time (radians),
!  in the range 0 to 2pi.
!
!  This revision:  2002 December 9
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 UTA, UTB, TTA, TTB, DPSI
!
!! LOCAL
  REAL*8, PARAMETER :: D2PI = 6.283185307179586476925287D0
!
!! FUNCTION CALLED
  REAL*8 GMST2000, EE2000

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Greenwich Sidereal Time, IAU 2000.
  GST2000 = GMST2000(UTA, UTB, TTA, TTB) + EE2000(TTA, TTB, DPSI)
  GST2000 = MOD(GST2000, D2PI)
  IF (GST2000 .LT. 0.D0) GST2000 = GST2000 + D2PI

  RETURN
END
