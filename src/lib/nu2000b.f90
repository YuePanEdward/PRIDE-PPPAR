!
SUBROUTINE NU2000B(DATE1, DATE2, DPSI, DEPS)
!+
!
!  Nutation, IAU 2000B (truncated) model.
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Given:
!     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
!
!  Returned:
!     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
!
!  This revision:  2002 November 25
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 DATE1, DATE2, DPSI, DEPS

!  Arcseconds to radians
  REAL*8, PARAMETER :: DAS2R = 4.848136811095359935899141D-6

!  Milliarcseconds to radians
  REAL*8, PARAMETER :: DMAS2R = DAS2R/1D3

!  Arc seconds in a full circle
  REAL*8, PARAMETER :: TURNAS = 1296000D0

!  2Pi
  REAL*8, PARAMETER :: D2PI = 6.283185307179586476925287D0

!  Units of 0.1 microarcsecond to radians
  REAL*8, PARAMETER :: U2R = DAS2R/1D7

!  Reference epoch (J2000), JD
  REAL*8, PARAMETER :: DJ0 = 2451545D0

!  Days per Julian century
  REAL*8, PARAMETER :: DJC = 36525D0

!  Miscellaneous
  REAL*8 T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG, DPSILS, DEPSLS, DPSIPL, DEPSPL
  INTEGER*4 I, J

!  -------------------------
!  Luni-Solar nutation model
!  -------------------------

!  Number of terms in the luni-solar nutation model
  INTEGER*4, PARAMETER :: NLS = 77

!  Coefficients for fundamental arguments
  INTEGER*4 NALS(5, NLS)

!  Longitude and obliquity coefficients
  REAL*8 CLS(6, NLS)

!  ------------------
!  Planetary nutation (radians)
!  ------------------

  REAL*8, PARAMETER :: DPPLAN = -0.135D0*DMAS2R, DEPLAN = +0.388D0*DMAS2R

!  n.b.  The above fixed terms account for the omission of the
!  long-period planetary terms in the truncated model.

!  ----------------------------------------
!  Tables of argument and term coefficients
!  ----------------------------------------

!
!  Luni-Solar argument multipliers:
!               L     L'    F     D     Om

  DATA((NALS(I, J), I=1, 5), J=1, 10)/ &
    0, 0, 0, 0, 1, &
    0, 0, 2, -2, 2, &
    0, 0, 2, 0, 2, &
    0, 0, 0, 0, 2, &
    0, 1, 0, 0, 0, &
    0, 1, 2, -2, 2, &
    1, 0, 0, 0, 0, &
    0, 0, 2, 0, 1, &
    1, 0, 2, 0, 2, &
    0, -1, 2, -2, 2/
  DATA((NALS(I, J), I=1, 5), J=11, 20)/ &
    0, 0, 2, -2, 1, &
    -1, 0, 2, 0, 2, &
    -1, 0, 0, 2, 0, &
    1, 0, 0, 0, 1, &
    -1, 0, 0, 0, 1, &
    -1, 0, 2, 2, 2, &
    1, 0, 2, 0, 1, &
    -2, 0, 2, 0, 1, &
    0, 0, 0, 2, 0, &
    0, 0, 2, 2, 2/
  DATA((NALS(I, J), I=1, 5), J=21, 30)/ &
    0, -2, 2, -2, 2, &
    -2, 0, 0, 2, 0, &
    2, 0, 2, 0, 2, &
    1, 0, 2, -2, 2, &
    -1, 0, 2, 0, 1, &
    2, 0, 0, 0, 0, &
    0, 0, 2, 0, 0, &
    0, 1, 0, 0, 1, &
    -1, 0, 0, 2, 1, &
    0, 2, 2, -2, 2/
  DATA((NALS(I, J), I=1, 5), J=31, 40)/ &
    0, 0, -2, 2, 0, &
    1, 0, 0, -2, 1, &
    0, -1, 0, 0, 1, &
    -1, 0, 2, 2, 1, &
    0, 2, 0, 0, 0, &
    1, 0, 2, 2, 2, &
    -2, 0, 2, 0, 0, &
    0, 1, 2, 0, 2, &
    0, 0, 2, 2, 1, &
    0, -1, 2, 0, 2/
  DATA((NALS(I, J), I=1, 5), J=41, 50)/ &
    0, 0, 0, 2, 1, &
    1, 0, 2, -2, 1, &
    2, 0, 2, -2, 2, &
    -2, 0, 0, 2, 1, &
    2, 0, 2, 0, 1, &
    0, -1, 2, -2, 1, &
    0, 0, 0, -2, 1, &
    -1, -1, 0, 2, 0, &
    2, 0, 0, -2, 1, &
    1, 0, 0, 2, 0/
  DATA((NALS(I, J), I=1, 5), J=51, 60)/ &
    0, 1, 2, -2, 1, &
    1, -1, 0, 0, 0, &
    -2, 0, 2, 0, 2, &
    3, 0, 2, 0, 2, &
    0, -1, 0, 2, 0, &
    1, -1, 2, 0, 2, &
    0, 0, 0, 1, 0, &
    -1, -1, 2, 2, 2, &
    -1, 0, 2, 0, 0, &
    0, -1, 2, 2, 2/
  DATA((NALS(I, J), I=1, 5), J=61, 70)/ &
    -2, 0, 0, 0, 1, &
    1, 1, 2, 0, 2, &
    2, 0, 0, 0, 1, &
    -1, 1, 0, 1, 0, &
    1, 1, 0, 0, 0, &
    1, 0, 2, 0, 0, &
    -1, 0, 2, -2, 1, &
    1, 0, 0, 0, 2, &
    -1, 0, 0, 1, 0, &
    0, 0, 2, 1, 2/
  DATA((NALS(I, J), I=1, 5), J=71, 77)/ &
    -1, 0, 2, 4, 2, &
    -1, 1, 0, 1, 1, &
    0, -2, 2, -2, 1, &
    1, 0, 2, 2, 1, &
    -2, 0, 2, 2, 2, &
    -1, 0, 0, 0, 2, &
    1, 1, 2, -2, 2/

!
!  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
!  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
!
!  Each row of coefficients in CLS belongs with the corresponding row of
!  fundamental-argument multipliers in NALS.
!

  DATA((CLS(I, J), I=1, 6), J=1, 10)/ &
    -172064161D0, -174666D0, 33386D0, 92052331D0, 9086D0, 15377D0, &
    -13170906D0, -1675D0, -13696D0, 5730336D0, -3015D0, -4587D0, &
    -2276413D0, -234D0, 2796D0, 978459D0, -485D0, 1374D0, &
    2074554D0, 207D0, -698D0, -897492D0, 470D0, -291D0, &
    1475877D0, -3633D0, 11817D0, 73871D0, -184D0, -1924D0, &
    -516821D0, 1226D0, -524D0, 224386D0, -677D0, -174D0, &
    711159D0, 73D0, -872D0, -6750D0, 0D0, 358D0, &
    -387298D0, -367D0, 380D0, 200728D0, 18D0, 318D0, &
    -301461D0, -36D0, 816D0, 129025D0, -63D0, 367D0, &
    215829D0, -494D0, 111D0, -95929D0, 299D0, 132D0/
  DATA((CLS(I, J), I=1, 6), J=11, 20)/ &
    128227D0, 137D0, 181D0, -68982D0, -9D0, 39D0, &
    123457D0, 11D0, 19D0, -53311D0, 32D0, -4D0, &
    156994D0, 10D0, -168D0, -1235D0, 0D0, 82D0, &
    63110D0, 63D0, 27D0, -33228D0, 0D0, -9D0, &
    -57976D0, -63D0, -189D0, 31429D0, 0D0, -75D0, &
    -59641D0, -11D0, 149D0, 25543D0, -11D0, 66D0, &
    -51613D0, -42D0, 129D0, 26366D0, 0D0, 78D0, &
    45893D0, 50D0, 31D0, -24236D0, -10D0, 20D0, &
    63384D0, 11D0, -150D0, -1220D0, 0D0, 29D0, &
    -38571D0, -1D0, 158D0, 16452D0, -11D0, 68D0/
  DATA((CLS(I, J), I=1, 6), J=21, 30)/ &
    32481D0, 0D0, 0D0, -13870D0, 0D0, 0D0, &
    -47722D0, 0D0, -18D0, 477D0, 0D0, -25D0, &
    -31046D0, -1D0, 131D0, 13238D0, -11D0, 59D0, &
    28593D0, 0D0, -1D0, -12338D0, 10D0, -3D0, &
    20441D0, 21D0, 10D0, -10758D0, 0D0, -3D0, &
    29243D0, 0D0, -74D0, -609D0, 0D0, 13D0, &
    25887D0, 0D0, -66D0, -550D0, 0D0, 11D0, &
    -14053D0, -25D0, 79D0, 8551D0, -2D0, -45D0, &
    15164D0, 10D0, 11D0, -8001D0, 0D0, -1D0, &
    -15794D0, 72D0, -16D0, 6850D0, -42D0, -5D0/
  DATA((CLS(I, J), I=1, 6), J=31, 40)/ &
    21783D0, 0D0, 13D0, -167D0, 0D0, 13D0, &
    -12873D0, -10D0, -37D0, 6953D0, 0D0, -14D0, &
    -12654D0, 11D0, 63D0, 6415D0, 0D0, 26D0, &
    -10204D0, 0D0, 25D0, 5222D0, 0D0, 15D0, &
    16707D0, -85D0, -10D0, 168D0, -1D0, 10D0, &
    -7691D0, 0D0, 44D0, 3268D0, 0D0, 19D0, &
    -11024D0, 0D0, -14D0, 104D0, 0D0, 2D0, &
    7566D0, -21D0, -11D0, -3250D0, 0D0, -5D0, &
    -6637D0, -11D0, 25D0, 3353D0, 0D0, 14D0, &
    -7141D0, 21D0, 8D0, 3070D0, 0D0, 4D0/
  DATA((CLS(I, J), I=1, 6), J=41, 50)/ &
    -6302D0, -11D0, 2D0, 3272D0, 0D0, 4D0, &
    5800D0, 10D0, 2D0, -3045D0, 0D0, -1D0, &
    6443D0, 0D0, -7D0, -2768D0, 0D0, -4D0, &
    -5774D0, -11D0, -15D0, 3041D0, 0D0, -5D0, &
    -5350D0, 0D0, 21D0, 2695D0, 0D0, 12D0, &
    -4752D0, -11D0, -3D0, 2719D0, 0D0, -3D0, &
    -4940D0, -11D0, -21D0, 2720D0, 0D0, -9D0, &
    7350D0, 0D0, -8D0, -51D0, 0D0, 4D0, &
    4065D0, 0D0, 6D0, -2206D0, 0D0, 1D0, &
    6579D0, 0D0, -24D0, -199D0, 0D0, 2D0/
  DATA((CLS(I, J), I=1, 6), J=51, 60)/ &
    3579D0, 0D0, 5D0, -1900D0, 0D0, 1D0, &
    4725D0, 0D0, -6D0, -41D0, 0D0, 3D0, &
    -3075D0, 0D0, -2D0, 1313D0, 0D0, -1D0, &
    -2904D0, 0D0, 15D0, 1233D0, 0D0, 7D0, &
    4348D0, 0D0, -10D0, -81D0, 0D0, 2D0, &
    -2878D0, 0D0, 8D0, 1232D0, 0D0, 4D0, &
    -4230D0, 0D0, 5D0, -20D0, 0D0, -2D0, &
    -2819D0, 0D0, 7D0, 1207D0, 0D0, 3D0, &
    -4056D0, 0D0, 5D0, 40D0, 0D0, -2D0, &
    -2647D0, 0D0, 11D0, 1129D0, 0D0, 5D0/
  DATA((CLS(I, J), I=1, 6), J=61, 70)/ &
    -2294D0, 0D0, -10D0, 1266D0, 0D0, -4D0, &
    2481D0, 0D0, -7D0, -1062D0, 0D0, -3D0, &
    2179D0, 0D0, -2D0, -1129D0, 0D0, -2D0, &
    3276D0, 0D0, 1D0, -9D0, 0D0, 0D0, &
    -3389D0, 0D0, 5D0, 35D0, 0D0, -2D0, &
    3339D0, 0D0, -13D0, -107D0, 0D0, 1D0, &
    -1987D0, 0D0, -6D0, 1073D0, 0D0, -2D0, &
    -1981D0, 0D0, 0D0, 854D0, 0D0, 0D0, &
    4026D0, 0D0, -353D0, -553D0, 0D0, -139D0, &
    1660D0, 0D0, -5D0, -710D0, 0D0, -2D0/
  DATA((CLS(I, J), I=1, 6), J=71, 77)/ &
    -1521D0, 0D0, 9D0, 647D0, 0D0, 4D0, &
    1314D0, 0D0, 0D0, -700D0, 0D0, 0D0, &
    -1283D0, 0D0, 0D0, 672D0, 0D0, 0D0, &
    -1331D0, 0D0, 8D0, 663D0, 0D0, 4D0, &
    1383D0, 0D0, -2D0, -594D0, 0D0, -2D0, &
    1405D0, 0D0, 4D0, -610D0, 0D0, 2D0, &
    1290D0, 0D0, 0D0, -556D0, 0D0, 0D0/

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and given date (JC).
  T = ((DATE1 - DJ0) + DATE2)/DJC

!  -------------------
!  LUNI-SOLAR NUTATION
!  -------------------

!
!  Fundamental (Delaunay) arguments from Simon et al. (1994)
!

!  Mean anomaly of the Moon.
  EL = MOD(485868.249036D0 + T*1717915923.2178D0, TURNAS)*DAS2R

!  Mean anomaly of the Sun.
  ELP = MOD(1287104.79305D0 + T*129596581.0481D0, TURNAS)*DAS2R

!  Mean argument of the latitude of the Moon.
  F = MOD(335779.526232D0 + T*1739527262.8478D0, TURNAS)*DAS2R

!  Mean elongation of the Moon from the Sun.
  D = MOD(1072260.70369D0 + T*1602961601.2090D0, TURNAS)*DAS2R

!  Mean longitude of the ascending node of the Moon.
  OM = MOD(450160.398036D0 - T*6962890.5431D0, TURNAS)*DAS2R

!  Initialize the nutation values.
  DP = 0D0
  DE = 0D0

!  Summation of luni-solar nutation series (in reverse order).
  DO I = NLS, 1, -1

!     Argument and functions.
    ARG = MOD(DBLE(NALS(1, I))*EL + &
              DBLE(NALS(2, I))*ELP + &
              DBLE(NALS(3, I))*F + &
              DBLE(NALS(4, I))*D + &
              DBLE(NALS(5, I))*OM, D2PI)
    SARG = SIN(ARG)
    CARG = COS(ARG)

!     Term.
    DP = DP + (CLS(1, I) + CLS(2, I)*T)*SARG &
         + CLS(3, I)*CARG
    DE = DE + (CLS(4, I) + CLS(5, I)*T)*CARG &
         + CLS(6, I)*SARG

  ENDDO

!  Convert from 0.1 microarcsec units to radians.
  DPSILS = DP*U2R
  DEPSLS = DE*U2R

!  ------------------
!  PLANETARY NUTATION
!  ------------------

!  Fixed terms to allow for long-period nutation.
  DPSIPL = DPPLAN
  DEPSPL = DEPLAN

!  -----
!  TOTAL
!  -----

!  Add planetary and luni-solar components.
  DPSI = DPSIPL + DPSILS
  DEPS = DEPSPL + DEPSLS

  RETURN
END
