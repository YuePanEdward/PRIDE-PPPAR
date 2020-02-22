REAL*8 FUNCTION EECT2000(DATE1, DATE2)
!+
!  - - - - - - - - -
!   E E C T 2 0 0 0
!  - - - - - - - - -
!
!  Equation of the equinoxes complementary terms, consistent with
!  IAU 2000 resolutions.
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Given:
!     DATE1,DATE2   d    TT date (JD = DATE1+DATE2)
!
!  Returned:
!     EECT2000      d    complementary terms (radians)
!
!  This revision:  2003 May 4
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 DATE1, DATE2

!  2Pi
  REAL*8, PARAMETER :: D2PI = 6.283185307179586476925287D0

!  Arcseconds to radians
  REAL*8, PARAMETER :: DAS2R = 4.848136811095359935899141D-6

!  Reference epoch (J2000), JD
  REAL*8, PARAMETER :: DJ0 = 2451545D0

!  Days per Julian century
  REAL*8, PARAMETER :: DJC = 36525D0

!  Time since J2000, in Julian centuries
  REAL*8 T

!  Miscellaneous
  INTEGER*4 I, J
  REAL*8 A, S0, S1

!  Fundamental arguments
  REAL*8 FA(14)

!  -----------------------------------------
!  The series for the EE complementary terms
!  -----------------------------------------

!  Number of terms in the series
  INTEGER*4, PARAMETER :: NE0 = 33, NE1 = 1

!  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
  INTEGER*4 KE0(14, NE0), KE1(14, NE1)

!  Sine and cosine coefficients
  REAL*8 SE0(2, NE0), SE1(2, NE1)
!
!! FUNCTION CALLED
  REAL*8 ANPM

!  Argument coefficients for t^0
  DATA((KE0(I, J), I=1, 14), J=1, 10)/ &
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  DATA((KE0(I, J), I=1, 14), J=11, 20)/ &
    1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  DATA((KE0(I, J), I=1, 14), J=21, 30)/ &
    0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, -2, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, -2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -1, &
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    2, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 1, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 4, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  DATA((KE0(I, J), I=1, 14), J=31, NE0)/ &
    0, 0, 2, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, -2, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0/

!  Argument coefficients for t^1
  DATA((KE1(I, J), I=1, 14), J=1, NE1)/ &
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/

!  Sine and cosine coefficients for t^0
  DATA((SE0(I, J), I=1, 2), J=1, 10)/ &
    +2640.96D-6, -0.39D-6, &
    +63.52D-6, -0.02D-6, &
    +11.75D-6, +0.01D-6, &
    +11.21D-6, +0.01D-6, &
    -4.55D-6, +0.00D-6, &
    +2.02D-6, +0.00D-6, &
    +1.98D-6, +0.00D-6, &
    -1.72D-6, +0.00D-6, &
    -1.41D-6, -0.01D-6, &
    -1.26D-6, -0.01D-6/
  DATA((SE0(I, J), I=1, 2), J=11, 20)/ &
    -0.63D-6, +0.00D-6, &
    -0.63D-6, +0.00D-6, &
    +0.46D-6, +0.00D-6, &
    +0.45D-6, +0.00D-6, &
    +0.36D-6, +0.00D-6, &
    -0.24D-6, -0.12D-6, &
    +0.32D-6, +0.00D-6, &
    +0.28D-6, +0.00D-6, &
    +0.27D-6, +0.00D-6, &
    +0.26D-6, +0.00D-6/
  DATA((SE0(I, J), I=1, 2), J=21, 30)/ &
    -0.21D-6, +0.00D-6, &
    +0.19D-6, +0.00D-6, &
    +0.18D-6, +0.00D-6, &
    -0.10D-6, +0.05D-6, &
    +0.15D-6, +0.00D-6, &
    -0.14D-6, +0.00D-6, &
    +0.14D-6, +0.00D-6, &
    -0.14D-6, +0.00D-6, &
    +0.14D-6, +0.00D-6, &
    +0.13D-6, +0.00D-6/
  DATA((SE0(I, J), I=1, 2), J=31, NE0)/ &
    -0.11D-6, +0.00D-6, &
    +0.11D-6, +0.00D-6, &
    +0.11D-6, +0.00D-6/

!  Sine and cosine coefficients for t^1
  DATA((SE1(I, J), I=1, 2), J=1, NE1)/ &
    -0.87D-6, +0.00D-6/

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and current date (JC).
  T = ((DATE1 - DJ0) + DATE2)/DJC

!  Fundamental Arguments (from IERS Conventions 2000)

!  Mean Anomaly of the Moon.
  FA(1) = ANPM((485868.249036D0 + &
                (715923.2178D0 + &
                 (31.8792D0 + &
                  (0.051635D0 + &
                   (-0.00024470D0) &
                   *T)*T)*T)*T)*DAS2R &
               + MOD(1325D0*T, 1D0)*D2PI)

!  Mean Anomaly of the Sun.
  FA(2) = ANPM((1287104.793048D0 + &
                (1292581.0481D0 + &
                 (-0.5532D0 + &
                  (+0.000136D0 + &
                   (-0.00001149D0) &
                   *T)*T)*T)*T)*DAS2R &
               + MOD(99D0*T, 1D0)*D2PI)

!  Mean Longitude of the Moon minus Mean Longitude of the Ascending
!  Node of the Moon.
  FA(3) = ANPM((335779.526232D0 + &
                (295262.8478D0 + &
                 (-12.7512D0 + &
                  (-0.001037D0 + &
                   (0.00000417D0) &
                   *T)*T)*T)*T)*DAS2R &
               + MOD(1342D0*T, 1D0)*D2PI)

!  Mean Elongation of the Moon from the Sun.
  FA(4) = ANPM((1072260.703692D0 + &
                (1105601.2090D0 + &
                 (-6.3706D0 + &
                  (0.006593D0 + &
                   (-0.00003169D0) &
                   *T)*T)*T)*T)*DAS2R &
               + MOD(1236D0*T, 1D0)*D2PI)

!  Mean Longitude of the Ascending Node of the Moon.
  FA(5) = ANPM((450160.398036D0 + &
                (-482890.5431D0 + &
                 (7.4722D0 + &
                  (0.007702D0 + &
                   (-0.00005939D0) &
                   *T)*T)*T)*T)*DAS2R &
               + MOD(-5D0*T, 1D0)*D2PI)

  FA(6) = ANPM(4.402608842D0 + 2608.7903141574D0*T)
  FA(7) = ANPM(3.176146697D0 + 1021.3285546211D0*T)
  FA(8) = ANPM(1.753470314D0 + 628.3075849991D0*T)
  FA(9) = ANPM(6.203480913D0 + 334.0612426700D0*T)
  FA(10) = ANPM(0.599546497D0 + 52.9690962641D0*T)
  FA(11) = ANPM(0.874016757D0 + 21.3299104960D0*T)
  FA(12) = ANPM(5.481293872D0 + 7.4781598567D0*T)
  FA(13) = ANPM(5.311886287D0 + 3.8133035638D0*T)
  FA(14) = (0.024381750D0 + 0.00000538691D0*T)*T

!  Evaluate the EE complementary terms.
  S0 = 0D0
  S1 = 0D0

  DO I = NE0, 1, -1
    A = 0D0
    DO J = 1, 14
      A = A + DBLE(KE0(J, I))*FA(J)
    END DO
    S0 = S0 + (SE0(1, I)*SIN(A) + SE0(2, I)*COS(A))
  END DO
  DO I = NE1, 1, -1
    A = 0D0
    DO J = 1, 14
      A = A + DBLE(KE1(J, I))*FA(J)
    END DO
    S1 = S1 + (SE1(1, I)*SIN(A) + SE1(2, I)*COS(A))
  END DO
  EECT2000 = (S0 + S1*T)*DAS2R

!  Finished.
  RETURN
END

REAL*8 FUNCTION ANPM(A)
  IMPLICIT NONE

  REAL*8 A
!
!! LOCAL
  REAL*8, PARAMETER :: DPI = 3.141592653589793238462643D0
  REAL*8, PARAMETER :: D2PI = 6.283185307179586476925287D0
  REAL*8 W

  W = MOD(A, D2PI)
  IF (DABS(W) .GE. DPI) W = W - SIGN(D2PI, A)
  ANPM = W

  RETURN
END
