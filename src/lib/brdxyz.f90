!
!! brdxyz.f90
!!
!!    Copyright (C) 2018 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! author: J. Geng, M. Ge
!! tester: Y. Pan, X. Chen, J. Zhou, S. Mao
!!
!!
!! purpose   : compute position and velocity of GPS satellite using broadcast ephemeris
!!
!! parameters:
!!     cmode  == string with 3 characters for pos,vel and clock, respectively.
!!             for example,   yyy -- all required, nny -- only clock
!!     lstick == use the last one (jeph) if it is not zero, to avoid discontinuity by two brdeph
!!     neph, ephem == total number of brdeph in 'ephem'
!!     iprn == PRN of the satellite its pos and vel to be computed
!!     nwk, sow ==  epoch time, (#week, sow) or (mjd,sod) GPS time
!!     xsat, dtsat == position , velocity and clock correction
!!     v ==
!!     dtref == t - toe/toc
!!
!
subroutine brdxyz(cmode, lstick, neph, ephem, iprn, jeph, nwk, sow, xsat, v, dtsat, dtref)
  implicit none
  include '../header/brdeph.h'

  character*3 cmode
  logical*1 lstick
  integer*4 neph
  type(brdeph) ephem(1:*)
  integer*4 nwk, jeph, iprn
  real*8 sow, xsat(1:*), dtsat, dtref, v
!
!! local variables
  integer*4 i, j, ieph
  real*8 gme, wearth
  real*8 pi, dt, a, xn, xm, ex, e, v0, vs, vc, phi, ccc, sss, du, dr, di, r, u, xi, xx, yy, &
    xnode, term, xpdot, ypdot, asc, xinc, xp, yp, asctrm
!
!! Angular velocity of the earth  and  Earth's gravitational constant
  data gme, wearth/3.986005d14, 7.2921151467d-5/
  pi = datan(1.d0)*4.d0
!
!! find out the nearest ephem.
  if (.not. lstick .or. jeph .eq. 0) then
    jeph = 0
    dtref = 12.d0    ! hours
    do ieph = 1, neph
      if (ephem(ieph)%svn .eq. iprn) then
        dt = (nwk - ephem(ieph)%week)*168.d0 + (sow - ephem(ieph)%toe)/3600.d0
        if (dabs(dt) .le. dtref) then
          dtref = dabs(dt)
          jeph = ieph
          if (dtref .le. 1.d0) exit
        endif
      endif
    enddo
  endif
  if (jeph .eq. 0) return
!
!! dt
  i = jeph
  dt = (nwk - ephem(jeph)%week)*604800.d0 + sow - ephem(jeph)%toe
!
!! satellite clock correction
  if (cmode(3:3) .eq. 'y') dtsat = ephem(i)%a0 + (ephem(i)%a1 + ephem(i)%a2*dt)*dt

  if (cmode(1:2) .eq. 'nn') return
!
!! compute sat. coordinate in wgs-84
  a = ephem(i)%roota**2
  xn = dsqrt(gme/a/a/a)
  xn = xn + ephem(i)%dn
!
!! iterate to solve EX
  xm = ephem(i)%m0 + xn*dt
  ex = xm
  e = ephem(i)%e
  do j = 1, 12
    ex = xm + e*dsin(ex)
  enddo
!
!! determination of V
  v0 = 1.d0 - e*dcos(ex)
  vs = dsqrt(1.d0 - e*e)*dsin(ex)/v0
  vc = (dcos(ex) - e)/v0
  v = datan2(vs, vc)
!..v=dabs(dasin(vs))
!..if(vc.ge.0.d0) then
!..  if(vs.lt.0.d0) v=2.d0*pi-v
!..else
!..  if(vs.le.0.d0) then
!..    v=pi+v
!..  else
!..    v=pi-v
!..  endif
!..endif
  phi = v + ephem(i)%omega

  ccc = dcos(2*phi)
  sss = dsin(2*phi)
  du = ephem(i)%cuc*ccc + ephem(i)%cus*sss
  dr = ephem(i)%crc*ccc + ephem(i)%crs*sss
  di = ephem(i)%cic*ccc + ephem(i)%cis*sss
  r = a*(1 - e*dcos(ex)) + dr
  u = phi + du

  xi = ephem(i)%i0 + di + ephem(i)%idot*dt
  xx = r*dcos(u)
  yy = r*dsin(u)
  xnode = ephem(i)%omega0 + (ephem(i)%omegadot - wearth)*dt
  xnode = xnode - wearth*ephem(i)%toe
  xsat(1) = xx*dcos(xnode) - yy*dcos(xi)*dsin(xnode)
  xsat(2) = xx*dsin(xnode) + yy*dcos(xi)*dcos(xnode)
  xsat(3) = yy*dsin(xi)

  if (cmode(2:2) .eq. 'n') return
!
!! velocity
  term = (xn*a)/dsqrt(1.d0 - e*e)
  xpdot = -dsin(u)*term
  ypdot = (e + dcos(u))*term
  asc = xnode
  xinc = xi
  xp = xx
  yp = yy
  asctrm = (ephem(i)%omegadot - wearth)
  xsat(4) = xpdot*dcos(asc) - ypdot*dcos(xinc)*dsin(asc) &
            - xp*dsin(asc)*asctrm - yp*dcos(xinc)*dcos(asc)*asctrm
  xsat(5) = xpdot*dsin(asc) + ypdot*dcos(xinc)*dcos(asc) &
            + xp*dcos(asc)*asctrm - yp*dcos(xinc)*dsin(asc)*asctrm
  xsat(6) = ypdot*dsin(xinc)
  return

end
