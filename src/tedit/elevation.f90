!
!! elevation.f90
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
!!!
!! purpose  : compute elevation angle of satellie
!! parameter:
!!   input  : neph,ephem ----- ephemeris number and ephemeris
!!            iprn ----------- satellite number
!!            jd0 ti --------- epoch time (julday)
!!            x,y,z ---------- receiver position
!!            elev ----------- the elevation angle
!!            dist ----------- the distance between satellite and reciver
!!            dtsat ---------- clock correction, units: meter
!!            v -------------- true anomaly
!!            jeph ----------- the prn of ephemris which is chosen
!!            lstick --------- use the last one (jeph) if it is not zero, to avoid discontinuity by two brdeph
!!   output :
!!
subroutine elevation(neph, ephem, iprn, jd0, ti, x, y, z, elev, dist, dtsat, v, jeph, lstick)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'

  type(brdeph) ephem(1:*)
  integer*4 neph, iprn, nepo, jeph, jd0
  real*8 ti, x, y, z, elev, v
  logical*1 lstick

! local
  integer*4 nwk, jd
  real*8 sow, dtsat, xsat(6), tsend, tsend1, dist, r12, r1, cosz, dtref

! function called
  integer*4 modified_julday

! nwk,sow -------- GPS week and second of week
  jd = jd0 + int((ti + 1.d-3)/86400.d0) - modified_julday(5, 1, 1980) - 1
  nwk = jd/7
  jd = jd - nwk*7
  sow = jd*86400.d0 + dmod(ti, 86400.d0)

  tsend = sow - 0.075d0
  do while (.true.)
!!     xsat, dtsat ---- position , velocity and clock correction
    call brdxyz('yyy', lstick, neph, ephem, iprn, jeph, nwk, tsend, xsat, v, dtsat, dtref)
!
!! not found
    if (dabs(dtref) .gt. 5.d0) then   !!!!!!!!!!!!!!!! can be modified
!    write(*,'(a,i3,f5.1)') '***WARNING(elevation): no ephemeris for PRN, &
!          dt(hours)',iprn,dtref
      dist = -1.d0
      return
    endif

    dist = dsqrt((xsat(1) - x)**2 + (xsat(2) - y)**2 + (xsat(3) - z)**2)
    tsend1 = sow - dist/VLIGHT
    if (dabs(tsend - tsend1) .lt. 1.d-9) exit
    tsend = tsend1
  enddo
!
! zenith-distance (approx.)
  r12 = x*(xsat(1) - x) + y*(xsat(2) - y) + z*(xsat(3) - z)
  r1 = dsqrt(x*x + y*y + z*z)
  cosz = r12/r1/dist
  elev = 90.d0 - dacos(cosz)/PI*180.d0
  dtsat = dtsat*VLIGHT

  return
end
