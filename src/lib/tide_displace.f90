!
!! tide_displace.f90
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
!! purpose   : compute all tide related station position deformation
!! parameters:
!!
!!        jd, t  -- epoch time
!!        xsite, xsun xlun -- station-, solar- and lunar-positions
!!        rot_f2j, rot_l2f -- rotation matrix from earth-fixed to inertial
!!                and from station to earth-fixed system
!!        sidtm, xpole,ypole -- sideral time, x and y pole positions
!!        olc  -- ocean loading coefficients
!!        dx  -- position correction
!! Author    : Ge Maorong
!
subroutine tide_displace(jd, t, xsit_j, xsit_f, xsun, xlun, rot_f2j, rot_l2f, lat, lon, sidtm, xpole, ypole, olc, disp)
  implicit none
  include '../header/const.h'

  integer*4 jd
  real*8 xsit_j(1:*), xsit_f(1:*), xsun(1:*), xlun(1:*), disp(1:*)  ! in J2000
  real*8 t, lat, lon, xpole, ypole, sidtm, rot_f2j(3, 3), rot_l2f(3, 3), olc(11, 6)
!
!! Fixed constants( Gravitational constants for Sun, Moon, Earth and
!! Love numbers h and l (see IERS Conversion 2003)
  real*8 h20, dh2, l20, dl2, h3, l3
  data h20, dh2, l20, dl2, h3, l3/0.6078d0, -0.0006d0, 0.0847d0, 0.0002d0, 0.292d0, 0.015d0/
!
!! local
  integer*4 i, j, jdutc, iyear, idoy
  real*8 xbod(3), ubod(3), usit(3), xpm, ypm, angle(11), tutc
  real*8 rbod, scal, usit_part, ubod_part, rsit, h2, l2, dummy
  real*8 dxi(3), colat, geoc(3), doy, gb, xs(3), xl(3)
!
!! function called
  real*8 dot, taiutc
!
!! initialization
  do i = 1, 3
    disp(i) = 0.d0
  enddo
  colat = PI/2.d0 - lat
  call timinc(jd, t, 19.d0 - taiutc(jd), jdutc, tutc)
!
!! 1. Displacement due to frequency-independent solid-Earth tide(in J2000)
  call matmpy(xsun, rot_f2j, xs, 1, 3, 3)
  call matmpy(xlun, rot_f2j, xl, 1, 3, 3)
  call solid_earth_tide(xsit_f(1:3)*1.d3, jdutc, tutc/3600.d0, xs*1.d3, xl*1.d3, dxi)
  call matmpy(rot_f2j, dxi, dxi, 3, 3, 1)
  do i = 1, 3
    disp(i) = disp(i) + dxi(i)*1.d-3
  enddo
!
!! 3. Displacement due to the pole tide
! displacement in east, south and radial  in mm, xpole,ypole in seconds of arc (equ.22, pp 67)
! dxi(3) east-north-radial
  doy = jdutc + tutc/86400.d0
  xpm = 0.054d0 + 0.00083d0*(doy - 51545.d0)/365.25d0
  ypm = 0.357d0 + 0.00395d0*(doy - 51545.d0)/365.25d0
  xpm = xpole - xpm
  ypm = -(ypole - ypm)
  dxi(1) = 9.d0*dcos(colat)*(xpm*dsin(lon) - ypm*dcos(lon))
  dxi(2) = 9.d0*dcos(2.d0*colat)*(xpm*dcos(lon) + ypm*dsin(lon))  ! to north
  dxi(3) = -32.d0*dsin(2.d0*colat)*(xpm*dcos(lon) + ypm*dsin(lon))
!
!! rotation matrix from east-north-radial to x-y-z, then to J2000
  call matmpy(rot_l2f, dxi, dxi, 3, 3, 1)
  call matmpy(rot_f2j, dxi, dxi, 3, 3, 1)
  do i = 1, 3
    disp(i) = disp(i) + dxi(i)*1.d-6
  enddo
!
!! 4. Displacement due to ocean-loading
  dxi(1:3) = 0.d0
  call ocean_tidal_loading(jdutc, tutc, olc, 1, 0.d0, dxi)
  call matmpy(rot_l2f, dxi, dxi, 3, 3, 1)
  call matmpy(rot_f2j, dxi, dxi, 3, 3, 1)
  do i = 1, 3
    disp(i) = disp(i) + dxi(i)*1.d-3
  enddo

  return
end
