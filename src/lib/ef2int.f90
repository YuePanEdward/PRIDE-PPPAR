!
!! ef2int.f90
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
!! purpose  : rotation from earth-fixed system to inertial system (IERS 1996)
!! parameter:
!!    input : erpfil -- ERP file
!!            jd,sod -- GPS time
!!    output: mate2j -- rotation matrix
!!            rmte2j -- first derivative of rotation matrix
!!            gast   -- GAST angle (radian)
!!            xpole,ypole -- pole movement (arcsec)
!
subroutine ef2int(erpfil, jd, sod, mate2j, rmte2j, gast, xpole, ypole)
  implicit none
  include '../header/const.h'

  integer*4 jd
  real*8 sod, gast, xpole, ypole, mate2j(3, 3), rmte2j(3, 3)
  character*(*) erpfil
!
!! local
  integer*4 jd_tdt, jd_ut1
  real*8 sod_tdt, sod_ut1, sec2rad, p1, p2, p3, epsilon0, dphi, depsi, eqe
  real*8 xhelp(2), taiut1r, dx, dy, dut1_long, dut1, dels
  real*8 pmat(3, 3), nmat(3, 3), gmat(3, 3), rmat(3, 3), xmat(3, 3), ymat(3, 3)
!
!! fuction called
  real*8 sp2000, gst2000

!
!! time second to arc radian
  sec2rad = datan(1.d0)/10800.d0
  call timinc(jd, sod, GPSTDT, jd_tdt, sod_tdt)
!
!! read IGS ERP
  call read_igserp(erpfil, jd, sod, taiut1r, xhelp)
  xpole = xhelp(1)
  ypole = xhelp(2)
!
!! long-term tide effect on ut1 (ut1=ut1r+dx)
  call ut1ut1r(jd_tdt + sod_tdt/86400.d0, dut1_long)
!
!! Tidal variations in Earth rotation
  call ortho_eop(jd_tdt + sod_tdt/86400.d0, dx, dy, dut1)
  xpole = xpole + dx
  ypole = ypole + dy
!
!! Quantity s'
  dels = sp2000(jd_tdt + sod_tdt/86400.d0, MJD2JD)
!
!! polar motion
  call pom2000(xpole/3600.d0/180.d0*PI, ypole/3600.d0/180.d0*PI, dels, xmat)
!
!! nutation
  call nu2000b(jd_tdt + sod_tdt/86400.d0, MJD2JD, dphi, depsi)
!
!! Greenwich
  call timinc(jd, sod, GPSTAI - taiut1r + dut1_long + dut1, jd_ut1, sod_ut1)
  gast = gst2000(jd_ut1 + sod_ut1/86400.d0, MJD2JD, jd_tdt + sod_tdt/86400.d0, MJD2JD, dphi)
!
!! Q matrix
  call cbpn2000(jd_tdt + sod_tdt/86400.d0, MJD2JD, dphi, depsi, pmat)
!
!! transformation
  call t2c2000(xmat, gast, pmat, rmte2j, mate2j)
  return
end
