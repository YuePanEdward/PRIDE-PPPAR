!
!! get_ant_corr.f90
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
!! purpose   : get antenna pointer in the whole table for receivers and satellites
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!             antnam,antnum   -- antenna name and serial number
!!    output : iptatx  -- pointer to global table
!!             enu     -- antenna offset
!
subroutine get_ant_ipt(fjd_beg, fjd_end, antnam, antnum, iptatx, enu)
  implicit none
  include '../header/const.h'
  include '../header/antatx.h'

  integer*4 iptatx, ript, sipt
  real*8 fjd_beg, fjd_end, zeni, azim, nadir, enu(1:*), var(2)
  character*20 antnam, antnum
!
!! local
  integer*4 i, j, natx, izen, iazi
  real*8 rad2deg, alpha, zen, azi, nad, x1, x2
  type(antatx) AX(MAXSIT + MAXSAT), ATX

  data natx/0/
  save natx, AX

  do i = 1, natx
    if (antnam .ne. AX(i)%antnam) cycle
    j = len_trim(antnum)
    if (j .ne. 0 .and. antnum(1:j) .ne. AX(i)%antnum(1:j)) cycle
    iptatx = i
    exit
  enddo
  if (iptatx .ne. 0) goto 5
!
!! if not in memory, read from file
  ATX%antnam = antnam
  ATX%antnum = antnum
  call rdatx(fjd_beg, fjd_end, ATX)
  natx = natx + 1
  AX(natx) = ATX
  antnam = ATX%antnam
  iptatx = natx
!
!! get antenna phase offset
5 continue
  if (antnam(1:5) .ne. 'BLOCK') then
    enu(1) = AX(iptatx)%neu(2, 1)
    enu(2) = AX(iptatx)%neu(1, 1)
    enu(3) = AX(iptatx)%neu(3, 1)
    enu(4) = AX(iptatx)%neu(2, 2)
    enu(5) = AX(iptatx)%neu(1, 2)
    enu(6) = AX(iptatx)%neu(3, 2)
  else
    do i = 1, 3
      enu(i) = AX(iptatx)%neu(i, 1)
      enu(i + 3) = AX(iptatx)%neu(i, 2)
    enddo
  endif

  return
!
!! purpose   : get antenna pcv for receivers and satellites
!! parameter :
!!    input  : ript,sipt -- pointer of receiver and satellite
!!             zeni  -- zenith angle in radian
!!             azim  -- azimuth angle in radian
!!             nadir -- nadir angle in radian
!!    output : var   -- phase center variation
!
  Entry get_ant_pcv(ript, sipt, zeni, azim, nadir, var)

  rad2deg = 180.d0/PI
  var(1) = 0.d0
  var(2) = 0.d0
!
!! receiver pcv
  zen = zeni*rad2deg
  azi = azim*rad2deg
  if (zen .gt. AX(ript)%zen2) zen = AX(ript)%zen2
  if (zen .lt. AX(ript)%zen1) zen = AX(ript)%zen1
  if (azi .lt. 0.d0) azi = azi + 360.d0
!
!! azimuth dependent
  iazi = 0
  if (AX(ript)%dazi .ne. 0.d0) iazi = int(azi/AX(ript)%dazi) + 1
!
!! zenith dependent
  izen = int((zen - AX(ript)%zen1)/AX(ript)%dzen) + 1
  do i = 1, 2
    x1 = AX(ript)%pcv(izen, iazi, i)
    x2 = AX(ript)%pcv(izen + 1, iazi, i)
    if (iazi .ne. 0) then
      alpha = azi/AX(ript)%dazi - iazi + 1
      x1 = x1 + (AX(ript)%pcv(izen, iazi + 1, i) - AX(ript)%pcv(izen, iazi, i))*alpha
      x2 = x2 + (AX(ript)%pcv(izen + 1, iazi + 1, i) - AX(ript)%pcv(izen + 1, iazi, i))*alpha
    endif
    alpha = (zen - AX(ript)%zen1)/AX(ript)%dzen - izen + 1
    var(i) = var(i) + x1 + (x2 - x1)*alpha
  enddo
!
!! satellite pcv
  if (AX(sipt)%dzen .eq. 0.d0) return
  nad = nadir*rad2deg
  if (nad .gt. AX(sipt)%zen2) nad = AX(sipt)%zen2
  if (nad .lt. AX(sipt)%zen1) nad = AX(sipt)%zen1
!
!! only zenith dependent
  izen = int((nad - AX(sipt)%zen1)/AX(sipt)%dzen) + 1
  do i = 1, 2
    x1 = AX(sipt)%pcv(izen, 0, i)
    x2 = AX(sipt)%pcv(izen + 1, 0, i)
    alpha = (nad - AX(sipt)%zen1)/AX(sipt)%dzen - izen + 1
    var(i) = var(i) + x1 + (x2 - x1)*alpha
  enddo

  return
end
