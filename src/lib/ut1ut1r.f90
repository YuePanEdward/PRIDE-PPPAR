!
!! ut1ut1r.f90
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
!! purpose  : compute tide variation on UT1 due to Earth rotation  (IERS1996, 2003)
!! parameter: rmjd -- mod. julian date
!!            dut1 -- ut1 correction
!!                 UT1 = UT1R + dut1
!! notice   : if ut1. table is UT1, dut1 must be subtracted for interpolation.
!!            dut1 must be added to the interpolated UT1R for coordinate trans.
!
subroutine ut1ut1r(rmjd, dut1)
  implicit none
  include 'ut1tid.h'

  real*8 rmjd, dut1
!
!! local
  integer*4 i, j
  real*8 argu(5), arg
!
!! get arguement
  call fund_arg_nutation(rmjd, argu)
  dut1 = 0.d0
  do i = 1, 62
    arg = 0.d0
    do j = 1, 5
      arg = arg + mi(j, i)*argu(j)
    enddo
    dut1 = dut1 + sc(1, i)*dsin(arg) + sc(2, i)*dcos(arg)
  enddo
!
!! 0.1 mimiseconds to rad.
  dut1 = dut1*1d-4
  return
end
