!
!! rot_scfix2j2000.f90
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
!! purpose   : rotation matrix from spacecraft(GPS)-fixed system to J2000
!!             The yaw-error correction should be implemented later.
!!
!! parameters:
!!        xsat -- satellite coordinates J2000
!!        xsun -- sun coordinates in J2000
!!        xscf,yscf,zscf -- unit vector of sc-fixed x-axes in J2000
!!                 They define the transformation from x(scf) to j2000 as following
!!                 x(j2000) =  [ xscf,yscf,zscf ] x(scf)
!!
!
subroutine rot_scfix2j2000(xsat, xsun, xscf, yscf, zscf)
  implicit none
  include '../header/const.h'

  real*8 xsat(1:*), xsun(1:*), xscf(1:*), yscf(1:*), zscf(1:*)
!
!! local variables
  integer*4 i, stat
  real*8 yaw_t, yaw_nominal, yaw_bias, cosa, sina, ax, ay, dump
  real*8 rad2deg, alpha, u_sc2sun(3)
!
!! functions used
  real*8 dot

  stat = 0 ! satellite status: shadow or not
  rad2deg = 45.d0/datan(1.d0)
!
!! sc-fixed z-axis, from sc to earth
  do i = 1, 3
    zscf(i) = -xsat(i)
  enddo
  call unit_vector(3, zscf, zscf, ax)
!
!! unint vector from sc to sun
  do i = 1, 3
    u_sc2sun(i) = xsun(i) - xsat(i)
  enddo
  call unit_vector(3, u_sc2sun, u_sc2sun, ax)
!
!! angle between u_sc2sun and zscf, if they near parallel yscf as cross product
!! of them can not be defined. Values from the last epoch should be used
  alpha = dot(3, zscf, u_sc2sun)
  if (dabs(alpha) .gt. 1.d0) alpha = nint(alpha)
  alpha = dacos(alpha)*rad2deg
  if (alpha .lt. 13.5d0) stat = 1

  if (alpha .ge. 1.d-6 .and. alpha .le. 180.d0 - 1.d-6) then
!
!! the sc-fixed y-axis.
    call cross(zscf, u_sc2sun, yscf)
    call unit_vector(3, yscf, yscf, ax)
!
!! the sc-fixed z-axis
    call cross(yscf, zscf, xscf)
  else
    write (*, '(a,f14.4)') '###WARNING(rot_scfix2j200): scx/y no definition, previous epoch used ', alpha
  endif

  return
!
!! yaw_nominal is defined in svnav.dat (GAMIT), yaw_t to be calculated and estimated
!! now set to zero
  yaw_bias = yaw_t - yaw_nominal
  yaw_bias = dmod(yaw_bias, 360.d0)
!
!! yaw_bias correction
  if (dabs(yaw_bias) .gt. 0.1d0) then
    yaw_bias = yaw_bias/rad2deg
    cosa = dcos(yaw_bias)
    sina = dsin(yaw_bias)
    do i = 1, 3
      ax = xscf(i)*cosa + yscf(i)*sina
      ay = -xscf(i)*sina + yscf(i)*cosa
      xscf(i) = ax
      yscf(i) = ay
    enddo
  endif

  return
end
