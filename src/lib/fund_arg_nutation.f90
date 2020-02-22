!
!! fund_arg_nutation.f90
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
!! purpose   : fundamental arguments of nutation
!!
!! parameters: rmjd_tdt -- mjd in TDT system
!!             arg      -- five nutation arguments
!!
!
subroutine fund_arg_nutation(rmjd_tdt, arg)
  implicit none

  real*8 rmjd_tdt, arg(1:*)
!
!! local
  logical*1 first
  integer*4 i, j
  real*8 coef(5, 5), deg2rad, sec2rad, ttc
!
!! IERS resolution 1996, pp23
  data coef/134.96340251d0, 1717915923.2178d0, 31.8792d0, 0.051635d0, -0.00024470d0, &
    357.52910918d0, 129596581.0481d0, -0.5532d0, 0.000136d0, -0.00001149d0, &
    93.27209062d0, 1739527262.8478d0, -12.7512d0, -0.001037d0, 0.00000417d0, &
    297.85019547d0, 1602961601.2090d0, -6.3706d0, 0.006593d0, -0.00003169d0, &
    125.04455501d0, -6962890.2665d0, 7.4722d0, 0.007702d0, -0.00005939d0/
  data first/.true./
  save coef, first

  if (first) then
    first = .false.
    coef(2, 5) = -6962890.5431d0
    deg2rad = datan(1.d0)/45.d0
    sec2rad = datan(1.d0)/162000.d0
    do i = 1, 5
      coef(1, i) = coef(1, i)*deg2rad
      do j = 2, 5
        coef(j, i) = coef(j, i)*sec2rad
      enddo
    enddo
  endif
!
!! calculate the angles
  ttc = (rmjd_tdt - 51544.5d0)/36525.d0
  do i = 1, 5
    arg(i) = coef(1, i) + (coef(2, i) + (coef(3, i) + (coef(4, i) + coef(5, i)*ttc)*ttc)*ttc)*ttc
  enddo

  return
end
