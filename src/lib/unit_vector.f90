!
!! unit_vector.f90
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
!! purpose  : normalize a vector
!!
!! parameter: n -- dimmension of  the vector
!!            v -- input vector
!!            u -- ouput unit vector
!!            length -- length of the vector
!!
!!
!! last mof.: 31-May-2003 by Maorong GE, CLEAN

subroutine unit_vector(n, v, u, length)
  implicit none
  integer*4 n
  real*8 v(1:*), u(1:*), length
!
!! local
  integer*4 i

  length = 0.d0
  do i = 1, n
    length = length + v(i)*v(i)
  enddo
  length = dsqrt(length)

  do i = 1, n
    u(i) = v(i)/length
  enddo

  return
end
