!
!! map_invnormal.f90
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
!! purpose  : map inversed normal matrix from one-way to multiple-difference ones
!! parameter:
!!    input : MD -- multiple difference struct (SD or DD)
!!    output: QN -- inversed normal matrix (lower triangular part)
!
subroutine map_invnormal(MD, QN, invx)
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include '../header/invnor.h'

  real*8 invx(1:*)
  type(difamb) MD(1:*)
  type(invm) QN
!
!! local
  integer*4 i, j, k, l, ir, ic, ierr
  real*8 op(4)
  real*8, pointer :: hp(:, :)

  data op/1.d0, -1.d0, -1.d0, 1.d0/
  save op

  allocate (hp(QN%indp, QN%indp + QN%nxyz), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(map_invnormal): memory allocation hp '
    call exit(1)
  endif
!
!! M * Q21
  do i = 1, QN%indp
    do j = 1, QN%nxyz
      hp(i, j) = 0.d0
      do k = 1, 4
        if (MD(i)%ipt(k) .eq. 0) exit
        hp(i, j) = hp(i, j) + op(k)*invx(QN%idq(j) + MD(i)%ipt(k))*137.d0/77.d0
      enddo
    enddo
  enddo
  do j = 1, QN%nxyz
    do i = 1, QN%ntot - QN%nxyz
      invx(QN%idq(j) + QN%nxyz + i) = 0.d0
      if (i .le. QN%indp) invx(QN%idq(j) + QN%nxyz + i) = hp(i, j)
    enddo
  enddo
!
!! M * Q22 * Mt
  do i = 1, QN%indp
    do j = 1, i
      hp(i, j) = 0.d0
      do k = 1, 4
        if (MD(i)%ipt(k) .eq. 0) exit
        do l = 1, 4
          if (MD(j)%ipt(l) .eq. 0) exit
          ir = max(MD(i)%ipt(k), MD(j)%ipt(l))
          ic = min(MD(i)%ipt(k), MD(j)%ipt(l))
          hp(i, j) = hp(i, j) + op(k)*invx(QN%idq(ic) + ir)*op(l)*(137.d0/77.d0)**2
        enddo
      enddo
    enddo
  enddo
  do j = QN%nxyz + 1, QN%ntot
    do i = j, QN%ntot
      invx(QN%idq(j) + i) = 0.d0
      if (i .le. QN%nxyz + QN%indp) invx(QN%idq(j) + i) = hp(i - QN%nxyz, j - QN%nxyz)
    enddo
  enddo
  deallocate (hp)

  return
end
