!
!! ambslv.f90
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
!! purpose  : resolve ambiguity
!! parameter:
!!    input : ncad   -- # of candidate ambiguities
!!            q22    -- cofactor matrix
!!    output: bias   -- float / fixed ambiguity estimates
!!            disall -- norm of optimum & suboptimum solutions
!
subroutine ambslv(ncad, q22, bias, disall)
  implicit none

  integer*4 ncad
  real*8 q22(1:*), bias(1:*), disall(1:*)
!
!! local
  real*8 dump

  if (ncad .gt. 1) then
    call lambda4(ncad, q22, bias, disall)
  else
    dump = bias(1)
    bias(1) = nint(bias(1))*1.d0
    dump = bias(1) - dump
    disall(1) = dump/q22(1)*dump
    dump = 1.d0 - dabs(dump)
    disall(2) = dump/q22(1)*dump
  endif

  return
end
