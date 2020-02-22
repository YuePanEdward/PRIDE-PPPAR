!
!! check_range.f90
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
! purpose   : (PC range - Brd range) to check receiver clock consistence (one epoch)
! instruction : flg --- 0 good
!                       1 problematic
!                       2 no data
!              nobs --- MAXSAT
subroutine check_range(iepo, nobs, mean, pmb, flg)
  implicit none
  include '../header/const.h'

  integer*4 nobs, flg(1:*), iepo
  real*8 mean, pmb(1:*)

  integer*4 i, j, k, imax_dif, mobs
  logical*1 true
  real*8 sig, mean1, sig1, maxdif, dif

  mean = 0.d0
  mobs = 0
! mobs ---------------- the numbers of observation satellite
  do i = 1, nobs
    if (flg(i) .eq. 0) mobs = mobs + 1
  enddo

! when only two obs existed
  if (mobs .le. 2) then
    j = 0
    k = 0
    do i = 1, nobs
      if (flg(i) .eq. 0 .and. j .eq. 0) j = i
      if (flg(i) .eq. 0 .and. j .ne. 0 .and. k .eq. 0) k = i
    enddo
    if (k .eq. 0) flg(j) = 1
    if (k .ne. 0 .and. dabs(pmb(j) - pmb(k)) .gt. 40.d0) then
      flg(j) = 1
      flg(k) = 1
    endif
    return
  endif

! mean value of recv clock
  j = 0
  do i = 1, nobs
    if (flg(i) .eq. 0) then
      mean = mean + pmb(i)/nobs
      j = j + 1
    endif
  enddo
  mean = mean*(dble(nobs)/j)

! sigma of recv clock
  sig = 0.d0
  do i = 1, nobs
    if (flg(i) .eq. 0) sig = sig + (pmb(i) - mean)*(pmb(i) - mean)
  enddo
  sig = dsqrt(sig/(j - 1))

100 continue ! find next largest recv clock

! find largest recv clock
  imax_dif = 0
  maxdif = 0.d0
  do i = 1, nobs
    if (flg(i) .eq. 0) then
      dif = pmb(i) - mean
      if (dabs(dif) .gt. maxdif) then
        maxdif = dabs(dif)
        imax_dif = i
      endif
    endif
  enddo

! new mean recv clock and its sigma when largest one is excluded
  mean1 = 0.d0
  do i = 1, nobs
    if (i .ne. imax_dif .and. flg(i) .eq. 0) mean1 = mean1 + pmb(i)/(j - 1)
  enddo

  sig1 = 0.d0
  do i = 1, nobs
    if (i .ne. imax_dif .and. flg(i) .eq. 0) then
      dif = pmb(i) - mean1
      sig1 = sig1 + dif*dif/(j - 2)
    endif
  enddo
  sig1 = dsqrt(sig1)

! whether remove or not
  true = sig/sig1 .gt. 2.5d0 .or. dabs(mean - mean1) .gt. 500.d0 .or. &
         dabs(pmb(imax_dif) - mean1) .gt. 3*sig1
  if (true .and. dabs(pmb(imax_dif) - mean1) .gt. 100.d0) then
    flg(imax_dif) = 1
    mean = mean1
    sig = sig1
    j = j - 1
  endif
!8/6/2006.Geng J.H.: add "true" in case of dead loop
  if (true .and. sig1 .gt. 1.d3 .and. j .ge. 3) goto 100

!do i=1,nobs
!  if(flg(i).eq.1) then
!    dif=pmb(i)-mean
!    sig=dif-nint(dif/VLIGHT*1.d3)/1.d3*VLIGHT
!  endif
!enddo

  return
end
