!
!! read_igserp.f90
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
!! purpose  : read IGS ERP & return TAI-UT1R
!! parameter:
!!    input : erpfil -- ERP file
!!            jd,sod -- requested time
!!    output: tmur   -- TAI-UT1R
!!            pole   -- pole movement
!! modified : May. 13, 2011: leap seconds in terms of igserp table time, not input time
!
subroutine read_igserp(erpfil, jd, sod, tmur, pole)
  implicit none
  include '../header/const.h'

  integer*4 jd
  real*8 sod, tmur, pole(2)
  character*(*) erpfil
!
!! local
  logical*1 lfirst
  integer*4 lfn, i, jd_tdt, ierr, kk
  real*8 lp, dut1, sod_tdt, alpha, fjd, mjdx(2), dat(3, 2)
  character*256 line

!
!! function called
  integer*4 get_valid_unit
  real*8 taiutc

  data lfirst/.true./
  save lfirst, lfn, mjdx, dat

  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file=erpfil, status='OLD', iostat=ierr)
    if (ierr .eq. 0) then
      write (*, '(2a)') '%%%MESSAGE(read_igserp): ERP read ', trim(erpfil)
    else
      write (*, '(2a)') '***ERROR(read_igserp): open file ', trim(erpfil)
      call exit(1)
    endif
    read (lfn, '(a)') line
    if (line(1:9) .ne. 'version 2' .and. line(1:9) .ne. 'VERSION 2') then
      write (*, '(2a)') '***ERROR(read_igserp): unknown version ', line(1:9)
      call exit(1)
    endif
    line = ' '
    do while (index(line, 'MJD') .eq. 0 .or. &
              index(line, 'UT1') .eq. 0 .or. &
              index(line, 'UTC') .eq. 0 .or. &
              index(line, 'LOD') .eq. 0)
      read (lfn, '(a)', end=100) line
    enddo
    read (lfn, '(a)', end=100) line
!
!! read two records
    mjdx = 0
    dat = 0.d0
    read (lfn, '(a)', end=100) line
    read (line, *, err=200) mjdx(1), (dat(i, 1), i=1, 3)
    read (lfn, '(a)', end=100) line
    read (line, *, err=200) mjdx(2), (dat(i, 2), i=1, 3)
    do i = 1, 2
      fjd = mjdx(i) + GPSTDT/86400.d0
      call ut1ut1r(fjd, dut1)
      lp = taiutc(int(mjdx(i)))
      dat(3, i) = lp + dut1 - dat(3, i)*1.d-7
    enddo
  endif
!
!! compare time tag
10 continue
  if (jd + sod/86400.d0 .ge. mjdx(1) .and. jd + sod/86400.d0 .le. mjdx(2)) then
    alpha = (jd + sod/86400.d0 - mjdx(1))/(mjdx(2) - mjdx(1))
    do i = 1, 2
      pole(i) = dat(i, 1) + alpha*(dat(i, 2) - dat(i, 1))
      pole(i) = pole(i)*1.d-6
    enddo
    tmur = dat(3, 1) + alpha*(dat(3, 2) - dat(3, 1))
  else if (jd + sod/86400.d0 .lt. mjdx(1)) then
    backspace lfn
    backspace lfn
    backspace lfn
    read (lfn, '(a)', end=100) line
    read (line, *, err=200) mjdx(1), (dat(i, 1), i=1, 3)
    read (lfn, '(a)', end=100) line
    read (line, *, err=200) mjdx(2), (dat(i, 2), i=1, 3)
    do i = 1, 2
      fjd = mjdx(i) + GPSTDT/86400.d0
      call ut1ut1r(fjd, dut1)
      lp = taiutc(int(mjdx(i)))
      dat(3, i) = lp + dut1 - dat(3, i)*1.d-7
    enddo
    goto 10
  else if (jd + sod/86400.d0 .gt. mjdx(2)) then
    mjdx(1) = mjdx(2)
    do i = 1, 3
      dat(i, 1) = dat(i, 2)
    enddo
    read (lfn, '(a)', end=100) line
    read (line, *, err=200) mjdx(2), (dat(i, 2), i=1, 3)
    fjd = mjdx(2) + GPSTDT/86400.d0
    call ut1ut1r(fjd, dut1)
    lp = taiutc(int(mjdx(2)))
    dat(3, 2) = lp + dut1 - dat(3, 2)*1.d-7
    goto 10
  endif

  return
100 write (*, '(a)') '***ERROR(read_igserp): end of file igserp'
  call exit(1)
200 write (*, '(2a)') '***ERROR(read_igserp): read file igserp ', trim(line)
  call exit(1)
!
!! reset
  Entry igserp_reset()
  close (lfn)
  lfirst = .true.
  return
end
