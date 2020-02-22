!
!! read_recclk.f90
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
!! purpose  : read receiver clock correction
!! parameter:
!!    input : recfil -- receiver clock file
!!            name   -- requested station
!!            jd,sod -- requested time
!!    output: x0     -- receiver clock correction (in meter)
!! note     : only read corresponding epoch without interpolation
!
subroutine read_recclk(recfil, name, jd, sod, x0)
  implicit none
  include '../header/const.h'

  integer*4 jd
  real*8 sod, x0
  character*(*) name, recfil
!
!! local
  logical*1 lfirst, lexist
  integer*4 i, lfn, nsit, jdf, jdx, iy, imon, id, ih, im, ierr
  real*8 dintv, a0(MAXSIT), sodf, sodx, sec, dt, dummy
  character*4 lname(MAXSIT)
  character*256 line
!
!! function called
  integer*4 get_valid_unit, modified_julday, pointer_string
  real*8 timdif

  data lfirst, lexist/.true., .false./
  save lfirst, lexist, lfn, jdf, sodf, nsit, lname, a0

  x0 = 0.d0
  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file=recfil, status='old', iostat=ierr)
    if (ierr .eq. 0) then
      write (*, '(2a)') '%%%MESSAGE(read_recclk): receiver clocks read ', trim(recfil)
    else
      return
    endif
    lexist = .true.
!
!! set interval
    dintv = 30.d0  ! default value in most cases
    line = ' '
    do while (index(line, 'END OF HEADER') .eq. 0)
      read (lfn, '(a)') line
      if (index(line, 'INTERVAL') .ne. 0) read (line, *) dintv
    enddo
!
!! read one epoch
    nsit = 0
    jdf = 0
    do while (.true.)
      read (lfn, '(a)') line
      nsit = nsit + 1
      if (nsit .gt. MAXSIT) then
        write (*, '(a)') '***ERROR(read_recclk): max site exceeded '
        call exit(1)
      endif
      read (line, '(a4,i5,4i3,f10.6,f17.6,f14.6)', err=100) lname(nsit), iy, imon, id, ih, im, sec, a0(nsit), dummy
      a0(nsit) = a0(nsit) + dummy
      jdx = modified_julday(id, imon, iy)
      sodx = ih*3600.d0 + im*60.d0 + sec
      if (jdf .eq. 0) then
        jdf = jdx
        sodf = sodx
      else if (timdif(jdx, sodx, jdf, sodf) .gt. MAXWND) then
        nsit = nsit - 1
        backspace lfn
        exit
      endif
    enddo
  endif
  if (.not. lexist) return
!
!! check time tag
10 dt = timdif(jd, sod, jdf, sodf)
  if (dabs(dt) .lt. MAXWND) then
    i = pointer_string(nsit, lname, name)
    if (i .eq. 0) then
      write (*, '(2a)') '###WARNING(read_recclk): station not found ', name(1:4)
      return
    endif
    x0 = a0(i)/VLIGHT
  else if (dt .lt. -MAXWND) then
    write (*, '(a,i7,f10.2)') '***ERROR(read_recclk): before ref. time ', jd, sod
    call exit(1)
  else if (dt .gt. MAXWND) then
    nsit = 0
    jdf = 0
    line = ' '
    do while (.true.)
      read (lfn, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
      nsit = nsit + 1
      if (nsit .gt. MAXSIT) then
        write (*, '(a)') '***ERROR(read_recclk): max site exceeded '
        call exit(1)
      endif
      read (line, '(a4,i5,4i3,f10.6,f17.6,f14.6)', err=100) lname(nsit), iy, imon, id, ih, im, sec, a0(nsit), dummy
      a0(nsit) = a0(nsit) + dummy
      jdx = modified_julday(id, imon, iy)
      sodx = ih*3600.d0 + im*60.d0 + sec
      if (jdf .eq. 0) then
        jdf = jdx
        sodf = sodx
      else if (timdif(jdx, sodx, jdf, sodf) .gt. MAXWND) then
        nsit = nsit - 1
        backspace lfn
        exit
      endif
    enddo
    goto 10
  endif

  return
100 write (*, '(2a)') '***ERROR(read_recclk): read file ', trim(line)
  call exit(1)
end
