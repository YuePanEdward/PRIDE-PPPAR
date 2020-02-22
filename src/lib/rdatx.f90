!
!! rdatx.f90
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
!! purpose   : read antenna atx file
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!    output : ATX%-- antenna correction
!
subroutine rdatx(fjd_beg, fjd_end, ATX)
  implicit none
  include '../header/const.h'
  include '../header/antatx.h'

  real*8 fjd_beg, fjd_end
  type(antatx) ATX
!
!! local
  logical*1 lfirst, lfound
  integer*4 i, j, k, lfn, iy, imon, id, ih, im, ncol, nrow, icol, irow, ierr
  real*8 sod, fjd0, fjd1
  character*1 atxtyp
  character*256 line
!
!! function called
  integer*4 get_valid_unit, modified_julday

  data lfirst/.true./
  save lfirst, lfn, atxtyp

  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file='abs_igs.atx', status='old', iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(a,i4)') '***ERROR(rdatx): open antenna file abs_igs.atx', ierr
      call exit(1)
    endif
!
!! read header
    line = ' '
    do while (index(line, 'END OF HEADER') .ne. 61)
      read (lfn, '(a)', end=100) line
      if (index(line, 'PCV TYPE / REFANT') .eq. 61) then
        read (line, '(a1)') atxtyp
      endif
    enddo
  endif
  if (atxtyp .ne. 'A') then
    write (*, '(a)') '***ERROR(rdatx): absolute antenna mode in abs_igs.atx'
    call exit(1)
  endif
!
!! find the right antenna
10 continue
  rewind lfn
  lfound = .false.
  do while (.true.)
    line = ' '
    do while (index(line, 'TYPE / SERIAL NO') .ne. 61)
      read (lfn, '(a)', end=100) line
    enddo
    if (ATX%antnam(1:5) .eq. 'BLOCK' .and. ATX%antnam(1:5) .ne. line(1:5) .or. &
        ATX%antnam(1:5) .ne. 'BLOCK' .and. ATX%antnam .ne. line(1:20) .or. &
        ATX%antnum .ne. line(21:40)) cycle
    ATX%dazi = 0.d0
    ATX%dzen = 0.d0
    ATX%nfreq = 0
    ATX%neu = 0.d0
    ATX%pcv(1:30, 0:80, 1:2) = 0.d0
    fjd0 = 0.d0
    fjd1 = 1.d10
    do while (index(line, 'START OF FREQUENCY') .ne. 61)
      read (lfn, '(a)', end=100) line
      if (index(line, 'DAZI') .eq. 61) then
        read (line, *, err=200) ATX%dazi
      else if (index(line, 'ZEN1 / ZEN2 / DZEN') .eq. 61) then
        read (line, *, err=200) ATX%zen1, ATX%zen2, ATX%dzen
      else if (index(line, 'VALID FROM') .eq. 61) then
        read (line, *, err=200) iy, imon, id, ih, im, sod
        fjd0 = modified_julday(id, imon, iy) + ih/24.d0 + im/1440.d0 + sod/86400.d0
      else if (index(line, 'VALID UNTIL') .eq. 61) then
        read (line, *, err=200) iy, imon, id, ih, im, sod
        fjd1 = modified_julday(id, imon, iy) + ih/24.d0 + im/1440.d0 + sod/86400.d0
      else if (index(line, '# OF FREQUENCIES') .eq. 61) then
        read (line, *, err=200) ATX%nfreq
        if (ATX%nfreq .gt. 2) then
          !write (*, '(a)') '###Warning(rdatx): num. of frequency greater than 2'
          ATX%nfreq = 2
        endif
      endif
    enddo
    if ((fjd0 - fjd_beg)*86400.d0 .gt. MAXWND .or. (fjd_end - fjd1)*86400.d0 .gt. MAXWND) cycle
    lfound = .true.
    exit
  enddo
100 continue
  backspace lfn
  if (.not. lfound) then
    if (ATX%antnam(1:5) .ne. 'BLOCK') then
      write (*, '(3a)') '###WARNING(rdatx): Antenna not found ', ATX%antnam, ATX%antnum
      if (ATX%antnam(17:20) .eq. 'NONE') then
        ATX%dazi = 5.d0
        ATX%zen1 = 0.d0
        ATX%zen2 = 90.d0
        ATX%dzen = 5.d0
        ATX%neu = 0.d0
        ATX%pcv(1:30, 0:80, 1:2) = 0.d0
        write (*, '(a)') '$$$MESSAGE(rdatx): Zero antenna applied instead'
      else
        ATX%antnam(17:20) = 'NONE'
        write (*, '(a)') '$$$MESSAGE(rdatx): dome NONE is tried. Search again ...'
        goto 10
      endif
      return
    else
      write (*, '(2a)') '***ERROR(rdatx): atx not found ', ATX%antnam
      call exit(1)
    endif
  endif
!
!! read offset and phase center variation
  ncol = 0
  nrow = 0
  if (ATX%dazi .ne. 0.d0) nrow = int(360.d0/ATX%dazi) + 1
  if (ATX%dzen .ne. 0.d0) ncol = int((ATX%zen2 - ATX%zen1)/ATX%dzen) + 1
  do i = 1, ATX%nfreq
!
!! find frequency
101 line = ' '
    do while (line(61:78) .ne. 'START OF FREQUENCY')
      read (lfn, '(a)') line
      if (index(line, 'END OF ANTENNA') .ne. 0) then
        write (*, '(a)') '***ERROR(rdatx): frquency not found '
        call exit(1)
      endif
    enddo
    if (line(4:4) .ne. 'G') goto 101
    read (line(5:6), *) k
    if (k .gt. 2) goto 101
    ATX%frq(k) = line(4:6)
!
!! reading
    irow = 0
    do while (.true.)
      read (lfn, '(a)') line
      if (line(61:76) .eq. 'END OF FREQUENCY') exit
      if (line(61:77) .eq. 'NORTH / EAST / UP') then
        read (line, *, err=200) (ATX%neu(j, k), j=1, 3)
!! convert unit from mm to m
        do j = 1, 3
          ATX%neu(j, k) = ATX%neu(j, k)*1.d-3
        enddo
      else if (line(4:8) .eq. 'NOAZI') then
        read (line(9:), *, err=200) (ATX%pcv(icol, 0, k), icol=1, ncol)
      else
        irow = irow + 1
        read (line(9:), *, err=200) (ATX%pcv(icol, irow, k), icol=1, ncol)
      endif
    enddo
    if (irow .ne. nrow) then
      write (*, '(a)') '***ERROR(rdatx): pcv lost '
      call exit(1)
    endif
!! convert unit from mm to m
    do irow = 0, nrow
      do icol = 1, ncol
        ATX%pcv(icol, irow, k) = ATX%pcv(icol, irow, k)*1.d-3
      enddo
    enddo
  enddo

  return
200 continue
  write (*, '(2a)') '***ERROR(rdatx): read file ', trim(line)
  call exit(1)
end
