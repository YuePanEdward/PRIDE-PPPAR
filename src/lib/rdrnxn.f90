!
!! rdrnxn.f90
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
!!! purpose   : read rinex navigation file.  choose only records
!!             within thime (fjd0,fjd1) and duplicated records are
!!             ignored.
!!
!! parameters: flneph -- rinex navigation file name
!!             fjd0,fjd1 -- start and stop time within the records are requested
!!             neph -- total number of records
!!             ephem -- broadcast ephemeris
!!
!! motified by: Xingyu Chen : add renix 3
!
subroutine rdrnxn(flneph, fjd0, fjd1, neph, ephem)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'

  type(brdeph) ephem(1:*)
  character*(*) flneph
  integer*4 neph, ierr
  real*8 fjd0, fjd1
!
!! local
  type(brdeph) eph
  integer*4 i1, i2, i3, i4, i5, j, iver, iunit, svn
  logical*1 already
  real*8 sec, dt1, dt2, a0, a1, a2
  character*80 line, fmt1
  character*1 satsys
!
!! function called
  integer*4 get_valid_unit

  neph = 0

  ierr = 0
  iunit = get_valid_unit(10)
  open (iunit, file=flneph, status='OLD', iostat=ierr)
  if (ierr .ne. 0) goto 100

  read (iunit, '(i6,a)', iostat=ierr) iver, line
  if (ierr .ne. 0) return

  line = adjustl(line)
  do while (index(line, 'END OF HEADER') .eq. 0 .and. line(1:1) .ne. ' ')
    read (iunit, '(a)', iostat=ierr) line
    if (ierr .ne. 0) goto 100
    line = adjustl(line)
  enddo

  if (iver .eq. 1) then
    fmt1 = '(i2,5i3,f5.1,3f19.12/5(3x,4d19.12/),3x,4d19.12)'
  elseif (iver .eq. 2) then
    fmt1 = '(i2,5i3,f5.1,3d19.12/5(3x,4d19.12/),3x,4d19.12/)'
  elseif (iver .eq. 3) then ! add rinex version 3
    fmt1 = '(6(4x,4d19.12/))'
  endif

  do while (.true.)
    if (iver .eq. 1 .or. iver .eq. 2) then
      read (iunit, fmt1, err=100, end=200) eph%svn, i1, i2, i3, i4, i5, sec, eph%a0, eph%a1, eph%a2, &
        eph%aode, eph%crs, eph%dn, eph%m0, eph%cuc, eph%e, eph%cus, eph%roota, eph%toe, eph%cic, &
        eph%omega0, eph%cis, eph%i0, eph%crc, eph%omega, eph%omegadot, eph%idot, eph%resvd0, &
        eph%week, eph%resvd1, eph%accu, eph%hlth, eph%tgd, eph%aodc
    elseif (iver .eq. 3) then ! add rinex version 3
      satsys = ""
      read (iunit, '(a1,i2,1X,i4,4(1x,i2),1x,f2.0,3d19.12)', iostat=ierr, end=200) satsys, svn, &
        i1, i2, i3, i4, i5, sec, a0, a1, a2
      if (ierr .ne. 0) cycle
      if (satsys .eq. 'G') then !GPS
        eph%svn = svn
        eph%a0 = a0
        eph%a1 = a1
        eph%a2 = a2
        read (iunit, fmt1, err=100, end=200) &
          eph%aode, eph%crs, eph%dn, eph%m0, &
          eph%cuc, eph%e, eph%cus, eph%roota, &
          eph%toe, eph%cic, eph%omega0, eph%cis, &
          eph%i0, eph%crc, eph%omega, eph%omegadot, &
          eph%idot, eph%resvd0, eph%week, eph%resvd1, &
          eph%accu, eph%hlth, eph%tgd, eph%aodc
      else
        cycle
      endif
    endif
!
!! check time
    dt1 = 0.d0
    dt2 = 0.d0
    if (fjd0 .ne. 0.d0) dt1 = (eph%week*7 + 44244) + eph%toe/86400.d0 - fjd0
    if (fjd1 .ne. 0.d0) dt2 = (eph%week*7 + 44244) + eph%toe/86400.d0 - fjd1
    if (dt1 .lt. -1.d0/24.d0 .or. dt2 .gt. 1.d0/24.d0) cycle
!
!! in case of repetition
    already = .false.
    do j = 1, neph
      if (ephem(j)%svn .eq. eph%svn .and. ephem(j)%week .eq. eph%week .and. ephem(j)%toe .eq. eph%toe) then
        already = .true.
        exit
      endif
    enddo
    if (.not. already .and. neph .lt. maxeph) then
      neph = neph + 1
      ephem(neph) = eph
    else if (neph .gt. maxeph) then
      write (*, *) '***WARNING(rdrnxn): too many ephemeris records(maxeph,neph),', maxeph, neph
      return
    endif
  enddo
  close (iunit)
  return

100 continue
  write (*, *) '***ERROR(rdrnxn): open/read nav. file,', trim(flneph)
  call exit(1)
200 continue
  return
end
