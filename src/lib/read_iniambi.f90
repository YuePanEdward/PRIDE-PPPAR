!
!! read_iniambi.f90
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
!! purpose  : read ambiguity estimates
!! parameter:
!!    input : flnamb -- ambiguity file
!!            name,iprn -- station name and satellite prn
!!            fjd0,fjd1 -- start and ending time of ambiguity
!!    output: ambini -- ambiguity estimate
!!            rwini  -- wide-lane ambiguity estimate
!
subroutine read_iniambi(flnamb, name, iprn, fjd0, fjd1, ambini, rwini, iflag)
  implicit none
  include '../header/const.h'

  integer*4 iprn, iflag
  real*8 fjd0, fjd1, ambini, rwini
  character*4 name, name1, name2
  character*20 flnamb
!
!! local
  logical*1 lfirst, lexist
  integer*4 lfn, i, ite, ierr
  real*8 t0, t1
  character*256 line
!
!! function used
  integer*4 get_valid_unit
  character*256 lower_string, upper_string

  data lfirst, lexist/.true., .true./
  save lfirst, lexist, lfn

  ambini = 0.d0
  rwini = 0.d0
  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file=flnamb, status='old', iostat=ierr)
    if (ierr .eq. 0) then
      write (*, '(2a)') '%%%MESSAGE(read_iniambi): ambiguity read ', trim(flnamb)
    else
      lexist = .false.
    endif
  endif
  iflag = 2                ! ambiguity file not exist
  if (.not. lexist) return
!
!! search for correct ambiguity estimate
  iflag = 1                ! ambiguity estimate not exist
  do ite = 1, 2
    do while (.true.)
      read (lfn, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
      read (line, '(4x,i4,2f22.6,2f18.10)') i, ambini, rwini, t0, t1
      name1 = lower_string(name)
      name2 = upper_string(name)
      if ((line(1:4) .ne. name2 .and. line(1:4) .ne. name1) .or. i .ne. iprn) cycle
      if ((t1 - fjd0)*86400.d0 .gt. -MAXWND .and. (fjd1 - t0)*86400.d0 .gt. -MAXWND) then
        iflag = 0          ! ambiguity estimate found
        exit
      else if ((t0 - fjd0)*86400.d0 .gt. MAXWND) then
        exit             ! no need to search following ambiguity estimates
      endif
    enddo
    if (iflag .eq. 0) exit
    if (iflag .eq. 1 .and. ite .eq. 1) rewind lfn
  enddo
!
!! check estimate
  if (iflag .eq. 1) then
    write (*, '(a,a4,i3,2f18.10)') '***ERROR(read_iniambi): ambiguity estimate not found ', name, iprn, fjd0, fjd1
    call exit(1)
  endif

  return
end
