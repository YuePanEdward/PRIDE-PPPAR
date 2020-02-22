!
!! get_lsq_args.f90
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
!! purpose  : get arguments and read options for lsq
!! parameter:
!!   output : LCF  -- lsq configure options
!!            SITE -- station infomation
!!            SAT  -- satellite information
!! tester: Y. Pan, X. Chen, J. Zhou, S. Mao
!
subroutine get_lsq_args(LCF, SITE, OB, SAT)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include '../header/satellite.h'
  include 'lsqcfg.h'

  type(lsqcfg) LCF
  type(station) SITE
  type(rnxobr) OB
  type(satellite) SAT(MAXSAT)
!
!! local
  logical*1 lexist
  integer*4 nargs, lfn, i, j, k, iunit, iprn, nfc, ierr
  integer*4 i1, i2, i3
  integer*4 iy, imon, id, ih, im
  real*8 is, seslen
  character*30 sesfil
  character*256 msg, key, rnxpath
  character*256 path
  integer*4 flag, rename
  type(orbhdr) OH
!
!! function called
  integer*4 get_valid_unit, modified_julday, pointer_int
  real*8 timdif
  character*256 findkey, lower_string
!
!! read arguments
  nargs = iargc()
  if (nargs .lt. 1) then
    write (*, '(a)') 'Usage: lsq sesfil'
    call exit(4)
  endif
  call getarg(1, sesfil)
  lfn = get_valid_unit(10)
  open (lfn, file=sesfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_lsq_args): open file ', trim(sesfil)
    call exit(1)
  endif
!
!! start & stop time
  msg = 'Session time'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) iy, imon, id, ih, im, is, seslen
  LCF%jd0 = modified_julday(id, imon, iy)
  LCF%sod0 = ih*3600.d0 + im*60.d0 + is
  call timinc(LCF%jd0, LCF%sod0, seslen, LCF%jd1, LCF%sod1)
!
!! define file names
  call file_name(.false., 'orb', ' ', iy, imon, id, ih, LCF%flnorb)
  call file_name(.false., 'pos', ' ', iy, imon, id, ih, LCF%flnpos)
  call file_name(.false., 'sck', ' ', iy, imon, id, ih, LCF%flnsck)
  call file_name(.false., 'rck', ' ', iy, imon, id, ih, LCF%flnrck)
  call file_name(.false., 'ztd', ' ', iy, imon, id, ih, LCF%flnztd)
  call file_name(.false., 'htg', ' ', iy, imon, id, ih, LCF%flnhtg)
  call file_name(.false., 'amb', ' ', iy, imon, id, ih, LCF%flnamb)
  call file_name(.false., 'res', ' ', iy, imon, id, ih, LCF%flnres)
  call file_name(.false., 'con', ' ', iy, imon, id, ih, LCF%flncon)
  call file_name(.false., 'neq', ' ', iy, imon, id, ih, LCF%flnneq)
  call file_name(.false., 'vmf', ' ', iy, imon, id, ih, LCF%flnvmf)
  call file_name(.false., 'fcb', ' ', iy, imon, id, ih, LCF%flnfcb)
!
!! erp filename
  LCF%flnerp = 'igserp'
!
!! sampling rate
  msg = 'Interval'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%dintv
!
!! pre-eliminate bias
  msg = 'Remove bias'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%lrmbias = .true.
  if (key(1:1) .ne. 'Y') LCF%lrmbias = .false.
!
!! ZTD model
  msg = 'ZTD model'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%ztdmod = trim(key)
!
!! Horizontal Troposphere Gradients
  msg = 'HTG model'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%htgmod = trim(key)
!
!! Rinex directory
  msg = 'Rinex directory'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  rnxpath = trim(key)
  i = len_trim(rnxpath)
  if (rnxpath(i:i) .ne. '/') then
    rnxpath(i + 1:i + 1) = '/'
  endif
!
!! GPS satellites
  rewind lfn
  msg = '+GPS satellites'
  key = ' '
  do while (key(1:15) .ne. msg(1:15))
    read (lfn, '(a)', end=100) key
  enddo
  LCF%nprn = 0
  call rdorbh(LCF%flnorb, iunit, OH)
!! check orbit span
  if (timdif(OH%jd0, OH%sod0, LCF%jd0, LCF%sod0) .gt. MAXWND) then
    LCF%jd0 = OH%jd0
    LCF%sod0 = OH%sod0
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data starting time revised due to orbits ', &
      LCF%jd0, LCF%sod0
  endif
  if (timdif(OH%jd1, OH%sod1, LCF%jd1, LCF%sod1) .lt. -MAXWND) then
    LCF%jd1 = OH%jd1
    LCF%sod1 = OH%sod1
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data ending time revised due to orbits ', &
      LCF%jd1, LCF%sod1
  endif
!! read satellites
  do while (key(1:15) .ne. '-GPS satellites')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    read (key(2:3), *, err=200) iprn
    if (pointer_int(OH%nprn, OH%prn, iprn) .eq. 0) cycle
    LCF%nprn = LCF%nprn + 1
    LCF%prn(LCF%nprn) = iprn
    SAT(LCF%nprn)%prn = iprn
    SAT(LCF%nprn)%iptatx = 0
    SAT(LCF%nprn)%sys = 'G'
    SAT(LCF%nprn)%typ = 'BLOCK'      ! GPS satellites only
    SAT(LCF%nprn)%xscf(1:3) = 0.d0
    SAT(LCF%nprn)%yscf(1:3) = 0.d0
    SAT(LCF%nprn)%zscf(1:3) = 0.d0
    SAT(LCF%nprn)%sclock = 0.d0
  enddo
  close (iunit)
!
!! Stations
  rewind lfn
  msg = '+Station used'
  key = ' '
  do while (key(1:13) .ne. msg(1:13))
    read (lfn, '(a)', end=100) key
  enddo
  do while (key(1:1) .ne. '-')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    SITE%name = key(2:5)
    SITE%skd = key(7:8)
    SITE%map = key(10:12)
    read (key(13:), *, err=200) SITE%dclk0, SITE%cutoff, SITE%dztd0, &
      SITE%qztd, SITE%dhtg0, SITE%qhtg, SITE%sigr, &
      SITE%sigp, (SITE%dx0(i), i=1, 3)
    SITE%cutoff = SITE%cutoff*PI/180.d0
    SITE%undu = 0.d0
    SITE%rclock = 0.d0
    do i = 1, LCF%nprn
      SITE%first(i) = .true.
      SITE%prephi(i) = 0.d0
    enddo
!! receiver antenna pointer to atx
    SITE%iptatx = 0
!! file names
    key = lower_string(SITE%name)
    OB%lfnrhd = 0
    call file_name(.false., 'rhd', 'SNAM='//key(1:4), iy, imon, id, ih, OB%rhdfil)
    if (index(SITE%skd, 'K') .ne. 0) then     ! for kinematic & pseudo-kinematic use
      SITE%ikin = 0
      call file_name(.false., 'kin', 'SNAM='//key(1:4), iy, imon, id, ih, SITE%kinfil)
    endif

    SITE%iunit = 0
    SITE%imet = 0
    i = len_trim(rnxpath)
    SITE%obsfil = rnxpath(1:i)
    call file_name(.false., 'rnxo', 'SNAM='//key(1:4), iy, imon, id, ih, SITE%obsfil(i + 1:))
    inquire (file=SITE%obsfil, exist=lexist)
    i = len_trim(SITE%obsfil)
    if (.not. lexist) then
      write (*, '(2a)') '###WARNING(get_lsq_args): observation not exist ', SITE%obsfil(1:i)
      call exit(1)
    endif
!! receiver clock jump file
    SITE%lfnjmp = 0
    path = '.'//SITE%name//'.jmp'
    inquire(file=path, exist=lexist)
    if (lexist) then
      SITE%lfnjmp = get_valid_unit(10)
      open(SITE%lfnjmp, file=path, status='old')
      write(*,'(2a)') '###INFO(get_lsq_args): read clock jump file: ', path
    end if
!! read position
    SITE%rlat = 0.d0
    SITE%rlon = 0.d0
    call read_position(LCF%flnpos, SITE%name, LCF%jd0, LCF%sod0, seslen, SITE%x, ierr)
    if (all(SITE%x(1:3) .eq. 1.d0)) then
      write (*, '(a,a4)') '###WARNING(get_lsq_args): no position ', SITE%name
      call exit(1)
    else
      call xyzblh(SITE%x(1:3)*1.d3, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
      SITE%geod(3) = SITE%geod(3)*1.d-3
      call rot_enu2xyz(SITE%geod(1), SITE%geod(2), SITE%rot_l2f)
      call oceanload_coef(SITE%name, &
                          SITE%geod(1), SITE%geod(2), SITE%rlat, &
                          SITE%rlon, SITE%olc)
    endif
    if(SITE%map(1:3) .eq. 'VM1') call vmf1_grid(LCF%flnvmf, SITE)
!
!! check availability of station
    call read_obsrhd(0, LCF%dintv, LCF%nprn, LCF%prn, OB)
    if (OB%ava_obs .lt. OB%rem_obs) then
      write (*, '(2a)') '###WARNING(get_lsq_args): bad observation quality ', SITE%obsfil(1:i)
      write (*, '(2a,i12)') '###WARNING(get_lsq_args): ', 'obs avaiable: ', OB%ava_obs
      !call exit(1)
    endif
    write (*, '(a4,1x,a4,1x,a,f7.1,a1,i3)') 'STA:', SITE%name, SITE%obsfil(1:i), &
      OB%ava_obs*1.d0/(OB%ava_obs + OB%rem_obs)*100.d0, '%', ierr
  enddo
  close (lfn)

  return
100 continue
  write (*, '(3a)') '***ERROR(get_lsq_args): find option ', trim(msg), trim(key)
  call exit(1)
200 continue
  write (*, '(3a)') '***ERROR(get_lsq_args): read option ', trim(msg), trim(key)
  call exit(1)
end
