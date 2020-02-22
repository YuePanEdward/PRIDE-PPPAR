!
!! lsq_add_ambcon.f90
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
!! purpose  : add ambiguity resolution constraint
!! parameter:
!!    input : jd,sod -- epoch time
!!            LCF    -- lsq control struct
!!            SITE   -- station struct
!!    output: NM,PM  -- normal matrix & parameter array
!
subroutine lsq_add_ambcon(jd, sod, LCF, SITE, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  integer*4 jd
  real*8 sod
  type(lsqcfg) LCF
  type(station) SITE
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  logical*1 lfirst, lexist
  integer*4 i, j, ir, ic, lfn, jdw, jdx, iy(2), imon(2), id(2), ih(2), im(2), isat, jsat, ipar, ix(2), ierr
  real*8 is(2), bc, iwl, rnl, fip, sodw, sodx, c(2), ini(4), op(4)
  character*2 ctyp
  character*4 cname, dname
  character*256 line
!
!! function called
  integer*4 get_valid_unit, modified_julday
  real*8 timdif

  data lfirst, lexist/.true., .false./, op/1.d0, -1.d0, -1.d0, 1.d0/
  save lfirst, lexist, lfn, c, ctyp

  if (lfirst) then
    lfirst = .false.
    inquire (file=LCF%flncon, exist=lexist)
    if (lexist) then
      write (*, '(a)') '%%%MESSAGE(lsq_add_ambcon): ambiguity constraint imposed'
      lfn = get_valid_unit(10)
      open (lfn, file=LCF%flncon, status='old', iostat=ierr)
      if (ierr .ne. 0) then
        write (*, '(2a)') '***ERROR(lsq_add_ambcon): open file ', trim(LCF%flncon)
        call exit(1)
      endif
      bc = 77.d0/60.d0
      c(1) = bc/(bc**2 - 1.d0)
      c(2) = 77.d0/137.d0
!
!! read header
      read (lfn, '(a)') line
      do while (index(line, 'END OF HEADER') .eq. 0)
        if (index(line, 'TYPE OF CONSTRAINT') .ne. 0) then
          ctyp = line(5:6)
        endif
        read (lfn, '(a)') line
      enddo
    endif
  endif
  if (.not. lexist) return

10 read (lfn, '(a)', end=100) line
  if (ctyp .eq. 'SD') then
    dname = ' '
    read (line, '(a4,1x,2(1x,i2.2,1x),2(i4,4i3,f10.6,1x),2f13.0,f8.3)', err=200) cname, isat, jsat, &
      iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2), iwl, rnl, fip
  endif
!
!! if station missing ...
  j = 0
  if (SITE%name .eq. cname) j = j + 1
  if (ctyp .eq. 'SD' .and. j .lt. 1) goto 10
!
!! time comaprison
  jdw = modified_julday(id(1), imon(1), iy(1))
  sodw = ih(1)*3600.d0 + im(1)*60.d0 + is(1)
  jdx = modified_julday(id(2), imon(2), iy(2))
  sodx = ih(2)*3600.d0 + im(2)*60.d0 + is(2)
!
!! if time is outside required duration
  if (timdif(LCF%jd0, LCF%sod0, jdw, sodw) .gt. 0.d0) then
    jdw = LCF%jd0
    sodw = LCF%sod0
  endif
  if (timdif(jdx, sodx, LCF%jd1, LCF%sod1) .gt. 0.d0) then
    jdx = LCF%jd1
    sodx = LCF%sod1
  endif
  if (timdif(jdx, sodx, jdw, sodw) .le. 0.d0) goto 10
!
!! check whether add new constraints
  if (timdif(jdx, sodx, jd, sod) .gt. MAXWND) then
    backspace lfn
    return
  endif
!
!! search corresponding one-way ambiguities in normal euation
  ix = 0
  do i = NM%nc + NM%np + 1, NM%imtx
    ipar = NM%iptp(i)
    if (PM(ipar)%pname(1:4) .ne. 'AMBC') cycle
    if (SITE%name .ne. cname) cycle
    if (LCF%prn(PM(ipar)%pcode(2)) .ne. isat .and. LCF%prn(PM(ipar)%pcode(2)) .ne. jsat) cycle
    if ((PM(ipar)%ptbeg - jdw)*86400.d0 - sodw .gt. MAXWND .or. (PM(ipar)%ptend - jdx)*86400.d0 - sodx .lt. -MAXWND) cycle
    if (SITE%name .eq. cname .and. LCF%prn(PM(ipar)%pcode(2)) .eq. isat) then
      ix(1) = PM(ipar)%ipt
      ini(1) = PM(ipar)%xini
    else if (SITE%name .eq. cname .and. LCF%prn(PM(ipar)%pcode(2)) .eq. jsat) then
      ix(2) = PM(ipar)%ipt
      ini(2) = PM(ipar)%xini
    endif
    if (ctyp .eq. 'SD' .and. ix(1) .ne. 0 .and. ix(2) .ne. 0) exit
  enddo
  if (ix(1) .eq. 0 .or. ix(2) .eq. 0) then
    write (*, '(a,a4,1x,a4,2i3,2(1x,i4,4i3,f10.6))') '***ERROR(lsq_add_ambcon): ambiguity not found ', cname, dname, &
      isat, jsat, iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2)
    call exit()
  endif
!
!! add constraint to normal equation
  if (ctyp .eq. 'SD') then
    bc = c(1)*iwl + c(2)*(rnl + fip) - ini(1) + ini(2)
  endif
  do i = 1, 2
    if (ix(i) .eq. 0) exit
    do j = i, 2
      if (ix(j) .eq. 0) exit
      ir = min(ix(i), ix(j))
      ic = max(ix(i), ix(j))
      NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*op(j)*1.d12
    enddo
    ir = ix(i)
    ic = NM%imtx + 1
    NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*1.d12*bc
  enddo
  NM%ltpl = NM%ltpl + 1.d12*bc**2
  NM%nobs = NM%nobs + 1

  goto 10

100 close (lfn)
  return
200 write (*, '(a,/,a)') '***ERROR(lsq_add_ambcon): read file ', trim(line)
  call exit()
end
