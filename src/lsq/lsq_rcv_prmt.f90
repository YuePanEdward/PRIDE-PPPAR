!
!! lsq_rcv_prmt.f90
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
!! purpose   : recover pre-eliminated parameters and residual from tmp files
!! parameter :
!!    input  : lfncid,lfnobs,lfnrem -- tmp files for recovery
!!             lfnres               -- residual file
!!             LCF                  -- least squares control
!!             SITE                 -- station control
!!             NM,PM                -- normal equation & parameters
!
subroutine lsq_rcv_prmt(lfncid, lfnobs, lfnrem, lfnres, LCF, SITE, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  integer*4 lfncid, lfnobs, lfnrem, lfnres
  type(lsqcfg) LCF
  type(station) SITE
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  character*2 cid
  integer*2 iflg
  integer*4 i, j, k, ipar, isit, isat, iobs, ntot, ipt(0:MAXPAR), jd, ierr
  integer*4 lfnpos, tmppos, lfnrck, tmprck, lfncck, tmpcck, lfnztd, lfnhtg, tmpztd, tmphtg, tmpkin, jdr, iy, imon, id, ih, im
  real*8 cof(MAXPAR), x(3), xc(6), sodr, phase, range, wphs, wrng, mw, elev, azim, dmap, wmap
  real*8 zdd, zwd, gni, gei, gnx, gex, xini, xcor, dummy
  character*4 nme
!zwx
  integer*4 openid
  logical*1 lopen
!
!! function called
  integer*4 get_valid_unit
  real*8 timdif
!
!! open necessary files
  tmprck = get_valid_unit(10)
  open (tmprck, file='tmp'//LCF%flnrck, form='unformatted')
  write (tmprck) '0000', 0, 0, 0, 0, 0, 0.d0, 0.d0, 0.d0
  if (LCF%ztdmod(1:3) .ne. 'FIX') then
    tmpztd = get_valid_unit(10)
    open (tmpztd, file='tmp'//LCF%flnztd, form='unformatted')
    write (tmpztd) '0000', 0, 0, 0, 0, 0, 0.d0, 0.d0, 0.d0, 0.d0
  endif
  if (LCF%htgmod(1:3) .ne. 'NON' .and. LCF%htgmod(1:3) .ne. 'FIX') then
    tmphtg = get_valid_unit(10)
    open (tmphtg, file='tmp'//LCF%flnhtg, form='unformatted')
    write (tmphtg) '0000', 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
  endif
  tmpkin = get_valid_unit(10)
  open (tmpkin, file='tmpkin', form='unformatted')
  write (tmpkin) 0, 0, 0.d0, 0.d0, 0.d0, 0.d0, 0
!
!! recover parameters and residuals
  jdr = 0; sodr = 0.d0
  do while (.true.)
    backspace lfncid
    read (lfncid, err=100) cid
    if (cid .eq. '00') exit
    if (cid .eq. 'ob') then
      backspace lfnobs
      read (lfnobs) jd, mw, isit, isat, ntot, ipt(0), (ipt(k), cof(k), k=1, ntot), &
                    phase, range, wphs, wrng, iflg, elev, azim, dmap, wmap
      dummy = 0.d0
      do i = 1, ntot
        dummy = dummy + PM(ipt(i))%xcor*cof(i)
      enddo
      if (ipt(0) .eq. 0) then
        range = range - dummy
      else
        phase = phase - dummy
        range = range - dummy + PM(ipt(ipt(0)))%xcor*cof(ipt(0))
        PM(ipt(ipt(0)))%xrms = PM(ipt(ipt(0)))%xrms + phase**2
      endif
      if (dabs((jdr - jd)*86400.d0 + sodr - mw) .gt. MAXWND) then
        jdr = jd
        sodr = mw
        call mjd2date(jd, mw, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jd, mw
      endif
      elev = elev/PI*180.d0
      azim = azim/PI*180.d0
      write (lfnres, '(2i3,2f10.3,2d16.8,i3,f8.3,f9.3)') isit, isat, phase, range, wphs, wrng, iflg, elev, azim
      backspace lfnobs
    else if (cid .eq. 'de') then
      backspace lfnrem
      read (lfnrem) isit, jd, mw, ntot, (ipt(k), k=1, ntot)
      if (dabs((jd - jdr)*86400.d0 + mw - sodr) .gt. MAXWND) then
        jdr = jd
        sodr = mw
        call mjd2date(jd, mw, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jd, mw
      endif
      do i = 1, ntot
        write (lfnres, '(2i3,2f10.3,2d16.8,i3,f8.3,f9.3)') isit, ipt(i), 0.d0, 0.d0, 0.d0, 0.d0, 4, 0.d0, 0.d0
      enddo
      backspace lfnrem
    else if (cid .eq. 'pc') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot), &
                    PM(ipar)%iobs, PM(ipar)%xini, PM(ipar)%ptime(1:2)
      dummy = cof(ipt(0))
      cof(ipt(0)) = -PM(ipar)%map*PM(ipar)%rw**2
      if (PM(ipar)%ipt .lt. 0) then
        PM(ipar)%ipt = 0
        cof(ipt(0)) = 0.d0
        NM%np = NM%np + 1
      endif
      do i = 1, ntot - 1
        cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
      enddo
      PM(ipar)%xcor = cof(ntot)/dummy
      if (PM(ipar)%iobs .gt. 0) then
        mw = PM(ipar)%ptime(1)
        if (PM(ipar)%pname(1:6) .eq. 'RECCLK') then
          call mjd2date(0, mw, iy, imon, id, ih, im, dummy)
          if (dabs(dummy-60.d0) .lt. maxwnd) then
            im = im + 1
            dummy = 0.d0
          endif
          write (tmprck) SITE%name, iy, imon, id, ih, im, dummy, PM(ipar)%xini, PM(ipar)%xcor
        else if (PM(ipar)%pname(1:6) .eq. 'SATCLK') then
          call mjd2date(0, mw, iy, imon, id, ih, im, dummy)
          if (dabs(dummy-60.d0) .lt. maxwnd) then
            im = im + 1
            dummy = 0.d0
          endif
          write (tmpcck) LCF%prn(PM(ipar)%pcode(2)), iy, imon, id, ih, im, dummy, (PM(ipar)%xini + PM(ipar)%xcor)/VLIGHT
        else if (PM(ipar)%pname(1:4) .eq. 'HTGC') then
          write (tmphtg) SITE%name, PM(ipar)%ptime(1:2), (PM(ipar + i)%xini, PM(ipar + i)%xcor, i=0, 1)
        else if (PM(ipar)%pname(1:5) .eq. 'STAPX') then
          write (tmpkin) PM(ipar)%pcode(1), int(mw), (mw - int(mw))*864.d2, &
                         (PM(ipar + i)%xini + PM(ipar + i)%xcor, i=0, 2), PM(ipar)%iobs
        endif
      endif
      backspace lfnrem
    else if (cid .eq. 'pz') then
      backspace lfnrem
      read (lfnrem) isit, ipar, zdd, zwd, dummy
      if ((dummy - PM(ipar)%ptime(1))*86400.d0 .gt. -MAXWND .and. &
          (PM(ipar)%ptime(2) - dummy)*86400.d0 .gt. -MAXWND .and. PM(ipar)%iobs .gt. 0) then
        call mjd2date(0, dummy, iy, imon, id, ih, im, mw)
        if (dabs(mw-60.d0) .lt. maxwnd) then
          im = im + 1
          mw = 0.d0
        endif
        write (tmpztd) SITE%name, iy, imon, id, ih, im, mw, zdd, zwd, PM(ipar)%xcor
      endif
      backspace lfnrem
    else if (cid .eq. 'am') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot)
      NM%ns = NM%ns + 1
      do i = 1, ntot - 1
        if (ipt(i) .ne. ipar) then
          cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
        endif
      enddo
      PM(ipar)%xcor = cof(ntot)/cof(ipt(0))
      backspace lfnrem
    endif
    backspace lfncid
  enddo
!
!! -------------WRITE FIELS--------------
!
!! write receiver clock file
  inquire (file=LCF%flnrck, opened=lopen, number=openid)   !add by zwx 20141031
  if (lopen) close (openid)
  lfnrck = get_valid_unit(10)
  open (lfnrck, file=LCF%flnrck)
  write (lfnrck, '(a14,46x,a)') 'Receiver Clock', 'COMMENT'
  write (lfnrck, '(f9.2,51x,a)') LCF%dintv, 'INTERVAL'
  write (lfnrck, '(60x,a)') 'END OF HEADER'
  do while (.true.)
    backspace tmprck
    read (tmprck) nme, iy, imon, id, ih, im, mw, xini, xcor
    if (nme .eq. '0000') exit
    write (lfnrck, '(a4,i5,4i3,f10.6,f17.6,f14.6)') nme, iy, imon, id, ih, im, mw, xini, xcor
    backspace tmprck
  enddo
  close (lfnrck); close (tmprck, status='delete')
!
!! write atm zenith delay file
  if (LCF%ztdmod(1:3) .ne. 'FIX') then
    inquire (file=LCF%flnztd, opened=lopen, number=openid) !add by zwx 20141031
    if (lopen) close (openid)
    lfnztd = get_valid_unit(10)
    open (lfnztd, file=LCF%flnztd)
    write (lfnztd, '(a25,35x,a)') 'Zenith Tropospheric Delay', 'COMMENT'
    write (lfnztd, '(f9.2,51x,a)') LCF%dintv, 'INTERVAL'
    write (lfnztd, '(60x,a)') 'END OF HEADER'
    do while (.true.)
      backspace tmpztd
      read (tmpztd) nme, iy, imon, id, ih, im, mw, zdd, zwd, xcor
      if (nme .eq. '0000') exit
      write (lfnztd, '(a4,i5,4i3,f10.6,3f11.6)') nme, iy, imon, id, ih, im, mw, zdd, zwd, xcor
      backspace tmpztd
    enddo
    close (lfnztd); close (tmpztd, status='delete')
  endif
!
!! write horizontal troposphere gradient
  if (LCF%htgmod(1:3) .ne. 'NON' .and. LCF%htgmod(1:3) .ne. 'FIX') then
    inquire (file=LCF%flnhtg, opened=lopen, number=openid) !add by zwx 20141031
    if (lopen) close (openid)
    lfnhtg = get_valid_unit(10)
    open (lfnhtg, file=LCF%flnhtg)
    write (lfnhtg, '(a32,28x,a)') 'Horizontal Troposphere Gradients', 'COMMENT'
    write (lfnhtg, '(f9.2,51x,a)') LCF%dintv, 'INTERVAL'
    write (lfnhtg, '(60x,a)') 'END OF HEADER'
    do while (.true.)
      backspace tmphtg
      read (tmphtg) nme, mw, dummy, gni, gnx, gei, gex
      if (nme .eq. '0000') exit
      call mjd2date(0, mw, iy, imon, id, ih, im, xini)
      write (lfnhtg, '(a4,i5,4i3,f10.6,$)') nme, iy, imon, id, ih, im, xini
      call mjd2date(0, dummy, iy, imon, id, ih, im, xini)
      write (lfnhtg, '(i5,4i3,f10.6,4f11.6)') iy, imon, id, ih, im, xini, gni, gnx, gei, gex
      backspace tmphtg
    enddo
    close (lfnhtg); close (tmphtg, status='delete')
  endif
!
!! write kinematic files
  if (index(SITE%skd, 'K') .ne. 0 .and. SITE%skd(2:2) .ne. 'F') then ! kinematic & pseudo-kinematic
    if (SITE%ikin .ne. 0) close (SITE%ikin)
    SITE%ikin = get_valid_unit(10)
    open (SITE%ikin, file=SITE%kinfil)
    write (SITE%ikin, '(a20,10x,a4,26x,a)') 'Kinematic Trajectory', SITE%name, 'COMMENT'
    write (SITE%ikin, '(f9.2,51x,a)') LCF%dintv, 'INTERVAL'
    write (SITE%ikin, '(60x,a)') 'END OF HEADER'
  endif
  do while (.true.)
    backspace tmpkin
    read (tmpkin) isit, jd, mw, (x(i), i=1, 3), iobs
    if (isit .eq. 0) exit
    cid = ' '
    if (iobs .lt. 5) cid = '* '
    write (SITE%ikin, '(i5,f9.2,a2,3f13.3)') jd, mw, cid, (x(i), i=1, 3)
    backspace tmpkin
  enddo
  close (tmpkin, status='delete')
  if (index(SITE%skd, 'K') .ne. 0) close (SITE%ikin)
  write (*, '(a,3i6)') 'Recovering: nc, np, ns ', NM%nc, NM%np, NM%ns

  return
100 write (*, '(a)') '***ERROR(lsq_rcv_prmt): back space or read file'
  call exit(1)
end
