!
!! write_diag_rpt.f90
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
!! purpose  : write rinex health diagnose report
!! parameter:
!!    input : flnrhd   -- diagnose file
!!            nepo     -- # of epochs
!!            nsat     -- # of satellites
!!            jd0,ti   -- time tag
!!            flagall  -- flags
!!            interval -- sampling interval
!
subroutine write_diag_rpt(flnrhd, nepo, nsat, jd0, ti, flagall, interval, sumepo)
  implicit none
  include '../header/const.h'

  integer*4 nepo, sumepo, nsat, jd0, flagall(MAXEPO, MAXSAT)
  real*8 ti(1:*), interval
  character*(*) flnrhd
!
!! local
  integer*4 i, j, ns, lfn, iepo, iprn, ie, jd, iy, imon, id, ih, im
  integer*4 maxamb, totamb, remobs, avaobs, xamb(MAXEPO), flag(MAXSAT, MAXEPO)
  real*8 sec, sod
  character*256 string, str(MAXSAT)
!
!! function called
  logical*1 istrue
  integer*4 get_valid_unit, find_ambd
!
!! initialization
  maxamb = 0
  totamb = 0
  remobs = 0
  avaobs = 0
  flag = transpose(flagall)
  do i = 1, nepo
    xamb(i) = 0
  enddo
!
!! statistics about rinex health
  do iprn = 1, nsat
    do iepo = 1, nepo
      if (istrue(flagall(iepo, iprn), 'nodata')) cycle
      if (istrue(flagall(iepo, iprn), 'ok')) then
        avaobs = avaobs + 1
        if (istrue(flagall(iepo, iprn), 'amb')) then
          totamb = totamb + 1
          j = find_ambd(nepo, flagall(1, iprn), iepo)
          do i = iepo, j
            xamb(i) = xamb(i) + 1
          enddo
        endif
      else
        remobs = remobs + 1
      endif
    enddo
  enddo
  maxamb = maxval(xamb, 1)
!
!! write header
  lfn = get_valid_unit(10)
  open (lfn, file=flnrhd)
  write (lfn, '(a21,39x,a)') 'Rinex Health Diagnose', 'COMMENT'
  write (lfn, '(2f10.2,40x,a)') interval, interval, 'INT AMB/DEL'
  write (lfn, '(3i10,30x,a)') maxamb, totamb, 0, 'AMB MAX/TOT/NEW'
  write (lfn, '(3i10,30x,a)') avaobs, remobs, 0, 'EPO AVA/REM/NEW'
  write (lfn, '(3i10,30x,a)') sumepo, nepo, 0, 'EFF EPO/SUM/NEW'
  write (lfn, '(60x,a)') 'END OF HEADER'
!
!! epoch diagnose
  do iepo = 1, nepo
    ns = 0
    do iprn = 1, nsat
      if (istrue(flag(iprn, iepo), 'nodata') .or. istrue(flag(iprn, iepo), 'good')) cycle
      if (.not. istrue(flag(iprn, iepo), 'ok')) then
        ns = ns + 1
        string = ' '
        if (istrue(flag(iprn, iepo), 'lowele')) then
          string = 'LOWELEVATION'
        else if (istrue(flag(iprn, iepo), 'shrt')) then
          string = 'SHORTPIECE'
        else if (istrue(flag(iprn, iepo), 'lwbad')) then
          string = 'BADWIDELANE'
        else if (istrue(flag(iprn, iepo), 'lgbad')) then
          string = 'BADIONOSPHERE'
        else if (istrue(flag(iprn, iepo), 'lccheck')) then
          string = 'CANNOTCHECKLC'
        else if (istrue(flag(iprn, iepo), 'no4')) then
          string = 'LESSTHAN4OBS'
        else if (istrue(flag(iprn, iepo), 'pcbad')) then
          string = 'BADRANGE'
        else if (istrue(flag(iprn, iepo), 'pc1ms')) then
          string = 'BAD1MS'
        endif
        write (str(ns), '(1x,i2,57x,a)') iprn, 'DEL_'//trim(string)
      else if (istrue(flag(iprn, iepo), 'amb')) then
        ns = ns + 1
        string = ' '
        if (istrue(flag(iprn, iepo), 'lli')) then
          string = 'FLAGINRINEX'
        else if (istrue(flag(iprn, iepo), 'bigsd')) then
          string = 'BIGLCJUMP'
        else if (istrue(flag(iprn, iepo), 'gap')) then
          string = 'BIGGAP'
        else if (istrue(flag(iprn, iepo), 'lwjump')) then
          string = 'WIDELANEJUMP'
        else if (istrue(flag(iprn, iepo), 'lwconn')) then
          string = 'IONO.FAILED'
        else if (istrue(flag(iprn, iepo), 'lgjump')) then
          string = 'IONO.JUMP'
        endif
        ie = find_ambd(nepo, flagall(1, iprn), iepo)
        jd = jd0 + int(ti(ie)/86400.d0)
        sod = ti(ie) - (jd - jd0)*86400.d0
        call mjd2date(jd, sod, iy, imon, id, ih, im, sec)
        write (str(ns), '(1x,i2,28x,i5,4i3,f11.7,1x,a)') iprn, iy, imon, id, ih, im, sec, &
          'AMB_'//trim(string)
      endif
    enddo
    if (ns .eq. 0) cycle
    jd = jd0 + int(ti(iepo)/86400.d0)
    sod = ti(iepo) - (jd - jd0)*86400.d0
    call mjd2date(jd, sod, iy, imon, id, ih, im, sec)
    write (lfn, '(a3,i5,4i3,f11.7)') 'TIM', iy, imon, id, ih, im, sec
    do i = 1, ns
      write (lfn, '(a)') str(i) (1:len_trim(str(i)))
    enddo
  enddo

  close (lfn)
  return
end
