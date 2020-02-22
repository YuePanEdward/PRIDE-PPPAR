!
!! read_rinex_file.f90
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
!! flnrnx ------------ rinex O_file
! tstart,sstart ----- start time: tstart(1):year, tstart(2):month, tstart(3):day, tstart(4):hour, tstart(5):minute, sstart:second.
! sesstion_length---- sesstion time length
! interval ---------- Epoch interval of processing. (rinex o_file Epoch interval: OH.intv)
! use_rinex_flag
! check_pc
! pclimit -----------
! cutoff_elevation -- cut off elevation
! use_brdeph -------- true or false
! neph -------------- the number of eph
! ephem -------------
! stanam ------------
! x,y,z -------------
! t_first_in_rinex --
! t_last_in_rinex ---
! v -----------------
subroutine read_rinex_file(flnrnx, tstart, sstart, session_length, interval, &
                           use_rinex_flag, check_pc, pclimit, cutoff_elevation, use_brdeph, &
                           neph, ephem, stanam, x, y, z, t_first_in_rinex, t_last_in_rinex, v, sumepo)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'
  include '../header/brdeph.h'

  logical*1 use_rinex_flag, check_pc, use_brdeph
  integer*4 neph, tstart(1:*), session_length
  real*4 v(MAXEPO, MAXSAT)
  real*8 sstart, x, y, z, interval, lglimit, nwlimit, cutoff_elevation, pclimit
  character*(*) flnrnx, stanam
  type(brdeph) EPHEM(1:*)

  integer*4 nepo, sumepo, nsat, flagall(MAXEPO, MAXSAT), jd0, nobs(MAXSAT)
  real*8 ti(MAXEPO), ts(MAXEPO), obs(MAXEPO, MAXSAT, 6), dcb(MAXSAT, 2), bias(MAXSAT, 4)
  common/OBSDAT/nepo, nsat, jd0, nobs, flagall, ti, ts, obs, dcb, bias

! local
  integer*4 nprn, prn(MAXSAT), iepo, iprn, j, k, iunit, ichn, ilast(MAXSAT), ieph, nused
  integer*4 ierr
  real*8 tobs, nwdif, lgdif, pgdif, elev, range, dtsat, dt, dwnd
  real*8 g, lambda1, lambda2, lambdaw, c1, c2
  real*8 t_first_in_rinex, t_last_in_rinex, vv
  type(rnxhdr) HD
  type(rnxobr) OB

! for range check
  integer*4 flg_pmb(MAXSAT)
  real*8 pmb(MAXSAT)

! receiver clock jump check
  real*8        jumpsum(2), jumpval(2), rangeval(2, MAXSAT), deltap, deltal
  real*8        tmp(2), ratio, sec, obsval(2, 4, MAXSAT), chkval(2, MAXSAT)
  integer*4     lfnjmp, nvalid, njump, cnt, jumpflag  ! 1:range, 2:phase

! funtion called
  logical*1 istrue
  integer*4 modified_julday, set_flag, get_valid_unit
  real*8 timdif

  g = 7.7d0/6.0d0
  c1 = g*g/(g*g - 1.d0)
  c2 = -1.d0/(g*g - 1.d0)
  lambda1 = VLIGHT/FREQ1
  lambda2 = VLIGHT/FREQ2
  lambdaw = VLIGHT/(FREQ1 - FREQ2)

  lglimit = 50.d0
  nwlimit = 50.d0

  lfnjmp = 0
  jumpsum = 0.d0
  obsval = 0.d0

  dcb = 0.d0
  bias = 0.d0

  dwnd = min(interval/10.d0, 0.3d0)
! open rinex observation file
  iunit = 10
  open (iunit, file=flnrnx, form='FORMATTED', status='OLD', iostat=ierr)
  if (ierr .ne. 0) then
!  write(*,'(a)') '***ERROR(read_rinex_file): open file error, '//flnrnx
    call exit(1)
  endif

! read head of rinex file
  call rdrnxoh(iunit, HD, ierr)
  if (ierr .ne. 0) then
!  write(*,'(a,i3)') '***ERROR(read_rinex_file): read rinex header error ',ierr
    call exit(1)
  endif

! coordinates of station
  if (use_brdeph) then
    if (x*y*z .eq. 0.d0) then
      if (HD%x*HD%y*HD%z .eq. 0.d0) then
!      write(*,'(a)') '***ERROR(read_rinex_file): no position when using brdeph '
        call exit(1)
      endif
      x = HD%x
      y = HD%y
      z = HD%z
    endif
  endif

! start time of first session
  if (tstart(2) .eq. 0) then
    do j = 1, 5
      tstart(j) = HD%t0(j)
    enddo
    sstart = HD%t0s
  endif
  call yr2year(tstart(1))

!! jd0:  julday
!! tobs: fractional julday
  jd0 = modified_julday(tstart(3), tstart(2), tstart(1))
  tobs = tstart(4)*3600.d0 + tstart(5)*60.d0 + sstart

! initialize
! nepo: the numbers of processing epoch
  nepo = min(nint(session_length/interval), MAXEPO)
  nsat = 0
  do iprn = 1, MAXSAT
    ilast(iprn) = 0
  enddo

! ti --- time of each processing epoch
! nepo - the number of epoch
  ti(1) = tobs
  do iepo = 1, nepo
    if (iepo .ge. 2) ti(iepo) = ti(iepo - 1) + interval
    ts(iepo) = ti(iepo)
    do iprn = 1, MAXSAT
      flagall(iepo, iprn) = set_flag(0, 'nodata')
      v(iepo, iprn) = 0.d0
      do j = 1, 6
        obs(iepo, iprn, j) = 0.d0
      enddo
    enddo
  enddo
!
! loop over session and epoch
  iepo = 1

  do while (iepo .le. nepo .and. ierr .eq. 0)
    nprn = 0
    prn = 0
! ***************************************************************************************** !
! purpose  : read one epoch data from a RINEX o-file
! parameter: jd0,tobs ----------- start time
!            dwnd --------------- window for time matching. If the obsersing time from the
!                                 rinex-file is close to the requested time within the
!                                 window, we take the data. Be careful with this parameter
!                                 when you are working sampling rate larger than 1Hz.
!            nprn,prn ----------- satellite prn number
!            HD ----------------- rinex header structure
!            OB ----------------- o_file body structure
! **************************************************************************************** !
    if (HD%ver .eq. 3) then
      call rdrnxoi3(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, dcb, bias, ierr)
    else
      call rdrnxoi(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, dcb, bias, ierr)
    endif

    if (ierr .eq. 0) then
      sumepo = sumepo + 1
      obsval(1,:,:) = obsval(2,:,:)
      obsval(2,:,:) = 0.d0
      ti(iepo) = tobs
      if (OB%nprn .ne. 0) ti(iepo) = ti(iepo) - timdif(jd0, tobs, OB%jd, OB%tsec)
      do ichn = 1, OB%nprn
        j = OB%prn(ichn)
! 8/5/2006.Geng.for GLONASS
        if (j .eq. 0) cycle
! nsat : largest PRN for a satellite
        if (nsat .lt. j) nsat = j
! 12/17/2006. Geng. remove satellite without broadcast information
        if (use_brdeph) then
          ieph = 0
! range --------- the distance between satellite and reciver
! dtsat --------- clock correction, units: meter
! vv    --------- true anomaly
          call elevation(neph, ephem, j, jd0, ti(iepo), x, y, z, elev, range, dtsat, vv, ieph, .false.)
          if (range .lt. 0.d0) then
            flagall(iepo, j) = set_flag(0, 'no4')
            cycle
          endif
        endif
        if (dabs(OB%obs(ichn, 1)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 2)) .gt. 1.d-3 .and. &
            dabs(OB%obs(ichn, 3)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 4)) .gt. 1.d-3) then
          flagall(iepo, j) = 0
          obsval(2,1,j) = OB%obs(ichn, 3) ! P1
          obsval(2,3,j) = OB%obs(ichn, 4) ! P2
          obsval(2,2,j) = OB%obs(ichn, 1)*lambda1 ! L1
          obsval(2,4,j) = OB%obs(ichn, 2)*lambda2 ! L2
          rangeval(1,j) = rangeval(2,j)
          rangeval(2,j) = range
          ! geometry-free : ionosphere observations
          obs(iepo, j, 1) = (lambda1*OB%obs(ichn, 1) - lambda2*OB%obs(ichn, 2))/(lambda2 - lambda1)
          ! **************************************************************** !
          !                       check geometry-free(lg)                    !
          ! **************************************************************** !
          if (ilast(j) .ne. 0) then
            lgdif = dabs(obs(iepo, j, 1) - obs(ilast(j), j, 1))/(ti(iepo) - ti(ilast(j)))
            if (lgdif .gt. lglimit) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lgjump')
!            write(*,'(a,i6,i4,f12.1)') ' bad ionosphere',iepo,j,lgdif
            endif
          endif
        else
          flagall(iepo, j) = set_flag(0, 'no4')
        endif
        if (use_rinex_flag .and. (mod(OB%lli(ichn, 1), 2) .ne. 0 .or. mod(OB%lli(ichn, 2), 2) .ne. 0)) then
          flagall(iepo, j) = set_flag(flagall(iepo, j), 'lli')
!        write(*,'(a,i6,i4,f12.1)') ' flag in rinex ',iepo,j
        endif
      enddo

      ! **************************************************************** !
      !               check & recover receiver clock jump                !
      !              Author:       Yuanxin Pan, 2019-07-15               !
      ! **************************************************************** !
      nvalid = 0 ! number of satellites with flag 'ok'
      njump = 0
      cnt = 0
      jumpval = 0.d0
      deltap = 0.d0
      deltal = 0.d0
      do ichn = 1, OB%nprn
        j = OB%prn(ichn)
        if (j .eq. 0) cycle
        if (.not.istrue(flagall(iepo,j), 'ok')) cycle
        if (dabs(obsval(1, 1, j)).lt.1.d-3 .or. dabs(obsval(2, 1, j)).lt.1.d-3) cycle
        nvalid = nvalid + 1
        chkval(1,j) = obsval(2,1,j)-obsval(1,1,j) - (obsval(2,2,j)-obsval(1,2,j)) ! f1
        chkval(2,j) = obsval(2,3,j)-obsval(1,3,j) - (obsval(2,4,j)-obsval(1,4,j)) ! f2
        if (dabs(chkval(1, j)).gt.20.d0 .and. dabs(chkval(2, j)).gt.20.d0) then ! 0.1 us jump: 30.d0
          njump = njump + 1
          if (rangeval(1, j).gt.1.d0 .and. rangeval(2, j).gt.1.d0) then
            cnt = cnt + 1
            tmp(1) = (obsval(2,1,j)-obsval(1,1,j) + obsval(2,3,j)-obsval(1,3,j))/2.d0 - (rangeval(2,j) - rangeval(1,j))
            tmp(2) = (obsval(2,2,j)-obsval(1,2,j) + obsval(2,4,j)-obsval(1,4,j))/2.d0 - (rangeval(2,j) - rangeval(1,j))
            deltap = deltap + tmp(1)
            deltal = deltal + tmp(2)
          endif
          jumpval(1) = jumpval(1) + chkval(1,j)
          jumpval(2) = jumpval(2) + chkval(2,j)
        endif
      enddo

      if (nvalid.ne.0 .and. nvalid.eq.njump) then
        jumpval = jumpval/(nvalid)
        deltap = deltap/cnt  ! P: range jump value
        deltal = deltal/cnt  ! L: phase jump value
        ratio = dabs(deltap/deltal)
        !sec = anint((jumpval(1)+jumpval(2))/2/vlight*1.d3) ! unit: ms
        sec = (jumpval(1)+jumpval(2))/2.d0  ! clk jmp unit: m
        ! inverse fix
        if(ratio < 1.d0) then ! phase jump
          jumpflag = 1
          jumpsum(1) = jumpsum(1) + sec!/1.d3*vlight ! range
        else
          jumpflag = 2
          jumpsum(2) = jumpsum(2) + sec!/1.d3*vlight ! phase
        end if
        !write(lfnjmp,'(f8.1,x,2f8.3,x,f20.1,x,i1,x,f28.14)') cd.ti(iepo), deltap/vlight*1.d3, deltal/vlig
        if (lfnjmp.eq.0) then
          lfnjmp = get_valid_unit(10)
          open(lfnjmp,file='.'//stanam//'.jmp',status='replace')
        end if
        write(lfnjmp,'(f8.1,2x,i1,x,f28.14)') ti(iepo), jumpflag, sec
      end if

      do ichn = 1, OB%nprn
        j = OB%prn(ichn)
        if (j .eq. 0) cycle
        if (nsat .lt. j) nsat = j
        if (use_brdeph) then
          ieph = 0
          call elevation(neph, ephem, j, jd0, ti(iepo), x, y, z, elev, range, dtsat, vv, ieph, .false.)
          if (range .lt. 0.d0) then
            flagall(iepo, j) = set_flag(0, 'no4')
            cycle
          endif
        endif
        if (dabs(OB%obs(ichn, 1)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 2)) .gt. 1.d-3 .and. &
            dabs(OB%obs(ichn, 3)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 4)) .gt. 1.d-3) then
          OB%obs(ichn, 3:4) = OB%obs(ichn, 3:4) + jumpsum(1) ! P1 & P2
          OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1
          OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2
          ! Melbourne-Wubbena (N1 - N2)
          obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g)/lambdaw
          ! ionosphere-free (pp37 LC)
          obs(iepo, j, 3) = c1*OB%obs(ichn, 1) + g*c2*OB%obs(ichn, 2)
          if (use_brdeph) then
            ieph = 0
            call elevation(neph, ephem, j, jd0, ti(iepo), x, y, z, elev, range, dtsat, vv, ieph, .false.)
            v(iepo, j) = vv/PI*180.0
            ! (sit-sat distance from PC)-(sit-sat distance from broadcast) to check recv clock
            obs(iepo, j, 6) = c1*OB%obs(ichn, 3) + c2*OB%obs(ichn, 4) - (range - dtsat)
            ! ******************************************* !
            !              check elevation                !
            ! ******************************************* !
            if (elev .lt. cutoff_elevation) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lowele')
            endif
          endif
          ! **************************************************************** !
          !                        check Melbourne-Wubbenamw                 !
          ! **************************************************************** !
          if (ilast(j) .ne. 0) then
            nwdif = dabs(obs(iepo, j, 2) - obs(ilast(j), j, 2))
            if (nwdif .gt. nwlimit) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lwjump')
  !            write(*,'(a,i6,i4,4f15.1)') ' bad widelane  ',iepo,j,nwdif
            endif
          endif
          ilast(j) = iepo
        else
          flagall(iepo, j) = set_flag(0, 'no4')
        endif
      enddo
      iepo = iepo + 1
      tobs = tobs + interval
    endif
  enddo
  close (10)
  if (lfnjmp .ne. 0) close(lfnjmp)
!
!! check receiver clock
  iepo = 1
  do while (iepo .le. nepo)
    nused = 0
    do j = 1, MAXSAT
      if (istrue(flagall(iepo, j), 'ok')) nused = nused + 1
    enddo
    if (nused .ne. 0) then
      if (check_pc) then
        do j = 1, MAXSAT
          flg_pmb(j) = 0
          pmb(j) = obs(iepo, j, 6)
          if (pmb(j) .eq. 0.d0) flg_pmb(j) = 2   ! '2' means no data for satellite
        enddo
!!    dt -------- mean value of recv clock
!     one epoch for check_range
        call check_range(iepo, MAXSAT, dt, pmb, flg_pmb)
        do j = 1, MAXSAT
          if (flg_pmb(j) .le. 1) then
            obs(iepo, j, 6) = pmb(j) - dt
          endif
! 245m for cleanresid 250m
          if (flg_pmb(j) .eq. 1 .or. dabs(obs(iepo, j, 6)) .gt. pclimit) then
            if (istrue(flagall(iepo, j), 'amb')) then
              call find_flag(iepo + 1, nepo, flagall(1, j), 'ok', k)
              if (k .gt. 0) flagall(k, j) = flagall(iepo, j)
            endif
            if (dabs(obs(iepo, j, 6)) .gt. 2*1.d5) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'pc1ms')
            else
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'pcbad')
            endif
!          write(*,'(a,i6,i4,3f12.1)') '*flag bad pc ',iepo,j,pmb(j)-dt,pmb(j),dt
          endif
        enddo
      endif
!! consider receiver clock
      ti(iepo) = ti(iepo) - dt/VLIGHT
    endif
    iepo = iepo + 1
  enddo

! time tag in second
  j = modified_julday(HD%t0(3), HD%t0(2), HD%t0(1))
  t_first_in_rinex = HD%t0(4)*3600.d0 + HD%t0(5)*60.d0 + HD%t0s + (j - jd0)*86400.d0
  t_last_in_rinex = ti(nepo)

  return
end
