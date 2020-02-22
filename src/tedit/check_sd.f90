!
!! check_sd.f90
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
!! purpose  : check single-difference LC
!! parameter:
!!    input : lfnsd --   0 : debug-sd = true
!!                      -1 : debug-sd = false
!!            ndgr  --   2
!!            niter --   3
!!            mepo  --   20 :arc length for checking
!!            nstep --   8  :step for next checking
!!          lclimit --
!!
!
subroutine check_sd(lfnsd, neph, ephem, ndgr, niter, mepo, nstep, x, y, z, interval, lclimit)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'
  include 'data_flag.h'
!
!! common part
  integer*4 nepo, nsat, flagall(MAXEPO, MAXSAT), jd0, nobs(MAXSAT)
  real*8 ti(MAXEPO), ts(MAXEPO), obs(MAXEPO, MAXSAT, 6)
  common/OBSDAT/nepo, nsat, jd0, nobs, flagall, ti, ts, obs
!
  integer*4 neph, lfnsd, ndgr, niter, mepo, nstep
  type(brdeph) EPHEM(MAXEPH)
  real*8 interval, lclimit, x, y, z
!
!! local
  logical*1 found, ref_ok, lwrite
  integer*4 i, j, k, i0, i1, istart, istop, ngood, iepo, kobs, ierr, jeph, jepo, nsd, iprn, iref, njump
  integer*4 ichecked(MAXEPO, MAXSAT), jump(MAXEPO*10, 4), nobs_epoch(MAXEPO), nobsprn(MAXSAT), nobsflg(MAXSAT)
  integer*4 flglg(MAXEPO), flglc(MAXEPO), flglgref(MAXEPO), ilc(MAXEPO), ilg(MAXEPO)
  real*8 lc(MAXEPO), lg(MAXEPO), range_ref(MAXEPO), range(MAXEPO)
  real*8 vlc(MAXEPO), vlg(MAXEPO), vlgref(MAXEPO)
  real*8 rmslg(MAXSAT), rmslc(MAXSAT)
  real*8 lambda1, dtsat, elev, dummy, a0
!
!! function called
  logical*1 istrue
  integer*4 set_flag
!
!! initialization
  lwrite = .true.
  if (lfnsd .eq. 0) lfnsd = 6
  if (lfnsd .lt. 0) lwrite = .false.
  lambda1 = VLIGHT/FREQ1
  ierr = 0
!! find out first and last epoch with enough data
!! nobs_epoch -- # of obs. per epoch, the number of available satellite
  istart = 0     ! start epoch
  istop = 0     ! end epoch
  ngood = 0     ! Effective number of epoch
  do iepo = 1, nepo
    nobs_epoch(iepo) = 0
    do iprn = 1, nsat
      ichecked(iepo, iprn) = 10
      if (istrue(flagall(iepo, iprn), 'ok')) then
        nobs_epoch(iepo) = nobs_epoch(iepo) + 1    ! the number of available satellite
      endif
    enddo
    if (nobs_epoch(iepo) .ge. 2) then
      if (istart .eq. 0) istart = iepo
      istop = iepo
      ngood = ngood + 1
    endif
  enddo
  if (istop - istart + 1 .le. ndgr + 2 .or. ngood .le. ndgr + 2) then
! remove all
    return
  endif

  njump = 0
  iepo = istart
  do while (iepo .le. istop)
    ngood = 0

! extract an epoch section from iepo
! mepo = 20; default value, the constant length of check arc is 20 epochs
    i1 = min(iepo + mepo - 1, istop)   ! set the start epoch of check arc
    i0 = max(istart, i1 - mepo + 1)    ! set the end epoch of check arc

    do k = i0, i1
      if (nobs_epoch(k) .ge. 2) ngood = ngood + 1
    enddo

! ndgr = 2; default value
    if (ngood .le. ndgr + 2) goto 100

! kobs ------ the number of epoch in check arc
    kobs = i1 - i0 + 1

! poly fit for each satellite in epoch section to choose reference satellite
!-nobsprn : # of observations for each satellite
!-nobsflg : # of flags (amb.) for each satellite
    do iprn = 1, nsat
      nobsprn(iprn) = 0
      nobsflg(iprn) = 0
      do jepo = 1, kobs
        k = i0 + jepo - 1
        if (istrue(flagall(k, iprn), 'ok')) then
          nobsprn(iprn) = nobsprn(iprn) + 1     ! the number of available epochs for each satellite
        endif
      enddo
      if (nobsprn(iprn) .gt. ndgr + 2) then
        do jepo = 1, kobs
          k = i0 + jepo - 1
          flglg(jepo) = 10
          if (istrue(flagall(k, iprn), 'ok')) then
            flglg(jepo) = 1
            lg(jepo) = obs(k, iprn, 1)          ! geometry-free obs value
          endif
          if (istrue(flagall(k, iprn), 'good')) then
            flglg(jepo) = 0
          endif
        enddo
        call check_for_jump(' LG ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                            2.0d0, 0.40d0, a0, rmslg(iprn), vlg, ilg, ierr, interval)
        if (ierr .ne. 0) then
!        write(*,'(a,2i6,i3)') ' Error check LG',i0,i1,iprn
        endif
        do jepo = 1, kobs
          k = i0 + jepo - 1
          if (flglg(jepo) .eq. 1) nobsflg(iprn) = nobsflg(iprn) + 1
        enddo
        if (ierr .ne. 0) rmslg(iprn) = 1.d3
      else
        nobsflg(iprn) = ngood
        rmslg(iprn) = 1.d3
      endif
    enddo
!!
!! select satisfactory referece satellte
    iref = 1
    do iprn = 1, nsat
      if (nobsprn(iprn) .gt. nobsprn(iref)) then
        iref = iprn
      else if (nobsprn(iprn) .eq. nobsprn(iref)) then
        if (nobsflg(iprn) .lt. nobsflg(iref)) then
          iref = iprn
        else if (nobsflg(iprn) .eq. nobsflg(iref)) then
          if (rmslg(iprn) .le. rmslg(iref)) iref = iprn
        endif
      endif
    enddo
!!
!! check referece satellite
!!-range_ref : distance from ref. satellites to static station in each epoch
    jeph = 0
    do jepo = 1, kobs
      k = i0 + jepo - 1
      range_ref(jepo) = 0.d0
      flglg(jepo) = 10
      if (istrue(flagall(k, iref), 'ok')) then
        call elevation(neph, ephem, iref, jd0, ti(k), x, y, z, elev, range_ref(jepo), dtsat, dummy, jeph, .true.)
        range_ref(jepo) = range_ref(jepo) - dtsat
        flglg(jepo) = 1
        if (istrue(flagall(k, iref), 'good')) flglg(jepo) = 0
        lg(jepo) = obs(k, iref, 1)
      endif
    enddo
    call check_for_jump(' LGR ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                        2.0d0, 0.40d0, a0, rmslg(iref), vlg, ilg, ierr, interval)
    if (ierr .ne. 0) then
!    write(*,'(a,i6,2i3)') ' Error check LGR ',i0,i1,iref
    endif
    do jepo = 1, kobs
      k = i0 + jepo - 1
      flglgref(jepo) = flglg(jepo)
      vlgref(jepo) = 0.d0
      if (ilg(jepo) .ne. 0) then
        vlgref(jepo) = vlgref(ilg(jepo)) + vlg(jepo)
      endif
    enddo
!!
!! check other satellites, form single diff LC obs
    do iprn = 1, nsat
      if (nobsprn(iprn) .gt. ndgr + 2 .and. iprn .ne. iref) then
        jeph = 0
        nsd = 0
        do jepo = 1, kobs
          range(jepo) = 0.d0
          k = i0 + jepo - 1
          flglc(jepo) = 10
          if (istrue(flagall(k, iprn), 'ok') .and. istrue(flagall(k, iref), 'ok')) then
            call elevation(neph, ephem, iprn, jd0, ti(k), x, y, z, elev, range(jepo), dtsat, dummy, jeph, .true.)
            range(jepo) = range(jepo) - dtsat
            ! LC single-difference = troposphere delay or ...
            lc(jepo) = obs(k, iprn, 3) - obs(k, iref, 3) - (range(jepo) - range_ref(jepo))/lambda1
            if (lwrite) write (lfnsd, '(i6,f15.3)') k, lc(jepo)
            flglc(jepo) = 1
            if (istrue(flagall(k, iprn), 'good') .and. istrue(flagall(k, iref), 'good')) then
              nsd = nsd + 1
              flglc(jepo) = 0
            endif
          endif
        enddo
        if (nsd .gt. ndgr + 2) then
          if (lwrite) write (lfnsd, *) ' Satellite ', iprn, iref
          call check_for_jump(' LC ', lfnsd, kobs, ti(i0), lc, flglc, ndgr, niter, &
                              0.2d0, 0.12d0, a0, rmslc(iprn), vlc, ilc, ierr, interval)
          if (ierr .eq. 0) then
            if (rmslc(iprn) .gt. lclimit .and. lwrite) &
              write (lfnsd, '(a,f10.6)') ' LC big rms ', rmslc(iprn)
            do jepo = 1, kobs
              if (ilc(jepo) .ne. 0) then
                k = i0 + jepo - 1
                if (ichecked(k, iprn) .ne. 0) ichecked(k, iprn) = 0
                if (ichecked(k, iref) .ne. 0) ichecked(k, iref) = 0
                if (ichecked(i0 + ilc(jepo) - 1, iprn) .eq. 10) &
                  ichecked(i0 + ilc(jepo) - 1, iprn) = 2
                if (ichecked(i0 + ilc(jepo) - 1, iref) .eq. 10) &
                  ichecked(i0 + ilc(jepo) - 1, iref) = 2
                if (dabs(vlc(jepo)) .gt. lclimit) then
                  if (lwrite) write (lfnsd, *) ' Found ', k, jepo, iprn, iref, vlc(jepo)
                  if (flglgref(jepo) .ne. 1 .and. rmslg(iref) .le. 0.10) ref_ok = .true.
                  found = .false.
                  do i = 1, njump
                    if (jump(i, 1) .eq. k .and. (jump(i, 2) .eq. iprn .or. &
                                                 jump(i, 3) .eq. iprn) .and. (jump(i, 2) .eq. iref .or. &
                                                                              jump(i, 3) .eq. iref)) then
                      found = .true.
                      cycle
                    endif
                  enddo
                  if (.not. found) then
                    njump = njump + 1
                    jump(njump, 1) = k
                    jump(njump, 2) = iprn
                    jump(njump, 3) = iref
                    jump(njump, 4) = 2
                    if (ref_ok .and. istrue(flagall(k - 1, iref), 'ok')) jump(njump, 4) = 1
                  endif
                endif
              endif
            enddo
          else
            if (lwrite) write (lfnsd, *) ' fit error, can not check this piece ', iprn, i0, i1, ierr, rmslc(iprn)
!          write(*,*) ' fit error, can not check this piece ',iprn,i0,i1,ierr,rmslc(iprn)
          endif
        endif
      endif
! next satellite
    enddo
100 continue
! next epoch
    iepo = iepo + nstep
  enddo

!!
!! count and remove unchecked data
  j = 0
  do iepo = 1, nepo
    do iprn = 1, nsat
      if (ichecked(iepo, iprn) .eq. 2) then
        if (lwrite) write (lfnsd, *) ' First Epoch ', iepo, iprn
        flagall(iepo, iprn) = set_flag(flagall(iepo, iprn), 'bigsd')
      else if (istrue(flagall(iepo, iprn), 'ok') .and. ichecked(iepo, iprn) .ne. 0) then
! shift flag to the next epoch
        if (.not. istrue(flagall(iepo, iprn), 'good')) then
          do k = iepo + 1, nepo
            if (istrue(flagall(k, iprn), 'ok')) then
              if (istrue(flagall(k, iprn), 'good')) flagall(k, iprn) = flagall(iepo, iprn)
              exit
            endif
          enddo
          j = j + 1
          flagall(iepo, iprn) = set_flag(flagall(iepo, iprn), 'lccheck')
          if (lwrite) write (lfnsd, *) ' Not checked ', iepo, iprn
        endif
      endif
    enddo
  enddo
!write(*,*) ' Not checked ',j
  if (lwrite) write (lfnsd, *) ' Not checked', j

! statistics of flags and set in 'flagall'
  j = 0
  do i = 1, njump
    if (istrue(flagall(jump(i, 1), jump(i, 2)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 2)), 'good')) then
        j = j + 1
      endif
!      write(*,'(a,2i6)') 'LC flagged ',jump(i,1),jump(i,2)
      flagall(jump(i, 1), jump(i, 2)) = set_flag(flagall(jump(i, 1), jump(i, 2)), 'bigsd')
    endif
    if (jump(i, 4) .eq. 2 .and. istrue(flagall(jump(i, 1), jump(i, 3)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 3)), 'good')) then
        j = j + 1
      endif
!      write(*,'(a,2i6)') 'LC flagged ',jump(i,1),jump(i,3)
      flagall(jump(i, 1), jump(i, 3)) = set_flag(flagall(jump(i, 1), jump(i, 3)), 'bigsd')
    endif
  enddo
  if (lwrite) write (lfnsd, *) ' TOTAL NEW FLAGs', j
!write(*,*) ' TOTAL NEW FLAGs',j

  return
end
