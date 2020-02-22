!
!! tedit.f90
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
!! Turboedit
!
program tedit
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'
  include 'data_flag.h'

! for the whole file, common
  integer*4 nepo, nsat, flagall(MAXEPO, MAXSAT), jd0, nobs(MAXSAT)
  real*8 ti(MAXEPO), ts(MAXEPO), obs(MAXEPO, MAXSAT, 6), dcb(MAXSAT, 2), bias(MAXSAT, 4)
  common/OBSDAT/nepo, nsat, jd0, nobs, flagall, ti, ts, obs, dcb, bias

! broadcast ephemeris
  integer*4 neph, sumepo
  type(brdeph) ephem(MAXEPH)

  integer*4 lfnsd
  logical*1 again, use_brdeph, use_rinex_flag, check_pc, check_lc, turbo_edit, &
    debug_sd, debug_tb, keep_end
! control parameters
  character*256 flnrnx, flneph, flnrhd, string*20, stanam*4
  integer*4 tstart(5), session_length, length_gap, length_short
  real*8 sstart, cutoff_elevation, interval, pclimit, lclimit, lglimit, lgrmslimit, &
    max_mean_namb, min_percent, min_mean_nprn
! local
  integer*4 iprn, iepo, nused
  real*8 x, y, z, t_first_in_rinex, t_last_in_rinex
  real*4 v(MAXEPO, MAXSAT)
! for one or single difference
!integer*4 flglc,flglg,flgnw
  real*8 pg(MAXEPO), sigpg(MAXEPO), respg(MAXEPO)
!! function used
  logical*1 istrue
!
!! get input from command line
  call get_control_parameter(flnrnx, flneph, flnrhd, check_lc, &
                             turbo_edit, use_brdeph, use_rinex_flag, check_pc, keep_end, &
                             tstart, sstart, session_length, length_gap, length_short, cutoff_elevation, &
                             max_mean_namb, min_percent, min_mean_nprn, interval, lclimit, pclimit, &
                             lglimit, lgrmslimit, stanam, x, y, z)
  debug_tb = .true.
  debug_sd = .false.

! read broadcast ephemeris
  if (use_brdeph) call rdrnxn(flneph, 0.d0, 0.d0, neph, ephem)

! read rinex observation file
  sumepo = 0
  dcb = 0.d0
  bias = 0.d0
  call read_rinex_file(flnrnx, tstart, sstart, session_length, interval, use_rinex_flag, &
                       check_pc, pclimit, cutoff_elevation, use_brdeph, neph, &
                       ephem, stanam, x, y, z, t_first_in_rinex, t_last_in_rinex, v, sumepo)

! remove short piece and flag gap
  do iprn = 1, nsat
    again = .true.
    do while (again)
! ********************************************************************** !
!               remove short piece and mark large gap                    !
! ********************************************************************** !
      call remove_short(keep_end, nepo, ti, flagall(1, iprn), length_short, length_gap, interval, flag_shrt, again)
    enddo
  enddo

! SD LC checking when necessary
  lfnsd = -1
  if (debug_sd) lfnsd = 0
  if (check_lc) then
! ********************************************************************** !
!                  check single-difference LC                            !
! ********************************************************************** !
    call check_sd(lfnsd, neph, ephem, 2, 3, 20, 8, x, y, z, interval, lclimit)
  endif
!
! widelane checking etc.
  do iprn = 1, nsat
!  write(*,'(/a,i4)') '## tedit for satellite, ',iprn
    nused = 0
    do iepo = 1, nepo
      obs(iepo, iprn, 4) = 0.d0
      obs(iepo, iprn, 5) = 0.d0
      if (istrue(flagall(iepo, iprn), 'ok')) then
        pg(iepo) = obs(iepo, iprn, 1)       ! geometry-free
        nused = nused + 1
      endif
    enddo
    if (nused .ne. 0 .and. turbo_edit) then
!! ****************************************************************** !!
!!                       Used MW obs to check epochs
!!           Nov. 1, 2007. "limit" should not be too small
!! ****************************************************************** !!
      call edit_widelane(nepo, ti, obs(1, iprn, 2), flagall(1, iprn), 3.0d0)
      if (check_lc) then
!! ****************************************************************** !!
!!                Check ionosphere observations: LG
!!           Nov. 1, 2007. "limit" should not be too small
!! ****************************************************************** !!
        call check_ionosphere(nepo, flagall(1, iprn), ti, pg, respg, sigpg, interval, lglimit, lgrmslimit)
        call lc_help(nepo, flagall(1, iprn))
      endif
    endif
    if (nused .ne. 0) then
      again = .true.
      do while (again)
        call remove_short(keep_end, nepo, ti, flagall(1, iprn), length_short, length_gap, interval, flag_shrt, again)
      enddo
      if (debug_tb) then
        do iepo = 1, nepo
          if (.not. istrue(flagall(iepo, iprn), 'nodata')) then
            string = 'ok'
            if (istrue(flagall(iepo, iprn), 'no4')) then
              string = 'No 4'
            else if (istrue(flagall(iepo, iprn), 'lowele')) then
              string = 'Low Elevation'
            else if (istrue(flagall(iepo, iprn), 'shrt')) then
              string = 'Short Piece'
            else if (istrue(flagall(iepo, iprn), 'lwbad')) then
              string = 'Bad Widelane'
            else if (istrue(flagall(iepo, iprn), 'lgbad')) then
              string = 'Bad Ionosphere'
            else if (istrue(flagall(iepo, iprn), 'lccheck')) then
              string = 'Can not check LC'
            else if (istrue(flagall(iepo, iprn), 'pcbad')) then
              string = 'Bad Range   '
            else if (istrue(flagall(iepo, iprn), 'pc1ms')) then
              string = 'Bad Range 1ms'
            else if (istrue(flagall(iepo, iprn), 'amb')) then
              string = 'AMB'
            endif
          endif
        enddo
      endif
    endif
  enddo

! write rhd file for further processing
  if (flnrhd(1:1) .ne. ' ') then
    write (*, *) 'Write RHD...'
    call write_diag_rpt(flnrhd, nepo, nsat, jd0, ts, flagall, interval, sumepo)
  endif

end
