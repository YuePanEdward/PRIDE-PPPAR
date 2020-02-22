!
!! rdrnxoi3.f90
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
!! purpose  : read one epoch data from a RINEX o-file
!!
!! parameter: lfn -- file unit
!!            jd0, sod0 --- julian day and second of day of the requested epoch
!!                          if they are zero, take the epoch the file pointer
!!                          points at.
!!            dwnd  --- window for time matching. If the obsersing time from the
!!                      rinex-file is close to the requested time within the
!!                      window, we take the data. Be careful with this parameter
!!                      when you are working sampling rate larger than 1Hz.
!!            nprn0,prn0 -- number of satellite and satellite PRNs are chosen
!!                          If nprn is zero, take all observation of the matched
!!                          epoch.
!!            HD -- rinex header structure
!!            OB -- observation structure
!!            ierr -- error code, end of file or read fil error
!!
!!
!! last mod.: 29-03-2003 by Xingyu Chen, rinex 3
!
subroutine rdrnxoi3(lfn, jd0, sod0, dwnd, nprn0, prn0, HD, OB, dcb, bias, ierr)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'

  integer*4 ierr, lfn, jd0, nprn0, prn0(1:*)
  real*8 sod0, dwnd
  type(rnxhdr) HD
  type(rnxobr) OB
!
!! local
  logical*1 prior_p1, prior_p2
  integer*4 ioerr, iy, im, id, ih, imi, nprn, prn(MAXSAT)
  integer*4 iflag, i, j, i0, nline, lli(MAXSAT), ssi(MAXSAT), ii, nobstype
  real*8 sec, ds, dt, c1, c2, p1, p2, obs(MAXTYP), dcb(MAXSAT, 2), bias(MAXSAT, 4)
  character*1 sysid(MAXSAT)
  character*80 line, cline, msg, name
!
!! RINEX-3 Signal Priority
  character*9 obs_prio
  integer*4 obs_prio_index(4)
  integer*4 prio_index
  character*1024 string
!
!! function used
  integer*4 modified_julday

  obs_prio = 'NMYXLSCWP'
  obs_prio_index = 0
  ierr = 0
  line = ' '
  prn = 0
10 continue  ! next record
  read (lfn, '(a)', end=200) line
  msg = ' '
!
!! start line
  if (line(1:1) .ne. '>') then   !!!! rinex 3.03
    goto 100
  endif
!
!! number of satellite
  read (line(33:35), '(i3)', iostat=ioerr) nprn
  if (ioerr .ne. 0) then
    msg = 'read satellite number error.'
  endif
  if (len_trim(msg) .ne. 0) goto 100
!
!! Check the RINEX 3 event flag
  read (line(30:32), '(i3)', iostat=ioerr) iflag
  if (ioerr .ne. 0) then
    msg = 'read event flag error.'
    goto 100
  else if (iflag .gt. 1) then
    msg = 'read internal antenna information error'
    do i = 1, nprn
      read (lfn, '(a80)', iostat=ioerr, end=200) line
      if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
        read (line, '(3f14.4)', err=100) HD%h, HD%e, HD%n
      endif
    enddo
    goto 10
  endif
!
!! initialization
  do i = 1, MAXSAT
    do j = 1, 6
      OB%obs(i, j) = 0.d0
    enddo
    OB%lli(i, 1:2) = 0
    OB%ssi(i, 1:2) = 0
  enddo
!
!! format of the time tag line
  msg = 'read time & svn error'
  read (line, '(1x,i5,4i3,f11.7,i3,i3)', err=100) iy, im, id, ih, imi, sec, iflag, nprn
  read (line(42:56), '(f15.12)', iostat=ioerr) dt
  if (ioerr .ne. 0) dt = 0.d0
  if (nprn .gt. MAXSAT + 100) then
    write (*, '(a,i3)') '***ERROR(rdrnxoi): nprn > maxsat ', nprn
    call exit()
  endif
!
!! check time
  if (im .le. 0 .or. im .gt. 12 .or. id .le. 0 .or. id .gt. 31 .or. ih .lt. 0 .or. ih .ge. 24 &
      .or. imi .lt. 0 .or. imi .gt. 60 .or. sec .lt. 0.d0 .or. sec .gt. 60.d0) then
    msg = 'epoch time incorrect'
    goto 100
  endif
!call yr2year(iy)
!
!! check on time tags. do not change the requested time if there is no data
  ds = 0.d0
  OB%jd = modified_julday(id, im, iy)
  OB%tsec = ih*3600.d0 + imi*60.d0 + sec
  if (jd0 .ne. 0) then
    ds = (jd0 - OB%jd)*86400.d0 + (sod0 - OB%tsec)
    if (ds .lt. -dwnd) then
      OB%jd = jd0
      OB%tsec = sod0
      backspace lfn
      OB%nprn = 0
      return
    else if (ds .gt. dwnd) then
      i = nprn
      do j = 1, i
        read (lfn, '(a)') line
      enddo
      line = ' '
      goto 10
    endif
  endif
!
!! in case of mult lines for PRN number, continous lines must be read here.
!ii=1
!do while(nprn.gt.12*ii)
!read(lfn,'(32x,12(a1,i2))',err=100) (sysid(i),prn(i),i=12*ii+1,min(12*ii+12,nprn))
!ii=ii+1
!end do
  ii = 0
!
!! read data. if more than 10 type 3 line should be merged to one
  do i = 1, nprn
!  if(sysid(i).eq.' ') sysid(i)='G'
    read (lfn, '(a)', err=100, end=200) string
    if (string(1:1) .eq. 'G') then
      ii = ii + 1
      !
      !! check if the sallite is requested
      read (string, '(a1,i2,20(f14.3,i1,i1))', err=100) sysid(ii), prn(ii), (obs(j), lli(j), ssi(j), j=1, min(HD%nobstyp, maxtyp))
      i0 = 0
      if (nprn0 .gt. 0) then
        do j = 1, nprn0
          if (prn0(j) .eq. prn(ii)) i0 = j
        enddo
      else
        i0 = ii
      endif
!
!! fill in the obs. structure
      if (i0 .ne. 0) then
        c1 = 0.d0; p1 = 0.d0; c2 = 0.d0; p2 = 0.d0
        obs_prio_index = 0
        prior_p1 = .false.
        prior_p2 = .false.
        do j = 1, HD%nobstyp
          prio_index = 0
          if (dabs(obs(j)) .lt. MAXWND) cycle
          if (HD%obstyp(j) (1:2) .eq. 'L1') then
            prio_index = index(obs_prio, HD%obstyp(j) (3:3))
            if (obs_prio_index(1) .lt. prio_index) then
              obs_prio_index(1) = prio_index
              OB%obs(i0, 1) = obs(j) - bias(prn0(i0), 1)*freq1/vlight
              OB%ssi(i0, 1) = ssi(j)
              OB%lli(i0, 1) = lli(j)
              if (dabs(obs(j)) .lt. 1.d-3) OB%ssi(i0, 1) = 1
            endif
          elseif (HD%obstyp(j) (1:2) .eq. 'L2') then
            prio_index = index(obs_prio, HD%obstyp(j) (3:3))
            if (obs_prio_index(2) .lt. prio_index) then
              obs_prio_index(2) = prio_index
              OB%obs(i0, 2) = obs(j) - bias(prn0(i0), 2)*freq2/vlight
              OB%ssi(i0, 2) = ssi(j)
              OB%lli(i0, 2) = lli(j)
              if (dabs(obs(j)) .lt. 1.d-3) OB%ssi(i0, 2) = 1
            endif
          elseif (HD%obstyp(j) (1:2) .eq. 'C1') then
            prio_index = index(obs_prio, HD%obstyp(j) (3:3))
            if (obs_prio_index(3) .lt. prio_index) then
              obs_prio_index(3) = prio_index
              OB%obs(i0, 3) = obs(j) - bias(prn0(i0), 3)
              if (obs_prio_index(3) .gt. index(obs_prio, 'C')) prior_p1 = .true.
            endif
          elseif (HD%obstyp(j) (1:2) .eq. 'C2') then
            prio_index = index(obs_prio, HD%obstyp(j) (3:3))
            if (obs_prio_index(4) .lt. prio_index) then
              obs_prio_index(4) = prio_index
              OB%obs(i0, 4) = obs(j) - bias(prn0(i0), 4)
              if (obs_prio_index(4) .gt. index(obs_prio, 'C')) prior_p1 = .true.
            endif
          endif
        enddo
!! if one of the phases is zero, or the data is removed before
        !write(317,'(i2,x,4f14.3)') prn(ii), OB%obs(i0,1:4)
        if (any(OB%obs(i0, 1:4) .eq. 0.d0)) then
          OB%obs(i0, 1:4) = 0.d0
        elseif (nprn0 .gt. 0) then  ! correct for P1C1 or P2C2 biases
          if (.not. prior_p1 .and. .not. HD%lc1p1) OB%obs(i0, 3) = OB%obs(i0, 3) + (dcb(prn0(i0), 1)/3.335641d0)
          if (.not. prior_p2 .and. .not. HD%lc2p2) OB%obs(i0, 4) = OB%obs(i0, 4) + (dcb(prn0(i0), 2)/3.335641d0)
        endif
      endif
    else
      cycle
    endif
  enddo
!
!! fill in PRNs
!! 8/5/2006.Geng J.H.: For GLONASS
  if (nprn0 .eq. 0) then
    OB%nprn = ii
    do i = 1, nprn
      if (sysid(i) .eq. 'G') then
        OB%prn(i) = prn(i)
      else
        OB%prn(i) = 0
      endif
    enddo
  else
    OB%nprn = nprn0
    do i = 1, nprn0
      OB%prn(i) = prn0(i)
    enddo
  endif
  OB%dtrcv = dt
!
!! normal ending
  return
!
!! error
100 continue
  ierr = 1
  inquire (unit=lfn, name=name)
  write (*, '(a/a/a)') '***ERROR(rdrnxoi): read file, '//trim(name), '   line :'//line(1:80), ' &
    msg  :'//trim(msg)
  call exit(1)
!
!! come here on end of file
200 continue
  ierr = 2
  inquire (unit=lfn, name=name)
  do i = 1, MAXSAT
    OB%obs(i, 1:4) = 0.d0
  enddo
  write (*, '(a/a)') '###WARNING(rdrnxoi): end of file, ', trim(name)
  return
end
