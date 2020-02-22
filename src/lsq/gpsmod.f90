!
!! gpsmod.f90
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
!! purpose   : modelling of GPS observations, observation equaitons are storaged in
!!             OB%amat and OB%omc etc.
!! parameters:
!!            jd,sod -- epoch time of the observations (GPS time)
!!            LCF    -- LSQ control parameters
!!            SITE   -- SITE information
!!            OB     -- observations
!!            SAT    -- satellite information
!! tester: Y. Pan, X. Chen, J. Zhou, S. Mao
!
subroutine gpsmod(jd, sod, LCF, SITE, OB, SAT)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include '../header/satellite.h'
  include 'lsqcfg.h'

  integer*4 jd
  real*8 sod
  type(lsqcfg) LCF
  type(station) SITE
  type(rnxobr) OB
  type(satellite) SAT(MAXSAT)
!
!! local
  integer*4 i, j, k, isat, jdc, jd_send, jd_recv, ierr
  real*8 sodc, sod_send, sod_recv, drecclk, dsatclk, &
    dsatclk0, dsatclkrat, drate, delay(2), ddelay, r1(3), &
    r2(3), r1leng, r2leng, phase1, phase2, range1, range2, &
    fjd, xsun(6), xlun(6), dx(3)
  real*8 gast, xpole, ypole, trpdel, trpart(3), dphwp, enu(6), pcv(2), sitrad, satrad, reldel, nadir
  real*8 xant_f(6, 2), xant_j(6, 2), xsat_j(6), dloudx(3), dump(3), &
    rot_f2j(3, 3), rot_rat(3, 3), xrot(3, 3)
!
!! function called
  real*8 dot
!
!! initial
  do i = 1, 6
    enu(i) = 0.d0
    xsun(i) = 0.d0
    xlun(i) = 0.d0
    xsat_j(i) = 0.d0
    xsat_j(i) = 0.d0
    do j = 1, 2
      xant_j(i, j) = 0.d0
      xant_f(i, j) = 0.d0
    enddo
  enddo
  do i = 1, 3
    trpart(i) = 0.d0
    dloudx(i) = 0.d0
    dump(i) = 0.d0
    dx(i) = 0.d0
    r1(i) = 0.d0
    r2(i) = 0.d0
    do j = 1, 3
      rot_f2j(i, j) = 0.d0
      rot_rat(i, j) = 0.d0
      xrot(i, j) = 0.d0
    enddo
  enddo
  delay(1) = 0.d0
  delay(2) = 0.d0
  pcv(1) = 0.d0
  pcv(2) = 0.d0
!
!! atmosphere delays
  OB%zdd = 0.d0; OB%zwd = 0.d0
  call troposphere_delay(jd, sod, SITE, OB%zdd, OB%zwd)
!
!! atmosphere gradients
  OB%nhtg = 0.d0; OB%ehtg = 0.d0
!
!! convert antenna offset from east-north-up to x-y-z in earth fixed system
  do j = 1, 2
    dump(1:3) = SITE%enu0(1:3) + SITE%enu(1:3, j)
    call matmpy(SITE%rot_l2f, dump, dump, 3, 3, 1)
    do i = 1, 3
      xant_f(i, j) = SITE%x(i) + dump(i)*1.d-3
      xant_f(i + 3, j) = 0.d0
    enddo
  enddo
!
!! receiver clock correction
  drecclk = 0.d0
  drecclk = SITE%rclock/VLIGHT
!
!! if receiver clock is too bad we repeat model starting from here
100 continue
  call timinc(OB%jd, OB%tsec, -drecclk, jd_recv, sod_recv)
!
!! Compute transformation matrix from earth-fixed to inertial system
  call ef2int(LCF%flnerp, jd_recv, sod_recv, rot_f2j, rot_rat, gast, xpole, ypole)
  do i = 1, 2
    call matmpy(rot_f2j, xant_f(1, i), xant_j(1, i), 3, 3, 1)
    call matmpy(rot_f2j, xant_f(4, i), dx, 3, 3, 1)
    call matmpy(rot_rat, xant_f(1, i), xant_j(4, i), 3, 3, 1)
    xant_j(4:6, i) = xant_j(4:6, i) + dx(1:3)
  enddo
!
!! read solar and lunar tables jd_send sod_send (AT)
  fjd = MJD2JD + jd_recv + (sod_recv - 0.075d0 + GPSTDT)/86400.d0
  call pleph(fjd, 11, 3, xsun)
  call pleph(fjd, 10, 3, xlun)
!
!! tide dpsplacement
  call tide_displace(jd_recv, sod_recv, xant_j(1, 1), xant_f(1, 1), xsun, &
                     xlun, rot_f2j, SITE%rot_l2f, SITE%geod(1), SITE%geod(2), gast, xpole, ypole, SITE%olc, dx)
  do j = 1, 2
    do i = 1, 3
      xant_j(i, j) = xant_j(i, j) + dx(i)
    enddo
  enddo
!
!! loop over all satellites
  do isat = 1, LCF%nprn
    if (OB%obs(isat, 1) .eq. 0.d0 .or. OB%obs(isat, 3) .eq. 0.d0) cycle
!
!! get satellite clock information
    call read_satclk(LCF%flnsck, OB%prn(isat), jd, sod, jdc, sodc, dsatclk0, dsatclkrat, ierr)
!
!### BLOCK 1
!### ITERATION OF SEND TIME and GET THE GEOMETRIC DISTANCE
    ddelay = 0.1d0
    do while (dabs(ddelay) .gt. 1.d-9)
      call timinc(OB%jd, OB%tsec, -drecclk - OB%delay(isat), jd_send, sod_send)
!
!! compute the satellite position and velocity coordinates
      call everett_interp_orbit(LCF%flnorb, .true., .true., jd_send, sod_send, LCF%prn(isat), xsat_j(1), xsat_j(4))
      if (all(xsat_j(1:6) .eq. 1.d15)) exit
!
!! the satellite unit vectors (should be corrected for yaw error)
      call rot_scfix2j2000(xsat_j, xsun, SAT(isat)%xscf, SAT(isat)%yscf, SAT(isat)%zscf)
      do i = 1, 3
        xsat_j(i) = xsat_j(i) + (SAT(isat)%xyz(1, 1)*SAT(isat)%xscf(i) + &
                                 SAT(isat)%xyz(2, 1)*SAT(isat)%yscf(i) + &
                                 SAT(isat)%xyz(3, 1)*SAT(isat)%zscf(i))*1.d-3
      enddo
!
!! geometric distance
      do i = 1, 3
        r1(i) = xsat_j(i) - xant_j(i, 1)
      enddo
      r1leng = dsqrt(dot(3, r1, r1))
      delay(1) = r1leng/VLIGHT*1.d3
      do i = 1, 3
        dloudx(i) = r1(i)/r1leng
      enddo
!
!! decide whether to iterate the delay calculation
      ddelay = delay(1) - OB%delay(isat)
      OB%delay(isat) = delay(1)
    enddo
    if (all(xsat_j(1:6) .eq. 1.d15)) cycle      ! satellite missing
    do i = 1, 3
      r2(i) = xsat_j(i) - xant_j(i, 2)
    enddo
    r2leng = dsqrt(dot(3, r2, r2))
    delay(2) = r2leng/VLIGHT*1.d3
!### END OF BLOCK 1
!###
!### BLOCK 2
!### COMPUTE CORRECTON
!
!! Compute the atmospheric corrections to the delay
!  compute the station-satellite elevation angle and azimuth from north
    call matmpy(r1, rot_f2j, r2, 1, 3, 3)
    call matmpy(r2, SITE%rot_l2f, dump, 1, 3, 3)
    OB%azim(isat) = datan2(dump(1), dump(2))
    OB%elev(isat) = datan(dump(3)/dsqrt(dump(1)**2 + dump(2)**2))
!
!! nadir angle for GPS satellite
    nadir = dot(3, xsat_j, r1)/dsqrt(dot(3, xsat_j, xsat_j))/dsqrt(dot(3, r1, r1))
    nadir = dacos(nadir)
!
!! get the nominal zenith delay, mapping function, and compute the path delay
    trpdel = 0.d0; trpart = 0.d0
    call troposphere_map(jd, sod, OB%elev(isat), SITE, OB%dmap(isat), OB%wmap(isat))
    trpdel = (OB%dmap(isat)*OB%zdd + OB%wmap(isat)*OB%zwd)/VLIGHT
    trpart(1) = OB%wmap(isat)
    if (LCF%htgmod(1:3) .ne. 'NON') then
      trpart(2) = OB%wmap(isat)/dtan(max(OB%elev(isat), 1.d-3))*dcos(OB%azim(isat))
      trpart(3) = OB%wmap(isat)/dtan(max(OB%elev(isat), 1.d-3))*dsin(OB%azim(isat))
      trpdel = trpdel + OB%nhtg*trpart(2)/VLIGHT + OB%ehtg*trpart(3)/VLIGHT
    endif
!
!! compute antenna / transmitter orientation dependent phase (cycle)
!! corrections for Right Circularly Polarized electro magnetic waves
    dphwp = 0.d0
    call phase_windup(SITE%first(isat), rot_f2j, SITE%rot_l2f, SAT(isat)%xscf, &
                      SAT(isat)%yscf, SAT(isat)%zscf, r1, SITE%prephi(isat), dphwp)
!
!! general relativistic time delay due to the Earth gravity (second)
    reldel = 2.d0*dot(3, xsat_j, xsat_j(4))/VLIGHT**2*1.d6
!
!! gravitional change effects on the satellite oscillators.
    sitrad = dsqrt(dot(3, xant_j(1, 1), xant_j(1, 1)))
    satrad = dsqrt(dot(3, xsat_j, xsat_j))
    reldel = reldel + 2.d0*GM/(VLIGHT/1.d3)**3*log((sitrad + satrad + r1leng)/(sitrad + satrad - r1leng))
!
!! pcv correction for satellite and receiver antenna
    call get_ant_pcv(SITE%iptatx, SAT(isat)%iptatx, PI/2.d0 - OB%elev(isat), OB%azim(isat), nadir, pcv)
    pcv(1) = pcv(1)/VLIGHT
    pcv(2) = pcv(2)/VLIGHT
!
!! delay
    delay(1) = delay(1) + trpdel + pcv(1) + reldel + dphwp/FREQ1
    delay(2) = delay(2) + trpdel + pcv(2) + reldel + dphwp/FREQ2
!
!! satellite clock correction, and save it for output clock estimates
    dsatclk = (jd - jdc)*86400.d0 + (sod - sodc)
    dsatclk = dsatclk0 + dsatclkrat*dsatclk
    SAT(isat)%sclock = dsatclk*VLIGHT
    phase1 = -FREQ1*(dsatclk - drecclk - delay(1))
    phase2 = -FREQ2*(dsatclk - drecclk - delay(2))
    range1 = phase1 - dphwp - FREQ1*pcv(1)
    range2 = phase2 - dphwp - FREQ2*pcv(2)
!
!! Finally, form observed minus calculated, in cycles.
    OB%omc(isat, 1) = OB%obs(isat, 1) - phase1
    OB%omc(isat, 2) = OB%obs(isat, 2) - phase2
    OB%omc(isat, 3) = FREQ1*OB%obs(isat, 3)/VLIGHT - range1
    OB%omc(isat, 4) = FREQ2*OB%obs(isat, 4)/VLIGHT - range2
!
!! weight of observations
    OB%var(isat, 1) = SITE%sigp**2
    OB%var(isat, 2) = SITE%sigp**2
    OB%var(isat, 3) = (SITE%sigr/0.19d0)**2
    OB%var(isat, 4) = (SITE%sigr/0.19d0)**2
    if (OB%elev(isat)/PI*180.d0 .le. 30.d0) then
      do i = 1, 4
        OB%var(isat, i) = OB%var(isat, i)/(2*dsin(OB%elev(isat)))**2
      enddo
    endif
!
!! Compute the delay rate and the predicted delay
    do i = 1, 3
      dump(i) = xsat_j(i + 3) - xant_j(i + 3, 1)
    enddo
    drate = dot(3, dump, r1)/(VLIGHT/1.d3*r1leng)
!  write(100,'(i2,f20.10)') isat,drate
    OB%delay(isat) = delay(1) + drate*LCF%dintv
!
!! observation equation
    xrot = rot_f2j
    call partial_gps(OB%npar, dloudx, drate, OB%pname, OB%ltog(1, isat), xrot, trpart, OB%amat(1, isat))
!
!! next satellite
  enddo
!
!! receiver clock offset
  k = count(OB%omc(1:LCF%nprn, 3) .ne. 0.d0)
  drecclk = sum(OB%omc(1:LCF%nprn, 3))
  if (k .ge. 1) drecclk = drecclk/(k*FREQ1)
  if (dabs(drecclk) .gt. 1.d-6 .or. (k .ge. 1 .and. SITE%rclock .eq. 0.d0)) then
    SITE%rclock = SITE%rclock + drecclk*VLIGHT
    drecclk = SITE%rclock/VLIGHT
    goto 100
  endif

  return
end
