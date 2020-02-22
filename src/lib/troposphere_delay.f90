!
!! troposphere_delay.f90
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
!! purpose  : troposhperic delay correction. Here we directly call several
!!            ztd model and mapping function code from GAMIT
!! parameter: jd,sod -- julian day (GPST)
!!            SITE%-- station information struct
!!            zdd,zwd --- dry and wet zenith delay, in meter
!
subroutine troposphere_delay(jd, sod, SITE, zdd, zwd)
  implicit none
  include '../header/const.h'
  include '../header/station.h'

  integer*4 jd
  real*8 sod, zdd, zwd
  type(station) SITE
!
!! local
  integer*4 i, jdutc
  real*8 al0, al1, tlen, sodutc
!
!! functions called
  real*8 saaszd, timdif, taiutc
!
!! transform time to UTC
  call timinc(jd, sod, 19.d0 - taiutc(jd), jdutc, sodutc)
!
!! a priori meteorological delays: hodrostatic & wet
  if (SITE%map(1:3) .eq. 'VM1') then    ! VMF1 troposphere model
    do i = 1, SITE%nvm - 1
      if (timdif(SITE%jdv(i), SITE%sodv(i), jdutc, sodutc) .le. 0.d0 .and. &
          timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), jdutc, sodutc) .ge. 0.d0) then
!! linear interpolation
        tlen = timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), SITE%jdv(i), SITE%sodv(i))
        al0 = timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), jdutc, sodutc)/tlen
        al1 = timdif(jdutc, sodutc, SITE%jdv(i), SITE%sodv(i))/tlen
        zdd = al0*SITE%vm1(3, i) + al1*SITE%vm1(3, i + 1)
        zwd = al0*SITE%vm1(4, i) + al1*SITE%vm1(4, i + 1)
        exit
      endif
    enddo
  else                               ! Saastamoinen model
!
!! hystrostatic zenith delay
    zdd = saaszd(SITE%p0, SITE%t0, 0.d0, 'R', SITE%geod(1), SITE%geod(3) - SITE%undu*1.d-3)
!
!! wet zenith delay
    zwd = saaszd(0.d0, SITE%t0, SITE%hr0, 'R', SITE%geod(1), SITE%geod(3) - SITE%undu*1.d-3)
  endif
!
!! check availability
  if (zdd .eq. 0.d0 .or. zwd .eq. 0.d0) then
    write (*, '(a,i5,f9.2)') '***ERROR(troposphere_delay): fail to compute a priori delays ', jd, sod
    call exit(1)
  endif

  return
end
