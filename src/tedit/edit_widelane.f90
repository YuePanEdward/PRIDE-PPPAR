!
!! edit_widelane.f90
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
!! purpose:
!!
!! input:   limit =3.0
!!          nw = obs(1,iprn,2)
!!
!!
!
subroutine edit_widelane(nepo, ti, nw, flagall, limit)
  implicit none
  integer*4 nepo, flagall(1:*)
  real*8 ti(1:*), nw(1:*), limit

! function called
  logical*1 ok, istrue
  integer*4 set_flag

! local
  integer*4 i, ilast, kobs, flg(nepo)
  real*8 mb, sigb, sigmb, dif, temp(nepo)

! simply check for jump
  ilast = 0
  do i = 1, nepo
    flg(i) = 2
    if (istrue(flagall(i), 'ok')) then
      flg(i) = 1
      if (.not. istrue(flagall(i), 'lwjump') .and. &
          .not. istrue(flagall(i), 'gap') .and. &
          .not. istrue(flagall(i), 'bigsd')) flg(i) = 0
!   if(ilast.ne.0) then
!     dif=nw(i)-nw(ilast)
!     if(dabs(dif).gt.limit*1.5) then
!       flg(i)=1
!       write(*,*) ' WM flagged by checking differenec ',i,dif
!     endif
!   endif
!   ilast=i
    endif
  enddo

! sequential mean of NW and flag possible cycle slips by testing
!          dabs(nwi-meani) > limit cycles
! if a data point is flagged, mean calculation is re-initialized
! by finding the next two points with a dif. less than limit.

  kobs = 0
  ilast = 0
  i = 1
  do while (i .le. nepo)
    if (flg(i) .lt. 2) then                      ! point with data
      if (flg(i) .eq. 1) kobs = 0
      if (kobs .eq. 0) then                      ! first point
        mb = nw(i)
        sigb = 0.d0
        ilast = i
        kobs = 1
        flg(i) = 1
      else if (kobs .eq. 1) then                 ! second point
        if (dabs(nw(i) - mb) .le. limit) then      ! two points are ok
          kobs = kobs + 1
          dif = nw(i) - mb
          mb = mb + dif/kobs
          sigb = sigb + (dif*dif - sigb)/kobs
        else                                  ! two points not ok
          dif = 0.d0                   !* only for output
!        write(*,'(a,i6,f13.2)') ' Widelane bad point,',i,nw(ilast)-nw(i)
!       flg(ilast)= 2                       ! flag the first one as bad. as short piece not bad
          flg(i) = 1                            ! piece start
          mb = nw(i)
        endif
      else                                    ! third or later point
        dif = nw(i) - mb
        if (dabs(dif) .le. limit .or. (dabs(dif) .le. 1.4*limit .and. &
                                       dabs(nw(i) - nw(ilast)) .le. 0.5)) then   ! ok. accumulate mean
          kobs = kobs + 1
          mb = mb + (nw(i) - mb)/kobs
          sigb = sigb + (dif*dif - sigb)/kobs
        else                                  !possible cycle slip
!        write(*,'(a,i6,f13.2)') ' Widelane cycle slip,',i,dif
          kobs = kobs + 1
          mb = nw(i)
          sigb = 0.d0
          kobs = 1
          flg(i) = 1
        endif
      endif
      ilast = i
      temp(i) = dif
    endif
    i = i + 1
  enddo

  do i = 1, nepo
    if (flg(i) .eq. 1 .and. .not. istrue(flagall(i), 'lwjump')) &
      flagall(i) = set_flag(flagall(i), 'lwjump')
  enddo
  return
end
