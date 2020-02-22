!
!! fix_ambiguity.f90
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
!! purpose  : fix ambiguity in one station
!! parameter:
!!    input : FCB  -- fractional cycle biases
!!    output: AS   -- ambiguity station struct
!!            ASD  -- single-difference struct
!
subroutine fix_ambiguity(FCB, AS, ASD)
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include 'abfcb.h'
  include 'ambsit.h'

  type(abfcb) FCB
  type(ambsit) AS
  type(difamb) ASD(1:*)
!
!! local
  logical*1 ldep
  integer*4 i, j, k, l, jd, ik, il, isat, jsat, isd, nwnfx, nwlfx, ntot, ndef
  real*8 prob, alpha, val, sodb, sode, fwl, swl
  character*5 fixed
  type(difamb) SDX, SD(MAXSD_SIT)

  nwnfx = 0
  nwlfx = 0
  do isd = 1, AS%nsd
    fixed = ' '
    ASD(isd)%dec = 0.d0
!
!! pointer to satellites
    isat = AS%isat(ASD(isd)%ipt(1))
    jsat = AS%isat(ASD(isd)%ipt(2))
!
!! try to fix
    call bdeci(ASD(isd)%rwl, ASD(isd)%swl, 1, FCB%wl_maxdev, FCB%wl_maxsig, prob, alpha)
    if (prob .gt. 0.d0 .and. alpha .gt. FCB%wl_alpha) then
      ASD(isd)%id = 1               ! wide-lane fixed
      fixed = 'WL'
      nwlfx = nwlfx + 1
!
!! compute narrow-lane ambiguities
      ASD(isd)%rnl = 137.d0/77.d0*(AS%xamb(ASD(isd)%ipt(1)) - AS%xamb(ASD(isd)%ipt(2))) - 60.d0/17.d0*nint(ASD(isd)%rwl)
      ASD(isd)%snl = 0.05d0 ! presumed statistics, see Ge et al. (2005)
!
!! try to fix
      call bdeci(ASD(isd)%rnl, ASD(isd)%snl, 1, FCB%nl_maxdev, FCB%nl_maxsig, prob, alpha)
      ASD(isd)%dec = alpha
      if (prob .gt. 0.d0 .and. alpha .gt. FCB%nl_alpha) then
        ASD(isd)%id = 0                     ! narrow-lane fixed already
        fixed = 'WL_NL'
        nwnfx = nwnfx + 1
      endif
!
!! output
      call timinc(FCB%jd0, FCB%sod0, (ASD(isd)%iepc(1) - 1)*FCB%dintv, jd, sodb)
      call timinc(FCB%jd0, FCB%sod0, (ASD(isd)%iepc(2) - 1)*FCB%dintv, jd, sode)
      write (*, '(i5,1x,a4,2i3,2(f14.3,f8.3),2f10.1,1x,a5)') isd, AS%name, FCB%prn(isat), FCB%prn(jsat), &
        ASD(isd)%rwl, ASD(isd)%swl, ASD(isd)%rnl, ASD(isd)%snl, sodb, sode, fixed
    endif
  enddo
  if (nwlfx .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'Wide/Narrow-lane FR(whole): ', AS%name, nwnfx, nwlfx, AS%nsd, &
      nwlfx*1.d2/AS%nsd, '%', nwnfx*1.d2/nwlfx, '%'
  endif
!
!! sort ASD
  do k = 1, 2
    do i = 1, AS%nsd - 1
      if (k .eq. 1 .and. ASD(i)%id .ne. 0) cycle
      if (k .eq. 2 .and. ASD(i)%id .ne. 1) cycle
      do j = i + 1, AS%nsd
        if (k .eq. 1 .and. ASD(j)%id .ne. 0) cycle
        if (k .eq. 2 .and. ASD(j)%id .ne. 1) cycle
        if (ASD(i)%dec .lt. ASD(j)%dec) then
          SDX = ASD(i)
          ASD(i) = ASD(j)
          ASD(j) = SDX
        endif
      enddo
    enddo
  enddo
!
!! select independent SD for this station
  ndef = 0
  do k = 1, 3
    do i = 1, AS%nsd
      if (k .eq. 1 .and. ASD(i)%id .ne. 0) cycle
      if (k .eq. 2 .and. ASD(i)%id .ne. 1) cycle
      if (k .eq. 3 .and. ASD(i)%id .ne. 2) cycle
      call check_amb_depend(AS%now, ndef, 2, ASD(i)%ipt, ldep)
      if (ldep) cycle
      if (k .le. 2) SD(ndef) = ASD(i)
    enddo
    if (k .eq. 1) nwnfx = ndef
    if (k .eq. 2) nwlfx = ndef
    if (k .eq. 3) ntot = ndef
  enddo
  call check_amb_depend(0, 0, 2, (/1, 1/), ldep)
  write (*, '(a,3i8,3(f6.1,a1,2x))') 'Station FR(indep): ', nwnfx, nwlfx, ntot, nwlfx*1.d2/ntot, '%', &
    nwnfx*1.d2/nwlfx, '%', nwnfx*1.d2/ntot, '%'
!! only save widelane fixed ones
  AS%nsd = nwlfx
  do i = 1, AS%nsd
    ASD(i) = SD(i)
  enddo
  return
end
