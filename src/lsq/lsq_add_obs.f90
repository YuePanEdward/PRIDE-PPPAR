!
!! lsq_add_obs.f90
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
!! purpose   : add observation equations to estimator LSQ
!! paraemters:
!!             lfncid,lfnobs -- tmp files for recovery
!!             jd,sod        -- float julian day
!!             isit          -- index of station
!!             LCF           -- LSQ control parameters
!!             OB            -- Observation struct
!!             PM,NM         -- normal matrix & PAR table
!
subroutine lsq_add_obs(lfncid, lfnobs, jd, sod, isit, LCF, OB, PM, NM, bias)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'
  include 'lsq.h'
  include 'lsqcfg.h'

  integer*4 isit, lfnobs, lfncid, jd
  real*8 sod
  type(lsqcfg) LCF
  type(rnxobr) OB
  type(prmt) PM(1:*)
  type(norm) NM
!
!! local
  logical*1 lfirst, lexist
  integer*4 ind, isat, ipar, i, k, ir, ic, nelem, ipt(0:MAXPAR)
  real*8 g1, g, rwl, lamdw, phase, range, wphs, wrng, wele, amat(MAXPAR), bias(MAXSAT, 4), wl

  data lfirst/.true./
  save lfirst, g, g1, lamdw

  if (lfirst) then
    lfirst = .false.
    g = 7.7d0/6.0d0
    g1 = g/(g*g - 1.d0)
    lamdw = VLIGHT/FREQ1/(1.d0 - 1.d0/g)
  endif
!
!! ambiguity file exist or not: whether initialize ambiguity estimates
  inquire (file=LCF%flnamb, exist=lexist)
  do isat = 1, LCF%nprn
    if (OB%omc(isat, 1) .eq. 0.d0 .or. OB%omc(isat, 3) .eq. 0.d0) cycle
    range = 0.d0; wrng = 0.d0
    phase = 0.d0; wphs = 0.d0
!
!! right hand side
    NM%nobs = NM%nobs + 1
    wrng = 1.d0/(OB%var(isat, 3) + OB%var(isat, 4))
    range = OB%omc(isat, 3) - g1*(OB%omc(isat, 4) - OB%omc(isat, 3)/g)

    NM%nobs = NM%nobs + 1
    wphs = 1.d0/(OB%var(isat, 1) + OB%var(isat, 2))
    phase = OB%omc(isat, 1) - g1*(OB%omc(isat, 2) - OB%omc(isat, 1)/g)

!
!! observation equations
    nelem = 0
    ipt(0) = 0
    do ipar = 1, OB%npar
      if (OB%ltog(ipar, isat) .eq. 0) cycle
!
!! extract non-zero ones
      nelem = nelem + 1
      ipt(nelem) = OB%ltog(ipar, isat)
      amat(nelem) = OB%amat(ipar, isat)
      ind = OB%ltog(ipar, isat)
      PM(ind)%iobs = PM(ind)%iobs + 1
!
!! initial value for the new ambiguity parameters
      if (OB%pname(ipar) (1:4) .eq. 'AMBC') then
        rwl = OB%obs(isat, 1) - OB%obs(isat, 2) - (g*OB%obs(isat, 3) + OB%obs(isat, 4))/(1.d0 + g)/lamdw
        wele = 1.d0
        !if(OB%elev(isat)*PI/180.d0.le.30.d0) wele=2.d0*dsin(OB%elev(isat))
        if (OB%elev(isat)*180.d0/PI .le. 30.d0) wele = 2.d0*dsin(OB%elev(isat))  !mod by zwx 20141031
        ipt(0) = nelem
        if (OB%flag(isat) .eq. 1 .and. .not. lexist) then
          PM(ind)%xini = phase - range
          PM(ind)%zw = rwl
        endif
        phase = phase - amat(nelem)*PM(ind)%xini
!
!! processing of wide-lane observation
        rwl = rwl - PM(ind)%zw
        PM(ind)%xrwl = PM(ind)%xrwl + rwl*wele
        PM(ind)%rw = PM(ind)%rw + wele
        PM(ind)%xswl = PM(ind)%xswl + wele*rwl**2
        PM(ind)%mele = PM(ind)%mele + OB%elev(isat)
      endif
!
!! next parameters
    enddo
!
!! save to recover residuals
    if (lfnobs .ne. 0) then
      write (lfncid) 'ob'
      write (lfnobs) jd, sod, isit, isat, nelem, ipt(0), (ipt(i), amat(i), i=1, nelem), &
        phase, range, wphs, wrng, OB%flag(isat), OB%elev(isat), OB%azim(isat), OB%dmap(isat), OB%wmap(isat)
    endif
!
!! transform to index of normal matrix
    do i = 1, nelem
      ipt(i) = PM(ipt(i))%ipt
    enddo
!
!! upper normal matrix is preferred
    do i = 1, nelem
      do k = i, nelem
        ir = min(ipt(i), ipt(k))
        ic = max(ipt(i), ipt(k))
        NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wphs
        if (i .ne. ipt(0) .and. k .ne. ipt(0)) NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wrng
      enddo
      ir = ipt(i)
      NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wphs*phase
      if (i .ne. ipt(0)) NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wrng*range
    enddo
    NM%ltpl = NM%ltpl + wphs*phase**2
    NM%ltpl = NM%ltpl + wrng*range**2
!
!! next satellite
  enddo

  return
end
