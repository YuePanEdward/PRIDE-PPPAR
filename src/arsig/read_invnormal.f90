!
!! read_invnormal.f90
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
!! purpose  : read inversed normal matrix (Q matrix)
!! parameter:
!!    input : FCB -- fractional part of initial phases
!!    output: QN  -- inversed normal matrix
!!            AS  -- ambiguity station struct
!
subroutine read_invnormal(FCB, PM, QN, AS)
  implicit none
  include '../header/const.h'
  include '../header/invnor.h'
  include 'abfcb.h'
  include 'ambsit.h'

  type(abfcb) FCB
  type(ambsit) AS
  type(pest) PM(1:*)
  type(invm) QN
!
!! local
  integer*4 lfn, i, j, isit, isat, nprn, prn(MAXSAT), ierr
  real*8 xrwl, xswl, elev
  character*20 anttyp
!
!! function called
  integer*4 get_valid_unit, pointer_int

  lfn = get_valid_unit(10)
  open (lfn, file=FCB%flnneq, form='unformatted', status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(read_invnormal): open file ', trim(FCB%flnneq)
    call exit(1)
  endif
  read (lfn) FCB%nsit, AS%name
  AS%now = 0
  read (lfn) nprn, (prn(i), i=1, nprn)
  read (lfn) QN%ntot, QN%vtpv, QN%frdm
  if (QN%ntot .gt. 3 + MAXOW_ST) then
    write (*, '(a)') '***ERROR(read_invnormal): too many parameters '
    call exit(1)
  endif
  allocate (QN%invx(1:QN%ntot, 1:QN%ntot), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(read_invnormal): memory allocation invx '
    call exit(1)
  endif
  allocate (QN%idq(QN%ntot), stat=ierr)
  do i = 1, QN%ntot
    QN%idq(i) = (i - 1)*QN%ntot
  enddo
!
!! read parameters
  QN%nxyz = 0
  do i = 1, QN%ntot
    read (lfn) PM(i)%pname, isit, isat, xrwl, PM(i)%xest, PM(i)%ptime(1:2), xswl, elev
    PM(i)%pcode(1) = isit
    PM(i)%pcode(2) = isat
    if (PM(i)%pname(1:4) .ne. 'AMBC') then
      QN%nxyz = QN%nxyz + 1
      cycle
    endif
    PM(i)%pcode(2) = pointer_int(FCB%nprn, FCB%prn, prn(isat))
    if (PM(i)%pcode(2) .eq. 0) then
      write (*, '(a,i2)') '***ERROR(read_invnormal): satellite not exist ', prn(isat)
      call exit(1)
    endif
    if (elev .le. FCB%cutoff .or. xswl*3.d0 .gt. 0.2d0 .or. (PM(i)%ptime(2) - PM(i)%ptime(1))*86400.d0 .lt. FCB%minsec_common) cycle
    AS%now = AS%now + 1
    if (AS%now .gt. MAXOW_ST) then
      write (*, '(a,a4)') '***ERROR(read_invnormal): too many one-way ambiguities ', AS%name
      call exit(1)
    endif
    AS%isat(AS%now) = PM(i)%pcode(2)
    AS%xamb(AS%now) = PM(i)%xest
    AS%xrwl(AS%now) = xrwl
    AS%xswl(AS%now) = xswl
    AS%iepc(1, AS%now) = nint(((PM(i)%ptime(1) - FCB%jd0)*86400.d0 - FCB%sod0)/FCB%dintv) + 1
    AS%iepc(2, AS%now) = nint(((PM(i)%ptime(2) - FCB%jd0)*86400.d0 - FCB%sod0)/FCB%dintv) + 1
    AS%ipt(AS%now) = i       ! pointer to Q matrix
  enddo
!
!! read inversed normal matrix
  QN%invx(1:QN%ntot, 1:QN%ntot) = 0.d0
  read (lfn) ((QN%invx(i, j), i=j, QN%ntot), j=1, QN%ntot)   ! lower part
  close (lfn)

  return
end
