!
!! rinex header information  
type rnxhdr 
  integer*4    ver,fact1,fact2,nobstyp
  integer*4    nprn,prn(maxsat)
  integer*4    t0(5),t1(5)
  real*8       intv,x,y,z,h,e,n,t0s,t1s
  character*1  sys
  character*3  obstyp(maxtyp)
  character*4  mark
  character*6  recnum,antnum
  character*12 file
  character*20 rectyp,anttyp
  logical*1 lc1p1,lc2p2
end type 
!
!! rinex observation records
type rnxobr 
  integer*4 jd
  integer*4 nprn,prn(maxsat)
  integer*4 lli(maxsat,2),ssi(maxsat,2)
  real*8    tsec,dtrcv
  real*8    obs(maxsat,6)
!
!! editing
  integer*4    lfnrhd,amb_epo,amb_tot,ava_obs,rem_obs,sum_obs,act_obs
  real*8       icyc(2,maxsat)
  character*20 rhdfil
!
!! epoch difference
  integer*4 iepd
  real*8    zmat(maxsat)
!
!! model
  integer*2 flag(maxsat)
  integer*4 npar,ltog(maxpar_sta,maxsat),lwlg(maxpar_sta,maxsat)
  real*8    ztdc,zage,zdd,zwd,nhtg,ehtg
  real*8    lifamb(maxsat,2),azim(maxsat),elev(maxsat),dmap(maxsat),wmap(maxsat)
  real*8    delay(maxsat),omc(maxsat,4),var(maxsat,4),amat(maxpar_sta,maxsat)
  real*8    sidereal(maxsat,2)
  character*15 pname(maxpar_sta)
!
!! wide-lane ambiguity
  integer*4         nien
  real*8            iage(maxsat),vard(maxsat),ida0(maxsat),ida1(maxsat)
  real*8,pointer :: iond(:,:)
!!  integer*4 abob(maxpp),ixab(maxpp)
!!  real*8 abwl(maxpp),iond(maxll,maxpp),iage(maxpp),vard(maxpp),ida0(maxpp),ida1(maxpp)
end type
