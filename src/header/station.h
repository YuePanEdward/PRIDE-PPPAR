type station
  character*4 name
  character*2 skd
  integer*4 iunit,lfnjmp,imet,ikin,iion,itrp,iptatx
  real*8 x(6),dx0(3),geod(3)
!
!! related file names
  character*256 obsfil,metfil
  character*20 kinfil,ionfil,trpfil
!
!! receiver & antenna
  real*8 enu0(3),enu(3,2),rot_l2f(3,3)
  character*20 rcvtyp,anttyp
!
!! clock correction
  real*8 rclock,dclk0
!
!! observation info
  real*8 sigr,sigp,cutoff
!
!! meteorology info
  integer*4 nmet
  character*2 mete(10)
!
!! meteorology
  character*3 map
  integer*4 nvm,jdv(13)
  real*8 ztdcor,dztd0,qztd,dhtg0,qhtg,p0,t0,hr0,undu,sodv(13),vm1(4,13) ! ah aw zhd zwd
!
!! phase wind-up
  logical*1 first(maxsat)
  real*8 prephi(maxsat)
!
!! ocean load
  real*8 olc(11,6),rlat,rlon
end type 
