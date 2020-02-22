type lsqcfg
!
!! start & stop time
  integer*4 jd0,jd1
  real*8 dintv,sod0,sod1
!
!! number of satellite & PRNs
  integer*4 nprn,prn(MAXSAT)
!
!! number of stations
  integer*4 nsit
!
!! ztd model
  character*15 ztdmod,htgmod
!
!! common file
  character*20 flnorb,flnsck,flnrck,flnpos,flnneq,flnerp
  character*20 flnztd,flnamb,flnres,flncon,flnhtg,flnvmf,flnfcb
!
!! keep bias
  logical*1 lrmbias
end type
