!
!! INVersed Normal Matrix & Parameter Information
type pest
  character*15 pname
  integer*4 pcode(2)
  real*8 xest,ptime(2)
end type

type invm
!! frdm : number of freedom
!! ntot : total number of parameters
!! nxyz : number of non-ambiguity parameters
!! indp : number of independent differential ambiguities
!! ndam : number of rest un-fixed differential ambiguities
!! ncad : number of to-be-fixed differential ambiguities
!! nfix : number of fixed differential ambiguities
  integer*4 frdm,nsiz,ntot,nxyz,indp,ndam,ncad,nfix
  real*8 vtpv
  integer*4,pointer :: idq(:)
  real*8,pointer :: invx(:,:)
end type
