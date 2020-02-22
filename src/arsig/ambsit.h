!
!! Ambiguity Information in One Station
type ambsit
  character*4 name
  integer*4 now,nsd
!! satellite index in global table
  integer*4 isat(MAXOW_ST)
!! epoch index
  integer*4 iepc(2,MAXOW_ST)
!! pointer to inversed normal matrix
  integer*4 ipt(MAXOW_ST)
!! useful value
  real*8 xamb(MAXOW_ST),xrwl(MAXOW_ST),xrms(MAXOW_ST),xswl(MAXOW_ST)
end type
