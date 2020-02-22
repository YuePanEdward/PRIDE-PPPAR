!
!! antenna atx
!
type antatx
  character*20 antnam
  character*20  antnum 
  integer*4    nfreq
  character*3  frq(2) 
  real*8 zen1,zen2,dzen,dazi
  real*8 neu(3,2)
  real*8 pcv(30,0:80,2)
end type
