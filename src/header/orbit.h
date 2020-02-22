!
!! Binary orbit header
!
type orbhdr
!
!! satelilte 
  character*10 :: sattyp = ''   ! satellite type  
  integer*4 nprn                ! # of satellites
  integer*4 prn(MAXSAT)         ! satellite number 
!
!! system tag
  character*80 :: iers = ''     ! IERS Conventions
!
!! time tag
  integer*4 jd0,jd1             !
  real*8 sod0,sod1              ! start & end time
!
!! sp3 interval (second)
  real*8 dintv
end type
