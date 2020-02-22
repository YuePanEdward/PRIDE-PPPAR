!
!! Constant parameters
integer*4,parameter :: oscr=6

real*8, parameter   :: vlight=299792458.d0,pi=3.141592654d0,   &
                       freq1=1.57542d9,freq2=freq1*60.d0/77.d0

real*8, parameter   :: gm=398600.4418d0,gms=1.32712442076d11,gmm=4.90280070054d3   ! km3/s-2
real*8, parameter   :: erad=6378.1366d0   ! km
real*8, parameter   :: omega=7.292115d-5  ! rad/s

real*8, parameter   :: gpstai=19.d0,taitdt=32.184d0,gpstdt=gpstai+taitdt,mjd2jd=2400000.5d0

integer*4,parameter :: maxsit=900,maxsat=64,maxepo=432000
!
!! maxtyp : observation types
integer*4,parameter :: maxeph=maxsat*48,maxtyp=32
!
integer*4,parameter :: maxpar_sta=10,maxow_st=maxsat*10,maxpar=maxsit*(maxpar_sta+6)+maxsit*maxow_st
integer*4,parameter :: maxsd_sit=maxsat*400
integer*4,parameter :: maxpp=maxsat*(maxsat-1)/2

real*8,parameter    :: maxwnd=0.01d0
