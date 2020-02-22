!!
!! purpose : GPS broadcast format
type brdeph
  integer*4 svn
  real*8 a0,a1,a2,aode,crs,dn,m0,cuc,e,cus,roota
  real*8 toe,cic,omega0,cis,i0,crc,omega,omegadot
  real*8 idot,resvd0,week,resvd1,accu,hlth,tgd,aodc
end type
