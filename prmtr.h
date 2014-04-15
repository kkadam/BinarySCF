!Needed for subroutines setbdy, bdygen, pot3, 
      parameter(jmax=254, jmax1=jmax+1, jmax2=jmax+2)
      parameter(kmax=254, kmax1=kmax+1, kmax2=kmax+2)
      parameter(lmax=512)

!            From subroutine bdygen...
      parameter(jkm1 = 1*jmax+kmax-1)
!            From subroutine pot3...
      parameter(lmax2 = 2*lmax, l2p1 = lmax/2 + 1,   &
               mwfw=3000, kkk = kmax/2 + 1,          &
               ijj = jmax-1, ikk = kmax-1)
      parameter(ifft=2*lmax)

