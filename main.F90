       program main
       implicit none       
       include "init.h" 

       print*,ra,rb,rc,type,rhom1,rhom2

       open(unit=11,file='output',form='formatted',status='unknown')
       open(unit=12,file='summary',form='formatted',status='unknown')
       
         call binary_scf(iam,ra,rb,rc,type,model_num,rhom1,rhom2)
       
       close(11)
       close(12)

       end program main
