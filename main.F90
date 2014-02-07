       program main
       implicit none       

       integer :: ra, rb, rc, type, iam, master, slave
       integer :: model_num, start_model, stop_model
       integer :: numnodes, counter
       integer :: ierror, done, I, N
       real :: rhom1, rhom2


! open and read initial conditions file
       open(unit=10,file='init',form='formatted',status='old')       
         read(10,*) ra, rb, rc, type, rhom1, rhom2
       write(*,*) 'Dispatched Model:'
       
       print*,N,ra,rb,rc,type,rhom1,rhom2

       open(unit=11,file='output',form='formatted',status='unknown')
       open(unit=12,file='summary',form='formatted',status='unknown')

         iam = 0
         model_num = 1000
         call binary_scf(iam,ra,rb,rc,type,model_num,rhom1,rhom2)

!         call flush_(11)
!         call flush_(12)

       close(11)
       close(12)
!       close(20)  wth??


       end program main
