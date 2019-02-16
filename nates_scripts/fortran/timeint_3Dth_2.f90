!     Code to generate 3D time history file (binary format) from interpolating in time from an old file 
!     Input: th.old (binary); 
!            timeint.in: 1st line: total # of days, nvrt (1 for elev; use 2*nvrt for hvel);
!                        2nd line: nope, (nond(i),i=1,nope) (nond: # of nodes on each segment that requires *3D.th);
!                        3rd line: new dt (old dt is read in from the 1st line of th.old).
!     Output: th.new (binary)

!     Compile: ifort -Bstatic -assume byterecl -O3 -o timeint_3Dth_2 timeint_3Dth_2.f90
      program riverforcing
      parameter(mnope=10)
      parameter(mnond=10000)
      parameter(mnv=160)
      parameter(nbyte=4)
      dimension th1(mnope,mnond,mnv),th2(mnope,mnond,mnv)
      dimension iglobal(mnope,mnond),th(mnope,mnond,mnv)
      dimension nond(mnope)

      open(21,file='timeint.in',status='old')
      read(21,*)rndays,nvrt
      if(nvrt>mnv) then
        print*, 'nvrt>mnv!'
        stop
      endif
      read(21,*)nope,(nond(i),i=1,nope)
      nnodes=0 !total # of open bnd nodes
      do i=1,nope
        if(nope>mnope.or.nond(i)>mnond) then
          print*, 'nope >mnope'
          stop
        endif
        nnodes=nnodes+nond(i)
      enddo !i
      read(21,*)dt
      close(21)

      irecl=nbyte*(1+nnodes*nvrt)
      open(17,file='th.old',access='direct',recl=irecl,status='old')
      open(18,file='th.new',access='direct',recl=irecl,status='replace')
      read(17,rec=1)dt0,(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
      print*, 'Old time step=',dt0
      irec_out=0
      tt0=rndays*86400
      nt0=tt0/dt0+1.e-5
      nt1=tt0/dt+1.e-5
!     Read first step in case dt<dt0
      timeout=dt
      ncount=0
      if(dt<dt0) then
!        do i=1,nope
!          do j=1,nond(i)
!            read(17,*)iglobal(i,j),(th2(i,j,l),l=1,nvrt)
!          enddo !j
!        enddo !i

        do
          ncount=ncount+1
!          write(18,'(e20.12)')timeout
!          do i=1,nope
!            do j=1,nond(i)
          write(18,rec=irec_out+1)timeout,(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
          irec_out=irec_out+1
!          do l=1,nvrt
!            if(abs(th2(i,j,l))>1.e8) stop 'Output too large'
!          enddo !l
!            enddo !j
!          enddo !i
          timeout=timeout+dt
          if(timeout>=dt0) exit
        enddo
      endif !dt.lt.dt0


!     timeout>= dt0
!     Interpolate
      do it=1,nt0
        read(17,rec=it)time2,(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
        if(abs(time2-it*dt0)>1.e-4) then
          print*, 'Time stamp wrong:',time2,it*dt0
          stop
        endif
!        do i=1,nope
!          do j=1,nond(i)
!            read(17,*)iglobal(i,j),(th2(i,j,l),l=1,nvrt)
!          enddo !j
!        enddo !i

        if(it>1.and.timeout>=time1.and.timeout<=time2) then
          do
            rat=(timeout-time1)/dt0
            if(rat<0.or.rat>1) then
              print*, 'ratio out of bound:',rat,it
              stop
            endif
            ncount=ncount+1
!            write(18,'(e20.12)')timeout
            do i=1,nope
              do j=1,nond(i)
                do l=1,nvrt
                  th(i,j,l)=th2(i,j,l)*rat+th1(i,j,l)*(1-rat)
                  if(abs(th(i,j,l))>1.e8) stop 'Output too large'
                enddo !l
              enddo !j
            enddo !i
            write(18,rec=irec_out+1)timeout,(((th(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
            irec_out=irec_out+1
            timeout=timeout+dt
            if(timeout>time2) exit
          enddo
        endif !it.gt.1 etc

        th1=th2
        time1=time2
      enddo !it=1,nt0
      
      if(ncount/=nt1) then
        print*, 'Miscount:',ncount,nt1
        stop
      endif

      stop
      end 
