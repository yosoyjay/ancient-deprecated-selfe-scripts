!     Generate fg.bp for interpolation (e.g., *3D.th)
!     Input: hgrid.gr3; open boundary segment # (iob)
!     Output: fg.bp
      implicit real*8(a-h,o-z)
      parameter(mnp=35000)
      parameter(mne=65000)
      parameter(mnv=70)
      parameter(mnope=9)
      parameter(mnond=1000)

      dimension x(mnp),y(mnp),dp(mnp),nm(mne,3)
      dimension nond(mnope),iond(mnope,mnond)

      open(13,file='gen_fg.in',status='old')
      read(13,*)iob
      close(13)

      write(*,*) iob

      open(14,file='hgrid.gr3',status='old')
      open(12,file='fg.bp')
      read(14,*)
      read(14,*) ne,np
      if(ne>mne.or.np>mnp) then
        write(11,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif
      
      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*) j,l,(nm(i,k),k=1,3)
      enddo !ii

!     Open bnds
      read(14,*) nope
      if(nope>mnope) then
        write(11,*) 'nope > mnope' 
        stop
      endif

      read(14,*) neta
      ntot=0
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) then
          write(11,*) 'nond(k) > mnond'
          stop
        endif
        do i=1,nond(k)
          read(14,*) iond(k,i)
        enddo
        if(iond(k,1)==iond(k,nond(k))) then
          write(11,*)'Looped open bnd:',k
          stop
        endif
        ntot=ntot+nond(k)
      enddo

      if(neta/=ntot) then
        write(11,*)'neta /= total # of open bnd nodes',neta,ntot
        stop
      endif
      close(14)

!     Output
      write(12,*)'fg.bp'
      write(12,*)nond(iob)
      do i=1,nond(iob)
        nd=iond(iob,i)
        write(12,'(i9,3e17.8)')nd,x(nd),y(nd),dp(nd)
      enddo !i
 
      stop
      end
