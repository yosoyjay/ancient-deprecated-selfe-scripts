      program ecp

      real*8 lambda_0,phi_1,x,y,lat,long,depth,pie
      integer i,np,ne,node,nelem,ntype,n1,n2,n3
      character*24 title

      open(11,file='fort.14.ll',status='old')
      open(12,file='fort.14.nos8',status='unknown')

      pie=3.14159265359
      lambda_0=-140.0*pie/180.0
      phi_1=50.0*pie/180.0

      read(11,'(a24)') title
      read(11,*) ne,np
      write(12,'(a24)') title
      write(12,*) ne,np
      do i=1,np
        read(11,*) node,long,lat,depth
        lat=lat*pie/180.0
        long=long*pie/180.0
        call cpp(x,y,long,lat,lambda_0,phi_1)
        write(12,10) node,x,y,depth
      end do
      do i=1,ne
        read(11,*) nelem,ntype,n1,n2,n3        
        write (12,*) nelem,ntype,n1,n2,n3
      end do
      
   10 format(i6,1x,f16.7,f16.7,1x,f16.7)
      end


      subroutine cpp(x,y,long,lat,lambda_0,phi_1)
      real*8 x,y,long,lat,lambda_0,phi_1
      r=6378206.4
      x=r*(long-lambda_0)*cos(phi_1)
      y=lat*r
      return
      end
