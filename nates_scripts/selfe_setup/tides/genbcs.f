      program genbcs

      implicit real*8(a-h,o-z)
      character*24 title
      dimension x(1000000),y(1000000)

      open(11,file='fort.14.nos8',status='old')
      open(12,file='fort.14.sta',status='unknown')

      read(11,'(a24)') title
      read(11,*) ne,np
      do i=1,np
        read(11,*) node, x(i), y(i), depth
      end do
      do i=1,ne
        read(11,*) nel,ntype,n1,n2,n3
      end do
      read(11,*) nbndo
      read(11,*) nonodes
      read(11,*) nonodes1
      write(12,*) 'edwilla_bndry'
      write(12,*) nonodes1
      do i=1,nonodes1
        read(11,*) nbnum
        write(12,*) i,x(nbnum), y(nbnum)
      end do

      end
