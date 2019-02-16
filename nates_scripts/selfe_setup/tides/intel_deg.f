c***********************************************************************
c
c       program intel
c
c       interpolates to specified points the amplitudes and phases of
c       elevation computed by teanl.  Puts results in a format compatibl
c       with GR
c
c station file format:
c
c alpha
c number of points
c n x y  > repeat for number of points
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      PARAMETER (nnodes = 46000, nels = 82000, npts = 2000)
      DIMENSION next(nnodes),xord(nnodes),yord(nnodes),depth(nnodes)
      DIMENSION n(nels),icon(nels,3),nn(npts),xp(npts),yp(npts),
     +          nelem(npts)
      DIMENSION elmod(nnodes),phase(nnodes)
      CHARACTER*14 unit3,unit14,unit15,unit16,unit12
      CHARACTER*40 alphid

      pi = acos(-1.0)
      fac = 180./pi
C
C read grid file name
C
      WRITE (*,9030)
   10 READ (*,9000) unit3
      OPEN (unit=3,file=unit3,status='old',iostat=ist)
      IF (ist.EQ.0) GO TO 20
      WRITE (*,9020)
      GO TO 10

C
C read .sta file name
C
   20 WRITE (*,9010)
   30 READ (*,9000) unit12
      OPEN (unit=12,file=unit12,status='old',iostat=ist)
      IF (ist.EQ.0) GO TO 60
      WRITE (*,9020)
      GO TO 30

C
C read .tct file name
C
   60 WRITE (*,9040)
   70 READ (*,9000) unit15
      OPEN (unit=15,file=unit15,status='old',iostat=ist)
      IF (ist.EQ.0) GO TO 80
      WRITE (*,9020)
      GO TO 70

C
C read grid
C
   80 READ (3,9000) alphid
      READ (3,*) nmel,nmnp
      DO 90 i = 1,nmnp
          READ (3,*) next(i),xord(i),yord(i),depth(i)
   90 CONTINUE
      DO 100 i = 1,nmel
          READ (3,*) n(i),ijk, (icon(i,j),j=1,ijk)
  100 CONTINUE
C
C read station file
C
      READ (12,9000) alphid
      READ (12,*) npi
C
C read points from station file
C
          DO 110 i = 1,npi
              READ (12,*) nn(i),xp(i),yp(i)
              print *, nn(i),xp(i),yp(i)
  110     CONTINUE
          CLOSE (unit=12)
C
C find the element containing each point
C
C note, might be more efficient to reverse the
C order of the loops, here
C
          ncnt = 0
          DO 130 j = 1,npi
              x = xp(j)
              y = yp(j)
C
C locate element containing this point
C
              DO 120 i = 1,nmel
                  nel = n(i)
                  n1 = icon(nel,1)
                  n2 = icon(nel,2)
                  n3 = icon(nel,3)
                  x1 = xord(n1)
                  x2 = xord(n2)
                  x3 = xord(n3)
                  y1 = yord(n1)
                  y2 = yord(n2)
                  y3 = yord(n3)
                  a1 = x2*y3 - x3*y2
                  b1 = y2 - y3
                  c1 = x3 - x2
                  a2 = x3*y1 - x1*y3
                  b2 = y3 - y1
                  c2 = x1 - x3
                  a3 = x1*y2 - x2*y1
                  b3 = y1 - y2
                  c3 = x2 - x1
                  a = a1 + a2 + a3
                  CALL belel(x,y,n1,n2,n3,x1,x2,x3,y1,y2,y3,c1,c2,c3,b1,
     +                       b2,b3,k)

                  if (i.eq.24061.and.j.eq.1) then
                  print *, i,j,((x-x1)*b3 + (y-y1)*c3),
     &                     (x-x1)*b3,(y-y1)*c3
                  print *, i,j,((x-x2)*b1 + (y-y2)*c1),
     &                     (x-x2)*b1,(y-y2)*c1
                  print *, i,j,((x-x3)*b2 + (y-y3)*c2),
     &                     (x-x3)*b2,(y-y3)*c2
                  endif

                  IF (k.EQ.1) THEN
                      nelem(j) = n(i)
                      ncnt = ncnt + 1
                      GO TO 130

                  END IF

  120         CONTINUE
              WRITE(*,*)'No element found for point #', j
  130     CONTINUE
          if (ncnt .eq. 0) then
              stop 'no elements found for any point'
          endif

C
C open results file
C
      WRITE (*,9060)
  160 READ (*,9000) unit14
      OPEN (unit=14,file=unit14,status='unknown',iostat=ist)
      IF (ist.EQ.0) GO TO 170
      WRITE (*,9020)
      GO TO 160

C
C read data from .tct file
C
  170 READ (15,9000) alphid
      READ (15,*) nfreq
      WRITE (14,9000) alphid
      WRITE (14,*) npi
      WRITE (14,*) nfreq
      DO 210 j = 1,nfreq
c
c itest - 1 -> frequency used, 0 -> not used so skip
c
          READ (15,*) omega,itest
          if (itest.eq.0) then
             read(15,*)
             goto 210
          endif
          WRITE (14,*) omega
          READ (15,9000) alphid
          WRITE (14,9000) alphid
C
C read elevations
C
          DO 180 i = 1,nmnp
              READ (15,*) ni,elmod(i),phase(i)
  180     CONTINUE
C
C skip velocities
C
cepm          DO 181 i = 1,nmnp
cepm              READ (15,*)
cepm  181     CONTINUE
C
C interpolate
C
          DO 200 i = 1,npi
              x = xp(i)
              y = yp(i)
              nel = nelem(i)
              n1 = icon(nel,1)
              n2 = icon(nel,2)
              n3 = icon(nel,3)
              x1 = xord(n1)
              x2 = xord(n2)
              x3 = xord(n3)
              y1 = yord(n1)
              y2 = yord(n2)
              y3 = yord(n3)
              a1 = x2*y3 - x3*y2
              b1 = y2 - y3
              c1 = x3 - x2
              a2 = x3*y1 - x1*y3
              b2 = y3 - y1
              c2 = x1 - x3
              a3 = x1*y2 - x2*y1
              b3 = y1 - y2
              c3 = x2 - x1
              a = a1 + a2 + a3
              rl1 = (a1+b1*x+c1*y)/a
              rl2 = (a2+b2*x+c2*y)/a
              rl3 = (a3+b3*x+c3*y)/a
c
c use A=elmod*cos(phase)
c and B=elmod*sin(phase) for interpolating amp and phase
c
              atmp = rl1*elmod(n1)*cos(phase(n1)) + 
     +               rl2*elmod(n2)*cos(phase(n2)) + 
     +               rl3*elmod(n3)*cos(phase(n3)) 
              btmp = rl1*elmod(n1)*sin(phase(n1)) + 
     +               rl2*elmod(n2)*sin(phase(n2)) + 
     +               rl3*elmod(n3)*sin(phase(n3)) 
              tmp = atan2(btmp,atmp)
              an = atmp/cos(tmp)
C degrees
C             phn = amod(tmp*fac+720.0,360.0)
C radians
              phn = tmp
cepm              WRITE (14,9090) nn(i),an,phn*180.0/pi
              if (phn.lt.0.0) phn=phn+2.0*pi
              WRITE (14,*) an,phn*180.0/pi
  200     CONTINUE
  210 CONTINUE

 9000 FORMAT (a)
 9010 FORMAT ('Name of file with location of stations: ',$)
 9020 FORMAT ('File does not exist, try again: ',$)
 9030 FORMAT ('Grid file : ',$)
 9040 FORMAT ('TEANL file : ',$)
 9050 FORMAT ('New station file : ',$)
 9060 FORMAT ('Output file: ',$)
 9070 FORMAT (a10)
 9080 FORMAT (f16.10)
 9090 FORMAT (i5,4f16.8)
 9100 FORMAT (' Station ',i2,' does not belong to any element ',i3)
 9110 FORMAT (
     +       ' Find numbers of elements containing the stations (y=1',
     +       ',n=0) ? ',$)
 9120 FORMAT (2f12.2,2i6)
      END

      SUBROUTINE belel(x,y,n1,n2,n3,x1,x2,x3,y1,y2,y3,c1,c2,c3,b1,b2,b3,
     +                 k)
      implicit real*8(a-h,o-z)
c
c       tests if the point of global coordinates (x,y) belongs to
c       non-isoparametric triangular element j.
c       uses vector product test, if the test is positive, ind=1
c       otherwise ind=0
c
      k = 1
      f = (x-x1)*b3 + (y-y1)*c3
      IF (f.LT.0.) GO TO 10
      f = (x-x2)*b1 + (y-y2)*c1
      IF (f.LT.0.) GO TO 10
      f = (x-x3)*b2 + (y-y3)*c2
      IF (f.LT.0.) GO TO 10
      RETURN

   10 k = 0
      RETURN

      END
