!===============================================================================
! Read in rank-specific hotstart outputs and combine them into hotstart.in.
! Gobal-local mappings are read in from separate files.

! Inputs:
!        rank-specific hotstart files; 
!        local_to_global_*; 
!        screen: nproc, ntracers; it_char (iteration #)
! Output: hotstart.in (unformatted binary). This format is different
!         between Intel and AMD!
!
!  ifort -Bstatic -O3 -assume byterecl -o combine_hotstart3 combine_hotstart3.f90

! Revisions: non-hydrostatic pressure added (v3.0a up).
!===============================================================================

program combine_hotstart1
!-------------------------------------------------------------------------------
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
  character(12) :: it_char
  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  allocatable ner(:),npr(:),nsr(:)
  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trel0(:,:,:),trel(:,:,:)
!-------------------------------------------------------------------------------
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

!  open(10,file='combine_hotstart1.in',status='old')
  write(*,*) 'Input nproc, ntracers:'
  read(*,*) nproc,ntracers
  write(*,*) 'Input iteration # of the hotstart file (before _00*_hotstart) :'
  read(*,'(a)')it_char 

  allocate(ner(0:nproc-1),npr(0:nproc-1),nsr(0:nproc-1),stat=istat) 
  if(istat/=0) stop 'Allocation error'

! Read mapping info
  fdb='local_to_global_0000'; fdb=adjustl(fdb)
  lfdb=len_trim(fdb)
  mxner=0; mxnpr=0; mxnsr=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*)ne_global,np_global,nvrt,ntmp
    if(ntmp/=nproc) then
      print*, '# of proc mistmatch!',ntmp,nproc,ne_global,np_global,nvrt
      stop
    endif
    read(10,*) 
    read(10,*)ner(irank); do i=1,ner(irank); read(10,*); enddo;
    read(10,*)npr(irank); do i=1,npr(irank); read(10,*); enddo;
    read(10,*)nsr(irank); do i=1,nsr(irank); read(10,*); enddo;
    mxner=max0(mxner,ner(irank))
    mxnpr=max0(mxnpr,npr(irank))
    mxnsr=max0(mxnsr,nsr(irank))
    close(10)
  enddo !irank

  allocate(ielg(0:nproc-1,mxner),iplg(0:nproc-1,mxnpr),islg(0:nproc-1,mxnsr), &
          stat=istat)
  if(istat/=0) stop 'Allocation error (3)'

  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*); read(10,*)
    read(10,*)ner(irank)
    do i=1,ner(irank)
      read(10,*)j,ielg(irank,i)
    enddo
    read(10,*)npr(irank)
    do i=1,npr(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)nsr(irank)
    do i=1,nsr(irank)
      read(10,*)j,islg(irank,i)
    enddo
    close(10)
  enddo !irank

! Compute global side #
!  ne_global=0
!  np_global=0
  ns_global=0
  do irank=0,nproc-1
!    do i=1,ner(irank)
!      ne_global=max0(ne_global,ielg(irank,i))
!    enddo !i 
!    do i=1,npr(irank)
!      np_global=max0(np_global,iplg(irank,i))
!    enddo !i 
    do i=1,nsr(irank)
      ns_global=max0(ns_global,islg(irank,i))
    enddo !i 
  enddo !irank
  print*, 'Global quantities:',ne_global,np_global,ns_global

  allocate(idry_e(ne_global),we(nvrt,ne_global),tsel(2,nvrt,ne_global), &
           idry_s(ns_global),su2(nvrt,ns_global),sv2(nvrt,ns_global), &
           tsd(nvrt,ns_global),ssd(nvrt,ns_global), &
           idry(np_global),eta2(np_global),tnd(nvrt,np_global),snd(nvrt,np_global), &
           tem0(nvrt,np_global),sal0(nvrt,np_global),q2(np_global,nvrt), &
           xl(np_global,nvrt),dfv(np_global,nvrt),dfh(np_global,nvrt), &
           dfq1(np_global,nvrt),dfq2(np_global,nvrt),qnon(nvrt,np_global),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'
  if(ntracers==0) then
    allocate(trel0(1,1,1),trel(1,1,1),stat=istat)
  else
    allocate(trel0(ntracers,nvrt,ne_global),trel(ntracers,nvrt,ne_global),stat=istat)
  endif
  if(istat/=0) stop 'Allocation error (3)'
!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------

  ! Open file
  it_char=adjustl(it_char)  !place blanks at end
  it_len=len_trim(it_char)  !length without trailing blanks
  fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb);

  !Gather all ranks
  do irank=0,nproc-1
    !Open input file
    fgb2=fgb
    write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
    ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ner(irank)+(8*nvrt+1)*nsr(irank)+(3+22*nvrt)*npr(irank))
    open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
    read(36,rec=1)time,it,ifile,(idry_e(ielg(irank,i)),(we(j,ielg(irank,i)),tsel(1:2,j,ielg(irank,i)), &
     &(trel0(l,j,ielg(irank,i)),trel(l,j,ielg(irank,i)),l=1,ntracers),j=1,nvrt),i=1,ner(irank)), &
     &(idry_s(islg(irank,i)),(su2(j,islg(irank,i)),sv2(j,islg(irank,i)), &
     &tsd(j,islg(irank,i)),ssd(j,islg(irank,i)),j=1,nvrt),i=1,nsr(irank)), &
     &(eta2(iplg(irank,i)),idry(iplg(irank,i)),(tnd(j,iplg(irank,i)),snd(j,iplg(irank,i)), &
     &tem0(j,iplg(irank,i)),sal0(j,iplg(irank,i)),q2(iplg(irank,i),j),xl(iplg(irank,i),j), &
     &dfv(iplg(irank,i),j),dfh(iplg(irank,i),j),dfq1(iplg(irank,i),j), &
     &dfq2(iplg(irank,i),j),qnon(j,iplg(irank,i)),j=1,nvrt),i=1,npr(irank)) 
    close(36)
  enddo !irank

  print*, 'time,it,ifile:',time,it,ifile

! Output
  open(36,file='hotstart.in',form='unformatted',status='replace')
  write(36) time,it,ifile
  do i=1,ne_global
    write(36) i,idry_e(i),(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
  enddo !i
  do i=1,ns_global
    write(36) i,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
  enddo !i
  do i=1,np_global
    write(36) i,eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
             dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
  enddo !i
  close(36)

  stop
end program combine_hotstart1
