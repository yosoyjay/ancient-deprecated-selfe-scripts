!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                	!
!                                ELFE:                                      		!
!         A Three-Dimensional Baroclinic Model for Unstructured Grids         		!
!		        Version 1.5e (July 2006)                     	      		!
!                                                                             		!
!                 Center for Coastal and Land-Margin Research                 		!
!             Department of Environmental Science and Engineering             		!
!                   OGI School of Science and Engineering,		      		!
!	              Oregon Health & Science University		      		!
!                       Beaverton, Oregon 97006, USA                            	!
!								 	      		!
!                   Scientific direction: Antonio Baptista                    		!
!                   Code development: Joseph Zhang                           		!
!											!
!	        Copyright 2003-2004 Oregon Health and Science University		!
!  		               All Rights Reserved					!
!    									      		!
! 	The heat exchange module makes use of the bulk aerodynamic surface flux 	!
!	algorithm introduced by Zeng et al (1998), and the polynomial fits to 		!
!	saturation vapor pressure of Flatau et al (1992):				! 
!	Zeng, X., M. Zhao, and R. E. Dickinson, 1998:  Intercomparison of bulk		!
!	aerodynamic algorithms for the computation of sea surface fluxes using		!
!	TOGA COARE and TAO data.  J. Clim., 11, 2628-2644.				!
!	Flatau, P. J., R. L. Walko and W. R. Cotton, 1992:  Polynomial fits to		!
!	saturation vapor pressure.  J. Appl. Meteor., 31, 1507-1513.			!
!											!
!	Attenuation of solar radiation (and solar heating) within the water column	!
!	is based upon the expression given by Paulson and Simpson (1977), for the	!
!	water types defined by Jerlov (1968):						!
!	Jerlov, N. G., Optical Oceanography, Elsevier, 1968.				!
!	Paulson, C. A., and J. J. Simpson, Irradiance measurements in the upper		!
!	ocean, J. Phys. Oceanogr., 7, 952-956, 1977.					!
!											!
!	In addition, the module must be linked with netcdf library.
!
!       The GOTM option was taken from gotm.net.
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             		!

!...  Data type consts
      module kind_par
        implicit none
        integer, parameter :: sng_kind1=4
        integer, parameter :: dbl_kind1=8
        real(kind=dbl_kind1), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in global
      end module kind_par


!...  definition of variables
!...
!
!************************************************************************
!     			mnp < mne < mns					*
!************************************************************************
!
      module global
        implicit none
        integer, parameter :: sng_kind=4
        integer, parameter :: dbl_kind=8

!...  	Dimensioning parameters
        integer, parameter :: mnp=30000
        integer, parameter :: mne=60000
        integer, parameter :: mns=90000
        integer, parameter :: mnv=54
!       user-defined tracer part
        integer, parameter :: ntracers=0
!       end user-defined tracer part
        integer, parameter :: mntr=max(2,ntracers)
        integer, parameter :: mnei=20 !neighbor
        integer, parameter :: mne_kr=100 !max. # of elements used in Kriging
        integer, parameter :: mnei_kr=6 !max. # of pts used in Kriging
        integer, parameter :: mnope=20 !# of open bnd segements
        integer, parameter :: mnond=1000 !max. # of open-bnd nodes on each segment
        integer, parameter :: mnland=220 !# of land bnd segements
        integer, parameter :: mnlnd=10000 !max. # of land nodes on each segment
        integer, parameter :: mnbfr=15 !# of forcing freqs.
        integer, parameter :: itmax=5000 !# of iteration for itpack solvers used for dimensioning
        integer, parameter :: nwksp=6*mnp+4*itmax !available work space for itpack solvers
        integer, parameter :: nbyte=4
        integer, parameter :: mnout=100 !max. # of output files
        integer, parameter :: mirec=1109000000 !max. record # to prevent output ~> 4GB
        real(kind=dbl_kind), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in kind_par

!...  	Important variables
        integer :: np,ne,ns,nvrt,ivcor,itheta_1,itheta_2,kz,nsig,imm,kr_co,indvel,ihconsv,isconsv
        real(kind=dbl_kind) :: h0,q2min,rho0,dt,pi,theta_b,theta_f,h_c,tempmin,tempmax,saltmin,saltmax,vis_coe1,vis_coe2
        real(kind=dbl_kind) :: h_s,s_con1 !hyperbolic functions used in ivcor=2
!	Consts. used in GLS closure
        character(len=2) :: mid,stab
        real(kind=dbl_kind) :: ubd0,ubd1,ubd2,ubd3,ubd4,ubd5,ubs0,ubs1,ubs2,ubs4,ubs5,ubs6, &
     &a2_cm03,schk,schpsi

!...    Output handles
        character(len=48) :: start_time,version,data_format='DataFormat v5.0'
!	' (ylz: for appearance)
        character(len=12) :: ifile_char
        character(len=48), dimension(mnout) :: outfile,variable_nm,variable_dim
        integer :: nrec,nspool,igmp,noutgm,ifile,noutput,ifort12(100)
        integer, dimension(mnout) :: ichan,irec,iof
!        real(kind=dbl_kind), dimension(mnout) :: vpos
        
! evm   character buffers for binary type conversions
        integer :: iwrite
        character(len=48) :: a_48
        character(len=16) :: a_16
        character(len=8)  :: a_8
        character(len=4)  :: a_4


!...    1D arrays
        integer :: kfp(mnp) !only for sflux routines
        integer :: kbp(mnp),kbs(mns),kbe(mne),kbp00(mnp)
        integer :: nne(mnp),nnp(mnp),idry(mnp),idry_s(mns),idry_e(mne),idry_e0(mne), &
     &isbnd(mnp),isbs(mns),iback(mnp),interpol(mne),ie_kr(mne),lqk(mne),krvel(mne)
        real(kind=dbl_kind), dimension(mnp) :: x,y,dp,hmod,eta1,eta2,xlmin2,xlon,ylat,bdef,bdef1,bdef2,dp00
        real(kind=dbl_kind), dimension(mne) :: area,radiel,xctr,yctr,dpe
        real(kind=dbl_kind), dimension(mns) :: snx,sny,distj,xcj,ycj,dps
        real(kind=dbl_kind), dimension(mnv) :: sigma,cs,dcs,ztot
        real(kind=dbl_kind) :: decorrel(mne_kr)

!...    2D and higher arrays
        integer :: nm(mne,3),nx(3,2),ic3(mne,3),ine(mnp,mnei),js(mne,3),is(mns,2), &
         &isidenode(mns,2),inp(mnp,mnei),iself(mnp,mnei),isidenei(mns,2),isidenei2(mns,4), &
         &itier_nd(mne_kr,0:mnei_kr)

        real(kind=dbl_kind) :: ssign(mne,3),z(mnv,mnp),zs(mnv,mns),ze(mnv,mne),su2(mnv,mns),sv2(mnv,mns), &
     &tem0(mnv,mnp),sal0(mnv,mnp),tnd(mnv,mnp),snd(mnv,mnp),tsd(mnv,mns),ssd(mnv,mns), &
     &prho(mnp,mnv),q2(mnp,mnv),xl(mnp,mnv),we(mnv,mne),sig_t(mnp,mnv),side_ac(mns,2,2),side_x(mns,2),dfh(mnp,mnv)
!       Note: q2 is TKE as in GLS model (=u_i*u_i/2)

        real(kind=dbl_kind), dimension(mnv,mnp) :: uu2,vv2,ww2
        real(kind=dbl_kind), dimension(mnv,mne,3) :: ufg,vfg
        
        real(kind=dbl_kind) :: dl(mne,3,2),akrmat_nd(mne_kr,mnei_kr+3,mnei_kr+3)
!       Arrays used in transport routine
        real(kind=dbl_kind) :: tsel(mnv,mne,2) !S,T at elements and half levels for upwind scheme
        real(kind=dbl_kind) :: tr_el(mnv,mne,mntr) !tracer converntration @ prism center; used as temp. storage
        real(kind=dbl_kind) :: bdy_frc(mnv,mne,mntr) !body force at prism center Q_{i,k}
        real(kind=dbl_kind) :: flx_sf(mne,mntr) !surface b.c. \kappa*dC/dz = flx_sf (at element center)
        real(kind=dbl_kind) :: flx_bt(mne,mntr) !bottom b.c.

      end module global

!...  Main program
      program elfe
      use global
#ifdef USE_GOTM
      use turbulence, only: init_turbulence, do_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
      use mtridiagonal, only: init_tridiagonal
      
#endif
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

!...  Output handles
      real(kind=sng_kind) :: floatout,floatout2,st,en
      character(len=12) :: it_char
      character(len=40) :: date,timestamp
      character(len=2) :: tvd_mid,flimiter,tvd_mid2,flimiter2
      logical :: up_tvd

!...  Geometry
      dimension xlon_e(mne),ylat_e(mne),cwidth(mnope)
      dimension sigmap(mnv,10),sigma_prod(mnv,mnv,-4:4)
      dimension icolor1(mnp),icolor2(mns),ifront(mnei_kr),ifront2(mnei_kr)
      dimension akr(mnei_kr+3,mnei_kr+3) !,bkr(mnei_kr+3,1)
      dimension akrp((mnei_kr+3)*(mnei_kr+4)/2),ipiv(mnei_kr+3),work4(mnei_kr+3)
!      dimension akrmat_nd(mne_kr,mnei_kr+3,mnei_kr+3),akrmat_sd(mne_kr,mnei_kr+3,mnei_kr+3)

!...  Boundary forcings
      dimension nond(mnope),iond(mnope,mnond),nlnd(mnland),ilnd(mnland,mnlnd)
      dimension iettype(mnope),ifltype(mnope),itetype(mnope),isatype(mnope),tobc(mnope),sobc(mnope),itrtype(mnope)
      dimension tamp(mnbfr),tnf(mnbfr),tfreq(mnbfr),jspc(mnbfr),tear(mnbfr)
      dimension amig(mnbfr),ff(mnbfr),face(mnbfr)
      dimension emo(mnope,mnond,mnbfr),efa(mnope,mnond,mnbfr)
      dimension vmo(mnope,mnbfr),vfa(mnope,mnbfr)
      dimension eth(mnope,mnond),tth(mnope,mnond,mnv),sth(mnope,mnond,mnv),qthcon(mnope),ath(mnope),trth(mnope,ntracers)
      dimension uth(mns,mnv),vth(mns,mnv),uthnd(mnope,mnond,mnv),vthnd(mnope,mnond,mnv)
      dimension eta_mean(mnp),atd(mnope),z_r(mnv),tem1(mnv),sal1(mnv)

!...  Flow arrays
      dimension ptbt(mns,mnv,4),sdbt(mns,mnv,4),bubt(mne,2) !ptbt dimension enlarged for flux limiters
!      dimension x3bt(mnp,mnv,3),nelvbt(mnp,mnv,2),x3bt2(3,mnv,3),nelvbt2(3,mnv,2)
      dimension out3(mnv,3),out2(12)
      dimension windx1(mnp),windy1(mnp),windx2(mnp),windy2(mnp)
      dimension windx(mnp),windy(mnp),tau(mnp,2),iadv(mnp),nsubd(mnp),windfactor(mnp)
      dimension pr1(mnp),airt1(mnp),shum1(mnp),pr2(mnp),airt2(mnp),shum2(mnp),pr(mnp)
      dimension sflux(mnp),srad(mnp),tauxz(mnp),tauyz(mnp)
      dimension fluxsu(mnp),fluxlu(mnp),hradu(mnp),hradd(mnp)
      dimension chi(mns),cori(mns),Cd(mns),Cdp(mnp),rough(mns),rough_p(mnp)
      dimension dfv(mnp,mnv),dfz(2:mnv),dzz(2:mnv) !,dz2(mnv)
      dimension hvis(mnv,mns),d2u(mnv,mns),d2v(mnv,mns),horcon(mns)
      dimension icoef(mnp+1),jcoef(mnp*(mnei+1)),e2coef(mnp*(mnei+1)),qel(mnp)
      dimension sparsem(mnp,0:mnei),elbc(mnp),imap(mnp),qel2(mnp),eta3(mnp)
      dimension hhat(mns),bigu(mns,2),ghat1(mne,2),sne(mnv,3),area_e(mnv)
      dimension bcc(mns,mnv,2),hp_int(mnv,mne,2),ctmp(0:mnv) !hp_int indices reversed
      dimension ibt_p(mnp),ibt_s(mns),t_nudge(mnp),s_nudge(mnp),dr_ds(mnp,mnv)
      dimension fun_lat(mnp,0:2),etp(mnp) !fun_lat_e(mne,0:2)
      dimension dav(mnp,2)
      dimension elevmax(mnp) !max. elev. at nodes for all steps for tsunami
      dimension fluxprc(mnp),fluxevp(mnp)
      real(kind=dbl_kind),dimension(0:mnv) :: h1d,SS1d,NN1d !,num1d,nuh1d
      real(kind=sng_kind), dimension(mnp,mnv) :: tnd_nu1,snd_nu1,tnd_nu2,snd_nu2,tnd_nu,snd_nu
      dimension trel0(mnv,mne,ntracers),trel(mnv,mne,ntracers),tr_nd(mnv,mnp,ntracers)

!...  Wild-card arrays
      dimension nwild(mne+12),nwild2(mne),swild(mnp+mnv+12+ntracers),swild2(mnv,10) !swild2 dimension must match that in vinter()
      dimension swild5(3,2),swild6(4,2),swild4(mns,mnv,2),swild7(2,2,2)

!...  Solver arrays for TRIDAG
      dimension alow(mnv),bdia(mnv),cupp(mnv),rrhs(mnv,100),soln(mnv,100),gam(mnv) !"100" in rrhs & soln must match tridag()

!     MY-G turbulence closure arrays
      dimension diffmax(mnp),diffmin(mnp),dfq1(mnp,mnv),dfq2(mnp,mnv),q2tmp(mnv),xltmp(mnv)
      dimension rzbt(mnv),shearbt(2:mnv),xlmax(mnv),cpsi3(2:mnv),cpsi2p(2:mnv),q2ha(2:mnv),xlha(2:mnv)
      dimension xlsc0(mnp)

!...  variables used by the itpack solvers
      dimension iwksp(3*mnp),wksp(nwksp),iparm(12),rparm(12)


!...
!...  First executible statement of Elfe  
!...
      if(mnp>=mne.or.mne>=mns) then
        write(*,*)'Make sure mnp < mne < mns'
        stop
      endif

!...  Tracer transport
      if(ntracers<0) then
        write(11,*)'Illegal ntracers:',ntracers
        stop
      endif

!...  Initialize arrays and variables
      isbnd=0
      isbs=0
      iback=0 !back-up flags for abnormal cases in S-coord.
      pr1=0; pr2=0; pr=0 !uniform pressure (the const. is unimportant)
      uth=-99; vth=-99; uthnd=-99; vthnd=-99; eta_mean=-99 !flags
      fluxsu00=0; srad00=0 !for nws/=3
      elevmax=-1.e34
!     tsel and trel for passing on to routine; use allocatable arrays later
      tsel=0; trel=0

!     for output
      airt1=0; shum1=0;  airt2=0; shum2=0; srad=0; fluxsu=0; fluxlu=0 
      hradu=0; hradd=0; sflux=0; windx=0; windy=0
      q2=0; xl=0 !for hotstart with itur/=3 only
!     Fort.12 flags
      ifort12=0

!...  define some constants and initial values
!...
      omega=7.29d-5 !angular freq. of earth rotation 
      rearth=6378206.4 !earth radius
      g=9.81
      rho0=1000. !ref. density for S=33 and T=10C
      pi=dacos(-1.0d0)
      shw=4184  !specific heat of pure water

      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(*,*)'nx wrong',i,j,nx(i,j)
            stop
          endif
        enddo !j
      enddo !i


!                                                                             *
!******************************************************************************
!                                                                             *
!			open input files				      *
!                                                                             *
!******************************************************************************
!                                                                             *

      open(14,file='hgrid.gr3',status='old')
      open(15,file='param.in',status='old')
      open(19,file='vgrid.in',status='old')

      open(11,file='fort.11') !fatal error message output
      open(12,file='fort.12') !non-fatal error message output
      open(16,file='mirror.out')
      open(10,file='total.dat') !output total mass etc.

      call date_and_time(date,timestamp)
      write(16,*)'Run begins at ',date,timestamp

!...  read the vertical layers information from vgrid.in
!...
      ivcor=2 !S only
      read(19,*) nvrt,kz,h_s !kz>=1
      if(nvrt>mnv.or.nvrt<3) then
        write(11,*)'nvrt > mnv or nvrt<4'
        stop
      endif
      if(kz<1.or.kz>nvrt-2) then
        write(11,*)'Wrong kz:',kz
        stop
      endif
      if(h_s<10) then
        write(11,*)'h_s needs to be larger:',h_s
        stop
      endif

!     # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz-1
       read(19,*)j,ztot(k)
       if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
         write(11,*)'z-level inverted:',k
         stop
        endif
      enddo !k
      read(19,*) !level kz       
!     In case kz=1, there is only 1 ztot(1)=-h_s
      ztot(kz)=-h_s

      nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
      read(19,*) !for adding comment "S levels"
      read(19,*)h_c,theta_b,theta_f
      if(h_c<5) then !large h_c to avoid 2nd type abnormaty
        write(11,*)'h_c needs to be larger:',h_c
        stop
      endif
      if(theta_b<0.or.theta_b>1) then
        write(11,*)'Wrong theta_b:',theta_b
        stop
      endif
      if(theta_f<=0) then 
        write(11,*)'Wrong theta_f:',theta_f 
        stop
      endif
!     Pre-compute constants
      s_con1=dsinh(theta_f)

      sigma(1)=-1 !bottom
      sigma(nsig)=0 !surface
      read(19,*) !level kz
      do k=kz+1,nvrt-1
        kin=k-kz+1
        read(19,*) j,sigma(kin)
        if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
          write(11,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
          stop
        endif
      enddo
      read(19,*) !level nvrt
      close(19)

!     Compute C(s) and C'(s)
      do k=1,nsig
        cs(k)=(1-theta_b)*dsinh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*(dtanh(theta_f*(sigma(k)+0.5))-dtanh(theta_f*0.5))/2/dtanh(theta_f*0.5)
        dcs(k)=(1-theta_b)*theta_f*dcosh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*theta_f/2/dtanh(theta_f*0.5)/dcosh(theta_f*(sigma(k)+0.5))**2
      enddo !k=1,nvrt

!...  Output some sample z-coordinates
      write(16,*)'---------------------------------------------'
      swild(1)=20; swild(2)=h_s/2; swild(3)=h_s; swild(4)=4000
      write(16,*)'h_c= ',h_c,' h_s=',h_s
      do i=1,4
        write(16,*)'Depth= ',swild(i)
        do k=kz,nvrt
          kin=k-kz+1
          hmod2=dmin1(swild(i),h_s)
          if(hmod2<=h_c) then
            zz=sigma(kin)*hmod2
          else
            zz=h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
          endif
          write(16,*)k,zz
        enddo !k
      enddo !i
      write(16,*)'---------------------------------------------'

!...  read unit 15 file
!...

      read(15,'(a48)') version
      read(15,'(a48)') start_time
      read(15,*) ipre !pre-processing flag (to output obs.out)
      read(15,*) nscreen
      read(15,*) iwrite !write mode
      read(15,*) imm !0: without bed deformation; 1: with bed deformation (e.g., tsunami)
!     For moving bed, the output is from original bottom to nvrt
      if(imm<0.or.imm>1) then
        write(11,*)'Unknown imm',imm
        stop
      endif
!     Initialize variables used in tsunami model (but bdef[1,2] and ibdef are available for all models)
      bdef=0 !total deformation
      ibdef=1 !# of time steps for deformation (deformation rate=0 when it>ibdef)
      if(imm==1) then !read in deformation at all nodes
        read(15,*) ibdef
        open(32,file='bdef.gr3',status='old') !connectivity part not used
        read(32,*)
        read(32,*)ntmp,np
        do i=1,np
          read(32,*)j,xtmp,ytmp,bdef(i) !total deformation
        enddo !i
        close(32)
      endif

      read(15,*) ihot
      if(ihot<0.or.ihot>2) then
        write(11,*)'Unknown ihot',ihot
        stop
      endif

      read(15,*) ics
      if(ics/=1.and.ics/=2) then
        write(11,*)'Unknown ics',ics
        stop
      endif

!...  Center of projection in degrees
      read(15,*) slam0,sfea0
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180

!...  Horizontal viscosity option
      read(15,*) ihorcon !=0 means all horcon=0 and no horcon.gr3 is needed
      
!...  Enough info to read unit 14 grid file
!...
      read(14,*) 
      read(14,*) ne,np
      if(ne>mne.or.np>mnp) then
        write(11,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif

      dpmax=0 !max. depth
      do i=1,np
        if(ics==1) then
          read(14,*) j,x(i),y(i),dp(i)
        else !=2
          read(14,*) j,xlon(i),ylat(i),dp(i)
          ylat(i)=ylat(i)*pi/180
          xlon(i)=xlon(i)*pi/180
          call cpp(x(i),y(i),xlon(i),ylat(i),slam0,sfea0)
        endif
        hmod(i)=dmin1(dp(i),h_s)
        if(dp(i)>dpmax) dpmax=dp(i)
      enddo !i=1,np
!     Save intial depth for bed deformation case
      dp00=dp

      if(ztot(1)>=-dpmax) then
        write(11,*)'1st z-level must be below max. depth:',dpmax
        stop
      endif

      do i=1,ne
        read(14,*) j,l,(nm(i,k),k=1,l)
        if(l/=3) then
          write(11,*)'SELFE cannot handle quads:',i
          stop
        endif
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(i)<=0) then
          write(11,*)'Negative area at',i
          stop
        endif
        radiel(i)=dsqrt(area(i)/pi)  !equivalent radius
!...    Derivatives of shape functions
        dl(i,1,1)=(y(n2)-y(n3))/2/area(i) !dL_1/dx
        dl(i,2,1)=(y(n3)-y(n1))/2/area(i) !dL_2/dx
        dl(i,3,1)=(y(n1)-y(n2))/2/area(i)
        dl(i,1,2)=(x(n3)-x(n2))/2/area(i) !dL_1/dy
        dl(i,2,2)=(x(n1)-x(n3))/2/area(i)
        dl(i,3,2)=(x(n2)-x(n1))/2/area(i)
      enddo !i=1,ne

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
          isbnd(iond(k,i))=k
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

!     Land bnds
      read(14,*) nland
      if(nland>mnland) then
        write(11,*) 'nland > mnland'
        stop
      endif

      read(14,*) nvel
      do k=1,nland
        read(14,*) nlnd(k)
        if(nlnd(k)>mnlnd) then
          write(11,*)'nlnd(k) > mnlnd',k,nlnd(k),mnlnd
          stop
        endif
        do i=1,nlnd(k)
          read(14,*) ilnd(k,i)
          if(isbnd(ilnd(k,i))==0) isbnd(ilnd(k,i))=-1 !overlap of open bnd
        enddo
      enddo !k=1,nland
      close(14)
!...  End fort.14

!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!     			Compute geometry 				      *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *

!...  compute the ball of elements and arrange in counter-clockwise fashion
      do i=1,np
        nne(i)=0
        nnp(i)=0
      enddo

      do i=1,ne
        do j=1,3
          nd=nm(i,j)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          ine(nd,nne(nd))=i
          iself(nd,nne(nd))=j !to be updated later
        enddo
      enddo

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,3
          ic3(i,j)=0 !index for bnd sides
          nd1=nm(i,nx(j,1))
          nd2=nm(i,nx(j,2))
          do k=1,nne(nd1)
            ie=ine(nd1,k)
            if(ie/=i.and.(nm(ie,1)==nd2.or.nm(ie,2)==nd2.or.nm(ie,3)==nd2)) ic3(i,j)=ie
          enddo !k
        enddo !j
      enddo !i

!     Re-arrange in counter-clockwise fashion
      do i=1,np
        if(isbnd(i)/=0) then !bnd ball
!         Look for starting bnd element
          icount=0
          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)
            if(ic3(ie,nx(id,2))==0) then
              icount=icount+1
              ine(i,1)=ie
              iself(i,1)=id
            endif
          enddo !j=1,nne(i)
          if(icount/=1) then
            write(11,*)'Illegal bnd node',i
            stop
          endif
        endif !bnd ball

!	Sequential search for the rest of elements
        nnp(i)=2
        inp(i,1)=nm(ine(i,1),nx(iself(i,1),1))
        inp(i,2)=nm(ine(i,1),nx(iself(i,1),2))
        do j=2,nne(i)
          new=ic3(ine(i,j-1),nx(iself(i,j-1),1))
          if(new==0) then
            write(11,*)'Incomplete ball',i
            stop
          endif
          ine(i,j)=new
          id=0
          do l=1,3
            if(nm(new,l)==i) id=l
          enddo !l
          if(id==0) then
            write(11,*)'Failed to find local index:',i,new
            stop
          endif
          iself(i,j)=id

          if(isbnd(i)==0.and.j==nne(i)) then !complete internal ball
!	    Check completeness
            if(nm(new,nx(id,2))/=inp(i,1)) then
              write(11,*)'Broken ball:',i
              stop
            endif
          else !one more node
            nnp(i)=nnp(i)+1
            if(nnp(i)>mnei) then
              write(11,*)'Too many neighbor nodes',i
              stop
            endif
            inp(i,nnp(i))=nm(new,nx(id,2))
          endif
        enddo !j=2,nne(i)
      enddo !i=1,np

!...  Check hanging nodes
!...
      ihang=0
      do i=1,np
        if(nne(i)==0) then
          ihang=1
          write(11,*)'Hanging node',i
        endif
      enddo
      if(ihang==1) then
        write(11,*)'Check fort.11 for hanging nodes'
        stop
      endif
      
!...  compute the sides information
!...
      ns=0 !# of sides
      do i=1,ne
        do j=1,3
          nd1=nm(i,nx(j,1))
          nd2=nm(i,nx(j,2))
          if(ic3(i,j)==0.or.i<ic3(i,j)) then !new sides
            ns=ns+1 
            if(ns>mns) then
              write(11,*)'Too many sides'
              stop
            endif
            js(i,j)=ns
            is(ns,1)=i
            isidenode(ns,1)=nd1
            isidenode(ns,2)=nd2
            xcj(ns)=(x(nd1)+x(nd2))/2
            ycj(ns)=(y(nd1)+y(nd2))/2
            dps(ns)=(dp(nd1)+dp(nd2))/2
            distj(ns)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
            if(distj(ns)==0) then
              write(11,*)'Zero side',ns
              stop
            endif
            thetan=datan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
            snx(ns)=dcos(thetan)
            sny(ns)=dsin(thetan)

            is(ns,2)=ic3(i,j) !bnd element => bnd side
!	    Corresponding side in element ic3(i,j)
            if(ic3(i,j)/=0) then !old internal side
              iel=ic3(i,j)
              index=0
              do k=1,3
                if(ic3(iel,k)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(*,*)'Wrong ball info',i,j
                stop
              endif
              js(iel,index)=ns
            endif !ic3(i,j).ne.0
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,3
      enddo !i=1,ne

      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

      do i=1,ne
        do j=1,3
          jsj=js(i,j)
          ssign(i,j)=(is(jsj,2)-2*i+is(jsj,1))/(is(jsj,2)-is(jsj,1))
        enddo
      enddo

      if(nscreen.eq.1) write(*,*) 'There are',ns,' sides in the grid...'
      write(16,*) 'There are',ns,' sides in the grid...'

!...  compute centers of each element & dpe
!...
      do i=1,ne
        xctr(i)=0
        yctr(i)=0
        dpe(i)=1.e10
        do j=1,3
          xctr(i)=xctr(i)+x(nm(i,j))/3
          yctr(i)=yctr(i)+y(nm(i,j))/3
          if(dpe(i)>dp(nm(i,j))) dpe(i)=dp(nm(i,j))
        enddo !j
      enddo !i=1,ne

!...  Compute open bnd sides and total lengths of each open bnd segments (for imposing discharge)
      if(ipre==1) open(32,file='obs.out')

      do i=1,nope
        cwidth(i)=0
        do j=1,nond(i)-1
          n1=iond(i,j)
          n2=iond(i,j+1)
          cwidth(i)=cwidth(i)+dsqrt((x(n2)-x(n1))**2+(y(n2)-y(n1))**2)
          if(inp(n1,1)==n2) then
            ie=ine(n1,1)
            id=iself(n1,1)
            isd=js(ie,nx(id,2))
          else if(inp(n1,nnp(n1))==n2) then
            ie=ine(n1,nne(n1))
            id=iself(n1,nne(n1))
            isd=js(ie,nx(id,1))
          else
            write(11,*)'Wrong bnd ball orientation:',n1,n2
            stop
          endif
          if(is(isd,2)/=0) then
            write(11,*)'Wrong bnd side:',n1,n2
            stop
          endif
          isbs(isd)=i

!	  Output obs.out
          if(ipre==1) write(32,*)i,isd,isidenode(isd,1),isidenode(isd,2)
        enddo !j
      enddo !i

!...  Compute neighborhood for internal sides for horizontal derivatives (viscosity)
!...  isidenei(ns,2): 2 neighboring elements of a side 
!...  isidenei2(ns,4): 4 neighboring sides of a side
      dxy_min=1.e34 !min. local x-coord.
      iabort=0 !abort flag
      loop18: do i=1,ns
        if(is(i,2)==0) cycle loop18

!       Internal sides
        node1=isidenode(i,1)
        node2=isidenode(i,2)
        do j=1,2
          ie=is(i,j)
          l0=lindex_s(i,ie)
          if(l0==0) then
            write(11,*)'Cannot find a side'
            stop
          endif
          nwild(2*j-1)=js(ie,nx(l0,1))
          nwild(2*j)=js(ie,nx(l0,2))
        enddo !j=1,2
        isidenei2(i,1:4)=nwild(1:4) !may be modified later

!       First part for hvis only
        if(ihorcon/=0) then !need to compute hvis later
!         Compute intersections 1 and 2
          x0=xcj(i); y0=ycj(i)
          x1=xcj(nwild(1)); y1=ycj(nwild(1))
          x2=xcj(nwild(2)); y2=ycj(nwild(2))
          x3=xcj(nwild(3)); y3=ycj(nwild(3))
          x4=xcj(nwild(4)); y4=ycj(nwild(4))
          call intersect4(x0,snx(i),x1,x2,y0,sny(i),y1,y2,iflag1,xin1,yin1,tt11,tt21) !pt 1
          call intersect4(x0,-snx(i),x3,x4,y0,-sny(i),y3,y4,iflag2,xin2,yin2,tt12,tt22) !pt 2
          if(iflag1/=1.or.iflag2/=1.or.tt11==0.or.tt12==0) then
            write(11,*)'Failed to find 1 or 2:',iflag1,iflag2
            stop
          endif

          if(tt21<0) then
            isidenei(i,1)=is(nwild(1),1)+is(nwild(1),2)-is(i,1)
          else if(tt21>1) then
            isidenei(i,1)=is(nwild(2),1)+is(nwild(2),2)-is(i,1)
          else !inside
            isidenei(i,1)=is(i,1)
          endif
          if(tt22<0) then
            isidenei(i,2)=is(nwild(3),1)+is(nwild(3),2)-is(i,2)
          else if(tt22>1) then
            isidenei(i,2)=is(nwild(4),1)+is(nwild(4),2)-is(i,2)
          else !inside
            isidenei(i,2)=is(i,2)
          endif
          if(isidenei(i,1)==0.or.isidenei(i,2)==0) then
            iabort=1
            write(11,*)'Bnd side reached:',node1,node2,isidenei(i,1:2)
            cycle loop18
          endif

          call area_coord2(isidenei(i,1),xin1,yin1,swild,ifl1)
          side_ac(i,1,1:2)=swild(1:2)
          call area_coord2(isidenei(i,2),xin2,yin2,swild,ifl2)
          side_ac(i,2,1:2)=swild(1:2)
          if(ifl1==1.or.ifl2==1) then
            iabort=1
            write(11,*)'Failed to extend stencil:',node1,node2,ifl1,ifl2,isidenei(i,1:2)
            write(11,*)xin1,yin1,xin2,yin2
            write(11,*)side_ac(i,1,1:2)
            write(11,*)'-----------------------------------------------'
            cycle loop18
          endif

!         local x-coord.
          side_x(i,1)=(xin1-x0)*snx(i)+(yin1-y0)*sny(i)
          side_x(i,2)=(xin2-x0)*snx(i)+(yin2-y0)*sny(i)
          if(side_x(i,1)>=0.or.side_x(i,2)<=0) then
            write(11,*)'x-coord. out of order:',side_x(i,1:2),node1,node2
            stop
          endif
          if(abs(side_x(i,1))<dxy_min) dxy_min=abs(side_x(i,1)) 
          if(abs(side_x(i,2))<dxy_min) dxy_min=abs(side_x(i,2)) 

!         Debug
!          if(isidenei(i,1)/=is(i,1).or.isidenei(i,2)/=is(i,2)) write(97,*)node1,node2,i
        endif !ihorcon/=0

!       For Shapiro filter
!       Check if pt "0" is inside
        x1=xcj(i); y1=ycj(i)
        x2=xcj(isidenei2(i,4)); y2=ycj(isidenei2(i,4))
        x3=xcj(isidenei2(i,1)); y3=ycj(isidenei2(i,1))
        rl14=dsqrt((x2-x3)**2+(y2-y3)**2)
        ar4=signa(x1,x2,x3,y1,y2,y3)
        x2=xcj(isidenei2(i,2)); y2=ycj(isidenei2(i,2))
        x3=xcj(isidenei2(i,3)); y3=ycj(isidenei2(i,3))
        rl23=dsqrt((x2-x3)**2+(y2-y3)**2)
        ar3=signa(x1,x2,x3,y1,y2,y3)
        if(ar3<=0.and.ar4<=0) then
          write(11,*)'Degenerate parallelogram'
          stop
        endif

!       Enlarge stencil if pt 0 is outside
        if(ar3<=0.or.ar4<=0) then
          if(ar3<=0) then
            if(is(isidenei2(i,2),2)==0.or.is(isidenei2(i,3),2)==0) then
              iabort=1
              write(11,*)node1,node2,', bnd side (3)'
              cycle loop18
            endif

            nwild(1)=2; nwild(2)=3
            do k=1,2
              id=isidenei2(i,nwild(k))
              ie2=is(id,1)+is(id,2)-is(i,k)
              l0=lindex_s(id,ie2)
              if(l0==0) then
                write(11,*)'Cannot find a side (9):',k
                stop
              endif
              isidenei2(i,nwild(k))=js(ie2,nx(l0,3-k))
            enddo !k
          endif !ar3
         
          if(ar4<=0) then
            if(is(isidenei2(i,1),2)==0.or.is(isidenei2(i,4),2)==0) then
              iabort=1
              write(11,*)node1,node2,', bnd side (4)'
              cycle loop18
            endif

            nwild(1)=1; nwild(2)=4
            do k=1,2
              id=isidenei2(i,nwild(k))
              ie2=is(id,1)+is(id,2)-is(i,k)
              l0=lindex_s(id,ie2)
              if(l0==0) then
                write(11,*)'Cannot find a side (8):',k
                stop
              endif
              isidenei2(i,nwild(k))=js(ie2,nx(l0,k))
            enddo !k
          endif !ar4
         
!         Check convexity of quad 1-4
          x1=xcj(isidenei2(i,1)); y1=ycj(isidenei2(i,1))
          x2=xcj(isidenei2(i,2)); y2=ycj(isidenei2(i,2))
          x3=xcj(isidenei2(i,3)); y3=ycj(isidenei2(i,3))
          x4=xcj(isidenei2(i,4)); y4=ycj(isidenei2(i,4))
          ar1=signa(x1,x2,x3,y1,y2,y3)
          ar2=signa(x1,x3,x4,y1,y3,y4)
          ar3=signa(x1,x2,x4,y1,y2,y4)
          ar4=signa(x2,x3,x4,y2,y3,y4)
          if(ar1<=0.or.ar2<=0.or.ar3<=0.or.ar4<=0) then
            iabort=1
            write(11,*)node1,node2,'  Concave quad '
            write(11,*)((isidenode(isidenei2(i,m),mm),mm=1,2),m=1,4)
            write(11,*)ar1,ar2,ar3,ar4
            write(11,*)'--------------------------------------------'
            cycle loop18
          endif

!         Check if pt "0" is inside
          x0=xcj(i); y0=ycj(i)
          ar1=signa(x1,x2,x0,y1,y2,y0)
          ar2=signa(x2,x3,x0,y2,y3,y0)
          ar3=signa(x3,x4,x0,y3,y4,y0)
          ar4=signa(x4,x1,x0,y4,y1,y0)
          if(ar1<=0.or.ar2<=0.or.ar3<=0.or.ar4<=0) then
            iabort=1
            write(11,*)node1,node2,'  pt outside quad '
            write(11,*)ar1,ar2,ar3,ar4
            write(11,*)((isidenode(isidenei2(i,m),mm),mm=1,2),m=1,4)
            write(11,*)'----------------------------------------'
            cycle loop18
          endif
        endif !pt 0 outside

      end do loop18 !i=1,ns

      if(iabort==1) then
        write(*,*)'Check fort.11 for problems in side neighborhood'
        stop
      endif
      write(*,*)'Min. local x-coord= ',dxy_min

!.... Output sidecenters.bp for ipre=1
      if(ipre==1) then
        close(32)
        open(32,file='sidecenters.bp')
        write(32,*) 'Sidegrid'
        write(32,*) ns
        do i=1,ns
          write(32,'(i10,2(1x,e20.12),1x,f12.4)') i,xcj(i),ycj(i),real(dps(i))
        enddo !i
        write(*,*)'Side info generated; bye'
        close(32)
        stop
      endif

      if(nscreen==1) write(*,*)'done computing geometry...'
      write(16,*)'done computing geometry...'

!...  Continue to read param.in
!...

!...  Implicitness for momentum
      read(15,*)thetai

!...  Baroclinic flags
      read(15,*) ibc,ibtp
      if(ibc/=0.and.ibc/=1) then
        write(11,*)'Unknown ibc'
        stop
      endif
      if(ibtp/=0.and.ibtp/=1) then
        write(11,*)'Unknown ibtp'
        stop
      endif

      if(ibc==0) then
        write(*,*)'You are using baroclinic model'
        write(16,*)'You are using baroclinic model'
        read(15,*) nrampbc,drampbc 
      else !ibc=1
        if(ibtp==0) then
          write(*,*)'Barotropic model without ST calculation'
          write(16,*)'Barotropic model without ST calculation'
        else !ibtp=1
          write(*,*)'Barotropic model with ST calculation'
          write(16,*)'Barotropic model with ST calculation'
        endif
      endif

!      read(15,*) !tempmin,tempmax,saltmin,saltmax
!      if(tempmin<0.or.tempmax>40.or.saltmin<0.or.saltmax>42) then
!        write(11,*)'Specified ST range invalid'
!        stop
!      endif
       tempmin=0; tempmax=40; saltmin=0; saltmax=42

      read(15,*) rnday

!...  dramp not used if nramp=0
      read(15,*) nramp,dramp
      if(nramp/=0.and.nramp/=1) then
        write(11,*)'Unknown nramp',nramp
        stop
      endif

      read(15,*) dt

!...  compute total number of time steps 
      ntime=rnday*86400/dt+0.5
      write(10,*)ntime
      write(10,'(a200)')'Time (hours), volume, mass, potential E, kinetic E, total E, friction loss (Joule), energy leak (Joule)'
!'

!...  input info on backtracking
!...
      read(15,*) !nsubfl !flag

!...  Advection flag for momentum eq.
      read(15,*) nadv !flag: 1-Euler; 2: R-K
      if(nadv<0.or.nadv>2) then
        write(11,*)'Unknown advection flag',nadv
        stop
      endif

      if(nadv==0) then
        open(42,file='adv.gr3',status='old')
        read(42,*)
        read(42,*) !ne,np
        do i=1,np
          read(42,*)j,xtmp,ytmp,tmp
          iadv(i)=tmp
          if(iadv(i)<0.or.iadv(i)>2) then
            write(11,*)'Unknown iadv',i
            stop
          endif
        enddo
        close(42)
      else !nadv/=0
        iadv=nadv
      endif

!...  Tracking step
      read(15,*) dtb_max1,dtb_max2 !min. or max. sub-step for Euler and R-K if nadv=0; otherwise only dtb_max1 is used

!...  Minimum depth allowed
      read(15,*) h0 
      if(h0<=0) then
        write(11,*)'h0 must be positive'
        stop
      endif

!...  Bottom friction
      read(15,*) nchi
      if(nchi==0) then !read in drag coefficients
        open(32,file='drag.gr3',status='old')
        read(32,*)
        read(32,*) 
        do i=1,np
          read(32,*)j,xtmp,ytmp,Cdp(i)
          if(Cdp(i)<0) then
            write(11,*)'Negative bottom drag',Cdp(i)
            stop
          endif
        enddo
        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          Cd(i)=(Cdp(n1)+Cdp(n2))/2
        enddo
        close(32)
      else if(nchi==1) then !read in roughness in meters
        read(15,*) Cdmax !max. Cd
        open(32,file='rough.gr3',status='old')
        read(32,*)
        read(32,*)
        do i=1,np
          read(32,*)j,xtmp,ytmp,rough_p(i)
!          if(rough_p(i)<0) then
!            write(11,*)'Negative bottom roughness',rough_p(i)
!            stop
!          endif
        enddo !i
        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          sm=dmin1(rough_p(n1),rough_p(n2))
          if(sm<0) then
            rough(i)=sm !<0
          else !both non-negative
            rough(i)=(rough_p(n1)+rough_p(n2))/2 !>=0
          endif
        enddo !i
        close(32)
      else
        write(11,*)'Unknown nchi', nchi
        stop
      endif

!     Coriolis
      read(15,*) ncor
      if(abs(ncor)>1) then
        write(11,*)'Unknown ncor',ncor
        stop
      endif
      if(ncor==-1) then !lattitude
        read(15,*) tmp
        coricoef=2*omega*dsin(tmp/180*pi)
        cori=coricoef
      else if(ncor==0) then
        read(15,*) coricoef
        cori=coricoef
      else !ncor=1
        write(*,*)'Check slam0 and sfea0 as variable Coriolis is used'
        write(16,*)'Check slam0 and sfea0 as variable Coriolis is used'
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np
          read(32,*)j,xlon(i),ylat(i)
          xlon(i)=xlon(i)*pi/180
          ylat(i)=ylat(i)*pi/180
        enddo !i
        close(32)

        open(31,file='coriolis.out')
        fc=2*omega*dsin(sfea0)
        beta=2*omega*dcos(sfea0)
        do i=1,ns
          id1=isidenode(i,1)
          id2=isidenode(i,2)
          sphi=(ylat(id1)+ylat(id2))/2
          cori(i)=fc+beta*(sphi-sfea0)
          if(iwrite==0) then
            write(31,*)i,xcj(i),ycj(i),cori(i)
          else !evm
            write(31,"(a,i6,a,f16.9,a,f16.9,a,es22.14e3,a)",advance="no") &
     &           " ",i," ",xcj(i)," ",ycj(i)," ",cori(i),"\n"
          endif
        enddo !i=1,ns
        close(31)
      endif

!     Wind (nws=3: for conservation check; otherwise same as nws=2)
      read(15,*) nws,wtiminc
      if(nws<0.or.nws>3) then
        write(11,*)'Unknown nws',nws
        stop
      endif
      if(nws>0.and.dt>wtiminc) then
        write(11,*)'wtiminc < dt'
        stop
      endif

      if(nws>=2) then !CORIE mode; read in hgrid.ll
#ifndef USE_SFLUX
        write(11,*)'COIRE mode needs sflux routines'
        stop
#endif
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np
          read(32,*)j,xlon(i),ylat(i)
          xlon(i)=xlon(i)*pi/180
          ylat(i)=ylat(i)*pi/180
        enddo !i
        close(32)
      endif

      windfactor=1 !intialize for default
      if(nws>0) then
        read(15,*) nrampwind,drampwind
        read(15,*) iwindoff
        if(iwindoff/=0) then
          open(32,file='windfactor.gr3',status='old')
          read(32,*)
          read(32,*)ntmp,np
          do i=1,np
            read(32,*)j,xtmp,ytmp,windfactor(i) 
            if(windfactor(i)<0) then
              write(11,*)'Wind scaling factor must be positive:',i,windfactor(i)
              stop
            endif
          enddo !i
          close(32)
        endif
      endif

!     Heat and salt conservation flags
      read(15,*) ihconsv,isconsv
      if(ihconsv<0.or.ihconsv>1.or.isconsv<0.or.isconsv>1) then
        write(11,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
        stop
      endif
      if(isconsv/=0.and.ihconsv==0) then
        write(11,*)'Evap/precip model must be used with heat exchnage model'
!'
        stop
      endif
      if(ihconsv/=0.and.nws<2) then
        write(11,*)'Heat budge model must have nws>=2'
        stop
      endif
      if(ihconsv/=0) then
        write(16,*)'Warning: you have chosen a heat conservation model'
        write(16,*)'which assumes start time at 0:00 PST!'

#ifndef USE_SFLUX
        write(11,*)'Pls enable USE_SFLUX to use heat budget model'
        stop
#endif
#ifndef USE_NETCDF
        write(11,*)'Pls enable USE_NETCDF to use heat budget model'
        stop
#endif
      endif
   
      if(isconsv/=0) then 
#ifndef PREC_EVAP
        write(11,*)'Pls enable PREC_EVAP:',isconsv
        stop
!       USE_SFLUX and USE_NETCDF are definitely enabled in Makefile when isconsv=1
#endif
      endif

!...  Turbulence closure options
      read(15,*) itur
      if(itur<-2.or.itur>4) then
        write(11,*)'Unknown turbulence closure model',itur
        stop
      endif

      if(itur==0) then
        read(15,*) dfv0,dfh0
        dfv=dfv0; dfh=dfh0
      else if(itur==-1) then !VVD
        open(43,file='vvd.dat',status='old')
        read(43,*) !nvrt
        do j=1,nvrt
          read(43,*)k,dfv0,dfh0
          do i=1,np
            dfv(i,j)=dfv0
            dfh(i,j)=dfh0
          enddo
        enddo !j
        close(43)
      else if(itur==-2) then !HVD
        open(43,file='hvd.mom',status='old')
        open(44,file='hvd.tran',status='old')
        read(43,*)
        read(43,*) !np
        read(44,*)
        read(44,*) !np
        do i=1,np
          read(43,*)k,xtmp,ytmp,dfv0
          read(44,*)k,xtmp,ytmp,dfh0
          do j=1,nvrt
            dfv(i,j)=dfv0
            dfh(i,j)=dfh0
          enddo
        enddo !i=1,np
        close(43)
        close(44)
      else if(itur==2) then !read in P&P coefficients
        read(15,*) h1_pp,vdmax_pp1,vdmin_pp1,tdmin_pp1,h2_pp,vdmax_pp2,vdmin_pp2,tdmin_pp2
        if(h1_pp>=h2_pp) then
          write(11,*)'h1_pp >= h2_pp in P&P'
          stop
        endif
        if(vdmax_pp1<vdmin_pp1.or.vdmax_pp2<vdmin_pp2) then
          write(11,*)'Wrong limits in P&P:',vdmax_pp1,vdmin_pp1,vdmax_pp2,vdmin_pp2
          stop
        endif
      else if(itur==3.or.itur==4) then !read in const. (cf. Umlauf and Burchard 2003)
!       Common variables for both models
        cmiu0=dsqrt(0.3d0)
!       read in mixing limits
        open(31,file='diffmin.gr3',status='old')
        open(32,file='diffmax.gr3',status='old')
        read(31,*)
        read(31,*)
        read(32,*)
        read(32,*)
        do i=1,np
          read(31,*)j,xtmp,ytmp,diffmin(i)
          read(32,*)j,xtmp,ytmp,diffmax(i)
          if(diffmax(i)<diffmin(i)) then
            write(11,*)'diffmin > diffmax:',i
            stop
          endif
        enddo !i
        close(31)
        close(32)

        if(itur==3) then
          read(15,*) mid,stab

!	  Constants used in GLS; cpsi3 later
          a2_cm03=2/cmiu0**3
          eps_min=1.e-12

          select case(mid)
            case('MY') 
              rpub=0; rmub=1; rnub=1; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
              if(stab.ne.'GA') then
                write(11,*)'MY must use Galperins ASM:',stab
                stop
              endif
            case('KL')
              rpub=0; rmub=1; rnub=1; schk=2.44; schpsi=2.44; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
            case('KE')
              rpub=3; rmub=1.5; rnub=-1; schk=1; schpsi=1.3; cpsi1=1.44; cpsi2=1.92
              q2min=1.0e-9; psimin=1.e-8
            case('KW')
              rpub=-1; rmub=0.5; rnub=-1; schk=2; schpsi=2; cpsi1=0.555; cpsi2=0.833
              q2min=1.0e-9; psimin=1.e-8 
            case('UB')
              rpub=2; rmub=1; rnub=-0.67; schk=0.8; schpsi=1.07; cpsi1=1; cpsi2=1.22
              q2min=1.0e-9; psimin=1.e-8 
            case default
              write(11,*)'Unknown closure:',mid
              stop
          end select
          if(rnub==0) then
            write(11,*)'Wrong input for rnub:',rnub
            stop
          endif

          if(stab.ne.'GA'.and.stab.ne.'KC') then
            write(11,*)'Unknown ASM:',stab
            stop
          endif

!	  Consts. used in Canuto's ASM (Model A)
          ubl1=0.1070
          ubl2=0.0032
          ubl3=0.0864
          ubl4=0.12
          ubl5=11.9
          ubl6=0.4
          ubl7=0
          ubl8=0.48
          ubs0=1.5*ubl1*ubl5**2
          ubs1=-ubl4*(ubl6+ubl7)+2*ubl4*ubl5*(ubl1-ubl2/3-ubl3)+1.5*ubl1*ubl5*ubl8
          ubs2=-0.375*ubl1*(ubl6**2-ubl7**2)
          ubs4=2*ubl5
          ubs5=2*ubl4
          ubs6=2*ubl5/3*(3*ubl3**2-ubl2**2)-0.5*ubl5*ubl1*(3*ubl3-ubl2)+0.75*ubl1*(ubl6-ubl7)
          ubd0=3*ubl5**2
          ubd1=ubl5*(7*ubl4+3*ubl8)
          ubd2=ubl5**2*(3*ubl3**2-ubl2**2)-0.75*(ubl6**2-ubl7**2)
          ubd3=ubl4*(4*ubl4+3*ubl8)
          ubd4=ubl4*(ubl2*ubl6-3*ubl3*ubl7-ubl5*(ubl2**2-ubl3**2))+ubl5*ubl8*(3*ubl3**2-ubl2**2)
          ubd5=0.25*(ubl2**2-3*ubl3**2)*(ubl6**2-ubl7**2)
!  	  print*, 'ubd2=',ubd2,',ubd4=',ubd4,',ubd2/ubd4=',ubd2/ubd4

!         Initialize k and l
          do i=1,np
            xlmin2(i)=2*q2min*0.1*dmax1(h0,dp(i)) !floor for non-surface layers
            do k=1,nvrt
              q2(i,k)=q2min
              xl(i,k)=xlmin2(i)
            enddo !k
          enddo !i
          dfv=0; dfh=0; dfq1=0; dfq2=0 !initialize for closure eqs.

        else !itur=4
#ifndef USE_GOTM
          write(11,*)'Compile with GOTM:',itur
          stop
#endif

        endif !itur=3 or 4
      endif !itur

!     i.c.
      read(15,*)icst
      if(icst/=1.and.icst/=2) then
        write(11,*)'Unknown i.c. flag',icst
        stop
      endif

!...  Earth tidal potential
      read(15,*) ntip,tip_dp !cut-off depth for applying tidal potential
      if(ntip>mnbfr) then
        write(11,*)'ntip > mnbfr',ntip,mnbfr
        stop
      endif
      
      if(ntip>0) then
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np
          read(32,*)j,xlon(i),ylat(i)
          xlon(i)=xlon(i)*pi/180
          ylat(i)=ylat(i)*pi/180
!...      Pre-compute species function to save time
          fun_lat(i,0)=3*dsin(ylat(i))**2-1
          fun_lat(i,1)=dsin(2*ylat(i))
          fun_lat(i,2)=dcos(ylat(i))**2
        enddo !i
        close(32)
      endif !ntip>0
      
      do i=1,ntip
        read(15,*) !tag
        read(15,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
        if(jspc(i).lt.0.or.jspc(i).gt.2) then
          write(11,*)'Illegal tidal species #',jspc(i)
          stop
        endif
        tear(i)=tear(i)*pi/180
      enddo !i

!...  Boundary forcing freqs.
      read(15,*) nbfr
      if(nbfr>mnbfr) then
        write(11,*)'nbfr > mnbfr',nbfr,mnbfr
        stop
      endif

      do i=1,nbfr
        read(15,*) !tag
        read(15,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
        face(i)=face(i)*pi/180
      enddo

      read(15,*) nope1
      if(nope1/=nope) then
        write(11,*)'Inconsistent # of open bnds',nope1,nope
        stop
      endif

      nettype=0 !total # of type I bnds
      nfltype=0
      ntetype=0
      nsatype=0
      nettype2=0 !total # of type IV bnds (3D input)
      nfltype2=0 
      ntetype2=0
      nsatype2=0
      do k=1,nope
        read(15,*) ntmp,iettype(k),ifltype(k),itetype(k),isatype(k)
        if(ntmp/=nond(k)) then
          write(11,*)'Inconsistent # of nodes at open boundary',k
          write(11,*)ntmp,nond(k)
          stop
        endif

        if(iettype(k)==1) then
          nettype=nettype+1
!	  Mock reading
          open(50,file='elev.th',status='old')
          do j=1,ntime
            read(50,*) ttt,et
          enddo !j
          rewind(50)
        else if(iettype(k)==2) then
          read(15,*) eth(k,1)
        else if(iettype(k)==3) then
          do i=1,nbfr
            read(15,*)  !freq. name
            do j=1,nond(k)
              read(15,*) emo(k,j,i),efa(k,j,i) !amp. and phase
              efa(k,j,i)=efa(k,j,i)*pi/180
            enddo
          enddo
        else if(iettype(k)==4) then
          nettype2=nettype2+1
          open(54,file='elev3D.th',status='old')
        else if(iettype(k)/=0) then
          write(11,*)'Invalid iettype'
          stop
        endif

        if(ifltype(k)==1) then
          nfltype=nfltype+1
          open(51,file='flux.th',status='old')
          do j=1,ntime
            read(51,*) ttt,qq
          enddo
          rewind(51)
        else if(ifltype(k)==2) then
          read(15,*) qthcon(k)
        else if(ifltype(k)==3) then
          do i=1,nbfr
            read(15,*)
            read(15,*) vmo(k,i),vfa(k,i) !uniform amp. and phase along each segment
            vfa(k,i)=vfa(k,i)*pi/180
          enddo
        else if(ifltype(k)==4) then
          nfltype2=nfltype2+1
          open(55,file='uv3D.th',status='old')
        else if(ifltype(k)==-1) then !Flather 1
          if(iettype(k)/=0) then
            write(11,*)'Flather obc requires iettype=0:',k
            stop
          endif
          read(15,*) eta_m0,qthcon(k)
          do j=1,nond(k)
            eta_mean(iond(k,j))=eta_m0
          enddo !j
        else if(ifltype(k)/=0) then
          write(11,*) 'Invalid ifltype:',ifltype(k)
          stop
        endif

        if(itetype(k)==1) then
          ntetype=ntetype+1
          open(52,file='temp.th',status='old')
          do j=1,ntime
            read(52,*) ttt,temp
          enddo
          rewind(52)
        else if(itetype(k)==2) then
          read(15,*) tth(k,1,1)
        else if(iabs(itetype(k))==4) then
          ntetype2=ntetype2+1
          open(56,file='temp3D.th',status='old')
          if(itetype(k)==-4) read(15,*) tobc(k) !nudging factor
        else if(itetype(k)==-1) then
          read(15,*) tobc(k) !nudging factor
          if(tobc(k)<0.or.tobc(k)>1) then
            write(11,*)'Temp. obc nudging factor wrong:',tobc(k),k
            stop
          endif
        else if(itetype(k)/=0.and.itetype(k)/=3) then
          write(11,*) 'INVALID VALUE FOR ITETYPE'
          stop
        endif

        if(isatype(k)==1) then
          nsatype=nsatype+1
          open(53,file='salt.th',status='old')
          do j=1,ntime
            read(53,*) ttt,sal
          enddo
          rewind(53) 
        else if(isatype(k)==2) then
          read(15,*) sth(k,1,1)
        else if(iabs(isatype(k))==4) then
          nsatype2=nsatype2+1
          open(57,file='salt3D.th',status='old')
          if(isatype(k)==-4) read(15,*) sobc(k) !nudging factor
        else if(isatype(k)==-1) then
          read(15,*) sobc(k) !nudging factor
          if(sobc(k)<0.or.sobc(k)>1) then
            write(11,*)'Salt. obc nudging factor wrong:',sobc(k),k
            stop
          endif
        else if(isatype(k)/=0.and.isatype(k)/=3) then
          write(11,*) 'INVALID VALUE FOR ISATYPE'
          stop
        endif
      enddo !k=1,nope

!     Global output parameters
      noutput=25+ntracers
      if(noutput>mnout) then
        write(11,*)'Increase mnout in the header to',noutput
        stop
      endif
      outfile(1)='elev.61'
      outfile(2)='pres.61'
      outfile(3)='airt.61'
      outfile(4)='shum.61'
      outfile(5)='srad.61'
      outfile(6)='flsu.61'
      outfile(7)='fllu.61'
      outfile(8)='radu.61'
      outfile(9)='radd.61'
      outfile(10)='flux.61'
      outfile(11)='evap.61'
      outfile(12)='prcp.61'
      outfile(13)='wind.62'
      outfile(14)='wist.62'
      outfile(15)='dahv.62'
      outfile(16)='vert.63'
      outfile(17)='temp.63'
      outfile(18)='salt.63'
      outfile(19)='conc.63'
      outfile(20)='tdff.63'
      outfile(21)='vdff.63'
      outfile(22)='kine.63'
      outfile(23)='mixl.63'
      outfile(24)='zcor.63'
      outfile(25)='hvel.64'
      variable_nm(1)='surface elevation'
      variable_nm(2)='atmopheric pressure'
      variable_nm(3)='air temperature'
      variable_nm(4)='specific humidity'
      variable_nm(5)='solar radiation'
      variable_nm(6)='fluxsu'
      variable_nm(7)='fluxlu'
      variable_nm(8)='hradu'
      variable_nm(9)='hradd'
      variable_nm(10)='total flux'
      variable_nm(11)='Evaporation rate (kg/m^2/s)'
      variable_nm(12)='Precipitation rate (kg/m^2/s)'
      variable_nm(13)='wind speed'
      variable_nm(14)='wind stress (m^2/s^2)'
      variable_nm(15)='Depth averaged horizontal velocity'
      variable_nm(16)='vertical velocity'
      variable_nm(17)='temperature in C'
      variable_nm(18)='salinity in psu'
      variable_nm(19)='density in kg/m^3'
      variable_nm(20)='eddy diffusivity in m^2/s'
      variable_nm(21)='eddy viscosity in m^2/s'
      variable_nm(22)='turbulent kinetic energy'
      variable_nm(23)='turbulent mixing length'
      variable_nm(24)='z coordinates'
      variable_nm(25)='horizontal velocity'

      variable_dim(1:12)='2D scalar'
      variable_dim(13:15)='2D vector'
      variable_dim(16:24)='3D scalar'
      variable_dim(25)='3D vector'
    
      do i=1,ntracers
        write(ifile_char,'(i03)')i
        outfile(25+i)='tracer_'//trim(ifile_char)//'.63' 
        variable_nm(25+i)='Tracer #'//trim(ifile_char)
        variable_dim(25+i)='3D scalar'
      enddo !i

      read(15,*) nspool,ihfskip !output and file spools
      if(nspool==0.or.ihfskip==0) then
        write(11,*)'Zero nspool'
        stop
      endif
      if(mod(ihfskip,nspool)/=0) then
        write(11,*)'ihfskip/nspool != integer'
        stop
      endif
      nrec=min(ntime,ihfskip)/nspool

      do i=1,noutput
        read(15,*) iof(i)
        if(iof(i)/=0.and.iof(i)/=1) then
          write(11,*)'Unknown output option',i,iof(i)
          stop
        endif
      enddo !i=1,noutput
      if(iof(24)==0) then
        write(16,*)'Reset zcor output flag'
        iof(24)=1
      endif

!...  Test output parameters
      read(15,*) noutgm 
      if(noutgm/=1.and.noutgm/=0) then
        write(11,*)'Unknown noutgm',noutgm
        stop
      endif
      
!...  input information about hot start output
!...
      read(15,*) nhot
      if(nhot/=0.and.nhot/=1) then
        write(11,*)'Unknown nhot',nhot
        stop
      endif

!...  Itpack solver info
      read(15,*) isolver,itmax1,iremove,zeta,tol
      if(itmax1>itmax) then
        write(11,*)'Increase itmax in header file'
        stop
      endif
      if(isolver<1.or.isolver>4) then
        write(11,*)'Unknown solver',isolver
        stop
      endif

!...  Compute flux flag
      read(15,*) iflux !,ihcheck

!...  Interpolation flag for S,T and vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. will use Kriging as the filter is not applied there
      read(15,*) lq,inter_mom
      if(lq<0.or.lq>2.or.inter_mom<-1.or.inter_mom>1) then
        write(11,*)'Unknown interpolation flag:',lq,inter_mom
        stop
      endif
      if(lq/=0) then
        lqk(1:mne)=lq
      else !lq==0
        open(32,file='lqk.gr3',status='old')
        read(32,*)
        read(32,*)
        do i=1,np
          read(32,*)j,xtmp,ytmp,swild(i)
          if(swild(i)<1.or.swild(i)>2) then
            write(11,*)'Unknown interpolation flag in lqk.gr3'
            stop
          endif
        enddo !i
        close(32)
        do i=1,ne
          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
          lqk(i)=min(swild(n1),swild(n2),swild(n3))
        enddo !i
      endif
      
      if(inter_mom/=-1) then
        krvel=inter_mom
      else !-1
        open(32,file='krvel.gr3',status='old')
        read(32,*)
        read(32,*)
        do i=1,np
          read(32,*)j,xtmp,ytmp,swild(i)
          if(swild(i)<0.or.swild(i)>1) then
            write(11,*)'Unknown interpolation flag in krvel.gr3'
            stop
          endif
        enddo !i
        close(32)
        do i=1,ne
          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
          krvel(i)=min(swild(n1),swild(n2),swild(n3))
        enddo !i
      endif

!...  Interpolation mode (1: along Z; 2: along S)
      open(32,file='interpol.gr3',status='old')
      read(32,*)
      read(32,*)
      do i=1,np
        read(32,*)j,xtmp,ytmp,swild(i)
        if(swild(i)/=1.and.swild(i)/=2) then
          write(11,*)'Unknown interpolation flag in interpol.gr3'
          stop
        endif
      enddo !i
      close(32)
      do i=1,ne
        n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
        interpol(i)=min(swild(n1),swild(n2),swild(n3))
      enddo !i

!     Make sure lqk=2 & interpol=2 are in pure S region
      do i=1,ne
        if(lqk(i)==2.or.interpol(i)==2) then
          if(dp(nm(i,1))>h_s.or.dp(nm(i,2))>h_s.or.dp(nm(i,3))>h_s) then
            write(11,*)'lqk or interpol=2 must be inside pure S region:',i,interpol(i)
            stop
          endif
        endif
      enddo !i

!...  Cut-off depth for BCC
      read(15,*) h_bcc1 !z- or sigma-

!...  Land b.c. option
      read(15,*) islip !0: free slip; otherwise no slip
      if(islip/=0.and.islip/=1) then
        write(11,*)'Unknow islip:',islip
        stop
      endif
      if(islip==1) read(15,*) hdrag0

!...  Nudging options
      read(15,*) inu_st,step_nu,vnh1,vnf1,vnh2,vnf2
      if(inu_st<0.or.inu_st>2.or.step_nu<dt) then
        write(11,*)'Check nudging inputs:',inu_st,step_nu
        stop
      endif
      if(inu_st/=0) then
        if(vnh1>=vnh2.or.vnf1<0.or.vnf1>1.or.vnf2<0.or.vnf2>1) then
          write(11,*)'Check vertical nudging limits:',vnh1,vnf1,vnh2,vnf2
          stop
        endif

        open(96,file='t_nudge.gr3',status='old')
        open(97,file='s_nudge.gr3',status='old')
        read(96,*)
        read(96,*) !ne,np
        read(97,*)
        read(97,*) !ne,np
        do i=1,np
          read(96,*)j,xtmp,ytmp,t_nudge(i)
          read(97,*)j,xtmp,ytmp,s_nudge(i)
          if(t_nudge(i)<0.or.t_nudge(i)>1.or.s_nudge(i)<0.or.s_nudge(i)>1) then
            write(11,*)'Wrong nudging factor at node:',i,t_nudge(i),s_nudge(i)
            stop
          endif
        enddo !i
        close(96)
        close(97)
        
        if(inu_st==2) then
          nrec_nu=nbyte*(1+np*nvrt) !single precision
          open(37,file='temp_nu.in',access='direct',recl=nrec_nu)
          open(35,file='salt_nu.in',access='direct',recl=nrec_nu)
        endif
      endif

!...  Surface min. mixing length for f.s. and max. for all; inactive 
      read(15,*) !xlmax00

      if(itur==3) then
        open(32,file='xlsc.gr3',status='old')
        read(32,*)
        read(32,*)
        do i=1,np
          read(32,*)j,xtmp,ytmp,xlsc0(i)
          if(xlsc0(i)<0.or.xlsc0(i)>1) then
            write(11,*)'Wroing xlsc0:',i,xlsc0(i)
           stop
          endif
        enddo !i
        close(32)
      endif

!...  Order of integration
      read(15,*) mmm
      if(mmm<0) then
        write(11,*)'mmm<0'
        stop
      endif
!     Pre-compute sigmap & sigma_prod for rint_lag()
      if(mmm>0) then
        do k=1,nsig
          do j=1,2*mmm+1
            if(j==1) then
              sigmap(k,j)=sigma(k)
            else
              sigmap(k,j)=sigmap(k,j-1)*sigma(k)
            endif
          enddo !j 

          if(k<nsig) then
            do l=1,k
              j1=max0(l,k-mmm)
              j2=min0(nsig,k+mmm)
              if(j1>=j2) then
                write(11,*)'Weird indices:',j1,j2,k,l
                stop
              endif

              do i=j1,j2
                if(abs(i-k)>4) then
                  write(11,*)'sigma_prod index out of bound'
                  stop
                endif
                sigma_prod(l,k,i-k)=1
                do j=j1,j2
                  if(j/=i) sigma_prod(l,k,i-k)=sigma_prod(l,k,i-k)*(sigma(i)-sigma(j))
                enddo !j
                if(sigma_prod(l,k,i-k)==0) then
                  write(11,*)'Impossible in sigma_prod'
                  stop
                endif
              enddo !i
            enddo !l
          endif !k<nsig
        enddo !k=1,nsig
      endif !mmm>0

!...  Drag formulation
      read(15,*) idrag
      if(idrag/=1.and.idrag/=2) then
        write(11,*)'Unknown idrag'
        stop
      endif
      if(idrag==1.and.itur>0) then
        write(11,*)'Linear drag requires itur<=0'
        stop
      endif
      if(idrag==1.and.nchi/=0) then
        write(11,*)'Linear drag requires nchi=0'
        stop
      endif

!...  ELAD correction option for heat exchange (inactive)
      read(15,*) !ielad 
!      if(ielad/=0.and.ielad/=1) then
!        write(11,*)'Unknown ielad:',ielad
!        stop
!      endif

!...  Option to limit \hat{H} to enhance stability for large friction in shallow area
      read(15,*) ihhat
      if(ihhat/=0.and.ihhat/=1) then
        write(11,*)'Unknown ihhat:',ihhat
        stop
      endif

!...  Transport options: ELM or upwind
      read(15,*) iupwind_t,iupwind_s !0: ELM; 1: upwind; 2: TVD
      if(iupwind_t<0.or.iupwind_t>2.or.iupwind_s<0.or.iupwind_s>2) then
        write(11,*)'Unknown iupwind:',iupwind_t,iupwind_s
        stop
      endif
      if(iupwind_t+iupwind_s==3) then
        write(11,*)'TVD cannot be combined with upwind:',iupwind_t,iupwind_s
        stop
      endif

!     tvd_mid: model AA (my own formulation); CC (Casulli's definition of upwind ratio)
      if(iupwind_t==2.or.iupwind_s==2) read(15,*) tvd_mid,flimiter

!.... Blending factor for vel. in btrack (1 for internal sides; 2 for bnd sides or nodes)
      read(15,*) vis_coe1,vis_coe2
      if(vis_coe1<0.or.vis_coe1>1.or.vis_coe2<0.or.vis_coe2>1) then
        write(11,*)'Illegal vis_coe:',vis_coe1,vis_coe2
        stop
      endif

      read(15,*) shapiro
      if(shapiro<0.or.shapiro>0.5) then
        write(11,*)'Illegal shapiro:',shapiro
        stop
      endif

!     Kriging option
!     Desired # of pts used in Kriging (smaller in reality), and choice of generalized covariance fucntion
!     For Gaussian, also a scale for decorrelation (decorrelation length=decorrel0*(local element size)
      read(15,*) nkrige,kr_co,decorrel0
      if(nkrige>mnei_kr) then
        write(11,*)'Increase mnei_kr:',nkrige
        stop
      endif
      if(kr_co<0.or.kr_co>5) then
        write(11,*)'Wrong kr_co:',kr_co
        stop
      endif

      ie_kr=0 !no Kriging; non-zero value points to local index in all Kriging elements
      ne_kr=0 !total # of elements in Kriging zone
      do i=1,ne
        if(krvel(i)==1) then
          ne_kr=ne_kr+1
          ie_kr(i)=ne_kr
        endif
      enddo !i
      if(ne_kr>mne_kr) then
        write(11,*)'Too many elements in Kriging zone:',ne_kr
        stop
      endif

!...  Max. for vel. magnitude
      read(15,*) rmaxvel
      if(rmaxvel<5) then
        write(11,*)'Illegal rmaxvel:',rmaxvel
        stop
      endif

!...  Inundation algorithm flag (1: better algorithm for fine resolution) 
      read(15,*) inunfl 
      if(inunfl/=0.and.inunfl/=1) then
        write(11,*)'Illegal inunfl:',inunfl
        stop
      endif

!...  Option for calculating nodal vel. (0: discontinous with Shapiro filter; 1: averaging
!...  w/o Shapiro filter: more diffusion)
      read(15,*) indvel
      if(indvel/=0.and.indvel/=1) then
        write(11,*)'Illegal indvel:',indvel
        stop
      endif

!...  Tracer transport
      read(15,*) itmp
      if(itmp/=ntracers) then
        write(11,*)'Mismatch in # of tracers:',itmp,ntracers
        stop
      endif

      if(ntracers>0) then
!        read(15,*) itr_ic !1: horizontal i.c.; 2: vertical
!        if(itr_ic/=1.and.itr_ic/=2) then
!          write(11,*)'Unknown tracer i.c. flag',itr_ic
!          stop
!        endif
        read(15,*) itr_met !=1: upwind; 2: TVD
        if(itr_met/=1.and.itr_met/=2) then
          write(11,*)'Unknown tracer method',itr_met
          stop
        endif
        if(itr_met==2) read(15,*) tvd_mid2,flimiter2

!       b.c.
        read(15,*) !nope
        do k=1,nope
          read(15,*) itrtype(k)
          if(itrtype(k)==2) then
            read(15,*) trth(k,1:ntracers)
          else if(itrtype(k)/=0.and.itrtype(k)/=3) then
            write(11,*)'Wrong itrtype:',k,itrtype(k)
            stop
          endif
        enddo !k
      endif !ntracers

!...  Check last parameter read in from fort.15
      write(*,*)'Last parameter in param.in is flimiter=',flimiter
      close(15)
!     End reading fort.15

!...  Compute neighborhood for Kriging
!     Construct itier_nd and itier_sd
      do i=1,ne
        if(ie_kr(i)==0) cycle

        ie=ie_kr(i) !local index
        itier_nd(ie,0)=3 !use 0 to store actual # of pts used in Kriging
        itier_nd(ie,1:3)=nm(i,1:3)
!       icolor[1,2]: for nodes or sides. 0: outside current ball; 1: inside the ball
        icolor1=0; icolor2=0
        icolor1(nm(i,1:3))=1
        nfront=3
        ifront(1:nfront)=nm(i,1:3) !new frontier nodes
        itr=0 !iteration #
        loop14: do
          itr=itr+1
          if(itr>1000000) stop 'Too many iterations in Kriging'

          nfront0=nfront
          ifront2(1:nfront)=ifront(1:nfront)
          nfront=0
          do j=1,nfront0
            nd=ifront2(j)
            do l=1,nnp(nd)
              nd2=inp(nd,l)
              if(icolor1(nd2)==0) then !new frontier node
                nfront=nfront+1
                if(nfront>mnei_kr) exit loop14 !itier_nd(ie,0) has not been updated; abort
                ifront(nfront)=nd2
                icolor1(nd2)=1
              endif
            enddo !l
          enddo !j=1,nfront0

          if(nfront==0) then !all nodes have been included
            exit loop14
          else
            nold=itier_nd(ie,0)
            if(nold+nfront>nkrige) exit loop14 !itier_nd(ie,0) has not been updated; abort
            do j=1,nfront
              itier_nd(ie,nold+j)=ifront(j)
            enddo !j
            itier_nd(ie,0)=itier_nd(ie,0)+nfront
          endif

!         Debug
!          if(i==14216) then
!            write(97,*)i,itr,nfront,(ifront(j),j=1,nfront)
!          endif

        end do loop14

!        if(i==14216) then
!          nwild=0 !color
!          do j=1,itier_nd(ie,0)
!            nd=itier_nd(ie,j)
!            nwild(nd)=1
!          enddo !j
!          write(98,*)itier_nd(ie,0)
!          write(98,*)np
!          do k=1,np
!            write(98,*)k,real(x(k)),real(y(k)),nwild(k)
!          enddo !k
!       endif
 
!       Check redundancy
!       Comment out after debugging
        do j1=1,itier_nd(ie,0)
          nd1=itier_nd(ie,j1)
          do j2=1,itier_nd(ie,0)
            nd2=itier_nd(ie,j2)
            if(j1/=j2.and.nd1==nd2) then
              write(11,*)'Redundant ball in node:',nd1
              stop
            endif
          enddo !j2
        enddo !j1

      enddo !i=1,ne

!...  Invert Kriging matrices
      akrmat_nd=-1.e34 !initialization for debugging
      err_max=0 !max. error in computing the inverse matices
      do k=1,ne
        if(ie_kr(k)==0) cycle

        ie=ie_kr(k) !local index
        decorrel(ie)=decorrel0*radiel(k) !used in Gaussian covar fucntion
        npp=itier_nd(ie,0)
        do i=1,npp
          n1=itier_nd(ie,i)
          do j=1,npp
            n2=itier_nd(ie,j)
            rr=dsqrt((x(n1)-x(n2))**2+(y(n1)-y(n2))**2)
            akr(i,j)=covar(kr_co,decorrel(ie),rr)
          enddo !j
          akr(i,npp+1)=1
          akr(i,npp+2)=x(n1)
          akr(i,npp+3)=y(n1)
        enddo !i=1,npp

        akr(npp+1,1:npp)=1
        akr(npp+2,1:npp)=x(itier_nd(ie,1:npp))
        akr(npp+3,1:npp)=y(itier_nd(ie,1:npp))
        akr((npp+1):(npp+3),(npp+1):(npp+3))=0
!        bkr(1:(npp+3),1)=0 !does not matter

!       Debug
        akrmat_nd(ie,1:(npp+3),1:(npp+3))=akr(1:(npp+3),1:(npp+3))

!        call gaussj(akr,npp+3,mnei_kr+3,bkr,1,1)

!       LAPACK routines for positive definite symmetric matrix below did not work
!       Note: in LAPACK, the matrix dimension is (LDA,*) so the dimensions will match
!        call dpotrf('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotrf:',info
!          stop
!        endif
!        call dpotri('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotri:',info
!          stop
!        endif
!        do i=1,npp+3
!          do j=i+1,npp+3
!            akr(j,i)=akr(i,j)
!          enddo !j
!        enddo !i

!       Pack symmetric matrix
        do j=1,npp+3
          do i=1,j
            akrp(i+j*(j-1)/2)=akr(i,j)
          enddo !i
        enddo !j
        call dsptrf('U',npp+3,akrp,ipiv,info)
        if(info/=0) then
          write(11,*)'Failed dsptrf:',info,k,decorrel(ie)
          write(11,*)(i,(j,akr(i,j),j=1,npp+3),i=1,npp+3)
          stop
        endif
        call dsptri('U',npp+3,akrp,ipiv,work4,info)
        if(info/=0) then
          write(11,*)'Failed dsptri:',info,k
          stop
        endif
!       Unpack
        do j=1,npp+3
          do i=1,j
            akr(i,j)=akrp(i+j*(j-1)/2)
          enddo !i
        enddo !j

        do i=1,npp+3
          do j=i+1,npp+3
            akr(j,i)=akr(i,j)
          enddo !j
        enddo !i

!       Check
        do i=1,npp+3
          do j=1,npp+3
            suma=0
            do l=1,npp+3
              suma=suma+akrmat_nd(ie,i,l)*akr(l,j)
            enddo !l
            if(i==j) suma=suma-1

            if(k==22910) then
              write(96,*)i,j,akrmat_nd(ie,i,j),akr(i,j),suma
            endif

            if(dabs(suma)>1.e-8) write(98,*)k,i,j,suma
            if(dabs(suma)>err_max) err_max=dabs(suma)
          enddo !j
        enddo !i

        akrmat_nd(ie,1:(npp+3),1:(npp+3))=akr(1:(npp+3),1:(npp+3))
      enddo !k=1,ne

      write(16,*)'Max. error in inverting Kriging maxtrice= ',err_max

!...  Read in horizontal viscosity 
      if(ihorcon/=0) then
        open(32,file='horcon.gr3',status='old')
        read(32,*)
        read(32,*) !ntmp,np
        do i=1,np
          read(32,*)j,xtmp,ytmp,swild(i)
        enddo !i
        close(32)
        do i=1,ns
          horcon(i)=(swild(isidenode(i,1))+swild(isidenode(i,2)))/2
          if(horcon(i)<0.or.horcon(i)>1) then
            write(11,*)'horcon out of bound:',horcon(i),i
            stop
          endif
        enddo !i
      endif !ihorcon/=0

      if(nscreen.eq.1) write(*,*)'done reading inputs...'
      write(16,*)'done reading inputs...'

!...  Initialize bottom index (may be updated for wet/dry or bed deformation) for icst=2
!...  These indices will be used for output only
      do i=1,np
        if(dp(i)<=h0) then
          kbp(i)=0 !dry
        else if(dp(i)<=h_s) then
          kbp(i)=kz
        else
          kbp(i)=0 !flag
          do k=1,kz-1
            if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
              kbp(i)=k
              exit
            endif
          enddo !k
          if(kbp(i)==0) then
            write(11,*)'Cannot find a bottom level for node:',i
            stop
          endif
        endif
      enddo !i
      kbp00=kbp

!								   *
!*******************************************************************
!								   *
!	Initialization for cold start alone			   *
!								   *
!*******************************************************************
!								   *

      if(ihot==0) then
!------------------------------------------------------------------
!...  initialize elevations and vel. 
!...  
      eta2=0
      su2=0; sv2=0
      we=0

!     Parabolic bowl problem 
!      bA=(1.25**2-1)/(1.25**2+1)
!      do i=1,np
!        eta2(i)=dsqrt(1-bA**2)/(1-bA)-1-(x(i)**2+y(i)**2)/1500/1500*((1-bA**2)/(1-bA)**2-1)
!      enddo

!...  read the initial salinity and temperature values from 
!...  salt.ic and temp.ic files. Initial S,T fields may vary
!...  either horizontally (and vertically homogeneous) or vertically 
!...  (horizontally homogeneous). For more general 3D case, use hot start.
!...
      tem0=10; sal0=0 !initialize for extreme cases of wet/dry
      if(ibc==1.and.ibtp==0) then
!	Reset icst
        icst=1
        tem0=10; sal0=0
      else !read in S,T
        if(icst==1) then
          open(24,file='temp.ic',status='old')
          open(25,file='salt.ic',status='old')
          read(24,*) 
          read(24,*) !np
          do i=1,np
            read(24,*) num,xtmp,ytmp,te
            if(te<tempmin.or.te>tempmax) then
              write(11,*)'Initial invalid T at',i,te
              stop
            endif
            do k=1,nvrt
              tem0(k,i)=te
            enddo !k
          enddo !i

          read(25,*) 
          read(25,*) !np
          do i=1,np
            read(25,*) num,xtmp,ytmp,sa
            if(sa<saltmin.or.sa>saltmax) then
              write(11,*)'Initial invalid S at',i,sa
              stop
            endif
            do k=1,nvrt
              sal0(k,i)=sa
            enddo !k
          enddo
          close(24)
          close(25)
        else !icst=2 
!         Read in intial mean S,T
          open(24,file='ts.ic',status='old')
          read(24,*)nz_r
          if(nz_r>mnv.or.nz_r<2) then
            write(11,*)'Change nz_r:',nz_r
            stop
          endif
          do k=1,nz_r
            read(24,*)j,z_r(k),tem1(k),sal1(k)
            if(tem1(k)<tempmin.or.tem1(k)>tempmax.or.sal1(k)<saltmin.or.sal1(k)>saltmax) then
              write(11,*)'Initial invalid S,T at',k,tem1(k),sal1(k)
              stop
            endif
            if(k>=2.and.z_r(k)<z_r(k-1)) then
              write(11,*)'Inverted z-level (0):',k
              stop
            endif
          enddo !k
          close(24)

          do i=1,np
            if(kbp(i)==0) cycle
!           Wet nodes
            do k=kbp(i),nvrt
              if(k>=kz) then !S levels
                kin=k-kz+1
                if(hmod(i)>h_c) then
                  zz=h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
                else
                  zz=hmod(i)*sigma(kin)
                endif
              else if(k==kbp(i)) then !z-levels
                zz=-dp(i)
              else
                zz=ztot(k)
              endif

              if(zz<=z_r(1)) then
                zrat=0; l0=1
              else if(zz>=z_r(nz_r)) then
                zrat=1; l0=nz_r-1
              else
                l0=0 !flag
                do l=1,nz_r-1
                  if(zz>z_r(l).and.zz<=z_r(l+1)) then
                    l0=l
                    exit
                  endif
                enddo !l
                if(l0==0) then
                  write(11,*)'Cannot find a level for S,T:',i,k,zz
                  stop
                endif
                zrat=(zz-z_r(l0))/(z_r(l0+1)-z_r(l0))
              endif
              tem0(k,i)=tem1(l0)+(tem1(l0+1)-tem1(l0))*zrat
              sal0(k,i)=sal1(l0)+(sal1(l0+1)-sal1(l0))*zrat
            enddo !k
!           Extend
            do k=1,kbp(i)-1
              tem0(k,i)=tem0(kbp(i),i)
              sal0(k,i)=sal0(kbp(i),i)
            enddo !k
          enddo !i
        endif
      endif !ibc.eq.1.and.ibtp.eq.0

!...  initialize S,T
      tnd=tem0; snd=sal0

      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        do k=1,nvrt
          tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
          ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
        enddo !k
      enddo !i

      do i=1,ne
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        do k=2,nvrt
          tsel(k,i,1)=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
          tsel(k,i,2)=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
        enddo !k
        tsel(1,i,1)=tsel(2,i,1) !mainly for hotstart format
        tsel(1,i,2)=tsel(2,i,2)
      enddo !i

!...  Tracers; user-defined tracer part
      if(ntracers>0) then
        trel0(1:nvrt,1:ne,1)=1 
        trel0(1:nvrt,1:ne,2)=0 
        trel=trel0
      endif 
!     end user-defined tracer part
      

!...  initialize wind for nws=1,2 (first two lines)
!...
      if(nws==1) then
        open(22,file='wind.th',status='old')
        read(22,*) wx1,wy1
        read(22,*) wx2,wy2
        do i=1,np
          windx1(i)=wx1
          windy1(i)=wy1
          windx2(i)=wx2
          windy2(i)=wy2
        enddo
        wtime1=0
        wtime2=wtiminc 
      endif

!	CORIE mode
      if(nws>=2) then
        wtime1=0
        wtime2=wtiminc 
#ifdef USE_SFLUX
        call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
        call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
#endif

!       Uncomment the following to overwrite wind (similary airt etc)
!       Search for "Overwrite wind with wind.th"; WARNING: wind.th has time step of dt not wtiminc, and 
!       starts from t=0,dt,2*dt ... (format: windu,windv; similar to nws=1).
!       Overwrite wind with wind.th
!        open(22,file='wind.th',status='old')
!        read(22,*) wx1,wy1
!        windx1=wx1; windx2=wx1
!        windy1=wy1; windy2=wy1
!       End

        if(nws==3) then
!         Open sflux.th (time interval is dt, not wtiminc!)
          open(23,file='sflux.th',status='old')
          read(23,*) !t=0
!         To heat up water, fluxsu00<0, srad00>0
          read(23,*) tmp,fluxsu00,srad00 !time, total surface flux, solar radiation
        endif !nws==3
      endif !nws>=2

!...  Read initial nudging S,T
      if(inu_st==2) then
        read(37,rec=1)floatout,((tnd_nu1(i,j),j=1,nvrt),i=1,np)
        read(35,rec=1)floatout,((snd_nu1(i,j),j=1,nvrt),i=1,np)
        read(37,rec=2)floatout,((tnd_nu2(i,j),j=1,nvrt),i=1,np)
        read(35,rec=2)floatout,((snd_nu2(i,j),j=1,nvrt),i=1,np)
        irec_nu=2
        time_nu=step_nu
      endif

!------------------------------------------------------------------
      endif !ihot=0

!...  Initialize GOTM for both cold and hot starts (for cde etc).
!...  For real hot start, q2, xl, dfv and dfh will use the values in hotstart.in;
!...  otherwise they will be assigned values below.
      if(itur==4) then
#ifdef USE_GOTM
          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
          call init_tridiagonal(nvrt-1)
#endif
      endif

      if(nscreen.eq.1) write(*,*)'done initializing cold start'
      write(16,*)'done initializing cold start'
      

!                                                                             
!******************************************************************************
!                                                                             *
!		hot start setup of the program				      *
!                                                                             *
!******************************************************************************
!
!     Record length for hot start files (double precision for all reals)
      ihot_len=nbyte*(3+(6*nvrt+1+4*nvrt*ntracers)*ne+(8*nvrt+1)*ns+3*np+20*np*nvrt+1)+12

      if(ihot/=0) then
        open(36,file='hotstart.in',access='direct',recl=ihot_len)
        read(36,rec=1)iths,time,(idry_e(i),(we(j,i),tsel(j,i,1),tsel(j,i,2), &
     &(trel0(j,i,l),trel(j,i,l),l=1,ntracers),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
     &dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),j=1,nvrt),i=1,np),ifile,ifile_char


!       Output hotstart.out for debugging
        open(32,file='hotstart.out')
        write(32,*)iths,time,ifile
        write(32,'(a12)')ifile_char
        write(32,*)'element,idry_e,level,we'
        do i=1,ne,100
          write(32,*)i,idry_e(i),(j,we(j,i),j=1,nvrt)
        enddo !i
        write(32,*)'side,idry_s,level,su2,sv2,tsd,ssd'
        do i=1,ns,100
          write(32,*)i,idry_s(i),(j,real(su2(j,i)),real(sv2(j,i)),real(tsd(j,i)),real(ssd(j,i)),j=1,nvrt)
        enddo !i
        write(32,*)'node,eta2,idry,level,tnd,snd,tem0,sal0,q2,xl'
        do i=1,np,100
          write(32,*)i,eta2(i),idry(i), &
     &(j,real(tnd(j,i)),real(snd(j,i)),real(tem0(j,i)),real(sal0(j,i)),real(q2(i,j)),real(xl(i,j)),dfv(i,j),j=1,nvrt)
        enddo !i
        close(32)

        if(itur==3) then
          do i=1,np
            do j=1,nvrt
              q2(i,j)=dmax1(q2min,q2(i,j))
              xl(i,j)=dmax1(xlmin2(i),xl(i,j))
            enddo
          enddo
        endif
        close(36)

!...    change time and iteration for forecast mode
!...    Causion: this affects all t.h. files (fort.5[0-3]) and wind files
        if(ihot==1) then
          time=0
          iths=0
        endif

        write(*,*)'hot start at time=',time,iths
        write(16,*)'hot start at time=',time,iths

!...  find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...
        if(nws==1) then
          open(22,file='wind.th',status='old')
          rewind(22)
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          do it=0,ninv
            read(22,*)wx1,wy1
          enddo
          read(22,*)wx2,wy2
          do i=1,np
            windx1(i)=wx1
            windy1(i)=wy1
            windx2(i)=wx2
            windy2(i)=wy2
          enddo
        endif

        if(nws>=2) then
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
#ifdef USE_SFLUX
          call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
#endif

!         Overwrite wind with wind.th
!          open(22,file='wind.th',status='old')
!          rewind(22)
!          do it=0,iths
!            read(22,*)wx1,wy1
!          enddo !it
!          windx1=wx1; windx2=wx1
!          windy1=wy1; windy2=wy1
!         End

          if(nws==3) then
!           Read sflux.th
            open(23,file='sflux.th',status='old')
            rewind(23)
            do it=0,iths
              read(23,*)
            enddo
            read(23,*)tmp,fluxsu00,srad00
          endif !nws==3
        endif !nws

!...    Nudging 
        if(inu_st==2) then
          irec_nu=time/step_nu+1
          time_nu=irec_nu*step_nu
          read(37,rec=irec_nu)floatout,((tnd_nu1(i,j),j=1,nvrt),i=1,np)
          read(35,rec=irec_nu)floatout,((snd_nu1(i,j),j=1,nvrt),i=1,np)
          read(37,rec=irec_nu+1)floatout,((tnd_nu2(i,j),j=1,nvrt),i=1,np)
          read(35,rec=irec_nu+1)floatout,((snd_nu2(i,j),j=1,nvrt),i=1,np)
          irec_nu=irec_nu+1
        endif

!...    Find positions in t.h. files 
        if(nettype>0) then
          do it=1,iths
            read(50,*) ttt,et
          enddo !it
        endif

        if(nfltype>0) then
          do it=1,iths
            read(51,*) ttt,qq
          enddo !it
        endif

        if(ntetype>0) then
          do it=1,iths
            read(52,*) ttt,te
          enddo !it
        endif

        if(nsatype>0) then
          do it=1,iths
            read(53,*) ttt,sal
          enddo !it
        endif

        if(nettype2>0) then
          do it=1,iths
            read(54,*) ttt
            do i=1,nope
              if(iettype(i)==4) then
                 do j=1,nond(i)
                   read(54,*)nd2,eth(i,j)
                 enddo !j
              endif
            enddo !i
          enddo !it
        endif

        if(nfltype2>0) then
          do it=1,iths
            read(55,*) ttt
            do i=1,nope
              if(ifltype(i)==4) then
                 do j=1,nond(i)
                   read(55,*)nd2,(uthnd(i,j,k),vthnd(i,j,k),k=1,nvrt)
                 enddo !j
              endif
            enddo !i
          enddo !it
        endif

        if(ntetype2>0) then
          do it=1,iths
            read(56,*) ttt
            do i=1,nope
              if(iabs(itetype(i))==4) then
                 do j=1,nond(i)
                   read(56,*)nd2,(tth(i,j,k),k=1,nvrt)
                 enddo !j
              endif
            enddo !i
          enddo !it
        endif

        if(nsatype2>0) then
          do it=1,iths
            read(57,*) ttt
            do i=1,nope
              if(iabs(isatype(i))==4) then
                 do j=1,nond(i)
                   read(57,*)nd2,(sth(i,j,k),k=1,nvrt)
                 enddo !j
              endif
            enddo !i
          enddo !it
        endif

!...  end hot start section
!...
      endif !ihot.ne.0

      if(nscreen.eq.1) write(*,*)'done initializing variables...'
      write(16,*)'done initializing variables...'

!                                                                             *
!******************************************************************************
!                                                                             *
!			open output files				      *
!                                                                             *
!******************************************************************************
!                                                                             *

!...  write global output headers
!...

      if(ihot<=1) then
        ifile=1 !output file #
!       Convert it to a string
        write(ifile_char,'(i12)') ifile
      endif
      call header

      if(nscreen.eq.1) write(*,*)'done initializing outputs'
      write(16,*)'done initializing outputs'

!                                                                             *
!******************************************************************************
!                                                                             *
!	    	          Time stepping 				      *
!                                                                             *
!******************************************************************************
!                                                                             *

      if(ihot==0) iths=0
      if(nscreen.eq.1) write(*,*)'time stepping begins...',iths+1,ntime
      write(16,*)'time stepping begins...',iths+1,ntime

!...  Compute initial bed deformation and update depths info
      do i=1,np
        bdef1(i)=bdef(i)/ibdef*min0(iths,ibdef)
        dp(i)=dp00(i)-bdef1(i)
        hmod(i)=dmin1(dp(i),h_s)
      enddo !i
      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
      do i=1,ne
        dpe(i)=1.e10
        do j=1,3
          if(dpe(i)>dp(nm(i,j))) dpe(i)=dp(nm(i,j))
        enddo !j
      enddo !i=1,ne

!...  Compute initial vgrid
      if(inunfl==0) then
        call levels0(iths,iths)
      else
        call levels1(iths,iths)
      endif
      if(nscreen.eq.1) write(*,*)'done computing initial vgrid...'
      write(16,*)'done computing initial vgrid...'

!...  Compute nodal vel. 
      call nodalvel(ifltype)
      if(nscreen.eq.1) write(*,*)'done computing initial nodal vel...'
      write(16,*)'done computing initial nodal vel...'

!...  Debug: test btrack alone
!      eta1=0; eta2=0; we=0
!      do i=1,ns
!        do k=1,nvrt
!          su2(k,i)=-ycj(i)*2*pi/3.0e3 
!          sv2(k,i)=xcj(i)*2*pi/3.0e3
!        enddo !k
!      enddo !i
!      do i=1,ne
!        do k=1,nvrt
!          do j=1,3
!            nd=nm(i,j)
!            ufg(k,i,j)=-y(nd)*2*pi/3.0e3
!            vfg(k,i,j)=x(nd)*2*pi/3.0e3
!          enddo !j
!        enddo !k
!      enddo !i
!      do i=1,np
!        do k=1,nvrt
!          uu2(k,i)=-y(i)*2*pi/3.0e3 !5.e-16/81*x(i)**4
!          vv2(k,i)=x(i)*2*pi/3.0e3
!          ww2(k,i)=0 !-1.e-4*z(k,i)*(50+z(k,i))
!        enddo !k
!      enddo !i
!...  End test btrack

!...  Compute initial density
      call eqstate
      if(nscreen.eq.1) write(*,*)'done computing initial density...'
      write(16,*)'done computing initial density...'

!...  Initialize heat budget model
      if(nws>=2.and.ihconsv/=0) then
#ifdef USE_SFLUX
        call surf_fluxes(wtime1,windx1,windy1,pr1,airt1,shum1,srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                   fluxprc,fluxevp, &
#endif
     &                   nws,fluxsu00,srad00)
#endif
        do i=1,np
          sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
        enddo
        if(nscreen.eq.1) write(*,*)'heat budge model completes...'
        write(16,*)'heat budge model completes...'
      endif

!...  Assign variables in GOTM for cold starts
      if(itur==4.and.(ihot==0.or.ihot==1.and.nramp==1)) then
#ifdef USE_GOTM
!          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
!          call init_tridiagonal(nvrt-1)
!         Debug
!          do k=0,nvrt-1
!            write(99,*)k,tke1d(k),L1d(k),num1d(k),nuh1d(k)
!          enddo !i
!          stop

          do j=1,np
            q2(j,1:nvrt) = tke1d(0:(nvrt-1))
            xl(j,1:nvrt) = L1d(0:(nvrt-1))
            dfv(j,1:nvrt) = dmin1(diffmax(j),dmax1(diffmin(j),num1d(0:(nvrt-1))))
            dfh(j,1:nvrt) = dmin1(diffmax(j),dmax1(diffmin(j),nuh1d(0:(nvrt-1))))
          enddo !j
#endif
      endif !itur==4 etc
!...
!...  Begin time stepping
!...
      do it=iths+1,ntime

      time=it*dt 

!...  define ramp function for boundary elevation forcing, wind and pressure
!...  forcing and tidal potential forcing
!...
      if(ibc==0) then
        if(nrampbc/=0) then
          rampbc=tanh(2*time/86400/drampbc)
        else
          rampbc=1
        endif
      endif

      if(nws>0.and.nrampwind/=0) then
        rampwind=tanh(2*time/86400/drampwind)
      else
        rampwind=1
      endif

      if(nramp==1) then
        ramp=tanh(2*time/86400/dramp)
      else
        ramp=1
      endif

!...  Compute new bed deformation
      do i=1,np
        bdef2(i)=bdef(i)/ibdef*min0(it,ibdef)
      enddo !i

!...  Bottom drag coefficients for nchi=1; Cd and Cdp for nchi=0 already read in
      if(nchi==1) then !idrag=2 
        Cdp=0; Cd=0 !for dry pts
        do i=1,ns
          if(idry_s(i)==1) cycle

!         Wet side
          htot=dps(i)+(eta2(isidenode(i,1))+eta2(isidenode(i,2)))/2
          if(rough(i)<=0) then !time-independent Cd
            Cd(i)=dabs(rough(i))
          else !roughness >0
            bthick=zs(kbs(i)+1,i)-zs(kbs(i),i) !thickness of bottom bnd layer
            if(bthick<=rough(i)) then
              if(ifort12(17)==0) then
                ifort12(17)=1
                write(12,*)'BL too fine:',i,bthick,rough(i),htot
              endif
              Cd(i)=Cdmax
            else
              Cd(i)=1/(2.5*dlog(bthick/rough(i)))**2 
              Cd(i)=dmin1(Cd(i),Cdmax)
            endif
          endif
        enddo !i=1,ns

!       Drag at nodes
        do i=1,np
          if(idry(i)==1) cycle
!         Wet node
          htot=dp(i)+eta2(i)
          if(rough_p(i)<=0) then !time-independent Cd
            Cdp(i)=dabs(rough_p(i))
          else !roughness >0
            bthick=z(kbp(i)+1,i)-z(kbp(i),i) !thickness of bottom bnd layer
            if(bthick<=rough_p(i)) then
              if(ifort12(5)==0) then
                ifort12(5)=1
                write(12,*)'BL too fine (2):',i,bthick,rough_p(i),htot
              endif
              Cdp(i)=Cdmax
            else
              Cdp(i)=1/(2.5*dlog(bthick/rough_p(i)))**2 
              Cdp(i)=dmin1(Cdp(i),Cdmax)
            endif
          endif
        enddo !i=1,np

!       Output Cd for first step
        if(it==iths+1) then
          open(32,file='Cd.out')
          write(32,*)'Drag coefficents for nchi=1'
          write(32,*)ns
          do i=1,ns
            write(32,'(i6,2e14.6,1x,e9.3)')i,xcj(i),ycj(i),Cd(i)
          enddo !i=1,ns
          close(32)
        endif
      endif !nchi==1

!     Check viscosity in linear drag formulation
!      if(idrag==1.and.it==iths+1) then
!        open(32,file='viscosity_linear.out')
!        write(32,*)'Compare viscosity for idrag=1'
!        write(32,*)ns,dfv0
!        do i=1,ns
!          if(idry_s(i)==0) then
!            bthick=zs(2,i)-zs(1,i)
!            write(32,'(i6,3e14.6)')i,xcj(i),ycj(i),bthick*Cd(i)
!          endif
!        enddo !i=1,ns
!        close(32)
!      endif

!...  Horizontal diffusion; compute d2u, d2v (incorporated hvis inside)
!     Bypss this section if ihorcon=0 (all horcon=0) and then isidenei() 
!     side_ac, and side_x are not used in the code 
      call system_clock(ist,icount_rate)

      d2u=0; d2v=0 !for dry sides
      if(ihorcon/=0) then !not all horcon(i)=0
!       Compute d[u,v]/dx (x being local normal) first (saved in sdbt) (incorporated hvis inside)
        sdbt=0 !for dry side etc

        do i=1,ns
          if(idry_s(i)==1) cycle

!         Wet side
          if(is(i,2)==0) then !bnd sides
            if(isbs(i)<=0.and.islip==1) then !no-slip land bnd 
              do k=kbs(i),nvrt
                ie=is(i,1)
                icount=0
                av_u=0
                av_v=0
                do j=1,3
                  id=js(ie,j)
                  if(is(id,2)==0) cycle
                  icount=icount+1
                  av_u=av_u+su2(k,id) !no vertical interpolation
                  av_v=av_v+sv2(k,id)
                enddo !j
                if(icount==0) then
                  write(11,*)'Isolated element'
                  stop
                endif
                av_u=av_u/icount
                av_v=av_v/icount
                av_mag=dsqrt(av_u**2+av_v**2)
                sdbt(i,k,1)=-hdrag0*av_mag*av_u
                sdbt(i,k,2)=-hdrag0*av_mag*av_v
              enddo !k
            endif
            cycle
          endif

          if(horcon(i)==0) cycle

!         Internal wet sides
          x_1=side_x(i,1)
          x_2=side_x(i,2)
          node1=isidenode(i,1)
          node2=isidenode(i,2)
          in11=lindex(node1,is(i,1))
          in12=lindex(node1,is(i,2))
          in21=lindex(node2,is(i,1))
          in22=lindex(node2,is(i,2))
          if(in11==0.or.in12==0.or.in21==0.or.in22==0) then
            write(11,*)'Wrong sides (9):',in11,in12,in21,in22,node1,node2,is(i,1:2)
            stop
          endif
          do k=kbs(i),nvrt
!           Do vertical interpolation
            swild7=-1.e34 !flag
            do l=1,4 !2 intersections + 2 adjacent elements of side i
              if(l<=2) then
                ie=isidenei(i,l) !parent element for 2 intersections
              else
                ie=is(i,l-2) !adjacent element
                if(ie==isidenei(i,l-2)) cycle
              endif

              do j=1,3 !sides
                id=js(ie,j)
                if(idry_s(id)==1) then
                  swild(1:2)=0
                else
                  kbb=kbs(id)
                  alow(kbb:nvrt)=zs(kbb:nvrt,id)
                  swild2(kbb:nvrt,1)=su2(kbb:nvrt,id)
                  swild2(kbb:nvrt,2)=sv2(kbb:nvrt,id)
                  call vinter(mnv,2,zs(k,i),kbb,nvrt,k,alow,swild2,swild,ibelow)
                endif
                swild6(j,1:2)=swild(1:2)

!               Debug
!                if(it==150) then
!                  write(99,*)i,l,j,swild(1),k,su2(kbb:nvrt,id)
!                  write(99,*)i,l,j,swild(2),k,sv2(kbb:nvrt,id)
!                endif

              enddo !j=1,3
              swild5(1,1:2)=swild6(2,1:2)+swild6(3,1:2)-swild6(1,1:2) !@ node 1
              swild5(2,1:2)=swild6(1,1:2)+swild6(3,1:2)-swild6(2,1:2)
              swild5(3,1:2)=swild6(1,1:2)+swild6(2,1:2)-swild6(3,1:2)
!             Save (u,v) at node[12] for y-derivatives
              if(ie==is(i,1)) then
                swild7(1,1,1:2)=swild5(in11,1:2)
                swild7(2,1,1:2)=swild5(in21,1:2)
              endif
              if(ie==is(i,2)) then
                swild7(1,2,1:2)=swild5(in12,1:2)
                swild7(2,2,1:2)=swild5(in22,1:2)
              endif

              if(l<=2) then
                swild(1:2)=side_ac(i,l,1:2)
                swild(3)=1-swild(1)-swild(2)
                soln(l,1:2)=0 !u[12], v[12]
                do j=1,3
                  soln(l,1)=soln(l,1)+swild(j)*swild5(j,1)
                  soln(l,2)=soln(l,2)+swild(j)*swild5(j,2)
                enddo !j

!               Debug
!                if(it==150) then
!                  write(98,*)i,k,l,soln(l,1),swild5(1:3,1)
!                  write(98,*)i,k,l,soln(l,2),swild5(1:3,2)
!                endif
              endif !l<=2

            enddo !l=1,4
            do i1=1,2
              do i2=1,2
                do i3=1,2
                  if(swild7(i1,i2,i3)<-1.e33) then
                    write(11,*)'swild7 not assigned:',i1,i2,i3
                    stop
                  endif
                enddo !i3
              enddo !i2
            enddo !i1

            u_1=soln(1,1); u_2=soln(2,1)
            v_1=soln(1,2); v_2=soln(2,2)
            dudx=(x_2**2*(u_1-su2(k,i))-x_1**2*(u_2-su2(k,i)))/x_1/x_2/(x_2-x_1)
            dvdx=(x_2**2*(v_1-sv2(k,i))-x_1**2*(v_2-sv2(k,i)))/x_1/x_2/(x_2-x_1)
            dudy=(swild7(2,1,1)+swild7(2,2,1)-swild7(1,1,1)-swild7(1,2,1))/2/distj(i)
            dvdy=(swild7(2,1,2)+swild7(2,2,2)-swild7(1,1,2)-swild7(1,2,2))/2/distj(i)

!           Transform to global coordinates
            dudx_g=dudx*snx(i)-dudy*sny(i)
            dudy_g=dudx*sny(i)+dudy*snx(i)
            dvdx_g=dvdx*snx(i)-dvdy*sny(i)
            dvdy_g=dvdx*sny(i)+dvdy*snx(i)

            hvis(k,i)=horcon(i)*(area(is(i,1))+area(is(i,2)))*dsqrt(dudx_g**2+dvdy_g**2+(dvdx_g+dudy_g)**2/2)
            sdbt(i,k,1)=hvis(k,i)*dudx
            sdbt(i,k,2)=hvis(k,i)*dvdx

!             Debug
!              if(it==500) then
!                write(98,*)node1,node2,k,hvis(k,i)
!              endif

          enddo !k=kbs(i),nvrt
        enddo !i=1,ns

!       Debug
!        if(it==500) stop

        do i=1,ns
          if(idry_s(i)==1) cycle

!         Wet sides
          do k=kbs(i),nvrt
            ta=0
            sumu=0 !integral
            sumv=0 !integral
            do j=1,2 !elements
              ie=is(i,j)
              if(ie==0) cycle

              ta=ta+area(ie)
              do l=1,3 !sides
                id=js(ie,l)
                if(is(i,2)/=0.and.id==i) cycle

!               Do vertical interpolation
                if(idry_s(id)==1) then
                  swild(1:2)=0
                else
                  kbb=kbs(id)
                  alow(kbb:nvrt)=zs(kbb:nvrt,id)
                  swild2(kbb:nvrt,1)=sdbt(id,kbb:nvrt,1)
                  swild2(kbb:nvrt,2)=sdbt(id,kbb:nvrt,2)
                  call vinter(mnv,2,zs(k,i),kbb,nvrt,k,alow,swild2,swild,ibelow)
                endif
                sumu=sumu+swild(1)*ssign(ie,l)*distj(id)
                sumv=sumv+swild(2)*ssign(ie,l)*distj(id)
              enddo !l=1,3
            enddo !j; 2 adjacent elements
         
            if(ta==0) then
              write(11,*)'Impossible 127'
              stop
            endif
            d2u(k,i)=sumu/ta
            d2v(k,i)=sumv/ta
          enddo !k=kbs(i),nvrt
        enddo !i=1,ns
      endif !ihorcon/=0 

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'done hvis and bottom fric: ',btimer,'seconds'
!'
      write(16,*)'done hvis and bottom fric: ',btimer,' seconds..'
      
!...  Earth tidal potential at nodes: pre-compute to save time
!...
      do i=1,np
        etp(i)=0
        do jf=1,ntip
          ncyc=int(tfreq(jf)*time/2/pi)
          arg=tfreq(jf)*time-ncyc*2*pi+jspc(jf)*xlon(i)+tear(jf)
          etp(i)=etp(i)+ramp*tamp(jf)*tnf(jf)*fun_lat(i,jspc(jf))*dcos(arg)
        enddo !jf
      enddo !i

!...  process new wind info 
!...
      if(nws==1) then
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          read(22,*) wx2,wy2
          do i=1,np
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
            windx2(i)=wx2
            windy2(i)=wy2
          enddo
        endif

        wtratio=(time-wtime1)/wtiminc
        do i=1,np
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        enddo !i
      endif !nws=1

!     CORIE mode
      if(nws>=2) then
        if(time>=wtime2) then
!...      Heat budget & wind stresses
          if(ihconsv/=0) then
#ifdef USE_SFLUX
            call surf_fluxes(wtime2,windx2,windy2,pr2,airt2,shum2,srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                       fluxprc,fluxevp, &
#endif
     &                       nws,fluxsu00,srad00)
#endif
            do i=1,np
              sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
            enddo
            if(nscreen.eq.1) write(*,*)'heat budge model completes...'
            write(16,*)'heat budge model completes...'
          endif !ihconsv.ne.0

          wtime1=wtime2
          wtime2=wtime2+wtiminc
          do i=1,np
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
            pr1(i)=pr2(i)
            airt1(i)=airt2(i)
            shum1(i)=shum2(i)
          enddo
#ifdef USE_SFLUX
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
#endif
        endif !time>=wtime2

        wtratio=(time-wtime1)/wtiminc
        do i=1,np
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i

!       Overwrite wind with wind.th
!        read(22,*)wx2,wy2
!        windx1=wx2; windy1=wy2
!        windx2=wx2; windy2=wy2
!        windx=wx2; windy=wy2
!       End

!       Read in new flux values for next step
        if(nws==3) read(23,*) tmp,fluxsu00,srad00
      endif !nws>=2

!...  Re-scale wind
      if(nws>0) then
        do i=1,np
          windx(i)=windx(i)*windfactor(i)
          windy(i)=windy(i)*windfactor(i)
        enddo !i
      endif

!...  compute wind stress components
      dragcmin=1.0d-3*(0.61+0.063*6)
      dragcmax=1.0d-3*(0.61+0.063*50)
      do i=1,np
        if(nws==0) then
          tau(i,1)=0
          tau(i,2)=0
        else if(nws==1.or.nws>=2.and.ihconsv==0) then
          wmag=dsqrt(windx(i)**2+windy(i)**2)
          dragcoef=1.0d-3*(0.61+0.063*wmag)
          dragcoef=dmin1(dmax1(dragcoef,dragcmin),dragcmax)
          tau(i,1)=dragcoef*0.001293*wmag*windx(i)*rampwind
          tau(i,2)=dragcoef*0.001293*wmag*windy(i)*rampwind
        else !nws>=2 and ihconsv !=0; tauxz and tauyz defined
          if(idry(i)==1) then
            tau(i,1)=0
            tau(i,2)=0
          else !rescale as well
            tau(i,1)=-tauxz(i)/rho0*rampwind*windfactor(i)**2 !sign and scale difference between stresses tauxz and tau
            tau(i,2)=-tauyz(i)/rho0*rampwind*windfactor(i)**2
          endif
        endif !nws
      enddo !i=1,np

      if(nscreen.eq.1) write(*,*)'done adjusting wind stress ...'
      write(16,*)'done adjusting wind stress ...'

!...  Read in temp. and salt for nudging
      if(inu_st==2) then
        if(time>time_nu) then
          irec_nu=irec_nu+1
          time_nu=time_nu+step_nu
          tnd_nu1=tnd_nu2
          snd_nu1=snd_nu2
          read(37,rec=irec_nu)floatout,((tnd_nu2(i,j),j=1,nvrt),i=1,np)
          read(35,rec=irec_nu)floatout,((snd_nu2(i,j),j=1,nvrt),i=1,np)
          if(floatout/=time_nu) then
            write(11,*)'Wrong nudging time:',floatout,time_nu
            stop
          endif
        endif !time>time_nu

!       Compute S,T
        rat=(time_nu-time)/step_nu
        if(rat<0.or.rat>1) then
          write(11,*)'Impossible 81:',rat
          stop
        endif
        tnd_nu=tnd_nu1+(1-rat)*(tnd_nu2-tnd_nu1)
        snd_nu=snd_nu1+(1-rat)*(snd_nu2-snd_nu1)
      endif !nudging

!...  Get new t.h. values *.th
!...
      if(nettype>0) then
        read(50,*) ttt,(ath(i),i=1,nettype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for eta',it,ttt
          stop
        endif
      
        icount=0
        do k=1,nope
          if(iettype(k)==1) then
            icount=icount+1
            if(icount>nettype) then
              write(11,*)'Wrong counting 1'
              stop
            endif
            eth(k,1)=ath(icount)
          endif
        enddo 
      endif

      if(nfltype>0) then
        read(51,*) ttt,(ath(i),i=1,nfltype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for flux',it,ttt,time
          stop
        endif

        icount=0
        do k=1,nope
          if(ifltype(k)==1) then
            icount=icount+1
            if(icount>nfltype) then
              write(11,*)'Wrong counting 2'
              stop
            endif
            qthcon(k)=ath(icount)
          endif
        enddo !k
      endif

      if(ntetype>0) then
        read(52,*) ttt,(ath(i),i=1,ntetype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for temp',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(itetype(k)==1) then
            icount=icount+1
            if(icount>ntetype) then
              write(11,*)'Wrong counting 3'
              stop
            endif
            tth(k,1,1)=ath(icount)
          endif
        enddo !k
      endif

      if(nsatype>0) then
        read(53,*) ttt,(ath(i),i=1,nsatype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for salt',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(isatype(k)==1) then
            icount=icount+1
            if(icount>nsatype) then
              write(11,*)'Wrong counting 4'
              stop
            endif
            sth(k,1,1)=ath(icount)
          endif
        enddo !k
      endif

      if(nettype2>0) then
        read(54,*) ttt
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for eta 2',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(iettype(k)==4) then
            icount=icount+1
            if(icount>nettype2) then
              write(11,*)'Wrong counting 7'
              stop
            endif
            do j=1,nond(k)
              nd=iond(k,j)
              read(54,*)nd2,eth(k,j)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in elev3D.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

      if(nfltype2>0) then
        read(55,*) ttt
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for flux 2',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(ifltype(k)==4) then
            icount=icount+1
            if(icount>nfltype2) then
              write(11,*)'Wrong counting 6'
              stop
            endif
            do j=1,nond(k)
              nd=iond(k,j)
              read(55,*)nd2,(uthnd(k,j,l),vthnd(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in uv.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

      if(ntetype2>0) then
        read(56,*) ttt
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for temp. 2',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(iabs(itetype(k))==4) then
            icount=icount+1
            if(icount>ntetype2) then
              write(11,*)'Wrong counting 8'
              stop
            endif
            do j=1,nond(k)
              nd=iond(k,j)
              read(56,*)nd2,(tth(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in temp3D.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

      if(nsatype2>0) then
        read(57,*) ttt
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(11,*)'Starting time wrong for salt 2',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(iabs(isatype(k))==4) then
            icount=icount+1
            if(icount>nsatype2) then
              write(11,*)'Wrong counting 9'
              stop
            endif
            do j=1,nond(k)
              nd=iond(k,j)
              read(57,*)nd2,(sth(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in salt3D.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

!...  Compute new vel. for flow b.c.
!     Average total depth for calcualtion of cross-section areas
      do k=1,nope
        if(ifltype(k)/=0) then
          atd(k)=0
          do i=1,nond(k)
            nd=iond(k,i)
            H2=dp(nd)+eta2(nd)
            if(H2<=h0) then
              write(11,*)'Dry bnd side:',H2,k,i
              stop
            endif
            atd(k)=atd(k)+H2/nond(k)
          enddo !i
        endif
      enddo !k

      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        ibnd=isbs(i)
        if(ibnd<=0) cycle

!       Open bnds
        if(iabs(ifltype(ibnd))==1.or.ifltype(ibnd)==2) then !including Flather 1
          if(atd(ibnd)<h0) then
            write(11,*)'Dry bnd side 2',ibnd,atd(ibnd)
            stop
          endif
          vnth0=qthcon(ibnd)*ramp/cwidth(ibnd)/atd(ibnd)
          do k=1,nvrt
            uth(i,k)=vnth0*snx(i)
            vth(i,k)=vnth0*sny(i)
          enddo !k
        else if(ifltype(ibnd)==3) then
          vnth0=0 !normal vel.
          do jfr=1,nbfr
            ncyc=int(amig(jfr)*time/2/pi)
            arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-vfa(ibnd,jfr)
            vnth0=vnth0+ramp*ff(jfr)*vmo(ibnd,jfr)*dcos(arg)
          enddo !jfr=1,nbfr
          do k=1,nvrt
            uth(i,k)=vnth0*snx(i)
            vth(i,k)=vnth0*sny(i)
          enddo !k
        else if(ifltype(ibnd)==4) then
!         Find bnd node indices for n1,n2
          j1=0; j2=0
          do j=1,nond(ibnd)
            nd=iond(ibnd,j)
            if(n1==nd) j1=j
            if(n2==nd) j2=j
            if(j1/=0.and.j2/=0) exit
          enddo !j
          if(j1==0.or.j2==0) then
            write(11,*)'Open bnd side has non-bnd node:',i,ibnd,n1,n2
            stop
          endif

          do k=1,nvrt
            if(uthnd(ibnd,j1,k)<-98.or.uthnd(ibnd,j2,k)<-98.or.vthnd(ibnd,j1,k)<-98.or.vthnd(ibnd,j2,k)<-98) then
              write(11,*)'Wrong time series of vel.'
              stop
            endif
            uth(i,k)=ramp*(uthnd(ibnd,j1,k)+uthnd(ibnd,j2,k))/2
            vth(i,k)=ramp*(vthnd(ibnd,j1,k)+vthnd(ibnd,j2,k))/2
          enddo !k
        endif
      enddo !i=1,ns

      if(nscreen.eq.1) write(*,*)'done flow b.c.'
      write(16,*)'done flow b.c.'

!
!************************************************************************
!									*
!			Backtracking 					*
!									*
!************************************************************************
!

!     Debug: test backtracking alone
!       eta1=0; eta2=0; we=0
!       do i=1,ns
!        do k=1,nvrt
!          su2(k,i)=-ycj(i)*2*pi/3.0e3
!          sv2(k,i)=xcj(i)*2*pi/3.0e3
!        enddo !k
!      enddo !i
!      do i=1,ne
!        do k=1,nvrt
!          do j=1,3
!            nd=nm(i,j)
!            ufg(k,i,j)=-y(nd)*2*pi/3.0e3
!            vfg(k,i,j)=x(nd)*2*pi/3.0e3
!          enddo !j
!        enddo !k
!      enddo !i
!      do i=1,np
!        do k=1,nvrt
!          uu2(k,i)=-y(i)*2*pi/3.e3 
!          vv2(k,i)=x(i)*2*pi/3.e3
!          ww2(k,i)=0 !-1.e-4*z(k,i)*(50+z(k,i))
!        enddo !k
!      enddo !i
!     End debug

      call system_clock(ist,icount_rate)

!...  From nodes and sidecenters, and whole levels
!...  ptbt, sdbt: interpolated values at whole levels
!     Pre-assign for dry sides etc.
      do j=1,nvrt
        do i=1,ns
          sdbt(i,j,1)=su2(j,i)
          sdbt(i,j,2)=sv2(j,i)
        enddo !i
!       ptbt(1:2) currently not used
        do i=1,np
          ptbt(i,j,1)=uu2(j,i)
          ptbt(i,j,2)=vv2(j,i)
        enddo !i
      enddo !j

      ibt_p=0
      ibt_s=0
      do i=1,ne
        if(idry_e(i)==1) cycle

!   	wet elements (nodes and sides are wet as well)
        ie0=i

!       nodes, sides
        do l=1,6
!         Bypass nodes if upwind scheme is used for both S,T
          if(iupwind_t/=0.and.iupwind_s/=0.and.l<=3) cycle

          if(l<=3.and.ibt_p(nm(i,l))==1.or.l>3.and.l<=6.and.ibt_s(js(i,l-3))==1) cycle

          if(l<=3) then
            jmin=kbp(nm(i,l)) 
          else
            jmin=kbs(js(i,l-3)) 
          endif
          do j=jmin,nvrt 
!           Initialize (x0,y0,z0),nnel and vel.
!	    Caution! nnel must be initialized inside this loop as it is updated inside.
            if(l<=3) then !nodes
              nd0=nm(i,l)
              iadvf=iadv(nd0)
              x0=x(nd0)
              y0=y(nd0)
              z0=z(j,nd0)
              uuint=uu2(j,nd0)
              vvint=vv2(j,nd0)
              wwint=ww2(j,nd0)
            else !sides
              isd0=js(i,l-3)
              n1=isidenode(isd0,1)
              n2=isidenode(isd0,2)
              iadvf=min(iadv(n1),iadv(n2))
              x0=xcj(isd0)
              y0=ycj(isd0)
              z0=zs(j,isd0)
              uuint=su2(j,isd0)
              vvint=sv2(j,isd0)
              wwint=(ww2(j,n1)+ww2(j,n2))/2
            endif
            vmag=dsqrt(uuint**2+vvint**2+wwint**2)
            nnel=ie0
            jlev=j
!            jlev=min(j+1,nvrt) !make sure j>=2 for division()

            if(vmag<=1.e-4) then !No activity 
              if(l<=3) then
                ptbt(nd0,j,3)=tnd(j,nd0)
                ptbt(nd0,j,4)=snd(j,nd0)
                ptbt(nd0,j,1)=uu2(j,nd0)
                ptbt(nd0,j,2)=vv2(j,nd0)
!                x3bt(nd0,j,1)=x0
!                x3bt(nd0,j,2)=y0
!                x3bt(nd0,j,3)=z0
!                nelvbt(nd0,j,1)=nnel
!                nelvbt(nd0,j,2)=jlev
              else !sides
                sdbt(isd0,j,3)=tsd(j,isd0)
                sdbt(isd0,j,4)=ssd(j,isd0)
                sdbt(isd0,j,1)=su2(j,isd0)
                sdbt(isd0,j,2)=sv2(j,isd0)
              endif
            else !do btrack
              if(nadv>0) then
                dtb_max=dtb_max1
              else if(iadvf<=1) then !nadv=0
                dtb_max=dtb_max1
              else
                dtb_max=dtb_max2
              endif
              call btrack(i,l,j,iadvf,dtb_max,uuint,vvint,wwint,x0,y0,z0,nnel,jlev,xt,yt,zt,ndiv,swild)
              if(l<=3) then
                ptbt(nd0,j,3)=swild(1)
                ptbt(nd0,j,4)=swild(2)
                if(iadvf==0) then
                  ptbt(nd0,j,1)=uu2(j,nd0)
                  ptbt(nd0,j,2)=vv2(j,nd0)
                else
                  ptbt(nd0,j,1)=uuint
                  ptbt(nd0,j,2)=vvint
                endif
!                x3bt(nd0,j,1)=xt
!                x3bt(nd0,j,2)=yt
!                x3bt(nd0,j,3)=zt
!                nelvbt(nd0,j,1)=nnel
!                nelvbt(nd0,j,2)=jlev
              else !sides
                sdbt(isd0,j,3)=swild(1)
                sdbt(isd0,j,4)=swild(2)
                if(iadvf==0) then
                  sdbt(isd0,j,1)=su2(j,isd0)
                  sdbt(isd0,j,2)=sv2(j,isd0)
                else
                  sdbt(isd0,j,1)=uuint
                  sdbt(isd0,j,2)=vvint
                endif
              endif !sides
            endif !do backtrack

          enddo !j=1,nvrt

          if(l<=3) then
            ibt_p(nd0)=1
          else if(l<=6) then
            ibt_s(isd0)=1
          endif
        enddo !l=1,6
      enddo !i=1,ne

!     Debug
!      do i=1,np
!        th=pi/2+2*pi/3000*time
!        x0=1.8e3*cos(th)
!        y0=1.8e3*sin(th)
!        do k=1,nvrt
!          prho(i,k)=dexp(-((x(i)-x0)**2+(y(i)-y0)**2)/2/600/600) !exact soln
!        enddo !k
!      enddo !i

!...  Compute division pts 
!...  bubt: total integrated value
      do i=1,ne
        bubt(i,1)=0; bubt(i,2)=0
        do j=1,3 !sides
          isd=js(i,j)
          if(idry_s(isd)==0) then
            do k=kbs(isd)+1,nvrt !layer
              bubt(i,1)=bubt(i,1)+(sdbt(isd,k,1)+sdbt(isd,k-1,1))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
              bubt(i,2)=bubt(i,2)+(sdbt(isd,k,2)+sdbt(isd,k-1,2))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
            enddo !k
          endif
        enddo !j
      enddo !i=1,ne

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'btrack took',btimer,'seconds...'
      write(16,*)'backtracking took',btimer,'seconds...'

!     Density gradient at nodes and whole levels using cubic spline
      dr_ds=0 !for sigma_t
      if(ibc==0) then
        do i=1,np
          if(idry(i)==1) cycle
          if(prho(i,1)<-98) then
            write(11,*)'Impossible 4'
            stop
          endif

          if(kbp(i)==kz) then
            drds_b=0
          else !kbp < kz
            drds_b=(sig_t(i,kz+1)-sig_t(i,kz))/(sigma(2)-sigma(1))
          endif
          do k=1,nsig
            klev=k-1+kz !kz<= klev <=nvrt
            if(k==1) then
              bdia(k)=(sigma(k+1)-sigma(k))/3
              cupp(k)=bdia(k)/2
              rrhs(k,1)=(sig_t(i,klev+1)-sig_t(i,klev))/(sigma(k+1)-sigma(k))-drds_b
            else if(k==nsig) then
              bdia(k)=(sigma(k)-sigma(k-1))/3
              alow(k)=bdia(k)/2
              rrhs(k,1)=-(sig_t(i,klev)-sig_t(i,klev-1))/(sigma(k)-sigma(k-1))
            else
              bdia(k)=(sigma(k+1)-sigma(k-1))/3
              alow(k)=(sigma(k)-sigma(k-1))/6
              cupp(k)=(sigma(k+1)-sigma(k))/6
              rrhs(k,1)=(sig_t(i,klev+1)-sig_t(i,klev))/(sigma(k+1)-sigma(k))- &
     &(sig_t(i,klev)-sig_t(i,klev-1))/(sigma(k)-sigma(k-1))
            endif
          enddo !k
          call tridag(mnv,nsig,1,alow,bdia,cupp,rrhs,soln,gam)

          do k=1,nsig
            klev=k-1+kz !kz<= klev <=nvrt
            if(k==1) then
              dr_ds(i,klev)=drds_b
            else if(k==nsig) then
              dr_ds(i,klev)=0
            else
              dr_ds(i,klev)=(sig_t(i,klev+1)-sig_t(i,klev))/(sigma(k+1)-sigma(k))- &
     &(sigma(k+1)-sigma(k))/6*(2*soln(k,1)+soln(k+1,1))
            endif
          enddo !k
        enddo !i=1,np
      endif !ibc

      if(nscreen.eq.1) write(*,*)'done density gradient...'
      write(16,*)'done density gradient...'

!
!************************************************************************
!                                                                       *
!               Turbulence closure schemes                              *
!       Compute turbulence diffusivities dfv, dfh,                      *
!       and in MY-G, also dfq[1,2].                                     *
!                                                                       *
!************************************************************************
!

!...  Scheme 2: Pacanowski and Philander (1981)
      if(itur==2) then
        dfv=0; dfh=0 !for dry nodes
        do i=1,np
          if(idry(i)==1) cycle
          if(prho(i,1)<-98) then
            write(11,*)'Impossible dry 1'
            stop
          endif

!         wet nodes
          if(dp(i)<=h1_pp) then
            vmax=vdmax_pp1
            vmin=vdmin_pp1
            tmin=tdmin_pp1
          else if(dp(i)<h2_pp) then
            vmax=vdmax_pp1+(vdmax_pp2-vdmax_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            vmin=vdmin_pp1+(vdmin_pp2-vdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            tmin=tdmin_pp1+(tdmin_pp2-tdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
          else !dps >= h2
            vmax=vdmax_pp2
            vmin=vdmin_pp2
            tmin=tdmin_pp2
          endif

          do k=kbp(i),nvrt
            if(k==kbp(i).or.k==nvrt) then
              drhodz=0
            else
              drhodz=(prho(i,k+1)-prho(i,k-1))/(z(k+1,i)-z(k-1,i))
            endif
            bvf=-g*(drhodz/rho0+g/1.5e3**2)
            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(i))
            dudz=(su2(k2,i)-su2(k1,i))/(z(k2,i)-z(k1,i))
            dvdz=(sv2(k2,i)-sv2(k1,i))/(z(k2,i)-z(k1,i))
            shear2=dmax1(dudz**2+dvdz**2,1.0d-10) 
            rich=dmax1(bvf/shear2,0.0d0)

!           vmax >= vmin
            dfv(i,k)=vmax/(1+5*rich)**2+vmin
            dfh(i,k)=dfv(i,k)/(1+5*rich)+tmin
          enddo !k      
        enddo !i=1,np

        if(nscreen==1) write(*,*) 'done turbulence closure (PP)...'
        write(16,*) 'done turbulence closure (PP)...'
      endif !itur=2

!... Scheme 4: GOTM
!    In GOTM, all turbulence variables are defined at whole levels from bottom to F.S.
!    and mean flow variables at half levels. So the bottom is at level 0 (our kbp), 
!    F.S. is at level nlev (out nvrt).

      if(itur==4) then
#ifdef USE_GOTM
!        if(abs(cde-cmiu0**3)>1.e-4) then
!          write(11,*)'Mismatch in GOTM call:',cde,cmiu0**3
!          stop
!        endif
         write(16,*)'cde, cmiu0**3 = ',cde,cmiu0**3

        do j=1,np
!         Dry nodes will have initial values
          if(idry(j)==1) cycle
      
!         Friction velocity: [\niu*|du/dz|]^0.5 (m/s)
          u_taus=sqrt(tau(j,1)**2+tau(j,2)**2)
          u_taub=sqrt(Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2))
          nlev=nvrt-kbp(j)
          do k=0,nlev 
            klev=k+kbp(j) !kbp <= klev <= nvrt
            if(k/=0) h1d(k)=z(klev,j)-z(klev-1,j)
!           Shear frequency squared (1/s^2): (du/dz)^2+(dv/dz)^2
!           Buoyancy frequency squared (1/s^2): -g/\rho0*(d\rho/dz))
            if(k==0.or.k==nlev) then
              if(dfv(j,klev)<=0) then
                write(11,*)'Negative viscosity:',dfv(j,klev),j,klev
                stop
              endif
              if(k==0) then
                SS1d(k)=u_taub**2/dfv(j,klev)
              else
                SS1d(k)=u_taus**2/dfv(j,klev)
              endif
              NN1d(k)=0
            else
              ztmp=z(klev+1,j)-z(klev-1,j)
              if(ztmp==0) then
                write(11,*)'Zero layer:',j,klev
                stop
              endif
              SS1d(k)=((uu2(klev+1,j)-uu2(klev-1,j))**2+(vv2(klev+1,j)-vv2(klev-1,j))**2)/ztmp**2
              NN1d(k)=-g/rho0*(prho(j,klev+1)-prho(j,klev-1))/ztmp
            endif
            tke1d(k)=q2(j,klev)
            L1d(k)=xl(j,klev)
            if(tke1d(k)<=0.or.L1d(k)<=0) then
              write(11,*)'Negative tke,mixl:',tke1d(k),L1d(k),j,klev
              stop
            endif
            eps1d(k)=cde*tke1d(k)**1.5/L1d(k) 
            num1d(k)=dfv(j,klev)
            nuh1d(k)=dfh(j,klev)

!           Debug
!            write(98,*)k,h1d(k),NN1d(k),SS1d(k)

          enddo !k=0,nlev
!          h1d(0)=h1d(1)
          toth=eta2(j)+dp(j)
!         surface and bottom roughness length (m)
          z0s=min(0.1d0,toth/10)
          if(Cdp(j)==0) then
            z0b=0
          else
            z0b=(z(kbp(j)+1,j)-z(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
          endif

!         Debug
!          write(99,*)j,'WOW1'
!          write(98,*)nlev,dt,toth,u_taus,u_taub,z0s,z0b,h1d(0)


          call do_turbulence(nlev,dt,toth,u_taus,u_taub,z0s,z0b,h1d,NN1d,SS1d)

!         Debug
!          write(99,*)j,'WOW2'

          q2(j,kbp(j):nvrt) = tke1d(0:nlev)
          xl(j,kbp(j):nvrt) = L1d(0:nlev)
!          eps(i,j,:) = eps1d
          do k=0,nlev
            klev=k+kbp(j)
            dfv(j,klev)=dmin1(diffmax(j),num1d(k)+diffmin(j)) 
            dfh(j,klev)=dmin1(diffmax(j),nuh1d(k)+diffmin(j))
          enddo !k
        enddo !j=1,np
#endif
      endif !itur==4
 
!... Scheme 3: Mellor-Yamada-Galperin & Umlauf-Burchard scheme
      if(itur==3) then
!------------------------------------------------------------
      call system_clock(ist,icount_rate)
      do j=1,np
        if(idry(j)==1) then
          do k=1,nvrt
            q2(j,k)=q2min; xl(j,k)=xlmin2(j)
            dfv(j,k)=0; dfh(j,k)=0; dfq1(j,k)=0; dfq2(j,k)=0
          enddo 
          cycle
        endif
        if(prho(j,1)<-98) then
          write(11,*)'Impossible dry 2'
          stop
        endif

!       Wet node; compute layer thickness etc.
!       Error: use ufg?
        do k=kbp(j)+1,nvrt
          dzz(k)=z(k,j)-z(k-1,j)
          dudz=(uu2(k,j)-uu2(k-1,j))/dzz(k)
          dvdz=(vv2(k,j)-vv2(k-1,j))/dzz(k)
          shearbt(k)=dudz**2+dvdz**2 !@ half levels
          rzbt(k)=g/rho0*(prho(j,k)-prho(j,k-1))/dzz(k)
          q2ha(k)=(q2(j,k)+q2(j,k-1))/2
          xlha(k)=(xl(j,k)+xl(j,k-1))/2

!         Compute c_psi_3
          if(mid.eq.'MY') then
            cpsi3(k)=0.9
          else !GLS models
            if(rzbt(k)>0) then !unstable
              cpsi3(k)=1
            else !stable
              select case(mid)
                case('KL')
                  cpsi3(k)=2.53
                case('KE')
                  cpsi3(k)=-0.52
                case('KW')
                  cpsi3(k)=-0.58
                case('UB')
                  cpsi3(k)=0.1
                case default
                  write(11,*)'Unknown closure model:',mid
                  stop
              end select
            endif
          endif

!         Wall proximity function      
          if(mid.eq.'MY'.or.mid.eq.'KL') then
            zctr=(z(k,j)+z(k-1,j))/2
            dists=eta2(j)-zctr
            distb=zctr+dp(j)
            if(dists==0.or.distb==0) then
              write(11,*)'Zero in proximity function:',j,k
              stop
            endif
            fwall=1+1.33*(xlha(k)/0.4/distb)**2+0.25*(xlha(k)/0.4/dists)**2
            cpsi2p(k)=fwall*cpsi2 !F_wall*cpsi2
          else !other GLS
            cpsi2p(k)=cpsi2
          endif
        enddo !k=kbp(j)+1,nvrt
        rzbt(kbp(j))=0 !for Galperin's clipping

!        write(90,*)'WOW1',it,j

!	Compute upper bound for xl 
        do k=kbp(j),nvrt
          dists=eta2(j)-z(k,j)
          distb=z(k,j)+dp(j)
          if(k==kbp(j)) then
            xlmax(k)=dmax1(xlmin2(j),dzz(k+1)*0.4)
          else if(k==nvrt) then
            xlmax(k)=dmax1(xlmin2(j),dzz(k)*0.4)
          else !internal layers
            xlmax(k)=0.4*dmin1(dists,distb)
          endif
!          xlmax(k)=dmax1(0.4*dmin1(dists,distb),xlmin2(j)) !can be very small
!          xlmax(k)=0.4*dists*distb/(dps(j)+etam)
!          xlmax(k)=0.4*dmin1(dp(j)+eta2(j),xlmax00)
          if(xlmax(k)<=0) then
            write(11,*)'Dist<0 in MY-G',j,k,eta2(j)+dp(j),dists,distb
            stop
          endif
        enddo !k

!	b.c. (computed using values from previous time except wind)
        q2fs=16.6**(2.0/3)*dsqrt(tau(j,1)**2+tau(j,2)**2)/2
        q2fs=dmax1(q2fs,q2min)
        q2bot=16.6**(2.0/3)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)/2
        q2bot=dmax1(q2bot,q2min)
        xlfs=dmax1(xlmin2(j),xlsc0(j)*dzz(nvrt)*0.4) 
        xlbot=dmax1(xlmin2(j),dmin1(2.5d0,xlsc0(j)*dzz(kbp(j)+1))*0.4) !"5" to prevent over-mixing

!        write(90,*)'WOW2',it,j

!	Matrix Q
        nqdim=nvrt-kbp(j)+1
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1 !row #
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(kin,1)=0
          if(k<nvrt) then
            tmp=(dfq1(j,k+1)+dfq1(j,k))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*q2(j,k)+q2(j,k+1))
            prod=(dfv(j,k+1)+dfv(j,k))/2*shearbt(k+1)
            buoy=(dfh(j,k+1)+dfh(j,k))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k+1)/2*(prod+buoy)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cmiu0**3*dsqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          endif

          if(k>kbp(j)) then
            tmp=(dfq1(j,k)+dfq1(j,k-1))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*q2(j,k)+q2(j,k-1))
            prod=(dfv(j,k)+dfv(j,k-1))/2*shearbt(k)
            buoy=(dfh(j,k)+dfh(j,k-1))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k)/2*(prod+buoy)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cmiu0**3*dsqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          endif
        enddo !k=1,nvrt

!	Soln for q2 at new level
        call tridag(mnv,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          if(k==nvrt) then
            q2tmp(k)=q2fs
          else if(k==kbp(j)) then
            q2tmp(k)=q2bot
          else
            q2tmp(k)=dmax1(soln(kin,1),q2min)
          endif
        enddo !k

!        write(90,*)'WOW4',it,j,(q2tmp(k),k=1,nvrt)
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Matrix QL
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(kin,1)=0
          if(k<nvrt) then
            tmp=(dfq2(j,k+1)+dfq2(j,k))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            psi_n=cmiu0**rpub*q2(j,k)**rmub*xl(j,k)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(j,k+1)**rmub*xl(j,k+1)**rnub !psi^n_{j,k+1}
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(j,k+1)+dfv(j,k))/2*shearbt(k+1)
            buoy=cpsi3(k+1)*(dfh(j,k+1)+dfh(j,k))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k+1)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k+1)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cpsi2p(k+1)*cmiu0**3*dsqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          else !k=nvrt
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(j,k)/xl(j,k)
          endif

          if(k>kbp(j)) then 
            tmp=(dfq2(j,k)+dfq2(j,k-1))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            psi_n=cmiu0**rpub*q2(j,k)**rmub*xl(j,k)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(j,k-1)**rmub*xl(j,k-1)**rnub !psi^n_{j,k-1}
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(j,k)+dfv(j,k-1))/2*shearbt(k)
            buoy=cpsi3(k)*(dfh(j,k)+dfh(j,k-1))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cpsi2p(k)*cmiu0**3*dsqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          else !k=kbp(j)
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(j,k)/xl(j,k)
          endif
        enddo !k=kbp(j),nvrt

!        write(90,*)'WOW5',it,j
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Soln for q2l and xl at new level
        call tridag(mnv,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)

!        write(90,*)'WOW6',it,j

        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          q2l=dmax1(soln(kin,1),psimin)
          if(k==nvrt) then
            xltmp(k)=xlfs
          else if(k==kbp(j)) then
            xltmp(k)=xlbot
          else
            xltmp(k)=(q2l*cmiu0**(-rpub)*q2tmp(k)**(-rmub))**(1/rnub)
          endif
!	  Galperin's clipping 
          if(rzbt(k)<0) then
            upper=dsqrt(-0.56*q2tmp(k)/rzbt(k))
            xltmp(k)=dmin1(xltmp(k),upper)
          endif
!	  Max. length based on dissipation; xlmin2 prevails
          xl_max=(cmiu0*dsqrt(q2tmp(k)))**3/eps_min
          xltmp(k)=dmax1(xlmin2(j),dmin1(xl_max,xltmp(k)))
!	  Impose max. depth limit
          xltmp(k)=dmax1(xlmin2(j),dmin1(xltmp(k),xlmax(k)))

          q2(j,k)=q2tmp(k)
          xl(j,k)=xltmp(k)
          if(q2(j,k)<0) then
            write(11,*)'Negative q2',q2(j,k),xl(j,k)
            stop
          endif

!         Compute vertical diffusivities at new time
          call asm(g,j,k,vd,td,qd1,qd2)
          dfv(j,k)=dmin1(diffmax(j),dmax1(diffmin(j),vd))
          dfh(j,k)=dmin1(diffmax(j),dmax1(diffmin(j),td))
          dfq1(j,k)=dmin1(diffmax(j),dmax1(diffmin(j),qd1))
          dfq2(j,k)=dmin1(diffmax(j),dmax1(diffmin(j),qd2))

!         Debug
!          write(90,*)'No. ',k,xl(j,k),dfh(j,k),dfv(j,k),dfq1(j,k),dfq2(j,k)
        enddo !k

!       Extend
        do k=1,kbp(j)-1
          q2(j,k)=q2(j,kbp(j))
          xl(j,k)=xl(j,kbp(j))
          dfv(j,k)=dfv(j,kbp(j))
          dfh(j,k)=dfh(j,kbp(j))
          dfq1(j,k)=dfq1(j,kbp(j))
          dfq2(j,k)=dfq2(j,kbp(j))
        enddo !k
      enddo !j=1,np

!      if(it.eq.1739) write(90,*)'WOW7',it

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'MYG-UB took',btimer,'seconds'
      write(16,*)'MYG-UB took',btimer,'seconds'

!------------------------------------------------------------
      endif !itur=3

!
!************************************************************************
!									*
!		Wave-continuity equations				*
!									*
!************************************************************************
!

      call system_clock(ist,icount_rate)

!...  compute elevation essential boundary conditions
!...
      elbc=-9999 !flags
      do i=1,nope
        do j=1,nond(i)
          nd=iond(i,j)
          if(iettype(i)==1.or.iettype(i)==2) then
            elbc(nd)=ramp*eth(i,1)
          else if(iettype(i)==3) then
            elbc(nd)=0 !initialize
            do jfr=1,nbfr
              ncyc=int(amig(jfr)*time/2/pi)
              arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-efa(i,j,jfr)
              elbc(nd)=elbc(nd)+ramp*ff(jfr)*emo(i,j,jfr)*dcos(arg)
            enddo !jfr=1,nbfr
          else if(iettype(i)==4) then
            elbc(nd)=ramp*eth(i,j)
          endif
        enddo !j=1,noe(i)
      enddo !i=1,nope

!...  Pre-compute some arrays: chi,hhat,bigu,ghat1
!...
      do i=1,ns
        if(idry_s(i)==1) then
          chi(i)=0
          hhat(i)=0
          bigu(i,1)=0
          bigu(i,2)=0
          cycle
        endif

!	Wet side
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        if(idrag==1) then
          chi(i)=Cd(i)
        else
          chi(i)=Cd(i)*dsqrt(sdbt(i,kbs(i)+1,1)**2+sdbt(i,kbs(i)+1,2)**2)
        endif
        hhat(i)=(eta2(n1)+eta2(n2))/2+dps(i)-chi(i)*dt
!	Enforce positivity
        if(ihhat==1) hhat(i)=dmax1(0.d0,hhat(i))

!	bigu1,2
        bigu(i,1)=0 !U^n_x
        bigu(i,2)=0 !U^n_y
        do k=kbs(i),nvrt-1
          bigu(i,1)=bigu(i,1)+(zs(k+1,i)-zs(k,i))*(su2(k,i)+su2(k+1,i))/2
          bigu(i,2)=bigu(i,2)+(zs(k+1,i)-zs(k,i))*(sv2(k,i)+sv2(k+1,i))/2
        enddo !k
      enddo !i=1,ns

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen==1) write(*,*)'1st preparation took ',btimer,'seconds'
      write(16,*)'1st preparation took ',btimer,'seconds'
        
      call system_clock(ist,icount_rate)
!     ghat1
      do i=1,ne
        if(idry_e(i)==1) then
          ghat1(i,1)=0
          ghat1(i,2)=0
          cycle
        endif

!	Wet elements
!	Excluding hvis, baroclinc force first
!       Warning: \hat{G}_1 must include all: Coriolis, atmo. pressure, tidal potential, horizontal difusion, and baroclinic
!              Remember to update both f (botf) and F (bigf)
        tau_x=0
        tau_y=0
        detadx=0
        detady=0
        dprdx=0
        dprdy=0
        detpdx=0
        detpdy=0
        chigamma=0
        ubstar=0 !bottom advection * \chi
        vbstar=0
        hhat_bar=0
        h_bar=0
        bigf1=0 !all in F except baroclinic and hvis
        bigf2=0
        botf1=0 !all in \chi*f_b except baroclinic and hvis
        botf2=0 
        do j=1,3 !node or side
          nd=nm(i,j)
          tau_x=tau_x+tau(nd,1)/3
          tau_y=tau_y+tau(nd,2)/3
!         idry_e(i) checked already
          detadx=detadx+eta2(nd)*dl(i,j,1)
          detady=detady+eta2(nd)*dl(i,j,2)
          dprdx=dprdx+pr(nd)*dl(i,j,1)
          dprdy=dprdy+pr(nd)*dl(i,j,2)
          if(dpe(i)>=tip_dp) then
            detpdx=detpdx+etp(nd)*dl(i,j,1)
            detpdy=detpdy+etp(nd)*dl(i,j,2)
          endif
          h_bar=h_bar+(dp(nd)+eta2(nd))/3

          isd=js(i,j)
          chigamma=chigamma+chi(isd)/3
          hhat_bar=hhat_bar+hhat(isd)/3
          ubstar=ubstar+chi(isd)*sdbt(isd,kbs(isd)+1,1)/3
          vbstar=vbstar+chi(isd)*sdbt(isd,kbs(isd)+1,2)/3
          bigf1=bigf1+cori(isd)*bigu(isd,2)/3
          bigf2=bigf2-cori(isd)*bigu(isd,1)/3
          botf1=botf1+chi(isd)*cori(isd)*sv2(kbs(isd)+1,isd)/3
          botf2=botf2-chi(isd)*cori(isd)*su2(kbs(isd)+1,isd)/3
        enddo !j=1,3
        bigf1=bigf1+h_bar*(0.69*g*detpdx-dprdx/rho0)
        bigf2=bigf2+h_bar*(0.69*g*detpdy-dprdy/rho0)
        botf1=botf1+chigamma*(0.69*g*detpdx-dprdx/rho0)
        botf2=botf2+chigamma*(0.69*g*detpdy-dprdy/rho0)

        ghat1(i,1)=bubt(i,1)+area(i)*dt*(bigf1+tau_x-ubstar-dt*botf1-g*(1-thetai)*hhat_bar*detadx)
        ghat1(i,2)=bubt(i,2)+area(i)*dt*(bigf2+tau_y-vbstar-dt*botf2-g*(1-thetai)*hhat_bar*detady)

!       Horizontal diffusion
        horx=0 
        hory=0
        do j=1,3 !side
          isd=js(i,j)
          do k=kbs(isd)+1,nvrt
            horx=horx+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2u(k,isd)+d2u(k-1,isd))/2
            hory=hory+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2v(k,isd)+d2v(k-1,isd))/2
          enddo !k
          horx=horx-dt*chigamma*area(i)/3*d2u(kbs(isd)+1,isd)
          hory=hory-dt*chigamma*area(i)/3*d2v(kbs(isd)+1,isd)
        enddo !j=1,3

        ghat1(i,1)=ghat1(i,1)+dt*horx
        ghat1(i,2)=ghat1(i,2)+dt*hory

!       Baroclinic force
        if(ibc==0) then
          if(prho(nm(i,1),1)<-98.or.prho(nm(i,2),1)<-98.or.prho(nm(i,3),1)<-98) then
            write(11,*)'Impossible dry 5'
            stop
          endif
!         Density is defined at 3 nodes
!         Area integrals evaluated using z- or sigma-
          node1=nm(i,1); node2=nm(i,2); node3=nm(i,3)
          hmax=dmax1(hmod(node1),hmod(node2),hmod(node3))
!         Initialize starting searching level for z- method
          do j=1,3
            nwild(2*j-1)=1 !starting level for j @ j+1
            nwild(2*j)=1 !j @ j+2
          enddo !j

          do k=kz,nvrt !S-levels first
            kin=k-kz+1 !S-index
            hbar=0 !average of dz_ds
            do j=1,3
              nd=nm(i,j)
              if(iback(nd)/=0) then !use traditional sigma
                dzds=eta2(nd)+hmod(nd)
              else
                dzds=eta2(nd)+h_c+(hmod(nd)-h_c)*dcs(kin)
              endif
              hbar=hbar+dzds/3
              out3(k,j)=dzds !save space
            enddo !j

            if(hmax<=h_bcc1) then !sigma-
              rrhs(k,1)=0 !dr_dx
              rrhs(k,2)=0 !dr_dy
              rrhs(k,3)=0 !dzdx
              rrhs(k,4)=0 !dzdy
              av=0 !average dr_ds
              sig_tm=(sig_t(node1,k)+sig_t(node2,k)+sig_t(node3,k))/3 !for removal of mean
              zmean=(z(k,node1)+z(k,node2)+z(k,node3))/3
              do j=1,3
                nd=nm(i,j)
                rrhs(k,1)=rrhs(k,1)+(sig_t(nd,k)-sig_tm)*dl(i,j,1)
                rrhs(k,2)=rrhs(k,2)+(sig_t(nd,k)-sig_tm)*dl(i,j,2)
                rrhs(k,3)=rrhs(k,3)+(z(k,nd)-zmean)*dl(i,j,1)
                rrhs(k,4)=rrhs(k,4)+(z(k,nd)-zmean)*dl(i,j,2)
                av=av+dr_ds(nd,k)/3
              enddo !j
              alow(kin)=-g/rho0*area(i)*(hbar*rrhs(k,1)-rrhs(k,3)*av) !integrand M * (-g)/rho0
              bdia(kin)=-g/rho0*area(i)*(hbar*rrhs(k,2)-rrhs(k,4)*av)

            else !z-
              rrhs(k,1)=0 !average dr_dx
              rrhs(k,2)=0 !average dr_dy
              icount=0 !valid nodes for drho
              do j=1,3 !node
                n1=nm(i,j)
                ifl=0 !flag to indicate valid computation

                do l=1,2 !other 2 nodes
                  nd=nm(i,nx(j,l))
                  isd=js(i,nx(j,3-l))
                  if(isidenode(isd,1)==n1) then
                    fac=1
                  else
                    fac=-1
                  endif
                  out2(2*l+1)=-fac*sny(isd) !l_x
                  out2(2*l+2)=fac*snx(isd) !l_y

                  if(z(k,n1)<z(kbp(nd),nd)) then
                    if(k==kbp(n1).or.k==nvrt) then
                      ifl=1
                      exit
                    else
                      zrat=(z(k,n1)-z(kbp(n1),n1))/(z(kbp(nd),nd)-z(kbp(n1),n1))
                      if(zrat<=0.or.zrat>=1) then
                        write(11,*)'Impossible 69:',zrat
                        stop
                      endif
                      rho_tmp=zrat*sig_t(nd,kbp(nd))+(1-zrat)*sig_t(n1,kbp(n1))
                      rl=zrat*distj(isd)
                      out2(l)=(rho_tmp-sig_t(n1,k))/rl
                    endif
                  else if(z(k,n1)>z(nvrt,nd)) then
                    if(k==kbp(n1).or.k==nvrt) then
                      ifl=1
                      exit
                    else
                      zrat=(z(nvrt,n1)-z(k,n1))/(z(nvrt,n1)-z(nvrt,nd))
                      if(zrat<=0.or.zrat>=1) then
                        write(11,*)'Impossible 69b:',zrat
                        stop
                      endif
                      rho_tmp=zrat*sig_t(nd,nvrt)+(1-zrat)*sig_t(n1,nvrt)
                      rl=zrat*distj(isd)
                      out2(l)=(rho_tmp-sig_t(n1,k))/rl
                    endif
                  else !must have a valid level
                    lev=0 !flag
                    nwild(2*j-2+l)=max0(nwild(2*j-2+l),kbp(nd))
                    do kk=nwild(2*j-2+l),nvrt-1
                      if(z(k,n1)>=z(kk,nd).and.z(k,n1)<=z(kk+1,nd)) then
                        lev=kk
                        zrat=(z(k,n1)-z(kk,nd))/(z(kk+1,nd)-z(kk,nd))
                        if(zrat<0.or.zrat>1) then
                          write(11,*)'Impossible 70:',zrat
                          stop
                        endif
                        exit
                      endif
                    enddo !kk
                    if(lev==0) then
                      write(11,*)'Failed to find a level in bcc:',n1,nd
                      stop
                    endif
                    nwild(2*j-2+l)=lev !for next k-iteration
                    rho_tmp=zrat*sig_t(nd,lev+1)+(1-zrat)*sig_t(nd,lev)
                    out2(l)=(rho_tmp-sig_t(n1,k))/distj(isd)
                  endif
                enddo !l=1,2

                if(ifl==0) then
                  icount=icount+1
                  delta=out2(3)*out2(6)-out2(5)*out2(4)
                  if(delta==0) then
                    write(11,*)'Ill formed element:',i
                    stop
                  endif
                  rrhs(k,1)=rrhs(k,1)+(out2(6)*out2(1)-out2(4)*out2(2))/delta
                  rrhs(k,2)=rrhs(k,2)+(out2(3)*out2(2)-out2(5)*out2(1))/delta
                endif !ifl==0
              enddo !j=1,3; nodes
              if(icount/=0) then
                rrhs(k,1)=rrhs(k,1)/icount
                rrhs(k,2)=rrhs(k,2)/icount
              endif

              alow(kin)=-g/rho0*area(i)*hbar*rrhs(k,1)
              bdia(kin)=-g/rho0*area(i)*hbar*rrhs(k,2)

            endif !sigma or z
          enddo !k=kz,nvrt

          do k=1,nsig-1 !S-index
            if(mmm==0) then !trapzoidal rule
              soln(k,1)=(alow(k+1)+alow(k))/2*(sigma(k+1)-sigma(k))
              soln(k,2)=(bdia(k+1)+bdia(k))/2*(sigma(k+1)-sigma(k))
            else !Lagrangian
              soln(k,1)=rint_lag(mnv,1,nsig,mmm,k,sigma,sigmap,sigma_prod,alow,gam,ctmp)
              soln(k,2)=rint_lag(mnv,1,nsig,mmm,k,sigma,sigmap,sigma_prod,bdia,gam,ctmp)
            endif
          enddo !k

          do k=kz,nvrt
            hp_int(k,i,1)=0 !\int f_{c} d\Omega
            hp_int(k,i,2)=0
            do kk=k,nvrt-1
              kin=kk-kz+1 !S-index
              hp_int(k,i,1)=hp_int(k,i,1)+soln(kin,1)
              hp_int(k,i,2)=hp_int(k,i,2)+soln(kin,2)
            enddo !kk
          enddo !k

!         Integrand for ghat1
          do l=1,nsig-1 !integrand=0 when l=nsig
            lev=l-1+kz
            do m=l,nsig
              mlev=m-1+kz
              sum21=0 !double sum 1
              sum22=0 !double sum 2
              do j=1,3
                do jj=1,3
                  if(j==jj) then
                    fac=2
                  else
                    fac=1
                  endif
                  sum21=sum21+fac*out3(lev,j)*out3(mlev,jj)
                  sum22=sum22+fac*out3(lev,j)*dr_ds(nm(i,jj),mlev)
                enddo !jj
              enddo !j

              if(hmax<=h_bcc1) then !sigma-
                alow(m)=-g/rho0*area(i)/12*(rrhs(mlev,1)*sum21-rrhs(mlev,3)*sum22) !integrand N
                bdia(m)=-g/rho0*area(i)/12*(rrhs(mlev,2)*sum21-rrhs(mlev,4)*sum22) !integrand N
              else !z
                alow(m)=-g/rho0*area(i)/12*sum21*rrhs(mlev,1)
                bdia(m)=-g/rho0*area(i)/12*sum21*rrhs(mlev,2)
              endif
            enddo !m
            cupp(l)=0 !outer integrand
            rrhs(l,5)=0 !outer integrand
            do m=l,nsig-1
              if(mmm==0) then !trapzoidal rule
                cupp(l)=cupp(l)+(alow(m+1)+alow(m))/2*(sigma(m+1)-sigma(m))
                rrhs(l,5)=rrhs(l,5)+(bdia(m+1)+bdia(m))/2*(sigma(m+1)-sigma(m))
              else !Lagrangian
                cupp(l)=cupp(l)+rint_lag(mnv,l,nsig,mmm,m,sigma,sigmap,sigma_prod,alow,gam,ctmp)
                rrhs(l,5)=rrhs(l,5)+rint_lag(mnv,l,nsig,mmm,m,sigma,sigmap,sigma_prod,bdia,gam,ctmp)
              endif
            enddo !m
          enddo !l=1,nsig-1
          cupp(nsig)=0; rrhs(nsig,5)=0

          bigfc1=0
          do l=1,nsig-1
            if(mmm==0) then !trapzoidal rule
              bigfc1=bigfc1+(cupp(l+1)+cupp(l))/2*(sigma(l+1)-sigma(l))
            else !Lagrangian
              bigfc1=bigfc1+rint_lag(mnv,1,nsig,mmm,l,sigma,sigmap,sigma_prod,cupp,gam,ctmp)
            endif
          enddo !l
          do l=1,nsig
            cupp(l)=rrhs(l,5)
          enddo !l
          bigfc2=0
          do l=1,nsig-1
            if(mmm==0) then !trapzoidal rule
              bigfc2=bigfc2+(cupp(l+1)+cupp(l))/2*(sigma(l+1)-sigma(l))
            else !Lagrangian
              bigfc2=bigfc2+rint_lag(mnv,1,nsig,mmm,l,sigma,sigmap,sigma_prod,cupp,gam,ctmp)
            endif
          enddo !l

!         z-levels
          if(kbe(i)<kz) then
!           hp_int
            do k=kbe(i),kz !all 3 nodes have z-levels
              if(k==kbe(i)+1) then
                soln(k,3)=dpe(i)+ztot(k) !thickness
              else if(k>kbe(i)+1) then
                soln(k,3)=ztot(k)-ztot(k-1) !thickness
              endif
              if(k/=kbe(i).and.soln(k,3)<=0) then
                write(11,*)'Thickness <=0:',i,k,soln(k,3)
                stop
              endif

              soln(k,1)=0 !dr_dx
              soln(k,2)=0 !dr_dy
              sig_tm=(sig_t(node1,k)+sig_t(node2,k)+sig_t(node3,k))/3 !for removal of mean
              do j=1,3 !nodes
                nd=nm(i,j)
                soln(k,1)=soln(k,1)+(sig_t(nd,k)-sig_tm)*dl(i,j,1)
                soln(k,2)=soln(k,2)+(sig_t(nd,k)-sig_tm)*dl(i,j,2)
!                if(k/=kbe(i)) then
!                  if(k-1<kbp(nd)) then
!                    write(11,*)'Out of bound (1):',i,j,nd
!                    stop
!                  endif
!                  soln(k,3)=soln(k,3)+(z(k,nd)-z(k-1,nd))/3
!                endif
              enddo !j
            enddo !k

            do k=kbe(i),kz
              soln(k,4)=0 !P_k
              soln(k,5)=0
              do m=k,kz-1
                soln(k,4)=soln(k,4)-g/rho0*area(i)*soln(m+1,3)*(soln(m,1)+soln(m+1,1))/2
                soln(k,5)=soln(k,5)-g/rho0*area(i)*soln(m+1,3)*(soln(m,2)+soln(m+1,2))/2
              enddo !m
            enddo !k

            do k=kbe(i),kz-1
              hp_int(k,i,1)=hp_int(kz,i,1)+soln(k,4)
              hp_int(k,i,2)=hp_int(kz,i,2)+soln(k,5)
            enddo !k

!           bigfc[1,2]
            do k=kbe(i),kz-1
              bigfc1=bigfc1+soln(k+1,3)/2*(hp_int(k+1,i,1)+hp_int(k,i,1))
              bigfc2=bigfc2+soln(k+1,3)/2*(hp_int(k+1,i,2)+hp_int(k,i,2))
            enddo !k
          endif !kbe(i)<kz

!         Debug
!          if(i==5159) then
!            write(93,*)bigfc1,bigfc2,chigamma,ghat1(i,1),ghat1(i,2)
!            do k=kbe(i),nvrt
!              write(93,*)k,hp_int(k,i,1),hp_int(k,i,2),z(k,node1)
!            enddo !k
!          endif

          ghat1(i,1)=ghat1(i,1)+rampbc*dt*(bigfc1-chigamma*dt*hp_int(kbe(i)+1,i,1))
          ghat1(i,2)=ghat1(i,2)+rampbc*dt*(bigfc2-chigamma*dt*hp_int(kbe(i)+1,i,2))
        endif !ibc==0
      enddo !i=1,ne

!...  Baroclinic force at side and whole levels
      if(ibc==0) then
        bcc=0 !x, y-component in global frame
        do i=1,ns
          if(idry_s(i)==1) cycle

!         Wet side
!         Caution: icase=2 only works for pure S region
          icase=1
          if(is(i,2)==0) then
            if(interpol(is(i,1))==2) icase=2
          else
            if(interpol(is(i,1))==2.and.interpol(is(i,2))==2) icase=2
          endif

          do k=kbs(i),nvrt
            ta=0
            do j=1,2
              ie=is(i,j)
              if(ie/=0.and.idry_e(ie)==0) then
                if(icase==1) then
                  kbb=kbe(ie)
                  swild2(kbb:nvrt,1)=hp_int(kbb:nvrt,ie,1)
                  swild2(kbb:nvrt,2)=hp_int(kbb:nvrt,ie,2)
                  alow(kbb:nvrt)=ze(kbb:nvrt,ie) 
                  call vinter(mnv,2,zs(k,i),kbe(ie),nvrt,k,alow,swild2,swild,ibelow)
                  bcc(i,k,1:2)=bcc(i,k,1:2)+swild(1:2)
                else !Error: in pure S region
                  km=max0(k,kbe(ie))
                  bcc(i,k,1:2)=bcc(i,k,1:2)+hp_int(km,ie,1:2)
                endif

                ta=ta+area(ie)
              endif
            enddo !j=1,2
            if(ta==0) then
              bcc(i,k,1:2)=0
            else
              bcc(i,k,1:2)=bcc(i,k,1:2)/ta*rampbc
            endif
          enddo !k
        enddo !i=1,ns
      endif !ibc==0

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen==1) write(*,*)'2nd preparation took ',btimer,'seconds'
      write(16,*)'2nd preparation took ',btimer,'seconds'

      call system_clock(ist,icount_rate)
!...  setup coefficient matrix, sparsem, for the wave equation
!...  No elevation essential b.c. are imposed yet
      do i=1,np
        do j=0,nnp(i)
          sparsem(i,j)=0
        enddo !j
        qel(i)=0

!	Area integrals I_{1,4}
        do j=1,nne(i)
          ie=ine(i,j)
          id=iself(i,j)

!	  I_1
          n2=nm(ie,nx(id,1))
          n3=nm(ie,nx(id,2))
          dot1=(x(n3)-x(n2))**2+(y(n3)-y(n2))**2
          dot2=(x(n3)-x(n2))*(x(i)-x(n3))+(y(n3)-y(n2))*(y(i)-y(n3))
          dot3=-dot1-dot2
          hhatb=(hhat(js(ie,1))+hhat(js(ie,2))+hhat(js(ie,3)))/3
          tmp0=area(ie)/6+g*thetai**2*dt**2/4/area(ie)*hhatb*dot1
          tmpj=area(ie)/12+g*thetai**2*dt**2/4/area(ie)*hhatb*dot2
          tmpj1=area(ie)/12+g*thetai**2*dt**2/4/area(ie)*hhatb*dot3
          sparsem(i,0)=sparsem(i,0)+tmp0
          sparsem(i,j)=sparsem(i,j)+tmpj
          if(isbnd(i)==0.and.j==nne(i)) then
            sparsem(i,1)=sparsem(i,1)+tmpj1
          else
            sparsem(i,j+1)=sparsem(i,j+1)+tmpj1
          endif
!	  Check dominance
          if(hhatb<0) then
            if(ihhat==0.and.ifort12(1)==0) then
              ifort12(1)=1
              write(12,*)'Modified depth < 0:',it,i,j,hhatb
            endif
            if(ihhat==1) then
              write(11,*)'Impossible hhat:',hhatb
              stop
            endif
          endif

!	  I_4
          isd1=js(ie,1)
          isd2=js(ie,2)
          isd3=js(ie,3)
          dot1=dl(ie,id,1)*(bigu(isd1,1)+bigu(isd2,1)+bigu(isd3,1))/3+ &
     &dl(ie,id,2)*(bigu(isd1,2)+bigu(isd2,2)+bigu(isd3,2))/3
          dot2=dl(ie,id,1)*ghat1(ie,1)+dl(ie,id,2)*ghat1(ie,2)
        
          qel(i)=qel(i)+(1-thetai)*dt*area(ie)*dot1+thetai*dt*dot2
          do l=1,3
            if(id==l) then
              fac=2
            else
              fac=1
            endif
            nd=nm(ie,l)
            qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+bdef2(nd)-bdef1(nd))
          enddo !l
        enddo !j=1,nne(i)

!	bnd integrals I_{2,3,5,6}; they all vanish at land bnds 
!	I_2,6 are not needed if essential b.c. are enforced by elminating rows and columns
        if(isbnd(i)>0) then !open bnd node
          ibnd=isbnd(i)
          do l=1,2 !side
            if(l==1) then
              ie=ine(i,1)
              id=iself(i,1)
              isd=js(ie,nx(id,2))
              nj=nm(ie,nx(id,1))
              ind=1
            else
              ie=ine(i,nne(i))
              id=iself(i,nne(i))
              isd=js(ie,nx(id,1))
              nj=nm(ie,nx(id,2))
              ind=nnp(i)
            endif

            nd=isidenode(isd,1)+isidenode(isd,2)-i
            if(nd/=nj) then
              write(11,*)'Impossible 79'
              stop
            endif

!	    I_3 
            if(isbs(isd)>0.and.ifltype(isbs(isd))/=0) then !natural or Flather 1 b.c.
              if(idry_s(isd)==1) then
                write(11,*)'Dry flow bnd:',isd,i,nd
                stop
              endif

              bigvn=0
              do k=kbs(isd),nvrt-1
                vn1=uth(isd,k)*snx(isd)+vth(isd,k)*sny(isd)
                vn2=uth(isd,k+1)*snx(isd)+vth(isd,k+1)*sny(isd)
                bigvn=bigvn+(zs(k+1,isd)-zs(k,isd))*(vn1+vn2)/2
              enddo !k
              ri3=distj(isd)*bigvn/2
              if(ifltype(isbs(isd))==-1) then !Flather 1
                if(eta_mean(i)<-98.or.eta_mean(nj)<-98) then
                  write(11,*)'Mismatch 1'
                  stop
                endif
                if(dps(isd)<=0) then
                  write(11,*)'Negative depth at Flather bnd:',i,dps(isd)
                  stop
                endif
                con0=distj(isd)/6*dsqrt(g*dps(isd)) !for coefficient matrix
                ri3=ri3-con0*(2*eta_mean(i)+eta_mean(nj))
                sparsem(i,0)=sparsem(i,0)+thetai*dt*con0*2
                sparsem(i,ind)=sparsem(i,ind)+thetai*dt*con0
              endif !Flather 1
              qel(i)=qel(i)-thetai*dt*ri3
            endif

!	    I_5
            if(isbs(isd)>0.and.idry_s(isd)==0) then
              Unbar=bigu(isd,1)*snx(isd)+bigu(isd,2)*sny(isd)
              qel(i)=qel(i)-(1-thetai)*dt*distj(isd)*Unbar/2
            endif
          enddo !l=1,2 sides
        endif !bnd node i
      enddo !i=1,np

!     Check symmetry (comment out afterwards)
      do i=1,np
        do j=1,nnp(i)
          nd=inp(i,j)
          index=0
          do l=1,nnp(nd)
            if(inp(nd,l)==i) index=l
          enddo !l
          if(index==0) then
            write(11,*)'Index not symmetric:',i,j,nd
            stop
          endif
          if(dabs(sparsem(i,j)-sparsem(nd,index))>1.e-5) then
            write(11,*)'Matrix not symmetric:',i,j,nd,sparsem(i,j),sparsem(nd,index)
            stop
          endif
        enddo !j
      enddo !i


!...  To impose elevation essential b.c., create a mapping between element index and actual eq. index 
      icount=0
      do i=1,np
        if(isbnd(i)>0.and.iettype(isbnd(i))/=0) then
          icount=icount+1
        else
          imap(i)=i-icount
        endif
      enddo !i

!...  Assemble sparse matrix format
!...
      neq=0 !final index of eqs.
      nnz=0 !!# of non-zero entries
      do i=1,np
        if(isbnd(i)>0.and.iettype(isbnd(i))/=0) cycle

        neq=neq+1
        nnz=nnz+1
        icoef(neq)=nnz
        jcoef(nnz)=imap(i)
        if(imap(i)/=neq) then
          write(11,*)'Impossible 300',i
          stop
        endif
        e2coef(nnz)=sparsem(i,0)
        qel2(neq)=qel(i)
        eta3(neq)=eta2(i) !initial guess
        do j=1,nnp(i)
          nd=inp(i,j)
          if(isbnd(nd)>0.and.iettype(isbnd(nd))/=0) then !essential b.c.
            if(elbc(nd)<-9998) then
              write(11,*)'Eta not assigned:',nd,isbnd(nd)
              stop
            endif
            qel2(neq)=qel2(neq)-sparsem(i,j)*elbc(nd)
          else if(nd>i) then !upper triangle only
            nnz=nnz+1
            jcoef(nnz)=imap(nd)
            e2coef(nnz)=sparsem(i,j)
          endif
        enddo !j
      enddo !i=1,np

      if(neq<=0) then
        write(11,*)'No eqs. to be solved!'
        stop
      endif
      icoef(neq+1)=nnz+1

!...  solve the wave equation for elevations at each element center
!...  jcg jacobi conjugate gradient solver from itpack2d, srcv2d.f
!...
      st=0 !secnds(0.0) !timing the process

!...  input information about solver
!...
      call dfault(iparm,rparm)
      iparm(1)=itmax1
      iparm(2)=1 !level of output msg
      iparm(4)=33 !output msg to fort.??
      iparm(5)=0 !symmetric system
!      iparm(10)=iremove
      iparm(11)=1 !no timing
      iparm(12)=1 !error analysis
      rparm(1)=zeta  !0.11102230E-13 !stopping criterion
      rparm(8)=tol  !1.0e-8

!jcg    call nspcg (jac1,basic,ndim,mdim,nor,4,
!jcg &        e2coef,jcoef,p,ip,eta2,
!jcg &        ubar,qel,wksp,iwksp,nwksp,inw,iparm,rparm,ier)

      if(isolver==1) then
        call jcg(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else if(isolver==2) then
        call jsi(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else if(isolver==3) then
        call ssorcg(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else !isolver=4
        call ssorsi(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      endif

!...  Save eta1
      do i=1,np
        eta1(i)=eta2(i)
      enddo

!...  Re-assemble new elevations
!...
      etatot=0
      do i=1,np
        if(isbnd(i)>0.and.iettype(isbnd(i))/=0) then
          eta2(i)=elbc(i)
        else
          eta2(i)=eta3(imap(i))
        endif
        etatot=etatot+dabs(eta2(i))
        if(eta2(i)>elevmax(i)) elevmax(i)=eta2(i)
      enddo !i

      if(nscreen.eq.1) write(*,*) 'etatot=',etatot
      write(16,*) 'etatot=',etatot

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'solver took',btimer,'seconds...',iparm(1),iparm(8)
      write(16,*)'solver took',btimer,'seconds...',iparm(1),iparm(8)

!
!************************************************************************
!									*
!		Momentum equations					*
!									*
!************************************************************************
!

      call system_clock(ist,icount_rate)

!...  Along each side
!...  (su2,sv2) are in global frame
      do j=1,ns
        if(idry_s(j)==1) then
          do k=1,nvrt
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
          cycle
        endif

!	Wet sides
        node1=isidenode(j,1)
        node2=isidenode(j,2)

!       Define layer thickness & diffusivities
        do k=kbs(j)+1,nvrt
          dzz(k)=zs(k,j)-zs(k-1,j)
          dfz(k)=(dfv(node1,k)+dfv(node2,k)+dfv(node1,k-1)+dfv(node2,k-1))/4
        enddo !k

!	Coefficient matrix
        ndim=nvrt-kbs(j)
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j) !eq. #
          alow(kin)=0 
          cupp(kin)=0
          bdia(kin)=0
          if(k<nvrt) then
            tmp=dt*dfz(k+1)/dzz(k+1)
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
          endif

          if(k>kbs(j)+1) then
            tmp=dt*dfz(k)/dzz(k)
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
          else !b.c.
            bdia(kin)=bdia(kin)+dt*chi(j)
          endif
        enddo !k

!	RHS 
!	Elevation gradient, atmo. pressure and earth tidal potential
        deta2dx=0
        deta2dy=0
        deta1dx=0
        deta1dy=0
        icount1=0
        icount2=0
        dprdx=0
        dprdy=0
        detpdx=0
        detpdy=0
        icount3=0
        do l=1,2
          ie=is(j,l)
          if(ie/=0) then
            itmp=0
            do m=1,3
             nd=nm(ie,m)
             if(eta2(nd)+dp(nd)<=h0) itmp=1
            enddo !m
            if(itmp==0) then !wet
              icount2=icount2+1
              do m=1,3
                deta2dx=deta2dx+eta2(nm(ie,m))*dl(ie,m,1)
                deta2dy=deta2dy+eta2(nm(ie,m))*dl(ie,m,2)
              enddo !m
            endif
            if(idry_e(ie)==0) then
              icount1=icount1+1
              if(dpe(ie)>=tip_dp) icount3=icount3+1
              do m=1,3
                deta1dx=deta1dx+eta1(nm(ie,m))*dl(ie,m,1)
                deta1dy=deta1dy+eta1(nm(ie,m))*dl(ie,m,2)
                dprdx=dprdx+pr(nm(ie,m))*dl(ie,m,1)
                dprdy=dprdy+pr(nm(ie,m))*dl(ie,m,2)
                if(dpe(ie)>=tip_dp) then
                  detpdx=detpdx+etp(nm(ie,m))*dl(ie,m,1)
                  detpdy=detpdy+etp(nm(ie,m))*dl(ie,m,2)
                endif
              enddo !m
            endif
          endif !ie/=0
        enddo !l
        if(icount1/=0) then
          deta1dx=deta1dx/icount1
          deta1dy=deta1dy/icount1
          dprdx=dprdx/icount1
          dprdy=dprdy/icount1
        endif
        if(icount3/=0) then
          detpdx=detpdx/icount3
          detpdy=detpdy/icount3
        endif
        if(icount2/=0) then
          deta2dx=deta2dx/icount2
          deta2dy=deta2dy/icount2
        endif

!	b.c. to be imposed at the end
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
          rrhs(kin,1)=0
          rrhs(kin,2)=0
!	  Elevation gradient, atmo. pressure and tidal potential
          if(k<nvrt) then
            rrhs(kin,1)=rrhs(kin,1)-dzz(k+1)/2*dt*(g*thetai*deta2dx+g*(1-thetai)*deta1dx+dprdx/rho0-0.69*g*detpdx)
            rrhs(kin,2)=rrhs(kin,2)-dzz(k+1)/2*dt*(g*thetai*deta2dy+g*(1-thetai)*deta1dy+dprdy/rho0-0.69*g*detpdy)
          endif
          if(k>kbs(j)+1) then 
            rrhs(kin,1)=rrhs(kin,1)-dzz(k)/2*dt*(g*thetai*deta2dx+g*(1-thetai)*deta1dx+dprdx/rho0-0.69*g*detpdx)
            rrhs(kin,2)=rrhs(kin,2)-dzz(k)/2*dt*(g*thetai*deta2dy+g*(1-thetai)*deta1dy+dprdy/rho0-0.69*g*detpdy)
          endif

!	  Coriolis, advection, wind stress, and horizontal viscosity
          if(k<nvrt) then
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*sdbt(j,k,1)+sdbt(j,k+1,1)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k+1,j))+dt*(2*d2u(k,j)+d2u(k+1,j)))
            rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/6*(2*sdbt(j,k,2)+sdbt(j,k+1,2)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k+1,j))+dt*(2*d2v(k,j)+d2v(k+1,j)))
          else !k=nvrt
            rrhs(kin,1)=rrhs(kin,1)+dt*(tau(node1,1)+tau(node2,1))/2
            rrhs(kin,2)=rrhs(kin,2)+dt*(tau(node1,2)+tau(node2,2))/2
          endif

          if(k>kbs(j)+1) then
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*sdbt(j,k,1)+sdbt(j,k-1,1)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k-1,j))+dt*(2*d2u(k,j)+d2u(k-1,j)))
            rrhs(kin,2)=rrhs(kin,2)+dzz(k)/6*(2*sdbt(j,k,2)+sdbt(j,k-1,2)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k-1,j))+dt*(2*d2v(k,j)+d2v(k-1,j)))
          endif 

!         Baroclinic
          if(ibc==0) then
            if(k<nvrt) then
               rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*dt*(2*bcc(j,k,1)+bcc(j,k+1,1))
               rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/6*dt*(2*bcc(j,k,2)+bcc(j,k+1,2))
            endif
            if(k>kbs(j)+1) then
               rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*dt*(2*bcc(j,k,1)+bcc(j,k-1,1))
               rrhs(kin,2)=rrhs(kin,2)+dzz(k)/6*dt*(2*bcc(j,k,2)+bcc(j,k-1,2))
            endif
          endif !ibc==0
        enddo !k

        call tridag(mnv,ndim,2,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
!         Impose limits
          su2(k,j)=dmax1(-rmaxvel,dmin1(rmaxvel,soln(kin,1)))
          sv2(k,j)=dmax1(-rmaxvel,dmin1(rmaxvel,soln(kin,2)))
        enddo !k
        if(Cd(j)==0) then
          su2(kbs(j),j)=su2(kbs(j)+1,j)
          sv2(kbs(j),j)=sv2(kbs(j)+1,j)
        else !no slip bottom
          su2(kbs(j),j)=0
          sv2(kbs(j),j)=0
        endif

!       Extend
        khh=0 !larger of the 2 element bottom indices
        do l=1,2 !element
          ie=is(j,l)
          if(ie/=0.and.idry_e(ie)==0.and.kbe(ie)>khh) khh=kbe(ie)
        enddo !l
        if(khh==0) then
          write(11,*)'Cannot find the higher bottom:',j,(is(j,l),l=1,2)
          stop
        endif
        if(kbs(j)>khh) then
          write(11,*)'Side index > elemnt:',kbs(j),khh
          stop
        endif

        do k=1,khh-1
          su2(k,j)=0 !su2(kbs(j),j)
          sv2(k,j)=0 !sv2(kbs(j),j)
        enddo !k

!	Impose b.c.
        do k=1,nvrt
          if(isbs(j)>0.and.ifltype(isbs(j))/=0) then !open bnd side
!            ibnd=isbs(j)
            if(uth(j,k)<-98.or.vth(j,k)<-98) then
              write(11,*)'Wrong vel. input:',uth(j,k),vth(j,k),node1,node2
              stop
            endif
            if(ifltype(isbs(j))==-1) then !Flather 1
              if(eta_mean(node1)<-98.or.eta_mean(node2)<-98) then
                write(11,*)'Flather bnd elevation not assigned:',isbs(j)
                stop
              endif
              if(dps(j)<=0) then
                write(11,*)'Flather bnd has negative depth:',isbs(j),dps(j)
                stop
              endif
              vnorm=uth(j,k)*snx(j)+vth(j,k)*sny(j)+dsqrt(g/dps(j))* &
     &(eta2(node1)+eta2(node2)-eta_mean(node1)-eta_mean(node2))/2
              su2(k,j)=vnorm*snx(j)
              sv2(k,j)=vnorm*sny(j)
            else !not Flather
              su2(k,j)=uth(j,k)
              sv2(k,j)=vth(j,k)
            endif !Flather or not
          endif !open bnd

          if(isbs(j)==0.and.is(j,2)==0) then !land bnd
            if(islip==0) then !free slip
              vtan=-su2(k,j)*sny(j)+sv2(k,j)*snx(j)
              su2(k,j)=-vtan*sny(j)
              sv2(k,j)=vtan*snx(j)
            else !no slip
              su2(k,j)=0
              sv2(k,j)=0
            endif
          endif
        enddo !k
      enddo !j=1,ns

!...  Shapiro filter (used only if indvel=0)
!     use bcc as temporary variable
      if(indvel==0) then
        bcc=0
        do i=1,ns
          if(is(i,2)==0.or.idry_s(i)==1) cycle

!         Internal wet sides
          do k=kbs(i)+1,nvrt
            suru=0
            surv=0
            do j=1,4
              id=isidenei2(i,j)
              if(idry_s(id)==1) then
                kin=k
              else
                kin=max(k,kbs(id)+1)
              endif
              suru=suru+su2(kin,id)
              surv=surv+sv2(kin,id)
            enddo !j
            bcc(i,k,1)=su2(k,i)+shapiro/4*(suru-4*su2(k,i))
            bcc(i,k,2)=sv2(k,i)+shapiro/4*(surv-4*sv2(k,i))
          enddo !k
        enddo !i

        do j=1,ns
          if(is(j,2)==0.or.idry_s(j)==1) cycle

          do k=kbs(j)+1,nvrt
            su2(k,j)=bcc(j,k,1)
            sv2(k,j)=bcc(j,k,2)
          enddo !k

!         Extend
          khh=0 !larger of the 2 element bottom indices
          do l=1,2 !element
            ie=is(j,l)
            if(ie/=0.and.idry_e(ie)==0.and.kbe(ie)>khh) khh=kbe(ie)
          enddo !l
          if(khh==0) then
            write(11,*)'Cannot find the higher bottom (2):',j,(is(j,l),l=1,2)
            stop
          endif
          if(kbs(j)>khh) then
            write(11,*)'Side index > elemnt (2):',kbs(j),khh
            stop
          endif

          do k=1,khh-1
            su2(k,j)=0 !su2(kbs(j),j)
            sv2(k,j)=0 !sv2(kbs(j),j)
          enddo !k
        enddo !i
      endif !indvel=0; Shapiro filter

      if(nscreen.eq.1) write(*,*)'done solving momentum eq...'
      write(16,*)'done solving momentum eq...'

!...  solve for vertical velocities
!...
      we=0 !for dry and below bottom levels
      do i=1,ne
        if(idry_e(i)==1) cycle

!	Wet elements with 3 wet nodes
!	Compute upward normals and areas @ all levels
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        av_bdef1=(bdef1(n1)+bdef1(n2)+bdef1(n3))/3 !average bed deformation
        av_bdef2=(bdef2(n1)+bdef2(n2)+bdef2(n3))/3
        if(kbe(i)==0) then
          write(11,*)'Impossible 95'
          stop
        endif
        do l=kbe(i),nvrt
          xcon=(y(n2)-y(n1))*(z(l,n3)-z(l,n1))-(y(n3)-y(n1))*(z(l,n2)-z(l,n1))
          ycon=(x(n3)-x(n1))*(z(l,n2)-z(l,n1))-(x(n2)-x(n1))*(z(l,n3)-z(l,n1))
          zcon=area(i)*2
          area_e(l)=dsqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0) then
            write(11,*)'Zero area:',i,l
            stop
          endif
          sne(l,1)=xcon/area_e(l)/2
          sne(l,2)=ycon/area_e(l)/2
          sne(l,3)=zcon/area_e(l)/2 !>0
        enddo !l

!       Bottom b.c.
!       Error: only we=0
!        dhdx=dp(n1)*dl(i,1,1)+dp(n2)*dl(i,2,1)+dp(n3)*dl(i,3,1)
!        dhdy=dp(n1)*dl(i,1,2)+dp(n2)*dl(i,2,2)+dp(n3)*dl(i,3,2)
!        ubar=(su2(1,js(i,1))+su2(1,js(i,2))+su2(1,js(i,3)))/3
!        vbar=(sv2(1,js(i,1))+sv2(1,js(i,2))+sv2(1,js(i,3)))/3
        we(kbe(i),i)=(av_bdef2-av_bdef1)/dt !-dhdx*ubar-dhdy*vbar

        do l=kbe(i),nvrt-1
          sum=0
          ubar=0
          vbar=0
          ubar1=0
          vbar1=0
          do j=1,3
            jsj=js(i,j)
            vnor1=su2(l,jsj)*snx(jsj)+sv2(l,jsj)*sny(jsj)
            vnor2=su2(l+1,jsj)*snx(jsj)+sv2(l+1,jsj)*sny(jsj)
            if(l<kbs(jsj).or.kbs(jsj)==0) then
              write(11,*)'Impossible 94'
              stop
            endif
            sum=sum+ssign(i,j)*(zs(l+1,jsj)-zs(l,jsj))*distj(jsj)*(vnor1+vnor2)/2

            ubar=ubar+su2(l,jsj)/3    
            ubar1=ubar1+su2(l+1,jsj)/3    
            vbar=vbar+sv2(l,jsj)/3    
            vbar1=vbar1+sv2(l+1,jsj)/3    
          enddo !j=1,3

!         Impose bottom no-flux b.c.
          if(l==kbe(i)) then
            bflux=we(kbe(i),i)
          else
            bflux=ubar*sne(l,1)+vbar*sne(l,2)+we(l,i)*sne(l,3)
          endif

          we(l+1,i)=(-sum-(ubar1*sne(l+1,1)+vbar1*sne(l+1,2))*area_e(l+1)+ &
     &bflux*area_e(l))/sne(l+1,3)/area_e(l+1)

!         Debug
!          tmp1=sum
!          tmp2=(ubar1*sne(l+1,1)+vbar1*sne(l+1,2)+we(l+1,i)*sne(l+1,3))*area_e(l+1)-bflux*area_e(l)
!          if(i==24044.and.it==2) write(97,*)l,tmp1,tmp2,tmp1+tmp2

        enddo !l
      enddo !i=1,ne

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'done solving w; momentum took ',btimer,' seconds'
      write(16,*)'done solving w; momentum took ',btimer,' seconds'

!     Debug: test transport alone
!      eta1=0; eta2=0; we=0
!      do i=1,ns
!        do k=1,nvrt
!          su2(k,i)=-ycj(i)*2*pi/3.0e3
!          sv2(k,i)=xcj(i)*2*pi/3.0e3
!        enddo !k
!      enddo !i
!      do i=1,ne
!        do k=1,nvrt
!          do j=1,3
!            nd=nm(i,j)
!            ufg(k,i,j)=-y(nd)*2*pi/3.0e3
!            vfg(k,i,j)=x(nd)*2*pi/3.0e3
!          enddo !j
!        enddo !k
!      enddo !i
!      do i=1,np
!        do k=1,nvrt
!          uu2(k,i)=-y(i)*2*pi/3.e3 
!          vv2(k,i)=x(i)*2*pi/3.e3
!          ww2(k,i)=0 !-1.e-4*z(k,i)*(50+z(k,i))
!        enddo !k
!      enddo !i
!     End debug

      call system_clock(ist,icount_rate)
      if(ibc==0.or.ibtp==1) then
!----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
! 			Transport equation						!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!...  Initialize S,T as flags
      tnd=-99; snd=-99; tsd=-99; ssd=-99 !flags

!*************************************************************************************
!        ELM option
!*************************************************************************************

      if(iupwind_t==0.or.iupwind_s==0) then
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!...  Along each node & side
!...      
      do l=1,2 !nodes&sides
        if(l==1) then
          limit=np
        else
          limit=ns
        endif
        do i=1,limit

          if(l==1.and.idry(i)==1.or.l==2.and.idry_s(i)==1) cycle

!         Define nodes/sides, and layer thickness & diffusivities
          if(l==1) then !nodes
            nd0=i
            kbb=kbp(i)
            depth=dp(nd0)
          else !sides
            isd0=i
            kbb=kbs(i)
            node1=isidenode(isd0,1)
            node2=isidenode(isd0,2)
            depth=dps(isd0)
          endif
!          dzmin=1.e25
          do k=kbb+1,nvrt
            if(l==1) then !nodes
              dzz(k)=z(k,nd0)-z(k-1,nd0)
              dfz(k)=(dfh(nd0,k)+dfh(nd0,k-1))/2
            else !sides
              dzz(k)=zs(k,isd0)-zs(k-1,isd0)
              dfz(k)=(dfh(node1,k)+dfh(node2,k)+dfh(node1,k-1)+dfh(node2,k-1))/4
            endif
!            if(dzz(k)<dzmin) dzmin=dzz(k)
          enddo !k
!          do k=1,nvrt
!            if(k==nvrt) then
!              dz2(k)=dzz(k)/2
!            else if(k==1) then
!              dz2(k)=dzz(k+1)/2
!            else
!              dz2(k)=(dzz(k+1)+dzz(k))/2
!            endif
!            if(dz2(k)<dzmin) dzmin=dz2(k)
!          enddo !k

!	  Coefficient matrix
          ndim=nvrt-kbb+1
          do k=kbb,nvrt
            kin=k-kbb+1
            alow(kin)=0
            cupp(kin)=0
            bdia(kin)=0
            if(k<nvrt) then
              tmp=dt*dfz(k+1)/dzz(k+1)
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+dzz(k+1)/2+tmp
            endif

            if(k>kbb) then
              tmp=dt*dfz(k)/dzz(k)
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+dzz(k)/2+tmp
            endif
          enddo !k

!	  RHS 
!	  b.c. to be imposed at the end
          do k=kbb,nvrt
            kin=k-kbb+1
            rrhs(kin,1)=0
            rrhs(kin,2)=0
            if(k<nvrt) then 
              if(l==1) then
                rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/2*ptbt(nd0,k,3)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/2*ptbt(nd0,k,4)
              else
                rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/2*sdbt(isd0,k,3)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/2*sdbt(isd0,k,4)
              endif
!            else if(ihconsv/=0) then !surface flux including solar
!              if(l==1) then
!                rrhs(k,1)=rrhs(k,1)+dt/rho0/shw*(sflux(nd0)+srad(nd0))
!              else
!                rrhs(k,1)=rrhs(k,1)+dt/rho0/shw*(sflux(node1)+sflux(node2)+srad(node1)+srad(node2))/2
!              endif
            endif

            if(k>kbb) then 
              if(l==1) then
                rrhs(kin,1)=rrhs(kin,1)+dzz(k)/2*ptbt(nd0,k,3)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k)/2*ptbt(nd0,k,4)
              else
                rrhs(kin,1)=rrhs(kin,1)+dzz(k)/2*sdbt(isd0,k,3)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k)/2*sdbt(isd0,k,4)
              endif
!            else if(ihconsv/=0) then !bottom solar
!              if(l==1) then
!                htot=dmin1(1.d2,eta1(nd0)+dp(nd0)) !to prevent underflow
!                if(htot<=h0) then
!                  write(11,*)'Total depth (1) < 0:',htot
!                  stop
!                endif
!                rrhs(k,1)=rrhs(k,1)-dt/rho0/shw*srad(nd0)*(0.8*dexp(-htot/0.9)+0.2*dexp(-htot/2.1))
!              else
!                htot=dmin1(1.d2,(eta1(node1)+eta1(node2))/2+dps(isd0))
!                if(htot<=h0) then
!                  write(11,*)'Total depth (2) < 0:',htot
!                  stop
!                endif
!                rrhs(k,1)=rrhs(k,1)-dt/rho0/shw*(srad(node1)+srad(node2))/2*(0.8*dexp(-htot/0.9)+0.2*dexp(-htot/2.1))
!              endif
            endif

!           Compute fucntion F() for heat exchange
!            if(ihconsv/=0) then !solar
!              if(l==1) then
!                zz1=dmax1(-1.d2,z(k,nd0)-eta1(nd0))
!                if(zz1>0) then
!                  write(11,*)'Above f.s.:',zz1
!                  stop
!                endif
!                rrhs(k,5)=srad(nd0)*(0.8*0.9*dexp(zz1/0.9)+0.2*2.1*dexp(zz1/2.1)) !F()
!              else
!                zz1=dmax1(-1.d2,zs(k,isd0)-(eta1(node1)+eta1(node2))/2)
!!               Error: turn ifort12(16) back on after debugging 
!                if(zz1>0) then !.and.ifort12(16)==0) then
!                  ifort12(16)=1
!                  write(12,*)'Above f.s. (2):',l,i,k,zz1
!                endif
!                zz1=dmin1(0.d0,zz1)
!                rrhs(k,5)=(srad(node1)+srad(node2))/2*(0.8*0.9*dexp(zz1/0.9)+0.2*2.1*dexp(zz1/2.1))
!              endif
!            endif !solar
          enddo !k=1,nvrt

!         Solar
!          if(ihconsv/=0) then
!            do k=1,nvrt
!              if(k<nvrt) rrhs(k,1)=rrhs(k,1)+dt/rho0/shw*(rrhs(k+1,5)-rrhs(k,5))/dzz(k+1)
!              if(k>1) rrhs(k,1)=rrhs(k,1)-dt/rho0/shw*(rrhs(k,5)-rrhs(k-1,5))/dzz(k)
!            enddo !k
!          endif

          call tridag(mnv,ndim,2,alow,bdia,cupp,rrhs,soln,gam)

!         Impose no flux condition at bottom B.L. for slipless bottom
          if(l==1.and.Cdp(nd0)/=0.or.l==2.and.Cd(isd0)/=0) soln(1,1:2)=soln(2,1:2)

!         Correct overshoots for S,T
!         Debug
!          if(l==1.and.i==23) then
!            write(98,*)totalflux
!            do k=1,nvrt
!              write(98,*)k,soln(k,1),ptbt(nd0,k,3),rrhs(k,5),rrhs(k,4)*dt
!            enddo
!          endif

          tmin=100; tmax=-100
          smin=100; smax=-100
          do k=kbb,nvrt
            if(l==1) then
              if(tmin>ptbt(nd0,k,3)) tmin=ptbt(nd0,k,3)
              if(tmax<ptbt(nd0,k,3)) tmax=ptbt(nd0,k,3)
              if(smin>ptbt(nd0,k,4)) smin=ptbt(nd0,k,4)
              if(smax<ptbt(nd0,k,4)) smax=ptbt(nd0,k,4)
              rrhs(k,3)=ptbt(nd0,k,3) !for debugging only
            else
              if(tmin>sdbt(isd0,k,3)) tmin=sdbt(isd0,k,3)
              if(tmax<sdbt(isd0,k,3)) tmax=sdbt(isd0,k,3)
              if(smin>sdbt(isd0,k,4)) smin=sdbt(isd0,k,4)
              if(smax<sdbt(isd0,k,4)) smax=sdbt(isd0,k,4)
              rrhs(k,3)=sdbt(isd0,k,3)
            endif
          enddo !k
!         Reset extrema for heat exchange
!          if(ihconsv/=0) then
!            tmin=tempmin
!            tmax=tempmax
!          endif

!          if(tmin>tmax.or.tmin<0.or.smin>smax.or.smin<0) then
!            write(11,*)'Illegal min/max:',tmin,tmax,smin,smax,l,i
!            stop
!          endif

!         Store S,T in swild2 temporarily
          do k=kbb,nvrt
            kin=k-kbb+1
            swild2(k,1:2)=soln(kin,1:2)
!            if(ihconsv/=0) swild2(k,1)=dmax1(tempmin,dmin1(tempmax,soln(kin,1)))
!            if(isconsv/=0) swild2(k,2)=dmax1(saltmin,dmin1(saltmax,soln(kin,2)))
          enddo !k

!	  Impose b.c. on bnd nodes & sides
          if(l==1) then !nodes
            if(isbnd(nd0)>0) then 
              ibnd=isbnd(nd0)
              ind=0
              do ll=1,nond(ibnd)
                ndo=iond(ibnd,ll)
                if(ndo==nd0) then
                  ind=ll
                  exit
                endif
              enddo !ll
              if(ind==0) then
                write(11,*)'Impossible 101'
                stop
              endif

              isd=0 !flag
              do ll=1,2 !side
                if(ll==1) then
                  ie=ine(nd0,1)
                  id=iself(nd0,1)
                  isd3=js(ie,nx(id,2))
                else
                  ie=ine(nd0,nne(nd0))
                  id=iself(nd0,nne(nd0))
                  isd3=js(ie,nx(id,1))
                endif
                if(isbs(isd3)==ibnd) then
                  isd=isd3
                  exit
                endif
              enddo !ll
              if(isd==0) then
                write(11,*)'Cannot find an open side:',nd0
                stop
              endif
              do k=1,nvrt
                if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  swild2(k,1)=tth(ibnd,1,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=tem0(k,nd0)
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=tth(ibnd,ind,k)
                else if(itetype(ibnd)==-4) then
                  vnn=su2(k,isd)*snx(isd)+sv2(k,isd)*sny(isd)
                  if(vnn<0) swild2(k,1)=tobc(ibnd)*tth(ibnd,ind,k)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==-1) then
                  vnn=su2(k,isd)*snx(isd)+sv2(k,isd)*sny(isd)
                  if(vnn<0) swild2(k,1)=tobc(ibnd)*tem0(k,nd0)+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sth(ibnd,1,1)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=sal0(k,nd0)
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=sth(ibnd,ind,k)
                else if(isatype(ibnd)==-4) then
                  vnn=su2(k,isd)*snx(isd)+sv2(k,isd)*sny(isd)
                  if(vnn<0) swild2(k,2)=sobc(ibnd)*sth(ibnd,ind,k)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==-1) then
                  vnn=su2(k,isd)*snx(isd)+sv2(k,isd)*sny(isd)
                  if(vnn<0) swild2(k,2)=sobc(ibnd)*sal0(k,nd0)+(1-sobc(ibnd))*swild2(k,2)
                endif
              enddo !k
            endif !isbnd>0
          else !sides
            if(isbs(isd0)>0) then
              ibnd=isbs(isd0)
              in1=0; in2=0
              do ll=1,nond(ibnd)
                ndo=iond(ibnd,ll)
                if(ndo==node1) in1=ll
                if(ndo==node2) in2=ll
                if(in1/=0.and.in2/=0) exit
              enddo !ll
              if(in1==0.or.in2==0) then
                write(11,*)'Impossible 102'
                stop
              endif

              do k=1,nvrt
                if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  swild2(k,1)=tth(ibnd,1,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=(tem0(k,node1)+tem0(k,node2))/2
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=(tth(ibnd,in1,k)+tth(ibnd,in2,k))/2
                else if(itetype(ibnd)==-4) then
                  vnn=su2(k,isd0)*snx(isd0)+sv2(k,isd0)*sny(isd0)
                  if(vnn<0) swild2(k,1)=tobc(ibnd)*(tth(ibnd,in1,k)+tth(ibnd,in2,k))/2+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==-1) then
                  vnn=su2(k,isd0)*snx(isd0)+sv2(k,isd0)*sny(isd0)
                  if(vnn<0) swild2(k,1)=tobc(ibnd)*(tem0(k,node1)+tem0(k,node2))/2+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sth(ibnd,1,1)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=(sal0(k,node1)+sal0(k,node2))/2
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=(sth(ibnd,in1,k)+sth(ibnd,in2,k))/2
                else if(isatype(ibnd)==-4) then
                  vnn=su2(k,isd0)*snx(isd0)+sv2(k,isd0)*sny(isd0)
                  if(vnn<0) swild2(k,2)=sobc(ibnd)*(sth(ibnd,in1,k)+sth(ibnd,in2,k))/2+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==-1) then
                  vnn=su2(k,isd0)*snx(isd0)+sv2(k,isd0)*sny(isd0)
                  if(vnn<0) swild2(k,2)=sobc(ibnd)*(sal0(k,node1)+sal0(k,node2))/2+(1-sobc(ibnd))*swild2(k,2)
                endif
              enddo !k
            endif !isbs(isd0)>0
          endif

!         Nudging
          if(inu_st/=0) then
            do k=kbb,nvrt
              if(l==1) then !nodes
                if(z(k,nd0)>=-vnh1) then
                  vnf=vnf1 !vertical nudging factor
                else if(z(k,nd0)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(z(k,nd0)+vnh1)/(-vnh2+vnh1) 
                else
                  vnf=vnf2
                endif
                tnu=t_nudge(nd0)+vnf
                snu=s_nudge(nd0)+vnf
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(11,*)'Nudging factor out of bound (1):',tnu,snu,vnf
                  stop
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tem0(k,nd0)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+sal0(k,nd0)*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tnd_nu(nd0,k)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+snd_nu(nd0,k)*snu
                endif
              else !sides
                if(zs(k,isd0)>=-vnh1) then
                  vnf=vnf1 !vertical nudging factor
                else if(zs(k,isd0)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(zs(k,isd0)+vnh1)/(-vnh2+vnh1) 
                else
                  vnf=vnf2
                endif
                tnu=(t_nudge(node1)+t_nudge(node2))/2+vnf
                snu=(s_nudge(node1)+s_nudge(node2))/2+vnf
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(11,*)'Nudging factor out of bound (2):',tnu,snu,vnf
                  stop
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tem0(k,node1)+tem0(k,node2))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(sal0(k,node1)+sal0(k,node2))/2*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tnd_nu(node1,k)+tnd_nu(node2,k))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(snd_nu(node1,k)+snd_nu(node2,k))/2*snu
                endif
              endif
            enddo !k
          endif !nudging

!         Extend
          do k=1,kbb-1
            swild2(k,1)=swild2(kbb,1)
            swild2(k,2)=swild2(kbb,2)
          enddo !k

!         Check bounds
!         Error: compute initial min. max, and check against them
          do k=1,nvrt
            if(swild2(k,1)<-98.or.swild2(k,2)<-98) then
              write(11,*)'Werid ST (1):',l,i,k,swild2(k,1),swild2(k,2)
              stop
            endif
          enddo !k

!         Output S,T
          do k=1,nvrt
            if(iupwind_t==0) then
              if(l==1) then
                tnd(k,nd0)=swild2(k,1)
              else
                tsd(k,isd0)=swild2(k,1)
              endif
            endif !iupwind_t

            if(iupwind_s==0) then
              if(l==1) then
                snd(k,nd0)=swild2(k,2)
              else
                ssd(k,isd0)=swild2(k,2)
              endif
            endif !iupwind_s
          enddo !k

        enddo !i=1,limit
      enddo !l=1,2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      endif
!     end of ELM option

!*************************************************************************************
!
!        Upwind and TVD option
!
!*************************************************************************************
      if(iupwind_t/=0.or.iupwind_s/=0) then
!       b.c. and body forces
        bdy_frc=0; flx_sf=0; flx_bt=0

!       Salt exchange
        if(isconsv/=0) then
          do i=1,ne
            if(idry_e(i)==1) cycle
       
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)
            evap=(fluxevp(n1)+fluxevp(n2)+fluxevp(n3))/3
            precip=(fluxprc(n1)+fluxprc(n2)+fluxprc(n3))/3
            flx_sf(i,2)=tsel(nvrt,i,2)*(evap-precip)/rho0
          enddo !i
        endif !isconsv/=0

!       Heat exchange
        if(ihconsv/=0) then
          do i=1,ne
            if(idry_e(i)==1) cycle

!           Wet element
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)

!           Surface flux
            sflux_e=(sflux(n1)+sflux(n2)+sflux(n3))/3
            flx_sf(i,1)=sflux_e/rho0/shw

!           Solar
            do k=kbe(i)+1,nvrt
              srad_e=(srad(n1)+srad(n2)+srad(n3))/3
!             Don't use eta2 as it has been updated but not z()
!              elv_e=(eta2(n1)+eta2(n2)+eta2(n3))/3
              dp1=dmin1(ze(nvrt,i)-ze(k-1,i),1.d2) !to prevent underflow
              dp2=dmin1(ze(nvrt,i)-ze(k,i),1.d2) !to prevent underflow
              if(dp2<0.or.dp2>dp1) then
                write(11,*)'Depth<0 in upwind transport:',i,k,dp1,dp2
                write(11,*)eta2(n1),eta2(n2),eta2(n3),ze(nvrt,i)
                do l=kbe(i),nvrt
                  write(11,*)l,z(l,n1),z(l,n2),z(l,n3)
                enddo !l
                stop
              endif

              if(k==kbe(i)+1) then
                srad1=0
              else
                srad1=srad_e*(0.8*dexp(-dp1/0.9)+0.2*dexp(-dp1/2.1))
              endif
              srad2=srad_e*(0.8*dexp(-dp2/0.9)+0.2*dexp(-dp2/2.1))
              if(srad2<srad1.and.ifort12(19)==0) then
                ifort12(19)=1
                write(12,*)'Reset negative solar hearting:',i,k,srad2,srad1,srad2-srad1
              endif
              bdy_frc(k,i,1)=dmax1(srad2-srad1,0.d0)/rho0/shw/(ze(k,i)-ze(k-1,i)) !Q
            enddo !k=kbe(i)+1,nvrt
          enddo !i=1,ne
        endif !heat exchange

        up_tvd=(iupwind_t==2.or.iupwind_s==2)
        tr_el(1:mnv,1:mne,1:2)=tsel(1:mnv,1:mne,1:2)
        call do_transport_tvd(0,up_tvd,tvd_mid,flimiter,2)
        tsel(1:mnv,1:mne,1:2)=tr_el(1:mnv,1:mne,1:2)

!       Impose b.c. & nudging
        do i=1,ne
          if(idry_e(i)==1) cycle

          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)

!         Impose b.c.
          ibnd=0 !flag
          do j=1,3
            isd=js(i,j)
            if(isbs(isd)>0) then
              ibnd=isbs(isd)
              isd00=isd
              exit !first open bnd counts
            endif
          enddo !j

          if(ibnd>0) then !open bnd
            node1=isidenode(isd00,1)
            node2=isidenode(isd00,2)
            ind1=0; ind2=0
            do ll=1,nond(ibnd)
              ndo=iond(ibnd,ll)
              if(ndo==node1) ind1=ll
              if(ndo==node2) ind2=ll
            enddo !ll
            if(ind1==0.or.ind2==0) then
              write(11,*)'Cannot find a local index'
              stop
            endif

            do k=kbe(i)+1,nvrt
              tsel01=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
              tsel02=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
              if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                tsel(k,i,1)=tth(ibnd,1,1)
              else if(itetype(ibnd)==3) then
                tsel(k,i,1)=tsel01
              else if(iabs(itetype(ibnd))==4) then
                tmp=(tth(ibnd,ind1,k)+tth(ibnd,ind1,k-1)+tth(ibnd,ind2,k)+tth(ibnd,ind2,k-1))/4
                if(itetype(ibnd)==4) then
                  tsel(k,i,1)=tmp
                else
                  vnn=su2(k,isd00)*snx(isd00)+sv2(k,isd00)*sny(isd00)
                  if(vnn<0) tsel(k,i,1)=tobc(ibnd)*tmp+(1-tobc(ibnd))*tsel(k,i,1)
                endif
              else if(itetype(ibnd)==-1) then
                vnn=su2(k,isd00)*snx(isd00)+sv2(k,isd00)*sny(isd00)
                if(vnn<0) tsel(k,i,1)=tobc(ibnd)*tsel01+(1-tobc(ibnd))*tsel(k,i,1)
              endif

              if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                tsel(k,i,2)=sth(ibnd,1,1)
              else if(isatype(ibnd)==3) then
                tsel(k,i,2)=tsel02
              else if(iabs(isatype(ibnd))==4) then
                tmp=(sth(ibnd,ind1,k)+sth(ibnd,ind1,k-1)+sth(ibnd,ind2,k)+sth(ibnd,ind2,k-1))/4
                if(isatype(ibnd)==4) then
                  tsel(k,i,2)=tmp
                else
                  vnn=su2(k,isd00)*snx(isd00)+sv2(k,isd00)*sny(isd00)
                  if(vnn<0) tsel(k,i,2)=sobc(ibnd)*tmp+(1-sobc(ibnd))*tsel(k,i,2)
                endif
              else if(isatype(ibnd)==-1) then
                vnn=su2(k,isd00)*snx(isd00)+sv2(k,isd00)*sny(isd00)
                if(vnn<0) tsel(k,i,2)=sobc(ibnd)*tsel02+(1-sobc(ibnd))*tsel(k,i,2)
              endif
            enddo !k
          endif !open bnd

          if(inu_st/=0) then
            do k=kbe(i)+1,nvrt
              if(ze(k,i)>=-vnh1) then
                vnf=vnf1 !vertical nudging factor
              else if(ze(k,i)>=-vnh2) then
                vnf=vnf1+(vnf2-vnf1)*(ze(k,i)+vnh1)/(-vnh2+vnh1)
              else
                vnf=vnf2
              endif
              tnu=(t_nudge(n1)+t_nudge(n2)+t_nudge(n3))/3+vnf
              snu=(s_nudge(n1)+s_nudge(n2)+s_nudge(n3))/3+vnf
              if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                write(11,*)'Nudging factor out of bound (1):',tnu,snu,vnf
                stop
              endif

              if(inu_st==1) then !to i.c.
                tsel01=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
                tsel02=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
                tsel(k,i,1)=tsel(k,i,1)*(1-tnu)+tsel01*tnu
                tsel(k,i,2)=tsel(k,i,2)*(1-snu)+tsel02*snu
              else if(inu_st==2) then
                tnd_nu_e=(tnd_nu(n1,k)+tnd_nu(n2,k)+tnd_nu(n3,k)+tnd_nu(n1,k-1)+tnd_nu(n2,k-1)+tnd_nu(n3,k-1))/6
                snd_nu_e=(snd_nu(n1,k)+snd_nu(n2,k)+snd_nu(n3,k)+snd_nu(n1,k-1)+snd_nu(n2,k-1)+snd_nu(n3,k-1))/6
                tsel(k,i,1)=tsel(k,i,1)*(1-tnu)+tnd_nu_e*tnu
                tsel(k,i,2)=tsel(k,i,2)*(1-snu)+snd_nu_e*snu
              endif
            enddo !k
          endif !inu_st/=0

!         Extend
          do k=1,kbe(i)
            tsel(k,i,1:2)=tsel(kbe(i)+1,i,1:2)
          enddo !k
        enddo !i=1,ne

!       Convert to S,T at nodes and sides and whole levels
!       Use hp_int to temporarily store values at elements and whole levels
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(11,*)'Ratio out of bound:',i,k,zrat
              stop
            endif
            hp_int(k,i,1:2)=(1-zrat)*tsel(k+1,i,1:2)+zrat*tsel(k,i,1:2)
          enddo !k
          hp_int(nvrt,i,1:2)=tsel(nvrt,i,1:2)
          hp_int(kbe(i),i,1:2)=tsel(kbe(i)+1,i,1:2)
        enddo !i=1,ne

        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            tt1=0; ss1=0
            ta=0
            do j=1,nne(i)
              ie=ine(i,j)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                tt1=tt1+hp_int(kin,ie,1)*area(ie)
                ss1=ss1+hp_int(kin,ie,2)*area(ie)
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet if and only if at least one surrounding element is wet
              write(11,*)'Isolated wet node (9):',i
              stop
            else
              if(iupwind_t/=0) tnd(k,i)=tt1/ta
              if(iupwind_s/=0) snd(k,i)=ss1/ta
            endif
          enddo !k
        enddo !i=1,np

        do i=1,ns
          if(idry_s(i)==1) cycle

          do k=1,nvrt
            tt1=0; ss1=0
            ta=0
            do j=1,2
              ie=is(i,j)
              if(ie/=0.and.idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                tt1=tt1+hp_int(kin,ie,1)*area(ie)
                ss1=ss1+hp_int(kin,ie,2)*area(ie)
              endif
            enddo !j
            if(ta==0) then 
              write(11,*)'Isolated wet side (9):',i,(is(i,j),j=1,2)
              stop
            else
              if(iupwind_t/=0) tsd(k,i)=tt1/ta
              if(iupwind_s/=0) ssd(k,i)=ss1/ta
            endif
          enddo !k
        enddo !i=1,ns

!...    Compute total mass of S,T
!...
        open(91,file='total_ST.dat')
        tot_heat=0
        tot_salt=0
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            tot_heat=tot_heat+vol*tsel(k,i,1)
            tot_salt=tot_salt+vol*tsel(k,i,2)
          enddo !k
        enddo !i=1,ne
        write(91,*)time/86400,tot_heat,tot_salt

      endif !upwind

      if(nscreen==1) write(*,*)'done solving transport equation'
      write(16,*)'done solving transport equation'

!----------------------------------------------------------------------
      endif !ibc.eq.0.or.ibtp.eq.1

!...  Tracer transport
      if(ntracers>0) then

!       user-defined tracer part
!       define bdy_frc, flx_sf, flx_bt
!       bdy_frc(kbe(i)+1:nvrt,1:ne,ntracers): body force at prism center Q_{i,k} (for all wet elements i)
!       flx_sf(mne,mntr): surface b.c. \kappa*dC/dz = flx_sf (at element center)
!       flx_bt(mne,mntr): bottom b.c.

!        bdy_frc=0; flx_sf=0; flx_bt=0
        rkk1=log(10.)/3/3600
        rkk2=log(10.)/30/3600
        flx_sf=0; flx_bt=0
        bdy_frc(1:nvrt,1:ne,1)=-rkk1*trel(1:nvrt,1:ne,1)
        bdy_frc(1:nvrt,1:ne,2)=-rkk2*trel(1:nvrt,1:ne,2)+0.5*rkk1*trel(1:nvrt,1:ne,1)
!       end user-defined tracer part

        up_tvd=itr_met==2
        tr_el(1:mnv,1:mne,1:ntracers)=trel(1:mnv,1:mne,1:ntracers)
        call do_transport_tvd(1,up_tvd,tvd_mid2,flimiter2,ntracers)
        trel(1:mnv,1:mne,1:ntracers)=tr_el(1:mnv,1:mne,1:ntracers)

!       Impose b.c. 
        do i=1,ne
          if(idry_e(i)==1) cycle

          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)

          ibnd=0 !flag
          do j=1,3
            isd=js(i,j)
            if(isbs(isd)>0) then
              ibnd=isbs(isd)
              isd00=isd
              exit !first open bnd counts
            endif
          enddo !j

          if(ibnd>0) then !open bnd
            node1=isidenode(isd00,1)
            node2=isidenode(isd00,2)
            ind1=0; ind2=0
            do ll=1,nond(ibnd)
              ndo=iond(ibnd,ll)
              if(ndo==node1) ind1=ll
              if(ndo==node2) ind2=ll
            enddo !ll
            if(ind1==0.or.ind2==0) then
              write(11,*)'Cannot find a local index'
              stop
            endif

            do k=kbe(i)+1,nvrt
              if(itrtype(ibnd)==2) then
                trel(k,i,1:ntracers)=trth(ibnd,1:ntracers)
              else if(itrtype(ibnd)==3) then
                trel(k,i,1:ntracers)=trel0(k,i,1:ntracers)
              endif
            enddo !k
          endif !open bnd

!         Extend
          do k=1,kbe(i)
            trel(k,i,1:ntracers)=trel(kbe(i)+1,i,1:ntracers)
          enddo !k
        enddo !i=1,ne

!       Convert to S,T at nodes and sides and whole levels
!       Use tr_el to temporarily store values at elements and whole levels
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(11,*)'Ratio out of bound:',i,k,zrat
              stop
            endif
            tr_el(k,i,1:ntracers)=(1-zrat)*trel(k+1,i,1:ntracers)+zrat*trel(k,i,1:ntracers)
          enddo !k
          tr_el(nvrt,i,1:ntracers)=trel(nvrt,i,1:ntracers)
          tr_el(kbe(i),i,1:ntracers)=trel(kbe(i)+1,i,1:ntracers)
        enddo !i=1,ne

        tr_nd=-99 !for dry nodes
        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            swild(1:ntracers)=0
            ta=0
            do j=1,nne(i)
              ie=ine(i,j)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                swild(1:ntracers)=swild(1:ntracers)+tr_el(kin,ie,1:ntracers)*area(ie)
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet if and only if at least one surrounding element is wet
              write(11,*)'Isolated wet node (9):',i
              stop
            else
              tr_nd(k,i,1:ntracers)=swild(1:ntracers)/ta
            endif
          enddo !k
        enddo !i=1,np

!...    Compute total mass 
!...
        open(92,file='total_TR.dat')
        swild(1:ntracers)=0
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            swild(1:ntracers)=swild(1:ntracers)+vol*trel(k,i,1:ntracers)
          enddo !k
        enddo !i=1,ne
        write(92,*)time/86400,swild(1:ntracers)
      endif
!...  End of tracer transport

!...  Update bed deformation and depth info
      do i=1,np
        bdef1(i)=bdef2(i)
        dp(i)=dp00(i)-bdef1(i)
        hmod(i)=dmin1(dp(i),h_s)
      enddo !i

      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
      do i=1,ne
        dpe(i)=1.e10
        do j=1,3
          if(dpe(i)>dp(nm(i,j))) dpe(i)=dp(nm(i,j))
        enddo !j
      enddo !i=1,ne

!...  Recompute vgrid and calculate rewetted pts
      if(inunfl==0) then
        call levels0(iths,it)
      else
        call levels1(iths,it)
      endif
      if(nscreen.eq.1) write(*,*) 'done recomputing levels...'
      write(16,*) 'done recomputing levels...'

!...  Compute nodal vel. for output and next backtracking
      call nodalvel(ifltype)

!...  Density (using new level indices)
      call eqstate

      call system_clock(ien,icount_rate)
      btimer=real(ien-ist)/icount_rate
      if(nscreen.eq.1) write(*,*)'transport eqs. took',btimer,'sec'
      write(16,*)'transport eqs. took',btimer,'sec'

!...  Compute depth averaged h-vel.
!...
      dav=0
      do i=1,np
        if(idry(i)==1) cycle
        do k=kbp(i),nvrt-1
          dav(i,1)=dav(i,1)+(uu2(k+1,i)+uu2(k,i))/2*(z(k+1,i)-z(k,i))
          dav(i,2)=dav(i,2)+(vv2(k+1,i)+vv2(k,i))/2*(z(k+1,i)-z(k,i))
        enddo !k
        htot=eta2(i)+dp(i)
        if(htot<=h0) then
          write(11,*)'Impossible 24'
          stop
        endif
        dav(i,1)=dav(i,1)/htot
        dav(i,2)=dav(i,2)/htot
      enddo !i=1,np

!...  Optional computation of fluxes and total volume etc.
      if(iflux/=0) then
!--------------------------------------------------
!     Compute total mass etc.
      tvol=0 !total volume
      tmass=0 !total mass
      tpe=0 !total potential energy
      tkne=0 !total kinetic energy (2D only)
      enerf=0 !energy loss due to bottom friction; only correct for 2D model
      ener_ob=0 !total wave enery out of open bnds; only correct for 0 mean flows!
      do i=1,ne
        if(idry_e(i)==1) cycle

        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        etam=(eta2(n1)+eta2(n2)+eta2(n3))/3
        tpe=tpe+0.5*rho0*g*area(i)*etam**2
        av_dep=etam+(dp(n1)+dp(n2)+dp(n3))/3
        tvol=tvol+area(i)*av_dep
        do k=kbe(i),nvrt-1
          ah=(z(k+1,n1)+z(k+1,n2)+z(k+1,n3)-z(k,n1)-z(k,n2)-z(k,n3))/3
        enddo !j

        do j=1,3
          nd=nm(i,j)
          do k=kbp(nd),nvrt-1
            tmass=tmass+area(i)*(prho(nd,k)+prho(nd,k+1))*(z(k+1,nd)-z(k,nd))/6
          enddo !k
          htot=eta2(nd)+dp(nd)
          if(htot<=h0) then
            write(11,*)'Impossible dry (9):',i,j,nd,htot
            stop
          endif

          isd=js(i,j)
          do k=kbs(isd),nvrt-1
            vmag1=su2(k,isd)**2+sv2(k,isd)**2
            vmag2=su2(k+1,isd)**2+sv2(k+1,isd)**2
            tkne=tkne+rho0*area(i)/6*(zs(k+1,isd)-zs(k,isd))*(vmag1+vmag2)/2
          enddo !k

!         enerf only correct for 2D model
          enerf=enerf+dt*area(i)/3*rho0*Cdp(nd)*dsqrt(dav(nd,1)**2+dav(nd,2)**2)**3

!         ener_ob
          isd=js(i,j)
          if(isbs(isd)>0) then !open bnd
            n1=isidenode(isd,1)
            n2=isidenode(isd,2)
            eta_m=(eta2(n1)+eta2(n2))/2
            vel_m1=(dav(n1,1)+dav(n2,1))/2
            vel_m2=(dav(n1,2)+dav(n2,2))/2
            ener_ob=ener_ob+rho0/2*dsqrt(g*dps(isd))*dt*(g*eta_m**2+dps(isd)*(vel_m1**2+vel_m2**2))*distj(isd)
          endif
        enddo !j=1,3

      enddo !i=1,ne
      write(10,*)time/3600,tvol,tmass,tpe,tkne,tpe+tkne,enerf,ener_ob

!     Fluxes
      open(13,file='fluxflag.gr3',status='old')
      read(13,*)
      read(13,*)ntmp,npflag
      if(npflag/=np) then
        write(11,*)'# of pts in fluxflag should = np'
        stop
      endif
      do i=1,np
        read(13,*)j,xtmp,ytmp,wild
        nwild2(i)=wild
      enddo
      close(13)
      do i=1,ne
        nwild(i)=0
        do j=1,3
          nd=nm(i,j)
          if(nwild(i)<nwild2(nd)) nwild(i)=nwild2(nd)
        enddo !j
      enddo !i

      open(9,file='flux.dat') !,status='append')
      tvol12=0 !total volume inside rgns 1 and 2
      fluxbnd=0 !total flux across natural bnds (positive outward)
      fluxchan=0  !total flux across unnatural bnds (fluxchan1+fluxchan2)
      fluxchan1=0  !flux at bnd 1 (positive outward)
      fluxchan2=0  !flux at bnd 2 (positive outward)
 
      tot_s=0 !total salt inside rgns 1 and 2; [PSU*m^3]
      flux_s=0 !flux out of  rgns 1 and 2 (positive out) [PSU*m^3/s]
      do i=1,ne
        if(nwild(i)==0) cycle

        do j=1,3 !nodes or sides
          if(idry_e(i)==0) then
            nd=nm(i,j)
            do k=kbp(nd),nvrt-1
              ah=z(k+1,nd)-z(k,nd)
              tvol12=tvol12+area(i)*ah/3
              tot_s=tot_s+(snd(k+1,nd)+snd(k,nd))/2*ah*area(i)/3
            enddo !k
          endif !wet element

          isd=js(i,j)
          if(idry_s(isd)==0) then
            do k=kbs(isd),nvrt-1
              if(ic3(i,j)==0) then !bnd side
                vnn=(su2(k+1,isd)+su2(k,isd))/2*snx(isd)+(sv2(k+1,isd)+sv2(k,isd))/2*sny(isd) 
                ftmp=ssign(i,j)*distj(isd)*(zs(k+1,isd)-zs(k,isd))*vnn
                fluxbnd=fluxbnd+ftmp
              else if(nwild(ic3(i,j))==0) then !channel side
                vnn=(su2(k+1,isd)+su2(k,isd))/2*snx(isd)+(sv2(k+1,isd)+sv2(k,isd))/2*sny(isd)
                ftmp=ssign(i,j)*distj(isd)*(zs(k+1,isd)-zs(k,isd))*vnn
                fluxchan=fluxchan+ftmp
                if(nwild(i)==1) then
                  fluxchan1=fluxchan1+ftmp
                else !nwild(i)=2
                  fluxchan2=fluxchan2+ftmp
                endif

                flux_s=flux_s+ftmp*(ssd(k+1,isd)+ssd(k,isd))/2 !*dt=total salt lost
              endif
            enddo !k
          endif !wet side
        enddo !j=1,3
      enddo !i=1,ne

      write(9,*)time/86400,fluxchan1,-fluxchan2,fluxbnd,tvol12,fluxchan,0.,0.
!      write(91,*)time/86400,tot_s,flux_s
      if(nscreen.eq.1) write(*,*)'done computing fluxes...'
      write(16,*)'done computing fluxes...'
!---------------------------------------------------------      
      endif !iflux ne 0
!...  end compute flux balance

!...  output global data 
!...

!     Debug
!      write(99,*)'Time step= ',it
!      do i=1,np
!        if(fluxprc(i)/=0) write(99,*)i,fluxprc(i)
!      enddo !i

      do j=1,noutput
        if(iof(j).eq.1.and.mod(it,nspool).eq.0) then
          if(iwrite==0) then
            write(ichan(j),rec=irec(j)+1) real(time)
            write(ichan(j),rec=irec(j)+2) it
            irec(j)=irec(j)+2
          else !evm
            a_4 = transfer(source=real(time),mold=a_4)
            write(ichan(j),"(a4)",advance="no") a_4
            a_4 = transfer(source=it,mold=a_4)
            write(ichan(j),"(a4)",advance="no") a_4
          endif
          do i=1,np
            if(iwrite.eq.0) then
              write(ichan(j),rec=irec(j)+1) real(eta2(i))
              irec(j)=irec(j)+1
            else !evm
              a_4 = transfer(source=real(eta2(i)),mold=a_4)
              write(ichan(j),"(a4)",advance="no") a_4
            endif
          enddo !i

          do i=1,np
            if(j<=12) then
              floatout=0 !in case some are not defined
              if(j.eq.1) then
                floatout=eta2(i)
              else if(j.eq.2.and.ihconsv.ne.0) then
                floatout=pr(i)
              else if(j.eq.3.and.ihconsv.ne.0) then
                floatout=airt1(i)
              else if(j.eq.4.and.ihconsv.ne.0) then
                floatout=shum1(i)
              else if(j.eq.5.and.ihconsv.ne.0) then
                floatout=srad(i)
              else if(j.eq.6.and.ihconsv.ne.0) then
                floatout=fluxsu(i)
              else if(j.eq.7.and.ihconsv.ne.0) then
                floatout=fluxlu(i)
              else if(j.eq.8.and.ihconsv.ne.0) then
                floatout=hradu(i)
              else if(j.eq.9.and.ihconsv.ne.0) then
                floatout=hradd(i)
              else if(j.eq.10.and.ihconsv.ne.0) then
                floatout=sflux(i)
              else if(j.eq.11.and.isconsv.ne.0) then
                floatout=fluxevp(i)
              else if(j.eq.12.and.isconsv.ne.0) then
                floatout=fluxprc(i)
              endif

              if(iwrite.eq.0) then
                write(ichan(j),rec=irec(j)+1) floatout
                irec(j)=irec(j)+1
              else !evm
                a_4 = transfer(source=floatout,mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
              endif
            else if(j<=15) then
              if(j==13) then
                if(nws==0) then
                  floatout=0
                  floatout2=0
                else
                  floatout=windx(i)
                  floatout2=windy(i)
                endif
              else if(j==14) then
                floatout=tau(i,1)
                floatout2=tau(i,2)
              else !j=15
                floatout=dav(i,1)
                floatout2=dav(i,2)
              endif
              if(iwrite.eq.0) then
                write(ichan(j),rec=irec(j)+1) floatout
                write(ichan(j),rec=irec(j)+2) floatout2
                irec(j)=irec(j)+2
              else !evm
                a_4 = transfer(source=floatout,mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
                a_4 = transfer(source=floatout2,mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
              endif
            else if(j<25) then
              do k=max0(1,kbp00(i)),nvrt
                floatout=0 !for some undefined variables
                if(j.eq.16) then
                  floatout=ww2(k,i)
                else
                  if(j.eq.17) then
                    floatout=tnd(k,i)
                  else if(j.eq.18) then
                    floatout=snd(k,i)
                  else if(j.eq.19) then
                    floatout=prho(i,k)
                  else if(j.eq.20) then
                    floatout=dfh(i,k)
                  else if(j.eq.21) then
                    floatout=dfv(i,k)
                  else if(j.eq.22) then
                    floatout=q2(i,k)
                  else if(j.eq.23) then
                    floatout=xl(i,k)
                  else if(j==24) then
                    if(idry(i)==1) then
                      floatout=0
                    else
                      floatout=z(max0(k,kbp(i)),i)
                    endif
                  endif
                endif

                if(iwrite.eq.0) then
                  write(ichan(j),rec=irec(j)+1) floatout
                  irec(j)=irec(j)+1
                else !evm
                  a_4 = transfer(source=floatout,mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
                endif
              enddo !k
            else if(j==25) then
              do k=max0(1,kbp00(i)),nvrt
                if(iwrite.eq.0) then
                  write(ichan(j),rec=irec(j)+1) real(uu2(k,i))
                  write(ichan(j),rec=irec(j)+2) real(vv2(k,i))
                  irec(j)=irec(j)+2
                else !evm
                  a_4 = transfer(source=real(uu2(k,i)),mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
                  a_4 = transfer(source=real(vv2(k,i)),mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
                endif
              enddo !k

            else !tracers; implies ntracers>0
              do k=max0(1,kbp00(i)),nvrt
                floatout=tr_nd(k,i,j-25) !warning: '25' may need to change
                if(iwrite.eq.0) then
                  write(ichan(j),rec=irec(j)+1) floatout
                  irec(j)=irec(j)+1
                else !evm
                  a_4 = transfer(source=floatout,mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
                endif
              enddo !k
            endif !j
          enddo !i=1,np

          if(nscreen.eq.1) write(*,'(a48)')'done outputting '//variable_nm(j)
          write(16,'(a48)')'done outputting '//variable_nm(j)
        endif !iof(j).eq.1.and.mod(it,nspool).eq.0
      enddo !j=1,noutput

!...  Test output 
      if(noutgm.eq.1.and.mod(it,nspool).eq.0) then
        write(100,rec=igmp+1) real(time)
        write(100,rec=igmp+2) it
        igmp=igmp+2
        do i=1,ns
          do k=1,nvrt
            write(100,rec=igmp+1) real(su2(k,i))
            write(100,rec=igmp+2) real(sv2(k,i))
            igmp=igmp+2
          enddo !k
        enddo !i
        if(nscreen.eq.1) write(*,*)'done writing test results...'
        write(16,*)'done writing test results...'
      endif !noutgm.eq.1.and.mod(it,nspool).eq.0

!...  Open new file
      if(it==ifile*ihfskip) then
        ifile=ifile+1
        write(ifile_char,'(i12)') ifile
        do i=1,noutput
          close(ichan(i))
        enddo
        close(100)
        call header
      endif

!...  End output section

!...  write out hot start information if nhot=1 and at correct time step
!...  note: the hot start files use a record length of 8 on both 32 bit
!...  workstations and the 64 bit cray.
!...
      if(nhot==1.and.mod(it,ihfskip)==0) then 
!       Convert it to a string 
        write(it_char,'(i12)')it

        open(36,file=it_char//'_hotstart',access='direct',recl=ihot_len)
        write(36,rec=1)it,time,(idry_e(i),(we(j,i),tsel(j,i,1),tsel(j,i,2), &
     &(trel0(j,i,l),trel(j,i,l),l=1,ntracers),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
     &dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),j=1,nvrt),i=1,np),ifile,ifile_char
        close(36)

        if(nscreen.eq.1) write(*,*) 'Hot start written',it,time
        write(16,*) 'Hot start written',it,time
      endif !end writing hot start file

      if(nscreen.eq.1) write(*,*)'Time step=', it,' Time=',time
      write(16,*)'Time step=', it,' Time=',time


!...  end time stepping loop
!...
      enddo !it=iths,ntime

!...  Output max. elevations
      open(32,file='maxelev.gr3')
      write(32,*)
      write(32,*)ne,np
      do i=1,np
        write(32,'(i10,2(1x,e20.14),1x,e9.3)')i,x(i),y(i),elevmax(i)
      enddo !i
      do i=1,ne
        write(32,*)i,3,(nm(i,j),j=1,3)
      enddo !i
      close(32)

      call date_and_time(date,timestamp)
      if(nscreen.eq.1) write(*,*)'Run completes successfully at ',date,timestamp
      write(16,*)'Run completes successfully at ',date,timestamp

!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *
!...  end of program
!...
      stop
      end

!
!********************************************************************************
!	This program solves a tridiagonal system. It was adapted		*
!	from "Numerical recipes in FORTRAN (pp.43 ). 				*
!										*
!	a,b,c,r: input vectors and are not modified by this program.		*
!	       b is the main diagonal, a below and c above. a(1) and c(n) are	*
!	       not used. r is the r.h.s.                        		*
!	nmax: dimension of a etc in the driving routine.			*
!	n: actual rank of the system.						*
!	nc: input; actual # of columns of rhs (<=100). 				*
!	u: output with nc columns (depending on input nc).			*
!	gam: a working array.							*
!********************************************************************************
!

      subroutine tridag(nmax,n,nc,a,b,c,r,u,gam)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      integer, intent(in) :: nmax,n,nc
      real(kind=dbl_kind1), dimension(nmax), intent(in) :: a,b,c
      real(kind=dbl_kind1), dimension(nmax,100), intent(in) :: r
      real(kind=dbl_kind1), dimension(nmax), intent(out) :: gam
      real(kind=dbl_kind1), dimension(nmax,100), intent(out) :: u
!      dimension a(nmax),b(nmax),c(nmax),r(nmax,5),u(nmax,5),gam(nmax)

      if(n.lt.1) then
        write(11,*)'Tridag assumes n>=1'
        stop
      endif
      if(nc.gt.100) then
        write(11,*)'Increase # of columns in tridag'
        stop
      endif
      if(b(1).eq.0.) then
        write(11,*)'tridag:  b(1)=0'
        stop
      endif

      bet=b(1)
      do i=1,nc
        u(1,i)=r(1,i)/bet
      enddo

      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0) then
          write(11,*)'tridag failed'
          stop
        endif
        do i=1,nc
          u(j,i)=(r(j,i)-a(j)*u(j-1,i))/bet
        enddo !i
      enddo !j

!     Backsubstitution
      do j=n-1,1,-1
        do i=1,nc
          u(j,i)=u(j,i)-gam(j+1)*u(j+1,i)
        enddo
      enddo
      
      return
      end

!
!********************************************************************
!								    
!	Routine to update z-coordinates and wetting and drying      
!       Use levels1() for better inundation if resolution is fine enough.
!								    
!********************************************************************
!
      subroutine levels0(iths,it)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iths,it

      dimension idry_new(mnp),idry_s_new(mns),out2(12+mnv)

!...  z-coor. for nodes
!...  
      do i=1,np
        if(dp(i)+eta2(i)<=h0) then !dry
          idry_new(i)=1 
          if(dp(i)>=h_s) then
            write(11,*)'Deep depth dry:',i
            stop
          endif
          kbp(i)=0
        else !wet
          idry_new(i)=0
!         S-levels
          do k=kz,nvrt
            kin=k-kz+1

            if(hmod(i)<=h_c) then
              if(ifort12(12)==0) then
                ifort12(12)=1
                write(12,*)'Initial depth too shallow for S:',i,hmod(i),h_c
              endif
              iback(i)=1
              z(k,i)=sigma(kin)*(hmod(i)+eta2(i))+eta2(i)
            else if(eta2(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
              write(11,*)'Pls choose a larger h_c (1):',eta2(i),h_c
              stop
            else
              z(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
            endif
          enddo !k=kz,nvrt

!         z-levels
          if(dp(i)<=h_s) then
            kbp(i)=kz
          else !bottom index 
            if(imm==1) then
              kbp(i)=0 !flag
              do k=1,kz-1
                if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                  kbp(i)=k
                  exit
                endif
              enddo !k
              if(kbp(i)==0) then
                write(11,*)'Cannot find a bottom level for node (3):',i
                stop
              endif
            endif !imm

            if(kbp(i)>=kz.or.kbp(i)<1) then
              write(11,*)'Impossible 92:',kbp(i),kz,i
              stop
            endif
            z(kbp(i),i)=-dp(i)
            do k=kbp(i)+1,kz-1
              z(k,i)=ztot(k)
            enddo !k
          endif

          do k=kbp(i)+1,nvrt
            if(z(k,i)-z(k-1,i)<=0) then
              write(11,*)'Inverted z-levels at:',i,k,z(k,i)-z(k-1,i),eta2(i),hmod(i)
              stop
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,np

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake wet nodes; a node is wet if and only if at least one surrounding element is wet
!...
      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()
      idry_s_new(1:np)=idry(1:np) !temporary save
      idry=1 !dry unless wet
      kbe=0
      do i=1,ne
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        idry_e(i)=max0(idry_new(n1),idry_new(n2),idry_new(n3))
        if(idry_e(i)==0) then
          idry(n1)=0; idry(n2)=0; idry(n3)=0
          kbe(i)=max0(kbp(n1),kbp(n2),kbp(n3))
          do k=kbe(i),nvrt
            ze(k,i)=(z(k,n1)+z(k,n2)+z(k,n3))/3
            if(k>=kbe(i)+1.and.ze(k,i)-ze(k-1,i)<=0) then
              write(11,*)'Weird element:'
              write(11,*)k,i,ze(k,i),ze(k-1,i)
              stop
            endif
          enddo !k
        endif
      enddo !i

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
!     Vel. as average; back-up S,T
      do i=1,np
        if(it/=iths.and.idry_s_new(i)==1.and.idry(i)==0) then
          do k=1,nvrt
            uu2(k,i)=0
            vv2(k,i)=0
            tnd(k,i)=0
            snd(k,i)=0
            icount=0
            do j=1,nnp(i)
              nd=inp(i,j)
              if(idry_s_new(nd)==0) then !all indices extended
                icount=icount+1
                uu2(k,i)=uu2(k,i)+uu2(k,nd)
                vv2(k,i)=vv2(k,i)+vv2(k,nd)
                tnd(k,i)=tnd(k,i)+tnd(k,nd)
                snd(k,i)=snd(k,i)+snd(k,nd)
              endif
            enddo !j
            if(icount==0) then
              if(ifort12(7)==0) then
                ifort12(7)=1
                write(12,*)'Isolated rewetted node:',i
              endif
              tnd(k,i)=tem0(k,i)
              snd(k,i)=sal0(k,i)
            else
              uu2(k,i)=uu2(k,i)/icount
              vv2(k,i)=vv2(k,i)/icount
              tnd(k,i)=tnd(k,i)/icount
              snd(k,i)=snd(k,i)/icount
            endif
          enddo !k=1,nvrt
        endif !rewetted
      enddo !i=1,np

      do i=1,np
        if(it/=iths.and.idry_s_new(i)==1.and.idry(i)==0) then
          do k=kbp(i),nvrt
            x0=x(i)
            y0=y(i)
            uuint=uu2(k,i)
            vvint=vv2(k,i)
            vmag=dsqrt(uuint**2+vvint**2)
            nnel=ine(i,1) !any element
            jlev=k
            if(idry_e0(nnel)==0) then
              write(11,*)'Starting element must be dry'
              stop
            endif
            if(vmag/=0) then
              uuint=uuint/vmag
              vvint=vvint/vmag
              call upwindtrack(i,jlev,nnel,x0,y0,uuint,vvint,out2,nfl)
              if(nfl==0) then
                tnd(k,i)=out2(1)
                snd(k,i)=out2(2)
              endif
            endif !vmag/=0
          enddo !k=kbp(i),nvrt
          do k=1,kbp(i)-1
            tnd(k,i)=tnd(kbp(i),i)
            snd(k,i)=snd(kbp(i),i)
          enddo !k
        endif !rewetted
      enddo !i=1,np

!...  z-coor. for sides
!...  A side is wet if and only if at least one of its elements is wet
      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        idry_s_new(i)=1
        do j=1,2 !elements
          ie=is(i,j)
          if(ie/=0.and.idry_e(ie)==0) idry_s_new(i)=0
        enddo !j

        kbs(i)=0 !dry
        if(idry_s_new(i)==0) then !wet side with 2 wet nodes
          if(dps(i)+(eta2(n1)+eta2(n2))/2<=h0) then
            write(11,*)'Weird side:',i,n1,n2,eta2(n1),eta2(n2)
            stop
          endif
          kbs(i)=max0(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(z(k,n1)+z(k,n2))/2
            if(k>=kbs(i)+1.and.zs(k,i)-zs(k-1,i)<=0) then
              write(11,*)'Weird side:'
              write(11,*)k,n1,n2,z(k,n1),z(k,n2),z(k-1,n1),z(k-1,n2)
              stop
            endif
          enddo !k
        endif
      enddo !i=1,ns

!     Compute vel., S,T for re-wetted sides 
!     Vel. as average; back-up S,T
      do i=1,ns
        if(it/=iths.and.idry_s(i)==1.and.idry_s_new(i)==0) then
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          do k=1,nvrt
            su2(k,i)=0
            sv2(k,i)=0
            tsd(k,i)=0
            ssd(k,i)=0
            icount=0
            do j=1,2
              ie=is(i,j)
              if(ie/=0) then
                do jj=1,3 !side
                  isd=js(ie,jj)
                  if(idry_s(isd)==0) then
                    icount=icount+1
                    su2(k,i)=su2(k,i)+su2(k,isd)
                    sv2(k,i)=sv2(k,i)+sv2(k,isd)
                    tsd(k,i)=tsd(k,i)+tsd(k,isd)
                    ssd(k,i)=ssd(k,i)+ssd(k,isd)
                  endif
                enddo !jj
              endif !ie/=0
            enddo !j
            if(icount==0) then
              if(ifort12(10)==0) then
                ifort12(10)=1
                write(12,*)'Isolated rewetted side:',i,n1,n2
              endif
              tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
              ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
            else
              su2(k,i)=su2(k,i)/icount
              sv2(k,i)=sv2(k,i)/icount
              tsd(k,i)=tsd(k,i)/icount
              ssd(k,i)=ssd(k,i)/icount
            endif
          enddo !k
        endif !rewetted
      enddo !i=1,ns

      do i=1,ns
        if(it/=iths.and.idry_s(i)==1.and.idry_s_new(i)==0) then
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          do k=kbs(i),nvrt
            x0=xcj(i)
            y0=ycj(i)
            uuint=su2(k,i)
            vvint=sv2(k,i)
            vmag=dsqrt(uuint**2+vvint**2)
            nnel=is(i,1) !any element
            jlev=k
            if(idry_e0(nnel)==0) then
              write(11,*)'Starting element must be dry (2)'
              stop
            endif
            if(vmag/=0) then
              uuint=uuint/vmag
              vvint=vvint/vmag
              call upwindtrack(i,jlev,nnel,x0,y0,uuint,vvint,out2,nfl)
              if(nfl==0) then
                tsd(k,i)=out2(1)
                ssd(k,i)=out2(2)
              endif
            endif
          enddo !k=kbs(i),nvrt
          do k=1,kbs(i)-1
            tsd(k,i)=tsd(kbs(i),i)
            ssd(k,i)=ssd(kbs(i),i)
          enddo !k
        endif !rewetted
      enddo !i=1,ns

      idry_s=idry_s_new

      return
      end

!
!********************************************************************
!								    
!	Inundation routine to update z-coordinates and wetting and drying.
!       Better for fine resolution.
!								    
!********************************************************************
!
      subroutine levels1(iths,it)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iths,it

      dimension idry2(mnp),idry_s2(mns),idry_e2(mne),out2(12+mnv),sutmp(mnv),svtmp(mnv)
      dimension isdf(mns),inew(mns),icolor(mnp),icolor2(mns)

!...  An element is wet if and only if depths at all nodes >h0 
!...  A node is wet if and only if at least one surrounding element is wet
!...  A side is wet if and only if at least one surrounding element is wet
!     Initialize element flags for first step
      if(it==iths) then
        idry_e=0
        do i=1,ne
          do j=1,3
            nd=nm(i,j)
            if(eta2(nd)+dp(nd)<=h0) then
              idry_e(i)=1
              exit
            endif
          enddo !j
        enddo !i
      endif !it

      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

!...  Wetting/drying algorithm
      idry_e2=idry_e !starting from step n's indices
      if(it/=iths) then
        istop=0 !stop iteration and go to extrapolation stage
        itr=0
        loop15: do
          itr=itr+1
          if(itr>200) then
            write(11,*)'Too many iterations in wet/dry'
            stop
          endif

!         Interface sides
          nsdf=0
          icolor=0 !nodes on the interface sides
          icolor2=0 !interface sides
          do i=1,ns
            if(is(i,2)/=0.and.idry_e2(is(i,1))+idry_e2(is(i,2))==1) then
              nsdf=nsdf+1
              isdf(nsdf)=i
              icolor(isidenode(i,1:2))=1
              icolor2(i)=1
            endif
          enddo !i
          if(nsdf==0) exit loop15

!         Final extrapolation
          if(istop==1) then
            icolor=0 !frontier nodes
            inew=0 !for initializing and counting su2 sv2
            do i=1,nsdf
              isd=isdf(i)
              if(idry_e2(is(isd,1))+idry_e2(is(isd,2))/=1) cycle

              if(idry_e2(is(isd,1))==1) then
                ie=is(isd,1)
              else 
                ie=is(isd,2)
              endif
              n1=isidenode(isd,1)
              n2=isidenode(isd,2)
              nodeA=nm(ie,1)+nm(ie,2)+nm(ie,3)-n1-n2

              if(icolor(nodeA)==1) cycle !this node is done

              icolor(nodeA)=1 !this node is done
              inun=0 !inundation flag
              do j=1,nne(nodeA)
                ie2=ine(nodeA,j)
                id=iself(nodeA,j)
                isd2=js(ie2,id)
                if(icolor2(isd2)==1) then
                  tmp=su2(nvrt,isd2)*snx(isd2)+sv2(nvrt,isd2)*sny(isd2)
                  flux_t=-tmp*ssign(ie2,id) !inward normal
                  if(flux_t>0) then
                    n1=isidenode(isd2,1)
                    n2=isidenode(isd2,2)
!                    avh=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/2
!                    vol=flux_t*dt*avh*distj(isd2) !inflow volume in one step
!                    avh3=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/3 !assume total depth at nodeA=0
!                    volmin=avh3*area(ie2)
                    etm=max(eta2(n1),eta2(n2))
                    if(etm+dp(nodeA)>h0) then
                      inun=1
                      exit
                    endif
                  endif
                endif
              enddo !j

              if(inun==1) then
                eta2(nodeA)=max(eta2(nodeA),-dp(nodeA)+2*h0)
                do j=1,nne(nodeA)
                  ie2=ine(nodeA,j)
                  id=iself(nodeA,j)
                  isd2=js(ie2,id)
                  if(icolor2(isd2)==1) then
                    do l=1,3
                      nd=nm(ie2,l)
                      if(eta2(nd)+dp(nd)<=h0) then 
                        write(11,*)'Failed to wet element:',l
                        stop
                      endif
                    enddo !l=1,3
                    idry_e2(ie2)=0
                    do l=1,2 !sides sharing nodeA
                      id1=js(ie2,nx(id,l))
                      if(inew(id1)==0) then
                        su2(1:nvrt,id1)=su2(1:nvrt,isd2)
                        sv2(1:nvrt,id1)=sv2(1:nvrt,isd2)
                        inew(id1)=1
                      else
                        su2(1:nvrt,id1)=su2(1:nvrt,id1)+su2(1:nvrt,isd2)
                        sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+sv2(1:nvrt,isd2)
                        inew(id1)=inew(id1)+1
                      endif
                    enddo !l=1,2
                  endif !icolor2(isd2)==1
                enddo !j=1,nne(nodeA)
              endif !inun==1
            enddo !i=1,nsdf

            do i=1,ns
              if(inew(i)/=0) then
                su2(1:nvrt,i)=su2(1:nvrt,i)/inew(i)
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/inew(i)
              endif
            enddo !i

!            exit loop15
            istop=2
            go to 991
          endif !istop=1; final extrapolation

          inew=0 !for initializing and counting su2 sv2
          istop=1 !stop iteration and go to extrapolation stage
          do i=1,nsdf
            isd=isdf(i)
            ifl=0 !flag for making dry
            do j=1,2
              nd=isidenode(isd,j)
              if(eta2(nd)+dp(nd)<=h0) then
                istop=0
                ifl=1
                do l=1,nne(nd)
                  idry_e2(ine(nd,l))=1
                enddo !l
              endif
            enddo !j=1,2 nodes
            if(ifl==1) cycle
            
!           2 end nodes have total depths > h0
            if(idry_e2(is(isd,1))+idry_e2(is(isd,2))/=1) cycle

            if(idry_e2(is(isd,1))==1) then
              ie=is(isd,1)
            else
              ie=is(isd,2)
            endif
            n1=isidenode(isd,1)
            n2=isidenode(isd,2)
            nodeA=nm(ie,1)+nm(ie,2)+nm(ie,3)-n1-n2
            l0=lindex(nodeA,ie)
!            if(l0==0.or.icolor(nodeA)==1.or.nodeA==n1.or.nodeA==n2) then
            if(l0==0.or.nodeA==n1.or.nodeA==n2) then
              write(11,*)'Frontier node outside, or on the interface:',l0,icolor(nodeA),nodeA,n1,n2,itr,it,iths
              do l=1,ns
                if(icolor2(l)==1) then
                  write(99,*)l,isidenode(l,1:2)
                  write(99,*)l,is(l,1:2),idry_e2(is(l,1:2)),idry_e(is(l,1:2))
                endif
              enddo !l
              do l=1,ne
                write(98,*)l,idry_e2(l),idry_e(l)
              enddo !l
              stop
            endif

            if(eta2(nodeA)+dp(nodeA)>h0) then !all 3 nodes have depths > h0
!             Check
              do j=1,3
                nd=nm(ie,j)
                if(eta2(nd)+dp(nd)<=h0) then
                  write(11,*)'Failed to wet element (11):',nd,nodeA
                  stop
                endif
              enddo !j

              istop=0
              idry_e2(ie)=0
              do j=1,2 !sides sharing nodeA
                id1=js(ie,nx(l0,j))
                if(icolor2(id1)==0) then
                  if(inew(id1)==0) then
                    su2(1:nvrt,id1)=su2(1:nvrt,isd)
                    sv2(1:nvrt,id1)=sv2(1:nvrt,isd)
                    inew(id1)=1
                  else
                    su2(1:nvrt,id1)=su2(1:nvrt,id1)+su2(1:nvrt,isd)
                    sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+sv2(1:nvrt,isd)
                    inew(id1)=inew(id1)+1
                  endif
                endif !icolor2(id)==0
              enddo !j
            endif
          enddo !i=1,nsdf; interfacial sides

!         Compute average vel. for rewetted sides
          do i=1,ns
            if(inew(i)/=0) then
              su2(1:nvrt,i)=su2(1:nvrt,i)/inew(i)
              sv2(1:nvrt,i)=sv2(1:nvrt,i)/inew(i)
            endif
          enddo !i

!         Enforce wet/dry flag consistency between nodes and elements due to added wet elements
991       continue
          idry2=1
          do i=1,ne
            if(idry_e2(i)==0) idry2(nm(i,1:3))=0
          enddo !i
          do i=1,ne
            itmp=0
            do j=1,3
              if(idry2(nm(i,j))==1) itmp=1
            enddo !j
!           Compute su2 sv2 for rewetted sides
            if(idry_e2(i)==1.and.itmp==0) then
              sutmp=0
              svtmp=0
              icount=0
              do j=1,3
                isd=js(i,j)
                if(idry_e2(is(isd,1))==0.or.is(isd,2)/=0.and.idry_e2(is(isd,2))==0) then !at least one wet element
                  icount=icount+1
                  sutmp(1:nvrt)=sutmp(1:nvrt)+su2(1:nvrt,isd)
                  svtmp(1:nvrt)=svtmp(1:nvrt)+sv2(1:nvrt,isd)
                endif
              enddo !j
              if(icount/=0) then
                do j=1,3
                  isd=js(i,j)
                  if(idry_e2(is(isd,1))==1) then
                    if(is(isd,2)==0.or.is(isd,2)/=0.and.idry_e2(is(isd,2))==1) then
                      su2(1:nvrt,isd)=sutmp(1:nvrt)/icount
                      sv2(1:nvrt,isd)=svtmp(1:nvrt)/icount
                    endif
                  endif
                enddo !j
              endif
            endif !rewetted sides

            idry_e2(i)=itmp
          enddo !i=1,ne
          if(istop==2) exit loop15

        end do loop15
      endif !it/=iths

!...  Isolated dry nodes (do nothing for isolated wet)
      do i=1,np
        if(dp(i)+eta2(i)<=h0) idry_e2(ine(i,1:nne(i)))=1
      enddo !i

!...  Wet/dry flags for nodes/sides
      idry2=1; idry_s2=1
      do i=1,ne
        if(idry_e2(i)==0) then
          idry2(nm(i,1:3))=0
          idry_s2(js(i,1:3))=0
        endif
      enddo !i

!...  Reset vel. at dry sides
      do i=1,ns
        if(idry_s2(i)==1) then
          su2(1:nvrt,i)=0
          sv2(1:nvrt,i)=0
        endif
      enddo !i

!...  Reset elevation at dry nodes
      do i=1,np
        if(idry2(i)==1) then
          eta2(i)=min(0.d0,-dp(i))
        endif
      enddo !i

!...  z-coor. for nodes
!...  
      do i=1,np
        if(eta2(i)<=h0-h_s) then
          write(11,*)'Deep depth dry:',i
          stop
        endif

        if(idry2(i)==1) then
          kbp(i)=0
        else !wet
          if(dp(i)+eta2(i)<=h0) then
            write(11,*)'Mismatch in node index (2):',i,dp(i)+eta2(i)
            stop
          endif
!         S-levels
          do k=kz,nvrt
            kin=k-kz+1

            if(hmod(i)<=h_c) then
              if(ifort12(12)==0) then
                ifort12(12)=1
                write(12,*)'Initial depth too shallow for S:',i,hmod(i),h_c
              endif
              iback(i)=1
              z(k,i)=sigma(kin)*(hmod(i)+eta2(i))+eta2(i)
            else if(eta2(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
              write(11,*)'Pls choose a larger h_c (1):',eta2(i),h_c
              stop
            else
              z(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
            endif
          enddo !k=kz,nvrt

!         z-levels
          if(dp(i)<=h_s) then
            kbp(i)=kz
          else !bottom index 
!            if(imm==1) then
!           Redo every step for wet/dry changes
            kbp(i)=0 !flag
            do k=1,kz-1
              if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                kbp(i)=k
                exit
              endif
            enddo !k
            if(kbp(i)==0) then
              write(11,*)'Cannot find a bottom level for node (3):',i
              stop
            endif
!            endif !imm

            if(kbp(i)>=kz.or.kbp(i)<1) then
              write(11,*)'Impossible 92:',kbp(i),kz,i
              stop
            endif
            z(kbp(i),i)=-dp(i)
            do k=kbp(i)+1,kz-1
              z(k,i)=ztot(k)
            enddo !k
          endif

          do k=kbp(i)+1,nvrt
            if(z(k,i)-z(k-1,i)<=0) then
              write(11,*)'Inverted z-levels at:',i,k,z(k,i)-z(k-1,i),eta2(i),hmod(i)
              stop
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,np

!...  Z-coord. for elements
      kbe=0
      do i=1,ne
        if(idry_e2(i)==0) then
          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)
          kbe(i)=max0(kbp(n1),kbp(n2),kbp(n3))
          do k=kbe(i),nvrt
            ze(k,i)=(z(k,n1)+z(k,n2)+z(k,n3))/3
            if(k>=kbe(i)+1.and.ze(k,i)-ze(k-1,i)<=0) then
              write(11,*)'Weird element:'
              write(11,*)k,i,ze(k,i),ze(k-1,i)
              stop
            endif
          enddo !k
        endif
      enddo !i

!...  z-coor. for sides
!...  A side is wet if and only if at least one of its elements is wet
      do i=1,ns
        kbs(i)=0 !dry
        if(idry_s2(i)==0) then !wet side with 2 wet nodes
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          if(dps(i)+(eta2(n1)+eta2(n2))/2<=h0) then
            write(11,*)'Weird side:',i,n1,n2,eta2(n1),eta2(n2)
            stop
          endif
          kbs(i)=max0(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(z(k,n1)+z(k,n2))/2
            if(k>=kbs(i)+1.and.zs(k,i)-zs(k-1,i)<=0) then
              write(11,*)'Weird side:'
              write(11,*)k,n1,n2,z(k,n1),z(k,n2),z(k-1,n1),z(k-1,n2)
              stop
            endif
          enddo !k
        endif
      enddo !i=1,ns

!     Compute S,T for re-wetted nodes (q2 and xl are fine)
!     back-up S,T
!     uu2 and vv2 are first estimates; will be overwritten by nodalvel later
      do i=1,np
        if(it/=iths.and.idry(i)==1.and.idry2(i)==0) then
          do k=1,nvrt
            uu2(k,i)=0
            vv2(k,i)=0
            tnd(k,i)=0
            snd(k,i)=0
            icount=0
            do j=1,nnp(i)
              nd=inp(i,j)
              if(idry(nd)==0) then !all indices extended
                icount=icount+1
                uu2(k,i)=uu2(k,i)+uu2(k,nd)
                vv2(k,i)=vv2(k,i)+vv2(k,nd)
                tnd(k,i)=tnd(k,i)+tnd(k,nd)
                snd(k,i)=snd(k,i)+snd(k,nd)
              endif
            enddo !j
            if(icount==0) then
              if(ifort12(7)==0) then
                ifort12(7)=1
                write(12,*)'Isolated rewetted node:',i
              endif
              tnd(k,i)=tem0(k,i)
              snd(k,i)=sal0(k,i)
            else
              uu2(k,i)=uu2(k,i)/icount
              vv2(k,i)=vv2(k,i)/icount
              tnd(k,i)=tnd(k,i)/icount
              snd(k,i)=snd(k,i)/icount
            endif
          enddo !k=1,nvrt
        endif !rewetted
      enddo !i=1,np

      do i=1,np
        if(it/=iths.and.idry(i)==1.and.idry2(i)==0) then
          do k=kbp(i),nvrt
            x0=x(i)
            y0=y(i)
            uuint=uu2(k,i)
            vvint=vv2(k,i)
            vmag=dsqrt(uuint**2+vvint**2)
            nnel=ine(i,1) !any element
            jlev=k
            if(idry_e0(nnel)==0) then
              write(11,*)'Starting element must be dry'
              stop
            endif
            if(vmag/=0) then
              uuint=uuint/vmag
              vvint=vvint/vmag
              call upwindtrack(i,jlev,nnel,x0,y0,uuint,vvint,out2,nfl)
              if(nfl==0) then
                tnd(k,i)=out2(1)
                snd(k,i)=out2(2)
              endif
            endif !vmag/=0
          enddo !k=kbp(i),nvrt
          do k=1,kbp(i)-1
            tnd(k,i)=tnd(kbp(i),i)
            snd(k,i)=snd(kbp(i),i)
          enddo !k
        endif !rewetted
      enddo !i=1,np

!     Compute vel., S,T for re-wetted sides 
!     back-up S,T
      do i=1,ns
        if(it/=iths.and.idry_s(i)==1.and.idry_s2(i)==0) then
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          do k=1,nvrt
!            su2(k,i)=0
!            sv2(k,i)=0
            tsd(k,i)=0
            ssd(k,i)=0
            icount=0
            do j=1,2
              ie=is(i,j)
              if(ie/=0) then
                do jj=1,3 !side
                  isd=js(ie,jj)
                  if(isd/=i.and.idry_s(isd)==0) then !wet at step n
                    icount=icount+1
!                    su2(k,i)=su2(k,i)+su2(k,isd)
!                    sv2(k,i)=sv2(k,i)+sv2(k,isd)
                    tsd(k,i)=tsd(k,i)+tsd(k,isd)
                    ssd(k,i)=ssd(k,i)+ssd(k,isd)
                  endif
                enddo !jj
              endif !ie/=0
            enddo !j
            if(icount==0) then
              if(ifort12(10)==0) then
                ifort12(10)=1
                write(12,*)'Isolated rewetted side:',i,n1,n2
              endif
              tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
              ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
            else
!              su2(k,i)=su2(k,i)/icount
!              sv2(k,i)=sv2(k,i)/icount
              tsd(k,i)=tsd(k,i)/icount
              ssd(k,i)=ssd(k,i)/icount
            endif
          enddo !k
        endif !rewetted
      enddo !i=1,ns

      do i=1,ns
        if(it/=iths.and.idry_s(i)==1.and.idry_s2(i)==0) then
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          do k=kbs(i),nvrt
            x0=xcj(i)
            y0=ycj(i)
            uuint=su2(k,i)
            vvint=sv2(k,i)
            vmag=dsqrt(uuint**2+vvint**2)
            nnel=is(i,1) !any element
            jlev=k
            if(idry_e0(nnel)==0) then
              write(11,*)'Starting element must be dry (2)'
              stop
            endif
            if(vmag/=0) then
              uuint=uuint/vmag
              vvint=vvint/vmag
              call upwindtrack(i,jlev,nnel,x0,y0,uuint,vvint,out2,nfl)
              if(nfl==0) then
                tsd(k,i)=out2(1)
                ssd(k,i)=out2(2)
              endif
            endif
          enddo !k=kbs(i),nvrt
          do k=1,kbs(i)-1
            tsd(k,i)=tsd(kbs(i),i)
            ssd(k,i)=ssd(kbs(i),i)
          enddo !k
        endif !rewetted
      enddo !i=1,ns

!...  Update wet/dry flags
      idry=idry2
      idry_s=idry_s2
      idry_e=idry_e2

      return
      end

!
!***************************************************************************
!									   
!     			Routine for backtracking. 			   
!       Input:
!             ielem: initial starting element;
!             l_ns: starting from nodes (l_ns<=3) or sides (3<l_ns<=6);
!             j0:   starting level;
!             iadvf: advection flag (0 or 1: use Euler tracking; 2: R-K tracking);
!             dtb_max: min. (for Euler) or max. (for R-K) tracking step;
!             x00,y00,z00: starting coordinates;
!             
!       Input/output:
!             uuint,vvint,wwint: starting and interpolated vel.;
!             nnel,jlev: initial and final element and level;
!
!       Output:
!             xt,yt,zt: final location;
!             out6(6): only first 2 are used for interpolated T,S;
!             idt: total # of steppings for R-K;
!									   
!***************************************************************************
!
      subroutine btrack(ielem,l_ns,j0,iadvf,dtb_max,uuint,vvint,wwint,x00,y00,z00,nnel,jlev,xt,yt,zt,idt,out6)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      real, parameter :: per1=1.e-3 !used to check error in tracking
      real, parameter :: safety=0.8 !safyty factor in estimating time step

      integer, intent(in) :: ielem,l_ns,j0,iadvf
      integer, intent(inout) :: nnel,jlev
      real(kind=dbl_kind), intent(in) :: dtb_max,x00,y00,z00
      real(kind=dbl_kind), intent(inout) :: uuint,vvint,wwint
      real(kind=dbl_kind), intent(out) :: xt,yt,zt,out6(6)
      integer, intent(out) :: idt

      real(kind=dbl_kind) :: vxl(3,2),vyl(3,2),vzl(3,2),vxn(3),vyn(3),vzn(3)
      real(kind=dbl_kind) :: arco(3),t_xi(6),s_xi(6),sig(3),subrat(4),ztmp(mnv), &
     &swild(10),swild2(mnv,10),swild3(mnv) !swild2's dimension must match vinter()
      real(kind=dbl_kind) :: tsdata(mnei_kr,2),al_beta(mnei_kr+3,4),uvdata(mnei_kr,2)
      real(kind=dbl_kind) :: a(7),b(7,6),dc(6),rk(6,3)

!     Constants used in 5th order R-K
      a(2)=0.2; a(3)=0.3; a(4)=0.6; a(5)=1; a(6)=0.875; a(7)=1
      b(2,1)=0.2; b(3,1)=0.075; b(3,2)=0.225; b(4,1)=0.3; b(4,2)=-0.9; b(4,3)=1.2
      b(5,1)=-11./54; b(5,2)=2.5; b(5,3)=-70./27; b(5,4)=35./27
      b(6,1)=1631./55296; b(6,2)=175./512; b(6,3)=575./13824; b(6,4)=44275./110592; b(6,5)=253./4096
!     b(7,*) are c(*)
      b(7,1)=37./378; b(7,2)=0; b(7,3)=250./621; b(7,4)=125./594; b(7,5)=0; b(7,6)=512./1771
      dc(1)=b(7,1)-2825./27648; dc(2)=0; dc(3)=b(7,3)-18575./48384
      dc(4)=b(7,4)-13525./55296; dc(5)=b(7,5)-277./14336; dc(6)=b(7,6)-0.25

      wwint00=wwint !save for quadratic interpolation

!...  Choose vis_coe
      if(indvel==1) then
        vis_coe=0
      else !indvel=0
        if(l_ns<=3) then
          vis_coe=vis_coe2
        else !sides
          isd=js(ielem,l_ns-3)
          if(is(isd,2)==0) then
            vis_coe=vis_coe2
          else
            vis_coe=vis_coe1
          endif
        endif      
      endif

!...  Euler tracking
      if(iadvf==1.or.iadvf==0) then
!------------------------------------------------------------------------------------------------
!     Compute # of sub-division based on local gradients
      ndelt_max=dt/dtb_max !dtb_max is min. dtb actualy in this case
      if(l_ns<=3) then !nodes
        nd=nm(ielem,l_ns)
        sum=0
        do i=1,nne(nd)
          ie=ine(nd,i)
          id=iself(nd,i)
          dudx=0; dudy=0; dvdx=0; dvdy=0
          do j=1,3
            dudx=dudx+ufg(j0,ie,j)*dl(ie,j,1) !not strictly along z
            dudy=dudy+ufg(j0,ie,j)*dl(ie,j,2) !not strictly along z
            dvdx=dvdx+vfg(j0,ie,j)*dl(ie,j,1) 
            dvdy=dvdy+vfg(j0,ie,j)*dl(ie,j,2) 
          enddo !j
          sum=sum+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)/nne(nd)
        enddo !i=1,nne(nd)
        ndelt=max0(1,min0(ndelt_max,int(sum)*4)) !>=1
      else !sides
        isd=js(ielem,l_ns-3)
        sum=0
        icount=0
        do i=1,2
          ie=is(isd,i)
          if(ie==0) cycle
          icount=icount+1

          dudx=0; dudy=0; dvdx=0; dvdy=0
          do j=1,3
            dudx=dudx+ufg(j0,ie,j)*dl(ie,j,1) !not strictly along z
            dudy=dudy+ufg(j0,ie,j)*dl(ie,j,2) !not strictly along z
            dvdx=dvdx+vfg(j0,ie,j)*dl(ie,j,1)
            dvdy=dvdy+vfg(j0,ie,j)*dl(ie,j,2)
          enddo !j
          sum=sum+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)
        enddo !i=1,2
        if(icount==0) then
          write(11,*)'Impossible 77'
          stop
        endif
        ndelt=max0(1,min0(ndelt_max,int(sum/icount)*4)) !>=1
      endif

      dtb=dt/ndelt
      nnel0=nnel
      jlev0=jlev
      x0=x00
      y0=y00
      z0=z00
      do idt=1,ndelt
        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint
        call quicksearch(1,idt,ielem,nnel0,jlev0,dtb,x0,y0,z0,xt,yt,zt,nnel,jlev,arco,zrat,ztmp,kbpl,iflqs1,ss)

!	nnel wet
!	No interpolate in time
        do j=1,3
          nd=nm(nnel,j)
          do l=1,2
            lev=jlev+l-2
            vxl(j,l)=(1-vis_coe)*uu2(lev,nd)+vis_coe*ufg(lev,nnel,j)
            vyl(j,l)=(1-vis_coe)*vv2(lev,nd)+vis_coe*vfg(lev,nnel,j)
            vzl(j,l)=ww2(lev,nd)
          enddo !l
        enddo !j

!	Interpolate in vertical 
        do j=1,3
          if(interpol(nnel)==1) then
            nd=nm(nnel,j)
            kbb=kbp(nd)
            swild3(kbb:nvrt)=z(kbb:nvrt,nd)
            swild2(kbb:nvrt,1)=(1-vis_coe)*uu2(kbb:nvrt,nd)+vis_coe*ufg(kbb:nvrt,nnel,j)
            swild2(kbb:nvrt,2)=(1-vis_coe)*vv2(kbb:nvrt,nd)+vis_coe*vfg(kbb:nvrt,nnel,j)
            swild2(kbb:nvrt,3)=ww2(kbb:nvrt,nd)

            call vinter(mnv,3,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
            vxn(j)=swild(1)
            vyn(j)=swild(2)
            vzn(j)=swild(3)

          else !pure S region
            vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
            vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
            vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
          endif
        enddo !j

!	Interpolate in horizontal
        uuint=vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
        vvint=vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
        wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

        if(iflqs1==1) exit 
        x0=xt
        y0=yt
        z0=zt
        nnel0=nnel
        jlev0=jlev
      enddo !idt=1,ndelt
!------------------------------------------------------------------------------------------------
      endif !Euler

!...  5th-order R-K tracking
      if(iadvf==2) then
!------------------------------------------------------------------------------------------------
      dtb0=dmin1(dtb_max,dt)
      t_m=dt !t^m
      idt=0
      rk(1,1)=-dtb0*uuint
      rk(1,2)=-dtb0*vvint
      rk(1,3)=-dtb0*wwint
!     Initialize nnel0 and jlev0 
      nnel0=nnel
      jlev0=jlev
      x0=x00
      y0=y00
      z0=z00
      uuint0=uuint
      vvint0=vvint
      wwint0=wwint

      icount=0      
      loop5: do 
        idt=idt+1
        icount=icount+1

!       k2-6 and k1 for the next step
        do k=2,7 !k=7 --> k1
          xt=x0
          yt=y0
          zt=z0
          do l=1,k-1
            xt=xt+rk(l,1)*b(k,l)
            yt=yt+rk(l,2)*b(k,l)
            zt=zt+rk(l,3)*b(k,l)
          enddo !l

          call quicksearch(1,idt,ielem,nnel0,jlev0,dtb0*a(k),x0,y0,z0,xt,yt,zt,nnel,jlev,arco,zrat,ztmp,kbpl,iflqs1,ss)

          if(k==7) then
            dx=0
            dy=0
            dz=0
            do l=1,k-1
              dx=dx+rk(l,1)*dc(l)
              dy=dy+rk(l,2)*dc(l)
              dz=dz+rk(l,3)*dc(l)
            enddo !l
            del_xy=dsqrt(dx*dx+dy*dy)
            del_z=dabs(dz)
            n1=nm(nnel,1) !wet node
            del0_xy=per1*radiel(nnel)
            jmin=max0(jlev,kbp(n1)+1)
            del0_z=per1*(z(jmin,n1)-z(jmin-1,n1))
            if(del0_z<=0) then
              write(11,*)'Negative layer:',del0_z
              stop
            endif
            if(del_xy==0) then
              dtb_xy=dtb0
            else 
              dtb_xy=safety*dtb0*(del0_xy/del_xy)**0.2
            endif
            if(del_z==0) then
              dtb_z=dtb0
            else 
              dtb_z=safety*dtb0*(del0_z/del_z)**0.2
            endif
!           Proposed time step for next iteration
            dtb=dmin1(dtb_xy,dtb_z,t_m,dtb_max) !dtb_xy,dtb_z,t_m >0
            if(del_xy>del0_xy.or.del_z>del0_z) then
!             Debug
!             write(90,*)ielem,icount,t_m,dtb0,dtb,idt

              if(icount<=5) then
!               Update dtb0 and k1, and redo current step
                dtb0=dtb
                rk(1,1)=-dtb0*uuint0
                rk(1,2)=-dtb0*vvint0
                rk(1,3)=-dtb0*wwint0
                cycle loop5
              else !trap reached
                if(ifort12(2)==0) then
                  ifort12(2)=1
                  write(12,*)'Adaptivity trap reached after ',icount
                endif
              endif
            endif

!           Successful
            t_m=t_m-dtb0
            if(t_m<0) then
              write(11,*)'Negative time level:',t_m
              stop
            endif
!           Update dtb0
            dtb0=dmin1(dtb,t_m)
            if(dtb0==0.and.t_m/=0) then
              write(11,*)'Weird btrack:',dtb0,t_m
              stop
            endif
            icount=0 !reset

!           Debug
!            close(90)
!            open(90,file='fort.90')
!            rewind(90)
          endif !k==7

!	  nnel wet
!	  Interpolate in time
          if(k==7) then
            trat=t_m/dt !t_m updated
          else
            trat=(t_m-a(k)*dtb0)/dt !dtb0 not updated for k=2-6
          endif
          do j=1,3
            nd=nm(nnel,j)
            do l=1,2
              lev=jlev+l-2
              vxl(j,l)=(1-vis_coe)*uu2(lev,nd)+vis_coe*ufg(lev,nnel,j)
              vyl(j,l)=(1-vis_coe)*vv2(lev,nd)+vis_coe*vfg(lev,nnel,j)
              vzl(j,l)=ww2(lev,nd)
            enddo !l
          enddo !j

!	  Interpolate in vertical 
          do j=1,3
            if(interpol(nnel)==1) then
              nd=nm(nnel,j)
              kbb=kbp(nd)
              swild3(kbb:nvrt)=z(kbb:nvrt,nd)
              swild2(kbb:nvrt,1)=(1-vis_coe)*uu2(kbb:nvrt,nd)+vis_coe*ufg(kbb:nvrt,nnel,j)
              swild2(kbb:nvrt,2)=(1-vis_coe)*vv2(kbb:nvrt,nd)+vis_coe*vfg(kbb:nvrt,nnel,j)
              swild2(kbb:nvrt,3)=ww2(kbb:nvrt,nd)
              call vinter(mnv,3,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              vxn(j)=swild(1)
              vyn(j)=swild(2)
              vzn(j)=swild(3)
            else !pure S region
              vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
              vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
              vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
            endif
          enddo !j

!	  Interpolate in horizontal
          uuint=vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
          vvint=vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
          wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

          if(k==7) then
            in=1
          else
            in=k
          endif
          rk(in,1)=-dtb0*uuint !dtb0 updated for k=7
          rk(in,2)=-dtb0*vvint
          rk(in,3)=-dtb0*wwint
!          if(iflqs1==1) exit loop5
        enddo !k=2,7
        x0=xt
        y0=yt
        z0=zt
        nnel0=nnel
        jlev0=jlev
        uuint0=uuint
        vvint0=vvint
        wwint0=wwint

        if(t_m<=0) exit loop5
      end do loop5
!------------------------------------------------------------------------------------------------
      endif !R-K

!     Kriging for vel. (excluding bnd nodes/sides)
      if(krvel(nnel)==1.and.(l_ns<=3.and.isbnd(nm(ielem,l_ns))==0.or.l_ns>3.and.is(js(ielem,l_ns-3),2)/=0)) then
!       Prepare data
        ie=ie_kr(nnel) !local index
        if(ie==0) then
          write(11,*)'Out of Kriging zone:',nnel
          stop
        endif
        npp=itier_nd(ie,0)
        do i=1,npp
          nd=itier_nd(ie,i)
          if(idry(nd)==1) then !i.c.
            uvdata(i,1)=0
            uvdata(i,2)=0
          else !wet
            if(interpol(nnel)==1) then
              kbb=kbp(nd)
              swild3(kbb:nvrt)=z(kbb:nvrt,nd)
              swild2(kbb:nvrt,1)=uu2(kbb:nvrt,nd)
              swild2(kbb:nvrt,2)=vv2(kbb:nvrt,nd)
              call vinter(mnv,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              uvdata(i,1:2)=swild(1:2)
            else !along S
              uvdata(i,1)=uu2(jlev,nd)*(1-zrat)+uu2(jlev-1,nd)*zrat
              uvdata(i,2)=vv2(jlev,nd)*(1-zrat)+vv2(jlev-1,nd)*zrat
            endif
          endif
        enddo !all ball nodes

        do i=1,npp+3
          al_beta(i,1:2)=0
          do j=1,npp
            al_beta(i,1:2)=al_beta(i,1:2)+akrmat_nd(ie,i,j)*uvdata(j,1:2)
          enddo !j
        enddo !i

        uuint=al_beta(npp+1,1)+al_beta(npp+2,1)*xt+al_beta(npp+3,1)*yt
        vvint=al_beta(npp+1,2)+al_beta(npp+2,2)*xt+al_beta(npp+3,2)*yt
        do i=1,npp
          nd=itier_nd(ie,i)
          rr=dsqrt((x(nd)-xt)**2+(y(nd)-yt)**2)
          covar2=covar(kr_co,decorrel(ie),rr)
          uuint=uuint+al_beta(i,1)*covar2
          vvint=vvint+al_beta(i,2)*covar2
        enddo !i
      endif !Kriging vel.

!...  Interpolation at the foot for S,T
!...
      out6=0 !initialize
!     nnel wet
      if(zrat<0.or.zrat>1) then
        write(11,*)'zrat wrong:',jlev
        stop
      endif
 
!     Split-linear, quadratic or Kriging
      if(lqk(nnel)==1) then
!-----------------------------------------------------------------------
!     Split-linear
      index=0
      do i=1,4
        if(i<=3) then
          n1=nm(nnel,i)
          n2=js(nnel,nx(i,2))
          n3=js(nnel,nx(i,1))
          aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
          aa2=signa(x(n1),xt,xcj(n3),y(n1),yt,ycj(n3))
          aa3=signa(x(n1),xcj(n2),xt,y(n1),ycj(n2),yt)
          aa=dabs(aa1)+dabs(aa2)+dabs(aa3)
          subrat(i)=dabs(aa-area(nnel)/4)*4/area(nnel)
          if(subrat(i)<100*small1) then
            index=1
            sig(1)=aa1*4/area(nnel)
            sig(2)=aa2*4/area(nnel)
            sig(1)=dmax1(0.0d0,dmin1(1.0d0,sig(1)))
            sig(2)=dmax1(0.0d0,dmin1(1.0d0,sig(2)))
            if(sig(1)+sig(2)>1) then
              sig(3)=0
              sig(2)=1-sig(1)
            else
              sig(3)=1-sig(1)-sig(2)
            endif

!           S,T extended
!            t_xi(1)=tnd(jlev,n1)*(1-zrat)+tnd(jlev-1,n1)*zrat
!            t_xi(2)=tsd(jlev,n2)*(1-zrat)+tsd(jlev-1,n2)*zrat

!            smax=-99; tmin=100; ibb=0 !flag
            do jj=1,3 !node or side
              if(jj==1) then
                kbb=kbp(n1)
                swild3(kbb:nvrt)=z(kbb:nvrt,n1)
                swild2(kbb:nvrt,1)=tnd(kbb:nvrt,n1)
                swild2(kbb:nvrt,2)=snd(kbb:nvrt,n1)
!                swild2(kbb:nvrt,3)=uu2(kbb:nvrt,n1)
!                swild2(kbb:nvrt,4)=vv2(kbb:nvrt,n1)
              else if(jj==2) then
                kbb=kbs(n2)
                swild3(kbb:nvrt)=zs(kbb:nvrt,n2)
                swild2(kbb:nvrt,1)=tsd(kbb:nvrt,n2)
                swild2(kbb:nvrt,2)=ssd(kbb:nvrt,n2)
!                swild2(kbb:nvrt,3)=su2(kbb:nvrt,n2)
!                swild2(kbb:nvrt,4)=sv2(kbb:nvrt,n2)
              else !=3
                kbb=kbs(n3)
                swild3(kbb:nvrt)=zs(kbb:nvrt,n3)
                swild2(kbb:nvrt,1)=tsd(kbb:nvrt,n3)
                swild2(kbb:nvrt,2)=ssd(kbb:nvrt,n3)
!                swild2(kbb:nvrt,3)=su2(kbb:nvrt,n3)
!                swild2(kbb:nvrt,4)=sv2(kbb:nvrt,n3)
              endif

              call vinter(mnv,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
!              if(ibelow==1) ibb=1
!              if(smax<swild(2)) smax=swild(2)
!              if(tmin>swild(1)) tmin=swild(1)
              t_xi(jj)=swild(1); s_xi(jj)=swild(2)

            enddo !jj=1,3

            out6(1)=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
            out6(2)=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)

            exit
          endif !subrat(i)<100*small1
        else !i=4
          n1=js(nnel,1)
          n2=js(nnel,2)
          n3=js(nnel,3)
          aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
          aa2=signa(xcj(n1),xt,xcj(n3),ycj(n1),yt,ycj(n3))
          aa3=signa(xcj(n1),xcj(n2),xt,ycj(n1),ycj(n2),yt)
          aa=dabs(aa1)+dabs(aa2)+dabs(aa3)
          subrat(i)=dabs(aa-area(nnel)/4)*4/area(nnel)
          if(subrat(i)<100*small1) then
            index=1
            sig(1)=aa1*4/area(nnel)
            sig(2)=aa2*4/area(nnel)
            sig(1)=dmax1(0.0d0,dmin1(1.0d0,sig(1)))
            sig(2)=dmax1(0.0d0,dmin1(1.0d0,sig(2)))
            if(sig(1)+sig(2)>1) then
              sig(3)=0
              sig(2)=1-sig(1)
            else
              sig(3)=1-sig(1)-sig(2)
            endif

!            t_xi(1)=tsd(jlev,n1)*(1-zrat)+tsd(jlev-1,n1)*zrat
!            t_xi(2)=tsd(jlev,n2)*(1-zrat)+tsd(jlev-1,n2)*zrat

!            smax=-99; tmin=100; ibb=0 !flag
            do jj=1,3 !side
              isd=js(nnel,jj)
              kbb=kbs(isd)
              swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
              swild2(kbb:nvrt,1)=tsd(kbb:nvrt,isd)
              swild2(kbb:nvrt,2)=ssd(kbb:nvrt,isd)
!              swild2(kbb:nvrt,3)=su2(kbb:nvrt,isd)
!              swild2(kbb:nvrt,4)=sv2(kbb:nvrt,isd)
              call vinter(mnv,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
!              if(ibelow==1) ibb=1
!              if(smax<swild(2)) smax=swild(2)
!              if(tmin>swild(1)) tmin=swild(1)
              t_xi(jj)=swild(1); s_xi(jj)=swild(2)

            enddo !jj=1,3

!            if(ibb==1) then
!              if(smax<-98) then
!                write(11,*)'Max. S failed to exist (2):',zt,nnel
!                stop
!              endif
!              t_xi(1:3)=tmin; s_xi(1:3)=smax
!            endif

            out6(1)=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
            out6(2)=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)

            exit
          endif !subrat(i)<100*small1
        endif !i<=3
      enddo !i=1,4

      if(index==0) then
        write(11,*)'Not in any sub-element',nnel,(subrat(i),i=1,4)
        write(11,*)xt,yt
        stop
      endif

!-----------------------------------------------------------------------
      else if(lqk(nnel)==2) then
!-----------------------------------------------------------------------
!     Quadratic interplation
!     For pure S only
      if(ss<-1.or.ss>0.or.kbpl/=kz) then
        write(11,*)'ss out of bound:',ss,kbpl
        stop
      endif

      if(wwint00>=0) then
        lin=-1 !lower interval
        if(jlev==kbpl+1) lin=-99
      else
        lin=1
        if(jlev==nvrt) lin=-98
      endif
      outq1=0; outq2=0 
      t_min=100; t_max=-100; s_min=100; s_max=-100

      do i=1,3 !nodes and sides
        nd=nm(nnel,i)
        isd=js(nnel,i)
        in1=nx(i,1)
        in2=nx(i,2)
!       check (range extended)
        if(tnd(jlev,nd)<-98.or.snd(jlev,nd)<-98.or.tsd(jlev,isd)<-98.or.ssd(jlev,isd)<-98) then
          write(11,*)'Wrong S,T:',i,nd,isd,tnd(jlev,nd),snd(jlev,nd),tsd(jlev,isd),ssd(jlev,isd)
          stop
        endif

        if(dabs(ss+1)<1.e-4.or.dabs(ss)<1.e-4) then !two surfaces
          if(dabs(ss+1)<1.e-4) then
            lev=kz
          else
            lev=nvrt
          endif
          t_n=tnd(lev,nd)
          s_n=snd(lev,nd)
          t_s=tsd(lev,isd)
          s_s=ssd(lev,isd)
          temp_min=dmin1(tnd(lev,nd),tsd(lev,isd))
          temp_max=dmax1(tnd(lev,nd),tsd(lev,isd))
          salt_min=dmin1(snd(lev,nd),ssd(lev,isd))
          salt_max=dmax1(snd(lev,nd),ssd(lev,isd))
        else if(lin<=-98) then !constrained bottom or surface
          if(lin==-99) then
!            zrat3=((zt-ztmp(kbpl))/(ztmp(kbpl+1)-ztmp(kbpl)))**2
            srat=((ss+1)/(sigma(2)+1))**2
          else 
!            zrat3=((zt-ztmp(nvrt))/(ztmp(nvrt-1)-ztmp(nvrt)))**2
!            zrat3=1-zrat3 !to put in same form
            srat=(ss/sigma(nsig-1))**2
            srat=1-srat !to put in same form
          endif
          if(srat<0.or.srat>1) then
            write(11,*)'Out of bound (9):',srat
            stop
          endif
          t_n=tnd(jlev,nd)*srat+tnd(jlev-1,nd)*(1-srat)
          s_n=snd(jlev,nd)*srat+snd(jlev-1,nd)*(1-srat)
          t_s=tsd(jlev,isd)*srat+tsd(jlev-1,isd)*(1-srat)
          s_s=ssd(jlev,isd)*srat+ssd(jlev-1,isd)*(1-srat)
          temp_min=dmin1(tnd(jlev,nd),tnd(jlev-1,nd),tsd(jlev,isd),tsd(jlev-1,isd))
          temp_max=dmax1(tnd(jlev,nd),tnd(jlev-1,nd),tsd(jlev,isd),tsd(jlev-1,isd))
          salt_min=dmin1(snd(jlev,nd),snd(jlev-1,nd),ssd(jlev,isd),ssd(jlev-1,isd))
          salt_max=dmax1(snd(jlev,nd),snd(jlev-1,nd),ssd(jlev,isd),ssd(jlev-1,isd))
        else !normal
          if(i==1) then !the following is indepdendent of loop i
            if(lin==1) then
              k1=jlev-1
            else !=-1
              k1=jlev-2
            endif           
            k2=k1+1; k3=k2+1
            if(k1<kbpl.or.k3>nvrt) then
              write(11,*)'Weird level:',k1,k2,k3
              stop
            endif
            k1s=k1-kz+1; k2s=k2-kz+1; k3s=k3-kz+1 !change to sigma indices
           
            denom=sigma(k3s)-2*sigma(k2s)+sigma(k1s)
            if(dabs(denom)<1.e-5) then !degenerate
              xi=2*(ss-sigma(k2s))/(sigma(k3s)-sigma(k1s))
            else
              del=(sigma(k3s)-sigma(k1s))**2-8*(sigma(k2s)-ss)*denom
              if(del<0) then
                write(11,*)'No inverse quadratic mapping:',del
                stop
              endif
              icount=0
              vzn(1)=(sigma(k1s)-sigma(k3s)+dsqrt(del))/2/denom !!temporary storage
              vzn(2)=(sigma(k1s)-sigma(k3s)-dsqrt(del))/2/denom
              xi_m=vzn(1) !for no root scenario
              if(dabs(vzn(2))<dabs(vzn(1))) xi_m=vzn(2)
              do l=1,2
                if(dabs(vzn(l))<=1+1.e-4) then
                  icount=icount+1
                  vxn(icount)=dmax1(-1.d0,dmin1(1.d0,vzn(l)))
                endif
              enddo !l
              if(icount==0) then
!                if(ifort12(15)==0) then
!                  ifort12(15)=1
                write(11,*)'No roots in inverse quadratic mapping:',(vzn(l),l=1,2)
                write(11,*)ss,sigma(k1s),sigma(k2s),sigma(k3s)
                stop
!                endif
!                xi=dmax1(-1.d0,dmin1(1.d0,xi_m))
              else if(icount==2) then
                if(ifort12(14)==0) then
                  ifort12(14)=1
                  write(12,*)'Warning: 2 roots:',(vxn(l),l=1,2)
                  write(12,*)ss,sigma(k1s),sigma(k2s),sigma(k3s),k1,k2,k3,kbpl
                endif
                xi=vxn(1)
              else !=1
                xi=vxn(1)
              endif
            endif !degenerate
           
            phi1=xi*(xi-1)/2; phi2=1-xi*xi; phi3=xi*(xi+1)/2
          endif !i==1

          t_n=tnd(k1,nd)*phi1+tnd(k2,nd)*phi2+tnd(k3,nd)*phi3
          s_n=snd(k1,nd)*phi1+snd(k2,nd)*phi2+snd(k3,nd)*phi3
          t_s=tsd(k1,isd)*phi1+tsd(k2,isd)*phi2+tsd(k3,isd)*phi3
          s_s=ssd(k1,isd)*phi1+ssd(k2,isd)*phi2+ssd(k3,isd)*phi3
          temp_min=dmin1(tnd(k1,nd),tnd(k2,nd),tnd(k3,nd),tsd(k1,isd),tsd(k2,isd),tsd(k3,isd))
          temp_max=dmax1(tnd(k1,nd),tnd(k2,nd),tnd(k3,nd),tsd(k1,isd),tsd(k2,isd),tsd(k3,isd))
          salt_min=dmin1(snd(k1,nd),snd(k2,nd),snd(k3,nd),ssd(k1,isd),ssd(k2,isd),ssd(k3,isd))
          salt_max=dmax1(snd(k1,nd),snd(k2,nd),snd(k3,nd),ssd(k1,isd),ssd(k2,isd),ssd(k3,isd))
        endif !normal

        outq1=outq1+t_n*(2*arco(i)*arco(i)-arco(i))+t_s*4*arco(in1)*arco(in2)
        outq2=outq2+s_n*(2*arco(i)*arco(i)-arco(i))+s_s*4*arco(in1)*arco(in2)
        if(temp_min<t_min) t_min=temp_min
        if(temp_max>t_max) t_max=temp_max
        if(salt_min<s_min) s_min=salt_min
        if(salt_max>s_max) s_max=salt_max
      enddo !i=1,3

      if(t_min>t_max) then
        write(11,*)'Illegal min/max for temp:',t_min,t_max,nnel
        stop
      endif
      if(s_min>s_max) then
        write(11,*)'Illegal min/max for salt:',s_min,s_max,nnel
        stop
      endif

      out6(1)=dmax1(t_min,dmin1(t_max,outq1))
      out6(2)=dmax1(s_min,dmin1(s_max,outq2))

!-----------------------------------------------------------------------
      endif !linear or quadratic


      return
      end

      function sums(x1,x2,x3,x4,y1,y2,y3,y4)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)

      sums=dabs((x4-x3)*(y2-y3)+(x2-x3)*(y3-y4))/2+&
     &     dabs((x4-x1)*(y3-y1)-(y4-y1)*(x3-x1))/2+&
     &     dabs((y4-y1)*(x2-x1)-(x4-x1)*(y2-y1))/2
      
      return
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      return
      end

      subroutine header
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      character(len=48) :: variable_out
!      integer :: ivs,i23d

      do i=1,noutput
        if(iwrite.eq.0) irec(i)=0 !record # for binary
        ichan(i)=100+i !output channel #
        if(i>=13.and.i<=15.or.i==25) then
          ivs=2
        else
          ivs=1
        endif
        if(i<=15) then
          i23d=2 !2 or 3D
        else
          i23d=3
        endif

        if(iof(i)==1) then
          if(iwrite.eq.0) then
            open(ichan(i),file=ifile_char//'_'//outfile(i),access='direct',recl=nbyte)
!	    ' (ylz)
          else !evm
            open(ichan(i),file=ifile_char//'_'//outfile(i))
          endif

          if(iwrite.eq.0) then
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) data_format(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) version(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) start_time(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            variable_out=variable_nm(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            variable_out=variable_dim(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte

            write(ichan(i),rec=irec(i)+1) nrec
            write(ichan(i),rec=irec(i)+2) real(dt*nspool)
            write(ichan(i),rec=irec(i)+3) nspool
            write(ichan(i),rec=irec(i)+4) ivs
            write(ichan(i),rec=irec(i)+5) i23d
            irec(i)=irec(i)+5

          else !evm
            write(ichan(i),'(a48)',advance="no") data_format
            write(ichan(i),'(a48)',advance="no") version
            write(ichan(i),'(a48)',advance="no") start_time
            write(ichan(i),'(a48)',advance="no") variable_nm(i)
            write(ichan(i),'(a48)',advance="no") variable_dim(i)
            a_4 = transfer(source=nrec,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=real(dt*nspool),mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=nspool,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=ivs,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=i23d,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
          endif

!	  Vertical grid
          if(iwrite.eq.0) then
            write(ichan(i),rec=irec(i)+1) nvrt
            write(ichan(i),rec=irec(i)+2) kz
            write(ichan(i),rec=irec(i)+3) real(h0)
            write(ichan(i),rec=irec(i)+4) real(h_s)
            write(ichan(i),rec=irec(i)+5) real(h_c)
            write(ichan(i),rec=irec(i)+6) real(theta_b)
            write(ichan(i),rec=irec(i)+7) real(theta_f)
            irec(i)=irec(i)+7
          else !evm
            a_4 = transfer(source=nvrt,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=kz,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(h0),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(h_s),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(h_c),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(theta_b),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(theta_f),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
          endif
          do k=1,kz-1
            if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+k) real(ztot(k))
            else !evm
              a_4 = transfer(source=real(ztot(k)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
           endif
          enddo
          do k=kz,nvrt
            kin=k-kz+1
            if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+k) real(sigma(kin))
            else !evm
              a_4 = transfer(source=real(sigma(kin)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
           endif
          enddo
          if(iwrite.eq.0) irec(i)=irec(i)+nvrt
          irecm=48/nbyte*5+5+7+nvrt !estimates of total record #

!	  Horizontal grid
          if(iwrite.eq.0) then
            write(ichan(i),rec=irec(i)+1) np
            write(ichan(i),rec=irec(i)+2) ne
            irec(i)=irec(i)+2
          else !evm
            a_4 = transfer(source=np,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=ne,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
          endif
          do m=1,np
            if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+1)real(x(m))
              write(ichan(i),rec=irec(i)+2)real(y(m))
              write(ichan(i),rec=irec(i)+3)real(dp00(m))
              write(ichan(i),rec=irec(i)+4)kbp00(m)
              irec(i)=irec(i)+4
            else !evm
              a_4 = transfer(source=real(x(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(y(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(dp00(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=kbp00(m),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
            endif
          enddo !m=1,np
          do m=1,ne
            if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+1) 3
              irec(i)=irec(i)+1
              do mm=1,3
                write(ichan(i),rec=irec(i)+mm)nm(m,mm)
              enddo !mm
              irec(i)=irec(i)+3
            else !evm
              a_4 = transfer(source=3,mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              do mm=1,3
                a_4 = transfer(source=nm(m,mm),mold=a_4)
                write(ichan(i),'(a4)',advance="no") a_4
              enddo !mm
            endif
          enddo !m

!	  Estimate total # of records
          irecm=irecm+3*np+4*ne
          if(i23d.eq.2) then !2D
            irecm=irecm+(2+np+np*ivs)*nrec
          else !3D
            irecm=irecm+(2+np+np*nvrt*ivs)*nrec
          endif

          if(irecm.gt.mirec) then
            write(11,*)'Output file too large',i,irecm
            stop
          endif
          
          if(iwrite.eq.0) then
            close(ichan(i)) ! do this to flush the write buffer
            open(ichan(i),file=ifile_char//'_'//outfile(i),access='direct',recl=nbyte)
!	    ' (ylz)
          endif
        endif !iof(i)=1
      enddo !i=1,noutput

!     Test output (old vis5 format)
      igmp=0
      if(noutgm.eq.1) then
        open(100,file=ifile_char//'_test.60',access='direct',recl=nbyte)
        igmp=(32+24+24)/nbyte
        write(100,rec=igmp+1) nrec
        write(100,rec=igmp+2) ns
        write(100,rec=igmp+3) real(dt*nspool)
        write(100,rec=igmp+4) nspool
        write(100,rec=igmp+5) 2
        igmp=igmp+5
        close(100)
        open(100,file=ifile_char//'_test.60',access='direct',recl=nbyte)
      endif

      return
      end



!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*dcos(phi0)
      y=phi*r

      return
      end

!
!********************************************************************************
!										*
!     Straightline search algorithm. Initially nnel0 is an element that 	*
!     encompasses (x0,y0). iloc=0: do not nudge initial pt; iloc=1: nudge.	* 
!     Input: iloc,nnel0,x0,y0,z0,xt,yt,zt,jlev0, time, and su2,sv2,ww2 for 	*
!	abnormal cases;								*
!     Output the updated end pt (xt,yt,zt) (if so), nnel1, jlev1, area          *
!       coordinates, vertical ratio and a flag nfl.				*
!     nfl=1 if a bnd or dry element is hit and vel. there is small,		* 
!	or death trap is reached.						*
!										*
!********************************************************************************
!
      subroutine quicksearch(iloc,idt,id0,nnel0,jlev0,time,x0,y0,z0,xt,yt,zt,nnel1,jlev1,arco,zrat,ztmp,kbpl,nfl,ss)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iloc,idt,id0,nnel0,jlev0
      real(kind=dbl_kind), intent(in) :: time,x0,y0,z0
      integer, intent(out) :: nnel1,jlev1,nfl
      real(kind=dbl_kind), intent(inout) :: xt,yt,zt
      integer, intent(out) :: kbpl
      real(kind=dbl_kind), intent(out) :: arco(3),zrat,ztmp(mnv),ss

      dimension out2(mnv)

      if(iloc>1) then
        write(11,*)'iloc > 1'
        stop
      endif
      if(idry_e(nnel0)==1) then
        write(11,*)'Starting element is dry'
        stop
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
!     Try area_coord
      nel=nnel0
      aa=0
      aa1=0
      do i=1,3
        n1=nm(nel,i)
        n2=nm(nel,nx(i,1))
        aa=aa+dabs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
        aa1=aa1+dabs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
      enddo !i
      ae=dabs(aa-area(nel))/area(nel)
      if(ae>small1) then
        write(11,*)'(x0,y0) not in nnel0 initially',ae,nnel0
        stop
      endif

      ae=dabs(aa1-area(nel))/area(nel)
      if(ae<small1) then
        nnel1=nel
        go to 400
      endif

!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!     An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
      if(iloc==0) then
        xcg=x0
        ycg=y0
      else if(iloc==1) then
!        weit=1./3; al=0; bet=0
!        xint=x(nm(nel,1))*(weit-al)+x(nm(nel,2))*(weit-bet)+x(nm(nel,3))*(weit+al+bet)
!        yint=y(nm(nel,1))*(weit-al)+y(nm(nel,2))*(weit-bet)+y(nm(nel,3))*(weit+al+bet)
!        xcg=(1-5.e-4)*x0+5.e-4*xint
!        ycg=(1-5.e-4)*y0+5.e-4*yint
        xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
        ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)
      endif

      pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
      if(pathl==0) then
        write(11,*)'Zero path',x0,y0,xt,yt,xcg,ycg
        stop
      endif

!     Starting edge nel_j
      nel_j=0
      do j=1,3
        jd1=nm(nel,nx(j,1))
        jd2=nm(nel,nx(j,2))
        call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag==1) then
          nel_j=j
          exit
        endif
      enddo !j=1,3
      if(nel_j==0) then
        write(11,*)'Found no intersecting edges I:',nel,xcg,ycg,xt,yt,ae
        stop
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it>1000) then
        if(ifort12(3)==0) then
          ifort12(3)=1
          write(12,*)'Death trap reached',idt,id0
        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel1=nel
        exit loop4
      endif
      md1=nm(nel,nx(nel_j,1))
      md2=nm(nel,nx(nel_j,2))
      
!     Compute z position 
      dist=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(dist/pathl.gt.1+1.0d-4) then
        write(11,*)'Path overshot'
        stop
      endif
      zin=zt-dist/pathl*(zt-zin)
      trm=trm*dist/pathl !time remaining
      
      pathl=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl==0.or.trm==0) then
        write(11,*)'Target reached'
        stop
      endif

      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
      if(ic3(nel,nel_j)==0.or.idry_e(ic3(nel,nel_j))==1) then
        lit=1
        isd=js(nel,nel_j)
        if(isidenode(isd,1)+isidenode(isd,2)/=md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

!       Nudge intersect (xin,yin), and update starting pt
        xin=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
        yin=(1-1.0d-4)*yin+1.0d-4*yctr(nel)
        xcg=xin
        ycg=yin

        vtan=-su2(jlev0,isd)*sny(isd)+sv2(jlev0,isd)*snx(isd)
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)
        zvel=(ww2(jlev0,md1)+ww2(jlev0,md2))/2
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=dsqrt(xvel**2+yvel**2)
        if(hvel<1.e-4) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel1=nel
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel,nel_j) !next front element
      aa=0
      do i=1,3
        k1=nm(nel,i)
        k2=nm(nel,nx(i,1))
        aa=aa+dabs(signa(x(k1),x(k2),xt,y(k1),y(k2),yt))
      enddo !i
      ae=dabs(aa-area(nel))/area(nel)
      if(ae<small1) then
        nnel1=nel
        exit loop4
      endif

!     Next intersecting edge
      do j=1,3
         jd1=nm(nel,nx(j,1))
         jd2=nm(nel,nx(j,2))
!        For abnormal case, same side (border side) cannot be hit again
         if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
         call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag==1) then
           nel_j=j !next front edge          
           cycle loop4
         endif
      enddo !j
      write(11,*)'Failed to find next edge I:',lit,xin,yin,xt,yt,nel,md1,md2,idt,id0,ae
      stop
!----------------------------------------------------------------------------------------
      end do loop4 

400   continue
!     No vertical exit from domain
      if(idry_e(nnel1)==1) then
        write(11,*)'Ending element is dry'
        stop
      endif

!     Compute area & sigma coord.
      call area_coord(nnel1,xt,yt,arco)
      n1=nm(nnel1,1)
      n2=nm(nnel1,2)
      n3=nm(nnel1,3)
      etal=eta2(n1)*arco(1)+eta2(n2)*arco(2)+eta2(n3)*arco(3)
      dep=dp(n1)*arco(1)+dp(n2)*arco(2)+dp(n3)*arco(3)
      if(etal+dep<h0) then
        write(11,*)'Weird wet element in quicksearch:',nnel1,eta2(n1),eta2(n2),eta2(n3)
        stop
      endif

!     Compute z-coordinates
      do k=kz,nvrt
        kin=k-kz+1
        hmod2=dmin1(dep,h_s)
        if(hmod2<=h_c) then
          ztmp(k)=sigma(kin)*(hmod2+etal)+etal
        else if(etal<=-h_c-(dep-h_c)*theta_f/s_con1) then
          write(11,*)'Pls choose a larger h_c (2):',etal,h_c
          stop
        else
          ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
        endif

!       Following to prevent underflow
        if(k==kz) ztmp(k)=-hmod2
        if(k==nvrt) ztmp(k)=etal
      enddo !k

      if(dep<=h_s) then
        kbpl=kz
      else !z levels
!       Find bottom index
        kbpl=0
        do k=1,kz-1
          if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
            kbpl=k
            exit
          endif
        enddo !k
        if(kbpl==0) then
          write(11,*)'Cannot find a bottom level at foot:',dep
          stop
        endif
        ztmp(kbpl)=-dep
        do k=kbpl+1,kz-1
          ztmp(k)=ztot(k)
        enddo !k
      endif

      do k=kbpl+1,nvrt
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(11,*)'Inverted z-level in quicksearch:',nnel1,etal,dep,ztmp(k)-ztmp(k-1)
          stop
        endif
      enddo !k

      if(zt<=ztmp(kbpl)) then
        zt=ztmp(kbpl)
        zrat=1
        jlev1=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0
        jlev1=nvrt
      else
        jlev1=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev1=k+1
            exit
          endif
        enddo !k
        if(jlev1==0) then
          write(11,*)'Cannot find a vert. level:',zt,etal,dep
          write(11,*)(ztmp(k),k=kbpl,nvrt)
          stop
        endif
        zrat=(ztmp(jlev1)-zt)/(ztmp(jlev1)-ztmp(jlev1-1))
      endif

      if(zrat<0.or.zrat>1) then
        write(11,*)'Sigma coord. wrong (4):',jlev1,zrat
        stop
      endif

      if(kbpl==kz) then !in pure S region
        ss=(1-zrat)*sigma(jlev1-kz+1)+zrat*sigma(jlev1-kz)
      else
        ss=-99
      endif


!      if(ss<sigma(jlev1-1).or.ss>sigma(jlev1)) then
!        write(11,*)'Sigma coord. wrong (5):',jlev1,ss,sigma(jlev1-1),sigma(jlev1)
!        stop
!      endif

      return
      end


!-----------------------------------------------------------------------------------------------!
!												!
!       Find first wet upwind element for rewetted points, and compute S,T. 			!
!       Use idry_e0 (from previous step) and S,T from previous step.				!
!       Abnormal exit: death trap or a horizontal bnd is reached and vel.=0 there.		!
!												!	
!-----------------------------------------------------------------------------------------------!
!
      subroutine upwindtrack(id0,jlev,nnel,x0,y0,uuint,vvint,out2,nfl)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: id0,jlev,nnel
      real(kind=dbl_kind), intent(in) :: x0,y0
      real(kind=dbl_kind), intent(inout) :: uuint,vvint
      integer, intent(out) :: nfl
      real(kind=dbl_kind), intent(out) :: out2(12)


      nfl=0 !trap flag
      
!     Starting element nel (must be dry)
      nel=nnel

!     An interior pt close to (x0,y0) to prevent underflow
      xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
      ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)

!     Starting edge nel_j
      nel_j=0
      do j=1,3 !side
        jd1=nm(nel,nx(j,1))
        jd2=nm(nel,nx(j,2))
        call intersect3(xcg,uuint,x(jd1),x(jd2),ycg,vvint,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag==1) then
          nel_j=j
          exit
        endif
      enddo !j=1,3
      if(nel_j==0) then
        write(11,*)'Found no intersecting edges III:',nel,xcg,ycg,uuint,vvint
        stop
      endif

      it=0
      loop7: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it>50) then
        if(ifort12(9)==0) then
          ifort12(9)=1
          write(12,*)'Death trap reached in upwindtrack',id0
        endif
        nfl=1
        exit loop7
      endif
      md1=nm(nel,nx(nel_j,1))
      md2=nm(nel,nx(nel_j,2))
      
      lit=0 !flag
!     For horizontal exit, compute tangential vel, update target uunit & vvint, and continue.
      if(ic3(nel,nel_j)==0) then !.or.idry_e(ic3(nel,nel_j))==1) then
        lit=1
        isd=js(nel,nel_j)
        if(isidenode(isd,1)+isidenode(isd,2)/=md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

!       Nudge intersect (xin,yin), and update starting pt
        xin=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
        yin=(1-1.0d-4)*yin+1.0d-4*yctr(nel)
        xcg=xin
        ycg=yin

        vtan=-su2(jlev,isd)*sny(isd)+sv2(jlev,isd)*snx(isd)
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)
        hvel=dsqrt(xvel**2+yvel**2)
        if(hvel==0) then
          nfl=1
          exit loop7
        endif
        uuint=xvel/hvel
        vvint=yvel/hvel
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel,nel_j) !next front element
      if(idry_e0(nel)==0) then
        isd0=0
        do j=1,3 !side
          isd=js(nel,j)
          n1=isidenode(isd,1)
          n2=isidenode(isd,2)
          if(n1==md1.and.n2==md2.or.n1==md2.and.n2==md1) isd0=isd
        enddo
        if(isd0==0) then
          write(11,*)'Wrong connectivity (3)'
          stop
        endif
        if(tsd(jlev,isd0)<-98.or.ssd(jlev,isd0)<-98) then
          write(11,*)'Impossible dry (9)'
          stop
        endif
        out2(1)=tsd(jlev,isd0)
        out2(2)=ssd(jlev,isd0)
        exit loop7
      endif

!     Next intersecting edge
      do j=1,3
         jd1=nm(nel,nx(j,1))
         jd2=nm(nel,nx(j,2))
         if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
         call intersect3(xcg,uuint,x(jd1),x(jd2),ycg,vvint,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag==1) then
           nel_j=j !next front edge          
           cycle loop7
         endif
      enddo !j
      write(11,*)'Failed to find next edge (4):',lit,xin,yin,nel,md1,md2,id0
      stop
!----------------------------------------------------------------------------------------
      end do loop7 

      return
      end

!
!********************************************************************************
!										*
!     Program to detect if an infinite line (3,4) and a semi-infinite line 	*
!     from pt 1 in (uuint,vvint) have common pts   				*
!     Assumption: the 3 pts are distinctive.					*
!     The eq. for the infinite line is: X=X3+(X4-X3)*tt2 (-\infty <tt2 <\infty).		*
!     The eq. for the semi-infinite line is: X=X1+(uuint,vvint)*tt1 (tt1<=0).	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection, tt1, tt2 can also be used.     *
!										*
!********************************************************************************
!

      subroutine intersect4(x1,uuint,x3,x4,y1,vvint,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: small2=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,uuint,x3,x4,y1,vvint,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=uuint*(y3-y4)-vvint*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=uuint*(y3-y1)-vvint*(x3-x1)

      if(delta/=0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1<=small2) then
          iflag=1
          xin=x3+(x4-x3)*tt2
          yin=y3+(y4-y3)*tt2
        endif
      endif

      return
      end

!
!********************************************************************************
!										*
!     Program to detect if a finite segment (3,4) and a semi-infinite line 	*
!     from pt 1 in (uuint,vvint) have common pts   				*
!     Assumption: the 3 pts are distinctive.					*
!     The eq. for the finite line is: X=X3+(X4-X3)*tt2 (0<=tt2<=1).		*
!     The eq. for the semi-infinite line is: X=X1+(uuint,vvint)*tt1 (tt1<=0).	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection, tt1, tt2 can also be used.     *
!										*
!********************************************************************************
!

      subroutine intersect3(x1,uuint,x3,x4,y1,vvint,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: small2=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,uuint,x3,x4,y1,vvint,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=uuint*(y3-y4)-vvint*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=uuint*(y3-y1)-vvint*(x3-x1)

      if(delta/=0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1<=small2.and.tt2>=-small2.and.tt2<=1+small2) then
          iflag=1
          xin=x3+(x4-x3)*tt2
          yin=y3+(y4-y3)*tt2
        endif
      endif

      return
      end

!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!										*
!********************************************************************************
!

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: small2=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta/=0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1>=-small2.and.tt1<=1+small2.and.tt2>=-small2.and.tt2<=1+small2) then
          iflag=1
          xin=x1+(x2-x1)*tt1
          yin=y1+(y2-y1)*tt1
        endif
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
!	Compute area coordinates of pt (xt,yt), which must be inside element nnel.	!
!	Impose bounds for area coordinates.						!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine area_coord(nnel,xt,yt,arco)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      real(kind=dbl_kind), intent(in) :: xt,yt
      real(kind=dbl_kind), intent(out) :: arco(3)

      n1=nm(nnel,1)
      n2=nm(nnel,2)
      n3=nm(nnel,3)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(1)=dmax1(0.0d0,dmin1(1.0d0,arco(1)))
      arco(2)=dmax1(0.0d0,dmin1(1.0d0,arco(2)))
      if(arco(1)+arco(2)>1) then
        arco(3)=0
        arco(1)=dmin1(1.d0,dmax1(0.d0,arco(1)))
        arco(2)=1-arco(1)
      else
        arco(3)=1-arco(1)-arco(2)
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!       Compute area coordinates of pt (xt,yt), which may not be inside element nnel.      !
!       ifl=0: inside; =1: outside.
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine area_coord2(nnel,xt,yt,arco,ifl)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      integer, intent(out) :: ifl
      real(kind=dbl_kind), intent(in) :: xt,yt
      real(kind=dbl_kind), intent(out) :: arco(3)
    
      ifl=0
      n1=nm(nnel,1)
      n2=nm(nnel,2)
      n3=nm(nnel,3)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(3)=1-arco(1)-arco(2)
      if(arco(1)<0.or.arco(1)>1.or.arco(2)<0.or.arco(2)>1.or.arco(3)<0.or.arco(3)>1) ifl=1

      return
      end

!
!***************************************************************************
!									   *
!     Solve for the density at nodes and sides.				   *
!     From Pond and Pickard's book.					   *
!     validity region: T: [0,40], S: [0:42]				   *
!									   *
!***************************************************************************
!   

      subroutine eqstate
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      den(t,s)=1000-0.157406+6.793952E-2*t-9.095290E-3*t**2 &
     &+1.001685E-4*t**3-1.120083E-6*t**4+6.536332E-9*t**5+ &
     &s*(0.824493-4.0899E-3*t+&
     &7.6438E-5*t**2-8.2467E-7*t**3+5.3875E-9*t**4)+&
     &dsqrt(s)**3*(-5.72466E-3+1.0227E-4*t-1.6546E-6*t**2)+&
     &4.8314E-4*s**2

      prho=-99 !flags
!      do l=1,2 !nodes & sides 
      do l=1,1 !nodes only
        if(l==1) then
          limit=np
        else 
          limit=ns
        endif
        do i=1,limit
          if(l==1.and.idry(i)==1.or.l==2.and.idry_s(i)==1) cycle

!         Valid S,T
          do k=1,nvrt
            if(l==1) then !node
              ttmp=tnd(k,i)
              stmp=snd(k,i)
            else !side
              ttmp=tsd(k,i)
              stmp=ssd(k,i)
            endif

            if(ttmp<-98.or.stmp<-98) then
              write(11,*)'Impossible dry (7):',l,i,k,ttmp,stmp
              stop
            endif
            if(ttmp<tempmin.or.ttmp>tempmax.or.stmp<saltmin.or.stmp>saltmax) then
              if(ifort12(6)==0) then
                ifort12(6)=1
                write(12,*)'Invalid temp. or salinity for density'
                write(12,*)ttmp,stmp,l,i,k
              endif
              ttmp=dmax1(tempmin,dmin1(ttmp,tempmax))
              stmp=dmax1(saltmin,dmin1(stmp,saltmax))
            endif
!	    Density at one standard atmosphere
            rho=den(ttmp,stmp)

            if(rho<980) then
              write(11,*)'Weird density at:',l,i,k,rho,ttmp,stmp
              stop
            endif
            if(l==1) then
              prho(i,k)=rho
              sig_t(i,k)=rho-rho0
            else if(l==2) then
!              srho(i,k)=rho
            endif
          enddo !k=1,nvrt
        enddo !i=1,limit
      enddo !l=1,1

      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
!    Generic routine to compute \int_{\sigma_k}^{\sigma_{k+1}} \psi(\sigma)d\sigma,	!
!    where Nmin<=k<=Nmax-1, \sigma & \psi(Nmin:Nmax), using Lagrangian  		!
!    interpolation of order 2*m (i.e., from k-m to k+m).				!
!    mnv: dimensioning parameter from driving routine (input);				!
!    Nmin, Nmax: limits of vertical levels (input);					!
!    m: order of Lagrangian polynormial (input);					!
!    k: input for limits;								!
!    sigma,sigmap,sigma_prod,psi: input (sigmap&sigma_prod are the pre-computed 	!
!                                  powers and products of sigma for speed)		!
!    gam, coef: working arrays (output).						!
!    WARNING: Nmax must =nsig, and 1<=Nmin<=nsig-1 for sigma_prod!!			!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function rint_lag(mnv,Nmin,Nmax,m,k,sigma,sigmap,sigma_prod,psi,gam,coef)
      implicit real*8(a-h,o-z)
      integer, intent(in) :: mnv,Nmin,Nmax,m,k
      real(kind=8), intent(in) :: sigma(mnv),sigmap(mnv,10),sigma_prod(mnv,mnv,-4:4),psi(mnv)
      real(kind=8), intent(out) :: gam(mnv),coef(0:mnv)

!     Sanity check
      if(Nmin>=Nmax.or.Nmax>mnv.or.Nmin<1) then
        write(11,*)'Check inputs in rint_lag:',Nmin,Nmax
        stop
      endif
      if(k>Nmax-1.or.k<Nmin) then
        write(11,*)'Wrong k:',k
        stop
      endif
      if(m<1) then
        write(11,*)'m<1',m
        stop
      endif
      if(m>3) then
        write(*,*)'m>3 not covered presently' 
        stop
      endif
      if(2*m+1>10) then
        write(11,*)'Re-dimension sigmap'
        stop
      endif

!     Compute J1,2
      j1=max0(Nmin,k-m)
      j2=min0(Nmax,k+m)
      if(j1>=j2) then
         write(11,*)'Weird indices:',j1,j2
         stop
      endif

!     Compute sum
      rint_lag=0
      do i=j1,j2
!       Denominator & assemble working array gam
!        prod=1
        id=0
        do j=j1,j2
          if(j/=i) then
            id=id+1
            gam(id)=-sigma(j)
          endif
        enddo !j
        if(id/=j2-j1.or.id>2*m) then
          write(11,*)'Miscount:',id,j2-j1,m
          stop
        endif

!       Inner sum
        if(id==1) then
          coef(0)=gam(1); coef(1)=1
        else if(id==2) then
          coef(0)=gam(1)*gam(2)
          coef(1)=gam(1)+gam(2)
          coef(2)=1
        else if(id==3) then
          coef(0)=gam(1)*gam(2)*gam(3)
          coef(1)=gam(1)*(gam(2)+gam(3))+gam(2)*gam(3)
          coef(2)=gam(1)+gam(2)+gam(3)
          coef(3)=1
        else if(id==4) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)
          coef(1)=gam(1)*gam(2)*(gam(3)+gam(4))+(gam(1)+gam(2))*gam(3)*gam(4)
          coef(2)=gam(1)*(gam(2)+gam(3))+(gam(1)+gam(3))*gam(4)+gam(2)*(gam(3)+gam(4))
!          coef(2)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(2)*gam(3)+gam(2)*gam(4)+gam(3)*gam(4)
          coef(3)=gam(1)+gam(2)+gam(3)+gam(4)
          coef(4)=1
        else if(id==5) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(4)*gam(5)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(5)
          coef(2)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+gam(1)*gam(3)*gam(4)+ &
     &gam(1)*gam(3)*gam(5)+gam(1)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)+gam(2)*gam(3)*gam(5)+ &
     &gam(2)*gam(4)*gam(5)+gam(3)*gam(4)*gam(5)
          coef(3)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(2)*gam(3)+ &
     &gam(2)*gam(4)+gam(2)*gam(5)+gam(3)*gam(4)+gam(3)*gam(5)+gam(4)*gam(5)
          coef(4)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)
          coef(5)=1
        else if(id==6) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)+gam(1)*gam(2)*gam(3)*gam(4)*gam(6)+&
     &gam(1)*gam(2)*gam(3)*gam(5)*gam(6)+gam(1)*gam(2)*gam(4)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(2)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(3)*gam(6)+ &
     &gam(1)*gam(2)*gam(4)*gam(5)+gam(1)*gam(2)*gam(4)*gam(6)+gam(1)*gam(2)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(1)*gam(3)*gam(4)*gam(6)+gam(1)*gam(3)*gam(5)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(6)+ &
     &gam(2)*gam(3)*gam(5)*gam(6)+gam(2)*gam(4)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)*gam(6)
           coef(3)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+ &
     &gam(1)*gam(2)*gam(6)+gam(1)*gam(3)*gam(4)+gam(1)*gam(3)*gam(5)+gam(1)*gam(3)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)+gam(1)*gam(4)*gam(6)+gam(1)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)+ &
     &gam(2)*gam(3)*gam(5)+gam(2)*gam(3)*gam(6)+gam(2)*gam(4)*gam(5)+gam(2)*gam(4)*gam(6)+ &
     &gam(2)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)+gam(3)*gam(4)*gam(6)+gam(3)*gam(5)*gam(6)+ &
     &gam(4)*gam(5)*gam(6)
           coef(4)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(1)*gam(6)+ &
     &gam(2)*gam(3)+gam(2)*gam(4)+gam(2)*gam(5)+gam(2)*gam(6)+gam(3)*gam(4)+gam(3)*gam(5)+ &
     &gam(3)*gam(6)+gam(4)*gam(5)+gam(4)*gam(6)+gam(5)*gam(6)
           coef(5)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)+gam(6)
           coef(6)=1
        else
          write(*,*)'Not covered:',id
          stop
        endif

        sum=0
        do l=0,id
          sum=sum+coef(l)/(l+1)*(sigmap(k+1,l+1)-sigmap(k,l+1))
        enddo !l

        if(abs(i-k)>4) then
          write(11,*)'sigma_prod index out of bound (2)'
          stop
        endif

        rint_lag=rint_lag+psi(i)/sigma_prod(Nmin,k,i-k)*sum
      enddo !i=j1,j2

      return
      end

!
!***************************************************************************
!								     	   *
!     Convert normal vel. to 3D nodal vel. at WHOLE levels.                *
!								     	   *
!***************************************************************************
!
      subroutine nodalvel(ifltype)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: ifltype(mnope)

      dimension swild(10),swild2(mnv,10),swild3(mnv) !swild2's dimension must match vinter()

!     Compute discontinuous hvel first (used in btrack)
      ufg=0; vfg=0
      do i=1,ne
        do k=1,nvrt
          do j=1,3
            nd=nm(i,j)
            isd0=js(i,j)
            isd1=js(i,nx(j,1))
            isd2=js(i,nx(j,2))
            ufg(k,i,j)=su2(k,isd1)+su2(k,isd2)-su2(k,isd0)
            vfg(k,i,j)=sv2(k,isd1)+sv2(k,isd2)-sv2(k,isd0)
!           Error: impose bounds for ufg, vfg?
          enddo !j
        enddo !k
      enddo !i=1,ne

      if(indvel==0) then
!-------------------------------------------------------------------------------
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
      do i=1,np
        if(idry(i)==1) cycle

!       Wet node
        do k=kbp(i),nvrt
          weit_w=0
          icount=0
          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)
            if(idry_e(ie)==0) then
              icount=icount+1
              uu2(k,i)=uu2(k,i)+ufg(k,ie,id)
              vv2(k,i)=vv2(k,i)+vfg(k,ie,id)
            endif
 
!              if(.not.(isbnd(i)>0.and.ifltype(isbnd(i))/=0.and.isbs(isd)/=isbnd(i))) then

            if(interpol(ie)==1) then !along Z
              if(idry_e(ie)==1) then
                swild(1)=0
              else !wet eleemnt; node i is also wet
                kbb=kbe(ie)
                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
                call vinter(mnv,1,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
              endif
            else !along S
              swild(1)=we(k,ie)
            endif !Z or S

            ww2(k,i)=ww2(k,i)+swild(1)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j
          if(icount==0) then
            write(11,*)'Isolated wet node (8):',i
            stop
          else
            uu2(k,i)=uu2(k,i)/icount
            vv2(k,i)=vv2(k,i)/icount
          endif
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np

!-------------------------------------------------------------------------------
      else !indvel=1: averaging vel.
!-------------------------------------------------------------------------------
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
      do i=1,np
        if(idry(i)==1) cycle

!       Wet node
        icase=2
        do j=1,nne(i)
          ie=ine(i,j)
          if(interpol(ie)==1) icase=1
        enddo !j

        do k=kbp(i),nvrt
          weit=0
          icount=0
          weit_w=0
          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)

            if(isbnd(i)/=0) then !bnd ball
              limit=1
            else !internal ball
              limit=2
            endif
            do l=2,limit,-1
              isd=js(ie,nx(id,l))
              if(.not.(isbnd(i)>0.and.ifltype(isbnd(i))/=0.and.isbs(isd)/=isbnd(i))) then
                if(icase==1) then !along Z
                  if(idry_s(isd)==1) then
                    swild(1:2)=0
                  else !wet side; node i is also wet
                    kbb=kbs(isd)
                    swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)
                    swild2(kbb:nvrt,2)=sv2(kbb:nvrt,isd)
                    swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
                    call vinter(mnv,2,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
                  endif
                else !along S
                  swild(1)=su2(k,isd)
                  swild(2)=sv2(k,isd)
                endif !Z or S

                icount=icount+1
                uu2(k,i)=uu2(k,i)+swild(1)/distj(isd)
                vv2(k,i)=vv2(k,i)+swild(2)/distj(isd)
                weit=weit+1/distj(isd)
              endif
            enddo !l

            if(interpol(ie)==1) then !along Z
              if(idry_e(ie)==1) then
                swild(1)=0
              else !wet eleemnt; node i is also wet
                kbb=kbe(ie)
                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
                call vinter(mnv,1,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
              endif
            else !along S
              swild(1)=we(k,ie)
            endif !Z or S

            ww2(k,i)=ww2(k,i)+swild(1)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j
          if(icount==0) then
            write(11,*)'Isolated open bnd node:',i,isbnd(i)
            stop
          endif
          uu2(k,i)=uu2(k,i)/weit
          vv2(k,i)=vv2(k,i)/weit
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np
!-------------------------------------------------------------------------------
      endif !discontinous or averaging vel.

!...  Compute discrepancy between avergaed and elemental vel. vectors 
!      do i=1,np
!	do k=1,nvrt
!	  testa(i,k)=0
!          do j=1,nne(i)
!	    iel=ine(i,j)
!	    index=0
!	    do l=1,3
!	      if(nm(iel,l).eq.i) index=l
!	    enddo !l
!	    if(index.eq.0) then
!	      write(*,*)'Wrong element ball'
!	      stop
!	    endif
!	    testa(i,k)=testa(i,k)+dsqrt((uuf(iel,index,k)-uu2(k,i))**2+
!     +(vvf(iel,index,k)-vv2(k,i))**2)/nne(i)
!	  enddo !j
!	enddo !k
!      enddo !i

      return
      end

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!	Algebraic Stress Models									!
!												!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine asm(g,i,j,vd,td,qd1,qd2)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: i,j
      real(kind=dbl_kind), intent(in) :: g
      real(kind=dbl_kind), intent(out) :: vd,td,qd1,qd2

      if(j<kbp(i).or.j>nvrt) then
        write(11,*)'Wrong input level:',j
        stop
      endif

!     Wet node i with prho defined; kbp(i)<=j<=nvrt
      if(j==kbp(i).or.j==nvrt) then
        drho_dz=0
      else
        drho_dz=(prho(i,j+1)-prho(i,j-1))/(z(j+1,i)-z(j-1,i))
      endif
      bvf=g/rho0*drho_dz
      Gh=xl(i,j)**2/2/q2(i,j)*bvf
      Gh=dmin1(dmax1(Gh,-0.28d0),0.0233d0)

      if(stab.eq.'GA') then
        sh=0.49393/(1-34.676*Gh)
        sm=(0.39327-3.0858*Gh)/(1-34.676*Gh)/(1-6.1272*Gh)
        cmiu=dsqrt(2.d0)*sm
        cmiup=dsqrt(2.d0)*sh
        cmiu1=dsqrt(2.d0)*0.2 !for k-eq
        cmiu2=dsqrt(2.d0)*0.2 !for psi-eq.
      else if(stab.eq.'KC') then !Kantha and Clayson
!       Error: Warner's paper wrong
!        Ghp=(Gh-(Gh-0.02)**2)/(Gh+0.0233-0.04) !smoothing
        Ghp=Gh
        sh=0.4939/(1-30.19*Ghp)
        sm=(0.392+17.07*sh*Ghp)/(1-6.127*Ghp)
        cmiu=dsqrt(2.d0)*sm
        cmiup=dsqrt(2.d0)*sh
        cmiu1=cmiu/schk
        cmiu2=cmiu/schpsi
      else
        write(11,*)'Unknown ASM:',mid
        stop
      endif

      vd=cmiu*xl(i,j)*dsqrt(q2(i,j))
      td=cmiup*xl(i,j)*dsqrt(q2(i,j))
      qd1=cmiu1*xl(i,j)*dsqrt(q2(i,j))
      qd2=cmiu2*xl(i,j)*dsqrt(q2(i,j))

      return
      end

!     Routine to do vertical interpolation in z
!     Inputs:
!       mnv: dimensioning paramter for za etc.
!       nc: actual # of variables
!       k1,k2: lower and upper limits for za, sint
!       k3: initial guess for level index (to speed up)
!       zt: desired interpolation level
!       za(k1:k2): z-cor for sint (must be in ascending order)
!       sint(k1:k2,1:nc): values to be interpolated from; dimensions must match driving program
!     Outputs:
!       sout(1:nc): interpolated value @ z=zt (bottom value if ibelow=1)
!       ibelow: flag indicating if zt is below za(k1)
!
      subroutine vinter(mnv,nc,zt,k1,k2,k3,za,sint,sout,ibelow)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)
      integer, intent(in) :: mnv,nc,k1,k2,k3
      integer, intent(out) :: ibelow
      real(kind=dbl_kind1), intent(in) :: zt,za(mnv),sint(mnv,10)
      real(kind=dbl_kind1), intent(out) :: sout(10)

      if(k1>k2.or.nc>10) then
        write(11,*)'k1>k2 in vinter()'
        stop
      endif

      if(zt<za(k1)) then
        ibelow=1
        sout(1:nc)=sint(k1,1:nc)
      else !normal
        ibelow=0
        if(zt==za(k1)) then
          sout(1:nc)=sint(k1,1:nc)
        else if(zt>=za(k2)) then
          sout(1:nc)=sint(k2,1:nc)
        else
          kout=0 !flag
          if(k3<k1.or.k3>k2) then
            l1=k1; l2=k2-1
          else
            if(zt<za(k3)) then
              l1=k1; l2=k3-1
            else
              l1=k3; l2=k2-1
            endif
          endif
          do k=l1,l2
            if(zt>=za(k).and.zt<=za(k+1)) then
              kout=k
              exit
            endif
          enddo !k
          if(kout==0.or.za(kout+1)-za(kout)==0) then
            write(11,*)'Failed to find a level in vinter():',kout,zt,(za(k),k=k1,k2)
            stop
          endif
          zrat=(zt-za(kout))/(za(kout+1)-za(kout))
          sout(1:nc)=sint(kout,1:nc)*(1-zrat)+sint(kout+1,1:nc)*zrat
        endif
      endif

      return
      end

!     Flux limiter functions used in TVD schemes
      function flux_lim(ss,flimiter)
      implicit real*8(a-h,o-z)
      character(len=2) :: flimiter

      if(flimiter.eq.'SB') then !Superbee
        flux_lim=max(0.d0,min(1.d0,2*ss),min(2.d0,ss))
      else if(flimiter.eq.'MM') then !MINMOD
        flux_lim=max(0.d0,min(1.d0,ss))
      else if(flimiter.eq.'OS') then !OSHER
        flux_lim=max(0.d0,min(2.d0,ss))
      else if(flimiter.eq.'VL') then !Van Leer
        flux_lim=(ss+abs(ss))/(1+abs(ss))
      else
        write(11,*)'Unknown limiter:',flimiter
        stop
      endif

      return
      end

!     Compute local index of a side (0 if not a local side)
      function lindex_s(i,ie)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: i,ie
 
      l0=0 !local index
      do l=1,3
        if(js(ie,l)==i) then
          l0=l
          exit
        endif
      enddo !l
      lindex_s=l0

      return
      end

!     Compute local index of a node (0 if not a local node)
      function lindex(i,ie)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: i,ie

      l0=0 !local index
      do l=1,3
        if(nm(ie,l)==i) then
          l0=l
          exit
        endif
      enddo !l
      lindex=l0

      return
      end

      function covar(kr_co,decor,hh)
      implicit real*8(a-h,o-z)

      if(hh<0) then
        write(11,*)'Negative hh in covar:',hh
        stop
      endif

      if(kr_co==1) then
        covar=-hh
      else if(kr_co==2) then
        if(hh==0) then
          covar=0
        else
          covar=hh*hh*dlog(hh)
        endif
      else if(kr_co==3) then
        covar=hh*hh*hh
      else if(kr_co==4) then
        h2=hh*hh
        covar=-h2*h2*hh
      else if(kr_co==5) then
        covar=dexp(-hh*hh/decor/decor)
      else
        write(11,*)'Unknown covariance function option:',kr_co
        stop
      endif

      return
      end

! -----------------------------------------------------------------------
! *** Gauss elimination routine with full pivoting
! ***
!     Input:
!          a(np,np) (in/out): original and inverted square matrix;
!          n: actual rank of a and b;
!          np: dimension of a and b (must match the driving routine);
!          b(np,mp) (in/out): RHS or solution vectors;
!          m: actual # of columns on the RHS 
!          mp: used in dimensioning of RHS b (must match the driving routine);
      subroutine gaussj(a,n,np,b,m,mp)
      implicit real*8(a-h,o-z)
      parameter (nmax=50000)
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)

!     Check dimension
      if(nmax<np) then
        write(11,*)'Increase nmax in gaussj:',np,nmax
        stop
      endif

         do 11 j=1,n
            ipiv(j)=0
 11      continue
         do 22 i=1,n
            big=0.d0
            do 13 j=1,n
               if(ipiv(j).ne.1) then
                  do 12 k=1,n
                     if (ipiv(k).eq.0) then
                        if (abs(a(j,k)).ge.big) then
                           big=abs(a(j,k))
                           irow=j
                           icol=k
                        endif
                     else if (ipiv(k).gt.1) then

                        stop 'singular matrix in gaussj'
                     endif
 12               continue
               endif
 13         continue
            ipiv(icol)=ipiv(icol)+1
                                       
            if (irow.ne.icol) then                                          
               do 14 l=1,n                                                   
                  dum=a(irow,l)                                               
                  a(irow,l)=a(icol,l)                                         
                  a(icol,l)=dum                                               
 14            continue                                                      
               do 15 l=1,m                                                   
                  dum=b(irow,l)                                               
                  b(irow,l)=b(icol,l)                                         
                  b(icol,l)=dum                                               
 15            continue                                                      
            endif                                                           
            indxr(i)=irow                                                   
            indxc(i)=icol 
                    
            if (a(icol,icol).eq.0.d0) stop 'singular matrix in gaussj'
               
            pivinv=1.d0/a(icol,icol)                                          
            a(icol,icol)=1.d0                                               
            do 16 l=1,n                                                     
               a(icol,l)=a(icol,l)*pivinv                                    
 16         continue                                                        
            do 17 l=1,m                                                     
               b(icol,l)=b(icol,l)*pivinv                                    
 17         continue                                                        
            do 21 ll=1,n                                                    
               if(ll.ne.icol)then                                            
                  dum=a(ll,icol)                                              
                  a(ll,icol)=0.d0                                               
                  do 18 l=1,n                                                 
                     a(ll,l)=a(ll,l)-a(icol,l)*dum                             
 18               continue                                                    
                  do 19 l=1,m                                                 
                     b(ll,l)=b(ll,l)-b(icol,l)*dum                             
 19               continue                                                    
               endif                                                         
 21         continue                                                        
 22      continue                                                          
         do 24 l=n,1,-1                                                    
            if(indxr(l).ne.indxc(l))then                                    
               do 23 k=1,n                                                   
                  dum=a(k,indxr(l))                                           
                  a(k,indxr(l))=a(k,indxc(l))                                 
                  a(k,indxc(l))=dum                                           
 23            continue                                                      
            endif                                                           
 24      continue                                                          
      return                                                            
      end


!*************************************************************************************
!
!        Do upwind and TVD transport
!
!*************************************************************************************
      subroutine do_transport_tvd(imod,up_tvd,tvd_mid,flimiter,ntr)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: imod !=0: ST equations; 1: other tracers
      logical, intent(in) :: up_tvd !'T' if TVD is used (must be for all tracers)
      character(len=2), intent(in) :: tvd_mid,flimiter
      integer, intent(in) :: ntr !# of tracers
!      real(kind=dbl_kind), intent(in) :: bdy_frc(mnv,mne,mntr) !body force at prism center Q_{i,k}
!      real(kind=dbl_kind), intent(in) :: flx_sf(mne,mntr) !surface b.c. \kappa*dC/dz = flx_sf (at element center)
!      real(kind=dbl_kind), intent(in) :: flx_bt(mne,mntr) !bottom b.c.
!      real(kind=dbl_kind), intent(inout) :: tr_el(mnv_in,mne_in,ntr) !tracer converntration @ prism center

!     Working temporary arrays in this routine
      real(kind=dbl_kind) :: trel_tmp(mnv,mne,mntr) !tracer @ elements and half levels
      real(kind=dbl_kind) :: flx_adv(mns,mnv,2) ! original horizontal flux (1:  the local x-driection) and vertical flux (2: positive upward) 
      real(kind=dbl_kind) :: flx_mod(mns,mnv,mntr,2) !limited advective fluxes (1: horizontal; 2: vertical)
      real(kind=dbl_kind) :: up_rat(mns,mnv,mntr,2) !upwind ratios (1: horizontal; 2: vertical)

      real(kind=dbl_kind) :: swild(mnv),sne(mnv,3),area_e(mnv),psumtr(mntr),delta_tr(mntr),adv_tr(mntr), &
     &alow(mnv),bdia(mnv),cupp(mnv),rrhs(mnv,100),soln(mnv,100),gam(mnv) !"100" in rrhs & soln must match tridag()

!     List of temporary variables used in this option (may be altered outside this module)
!     hp_int: S,T @ elements and half levels;
!     bcc: save solar heat;
!     bubt: save surface flux (ns >= ne);
!     swild4(:,:,1:2): original horizontal flux (1:  the local x-driection) and vertical flux (2: positive upward)
!     ptbt(:,:,1:4): upwind ratios for T,S (order follows sdbt)
!     sdbt(:,:,1:4): limited horizontal flux (1 & 3 (for T,S); in the local x-driection) and vertical flux (2 & 4: positive upward)

!     Check
      if(ntr>100) then
        write(11,*)'Too many tracers:',ntr
        stop
      endif

!     For rewetted elements, tr_el takes the value from last wet step

!     Compute (pre-limiting) fluxes at all faces 
      flx_adv=-1.e34 !flags
      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet element with 3 wet nodes
!       Compute upward normals and areas @ all levels
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        isd1=js(i,1)
        isd2=js(i,2)
        isd3=js(i,3)
        if(kbe(i)==0) then
          write(11,*)'Impossible 95 (2)'
          stop
        endif
        do l=kbe(i),nvrt
          xcon=(y(n2)-y(n1))*(z(l,n3)-z(l,n1))-(y(n3)-y(n1))*(z(l,n2)-z(l,n1))
          ycon=(x(n3)-x(n1))*(z(l,n2)-z(l,n1))-(x(n2)-x(n1))*(z(l,n3)-z(l,n1))
          zcon=area(i)*2
          area_e(l)=dsqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0) then
            write(11,*)'Zero area (2):',i,l
            stop
          endif
          sne(l,1)=xcon/area_e(l)/2
          sne(l,2)=ycon/area_e(l)/2
          sne(l,3)=zcon/area_e(l)/2 !>0
        enddo !l

!       Compute vertical fluxes first
        do k=kbe(i),nvrt
          if(k==kbe(i)) then !bottom normal vel. is we(kbe(i),i)
            dot1=we(kbe(i),i)
          else
            dot1=(su2(k,isd1)+su2(k,isd2)+su2(k,isd3))/3*sne(k,1)+ & !upward normal vel.
     &           (sv2(k,isd1)+sv2(k,isd2)+sv2(k,isd3))/3*sne(k,2)+we(k,i)*sne(k,3)
          endif
          flx_adv(i,k,2)=dot1*area_e(k) !vertical flux (positive upward)

!               Debug
!                if(it==46.and.i==58422) write(99,*)k,we(k,i),dot1,area_e(k),flx_adv(i,k,2)

          if(k/=kbe(i)) swild(k)=flx_adv(i,k,2)-flx_adv(i,k-1,2) !local volume conservation
!         Limit flux later
        enddo !k=kbe(i),nvrt

!       Horizontal fluxes
        do j=1,3 !side
          jsj=js(i,j)
          iel=ic3(i,j)

          do k=kbe(i)+1,nvrt
            if(flx_adv(jsj,k,1)>-1.e33) cycle !already computed
            if(iel==0.and.isbs(jsj)==0) then !land
              flx_adv(jsj,k,1)=0 
              cycle
            endif
            
!           Check near bottom vel.
            tmp=dabs(su2(k,jsj))+dabs(sv2(k,jsj))+dabs(su2(k-1,jsj))+dabs(sv2(k-1,jsj))
            if(iel/=0.and.k<kbe(iel).and.tmp/=0) then
              write(11,*)'Non-zero hvel below element bottom:',i,j,iel,k,tmp
              stop
            endif

            vnor1=su2(k,jsj)*snx(jsj)+sv2(k,jsj)*sny(jsj)
            vnor2=su2(k-1,jsj)*snx(jsj)+sv2(k-1,jsj)*sny(jsj)
            if(k<kbs(jsj)+1) then
              write(11,*)'Side/element indices conflict:',i,k
              stop
            endif
            flx_adv(jsj,k,1)=(zs(k,jsj)-zs(k-1,jsj))*distj(jsj)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)
            tmp=ssign(i,j)*flx_adv(jsj,k,1) !local outward flux
            swild(k)=swild(k)+tmp !volume conservation metric

!               Debug
!                if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flx_adv(jsj,k,1)

          enddo !k=kbe(i)+1,nvrt
        enddo !j=1,3
      enddo !i=1,ne

!     Initialize limited flux for upwind
      do i=1,ntr
        flx_mod(1:mns,1:mnv,i,1:2)=flx_adv(1:mns,1:mnv,1:2)
      enddo !i
   
      it_sub=0
      time_r=dt !time remaining
      loop11: do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      it_sub=it_sub+1

!     Compute flux limiters and modify fluxes
      if(up_tvd) then !neither can be upwind any more
        up_rat=-1.e34 !flags
!       Vertical limiters
        ntot_v=0 !total # of vertical faces that have large limiters (for first tracer)
        do i=1,ne
          if(idry_e(i)==1) cycle

          up_rat(i,1:mnv,1:ntr,2)=-1 !initialize upwind ratio
          do k=kbe(i)+1,nvrt-1 !bottom and surface flux unchanged at -1
            if(flx_adv(i,k,2)<-1.e33) then
              write(11,*)'Left out vertical flux (3):',i,k
              stop
            else if(flx_adv(i,k,2)>0) then
              kup=k !upwind prism
              kdo=k+1 !downwind prism
            else
              kup=k+1 
              kdo=k
            endif

            psum=0 !sum of original fluxes
            psumtr(1:ntr)=0 !sum of products
            if(flx_adv(i,kup,2)<-1.e33.or.flx_adv(i,kup-1,2)<-1.e33) then
              write(11,*)'Left out vertical flux (4):',i,kup
              stop
            endif
            if(flx_adv(i,kup,2)<0.and.kup/=nvrt) then
              psum=psum+abs(flx_adv(i,kup,2))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(i,kup,2))*(tr_el(kup+1,i,1:ntr)-tr_el(kup,i,1:ntr))
            endif
            if(flx_adv(i,kup-1,2)>0.and.kup/=kbe(i)+1) then
              psum=psum+abs(flx_adv(i,kup-1,2))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(i,kup-1,2))*(tr_el(kup-1,i,1:ntr)-tr_el(kup,i,1:ntr))
            endif
            do j=1,3
              jsj=js(i,j)
              ie=ic3(i,j)
              if(flx_adv(jsj,kup,1)<-1.e33) then
                write(11,*)'Left out horizontal flux (5):',jsj,kup
                stop
              endif
              if(ie/=0.and.idry_e(ie)==0.and.ssign(i,j)*flx_adv(jsj,kup,1)<0) then
                psum=psum+abs(flx_adv(jsj,kup,1))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(jsj,kup,1))*(tr_el(kup,ie,1:ntr)-tr_el(kup,i,1:ntr))
              endif
            enddo !j

            if(tvd_mid.eq.'AA') then !my formulation
              do j=1,ntr
                tmp=(tr_el(kup,i,j)-tr_el(kdo,i,j))*abs(flx_adv(i,k,2))
                if(abs(tmp)>1.e-20) up_rat(i,k,j,2)=psumtr(j)/tmp
              enddo !j
            else if(tvd_mid.eq.'CC') then !Casulli's
              do j=1,ntr
                tmp=(tr_el(kup,i,j)-tr_el(kdo,i,j))*psum
                if(abs(tmp)>1.e-20) up_rat(i,k,j,2)=psumtr(j)/tmp
              enddo !j
            else
              write(11,*)'Unknown tvd_mid:',tvd_mid
              stop
            endif

            if(flux_lim(up_rat(i,k,1,2),flimiter)>0.1) ntot_v=ntot_v+1
            
          enddo !k=kbe(i)+1,nvrt-1
        enddo !i=1,ne

!       Horizontal limiters
        ntot_h=0 !total # of horizontal faces that have large limiters (for 1st tracer)
        do i=1,ns
          if(idry_s(i)==1) cycle

!         At least one element is wet
          up_rat(i,1:mnv,1:ntr,1)=-1 !initialize (for below bottom etc)
          if(is(i,2)==0.or.idry_e(is(i,2))==1.or.idry_e(is(i,1))==1) cycle

!         Not bnd face; 2 elements are wet
          kb1=min(kbe(is(i,1)),kbe(is(i,2)))
          kb=max(kbe(is(i,1)),kbe(is(i,2)))
          do k=kb1+1,kb-1
            if(flx_adv(i,k,1)/=0) then
              write(11,*)'Pls zero out the excess layers:',flx_adv(i,k,1),i,is(i,1),is(i,2),k,kb1,kb
              stop
            endif
          enddo !k
 
!         Leave k=kb unchanged
          do k=kb+1,nvrt !prisms
            if(flx_adv(i,k,1)<-1.e33) then
              write(11,*)'Left out horizontal flux (3):',i,k
              stop
            else if(flx_adv(i,k,1)>0) then
              iup=is(i,1); ido=is(i,2)
            else
              iup=is(i,2); ido=is(i,1)
            endif

            psum=0
            psumtr(1:ntr)=0
            if(flx_adv(iup,k,2)<-1.e33.or.flx_adv(iup,k-1,2)<-1.e33) then
              write(11,*)'Left out vertical flux (6):',iup,k
              stop
            endif
            if(flx_adv(iup,k,2)<0.and.k/=nvrt) then
              psum=psum+abs(flx_adv(iup,k,2))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(iup,k,2))*(tr_el(k+1,iup,1:ntr)-tr_el(k,iup,1:ntr))
            endif
            if(flx_adv(iup,k-1,2)>0.and.k>kbe(iup)+1) then
              psum=psum+abs(flx_adv(iup,k-1,2))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(iup,k-1,2))*(tr_el(k-1,iup,1:ntr)-tr_el(k,iup,1:ntr))
            endif
            do j=1,3
              jsj=js(iup,j)
              ie=ic3(iup,j)
              if(flx_adv(jsj,k,1)<-1.e33) then
                write(11,*)'Left out horizontal flux (6):',jsj,k
                stop
              endif
              if(ie/=0.and.idry_e(ie)==0.and.ssign(iup,j)*flx_adv(jsj,k,1)<0) then
                psum=psum+abs(flx_adv(jsj,k,1))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(jsj,k,1))*(tr_el(k,ie,1:ntr)-tr_el(k,iup,1:ntr))
              endif
            enddo !j
     
            if(tvd_mid.eq.'AA') then
              do j=1,ntr
                tmp=(tr_el(k,iup,j)-tr_el(k,ido,j))*abs(flx_adv(i,k,1))
                if(abs(tmp)>1.e-20) up_rat(i,k,j,1)=psumtr(j)/tmp
              enddo !j
            else !model CC
              do j=1,ntr
                tmp=(tr_el(k,iup,j)-tr_el(k,ido,j))*psum
                if(abs(tmp)>1.e-20) up_rat(i,k,j,1)=psumtr(j)/tmp
              enddo !j
            endif

            if(flux_lim(up_rat(i,k,1,1),flimiter)>0.1) ntot_h=ntot_h+1
          enddo !k=kb+1,nvrt
        enddo !i=1,ns

!       Debug
!        if(it==1.and.it_sub==1) then
!          do i=1,ne
!            do j=1,3
!              jsj=js(i,j)
!              write(99,*)is(jsj,1),is(jsj,2),up_rat(jsj,nvrt,1)
!            enddo !j
!          enddo !i
!          stop
!        endif

!       Modifed fluxes flx_mod (their signs do not change) 
!       Vertical fluxes
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1 !leave out the bnd
!           Compute \delta_i
            if(flx_adv(i,k,2)>0) then
              kup=k !upwind prism
            else
              kup=k+1
            endif

            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces
              if(flx_adv(i,kup-l,2)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat(i,kup-l,j,2)
                  if(rat<-1.e33) then
                    write(11,*)'Left out (1):',i,kup-l,rat,it_sub,j
                    stop
                  else if(rat/=0) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(11,*)'Flux limiting failed (1):',tmp,rat,flx_adv(i,kup-l,2),l,kup
                      stop
                    endif 
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=js(i,j)
              ie=ic3(i,j)
              if(ssign(i,j)*flx_adv(jsj,kup,1)>0) then
                do jj=1,ntr
                  rat=up_rat(jsj,kup,jj,1)
                  if(rat<-1.e33) then
                    write(11,*)'Left out (3):',i,j,kup,rat,jj
                    stop
                  else if(rat/=0) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(11,*)'Flux limiting failed (3):',tmp,rat,jj
                      stop
                    endif 
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif
            enddo !j=1,3

            do j=1,ntr
              flx_mod(i,k,j,2)=flx_adv(i,k,2)*(1-flux_lim(up_rat(i,k,j,2),flimiter)/2+delta_tr(j))
            enddo !j
          enddo !k=kbe(i)+1,nvrt-1  
        enddo !i=1,ne

!       Horizontal fluxes
        do i=1,ns
          if(idry_s(i)==1.or.is(i,2)==0.or.idry_e(is(i,1))==1.or.idry_e(is(i,2))==1) cycle

!         Both elements are wet
          kb=max(kbe(is(i,1)),kbe(is(i,2)))
          do k=kb+1,nvrt
            if(flx_adv(i,k,1)>0) then
              iup=is(i,1)
            else
              iup=is(i,2)
            endif
 
            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces
              if(flx_adv(iup,k-l,2)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat(iup,k-l,j,2)
                  if(rat<-1.e33) then
                    write(11,*)'Left out (5):',iup,k-l,rat,j
                    stop
                  else if(rat/=0) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(11,*)'Flux limiting failed (5):',tmp,rat,j
                      stop
                    endif
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=js(iup,j)
              ie=ic3(iup,j)
              if(ssign(iup,j)*flx_adv(jsj,k,1)>0) then !outflow
                do jj=1,ntr
                  rat=up_rat(jsj,k,jj,1)
                  if(rat<-1.e33) then
                    write(11,*)'Left out (7):',iup,ie,k,rat,jj
                    stop
                  else if(rat/=0) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(11,*)'Flux limiting failed (7):',tmp,rat,jj
                      stop
                    endif
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif !outflow
            enddo !j=1,3

            do j=1,ntr
              flx_mod(i,k,j,1)=flx_adv(i,k,1)*(1-flux_lim(up_rat(i,k,j,1),flimiter)/2+delta_tr(j)) 
            enddo !j
          enddo !k=kb+1,nvrt
        enddo !i=1,ns

      endif !flux limiter

!     Compute sub time step
!     Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(i,j) is dry)
!     Causion: \hat{S}^- conditions must be consistent later in the advective flux part!!!!!!
!     Implicit vertical flux for upwind; explicit for TVD

      if(up_tvd.or.it_sub==1) then !for upwind, only compute dtb for the first step
        dtb=time_r
        dtb_alt=time_r !alternative for TVD (more restrictive)
        ie01=0 !element # where the exteme is attained
        lev01=0 !level #
        in_st=0 !tracer #
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt !prism
            qj=0 !sum of original fluxes for all inflow bnds
            psumtr(1:ntr)=0 !sum of modified fluxes for all inflow bnds
            nplus=0 !# of outflow bnds
   
            if(up_tvd) then !neither can be upwind any more
              if(k/=nvrt.and.flx_mod(i,k,1,2)<0) then !flx_mod and flx_adv same sign
                qj=qj+abs(flx_adv(i,k,2))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_mod(i,k,1:ntr,2))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flx_adv(i,k,2),flx_mod(i,k,2)
              endif
              if(k-1/=kbe(i).and.flx_mod(i,k-1,1,2)>0) then
                qj=qj+abs(flx_adv(i,k-1,2))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_mod(i,k-1,1:ntr,2))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flx_adv(i,k-1,2),flx_mod(i,k-1,2)
              endif
            endif !flux limiter

            if(k/=nvrt.and.flx_adv(i,k,2)>0) nplus=nplus+1
            if(k-1/=kbe(i).and.flx_adv(i,k-1,2)<0) nplus=nplus+1 

            do j=1,3
              jsj=js(i,j)
              ie=ic3(i,j)
              do jj=1,ntr
                if(flx_mod(jsj,k,jj,1)<-1.e33) then
                  write(11,*)'Left out horizontal flux (10):',i,k,j,jj
                  stop
                endif
              enddo !jj=1,ntr

              if(ie/=0.and.idry_e(ie)==0) then
                if(ssign(i,j)*flx_mod(jsj,k,1,1)<0) then !flx_mod(:) same sign as flx_adv
                  qj=qj+abs(flx_adv(jsj,k,1))
                  psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_mod(jsj,k,1:ntr,1))
!                 Debug
!                   if(it==46.and.it_sub==1.and.i==58422) write(99,*)j,k,flx_adv(jsj,k,1),flx_mod(jsj,k,1)
                endif
                if(ssign(i,j)*flx_adv(jsj,k,1)>0) nplus=nplus+1
              endif
            enddo !j

            vj=area(i)*(ze(k,i)-ze(k-1,i))

!               Debug
!                if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,nplus,vj

            do jj=1,ntr
              if(psumtr(jj)/=0) then
                tmp=vj/psumtr(jj)*(1-1.e-6) !safety factor included
                if(tmp<dtb) then
                  dtb=tmp 
                  ie01=i; lev01=k; in_st=jj
                endif
              endif
            enddo !jj

            if(qj/=0) dtb_alt=min(dtb_alt,vj/(1+nplus)/qj*(1-1.e-10)) !safety factor included
          enddo !k=kbe(i)+1,nvrt
        enddo !i=1,ne

        if(dtb<=0.or.dtb>time_r) then
          write(11,*)'Illegal sub step:',dtb,time_r
          stop
        endif

!       Output time step
        if(up_tvd) write(18,*)it,it_sub,dtb,dtb_alt,ie01,lev01,in_st
      endif !compute dtb

      dtb=min(dtb,time_r) !for upwind
      time_r=time_r-dtb

!     Store last step's S,T
      trel_tmp=tr_el

      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with 3 wet nodes
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)

!       Check if having a dry neighbor (interface)
!       Interface element will not be subject to ELAD because the vel. there may be altered by levels()?
!        iedge=0 !flag
!        do j=1,3
!          if(ic3(i,j)/=0.and.idry_e(ic3(i,j))==1) iedge=1
!        enddo !j 

!       Matrix
!        tmin=tempmax+1; tmax=tempmin-1 !extrema for ELAD
!        smin=saltmax+1; smax=saltmin-1
        ndim=nvrt-kbe(i)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i) 
          alow(kin)=0
          cupp(kin)=0
          bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
          if(bigv<=0) then
            write(11,*)'Negative volume: ',bigv,i,k
            stop
          endif
          bdia(kin)=1
          if(k<nvrt) then
            av_df=(dfh(n1,k)+dfh(n2,k)+dfh(n3,k))/3
            av_dz=(ze(k+1,i)-ze(k-1,i))/2
            if(av_dz<=0) then
              write(11,*)'Impossible 111'
              stop
            endif
            tmp=area(i)*dtb*av_df/av_dz/bigv
            cupp(kin)=cupp(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif

          if(k>kbe(i)+1) then
            av_df=(dfh(n1,k-1)+dfh(n2,k-1)+dfh(n3,k-1))/3
            av_dz=(ze(k,i)-ze(k-2,i))/2
            if(av_dz<=0) then
              write(11,*)'Impossible 112'
              stop
            endif
            tmp=area(i)*dtb*av_df/av_dz/bigv
            alow(kin)=alow(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif

!         b.c. to be imposed at the end
!         Advective flux
!         Strike out \hat{S}^- (see above)
          psumtr(1:ntr)=0 !sum of modified fluxes at all inflow bnds 
          delta_tr(1:ntr)=0 !sum of tracer fluxes at all inflow bnds
          adv_tr(1:ntr)=trel_tmp(k,i,1:ntr) !alternative mass conservative form for the advection part
          if(ntr>1.and.flx_mod(i,k,1,2)*flx_mod(i,k,2,2)<0) then
            write(11,*)'Left out vertical flux (0):',i,k,flx_mod(i,k,1:2,2)
            stop
          endif
          do jj=1,ntr
            if(flx_mod(i,k,jj,2)<-1.e33) then
              write(11,*)'Left out vertical flux:',i,k,flx_mod(i,k,jj,2),jj
              stop
            endif
          enddo !jj

          if(k/=nvrt.and.flx_mod(i,k,1,2)<0) then !all flx_mod(:) same sign
            if(up_tvd) then !neither can be upwind any more
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(i,k,jj,2))
                delta_tr(jj)=delta_tr(jj)+abs(flx_mod(i,k,jj,2))*trel_tmp(k+1,i,jj)
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(i,k,2))*(trel_tmp(k+1,i,jj)-trel_tmp(k,i,jj))
              enddo !jj
            else !upwind
              tmp=abs(flx_mod(i,k,1,2))*dtb/bigv !flx_mod(:) all same for upwind
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif
          if(k-1/=kbe(i).and.flx_mod(i,k-1,1,2)>0) then
            if(up_tvd) then !neither can be upwind any more
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(i,k-1,jj,2))
                delta_tr(jj)=delta_tr(jj)+abs(flx_mod(i,k-1,jj,2))*trel_tmp(k-1,i,jj)
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(i,k-1,2))*(trel_tmp(k-1,i,jj)-trel_tmp(k,i,jj))
              enddo !jj
            else !upwind
              tmp=abs(flx_mod(i,k-1,1,2))*dtb/bigv
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif

!         Additional terms in adv_tr
          if(up_tvd) then
            if(k/=nvrt) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(i,k,2))*(trel_tmp(k,i,jj)-trel_tmp(k+1,i,jj))* &
     &flux_lim(up_rat(i,k,jj,2),flimiter)/2
              enddo !jj
            endif
            if(k-1/=kbe(i)) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(i,k-1,2))*(trel_tmp(k,i,jj)-trel_tmp(k-1,i,jj))* &
     &flux_lim(up_rat(i,k-1,jj,2),flimiter)/2
              enddo !jj
            endif
          endif !TVD

          do j=1,3
            jsj=js(i,j)
            iel=ic3(i,j)
            if(iel==0.or.idry_e(iel)==1) cycle 

            if(ntr>1.and.flx_mod(jsj,k,1,1)*flx_mod(jsj,k,2,1)<0) then
              write(11,*)'Left out horizontal flux (0):',i,j,k,flx_mod(jsj,k,1:2,1)
              stop
            endif
            do jj=1,ntr
              if(flx_mod(jsj,k,jj,1)<-1.e33) then
                write(11,*)'Left out horizontal flux:',i,j,k,flx_mod(jsj,k,jj,1),jj
                stop
              endif
            enddo !jj

            if(ssign(i,j)*flx_mod(jsj,k,1,1)<0) then
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(jsj,k,jj,1))
                delta_tr(jj)=delta_tr(jj)+abs(flx_mod(jsj,k,jj,1))*trel_tmp(k,iel,jj)
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(jsj,k,1))*(trel_tmp(k,iel,jj)-trel_tmp(k,i,jj))
              enddo !jj
            endif

            if(up_tvd) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(jsj,k,1))*(trel_tmp(k,i,jj)-trel_tmp(k,iel,jj))* &
     &flux_lim(up_rat(jsj,k,jj,1),flimiter)/2
              enddo !jj
            endif
          enddo !j=1,3

!         Check Courant number
          do jj=1,ntr
            if(1-dtb/bigv*psumtr(jj)<0) then
              write(11,*)'Courant # condition violated:',i,k,1-dtb/bigv**psumtr(jj),jj
              stop
            endif
          enddo !jj

          rrhs(kin,1:ntr)=adv_tr(1:ntr) !(1-dtb/bigv*qj1)*trel_tmp(k,i,1)+dtb/bigv*qjt

!         Check consistency between 2 formulations in TVD
!            if(up_tvd) then 
!              if(abs(adv_t-rrhs(kin,1))>1.e-4.or.abs(adv_s-rrhs(kin,2))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(kin,1),adv_s,rrhs(kin,2)
!                stop
!              endif
!            endif !TVD

!         Body source
          rrhs(kin,1:ntr)=rrhs(kin,1:ntr)+dtb*bdy_frc(k,i,1:ntr)

!         b.c.
          if(k==nvrt) rrhs(kin,1:ntr)=rrhs(kin,1:ntr)+area(i)*dtb*flx_sf(i,1:ntr)/bigv
          if(k==kbe(i)+1) rrhs(kin,1:ntr)=rrhs(kin,1:ntr)-area(i)*dtb*flx_bt(i,1:ntr)/bigv
        enddo !k=kbe(i)+1,nvrt

!       if(tmin>tmax.or.tmin<tempmin.or.smin>smax.or.smin<saltmin) then
!         write(11,*)'Illegal min/max:',tmin,tmax,smin,smax,i
!         stop
!       endif

        call tridag(mnv,ndim,ntr,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i)
          tr_el(k,i,1:ntr)=soln(kin,1:ntr)
       
          if(imod==0) then !ST
            if(ihconsv/=0) then
              tr_el(k,i,1)=dmax1(tempmin,dmin1(tempmax,soln(kin,1)))
            endif
            if(isconsv/=0) then
              tr_el(k,i,2)=dmax1(saltmin,dmin1(saltmax,soln(kin,2)))
            endif
          endif !imod
        enddo !k

!       Extend
        do k=1,kbe(i)
          tr_el(k,i,1:ntr)=tr_el(kbe(i)+1,i,1:ntr)
        enddo !k
      enddo !i=1,ne

      if(time_r<1.e-8) exit loop11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       end do loop11

!     Output # of divisions etc.
      if(up_tvd) then 
        write(17,*)'Total # of vertical and S faces limited = ',ntot_h,ntot_v
      endif
      write(17,*)it,it_sub
      
      return
      end
