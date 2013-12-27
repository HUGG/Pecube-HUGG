module definitions

  type faulttype

! A fault vertical cross-sectional (2D) geometry is defined by a set of points (x,y)
! Its 3D location is defined by two points in the horizontal plane (x1,y1),(x2,y2)
! Those points limit its lateral extent
! In the (x,y) coordinate system (in which the faul vertical geometry is defined)
! x lies to the right of the (x1,y1),(x2,y2) line and z is vertical upwards
! and zero at the surface

! n is the number of points used to define the fault
! x(n) and y(n) are the coordinates of the points
! samef(n), samex(n) and samey(n) can be used to set the x or y coordinates of a
! given fault node to be equal to those of another fault node (on fault samef)
! (x1,y1) and (x2,y2) are the coordinates of the trace of the fault at the surface
! nstep is the number of time intevals used to define the motion story of the fault
! per interval: timestart,timeend a,d velo
! timestart and timeend are the start and end time of the interval
! velo is the velocity across the fault (along the plane of the fault) for the interval
! conv_rate is used to define a convergence rate across a pair of thrust faults
! with the same geometry
! conv_factor is used to partition conv_rate between over- and under-thrusting

! By convention as one goes along the fault from point 1 to point n
! it is the block to the right that moves

! The sign of velo determines whether the sign of the x-component of velocity

! from that information we compute xs(n-1),ys(n-1), the direction
! of each of the (n-1) segments connecting the n points

! as well as xn,yn the normal to the (x1,y1)-(x2,y2) line in the z=0 plane

    integer n,nstep,faultmodel
    integer,dimension(:),pointer::samef,samex,samey
    double precision,dimension(:),pointer::x,y
    double precision x1,y1,x2,y2
    double precision,dimension(:),pointer::timestart,timeend,velo
    double precision,dimension(:),pointer::xs,ys,alpha,kbm
    double precision,dimension(:),pointer::conv_rate,conv_factor
    double precision xn,yn
    ! Duplex parameters
    double precision dstart,dend,dinner,douter,dvelo

  end type faulttype

!=====[EDGE]====================================================================
!this type is to store edges in a trianglulation
! it is used to update (in a generalized Delaunay sense)
! the triangulation of the 3D points on the surfaces
! for each edge:
! n1, n2 are the node numbers defining the edge
! t1, t2 are the triangle numbers on either side of the edge
! going from n1 to n2, t1 is to the left and t2 is to the right
! m1, m2 are the node numbers of the two other nodes making t1 and t2

  type edge
    integer n1,n2,m1,m2,t1,t2
  end type edge

end module definitions

!---------------------------

      program Pecube

      use definitions

      type (faulttype) fault(1024)
      real*4 param(1024),range(2,1024),misfit
      integer nd,nfault
      integer nproc,iproc,ierr
! MOD Jean 31/10/08
      character iq*1
! end MOD

      iproc=0
      nproc=1

      if (iproc.eq.0) then
      nd=0
      if (iproc.eq.0) then
      open (77,file='input/fault_parameters.txt',status='old')
1     read (77,'(a1)') iq
      if (iq.eq.'$' .or. iq.eq.' ') goto 1
      backspace (77)
      read (77,*) nfault
      if (nfault.eq.0) nfault=1
      close (77)
      endif
      call create_pecube_in (fault,nfault,nd,range,param)
      close (7)
      endif

      if (nd.eq.0) then
      call forward (nd,param,misfit)
      else
      if (iproc.ne.0) stop 'cannot run the parallel version of Pecube in NA mode'
      call NA (nd,range)
      endif

      end

!---------------------------

      subroutine user_init (nd,range,scales)

      real*4 range(2,*)
      real*4 scales(*)
      integer nd

      scales(1:nd)=-1.

      return
      end
      
!---------------------------

      subroutine writemodels (nd,ntot,models,misfit,ns1,ns2,itmax, &
                              nh_max,nh,header)

      real*4 models (nd,*)
      real*4 misfit (ntot)
      character*(*) header

      open (87,file='NA/NA_results.txt',status='unknown')
        do i=1,ntot
        write (87,'(100e15.7)') misfit(i),(models(k,i),k=1,nd)
        enddo
      close (87)

      return
      end

!---------------------------

      subroutine forward (nd,param,misfit)

      use definitions

!      implicit real*8 (a-h,o-z)
      implicit none

      type (faulttype),dimension(:),allocatable  :: fault,faultp

      double precision,dimension(:),allocatable :: x,y,z,xp,yp,zp
      double precision,dimension(:),allocatable :: u,up,temp,tempp !by VKP
      double precision,dimension(:),allocatable :: t,tp,f
      integer,dimension(:,:),allocatable :: icon
      double precision,dimension(:,:,:),allocatable :: ael
      double precision,dimension(:,:),allocatable :: bel
      integer, dimension(:),allocatable :: kfix,ielsurf,ieo
      double precision,dimension(:),allocatable :: xsurf,ysurf,zsurf,zsurfp
      double precision,dimension(:),allocatable :: usurf,usurfp,tempsurf,tempsurfp !by VKP
      double precision,dimension(:),allocatable :: topoa,topob,rsurf,rebound
      integer,dimension(:,:),allocatable :: iconsurf,neighbour
      double precision,dimension(:,:),allocatable :: xdepth,ydepth,zdepth
      double precision,dimension(:,:),allocatable::ftdist,ftdisto
      double precision,dimension(:),allocatable :: xdeptho,ydeptho,zdeptho,zsurfo,zsurfpo
      double precision,dimension(:,:),allocatable :: tprev
      double precision,dimension(:),allocatable :: tprevo
      double precision, dimension(:,:),allocatable:: xexhumation,yexhumation,zexhumation
      double precision, dimension(:),allocatable:: xexhumationo,yexhumationo,zexhumationo
      double precision,dimension(:),allocatable::we1,we2,we3,we4
      integer,dimension(:),allocatable :: proc,eproc,ice
      real*4,dimension(:),allocatable :: t2,t2p,t3,t4,t5
      real*4,dimension(:),allocatable::vxx,vyy,vzz,exhumation

      character run*5,cstep*3,iq*1

      real*4 times(1000)
      real*4,dimension(:),allocatable :: age1,age2,age3,age4,age5,age6,age7,age8
      real*4,dimension(:),allocatable :: raman
      real*4 t1

      integer nproc,iproc,nfault,npe,mpe,nsurf,nz,nelemsurf,nstep,ilog,iterative,interpol,ierr
      integer isoflag,nx,ny,nxiso,nyiso,niteradvec,istep,kstep
      integer i1,i2,i3,i4,iout,istatic,ntime,nnode,nelem,ie,iesurf,i,j,k,kk,ido,nelemloc,je
      integer irec,jrec,in,itime,ij,iu,ic,ieobs,inp,iobs,krec,ktime,niter,nobs !VKP
      integer ftlflag,mftflag,fltflag
      double precision x1,x2,x3,x4,x5,dt,time,alpha,zsurfmax,fact,tsurf,zh,ftime,ftimep,tprevious,tfake
      double precision eps,zl,diffusivity,heatproduction,zsurf1,zsurf2
      double precision tempsurf1,tempsurf2 !by VKP
      double precision tmax,tmsl,tlapse,tau,rhoc,rhom,xstep,ystep,young,poisson,thickness
      double precision xlonmin,xlonmax,xlatmin,xlatmax,xxx,yyy,zzz,vx,vy,vz,dxxx,dyyy,dzzz
      double precision xmin,xmax,ymin,ymax,timesurf,timesurfp,tfinal,dtimesurf
      double precision xlonobs,xlatobs,tnow,vzobs
      double precision ageheobs,dageheobs
      double precision ageftobs,dageftobs
      double precision ageheZobs,dageheZobs
      double precision ageftZobs,dageftZobs
      double precision agearKobs,dagearKobs
      double precision agearBobs,dagearBobs
      double precision agearMobs,dagearMobs
      double precision agearHobs,dagearHobs
      double precision tlobs(17),tlZobs(17),length(17)
      double precision Tramanobs,dTramanobs
      real*4,dimension(:),allocatable::ageft,agehe,ageheZ,ageftZ,agearK,agearB
      real*4,dimension(:),allocatable::agearM,agearH,hei,Traman
      double precision wobs1,wobs2,wobs3,wobs4,x1f,x2f,y1f,y2f
      double precision def,dif,heightobs,friction
      double precision,dimension(:),allocatable::aheoreg,hheoreg,ahepreg
      double precision,dimension(:),allocatable::aftoreg,hftoreg,aftpreg
      double precision,dimension(:),allocatable::aheZoreg,hheZoreg,aheZpreg
      double precision,dimension(:),allocatable::aftZoreg,hftZoreg,aftZpreg
      double precision,dimension(:),allocatable::aarKoreg,harKoreg,aarKpreg
      double precision,dimension(:),allocatable::aarBoreg,harBoreg,aarBpreg
      double precision,dimension(:),allocatable::aarMoreg,harMoreg,aarMpreg
      double precision,dimension(:),allocatable::aarHoreg,harHoreg,aarHpreg
      double precision,dimension(:),allocatable::aramanoreg,hramanoreg,aramanpreg
      double precision,dimension(:),allocatable::hpreg
      integer nhe,nft,nheZ,nftZ,narK,narB,narM,narH,nraman,ageflag(10)
      integer advect_topo
      double precision aheom,hheom,heoi,heos
      double precision aftom,hftom,ftoi,ftos
      double precision aheZom,hheZom,heZoi,heZos
      double precision aftZom,hftZom,ftZoi,ftZos
      double precision aarKom,harKom,arKoi,arKos
      double precision aarBom,harBom,arBoi,arBos
      double precision aarMom,harMom,arMoi,arMos
      double precision aarHom,harHom,arHoi,arHos
      double precision ahepm,hhepm,hepi,heps
      double precision aftpm,hftpm,ftpi,ftps
      double precision aheZpm,hheZpm,heZpi,heZps
      double precision aftZpm,hftZpm,ftZpi,ftZps
      double precision aarKpm,harKpm,arKpi,arKps
      double precision aarBpm,harBpm,arBpi,arBps
      double precision aarMpm,harMpm,arMpi,arMps
      double precision aarHpm,harHpm,arHpi,arHps
      double precision vtopo,vxtopo,vytopo,vxtopo0

      logical interpol_success,misfit_slope,vivi,l1,calc_surface_ages,PTt_output

      integer nd,nmisfit,ndmisfit,misfit_type
      real*4 param(nd)
      real*4 misfit
      real*4 range(2,1024)

      double precision :: cur_depth,xy_mean,zmax,zmin,zmin_orig,scale_fact
      integer :: num_left,nzin,plane_store1,plane_store2
      logical :: depth_flag

      if (nd.ne.0) print*,nd,param(1:nd)

      call cpu_time (times(1))
      nproc=1
      iproc=0

      eps=tiny(eps)

        do i=1,17
        length(i)=i
        enddo

! Pecube is a Finite Element solver of the 3D, transient heat transfer
! equation that allows for conduction, vertical advection and production of
! heat. Pecube also allows for the surface geometry to vary with time.

! This version is an improvement (we hope) of the original version
! available by download from the author's website in that it reads an
! input (topography) file as well as a set of cooling ages (He and FT in
! apatite of known location and calculates synthetic ages at the same
! locations for comparison.

! A flexural isostatic model has also been included to
! calculate the effect of isostasy in amplifying the exhumation caused
! by erosional unloading.

! Two dimensional faults have also been included; they can be orientated
! in any direction; more than one can be active at any given time and
! they advect each other (this only works for simple scenario)

! To understand how this version works, read the information in the file
! README

! comment this line out if you want Pecube to read your own input file
! (named Pecube.in)

      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Reading input'

      if (iproc.eq.0) then
      open (77,file='input/fault_parameters.txt',status='old')
1     read (77,'(a1)') iq
      if (iq.eq.'$' .or. iq.eq.' ') goto 1
      backspace (77)
      read (77,*) nfault
      if (nfault.eq.0) nfault=1
      close (77)
      endif

      if (nfault.gt.0) allocate (fault(nfault),faultp(nfault))
      call create_pecube_in (fault,nfault,nd,range,param)
      call cpu_time (times(2))

! opens input files

      rewind (7)

! read in general information

! first line:
! run: 5 character string that will determine the name of the folder where the input file
!      will be copied to and where the output files (Pecube.ou and Pecube.ptt) will be stored.

      read (7,'(a)') run

! second line
! npe: number of nodes per surface (2D) elements (3 = triangular elements - 4 = rectangular elements)
! nsurf: number of nodes defining the surface geometry
! nzin: number of nodes in the vertical direction
! nelemsurf: number of surface 2D elements
! zl: thickness of crustal layer (= depth at which t is fixed) (in km)
! diffusivity: thermal diffusivity (in km^2/Myr)
! heatproduction: heat production (in degC/Myr)

      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction,friction
      vivi=.FALSE.
        if (nsurf.lt.0) then
        nsurf=-nsurf
        vivi=.TRUE.
        endif
      mpe=2*npe

! second line:
! tmax: basal temperature (at z=-zl) (in degC)
! tmsl: temperature at mean sea level (z=0) (in degC)
! tlapse: lapse rate (indegC/km) 
! nstep: number of stages (the surface topo and exhumatoin rate can be specified at each time stage)
! ilog: redundant (dont use)
! iterative: iterative solver method (1 = Gauss Siedel with overrelaxation - 2 = Conjugate gradient)
! interpol: interpolation used (1 = linear - 2 = quadratic, not recommended)

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso
      read (7,*) xstep,ystep,young,poisson,thickness

      allocate (iconsurf(npe,nelemsurf),ielsurf(nsurf))
      allocate (neighbour(npe,nelemsurf))
      allocate (xsurf(nsurf),ysurf(nsurf),zsurf(nsurf),zsurfp(nsurf))
      allocate (usurf(nsurf),usurfp(nsurf),tempsurf(nsurf),tempsurfp(nsurf)) !VKP
      allocate (topoa(nsurf),topob(nsurf),rsurf(nsurf))
      allocate (tprev(nsurf,nstep))
      allocate (xdepth(nsurf,nstep),ydepth(nsurf,nstep),zdepth(nsurf,nstep))
      allocate (xexhumation(nsurf,nstep),yexhumation(nsurf,nstep),zexhumation(nsurf,nstep))

      tprev=0.d0
      usurf=0.d0 !VKP
      usurfp=0.d0 !VKP
      tempsurf=0.d0 !VKP
      tempsurfp=0.d0 !VKP

      if (ilog.eq.1) open (9,file='Pecube.log',status='unknown')

! read in nodal geometry

      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax

! read in topo advection info

      read (7,*) vtopo,advect_topo

! read in misfit type and # of parameters

      read (7,*) misfit_type,ndmisfit

! xsurf, ysurf: x- and y-locations of the surface nodes
! iconsurf: connectivity matrix describing the 2D surface elements

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

          do kstep=0,nstep
          read (7,*) timesurfp,iout
          read (7,*) (zsurfp(i),i=1,nsurf)
          if (vivi) read (7,*) (usurfp(i),i=1,nsurf)  !VKP
          if (vivi) read (7,*) (tempsurfp(i),i=1,nsurf)  !VKP
          enddo

      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Advecting rocks'

      xmin=minval(xsurf)
      xmax=maxval(xsurf)
      ymin=minval(ysurf)
      ymax=maxval(ysurf)

! go through the input file a first time to determine the depth of the
! points that will end up at the surface at the end of the model run

      niteradvec=3
      vxtopo=0.d0
      vytopo=0.d0
      vxtopo0=0.d0
      do istep=1,nstep
      xexhumation(1:nsurf,istep)=xsurf
      yexhumation(1:nsurf,istep)=ysurf
      zexhumation(1:nsurf,istep)=zl
      enddo
        do istep=nstep,1,-1
        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Time step ',istep
        timesurfp=0.
        rewind (7)
        read (7,*)
        read (7,*) i1,i2,i3,i4,x1,x2,x3
        read (7,*) x1,x2,x3,i1,i2,i3,i4
        read (7,*) x1,x2,x3
        read (7,*) i1,i2,i3,i4
        read (7,*) x1,x2,x3,x4,x5
        read (7,*) x1,x2,x3,x4
        read (7,*) x1,i1
        read (7,*) x1,x2
        read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
        read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)
          do kstep=0,istep-1
          read (7,*) timesurfp,iout
          read (7,*) (zsurfp(i),i=1,nsurf)
          if (vivi) read (7,*) (usurfp(i),i=1,nsurf)  !VKP
          if (vivi) read (7,*) (tempsurfp(i),i=1,nsurf)  !VKP
          enddo
        read (7,*) timesurf,iout
        read (7,*) (zsurf(i),i=1,nsurf)
        if (vivi) read (7,*) (usurf(i),i=1,nsurf) !VKP
        if (vivi) read (7,*) (tempsurf(i),i=1,nsurf) !VKP

! Calculate proper z-node spacing if nz is zero
! cspath and dwhipp 11/07
        if (nzin.eq.0) then                                                     ! If number of input z levels is zero, then code will calculate nz and spacing
          zmin=minval(zsurf)                                                    !   of z levels.  z levels will be spaced 1:1 with the input x and y spacing
          zmin_orig=zmin                                                        !   down to ~5 km below the model surface, 3:1 down to ~15 km below the surface
          zmax=maxval(zsurf)                                                    !   and ~9:1 for the rest of the model.
          xy_mean=(xstep+ystep)/2                                               ! This first portion of the variable z spacing code determines the number of
          cur_depth=zmin+zl                                                     !   z levels for the new geometry.
          nz=1
          depth_flag=.true.
          do while (cur_depth.gt.0.)                                            ! Work down from min elevation to base determining number of z levels needed
            nz=nz+1                                                             ! While still in this loop, increment nz
            if (cur_depth.gt.(zmin+zl)-5.) then                                 ! If within 5 km of model top surface, space node planes at xy_mean (1:1)
              cur_depth=cur_depth-(xy_mean/1000.)                               ! Subtract off new node plane spacing from remaining depth range
              plane_store1=nz-1                                                 ! Store number of planes used at 1:1 spacing
            else if (cur_depth.gt.(zmin+zl)-15.) then                           ! If within 15 km of model surface, space node planes at 3*xy_mean (3:1)
              cur_depth=cur_depth-3.*(xy_mean/1000.)                            ! Subtract off new node plane spacing from remaining depth range
              plane_store2=nz-1                                                 ! Store number of planes used at 3:1 spacing
            else                                                                ! If greater than 15 km from model surface, use ~9:1 node plane spacing
              if (depth_flag) then                                              ! Calculate node plane spacing on first click through this condition
                depth_flag=.false.                                              ! Set depth_flag to false to avoid repeating this calculation
                num_left=int((cur_depth)/(9.*xy_mean/1000.))+1                  ! Number of remaining node planes is equal to the remaining model depth over
                fact=cur_depth/(num_left*(xy_mean/1000.))                       !   the 9:1 spacing increment, plus one.  This yields spacing of <=9:1.
              endif
              if (int(cur_depth-fact*(xy_mean/1000.)).eq.0.) then               ! If near the base of the model, set the cur_depth to zero
                cur_depth=0.
              else                                                              ! Subtract off new node spacing from remaining depth range
                cur_depth=cur_depth-fact*(xy_mean/1000.)
              endif
            endif
          enddo
        else if (nzin.gt.0) then                                                ! If number of input z levels is positive, then use that number for nz
          nz=nzin
        else                                                                    ! Stop if input nz value is negative
          print *,'Error in Pecube.in: nz must be zero or a positive integer'
          stop
        endif

        if (istep.eq.nstep) then
          read (7,*) nobs
          allocate (xdeptho(nobs),ydeptho(nobs),zdeptho(nobs))
          allocate (zsurfo(nobs),zsurfpo(nobs))
          allocate (xexhumationo(nobs),yexhumationo(nobs),zexhumationo(nobs))
          allocate (tprevo(nobs),ieo(nobs),we1(nobs),we2(nobs),we3(nobs),we4(nobs))
          do i=1,nobs
            read (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, & ! By Xav
                       ageheZobs,dageheZobs,ageftZobs,dageftZobs,agearKobs,dagearKobs, &
                       agearBobs,dagearBobs,agearMobs,dagearMobs,agearHobs,dagearHobs, &
                       tlobs,Tramanobs,dTramanobs,heightobs,ieo(i),we1(i),we2(i),we3(i),we4(i)   ! By Xav

            xexhumationo(i)=we1(i)*xsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*xsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*xsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*xsurf(iconsurf(npe,ieo(i)))
            yexhumationo(i)=we1(i)*ysurf(iconsurf(1,ieo(i))) &
                           +we2(i)*ysurf(iconsurf(2,ieo(i))) &
                           +we3(i)*ysurf(iconsurf(3,ieo(i))) &
                           +we4(i)*ysurf(iconsurf(npe,ieo(i)))
            zexhumationo(i)=we1(i)*zsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*zsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*zsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*zsurf(iconsurf(npe,ieo(i)))+zl+min(0.d0,heightobs/1.e3)
          enddo
          read (7,*) ageflag,calc_surface_ages,PTt_output
        endif
! isostatic rebound
        topoa=zsurfp
        topob=zsurf
          if (isoflag.eq.1) then
          call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf, &
                                  rhoc,rhom,nx,ny,nxiso,nyiso, &
                                  xstep,ystep,young,poisson,thickness)
          else
          rsurf=0.
          endif
        zexhumation(1:nsurf,istep)=zexhumation(1:nsurf,istep)+zsurf
        if (istep.eq.nstep) tfinal=timesurf
        call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                      nz,timesurf,timesurfp,istep,eps,ilog, &
                      dt,ntime,istatic,fault,nfault,usurf)
        if (nd.eq.0) print*,ntime,'time steps required'
          do ktime=1,ntime
            time=timesurf-dt*ktime
            if (advect_topo == 1 .and. ktime == 1 .and. istep == nstep) then
              call calculate_topo_advection (time,fault,nfault,1,zl,         &
                                             advect_topo,vtopo,vxtopo,vytopo,&
                                             vxtopo0,iproc,nd)
            endif

            ! All surface particles loop - Track all surface particles if 
            ! calc_surface_ages flag is true - dwhipp 08/12
            if (calc_surface_ages) then
              do kstep=istep,nstep
                do i=1,nsurf
                  xxx=xexhumation(i,kstep)
                  yyy=yexhumation(i,kstep)
                  zzz=zexhumation(i,kstep)
                  do k=1,niteradvec
                    call find_velo (xxx,yyy,zzz,time,-dt,fault,nfault,vx,vy,vz,  &
                                    xmin,xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
                    dxxx=xexhumation(i,kstep)-dt*vx/2.-xxx
                    dyyy=yexhumation(i,kstep)-dt*vy/2.-yyy
                    dzzz=zexhumation(i,kstep)-dt*vz/2.-zzz
                    xxx=xxx+dxxx
                    yyy=yyy+dyyy
                    zzz=zzz+dzzz
                  enddo
                  xexhumation(i,kstep)=xexhumation(i,kstep)-dt*vx
                  yexhumation(i,kstep)=yexhumation(i,kstep)-dt*vy
                  zexhumation(i,kstep)=zexhumation(i,kstep)-dt*vz-rsurf(i)/ntime-dt*usurf(i) !VKP
                enddo
              enddo
            endif

            ! Tracked particles loop
            do i=1,nobs
              xxx=xexhumationo(i)
              yyy=yexhumationo(i)
              zzz=zexhumationo(i)
              do k=1,niteradvec
                call find_velo (xxx,yyy,zzz,time,-dt,fault,nfault,vx,vy,vz,    &
                                xmin,xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
                dxxx=xexhumationo(i)-dt*vx/2.-xxx
                dyyy=yexhumationo(i)-dt*vy/2.-yyy
                dzzz=zexhumationo(i)-dt*vz/2.-zzz
                xxx=xxx+dxxx
                yyy=yyy+dyyy
                zzz=zzz+dzzz
              enddo
              xexhumationo(i)=xexhumationo(i)-dt*vx
              yexhumationo(i)=yexhumationo(i)-dt*vy
              zexhumationo(i)=zexhumationo(i)-dt*vz-we1(i)*(rsurf(iconsurf(1,ieo(i)))/ntime+dt*usurf(iconsurf(1,ieo(i)))) &
                                                   -we2(i)*(rsurf(iconsurf(2,ieo(i)))/ntime+dt*usurf(iconsurf(2,ieo(i)))) &
                                                   -we3(i)*(rsurf(iconsurf(3,ieo(i)))/ntime+dt*usurf(iconsurf(3,ieo(i)))) &
                                                   -we4(i)*(rsurf(iconsurf(npe,ieo(i)))/ntime+dt*usurf(iconsurf(npe,ieo(i))))
            enddo

          ! Update fault positions if using topographic advection
          if (abs(advect_topo) == 1) then
            do i=1,nfault
              do k=1,fault(i)%n
                fault(i)%x(k)=fault(i)%x(k)+dt*vxtopo0
              enddo
            enddo
          endif

          ! Fault advection call
          if (fltflag.eq.1) call move_fault (fault,faultp,nfault,-dt,time,xmin,&
                                             xmax,ymin,ymax,vxtopo,vytopo,     &
                                             vxtopo0)
          enddo

        zsurfp=zsurf
        timesurfp=timesurf
        enddo

! depth is the initial depth of the rocks that end up at the location of the
! surface nodes at the end of the model experiment

      xdepth=xexhumation
      ydepth=yexhumation
      zdepth=zexhumation
      xdeptho=xexhumationo
      ydeptho=yexhumationo
      zdeptho=zexhumationo

! reset the input file to its proper position

!      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'initializing'
      if (iproc.eq.0 .and. nd.eq.0) print*,''

      rewind (7)
      read (7,'(a)') run
      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction
      nsurf=iabs(nsurf)
      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      nstep=abs(nstep)
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso
      read (7,*) xstep,ystep,young,poisson,thickness
      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax
      read (7,*) vtopo,advect_topo
      read (7,*) misfit_type,ndmisfit
      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      nnode=nsurf*nz
      nelem=nelemsurf*(nz-1)
      if (ilog.eq.1) write (9,*) 'nnode/nelem= ',nnode,nelem

! opens output files

! Pecube.out contains the temperature field at the end of each stage

      if (iproc.eq.0.and.nd.eq.0) open (8,file=run//'/Pecube.out',status='unknown',access='direct', &
            recl=4*(4+7*nnode+mpe*nelem))

! Ages.out contains the ages at the end of each stage

      if (iproc.eq.0.and.nd.eq.0) open (11,file=run//'/Ages.out',status=   &
            'unknown',access='direct', recl=4*(4+14*nsurf+npe*nelemsurf))

! Pecube.ptt contains the depth-temperture-paths of all surface nodes
! Files created only if calc_surface_ages is true - dwhipp 08/12

      if (calc_surface_ages) then
        do istep=1,nstep
          open (100+istep,status='scratch',access='direct',recl=4*(1+nsurf*4))
        enddo
      endif
      open (200,status='scratch',access='direct',recl=4*(1+nobs*4))

      allocate (x(nnode),y(nnode),z(nnode),t(nnode))
      allocate (u(nnode)) !VKP
      allocate (xp(nnode),yp(nnode),zp(nnode),tp(nnode))
      allocate (up(nnode)) !VKP
      allocate (icon(mpe,nelem))
      allocate (kfix(nnode))
      allocate (f(nnode),rebound(nnode))

      u=0.d0
      up=0.d0

! build 3D element connectivity

      ie=0
        do iesurf=1,nelemsurf
          do k=1,(nz-1)
          ie=ie+1
            do kk=1,npe
            icon(kk,ie)=(iconsurf(kk,iesurf)-1)*nz+k
            icon(kk+npe,ie)=icon(kk,ie)+1
            enddo
          enddo
        enddo

      if (ie.ne.nelem) then
      stop 'nelem mismatch'
      endif

! finds processor topology
      
      allocate (proc(nnode),eproc(nelem))
      call define_proc (nproc,nsurf,nz,proc)
        do ie=1,nelem
        ido=0
          do k=1,mpe
          if (proc(icon(k,ie)).eq.iproc) ido=1
          enddo
        eproc(ie)=ido
        enddo

! allocate reduced number of elemental matrices

      nelemloc=0
        do ie=1,nelem
          if (eproc(ie).eq.1) then
          nelemloc=nelemloc+1
          endif
        enddo

      allocate (ael(mpe,mpe,nelemloc),bel(mpe,nelemloc),ice(nelemloc))

      je=0
        do ie=1,nelem
          if (eproc(ie).eq.1) then
          je=je+1
          ice(je)=ie
          endif
        enddo

      if (ilog.eq.1) then
      write (9,*) 'icon'
        do ie=1,nelem
        write (9,*) (icon(k,ie),k=1,mpe)
        enddo
      endif

! finds neighbour conectivity matrix

      call find_neighbours (iconsurf,neighbour,npe,nelemsurf,nsurf)
      ielsurf=1

! initialize global parameters
! alpha is the time integration parameter

      alpha=0.5
      time=0.
      timesurf=0.
      irec=0
      jrec=0

      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Start of time stepping'

! begining of surface stepping (stageing)

      do istep=0,nstep

        if (ilog.eq.1) write (9,*) 'istep= ',istep

        if (iproc.eq.0 .and. nd.eq.0) write (6,*) 'Stage : ',istep,' of ',nstep

        if (istep.ne.0) zsurfp=zsurf
        if (istep.ne.0) usurfp=usurf !VKP
        if (istep.ne.0) tempsurfp=tempsurf !VKP
        timesurfp=timesurf

! read the step information

! timesurf is the time since the begining of the experiment
! Peclet is the exhumation velocity since the end of the last stage
! iout indicates if the T information is to be saved at the end of the stage

        read (7,*) timesurf,iout
        if (istep.eq.0) then
          if (timesurf.gt.eps.and.iproc.eq.0) then
            write (6,*) 'timesurf= ',timesurf
            write (6,*) 'first topography record must be at time zero ...'
            stop
          endif
        endif
        read (7,*) (zsurf(i),i=1,nsurf)
        if (vivi) read (7,*) (usurf(i),i=1,nsurf) !VKP
        if (vivi) read (7,*) (tempsurf(i),i=1,nsurf) !VKP

! initial (or zeroth) step

        if (istep.eq.0) then
          zsurfp=zsurf
          usurfp=usurf !VKP
          tempsurfp=tempsurf !VKP
          zsurfmax=maxval(zsurf)
          in=0
          if (nzin.gt.0) then                                                     ! If the input nz value is positive, build geometry with constant z node plane
            do i=1,nsurf                                                          !   spacing equal to the total model thickness divided by the input nz
              do k=1,nz
                in=in+1
                fact=float(k-1)/float(nz-1)
                if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
                xp(in)=xsurf(i)
                x(in)=xp(in)
                yp(in)=ysurf(i)
                y(in)=yp(in)
                zh=zsurf(i)+zl
                zp(in)=zh*fact
                up(in)=usurf(i)
                u(in)=up(in)
                kfix(in)=0
                ! Modified to identify top and bottom fixed T B/Cs (original
                ! version below) - dwhipp (09/07)
                if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                ! Commented out for now - dwhipp 08/12
                !if (k.eq.1) kfix(in)=1
                !if (k.eq.nz) kfix(in)=2
              enddo
            enddo
          else                                                                    ! If the input nz value is zero, build geometry with variable z node plane
            zmin=minval(zsurf)                                                    !   spacing.  The spacing is such that from the top surface of the model down
            zmax=maxval(zsurf)                                                    !   to 5 km below it, the spacing is 1:1 with the x and y node spacing (xy_mean)
            do i=1,nsurf                                                          !   Below that, down to 15 km below the surface, the spacing is ~3:1 (3*xy_mean)
              do k=1,nz                                                           !   Below that, the spacing is ~9:1 down to the base of the model
                in=in+1                                                           ! Set the topography scaling factor to zero if the top surface is flat
                if ((zmax-zmin).eq.0.) then
                  scale_fact=0.
                else                                                              ! Define scaling factor to shift nodes beneath the topography
                  scale_fact=(zsurf(i)-zmin)*(1/(zsurf(i)+zl))
                endif
                if (k.gt.nz-(plane_store1)) then
                  if (k.eq.nz) then
                    scale_fact=0.
                    cur_depth=zsurf(i)+zl
                  else
                    cur_depth=cur_depth+(xy_mean/1000.)
                  endif
                else if (k.gt.nz-(plane_store2)) then
                  cur_depth=cur_depth+3*(xy_mean/1000.)
                else
                  if (k.eq.1) then
                    cur_depth=0.
                  else
                    cur_depth=cur_depth+fact*(xy_mean/1000.)
                  endif
                endif
                xp(in)=xsurf(i)
                x(in)=xp(in)
                yp(in)=ysurf(i)
                y(in)=yp(in)
                zp(in)=cur_depth+cur_depth*scale_fact
                up(in)=usurf(i)
                u(in)=up(in)
                kfix(in)=0
                ! Modified to identify top and bottom fixed T B/Cs (original
                !   version below) - dwhipp (09/07)
                if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                ! Commented out for now - dwhipp 08/12
                !if (k.eq.1) kfix(in)=1
                !if (k.eq.nz) kfix(in)=2
              enddo
            enddo
          endif

! calculates initial temperature

          do i=1,nsurf
            tsurf=-zsurf(i)*tlapse+tmsl
            if (vivi) tsurf = tempsurf(i)
            do k=1,nz
              in=(i-1)*nz+k
              zh=zsurf(i)+zl
              tp(in)=tsurf+(tmax-tsurf)*(zh-zp(in))/zh
            enddo
          enddo

          t=tp

          if (ilog.eq.1) then
          write (9,*) 'Nodal geometry'
            do i=1,nnode
              write (9,*) i,xp(i),yp(i),zp(i),tp(i),kfix(i)
            enddo
          endif

        endif

        call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp,nz,timesurf,timesurfp, &
                     istep,eps,ilog,dt,ntime,istatic,fault,nfault,usurf)

! beginning of time stepping

        do itime=1,ntime
          call cpu_time (times(3))
          time=time+dt
          ftime=float(itime)/float(ntime)
          ftimep=float(itime-1)/float(ntime)
! Mod made by Jean on 15/3/2010
! now one should use tau=0 to indicate a linear change in topography
! a positive value for a change late in the time step
! a negative value for a change early in the time step (most meaningful geomorphologically)
          if (tau.ne.0.d0) then
            ftime=(1.-exp(-ftime*tau/tfinal))/(1.-exp(-tau/tfinal))
            ftimep=(1.-exp(-ftimep*tau/tfinal))/(1.-exp(-tau/tfinal))
		  endif

!        if (iproc.eq.0.and.nd.eq.0) write (6,'(a,i4,a,i4,a,i3,a,i3,a,i6,a)') &
!         'Doing time step ',itime,' of ',ntime, &
!         ' in stage ',istep,' of ',nstep

          if (ilog.eq.1) write (9,*) 'itime,time= ',itime,time

! build new node geometry

          z=zp

          do i=1,nsurf
            in=i*nz
            inp=i*nz-1
            zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime+zl
            zsurf2=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep+zl
            topoa(i)=zsurf2
            topob(i)=zsurf1
            z(in)=z(in)+zsurf1-zsurf2
            z(inp)=z(inp)+(zsurf1-zsurf2)/2.
! Mod made by Jean on 23/3/2010
! This insures that the surface temperature follows the lapse rate at every time step
! not just at t=0
            tp(in)=-(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep)*tlapse+tmsl
            t(in)=-(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)*tlapse+tmsl
            if (vivi) then !VKP
              tempsurf1 = tempsurfp(i) + (tempsurf(i)-tempsurfp(i))*ftime !VKP
              tempsurf2 = tempsurfp(i) + (tempsurf(i)-tempsurfp(i))*ftimep !VKP
!              t(in) = t(in) + tempsurf1 - tempsurf2 !VKP
              tp(in)=tempsurf2
              t(in)=tempsurf1
            endif !VKP
          enddo

          in=0
          do i=1,nsurf !VKP
            do k=1,nz !VKP
              in=in+1 !VKP
              u(in)=usurf(i) !VKP
            enddo !VKP
         enddo !VKP

         if (isoflag.eq.1) then
           call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf,rhoc,   &
                                  rhom,nx,ny,nxiso,nyiso,xstep,ystep,young,    &
                                  poisson,thickness)
         else
           rsurf=0.
         endif

         do i=1,nsurf
           do k=1,nz
             rebound(k+(i-1)*nz)=rsurf(i)
           enddo
         enddo

        if (in.ne.nnode) then
          stop 'nnode mismatch'
        endif

! build local FE matrices

!      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Building matrix '

        call cpu_time (times(7))
        je=0
        do je=1,nelemloc
          ie=ice(je)
          call make_matrix (mpe,ael(1,1,je),bel(1,je),icon(1,ie),x,y,z,u,xp,yp,&
                            zp,up,kfix,diffusivity,heatproduction,x1f,y1f,x2f, &! VKP
                            y2f,def,dif,alpha,dt,time,tp,nnode,istatic,fault,  &
                            nfault,rebound,xmin,xmax,ymin,ymax,friction,vxtopo,&
                            vytopo,vxtopo0)
        enddo
        call cpu_time (times(8))

!        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Building time : ',times(8)-times(7)

! build global RHS vector

        f=0.d0
        do je=1,nelemloc
          ie=ice(je)
          do k=1,mpe
            ic=icon(k,ie)
            f(ic)=f(ic)+bel(k,je)
          enddo
        enddo

! solve global FE equations

!        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Solving matrix '

        call cpu_time (times(5))
        if (iterative.eq.1) then
          call solve_iterative (mpe,1,ael,f,t,kfix,icon,nnode,nelemloc,niter,  &
                               proc,ice)
          if (ilog.eq.1) write (9,*) niter,' iterations'
        elseif (iterative.eq.2) then
          stop 'solution strategy not implemented'
        else
          stop 'solution strategy not implemented'
        endif
        call cpu_time (times(6))

!        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Solution time : ',times(6)-times(5),niter

        if (iproc.eq.0.and.nd.eq.0) call screen_counter (itime,ntime,niter)

! stretch old grid - now with variable z-node spacing
        if (nzin.gt.0) then                                                       ! Linearly shift nodes if the input nz value is positive
          do i=1,nsurf
            do k=1,nz
              in=(i-1)*nz+k
              fact=float(k-1)/float(nz-1)
              if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
              zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime
              zh=zsurf1+zl
              zp(in)=zh*fact
            enddo
          enddo
        else                                                                      ! Shift variably spaced z node planes for new topography
          zmin=minval(zsurfp+(zsurf-zsurfp)*ftime)                                ! Calculate new zmin and zmax for new topo
          zmax=maxval(zsurfp+(zsurf-zsurfp)*ftime)
          do i=1,nsurf
            do k=1,nz
              in=(i-1)*nz+k
              if ((zmax-zmin).eq.0.) then                                         ! Set scale_fact to zero if surface topo is flat
                if (zmin.eq.0.) then
                  scale_fact=0.
                else
                  scale_fact=(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)*(1/zl)        ! Set scale factor to allow for topography shifts
                endif
              else                                                                ! Set scale factor to allow for topography shifts and relief scaling
                scale_fact=((zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)-zmin_orig)*(1/zl)
              endif
              if (k.gt.nz-(plane_store1)) then
                if (k.eq.nz) then
                  cur_depth=(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)+zl
                  scale_fact=0.
                else
                  cur_depth=cur_depth+(xy_mean/1000.)
                endif
              else if (k.gt.nz-(plane_store2)) then
                cur_depth=cur_depth+3.*(xy_mean/1000.)
              else
                if (k.eq.1) then
                  cur_depth=0.
                else
                  cur_depth=cur_depth+fact*(xy_mean/1000.)
                endif
              endif
              zp(in)=cur_depth+cur_depth*scale_fact
            enddo
          enddo
        endif

! Update fault positions if using topographic advection
        if (abs(advect_topo) == 1) then
          do i=1,nfault
            do k=1,fault(i)%n
              fault(i)%x(k)=fault(i)%x(k)-dt*vxtopo0
            enddo
          enddo
        endif

! update fault geometry
        if (fltflag.eq.1) call move_fault (fault,faultp,nfault,dt,time,xmin,   &
                                          xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)

! interpolate result onto undeformed mesh
        do i=1,nsurf
          ij=(i-1)*nz+1
          call interpolate (t(ij),tp(ij),z(ij),zp(ij),nz)
        enddo

! update ptt

        zsurfo=we1*zsurf(iconsurf(1,ieo))+we2*zsurf(iconsurf(2,ieo))+ &
               we3*zsurf(iconsurf(4,ieo))+we4*zsurf(iconsurf(4,ieo))
        zsurfpo=we1*zsurfp(iconsurf(1,ieo))+we2*zsurfp(iconsurf(2,ieo))+ &
                we3*zsurfp(iconsurf(4,ieo))+we4*zsurfp(iconsurf(4,ieo))

        if (calc_surface_ages) then
          do kstep=max(1,istep),nstep
            do i=1,nsurf
              xxx=xdepth(i,kstep)
              yyy=ydepth(i,kstep)
              zzz=zdepth(i,kstep)
              do k=1,niteradvec
                call find_velo (xxx,yyy,zzz,time,dt,fault,nfault,vx,vy,vz,xmin,&
                               xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
                dxxx=xdepth(i,kstep)+dt*vx/2.-xxx
                dyyy=ydepth(i,kstep)+dt*vy/2.-yyy
                dzzz=zdepth(i,kstep)+dt*vz/2.-zzz
                xxx=xxx+dxxx
                yyy=yyy+dyyy
                zzz=zzz+dzzz
              enddo
              xdepth(i,kstep)=xdepth(i,kstep)+vx*dt
              ydepth(i,kstep)=ydepth(i,kstep)+vy*dt
              zdepth(i,kstep)=zdepth(i,kstep)+vz*dt+rsurf(i)+dt*usurf(i) ! VKP
              call find_element (xdepth(i,kstep),ydepth(i,kstep),              &
                                zdepth(i,kstep),tnow,x,y,z,t,xsurf,ysurf,zsurf,&
                                ielsurf(i),neighbour,iconsurf,icon,nelemsurf,  &
                                nelem,nsurf,nz,nnode,npe,mpe,interpol_success)
              if (interpol_success) then
                tprev(i,kstep)=tnow
              else
                if (zdepth(i,kstep).gt.zl) then
                  tprev(i,kstep)=0.d0
                elseif (zdepth(i,kstep).lt.0.d0) then
                  tprev(i,kstep)=tmax
                else
!                  tprev(i,kstep)=tmax-tmax*zdepth(i,kstep)/zl
                  tprev(i,kstep)=tmax-tmax*zdepth(i,kstep)/(zl+zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)
                endif
              endif
            enddo
          enddo
        endif
        do i=1,nobs
          xxx=xdeptho(i)
          yyy=ydeptho(i)
          zzz=zdeptho(i)
          do k=1,niteradvec
            call find_velo (xxx,yyy,zzz,time,dt,fault,nfault,vx,vy,vz,xmin,    &
                            xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
            dxxx=xdeptho(i)+dt*vx/2.-xxx
            dyyy=ydeptho(i)+dt*vy/2.-yyy
            dzzz=zdeptho(i)+dt*vz/2.-zzz
            xxx=xxx+dxxx
            yyy=yyy+dyyy
            zzz=zzz+dzzz
          enddo
          xdeptho(i)=xdeptho(i)+vx*dt
          ydeptho(i)=ydeptho(i)+vy*dt
          zdeptho(i)=zdeptho(i)+vz*dt+we1(i)*(rsurf(iconsurf(1,ieo(i)))/ntime+dt*usurf(iconsurf(1,ieo(i)))) &
                                     +we2(i)*(rsurf(iconsurf(2,ieo(i)))/ntime+dt*usurf(iconsurf(2,ieo(i)))) &
                                     +we3(i)*(rsurf(iconsurf(3,ieo(i)))/ntime+dt*usurf(iconsurf(3,ieo(i)))) &
                                     +we4(i)*(rsurf(iconsurf(npe,ieo(i)))/ntime+dt*usurf(iconsurf(npe,ieo(i))))

          call find_element (xdeptho(i),ydeptho(i),zdeptho(i),tnow,x,y,z,t,    &
                             xsurf,ysurf,zsurf,ielsurf(i),neighbour,iconsurf,  &
                             icon,nelemsurf,nelem,nsurf,nz,nnode,npe,mpe,      &
                             interpol_success)
          if (interpol_success) then
            tprevo(i)=tnow
          else
            if (zdeptho(i).gt.zl) then
              tprevo(i)=0.d0
            elseif (zdeptho(i).lt.0.d0) then
              tprevo(i)=tmax
            else
!              tprevo(i)=tmax-tmax*zdeptho(i)/zl
              tprevo(i)=tmax-tmax*zdeptho(i)/(zl+zsurfpo(i)+(zsurfo(i)-zsurfpo(i))*ftime)
            endif
          endif
        enddo

!      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'PTt Updated'

! saves tt and pt

        jrec=jrec+1
        tfake=tfinal-time
        if (jrec.eq.1) tfake=tprevious-time
        if (calc_surface_ages) then
          do kstep=max(1,istep),nstep
            write (100+kstep,rec=jrec) sngl(tfake),sngl(tprev(1:nsurf,kstep)), &
                  sngl(zl+zsurfp(1:nsurf)+(zsurf(1:nsurf)-zsurfp(1:nsurf))*ftime-zdepth(1:nsurf,kstep)), &
                  sngl(xdepth(1:nsurf,kstep)),sngl(ydepth(1:nsurf,kstep))
          enddo
        endif
        write (200,rec=jrec) sngl(tfake),sngl(tprevo(1:nobs)), &
              sngl(zl+zsurfpo(1:nobs)+(zsurfo(1:nobs)-zsurfpo(1:nobs))*ftime-zdeptho(1:nobs)), &
              sngl(xdeptho(1:nobs)),sngl(ydeptho(1:nobs))

! end of time stepping

        call cpu_time (times(4))

!        if (iproc.eq.0.and.nd.eq.0) write(6,*) 'This time step : ',times(4)-times(3)

        enddo

      if (iproc.eq.0 .and. nd.eq.0) print*,''

! clean up

      if (calc_surface_ages) then
        if (istep.ne.0) then
          allocate(t2(nsurf),t2p(nsurf),t3(nsurf),t4(nsurf),t5(nsurf))
          t2p=0.
          do krec=jrec,1,-1
            read (100+istep,rec=krec) t1,t2,t3,t4,t5
            do i=1,nsurf
              if (t2(i).lt.0.) t2(i)=t2p(i)
              t2p(i)=t2(i)
            enddo
            write (100+istep,rec=krec) t1,t2,t3,t4,t5
          enddo
          deallocate (t2,t2p,t3,t4,t5)
  
  !      jrec=jrec+1
          write (100+istep,rec=jrec+1) sngl(-1.d0),sngl(tprev(1:nsurf,istep)), &
               sngl(zl+zsurfp+(zsurf-zsurfp)*ftime-zdepth(1:nsurf,istep)), &
               sngl(xdepth(1:nsurf,istep)),sngl(ydepth(1:nsurf,istep))
  
        endif

! calculates ages

        allocate (age1(nsurf),age2(nsurf),age3(nsurf),age4(nsurf),age5(nsurf))
        allocate (age6(nsurf),age7(nsurf),age8(nsurf),ftdist(17,nsurf))
        allocate (raman(nsurf))

        if (istep.eq.0.or.iout.eq.0) then

          age1=tprevious;age2=tprevious;age3=tprevious;age4=tprevious;age5=tprevious
          age6=tprevious;age7=tprevious;age8=tprevious;ftdist=0.;raman=0.

        else

          if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Calculating ages'

          call calculate_ages (jrec,nsurf,nz,istep,age1,age2,age3,age4,age5,   &
                               age6,age7,age8,ftdist,raman,tprevious,ftlflag,  &
                               ageflag,iproc,nd,calc_surface_ages,PTt_output)

        endif
      endif

! write output

      if (iout.eq.1) then
        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Saving at step ',istep
        irec=irec+1

        allocate (vxx(nnode),vyy(nnode),vzz(nnode))
        ij=0
        iu=0 !VKP
        do i=1,nsurf
          iu=iu+1 !VKP
          do k=1,nz
            ij=ij+1
            call find_velo (x(ij),y(ij),z(ij),time,dt,fault,nfault,vx,vy,vz, &
                            xmin,xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
            vxx(ij)=vx;vyy(ij)=vy;vzz(ij)=vz+usurf(iu) !VKP
          enddo
        enddo
        if (nd.eq.0) write (8,rec=irec) nnode,nelem,mpe,sngl(tfinal-time),sngl(x),sngl(y),sngl(z),sngl(t), &
                           vxx,vyy,vzz,((icon(j,i),j=1,mpe),i=1,nelem)
        deallocate (vxx,vyy,vzz)

        if (calc_surface_ages) then
          allocate (exhumation(nsurf))
          exhumation=0.d0
          do i=1,nsurf
            dtimesurf=timesurf-timesurfp
            if (dtimesurf.gt.0.d0) exhumation(i)=-(zsurf(i)-zsurfp(i))/dtimesurf
            call find_velo (xsurf(i),ysurf(i),zl+zsurf(i),time,dt,fault,nfault,  &
                            vx,vy,vz,xmin,xmax,ymin,ymax,vxtopo,vytopo,vxtopo0)
            exhumation(i)=exhumation(i)+vz+usurf(i) !VKP
          enddo
          if (nd.eq.0) write (11,rec=irec) nsurf,nelemsurf,npe,sngl(tfinal-time),sngl(xsurf),sngl(ysurf),sngl(zl+zsurf), &
                              exhumation,age1,age2,age3,age4,age5,age6,age7,age8, &
                              (sngl(sum(length*ftdist(:,j))),j=1,nsurf),raman,iconsurf
          deallocate (exhumation)

          write (cstep,'(i3)') istep
          if (istep.lt.10) cstep(1:2)='00'
          if (istep.lt.100) cstep(1:1)='0'
          if (nd.eq.0) then
            open (12,file=run//'/Ages'//cstep//'.txt',status='unknown')
            write (12,'(a156)') 'Longitude   Latitude    Height      HeApatite   HeZircon    FTApatite   ' &
                              //'FTZircon    ArKFeldspar ArBiotite   ArMuscovite ArHornblend ' &
                              //'FTApMeanTL  RamanT      '
            do i=1,nsurf
              xxx=xlonmin+(xlonmax-xlonmin)*(xsurf(i)-xmin)/(xmax-xmin)
              yyy=xlatmin+(xlatmax-xlatmin)*(ysurf(i)-ymin)/(ymax-ymin)
              write (12,'(21f12.4)') xxx,yyy,zsurf(i),age1(i),age2(i),age3(i),age4(i),age5(i), &
                                           age6(i),age7(i),age8(i),sngl(sum(length*ftdist(:,i))), &
                                           raman(i)
            enddo
            close (12)
          endif
        endif
      endif

      if (calc_surface_ages) then
        if (istep.ne.nstep) deallocate (age1,age2,age3,age4,age5,age6,age7,age8,ftdist,raman)
        if (istep.ne.0) close (100+istep)
      else
        if (nd.eq.0) write (11,rec=irec) nsurf,nelemsurf,npe,                &
        sngl(tfinal-time),sngl(xsurf),sngl(ysurf),sngl(zl+zsurf),            &
        (0.0,i=1,nsurf),(0.0,i=1,nsurf),(0.0,i=1,nsurf),(0.0,i=1,nsurf),     &
        (0.0,i=1,nsurf),(0.0,i=1,nsurf),(0.0,i=1,nsurf),(0.0,i=1,nsurf),     &
        (0.0,i=1,nsurf),(0.0,i=1,nsurf),(0.0,i=1,nsurf),iconsurf
      endif

! end of surface stepping

      enddo

    zsurfo=we1*zsurf(iconsurf(1,ieo))+we2*zsurf(iconsurf(2,ieo))+ &
           we3*zsurf(iconsurf(4,ieo))+we4*zsurf(iconsurf(4,ieo))
    zsurfpo=we1*zsurfp(iconsurf(1,ieo))+we2*zsurfp(iconsurf(2,ieo))+ &
            we3*zsurfp(iconsurf(4,ieo))+we4*zsurfp(iconsurf(4,ieo))
    write (200,rec=jrec+1) sngl(-1.d0),sngl(tprevo(1:nobs)), &
           sngl(zl+zsurfpo(1:nobs)+(zsurfo(1:nobs)-zsurfpo(1:nobs))*ftime-zdeptho(1:nobs)), &
           sngl(xdeptho(1:nobs)),sngl(ydeptho(1:nobs))

      if (iproc.eq.0.and.nd.eq.0) then
      close (8)
      open (8,file=run//'/Pecube.out',status='unknown',access='direct', &
            recl=4)
      write (8,rec=irec*(4+7*nnode+mpe*nelem)+1) -1
      close (8)
      endif

      if (iproc.eq.0.and.nd.eq.0) then
        close (11)
        open (11,file=run//'/Ages.out',status='unknown',access='direct', &
              recl=4)
        write (11,rec=irec*(4+14*nsurf+npe*nelemsurf)+1) -1
        close (11)
      endif

      if (ilog.eq.1) close (9)

! calculate misfit

      if (iproc.eq.0) then
        ! added by Jean to change misfit to slope rather than exact ages
        ! (4/6/2008)
        open (13,file=run//'/Comparison.txt',status='unknown')
        misfit=0.
        nmisfit=0
        nhe=0
        nft=0
        nheZ=0
        nftZ=0
        narK=0
        narB=0
        narM=0
        narH=0
        nraman=0
        read (7,*) nobs
        allocate (aheoreg(nobs),ahepreg(nobs),hheoreg(nobs))
        allocate (aftoreg(nobs),aftpreg(nobs),hftoreg(nobs))
        allocate (aheZoreg(nobs),aheZpreg(nobs),hheZoreg(nobs))
        allocate (aftZoreg(nobs),aftZpreg(nobs),hftZoreg(nobs))
        allocate (aarKoreg(nobs),aarKpreg(nobs),harKoreg(nobs))
        allocate (aarBoreg(nobs),aarBpreg(nobs),harBoreg(nobs))
        allocate (aarMoreg(nobs),aarMpreg(nobs),harMoreg(nobs))
        allocate (aarHoreg(nobs),aarHpreg(nobs),harHoreg(nobs))
        allocate (aramanoreg(nobs),aramanpreg(nobs),hramanoreg(nobs))
        allocate (hpreg(nobs))
        allocate (hei(nobs),agehe(nobs),ageheZ(nobs),ageft(nobs),ageftZ(nobs),   &
                  agearB(nobs),agearK(nobs),agearM(nobs),agearH(nobs),           &
                  ftdisto(17,nobs),Traman(nobs))

        if (nd.eq.0) write (13,*) nobs
        call calculate_ages (jrec,nobs,nz,100,agehe,ageheZ,ageft,ageftZ,agearK,&
                             agearB,agearM,agearH,ftdisto,Traman,tprevious,    &
                             ftlflag,ageflag,iproc,nd,calc_surface_ages,       &
                             PTt_output)
        do iobs=1,nobs
          !read (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, &
          !           heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
          read (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, & ! By Xav
                     ageheZobs,dageheZobs,ageftZobs,dageftZobs,agearKobs,dagearKobs, &
                     agearBobs,dagearBobs,agearMobs,dagearMobs,agearHobs,dagearHobs, &
                     tlobs,Tramanobs,dTramanobs,heightobs,ieobs,wobs1,wobs2,wobs3,wobs4   ! By Xav
          hei(iobs)=wobs1*zsurf(iconsurf(1,ieobs)) &
                   +wobs2*zsurf(iconsurf(2,ieobs)) &
                   +wobs3*zsurf(iconsurf(3,ieobs)) &
                   +wobs4*zsurf(iconsurf(npe,ieobs)) 
          hei(iobs)=hei(iobs)*1000.+min(0.d0,heightobs)
          if (ageheobs.gt.0.) then
            nhe=nhe+1
            aheoreg(nhe)=ageheobs
            hheoreg(nhe)=hei(iobs)
          endif
          if (ageftobs.gt.0.) then
            nft=nft+1
            aftoreg(nft)=ageftobs
            hftoreg(nft)=hei(iobs)
          endif
          if (ageheZobs.gt.0.) then
            nheZ=nheZ+1
            aheZoreg(nheZ)=ageheZobs
            hheZoreg(nheZ)=hei(iobs)
          endif
          if (ageftZobs.gt.0.) then
            nftZ=nftZ+1
            aftZoreg(nftZ)=ageftZobs
            hftZoreg(nftZ)=hei(iobs)
          endif
          if (agearKobs.gt.0.) then
            narK=narK+1
            aarKoreg(narK)=agearKobs
            harKoreg(narK)=hei(iobs)
          endif
          if (agearBobs.gt.0.) then
            narB=narB+1
            aarBoreg(narB)=agearBobs
            harBoreg(narB)=hei(iobs)
          endif
          if (agearMobs.gt.0.) then
            narM=narM+1
            aarMoreg(narM)=agearMobs
            harMoreg(narM)=hei(iobs)
          endif
          if (agearHobs.gt.0.) then
            narH=narH+1
            aarHoreg(narH)=agearHobs
            harHoreg(narH)=hei(iobs)
          endif
          if (Tramanobs.gt.0.) then
            nraman=nraman+1
            aramanoreg(nraman)=Tramanobs
            hramanoreg(nraman)=hei(iobs)
          endif
          ahepreg(iobs)=agehe(iobs)
          aftpreg(iobs)=ageft(iobs)
          aheZpreg(iobs)=ageheZ(iobs)
          aftZpreg(iobs)=ageftZ(iobs)
          aarKpreg(iobs)=agearK(iobs)
          aarBpreg(iobs)=agearB(iobs)
          aarMpreg(iobs)=agearM(iobs)
          aarHpreg(iobs)=agearH(iobs)
          aramanpreg(iobs)=Traman(iobs)
          hpreg(iobs)=hei(iobs)

          call find_velo (xdeptho(iobs),ydeptho(iobs),zdeptho(iobs),time,dt,   &
                          fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,vxtopo,    &
                          vytopo,vxtopo0)

          vzobs=vz+we1(iobs)*(rsurf(iconsurf(1,ieo(iobs)))/ntime+dt*usurf(iconsurf(1,ieo(iobs)))) &
                  +we2(iobs)*(rsurf(iconsurf(2,ieo(iobs)))/ntime+dt*usurf(iconsurf(2,ieo(iobs)))) &
                  +we3(iobs)*(rsurf(iconsurf(3,ieo(iobs)))/ntime+dt*usurf(iconsurf(3,ieo(iobs)))) &
                  +we4(iobs)*(rsurf(iconsurf(npe,ieo(iobs)))/ntime+dt*usurf(iconsurf(npe,ieo(iobs))))

          if (nd.eq.0) write (13,'(128f12.4)') xlonobs,xlatobs,heightobs,      &
                              hei(iobs),vzobs,ageheobs,agehe(iobs),ageftobs,   &
                              ageft(iobs),ageheZobs,ageheZ(iobs),ageftZobs,    &
                              ageftZ(iobs),agearKobs,agearK(iobs),agearBobs,   &
                              agearB(iobs),agearMobs,agearM(iobs),agearHobs,   &
                              agearH(iobs),tlobs,ftdisto(1:17,iobs),Tramanobs, &
                              Traman(iobs)
          if (ageheobs.gt.0.) then
            !misfit=misfit+(ageheobs-agehe(iobs))**2/dageheobs**2
            misfit=misfit+age_misfit(ageheobs,agehe(iobs),dageheobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (ageftobs.gt.0.) then
            !misfit=misfit+(ageftobs-ageft(iobs))**2/dageftobs**2
            misfit=misfit+age_misfit(ageftobs,ageft(iobs),dageftobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (ageheZobs.gt.0.) then
            !misfit=misfit+(ageheZobs-ageheZ(iobs))**2/dageheZobs**2
            misfit=misfit+age_misfit(ageheZobs,ageheZ(iobs),dageheZobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (ageftZobs.gt.0.) then
            !misfit=misfit+(ageftZobs-ageftZ(iobs))**2/dageftZobs**2
            misfit=misfit+age_misfit(ageftZobs,ageftZ(iobs),dageftZobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (agearKobs.gt.0.) then
            !misfit=misfit+((agearKobs-agearK(iobs))/dagearKobs)**2
            misfit=misfit+age_misfit(agearKobs,agearK(iobs),dagearKobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (agearBobs.gt.0.) then
            !misfit=misfit+(agearBobs-agearB(iobs))**2/dagearBobs**2
            misfit=misfit+age_misfit(agearBobs,agearB(iobs),dagearBobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (agearMobs.gt.0.) then
            !misfit=misfit+(agearMobs-agearM(iobs))**2/dagearMobs**2
            misfit=misfit+age_misfit(agearMobs,agearM(iobs),dagearMobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (agearHobs.gt.0.) then
            !misfit=misfit+(agearHobs-agearH(iobs))**2/dagearHobs**2
            misfit=misfit+age_misfit(agearHobs,agearH(iobs),dagearHobs,misfit_type)
            nmisfit=nmisfit+1
          endif
          if (tlobs(1).gt.-1.d-3) then
            do j=1,17
              ! NOTE: The FT distribution misfit calculation does not use the
              ! misfit types given in the topo_parameters file.
              misfit=misfit+(tlobs(j)-ftdisto(j,iobs))**2
              !misfit=misfit+age_misfit(tlobs(j),real(ftdisto(j,iobs)),1.d0,misfit_type)
              nmisfit=nmisfit+1
            enddo
          endif
          if (Tramanobs.gt.0) then
            !misfit=misfit+(Tramanobs-Traman(iobs))**2/dTramanobs**2
            misfit=misfit+age_misfit(Tramanobs,Traman(iobs),dTramanobs,misfit_type)
            nmisfit=nmisfit+1
          endif
        enddo

        if (nmisfit.ne.0) then
          if (nd.ne.0) ndmisfit=nd
          select case(misfit_type)
          case(1)
            misfit=sqrt(misfit)/nmisfit
          case(2)
            misfit=misfit/(nmisfit-ndmisfit-1)
          end select
        endif

        !if (nmisfit.ne.0) misfit=sqrt(misfit)/nmisfit

        ! set misfit_slope to .TRUE. to use the slopes for the misfit
        misfit_slope=.FALSE.
        if (mftflag.eq.1) misfit_slope=.TRUE.
        if (misfit_slope) then
          if (nhe.gt.1) call regression (aheoreg,hheoreg,nhe,aheom,hheom,heoi,heos)
          if (nft.gt.1) call regression (aftoreg,hftoreg,nft,aftom,hftom,ftoi,ftos)
          if (nheZ.gt.1) call regression (aheZoreg,hheZoreg,nheZ,aheZom,hheZom,heZoi,heZos)
          if (nftZ.gt.1) call regression (aftZoreg,hftZoreg,nftZ,aftZom,hftZom,ftZoi,ftZos)
          if (narK.gt.1) call regression (aarKoreg,harKoreg,narK,aarKom,harKom,arKoi,arKos)
          if (narB.gt.1) call regression (aarBoreg,harBoreg,narB,aarBom,harBom,arBoi,arBos)
          if (narM.gt.1) call regression (aarMoreg,harMoreg,narM,aarMom,harMom,arMoi,arMos)
          if (narH.gt.1) call regression (aarHoreg,harHoreg,narH,aarHom,harHom,arHoi,arHos)
          if (nobs.gt.1) call regression (ahepreg,hpreg,nobs,ahepm,hhepm,hepi,heps)
          if (nobs.gt.1) call regression (aftpreg,hpreg,nobs,aftpm,hftpm,ftpi,ftps)
          if (nobs.gt.1) call regression (aheZpreg,hpreg,nobs,aheZpm,hheZpm,heZpi,heZps)
          if (nobs.gt.1) call regression (aftZpreg,hpreg,nobs,aftZpm,hftZpm,ftZpi,ftZps)
          if (nobs.gt.1) call regression (aarKpreg,hpreg,nobs,aarKpm,harKpm,arKpi,arKps)
          if (nobs.gt.1) call regression (aarBpreg,hpreg,nobs,aarBpm,harBpm,arBpi,arBps)
          if (nobs.gt.1) call regression (aarMpreg,hpreg,nobs,aarMpm,harMpm,arMpi,arMps)
          if (nobs.gt.1) call regression (aarHpreg,hpreg,nobs,aarHpm,harHpm,arHpi,arHps)
          misfit=0.d0
          if (nhe.gt.1) misfit=misfit+abs(heos-heps)/abs(heos)+abs(aheom-ahepm)/abs(aheom)
          if (nft.gt.1) misfit=misfit+abs(ftos-ftps)/abs(ftos)+abs(aftom-aftpm)/abs(aftom)
          if (nheZ.gt.1) misfit=misfit+abs(heZos-heZps)/abs(heZos)+abs(aheZom-aheZpm)/abs(aheZom)
          if (nftZ.gt.1) misfit=misfit+abs(ftZos-ftZps)/abs(ftZos)+abs(aftZom-aftZpm)/abs(aftZom)
          if (narK.gt.1) misfit=misfit+abs(arKos-arKps)/abs(arKos)+abs(aarKom-aarKpm)/abs(aarKom)
          if (narB.gt.1) misfit=misfit+abs(arBos-arBps)/abs(arBos)+abs(aarBom-aarBpm)/abs(aarBom)
          if (narM.gt.1) misfit=misfit+abs(arMos-arMps)/abs(arMos)+abs(aarMom-aarMpm)/abs(aarMom)
          if (narH.gt.1) misfit=misfit+abs(arHos-arHps)/abs(arHos)+abs(aarHom-aarHpm)/abs(aarHom)
          if (nd.eq.0) then
            if (nhe.gt.0) print*,'HeAp s o/p',nhe,heos,heps
            if (nft.gt.0) print*,'FTAp s o/p',nft,ftos,ftps
            if (nheZ.gt.0) print*,'HeZr s o/p',nheZ,heZos,heZps
            if (nftZ.gt.0) print*,'HeZr s o/p',nftZ,ftZos,ftZps
            if (narK.gt.0) print*,'ArKF s o/p',narK,arKos,arKps
            if (narB.gt.0) print*,'HeBi s o/p',narB,arBos,arBps
            if (narM.gt.0) print*,'HeMu s o/p',narM,arMos,arMps
            if (narH.gt.0) print*,'HeHo s o/p',narH,arHos,arHps
            if (nhe.gt.0) print*,'HeAp m o/p',nhe,aheom,ahepm
            if (nft.gt.0) print*,'FTAp m o/p',nft,aftom,aftpm
            if (nheZ.gt.0) print*,'HeZr m o/p',nheZ,aheZom,aheZpm
            if (nftZ.gt.0) print*,'HeZr m o/p',nftZ,aftZom,aftZpm
            if (narK.gt.0) print*,'HeKF m o/p',narK,aarKom,aarKpm
            if (narB.gt.0) print*,'HeBi m o/p',narB,aarBom,aarBpm
            if (narM.gt.0) print*,'HeMu m o/p',narM,aarMom,aarMpm
            if (narH.gt.0) print*,'HeHo m o/p',narH,aarHom,aarHpm
          endif
        endif
        if (nd.ne.0) then
          open (71,file='NA/NA_int_res.txt',status='unknown',position='append')
          write (71,'(100e14.5)') misfit,(param(i),i=1,nd)
          close (71)
        endif
        if (nd.eq.0) write (6,'(a,f12.4,a,i3,a,i1,a)') ' Misfit : ',misfit,' (',nmisfit,' ages, misfit type ',misfit_type,')'
        close (13)

        ! added by Jean to change misfit to slope rather than exact ages
        ! (4/6/2008)
        deallocate (aheoreg,ahepreg,hheoreg)
        deallocate (aftoreg,aftpreg,hftoreg)
        deallocate (aheZoreg,aheZpreg,hheZoreg)
        deallocate (aftZoreg,aftZpreg,hftZoreg)
        deallocate (aarKoreg,aarKpreg,harKoreg)
        deallocate (aarBoreg,aarBpreg,harBoreg)
        deallocate (aarMoreg,aarMpreg,harMoreg)
        deallocate (aarHoreg,aarHpreg,harHoreg)
        deallocate (aramanoreg,aramanpreg,hramanoreg)
        deallocate (hei,agehe,ageheZ,ageft,ageftZ,agearB,agearK,agearM,agearH,ftdisto,Traman)
        deallocate (hpreg)
      endif

      if (calc_surface_ages) deallocate (age1,age2,age3,age4,age5,age6,age7,age8,raman)
      !deallocate (age1,age2,age3,age4,age5,age6,age7,age8,raman)

      deallocate (iconsurf,ielsurf,neighbour)
      deallocate (xsurf,ysurf,zsurf,zsurfp)
      deallocate (topoa,topob,rsurf,rebound)
      deallocate (xdepth,ydepth,zdepth,tprev)
      deallocate (xexhumation,yexhumation,zexhumation)

      deallocate (x,y,z,t)
      deallocate (xp,yp,zp,tp)
      deallocate (icon,ael,bel)
      deallocate (kfix)
      deallocate (f)
      deallocate (proc,eproc,ice)

      call cpu_time (times(9))

      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Total times : ',times(9)-times(1)

      close (7)

      return

      contains
        function age_misfit(obs_age,pred_age,sigma_obs_age,misfit_type)
          implicit none
        
          double precision, intent(in) :: obs_age,sigma_obs_age
          integer, intent(in) :: misfit_type
          real*4 :: pred_age
          real*4 :: age_misfit
        
          select case(misfit_type)
          case(1)
            age_misfit=(pred_age-obs_age)**2/sigma_obs_age**2
          case(2)
            age_misfit=(pred_age-obs_age)**2/sigma_obs_age**2
          case(3)
            age_misfit=((pred_age-obs_age)/sigma_obs_age)**2
          case default
            write (*,*) 'Error: Unsupported value for misfit_type.'
            write (*,*) 'misfit_type must be an integer between 1 and 3. Stopping.'
            stop
          end select
        end function age_misfit

      end

!-----------
      subroutine regression (x,y,n,xmean,ymean,intercept,slope)

! this routine is used to calculate the means, intercept and slope
! of the regression between two datasets of length n stored in x and y

      implicit none
      
      integer n
      double precision x(n),y(n),xmean,ymean,intercept,slope

      xmean=sum(x)/n
      ymean=sum(y)/n
      slope=sum((x-xmean)*(y-ymean))/sum((x-xmean)**2)
      intercept=ymean-xmean*slope

      return
      end
