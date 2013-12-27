      subroutine create_pecube_in (fault,nfault,nd,range,param)

      use definitions

      implicit none

! Subroutine to create the input file Pecube.in to be used by Pecube

! The user can modify this file at will but the output it produces (Pecube.in)
! must obey rules that are described in the user guide and in Pecube.f90

      integer i1,j1,nfault,ieobs,ilog,icol,jd,jcol

      type (faulttype) fault(nfault)

      double precision,dimension(:,:),allocatable::zNZ
      double precision,dimension(:),allocatable::z,xz,yz
      double precision,dimension(:),allocatable::timek,topomag,topooffset
      integer,dimension(:),allocatable::iout
      integer,dimension(:,:),allocatable::iconz,neighz
      integer i,j,k,nfnme,nx0,ny0,nskip,nstep,istep,isoflag,nxiso,nyiso,nzin,ij,it
      integer nx,ny,nsurf,nelem,nne,nobsfile,iterative,interpol,icon1,icon2,icon3,icon4,nobs
      integer ftlflag,mftflag,fltflag,ageflag(10),calc_surface_ages_in,advect_topo
      integer misfit_type,ndmisfit,PTt_output_in
      double precision dx,dy,ddx,ddy,xlon,xlat,tau,rhoc,rhom,young,poisson,thickness,tprevious
      double precision crustal_thickness,diffusivity,tmax,tmsl,tlapse,heatproduction,friction
      double precision xl,yl,zl,xlon1,xlat1,xlon2,xlat2,x,y,r,s
      double precision xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs
      double precision wobs1,wobs2,wobs3,wobs4,a1,b1,c1,a2,b2,c2,a3,b3,c3,surf2
      double precision ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs
      double precision ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs
      double precision ftdist(17),Tramanobs,dTramanobs
      double precision vtopo

      character run*5,fnme*300,obsfile*300,line*1024,c5*5

      real*4 range(2,*),param(*)
      integer nd,nd0,nc5
      logical vivi !VKP
      logical calc_surface_ages,PTt_output

      nd0=nd
      nd=0

      open (54,file='input/topo_parameters.txt',status='old')
      open (55,status='scratch')
    1 read (54,'(a1024)',end=2) line
        if (line(1:1).ne.'$'.and. line(1:1).ne.' ') then
          if (scan(line,'$').ne.0) then
            do i=scan(line,'$'),1024
            line(i:i)=' '
            enddo
          endif
        k=1
          do j=1,1024
            if (line(j:j).eq.' '.or.line(j:j).eq.',') then
            if (j.ne.k) write (55,'(a)') line(k:j-1)
            k=j+1
            endif
          enddo
        endif
      goto 1
    2 close (54)
      rewind (55)

! run is the name of the run (assumes a directory of that name exists)
! should be 5 character long
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'run cannot be specified as a range'
      backspace (55)
      read (55,'(a5)') run

! fnme is the name of the input topographic file 
! nfnme is the length of the file name

        do i=1,300
        fnme(i:i)=' '
        enddo
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'fnme cannot be specified as a range'
      backspace (55)
      read (55,'(a)') fnme
        do i=1,300
        if (fnme(i:i).ne.' ') nfnme=i
        enddo
      vivi=.FALSE. !VKP
      if (fnme(nfnme:nfnme).eq.'/') vivi=.TRUE. !VKP
      if (vivi) nfnme=nfnme-1 !VKP

! nx0 and ny0 are the number of points in the input topographic file
! dx and dy are the data spacings in the x- and y-direction in the input topographic
! file (in degrees)
! nskip is the number of points that are skipped when sampling the initial
! topo file
! xlat and xlon are the latitude/longitude of the bottom left corner of the
! data set (in deg.)
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nx0 cannot be specified as a range'
      backspace (55)
      read (55,*) nx0
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'ny0 cannot be specified as a range'
      backspace (55)
      read (55,*) ny0
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'dx cannot be specified as a range'
      backspace (55)
      read (55,*) dx
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'dy cannot be specified as a range'
      backspace (55)
      read (55,*) dy
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nskip cannot be specified as a range'
      backspace (55)
      read (55,*) nskip
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'xlon cannot be specified as a range'
      backspace (55)
      read (55,*) xlon
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'xlat cannot be specified as a range'
      backspace (55)
      read (55,*) xlat

! nstep is the number of time step
! tau is the erosion time scale (assuming an exponential decvrease in topography) (in Myr)
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nstep cannot be specified as a range'
      backspace (55)
      read (55,*) nstep
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      tau=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      tau=param(jd)
      else
      read (line,*) tau
      endif

! for each time step + 1
! timek is the time (in Myr)
! topomag is the topographic amplification factor at this step
! topooffset is the topographic offset
      allocate (timek(nstep+1),topomag(nstep+1),topooffset(nstep+1),iout(nstep+1))

        do istep=1,nstep+1
        read (55,'(a1024)') line
        icol=scan(line,':')
        jcol=scan(line,'#')
        if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        timek(istep)=param(nd)
        elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        timek(istep)=param(jd)
        else
        read (line,*) timek(istep)
        endif
        read (55,'(a1024)') line
        icol=scan(line,':')
        jcol=scan(line,'#')
        if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        topomag(istep)=param(nd)
        elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        topomag(istep)=param(jd)
        else
        read (line,*) topomag(istep)
        endif
        read (55,'(a1024)') line
        icol=scan(line,':')
        jcol=scan(line,'#')
        if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        topooffset(istep)=param(nd)
        elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        topooffset(istep)=param(jd)
        else
        read (line,*) topooffset(istep)
        endif
        read (55,'(a1024)') line
        if (scan(line,':').ne.0) stop 'iout cannot be specified as a range'
        backspace (55)
        read (55,*) iout(istep)
        enddo

! converts geological time into model time
        do istep=nstep+1,1,-1
        timek(istep)=timek(1)-timek(istep)
        enddo

! isostasy flag (0 no isostasy, 1 isostasy on)
! rhoc and rhom are the densities for the crust and mantle, respectively (in kg/m3)
! these values are used in the isostatic calculations
! young is the elastic plate young modulus (in Pa)
! poisson is poisson's ratio (dimensionless)
! thickness is the elastic thickness of the plate (in km)
! nxiso and nyiso are the resolutions in the x- and y-directions of the grid on
! which the isostatic (flexural) calculations are performed (including the FFT)
! note that these numbers must be powers of two.
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'isoflag cannot be specified as a range'
      backspace (55)
      read (55,*) isoflag
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      rhoc=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      rhoc=param(jd)
      else
      read (line,*) rhoc
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      rhom=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      rhom=param(jd)
      else
      read (line,*) rhom
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      young=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      young=param(jd)
      else
      read (line,*) young
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      poisson=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      poisson=param(jd)
      else
      read (line,*) poisson
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      thickness=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      thickness=param(jd)
      else
      read (line,*) thickness
      endif
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nxiso cannot be specified as a range'
      backspace (55)
      read (55,*) nxiso
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nyiso cannot be specified as a range'
      backspace (55)
      read (55,*) nyiso

! crustal thickness is the averaged crustal thickness (i.e. the depth at which the
! temperature is assumed to be constant) (in km)
! nzin is the number of points in the z-direction
!diffusivity is the heat diffusivity (in km2/Myr)
! tmax is the basal temperature (in C)
! tmsl is the temperature at the top of the model (at z=0)
! tlapse is the lapse rate (or change of temperature with height in the atmosphere)
! (in C/km)
! heatproduction is the rate of heat production (in C/Myr)
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      crustal_thickness=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      crustal_thickness=param(jd)
      else
      read (line,*) crustal_thickness
      endif
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'nz cannot be specified as a range'
      backspace (55)
      read (55,*) nzin
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      diffusivity=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      diffusivity=param(jd)
      else
      read (line,*) diffusivity
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      tmax=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      tmax=param(jd)
      else
      read (line,*) tmax
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      tmsl=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      tmsl=param(jd)
      else
      read (line,*) tmsl
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      tlapse=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      tlapse=param(jd)
      else
      read (line,*) tlapse
      endif
      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      heatproduction=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      heatproduction=param(jd)
      else
      read (line,*) heatproduction
      endif

      crustal_thickness=crustal_thickness*1.d3

! obsfile is the name of the observation file
        do i=1,300
        obsfile(i:i)=' '
        enddo
      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'obsfile cannot be specified as a range'
      backspace (55)
      read (55,'(a)') obsfile
        do i=1,300
        if (obsfile(i:i).ne.' ') nobsfile=i
        enddo

! new parameters added in Pecube2

      tprevious=timek(nstep+1)
      ftlflag=0
      mftflag=0
      fltflag=0
      friction=0.d0
      ageflag=1

      read (55,'(a1024)',end=999) line
      if (scan(line,':').ne.0) stop 'default age cannot be specified as a range'
      backspace (55)
      read (55,*) tprevious
      if (tprevious.eq.0.) tprevious=timek(nstep+1)

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'FT model flag cannot be specified as a range'
      backspace (55)
      read (55,*) ftlflag

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'misfit flag cannot be specified as a range'
      backspace (55)
      read (55,*) mftflag

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'default age cannot be specified as a range'
      backspace (55)
      read (55,*) fltflag

      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      friction=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      friction=param(jd)
      else
      read (line,*) friction
      endif

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'topographic advection flag cannot be specified as a range'
      backspace (55)
      read (55,*) advect_topo
      if (advect_topo /= 0 .and. advect_topo /= 1 .and. advect_topo /= -1) then
        write (*,*) 'Bad value, ',advect_topo,' for topographic advection flag. Input value must be'
        write (*,*) '"0", "1" or "-1". Exiting.'
        stop
      endif

      read (55,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
      nd=nd+1
      read (line(1:icol-1),*) range(1,nd)
      read (line(icol+1:1024),*) range(2,nd)
      vtopo=param(nd)
      elseif (jcol.ne.0) then
      read (line(jcol+1:1024),*) jd
      vtopo=param(jd)
      else
      read (line,*) vtopo
      endif

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'Misfit type cannot be specified as a range'
      backspace (55)
      read (55,*) misfit_type

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'Number of free parameters cannot be specified as a range'
      backspace (55)
      read (55,*) ndmisfit

      do i=1,10
        read (55,'(a1024)') line
        if (scan(line,':').ne.0) stop 'age flags cannot be specified as a range'
        backspace (55)
        read (55,*) ageflag(i)
      enddo

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'Surface age calculation flag cannot be specified as a range'
      backspace (55)
      read (55,*) calc_surface_ages_in

      if (calc_surface_ages_in == 0) then
        calc_surface_ages=.false.
      elseif (calc_surface_ages_in == 1) then
        calc_surface_ages=.true.
      else
        write (*,*) 'Bad value for calculate surface ages flag. Input value must be'
        write (*,*) '"0" or "1". Exiting.'
        stop
      endif

      read (55,'(a1024)') line
      if (scan(line,':').ne.0) stop 'PTt output flag cannot be specified as a range'
      backspace (55)
      read (55,*) PTt_output_in

      if (PTt_output_in == 0) then
        PTt_output=.false.
      elseif (PTt_output_in == 1) then
        PTt_output=.true.
      else
        write (*,*) 'Bad value for PTt output flag. Input value must be'
        write (*,*) '"0" or "1". Exiting.'
        stop
      endif

 999  close (55)

! reads in topography

if (.not.vivi) then !VKP

      if (nx0.gt.0) then

      allocate (zNZ(nx0,ny0))
        if (fnme(1:nfnme).eq.'Nil') then
        zNZ=0.d0
        else
        open (8,file='data/'//fnme(1:nfnme),status='old')
        read (8,*) zNZ
        close (8)
        endif

      nx=(nx0-1)/nskip+1
      ny=(ny0-1)/nskip+1
      allocate (z(nx*ny))
      ij=0
        do j=1,ny0,nskip
          do i=1,nx0,nskip
          ij=ij+1
          z(ij)=zNZ(i,j)
          enddo
        enddo

      deallocate (zNZ)
	  
      xlon1=xlon
      xlat1=xlat
      xlon2=xlon+(nx-1)*dx*nskip
      xlat2=xlat+(ny-1)*dy*nskip

      else

      allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0))
        open (8,file='data/'//fnme(1:nfnme),status='old')
          do i=1,-nx0
          read (8,*) z(i)
          enddo
        close (8)
        open (8,file='data/'//fnme(1:nfnme)//'.geometry',status='old')
          do i=1,-nx0
          read (8,*) xz(i),yz(i)
          enddo
          do i=1,-ny0
          read (8,*) (iconz(k,i),k=1,3)
          enddo
        close (8)

      xlon1=minval(xz)
      xlat1=minval(yz)
      xlon2=maxval(xz)
      xlat2=maxval(yz)
      xlon=xlon1
      xlat=xlat1

      endif

else !VKP

      if (nx0.gt.0) then

      nx=(nx0-1)/nskip+1 !VKP
      ny=(ny0-1)/nskip+1 !VKP
      allocate (z(nx*ny)) !VKP
      topomag=1.d0 !VKP
      topooffset=0.d0 !VKP

      else

      allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0))
        open (8,file='data/'//fnme(1:nfnme)//'/geometry',status='old')
          do i=1,-nx0
          read (8,*) xz(i),yz(i)
          enddo
          do i=1,-ny0
          read (8,*) (iconz(k,i),k=1,3)
          enddo
        close (8)
      topomag=1.d0 !VKP
      topooffset=0.d0 !VKP
      xlon1=minval(xz)
      xlat1=minval(yz)
      xlon2=maxval(xz)
      xlat2=maxval(yz)
      xlon=xlon1
      xlat=xlat1

      endif

endif !VKP

      if (nx0.gt.0) then

      xl=dx*(nx-1)*nskip*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
      yl=dy*(ny-1)*nskip*111.11
      zl=crustal_thickness/1.e3
      z=z/crustal_thickness*zl

      else

      zl=crustal_thickness/1.e3
      z=z/crustal_thickness*zl
      xl=(xlon2-xlon1)*111.11*cos((xlat+(xlat2-xlat1)/2.)*3.141592654/180.)
      yl=(xlat2-xlat1)*111.11

      endif

! reads in fault definitions

      call read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,yl,zl,timek(nstep+1), &
                                     nd,range,param)

      if (nd0.gt.0) nd=nd0

      open (7,status='scratch')

      ilog=0
      iterative=1
      interpol=1
      nsurf=nx*ny
      nelem=(nx-1)*(ny-1)
      nne=4
      ddx=xl/(nx-1)*1.d3
      ddy=yl/(ny-1)*1.d3
        if (nx0.lt.0) then
        nsurf=-nx0
        nelem=-ny0
        nne=3
        ddx=xl/sqrt(dble(nsurf))*1.d3
        ddy=yl/sqrt(dble(nsurf))*1.d3
	endif

      write (7,'(a)') run
        if (vivi) then !VKP
        write (7,*) nne,-nsurf,nzin,nelem,zl,diffusivity,heatproduction,friction !VKP
        else !VKP
        write (7,*) nne,nsurf,nzin,nelem,zl,diffusivity,heatproduction,friction
        endif !VKP
      write (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      write (7,*) isoflag,tau,rhoc,rhom
      write (7,*) nx,ny,nxiso,nyiso
      write (7,*) ddx,ddy,young,poisson,thickness*1.d3
      write (7,*) xlon1,xlon2,xlat1,xlat2
      write (7,*) vtopo,advect_topo
      write (7,*) misfit_type,ndmisfit

      if (nx0.gt.0) then

        do j=1,ny
          do i=1,nx
          x=xl*float(i-1)/float(nx-1)
          y=yl*float(j-1)/float(ny-1)
          write (7,*) x,y
          enddo
        enddo

        do j=1,ny-1
          do i=1,nx-1
          icon1=(j-1)*nx+i
          icon2=icon1+1
          icon3=icon1+nx+1
          icon4=icon1+nx
          write (7,*) icon1,icon2,icon3,icon4
          enddo
        enddo

      else

        do i=1,nsurf
        x=(xz(i)-xlon1)*111.11*cos((xlat+(xlat2-xlat1)/2.)*3.141592654/180.)
        y=(yz(i)-xlat1)*111.11
        write (7,*) x,y
        enddo

        do i=1,nelem
        write (7,*) (iconz(k,i),k=1,3)
        enddo

      endif

        do istep=0,nstep
        write (7,*) timek(istep+1),iout(istep+1)
          if (vivi) then !VKP
            if (istep.lt.10) then !VKP
            write (c5(1:1),'(i1)') istep !VKP
            nc5=1 !VKP
            elseif (istep.lt.100) then !VKP
            write (c5(1:2),'(i2)') istep !VKP
            nc5=2 !VKP
            elseif (istep.lt.1000) then !VKP
            write (c5(1:3),'(i3)') istep !VKP
            nc5=3 !VKP
            elseif (istep.lt.10000) then !VKP
            write (c5(1:4),'(i4)') istep !VKP
            nc5=4 !VKP
            else !VKP
            write (c5(1:5),'(i5)') istep !VKP
            nc5=5 !VKP
            endif !VKP
            if (nx0.gt.0) then
            allocate (zNZ(nx0,ny0))
            open (67,file='data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip  !VKP
                do i=1,nx0,nskip  !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            z=z/1.e3 !VKP
            close (67) !VKP
            else
            open (67,file='data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            z=z/1.e3
            close (67)
            endif
          endif !VKP
        write (7,*) (z(k)*topomag(istep+1)+topooffset(istep+1),k=1,nsurf)
          if (vivi) then !VKP
            if (nx0.gt.0) then
            open (67,file='data/'//fnme(1:nfnme)//'/uplift'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip !VKP
                do i=1,nx0,nskip !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            close (67) !VKP
            write (7,*) (z(k),k=1,nsurf) !VKP
            open (67,file='data/'//fnme(1:nfnme)//'/temp'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip !VKP
                do i=1,nx0,nskip !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            close (67) !VKP
            write (7,*) (z(k),k=1,nsurf) !VKP
            deallocate (zNZ)
            else
            open (67,file='data/'//fnme(1:nfnme)//'/uplift'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            close (67)
            write (7,*) (z(k),k=1,nsurf) !VKP
            open (67,file='data/'//fnme(1:nfnme)//'/temp'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            close (67)
            write (7,*) (z(k),k=1,nsurf)
            endif
          endif !VKP
        enddo

! observations

      if (obsfile(1:nobsfile).eq.'Nil') then
      nobs=0
      write (7,*) nobs
      else
        if (nx0.lt.0) then
        allocate (neighz(3,nelem))
        call neighbours (iconz,neighz,nelem)
        it=1
        endif
      open (8,file='data/'//obsfile(1:nobsfile),status='old')
      read (8,*) nobs
      write (7,*) nobs
        do i=1,nobs
!        read (8,*) xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs
        read (8,*) xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                   ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                   ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                   ftdist,Tramanobs,dTramanobs
        if (nx0.gt.0) then
        i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
        if (i1.eq.nx) i1=nx-1
        j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
        if (j1.eq.ny) j1=ny-1
        ieobs=i1+(j1-1)*(nx-1)
        r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
        r=-1.+2.*r
        s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
        s=-1.+2.*s
        wobs1=(1.-r)*(1.-s)/4.
        wobs2=(1.+r)*(1.-s)/4.
        wobs3=(1.+r)*(1.+s)/4.
        wobs4=(1.-r)*(1.+s)/4.
        else
        call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
        ieobs=it
        surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
             -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
        a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
        b1=yz(iconz(2,it))-yz(iconz(3,it))
        c1=xz(iconz(3,it))-xz(iconz(2,it))
        a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
        b2=yz(iconz(3,it))-yz(iconz(1,it))
        c2=xz(iconz(1,it))-xz(iconz(3,it))
        a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
        b3=yz(iconz(1,it))-yz(iconz(2,it))
        c3=xz(iconz(2,it))-xz(iconz(1,it))
        wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
        wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
        wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
        wobs4=0.
        endif
        write (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                    ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                    ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                    ftdist,Tramanobs,dTramanobs,heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
      close (8)
      endif

      write (7,*) ageflag,calc_surface_ages,PTt_output

      deallocate (z)
      if (nx0.lt.0) deallocate (xz,yz,iconz)
      deallocate (timek,topomag,topooffset)

      return
      end