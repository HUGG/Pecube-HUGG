subroutine testTopo (nx,ny,nz,xl,yl,zl,timesteps,nstep,xlon1,xlat1,xlon2,xlat2,&
                     fltflag,advect_topo,vtopo)

implicit none

character line*1024,run*5,fnme*300,obsfile*300,cs*3
integer nx0,ny0,nskip,nstep,istep,isoflag,nxiso,nyiso,nz,i,j,k,ii,jj,nx,ny,nobs,ij
integer nfnme,nobsfile,iunit,icon1,icon2,icon3,icon4,itime,fltflag,mftflag,ftlflag
double precision dx,dy,xlon,xlat,tau,rhoc,rhom,young,poisson,thickness,dum,tprevious
double precision crustal_thickness,diffusivity,tmax,tmsl,tlapse,heatproduction
double precision zl,xl,yl,x,y,h
double precision,dimension(:),allocatable::timek,topomag,topooffset,a1,a2,a3,a4,a5,a6,a7,a8
double precision,dimension(:,:),allocatable::zNZ
double precision,dimension(:),allocatable::z,xz,yz
double precision timesteps(1000)
double precision xlon1,xlat1,xlon2,xlat2,vtopo,friction
integer,dimension(:,:),allocatable::iconz
integer,dimension(:),allocatable :: iout

logical vivi
character c5*5
integer nc5,advect_topo

open (77,status='scratch')

  open (76,file='input/topo_parameters.txt',status='old')
1 read (76,'(a1024)',end=2) line
  if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (77,'(a)') line
  goto 1
2 close (76)

rewind (77)

read (77,'(a5)') run

  do i=1,300
  fnme(i:i)=' '
  enddo
read (77,'(a)') fnme
  do i=1,300
  if (fnme(i:i).ne.' ') nfnme=i
  enddo
vivi=.FALSE.
if (fnme(nfnme:nfnme).eq.'/') vivi=.TRUE.
if (vivi) nfnme=nfnme-1

read (77,*) nx0,ny0,dx,dy,nskip,xlon,xlat

read (77,*) nstep,tau

allocate (timek(nstep+1),topomag(nstep+1),topooffset(nstep+1),iout(nstep+1))

do istep=1,nstep+1
  read (77,*) timek(istep),topomag(istep),topooffset(istep),iout(istep)
enddo

! converts geological time into model time
  do istep=nstep+1,1,-1
  timek(istep)=timek(1)-timek(istep)
  timesteps(istep)=timek(istep)
  enddo

read (77,*) isoflag,rhoc,rhom,young,poisson,thickness,nxiso,nyiso

read (77,*) crustal_thickness,nz,diffusivity,tmax,tmsl,tlapse,heatproduction

crustal_thickness=crustal_thickness*1.d3

  do i=1,300
  obsfile(i:i)=' '
  enddo
read (77,'(a)') obsfile
  do i=1,300
  if (obsfile(i:i).ne.' ') nobsfile=i
  enddo
print*,obsfile(1:nobsfile)

read (77,*,end=999) tprevious,ftlflag,mftflag,fltflag,friction,advect_topo,vtopo

!read (77,*,end=999) tprevious
!read (77,*,end=999) ftlflag
!read (77,*,end=999) mftflag
!read (77,*,end=999) fltflag

999 close (77)

if (advect_topo /= 0 .and. advect_topo /= 1 .and. advect_topo /= -1) then
  write (*,*) 'Bad value, ',advect_topo,' for topographic advection flag. Input value must be'
  write (*,*) '"0", "1" or "-1". Exiting.'
  stop
endif

if (.not.vivi) then

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
  xlon2=xlon+(nx0-1)*dx
  xlat2=xlat+(ny0-1)*dy
  
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

else

  if (nx0.gt.0) then

  nx=(nx0-1)/nskip+1 !VKP
  ny=(ny0-1)/nskip+1 !VKP
  allocate (zNZ(nx0,ny0)) !VKP
  allocate (z(nx*ny)) !VKP
  topomag=1.d0 !VKP
  topooffset=0.d0 !VKP

  ! Calculate latitude/longitude ranges - dwhipp 11.14
  xlon1=xlon
  xlat1=xlat
  xlon2=xlon+(nx0-1)*dx
  xlat2=xlat+(ny0-1)*dy

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

endif

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

iunit=30
do istep=0,nstep
write(cs,'(i3)') istep
if (istep.lt.10) cs(1:2)='00'
if (istep.lt.100) cs(1:1)='0'
itime=int(timesteps(istep+1))

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
    open (67,file='data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
    read (67,*) zNZ !VKP
    ij=0
      do j=1,ny0,nskip  !VKP
        do i=1,nx0,nskip  !VKP
        ij=ij+1
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

open(unit=iunit,file='VTK/Topo'//cs//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'surface'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
if (nx0.gt.0) then
write(iunit,'(a7,i10,a6)')'POINTS ',nx*ny,' float'
ij=0
  do j=1,ny
  do i=1,nx
  x=xl*float(i-1)/float(nx-1)
  y=yl*float(j-1)/float(ny-1)
  ij=ij+1
  write(iunit,'(3f16.11)') x,y,z(ij)*topomag(istep+1)+topooffset(istep+1)+zl
  enddo
  enddo
write(iunit,'(A6, 2I10)') 'CELLS ',(nx-1)*(ny-1),(1+4)*(nx-1)*(ny-1)
  do j=1,ny-1
    do i=1,nx-1
    icon1=(j-1)*nx+i-1
    icon2=icon1+1
    icon3=icon1+nx+1
    icon4=icon1+nx
    write (iunit,'(9I10)') 4,icon1,icon2,icon3,icon4
    enddo
  enddo
write(iunit,'(A11, I10)') 'CELL_TYPES ',(nx-1)*(ny-1)
  do k=1,(nx-1)*(ny-1)
  write(iunit,'(I2)')9 ! rectangles
  enddo
write(iunit,'(a11,i10)')'POINT_DATA ',nx*ny
write(iunit,'(a)')'SCALARS Topo float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
ij=0
  do j=1,ny
    do i=1,nx
    ij=ij+1
    write(iunit,'(f18.13)') z(ij)*topomag(istep+1)+topooffset(istep+1)
    enddo
  enddo
else
write(iunit,'(a7,i10,a6)')'POINTS ',-nx0,' float'
  do i=1,-nx0
  write(iunit,'(3f16.11)') (xz(i)-xlon1)/(xlon2-xlon1)*xl, &
                           (yz(i)-xlat1)/(xlat2-xlat1)*yl, &
                           z(i)*topomag(istep+1)+topooffset(istep+1)+zl
  enddo
write(iunit,'(A6, 2I10)') 'CELLS ',-ny0,(1+3)*(-ny0)
    do i=1,-ny0
    write (iunit,'(9I10)') 3,iconz(:,i)-1
    enddo
write(iunit,'(A11, I10)') 'CELL_TYPES ',-ny0
  do k=1,-ny0
  write(iunit,'(I2)')5 ! rectangles
  enddo
write(iunit,'(a11,i10)')'POINT_DATA ',-nx0
write(iunit,'(a)')'SCALARS Topo float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,-nx0
    write(iunit,'(f18.13)') z(i)*topomag(istep+1)+topooffset(istep+1)
    enddo
endif
close(iunit)
print*,'Step',istep,'done'
enddo

if (obsfile(1:nobsfile).eq.'Nil') return

open (88,file='data/'//obsfile(1:nobsfile),status='old')
read (88,*) nobs
allocate (a1(nobs),a2(nobs),a3(nobs),a4(nobs),a5(nobs),a6(nobs),a7(nobs),a8(nobs))

open(unit=iunit,file='VTK/Data.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'AgeData'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nobs,' float'
  do i=1,nobs
  read (88,*) x,y,h,a1(i),dum,a2(i),dum,a3(i),dum,a4(i),dum,a5(i),dum,a6(i),dum,a7(i),dum,a8(i),dum
  x=(x-xlon1)/(xlon2-xlon1)*xl
  y=(y-xlat1)/(xlat2-xlat1)*yl
  write(iunit,'(3f16.11)') x,y,h/1.e3+zl
  enddo
a1=max(a1,-1.d0)
a2=max(a2,-1.d0)
a3=max(a3,-1.d0)
a4=max(a4,-1.d0)
a5=max(a5,-1.d0)
a6=max(a6,-1.d0)
a7=max(a7,-1.d0)
a8=max(a8,-1.d0)

write(iunit,'(A6, 2I10)') 'CELLS ',1,nobs+1
write (iunit,'(256I10)') nobs,(i-1,i=1,nobs)
write(iunit,'(A11, I10)') 'CELL_TYPES ',1
write (iunit,'(I10)') 2
write(iunit,'(a11,i10)')'POINT_DATA ',nobs
write(iunit,'(a)')'SCALARS ApatiteHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a1(i)
  enddo
write(iunit,'(a)')'SCALARS ApatiteFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a2(i)
  enddo
write(iunit,'(a)')'SCALARS ZirconHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a3(i)
  enddo
write(iunit,'(a)')'SCALARS ZirconFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a4(i)
  enddo
write(iunit,'(a)')'SCALARS KSparArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a5(i)
  enddo
write(iunit,'(a)')'SCALARS BiotiteArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a6(i)
  enddo
write(iunit,'(a)')'SCALARS MuscoviteArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a7(i)
  enddo
write(iunit,'(a)')'SCALARS HornblendeArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a8(i)
  enddo
close(iunit)

end