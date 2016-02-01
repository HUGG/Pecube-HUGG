!----------------------------

subroutine read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,  &
                                     yl,zl,timeend,nd,range,param)

! this subroutines reads the fault information from the fault_parameters.txt file
! and generates the necessary complementary information to be used for fault
! kinematics and movement

use definitions

implicit none

integer nfault,i,j,k,icol,jcol,kcol,jd,faultmodel,advect_duplex
type(faulttype) fault(nfault)
double precision xn,yn,xyn,x1,y1,x2,y2,timeend
double precision xlon1,xlat1,xlon2,xlat2,xl,yl,zl
double precision dstart,dend,dvelo,dinner,douter,eps,d1,vduplexadv
character line*1024
logical partition

real*4 range(2,*),param(*)
integer nd

if (nfault.eq.0) return

      open (76,file='input/fault_parameters.txt',status='old')
      open (77,status='scratch')
    1 read (76,'(a1024)',end=2) line
      if (line(1:1).ne.'$'.and. line(1:1).ne.' ') then
        if (scan(line,'$').ne.0) then
          do i=scan(line,'$'),1024
            line(i:i)=' '
          enddo
        endif
        k=1
        do j=1,1024
          if (line(j:j).eq.' '.or.line(j:j).eq.',') then
            if (j.ne.k) write (77,'(a)') line(k:j-1)
            k=j+1
          endif
        enddo
      endif
      goto 1
    2 close (76)
      rewind (77)

partition=.false.
eps=1.d-10

read (77,*)
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  x1=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  x1=param(jd)
else
  read (line,*) x1
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  y1=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  y1=param(jd)
else
  read (line,*) y1
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  x2=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  x2=param(jd)
else
  read (line,*) x2
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  y2=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  y2=param(jd)
else
  read (line,*) y2
endif
x1=(x1-xlon1)/(xlon2-xlon1)*xl
y1=(y1-xlat1)/(xlat2-xlat1)*yl
x2=(x2-xlon1)/(xlon2-xlon1)*xl
y2=(y2-xlat1)/(xlat2-xlat1)*yl

read (77,'(a1024)') line
if (scan(line,':').ne.0) stop 'faultmodel cannot be specified as a range'
backspace (77)
read (77,*) faultmodel
if (faultmodel > 2 .or. faultmodel < 1) stop 'faultmodel must have a value of either 1 or 2'

do i=1,nfault
  read (77,*) fault(i)%n
  if (fault(i)%n.lt.0) then
    allocate (fault(i)%x(4),fault(i)%y(4))
    do k=1,4
      read (77,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%x(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%x(k)=param(jd)
      else
        read (line,*) fault(i)%x(k)
      endif
    enddo
  else
    allocate (fault(i)%x(fault(i)%n),fault(i)%y(fault(i)%n))
    allocate (fault(i)%samef(fault(i)%n),fault(i)%samex(fault(i)%n))
    allocate (fault(i)%samey(fault(i)%n))
    do k=1,fault(i)%n
      read (77,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%x(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%x(k)=param(jd)
      else
        read (line,*) fault(i)%x(k)
      endif
      read (77,'(a1024)') line
      icol=scan(line,':')
      jcol=scan(line,'#')
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%y(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%y(k)=param(jd)
      else
        read (line,*) fault(i)%y(k)
      endif
      !read (77,'(a1024)') line
      !if (scan(line,':').ne.0) stop 'samef cannot be specified as a range'
      !backspace (77)
      read (77,*) fault(i)%samef(k)
      if (fault(i)%samef(k) > nfault) then
        write (*,*) ''
        write (*,'(a)') '*** Error in fault_parameters.txt file ***'
        write (*,'(a,i3,a,i3,a)') 'Fault ',i,' is trying to use the same geometry as fault ',&
                                  fault(i)%samef(k),', but'
        write (*,'(a,i3,a)') 'there are only ',nfault,' defined faults.'
        write (*,*) ''
        stop
      endif
      !read (77,'(a1024)') line
      !if (scan(line,':').ne.0) stop 'samex cannot be specified as a range'
      !backspace (77)
      read (77,*) fault(i)%samex(k)
      if (fault(i)%samex(k) > 0) then
        if (fault(i)%samex(k) > fault(fault(i)%samef(k))%n) then
          write (*,*) ''
          write (*,'(a)') '*** Error in fault_parameters.txt file ***'
          write (*,'(a,i3,a,i3,a,i3,a)') 'Fault ',i,' is trying to use the x position of particle ',&
                                       fault(i)%samex(k),' on fault',fault(i)%samef(k),','
          write (*,'(a,i3,a,i3)') 'but there are only ',fault(fault(i)%samef(k))%n,&
                                  ' points defining fault ',fault(i)%samef(k)
          write (*,*) ''
          stop
        endif
      endif
      !read (77,'(a1024)') line
      !if (scan(line,':').ne.0) stop 'samey cannot be specified as a range'
      !backspace (77)
      read (77,*) fault(i)%samey(k)
      if (fault(i)%samey(k) > 0) then
        if (fault(i)%samey(k) > fault(fault(i)%samef(k))%n) then
          write (*,*) ''
          write (*,'(a)') '*** Error in fault_parameters.txt file ***'
          write (*,'(a,i3,a,i3,a,i3,a)') 'Fault ',i,' is trying to use the y position of particle ',&
                                       fault(i)%samey(k),' on fault',fault(i)%samef(k),','
          write (*,'(a,i3,a,i3)') 'but there are only ',fault(fault(i)%samef(k))%n,&
                                  ' points defining fault ',fault(i)%samef(k)
          write (*,*) ''
          stop
        endif
      endif
      fault(i)%y(k)=fault(i)%y(k)+zl
    enddo
    fault(i)%x1=x1;fault(i)%y1=y1;fault(i)%x2=x2;fault(i)%y2=y2
  endif
  read (77,*) fault(i)%nstep
  allocate (fault(i)%timestart(fault(i)%nstep),fault(i)%timeend(fault(i)%nstep), &
            fault(i)%velo(fault(i)%nstep),fault(i)%conv_rate(fault(i)%nstep),  &
            fault(i)%conv_factor(fault(i)%nstep))
  do k=1,fault(i)%nstep
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    kcol=scan(line,'*')
    if (partition .and. i == 2) then
      fault(i)%timestart(k)=fault(1)%timestart(k)
      ! Dummy reads below
      if (icol.ne.0) then
        read (line(1:icol-1),*) d1
        read (line(icol+1:1024),*) d1
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) d1
      elseif (kcol.ne.0) then
        fault(i)%timestart(k) = fault(1)%timestart(k)
      else
        read (line,*) d1
      endif
    else
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%timestart(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%timestart(k)=param(jd)
      elseif (kcol.ne.0) then
        if (k.eq.1) stop 'cannot use wild card for first time step'
        fault(i)%timestart(k)=fault(i)%timeend(k-1)
      else
        read (line,*) fault(i)%timestart(k)
      endif
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (partition .and. i == 2) then
      fault(i)%timeend(k)=fault(1)%timeend(k)
      ! Dummy reads below
      if (icol.ne.0) then
        read (line(1:icol-1),*) d1
        read (line(icol+1:1024),*) d1
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) d1
      else
        read (line,*) d1
      endif
    else
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%timeend(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%timeend(k)=param(jd)
      else
        read (line,*) fault(i)%timeend(k)
      endif
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (partition .and. i == 2) then
      fault(i)%velo(k)=fault(1)%velo(k)
      ! Dummy reads below
      if (icol.ne.0) then
        read (line(1:icol-1),*) d1
        read (line(icol+1:1024),*) d1
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) d1
      else
        read (line,*) d1
      endif
    else
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%velo(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%velo(k)=param(jd)
      else
        read (line,*) fault(i)%velo(k)
      endif
    endif
    ! New inputs for constant convergence rate models
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    kcol=scan(line,'*')
    if (partition .and. i == 2) then
      fault(i)%conv_factor(k)=fault(1)%conv_factor(k)
      ! Dummy reads below
      if (icol.ne.0) then
        read (line(1:icol-1),*) d1
        read (line(icol+1:1024),*) d1
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) d1
      elseif (kcol.ne.0) then
        fault(i)%conv_factor(k)=fault(1)%conv_factor(k)
      else
        read (line,*) d1
      endif
    else
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%conv_factor(k)=param(nd)
        if (i == 1 .and. range(1,nd) > -eps .and. range(2,nd) > -eps) partition=.true.
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%conv_factor(k)=param(jd)
      elseif (kcol.ne.0) then
        if (k.eq.1) stop 'Cannot use wild card for first time step'
        fault(i)%conv_factor(k)=fault(i)%conv_factor(k-1)
      else
        read (line,*) fault(i)%conv_factor(k)
        if (i == 1 .and. fault(i)%conv_factor(k) > -eps) partition=.true.
      endif
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    kcol=scan(line,'*')
    if (partition .and. i == 2) then
      fault(i)%conv_rate(k)=fault(1)%conv_rate(k)
      ! Dummy reads below
      if (icol.ne.0) then
        read (line(1:icol-1),*) d1
        read (line(icol+1:1024),*) d1
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) d1
      elseif (kcol.ne.0) then
        fault(i)%conv_rate(k)=fault(1)%conv_rate(k)
      else
        read (line,*) d1
      endif
    else
      if (icol.ne.0) then
        nd=nd+1
        read (line(1:icol-1),*) range(1,nd)
        read (line(icol+1:1024),*) range(2,nd)
        fault(i)%conv_rate(k)=param(nd)
      elseif (jcol.ne.0) then
        read (line(jcol+1:1024),*) jd
        fault(i)%conv_rate(k)=param(jd)
      elseif (kcol.ne.0) then
        if (k.eq.1) stop 'Cannot use wild card for first time step'
        fault(i)%conv_rate(k)=fault(i)%conv_rate(k-1)
      else
        read (line,*) fault(i)%conv_rate(k)
      endif
    endif
    if (partition .and. i == 2) then
      fault(i)%timestart(k)=fault(1)%timestart(k)
      fault(i)%timeend(k)=fault(1)%timeend(k)
    else
      fault(i)%timestart(k)=timeend-fault(i)%timestart(k)
      fault(i)%timeend(k)=timeend-fault(i)%timeend(k)
    endif
  enddo
enddo

read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  dstart=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  dstart=param(jd)
else
  read (line,*) dstart
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  dend=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  dend=param(jd)
else
  read (line,*) dend
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  dvelo=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  dvelo=param(jd)
else
  read (line,*) dvelo
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  douter=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  douter=param(jd)
else
  read (line,*) douter
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  dinner=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  dinner=param(jd)
else
  read (line,*) dinner
endif
dstart=timeend-dstart
dend=timeend-dend

read (77,'(a1024)') line
if (scan(line,':').ne.0) stop 'duplex advection flag cannot be specified as a range'
backspace (77)
read (77,*) advect_duplex
if (advect_duplex /= 0 .and. advect_duplex /= 1 .and. advect_duplex /= -1) then
  write (*,*) 'Bad value, ',advect_duplex,' for duplex advection flag. Input value must be'
  write (*,*) '"0", "1" or "-1". Exiting.'
  stop
endif

read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
  nd=nd+1
  read (line(1:icol-1),*) range(1,nd)
  read (line(icol+1:1024),*) range(2,nd)
  vduplexadv=param(nd)
elseif (jcol.ne.0) then
  read (line(jcol+1:1024),*) jd
  vduplexadv=param(jd)
else
  read (line,*) vduplexadv
endif
close (77)

do i=1,nfault
  if (fault(i)%n.gt.0) then
    allocate (fault(i)%xs(fault(i)%n-1),fault(i)%ys(fault(i)%n-1))
    allocate (fault(i)%alpha(fault(i)%n-1),fault(i)%kbm(fault(i)%n-2))
  endif
enddo

call calculate_fault_parameters (fault,nfault)

do i=1,nfault
  if (fault(i)%n.gt.0) then
    if (fault(i)%x(fault(i)%n).gt.fault(i)%x(1) .and. &
        fault(i)%y(fault(i)%n).lt.fault(i)%y(1)) &
        fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
    if (fault(i)%x(fault(i)%n).lt.fault(i)%x(1) .and. &
        fault(i)%y(fault(i)%n).gt.fault(i)%y(1)) &
        fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
  endif
enddo

do i=1,nfault
  xn=fault(i)%y2-fault(i)%y1
  yn=fault(i)%x1-fault(i)%x2
  xyn=sqrt(xn**2+yn**2)
  fault(i)%xn=xn/xyn
  fault(i)%yn=yn/xyn
  fault(i)%faultmodel=faultmodel
  fault(i)%dstart=dstart
  fault(i)%dend=dend
  fault(i)%dvelo=dvelo
  fault(i)%dinner=dinner
  fault(i)%douter=douter
  fault(i)%advect_duplex=advect_duplex
  fault(i)%vduplexadv=vduplexadv
enddo

return

end subroutine read_in_fault_parameters
