!---------------------

subroutine find_velocity (x,y,z,vx,vy,vz,fault,nfault,time,iflag,xmin,xmax,    &
                          ymin,ymax,vxtopo,vytopo,vxtopo0)

! this subroutine returns a velocity vector at a given location
! and time from a set of faults; note that the fault gemetry
! must be compatible with the time (it is supposed to be updated
! somewhere else)
! when iflag=1, global coordinates are used
! when iflag=0, local "fault" coordinates are used

use definitions

implicit none

type(faulttype) fault(nfault)
double precision :: x,y,z,vx,vy,vz,eps,vxf,vyf,vzf,velo,time,xp,yp
double precision :: r,s,h1,h2,h3,h4,xmin,xmax,ymin,ymax,pi
double precision :: vxtopo,vytopo,vxtopo0,velohw,velofw
integer :: nv,nfault,i,k,iflag
logical :: velo_assigned

! Note (x-x0)*xn+(y-y0)*yn is positive when x,y is on top of the plane
! passing through (x0,y0) and having (xn,yn) as a normal

eps=1.d-10
pi=4.d0*atan(1.d0)

vx=0.d0
vy=0.d0
vz=0.d0

do i=1,nfault
  ! first we rotate/translate the coordinates of the point, if needed
  if (iflag.eq.1.and.fault(i)%n.gt.0) then
    xp=(x-fault(i)%x1)*fault(i)%xn+(y-fault(i)%y1)*fault(i)%yn
    yp=y-fault(i)%yn*xp
    if ((yp-fault(i)%y1)*(yp-fault(i)%y2).gt.eps) goto 333
    !yp=-(x-fault(i)%x1)*fault(i)%yn+(y-fault(i)%y1)*fault(i)%xn
    !if (yp*(yp-sqrt((fault(i)%x1-fault(i)%x2)**2+(fault(i)%y1-fault(i)%y2)**2)).gt.0.d0) goto 333
    yp=z
  else
    xp=x
    yp=z
  endif

  velo=0.d0
  do k=1,fault(i)%nstep
    if ((time-fault(i)%timestart(k))*(time-fault(i)%timeend(k)).le.eps) then
      if (fault(i)%conv_factor(k).lt.-eps) then
        velo=fault(i)%velo(k)
      else
        if (i.eq.1) then
          velo=fault(i)%conv_rate(k)*(fault(i)%conv_factor(k))
          velofw = velo
          velohw = fault(i)%conv_rate(k) - velofw
        elseif (i.eq.2) then
          velo=fault(i)%conv_rate(k)-fault(1)%conv_rate(k)*(fault(1)%conv_factor(k))
        else
          write (*,*) 'Use of constant convergence model with more than 2 fault'
          write (*,*) 'segments is not supported. Exiting.'
          stop
        endif
      endif
    endif
  enddo
  nv=0
  vxf=0.d0
  vyf=0.d0
  vzf=0.d0

  if (fault(i)%n.lt.0) then
    nv=nv+1
    r=(x-xmin)/(xmax-xmin)*2.d0-1.d0
    s=(y-ymin)/(ymax-ymin)*2.d0-1.d0
    h1=(1.d0-r)*(1.d0-s)/4.d0
    h2=(1.d0+r)*(1.d0-s)/4.d0
    h3=(1.d0-r)*(1.d0+s)/4.d0
    h4=(1.d0+r)*(1.d0+s)/4.d0
    vzf=vzf+velo*(h1*fault(i)%x(1)+h2*fault(i)%x(2)+h3*fault(i)%x(3)+h4*fault(i)%x(4))
  else
    if (fault(i)%faultmodel == 2) then
      if (i == 1) then
        velo_assigned=.false.
        if (fault(i)%x(2) - fault(i)%x(1) >= 0.0) then
          ! Do this part of the loop if the fault points are given in order of
          ! decreasing x value
          !
          ! This could probably be done in a smarter way, but it works...
          do k = 1,fault(i)%n-2
            if (xp-fault(i)%x(k+1) < (1.d0/fault(i)%kbm(k))*(yp-fault(i)%y(k+1))) then
              if (yp-fault(i)%y(k) >= fault(i)%alpha(k)*(xp-fault(i)%x(k))) then
                vxf=-velohw*cos(atan(fault(i)%alpha(k)))
                vzf=-velohw*sin(atan(fault(i)%alpha(k)))
                velo_assigned=.true.
                exit
              else
                vxf=velofw*cos(atan(fault(i)%alpha(k)))
                vzf=velofw*sin(atan(fault(i)%alpha(k)))
                velo_assigned=.true.
                exit
              endif
            endif
          enddo
        else
          ! Do this part of the loop if the fault points are given in order of
          ! increasing x value
          do k = 1,fault(i)%n-2
            if (xp-fault(i)%x(k+1) >= (1.d0/fault(i)%kbm(k))*(yp-fault(i)%y(k+1))) then
              if (yp-fault(i)%y(k) >= fault(i)%alpha(k)*(xp-fault(i)%x(k))) then
                vxf=-velohw*cos(atan(fault(i)%alpha(k)))
                vzf=-velohw*sin(atan(fault(i)%alpha(k)))
                velo_assigned=.true.
                exit
              else
                vxf=velofw*cos(atan(fault(i)%alpha(k)))
                vzf=velofw*sin(atan(fault(i)%alpha(k)))
                velo_assigned=.true.
                exit
              endif
            endif
          enddo
        endif
        if (.not.velo_assigned) then
          if (yp-fault(i)%y(fault(i)%n-1) >= fault(i)%alpha(fault(i)%n-1)*(xp-fault(i)%x(fault(i)%n-1))) then
            vxf=-velohw*cos(fault(i)%alpha(fault(i)%n-1))
            vzf=-velohw*sin(fault(i)%alpha(fault(i)%n-1))
            velo_assigned=.true.
          else
            vxf=velofw*cos(fault(i)%alpha(fault(i)%n-1))
            vzf=velofw*sin(fault(i)%alpha(fault(i)%n-1))
            velo_assigned=.true.
          endif
        endif
      endif
    else          
      ! Use original fault model
      if (fault(i)%faultmodel /= 1) then
        write (*,'(a,i2,a)') 'Unexpected fault model type (',fault(i)%faultmodel,'). Using original fault model'
      endif
      ! first check segment by segment
      do k=1,fault(i)%n-1
        ! check if we are in the part of the plane that is on the top side of the segment
        if ((xp-fault(i)%x(k))*fault(i)%xs(k)+(yp-fault(i)%y(k))*fault(i)%ys(k).lt.eps) goto 111
        if ((xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k).lt.eps) goto 111
        if ((xp-fault(i)%x(k+1))*fault(i)%xs(k)+(yp-fault(i)%y(k+1))*fault(i)%ys(k).gt.eps) goto 111
        nv=nv+1
        vxf=vxf+velo*fault(i)%xs(k)
        vzf=vzf+velo*fault(i)%ys(k)
111 continue
      enddo
      ! then point by point
      do k=1,fault(i)%n
        ! now check if we are in a sector between two "diverging" segments
        if (k.ne.1) then
          if ((xp-fault(i)%x(k))*fault(i)%ys(k-1)-(yp-fault(i)%y(k))*fault(i)%xs(k-1).lt.eps) goto 222
          if ((xp-fault(i)%x(k))*fault(i)%xs(k-1)+(yp-fault(i)%y(k))*fault(i)%ys(k-1).lt.eps) goto 222
        endif
        if (k.ne.fault(i)%n) then
          if ((xp-fault(i)%x(k))*fault(i)%xs(k)+(yp-fault(i)%y(k))*fault(i)%ys(k).gt.eps) goto 222
          if ((xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k).lt.eps) goto 222
        endif
        if (k.ne.1) then
          nv=nv+1
          vxf=vxf+velo*fault(i)%xs(k-1)
          vzf=vzf+velo*fault(i)%ys(k-1)
        endif
        if (k.ne.fault(i)%n) then
          nv=nv+1
          vxf=vxf+velo*fault(i)%xs(k)
          vzf=vzf+velo*fault(i)%ys(k)
        endif
222 continue
      enddo
    endif
    ! here we reorientate the velocity if needed
    if (iflag.eq.1) then
      vyf=vxf*fault(i)%yn
      vxf=vxf*fault(i)%xn
    endif
  endif
  if (fault(i)%faultmodel == 2) then
    if (i == 1) then
      vx=vxf
      vy=vyf
      vz=vzf
    endif
  else
    if (nv.gt.0) then
      vx=vx+vxf/nv
      vy=vy+vyf/nv
      vz=vz+vzf/nv
    endif
  endif
333 continue
enddo

! Optional zone of accretion from F. Herman
! Applied regardless of fault activity
! Rectangular erosion function
if (time.gt.fault(1)%dstart .and. time.le.fault(1)%dend) then
  if (-xp.lt.fault(1)%dinner) then
!  if (-xp.lt.fault(1)%dinner .and. -xp.gt.fault(1)%douter) then
!    vz=vz+fault(1)%dvelo
    vz=vz+fault(1)%dvelo*(atan((-xp-fault(1)%douter))+pi/2.)/pi
  endif
endif

if (iflag == 1) then
  vx=vx-vxtopo
  vy=vy-vytopo
else
  vx=vx-vxtopo0
endif

return

end subroutine find_velocity
