subroutine calculate_topo_advection (time,fault,nfault,iflag,zl,advect_topo,   &
                                     vtopo,vxtopo,vytopo,vxtopo0,iproc,nd)

! Calculates the advection velocity for topographic advection, if enabled

use definitions

implicit none

type(faulttype) fault(nfault)
double precision :: time,zl,vtopo,vxtopo,vytopo,vxtopo0
integer nfault,iflag,advect_topo,iproc,nd

integer :: i,k
double precision :: maxxs(nfault),velo(nfault)
double precision :: eps,xn,yn

eps    = 1.d-10

maxxs  = 0.d0
velo   = 0.d0
vxtopo = 0.d0

! Find optimal topo advection velocity
do i=1,nfault
  ! Find velocity for all faults
  do k=1,fault(i)%nstep
    if (advect_topo == 1) then                                                  ! Calculate optimal topo advection rate for all times
      if (fault(i)%conv_factor(k).lt.-eps) then
  	    velo(i)=fault(i)%velo(k)
      else
        ! NOTE: This section should probably be adjusted to handle the original
        ! Pecube fault model the same way as the kink-band model
        !
        ! dwhipp - 23.04.2015
        if (fault(i)%x(2) - fault(i)%x(1) >= 0.0) then
          velo(i)=fault(i)%conv_rate(k)*(fault(i)%conv_factor(k))
        elseif (fault(i)%x(2) - fault(i)%x(1) < 0.0) then
          velo(i)=fault(i)%conv_rate(k)-fault(1)%conv_rate(k)*(fault(1)%conv_factor(k))
        endif
        if (i > 2) then
          write (*,*) 'Use of constant convergence model with more than 2 fault'
          write (*,*) 'segments is not supported. Exiting.'
          stop
        endif
      endif
    elseif (advect_topo == -1) then                                             ! Calculate optimal topo advection rate by time step
      write (*,*) 'The per-step optimal topographic advection calculation advect_topo=-1'
      write (*,*) 'is currently disabled. Consider using the standard optimal advection'
      write (*,*) 'velocity calculation (advect_topo=1) instead.'
      write (*,*) 'Stopping.'
      stop
    !  if ((time-fault(i)%timestart(k))*(time-fault(i)%timeend(k)).le.eps) then
    !    if (fault(i)%conv_factor(k).lt.-eps) then
    !	    velo(i)=fault(i)%velo(k)
    !    else
    !      if (i.eq.1) then
    !        velo(i)=fault(i)%conv_rate(k)*(fault(i)%conv_factor(k))
    !      elseif (i.eq.2) then
    !        velo(i)=fault(i)%conv_rate(k)-fault(1)%conv_rate(k)*(fault(1)%conv_factor(k))
    !      else
    !        write (*,*) 'Use of constant convergence model with more than 2 fault'
    !        write (*,*) 'segments is not supported. Exiting.'
    !        stop
    !      endif
    !    endif
    !  endif
    else
      write (*,*) 'Bad value, ',advect_topo,' for topographic advection flag. Input value must be'
      write (*,*) '"0", "1" or "-1". Exiting.'
      stop
    endif
  enddo

  ! Find maximum horizontal velocity component from all faults
  do k=1,fault(i)%n-1
    if (.not.(fault(i)%y(k+1)-zl >= 0.d0 .and. fault(i)%y(k)-zl >= 0.d0)) then
      if (abs(fault(i)%xs(k)) > abs(maxxs(i))) maxxs(i) = fault(i)%xs(k)
    endif
  
    ! Calculate advection velocity from velo and maxxs, and update topo
    ! advection velocity if this fault has thrust-sense motion
    if (velo(i) < 0.d0) then
      if (fault(i)%x(k+1)-fault(i)%x(k) < 0.d0) then
        if (abs(velo(i)*maxxs(i)) > vxtopo) then
          vxtopo=velo(i)*maxxs(i)
          xn=fault(i)%xn
          yn=fault(i)%yn
        endif
      endif
    elseif (velo(i) > 0.d0) then
      write (*,*) 'Calculation of the topographic advection velocity is not supported'
      write (*,*) 'for normal-sense faults. Either specify an advection velocity in'
      write (*,*) 'the input file or disable the use of topographic advection.'
      write (*,*) 'Exiting.'
      stop
    endif
  enddo
enddo

! If user provided an advection velocity, use that with best normals from above
if (vtopo >= 0.d0) then
  vxtopo0=vtopo
  ! Here we reorientate the velocity if needed
  if (iflag.eq.1) then
    vytopo=vtopo*yn
    vxtopo=vtopo*xn
  endif
! Otherwise, use optimal topo advection velocity and normals
else
  vxtopo0=vxtopo
  ! Here we reorientate the velocity if needed
  if (iflag.eq.1) then
    vytopo=vxtopo*yn
    vxtopo=vxtopo*xn
  endif
endif

if (iproc.eq.0.and.nd.eq.0) then
  if (iflag == 1) then
    write (*,'(a,f8.4,a)') ' Topographic advection x-velocity: ',vxtopo,' mm/yr'
    write (*,'(a,f8.4,a)') ' Topographic advection y-velocity: ',vytopo,' mm/yr'
  else
    write (*,'(a,f8.4,a)') ' Topographic advection velocity normal to fault trace: ',vxtopo0,' mm/yr'
  endif
endif

end subroutine calculate_topo_advection