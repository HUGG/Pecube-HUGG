      subroutine find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                          nz,timesurf,timesurfp,istep,eps,ilog, &
                          dt,ntime,istatic,fault,nfault,usurf)

      use definitions

      implicit none

!      implicit real*8 (a-h,o-z)

      type (faulttype) fault(nfault)

      double precision zsurf(nsurf),zsurfp(nsurf),usurf(nsurf)
      double precision zl,diffusivity,timesurf,timesurfp,eps,dt
      double precision dt1,dt2,dt3,Pecletz,Pesurf,dzmax
      integer nsurf,nz,istep,ilog,ntime,istatic,nfault,i,k

! first constrain on time step from conduction

      dt1=zl**2/diffusivity/100.

! second constrain on time step from advection

      dt2=dt1
      Pecletz=0.d0
        do i=1,nfault
          do k=1,fault(i)%nstep
            if ((timesurf-fault(i)%timestart(k))*(timesurf-fault(i)%timeend(k)).le.0.d0) then
              if (fault(i)%conv_factor(k) > eps) then
                Pecletz=max(Pecletz,abs(fault(i)%conv_rate(k)))
              else
                Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
              endif
            endif
            if ((timesurfp-fault(i)%timestart(k))*(timesurfp-fault(i)%timeend(k)).le.0.d0) then
              if (fault(i)%conv_factor(k) > eps) then
                Pecletz=max(Pecletz,abs(fault(i)%conv_rate(k)))
              else
                Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
              endif
            endif
            if ((fault(i)%timestart(k)-timesurf)*(fault(i)%timestart(k)-timesurfp).le.0.d0) then
              if (fault(i)%conv_factor(k) > eps) then
                Pecletz=max(Pecletz,abs(fault(i)%conv_rate(k)))
              else
                Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
              endif
            endif
          enddo
        enddo
      Pecletz=max(Pecletz,maxval(abs(usurf)))
      if (abs(Pecletz).gt.eps) dt2=zl/abs(Pecletz)/100.

! third constrain on time step from surface lowering

      dt3=dt1
        if (istep.ne.0) then
        dzmax=0.
          do i=1,nsurf
          dzmax=max(dzmax,zsurfp(i)-zsurf(i))
          enddo
        Pesurf=dzmax/(timesurf-timesurfp)
!        if (abs(Pesurf).gt.eps) dt3=dzmin/Pesurf/5.
        if (abs(Pesurf).gt.eps) dt3=zl/Pesurf/5.
        endif

! find optimum time step and number of steps

      dt=min(dt1,dt2)
      dt=min(dt,dt3)

      if (istep.ne.0) dt=min(dt,timesurf-timesurfp)
      ntime=int((timesurf-timesurfp)/dt)
      ntime=ntime+1
      dt=(timesurf-timesurfp)/ntime
      istatic=0

        if (istep.eq.0) then
        ntime=1
        dt=0.
        istatic=1
        endif

      if (ilog.eq.1) write (9,*) 'ntime/dt/dt1/dt2/dt3= ',ntime,dt,dt1,dt2,dt3

      return
      end
