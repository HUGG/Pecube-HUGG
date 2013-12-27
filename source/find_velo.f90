      subroutine find_velo (x,y,z,time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,&
                            ymax,vxtopo,vytopo,vxtopo0)

! this extra layer between Pecube and the velocity function geometry ensures
! that the advection of the "lagrangian particles" for which we track the thermal
! history (to calculate the age) is second order accurate; it uses a mid-point
! algorithm.

! Note that this utility has been turned off in this new version of Pecube
! this routine is useless but has been kept for compatibility with ancient versions

      use definitions

      implicit none

      double precision :: x,y,z,time,dt,vx,vy,vz,xmin,xmax,ymin,ymax
      double precision :: vxtopo,vytopo,vxtopo0
      integer :: nfault

      type (faulttype) fault(nfault)

      call find_velocity (x,y,z,vx,vy,vz,fault,nfault,time,1,xmin,xmax,ymin,   &
                          ymax,vxtopo,vytopo,vxtopo0)

      return

      end
