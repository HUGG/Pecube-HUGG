!---------------

subroutine calculate_fault_parameters (fault,nfault)

! given the fault geometry recalculates fault parameters needed by other routines

use definitions

implicit none

type(faulttype) fault(nfault)
integer nfault,i,k
double precision xn,yn,xyn,pi

pi=4.d0*atan(1.d0)

! If any fault nodes should have the same x or y coordinates as another fault,
! modify the geometry here
! dwhipp - 04/12
do i=1,nfault
  do k=1,fault(i)%n
    if (fault(i)%samef(k) > 0) then
      if (fault(i)%samex(k) > 0) fault(i)%x(k)=fault(fault(i)%samef(k))%x(fault(i)%samex(k))
      if (fault(i)%samey(k) > 0) fault(i)%y(k)=fault(fault(i)%samef(k))%y(fault(i)%samey(k))
    endif
  enddo
enddo

do i=1,nfault
  do k=1,fault(i)%n-1
    xn=fault(i)%x(k+1)-fault(i)%x(k)
    yn=fault(i)%y(k+1)-fault(i)%y(k)
    xyn=sqrt(xn**2+yn**2)
    xn=xn/xyn
    yn=yn/xyn
    fault(i)%xs(k)=xn
    fault(i)%ys(k)=yn
    fault(i)%alpha(k)=fault(i)%ys(k)/fault(i)%xs(k)
  enddo

! Define kink-band planes for kink-band fault kinematic model
  do k=1,fault(i)%n-2
    !fault(i)%kbm(k)=tan((atan(fault(i)%alpha(k))+atan(fault(i)%alpha(k+1)))/2.d0 + pi/2.d0)
    fault(i)%kbm(k)=(fault(i)%alpha(k)+fault(i)%alpha(k+1))/2.d0 + tan(pi/2.d0)
  enddo    
enddo

return

end subroutine calculate_fault_parameters
