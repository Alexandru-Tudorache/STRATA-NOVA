module observer_projection
  implicit none
  private
  public :: equal_arrival_time_curve, linear_interp
  real*8, parameter :: c_light = 2.99792458d10
  real*8, parameter :: pi_mgr = 3.1415926535897932384626433832795d0
contains

  subroutine equal_arrival_time_curve(t,lum,radius,n,nmu,observer)
    integer, intent(in) :: n,nmu
    real*8, intent(in) :: t(n),lum(n),radius(n)
    real*8, intent(out) :: observer(n)
    real*8, allocatable :: mu(:),w(:)
    real*8 :: delay,target,val,sumv,min_target
    integer :: k,j

    allocate(mu(nmu),w(nmu))
    call gauss_legendre_01(nmu,mu,w)
    ! integrate delayed luminosity over the visible hemisphere
    do k=1,n
       delay=max(radius(k),0.0d0)/c_light
       min_target=t(k)-delay
       ! retain the input sample where the earlier history is unavailable
       if (min_target.lt.t(1)) then
          observer(k)=lum(k)
       else
          sumv=0.0d0
          do j=1,nmu
             target=t(k)-delay*(1.0d0-mu(j))
             val=linear_interp(t,lum,target)
             sumv=sumv+2.0d0*mu(j)*w(j)*val
          enddo
          observer(k)=sumv
       endif
    enddo
    deallocate(mu,w)
  end subroutine equal_arrival_time_curve


  subroutine gauss_legendre_01(n,x,w)
    integer, intent(in) :: n
    real*8, intent(out) :: x(n),w(n)
    integer :: i,j,m
    real*8 :: z,z1,p1,p2,p3,pp,eps,xm,xl

    eps=3.0d-14
    m=(n+1)/2
    xm=0.5d0
    xl=0.5d0
    ! solve for gauss-legendre nodes and weights on the unit interval
    do i=1,m
       z=cos(pi_mgr*(dble(i)-0.25d0)/(dble(n)+0.5d0))
       do
          p1=1.0d0; p2=0.0d0
          do j=1,n
             p3=p2; p2=p1
             p1=((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/dble(j)
          enddo
          pp=dble(n)*(z*p1-p2)/(z*z-1.0d0)
          z1=z
          z=z1-p1/pp
          if (abs(z-z1).le.eps) exit
       enddo
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)*xl*2.0d0
       w(n+1-i)=w(i)
    enddo
  end subroutine gauss_legendre_01


  ! interpolate between bracketing samples with endpoint extension
  real*8 function linear_interp(x,y,xq) result(v)
    real*8, intent(in) :: x(:),y(:),xq
    integer :: lo,hi,mid,n
    real*8 :: f
    n=size(x)
    if (xq.le.x(1)) then
       v=y(1); return
    else if (xq.ge.x(n)) then
       v=y(n); return
    endif
    ! locate the bracketing samples before linear interpolation
    lo=1; hi=n
    do while(hi-lo.gt.1)
       mid=(lo+hi)/2
       if (x(mid).le.xq) then
          lo=mid
       else
          hi=mid
       endif
    enddo
    if (x(hi).le.x(lo)) then
       v=y(lo)
    else
       f=(xq-x(lo))/(x(hi)-x(lo))
       v=(1.0d0-f)*y(lo)+f*y(hi)
    endif
  end function linear_interp
end module observer_projection
