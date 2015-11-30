      subroutine qpot(am,x,ntraj,ndim)

      implicit real*8(a-h,o-z)
      
! --- linear quantum force subroutine; fit {r} with linear basis {1,x1,x2,...}            
      integer*4  :: Ntraj,ndim
      real*8     :: w(Ntraj),x(ndim,Ntraj),am(ndim)
      real*8     :: qf(ndim,Ntraj),qp(Ntraj)
      real*8 :: f(ndim+1),r(ndim,Ntraj)
      real*8 :: s(ndim+1,ndim+1),c(ndim+1,ndim)
      
!---- define c matrix
      c=0d0
      do j=1,ndim
        c(j,j)=-0.5d0
      enddo

! --- basis vector f = (1,x(1),x(2),...,x(ndim)) for each trajectory
      f = 0d0
      s = 0d0 
      do i=1,Ntraj
        do j=1,ndim
          f(j)=x(j,i)
        enddo
        f(ndim+1)=1d0

! --- matrix S=f X f [ndim+1,ndim+1]
         do k1=1,ndim+1
          do k2=1,k2
            s(k1,k2)=w(i)*f(k1)*f(k2)+s(k1,k2)
          enddo
        enddo

      enddo 
     
      do k1=1,ndim+1
        do k2=1,ndim+1
          s(k2,k1) = s(k1,k2)
        enddo 
      enddo 

!-------- calculate matrix c(t) --------------------
      call DPOSV('U',ndim+1,ndim,s,ndim+1,c,ndim+1,INFO)
      if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
      end if

!----- the momentum operator r=cf-------------
      r=0d0
      qp=0d0
!	qf=0d0

! --- calculate quantum potential
      do i=1,Ntraj
        do j=1,ndim
          f(j)=x(j,i)
        enddo
        f(ndim+1)=1d0

        do j=1,ndim
          r(j,i) = 0d0
          do n=1,ndim+1
              r(j,i)=c(n,j)*f(n)+r(j,i)
        !r(2,i)=c(1,2)*x(1,i)+c(2,2)*x(2,i)+c(3,2)
          enddo 
          qp(i)=-r(j,i)**2/(2d0*am(j))-c(j,j)/(2d0*am(j))+qp(i)
        enddo
      enddo

! --- quantum force
      do i=1,Ntraj
        do j=1,ndim
          qf(j,i) = 0d0
	    do n=1,ndim
        !qfx(i)=rx(i)*c(1,1)/m1+ry(i)*c(1,2)/m2
        !qfy(i)=rx(i)*c(2,1)/m1+ry(i)*c(2,2)/m2
	      qf(j,i)=r(n,i)*c(j,n)/am(n)+qf(j,i)
	!qf(2,i)=r(2,i)*c(2,1)/m(1)+r(2,i)*c(2,2)/m(2)
          enddo
        enddo
      enddo

      return
      end subroutine
