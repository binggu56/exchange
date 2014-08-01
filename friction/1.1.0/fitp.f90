      subroutine fitp(am1,am2,x,y,px,py,rx,ry,w,Ntraj,qp,qfx,qfy,apx,apy,frx,fry,c,err)
!     ----------------------------------------------------------------------
!     subroutine to do least-square fitting to momemtum and non-classical momemtum 
!     full-dimenional fitting with basis f=(1,x,x^2,x^3,y,y^2,y^3,xy)
!     ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      integer*4,intent(IN) :: Ntraj
      real*8, intent(in)   :: am1,am2
      real*8, intent(IN)   :: px(Ntraj),py(Ntraj),w(Ntraj),x(Ntraj),y(Ntraj)
      real*8,dimension(ntraj) :: dxpx,dypx,dxpy,dypy
      integer*4, parameter :: nb = 8, ndim = 2

      real*8 :: f(nb),s(nb,nb),c(nb,2*ndim),dxf(nb),dyf(nb)
      real*8 :: mat(nb,nb),err(2)
      real*8,dimension(ntraj) :: dxrx,dxry,dyrx,dyry,rx,ry,qfx,qfy,qp,apx,apy,arx,ary, &
                                 dxxpx,dxypy,dyypy,dyxpx,frx,fry
      real*8,dimension(nb)    :: dxxf,dxyf,dyxf,dyyf

      c = 0d0
      s = 0d0

      do j=1,nb
        do i=1,Ntraj
          f = (/x(i),x(i)**2,x(i)**3,y(i),y(i)**2,y(i)**3,x(i)*y(i),1d0/)
!          dxf = (/1d0,2d0*x(i),3d0*x(i)**2,0d0,0d0,0d0,y(i),0d0/)
!          dyf = (/0d0,0d0,0d0,1d0,2d0*y(i),3d0*y(i)**2,x(i),0d0/)
          c(j,1) = c(j,1)+px(i)*f(j)*w(i)
          c(j,2) = c(j,2)+py(i)*f(j)*w(i)
          c(j,ndim+1) = c(j,ndim+1)+w(i)*f(j)*rx(i)
          c(j,ndim+2) = c(j,ndim+2)+w(i)*f(j)*ry(i)
        enddo
      enddo

!      cq = -0.5*cq


      do m=1,nb
        do n=1,nb
          do i=1,Ntraj
            f=(/x(i),x(i)**2,x(i)**3,y(i),y(i)**2,y(i)**3,x(i)*y(i),1d0/)
            s(m,n)=w(i)*f(m)*f(n)+s(m,n)
          end do
        end do
      end do

      mat = s
! calculate matrix c(t)
        call DPOSV('U',nb,4,mat,nb,c,nb,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "polynomial fitting : matrix fails"
        stop
        end if
     
!      mat = s
!      call DPOSV('U',nb,2,mat,nb,cq,nb,INFO)
!      if(INFO/=0) then
!      print *, "info=",info
!      print *, "quantum potential : matrix fails"
!      stop
!      end if


!      do i=1,Ntraj
!c          dpx(i) = c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2
!        f=(/x(i),x(i)**2,x(i)**3,y(i),y(i)**2,y(i)**3,x(i)*y(i),1d0/)
!        dxf = (/1d0,2d0*x(i),3d0*x(i)**2,0d0,0d0,0d0,y(i),0d0/)
!        dyf = (/0d0,0d0,0d0,1d0,2d0*y(i),3d0*y(i)**2,x(i),0d0/)
       

      
!c     &             py(i)*(c(1,2)+2d0*c(3,2)*x(i)+c(5,2)*y(i))/am2
!c          dpy(i) = px(i)*(c(2,1)+2d0*c(4,1)*y(i)+c(5,1)*x(i))/am1+
!c     &             py(i)*(c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/am2
!! calculate quantum potential
!c          qp(i) = (c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2)/2d0/am1
!c     &            (c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/2d0/am2
!c          qfx(i) = -(2d0*c(2,1)+6d0*x(i)*c(3,1))/2d0/am1
!c     qfy(i) = -(c(5,1)/2d0/am1+c(4,2)/am2)
      apx = 0d0
      apy = 0d0
      dxpx = 0d0
      dxpy = 0d0
      dypx = 0d0
      dypy = 0d0
      arx = 0d0
      ary = 0d0
      dxrx = 0d0
      dxry = 0d0
      dyrx = 0d0
      dyry = 0d0
      dxxpx = 0d0
      dxypy = 0d0
      dyypy = 0d0
      dyxpx = 0d0
      qp = 0d0
      frx = 0d0
      fry = 0d0

      traj_loop:do i=1,ntraj
        f=(/x(i),x(i)**2,x(i)**3,y(i),y(i)**2,y(i)**3,x(i)*y(i),1d0/)
        dxf = (/1d0,2d0*x(i),3d0*x(i)**2,0d0,0d0,0d0,y(i),0d0/)
        dxxf = (/0d0,2d0,6d0*x(i),0d0,0d0,0d0,0d0,0d0/)
        dyxf = (/0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0/)
        dxyf = (/0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0/)
        dyf = (/0d0,0d0,0d0,1d0,2d0*y(i),3d0*y(i)**2,x(i),0d0/)
        dyyf = (/0d0,0d0,0d0,0d0,2d0,6d0*y(i),0d0,0d0/)

        do j=1,nb
          apx(i) = apx(i)+c(j,1)*f(j)
          apy(i) = apy(i)+c(j,2)*f(j)
          dxpx(i) = dxpx(i)+dxf(j)*c(j,1)
          dxpy(i) = dxpy(i)+dxf(j)*c(j,2)
          dypx(i) = dypx(i)+dyf(j)*c(j,1)
          dypy(i) = dypy(i)+dyf(j)*c(j,2)
          dxxpx(i) = dxxpx(i)+dxxf(j)*c(j,1)
          dxypy(i) = dxypy(i)+dxyf(j)*c(j,2)
          dyypy(i) = dyypy(i)+dyypy(j)*c(j,2)
          dyxpx(i) = dyxpx(i)+dyxf(j)*c(j,1)
        enddo

        do j=1,nb
          arx(i)  = arx(i)   +c(j,1+ndim)*f(j)
          ary(i)  = ary(i)   +c(j,2+ndim)*f(j)
          dxrx(i) = dxrx(i) +dxf(j)*c(j,1+ndim)
          dyrx(i) = dyrx(i) +dyf(j)*c(j,1+ndim)
          dyry(i) = dyry(i) +dyf(j)*c(j,2+ndim)
          dxry(i) = dxry(i) +dxf(j)*c(j,2+ndim)
        enddo
        
        qp(i) = -1d0/2d0/am1*(arx(i)**2+dxrx(i))-1d0/2d0/am2*(ary(i)**2+dyry(i))
        qfx(i) = 1d0/2d0/am1*(2d0*arx(i)*dxrx(i)+2d0*c(2,1+ndim)+6d0*x(i)*c(3,1+ndim))+ &
                 1d0/2d0/am2*(2d0*ary(i)*dxry(i)+c(7,2+ndim))
        qfy(i) = 1d0/2d0/am1*(2d0*arx(i)*dyrx(i)+c(7,1+ndim))+                      &
                1d0/2d0/am2*(2d0*ary(i)*dyry(i)+2d0*c(5,2+ndim)+6d0*y(i)*c(6,2+ndim))
!      enddo

        frx(i) = -((dxpx(i)*arx(i)/am1+dxpy(i)*ary(i)/am2)+   &
                 (dxxpx(i)/2d0/am1+dxypy(i)/2d0/am2))
        fry(i) = -((dypx(i)*arx(i)/am1+dypy(i)*ary(i)/am2)+   &
                 (dyxpx(i)/2d0/am1+dyypy(i)/2d0/am2))
      enddo traj_loop

      call err_fit(ntraj,rx,arx,err_rx)
      call err_fit(ntraj,ry,ary,err_ry)

      err(1) = err_rx
      err(2) = err_ry

!      write(*,*) 'rx',arx(1),arx(2),arx(3)
!      write(*,*) 'mass', am1,am2
!      write(*,*) 'quantum potential',qp(1),qp(2),qp(3)
      return
      end subroutine

!c     --------------------------------
!c     check the error of poly-fit of p
!c     --------------------------------
!      subroutine errp(nt,p,x,c,err)
!      implicit real*8(a-h,o-z)
!
!      integer*4, intent(in)  :: nt
!      real*8,    intent(in)  :: p(nt),x(nt),c(4,1)
!      real*8,    intent(out) :: err
!      
!      common/smaple/xmin,xmax,dx
!     
!      err = 0d0
!      do i=1,nt
!        pfit = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)*x(i)**3+c(4,1)
!        err = err+(p(i)-pfit)**2*dx
!      enddo
!
!      return
!      end subroutine



      subroutine err_fit(nt,r,ar,err)
      implicit real*8(a-h,o-z)

      integer*4, intent(in)  :: nt
      real*8,    intent(in)  :: r(nt),ar(nt)
      real*8,    intent(out) :: err
     
      err = 0d0
      do i=1,nt
!         = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)*x(i)**3+c(4,1)
        err = err+(r(i)-ar(i))**2
      enddo

      return
      end subroutine
