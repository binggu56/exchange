      program main
c     ---------DESCRIPTION-----------------------------------------------
c     friction 2-D
c     lagrangian frame
c     qubic basis with linear coupling decomposing into two-step fitting
c     apply for harmonic system (2-dim with linear coupling)
c     completed 
c     ---------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8 :: ki,err(2),cr(8,4)
      integer*4, parameter :: ndim = 2
     	real*8,allocatable   :: ke(:),pe(:),x(:),y(:),px(:),py(:),w(:)
      real*8,allocatable   ::  qfx(:),qfy(:),qp(:),fx(:),c0(:),
     +                       fy(:),s(:),dxpx(:),dypx(:),dxpy(:),
     +                       dypy(:),rx(:),ry(:),apx(:),apy(:),
     +                       frx(:),fry(:),arx(:),ary(:)
!     	real*8 :: f(3)
	complex*16 :: im
	real :: gasdev

      open(100,file = 'en.dat')
      open(101,file = 'w.dat')
      open(102,file = 'traj.dat')  
      open(103,file = 'wft.dat')
      open(104,file = 'err.dat')
      open(105,file = 'fit.dat')

      open(5,file='IN')
      read(5,*) Ntraj
      read(5,*) kmax,dt
      read(5,*) am1,am2
      read(5,*) idum1,idum2
      read(5,*) ax,ay
      read(5,*) qx0,qy0
      read(5,*) px0,py0
      read(5,*) fc
      close(5)
      write(*,1011) 
1011  format('Quantum Trajectories Code with Friction'/,
     +       'version : 1.1'/,
     +       'date    : 1/17/2014'/)  

	allocate(ke(Ntraj),pe(Ntraj),x(Ntraj),
     &         y(Ntraj),px(Ntraj),py(Ntraj),s(Ntraj),
     &         w(Ntraj),qfx(Ntraj),qfy(Ntraj),
     &         qp(Ntraj),fx(Ntraj),fy(Ntraj) ,
     +         rx(ntraj),ry(ntraj),c0(ntraj))

      allocate(dxpx(ntraj),dxpy(ntraj),dypx(ntraj),dypy(ntraj),
     +         apx(Ntraj),apy(Ntraj),
     +         arx(ntraj),ary(Ntraj),
     +         frx(ntraj),fry(ntraj))

      dt2=dt/2d0
      t=0d0
      im=(0d0,1d0)

      pi = 4.d0*atan(1.d0)

      anrm = dsqrt(dsqrt(4d0*ax*ay)/pi)

      ww = 0d0
      do i=1,Ntraj
        x(i)=gasdev(idum1)
        y(i)=gasdev(idum2)
        x(i)=x(i)/sqrt(4d0*ax)+qx0
        y(i)=y(i)/sqrt(4d0*ay)+qy0
!        c0(i) = -ax*(x(i)-qx0)**2-ay*(y(i)-qy0)**2+log(anrm)
        w(i) = 1d0/ntraj
        rx(i) = -2d0*ax*(x(i)-qx0)
        ry(i) = -2d0*ay*(y(i)-qy0)
      enddo 

      write(*,*) 'weights =',ww
       

      px = 0d0
      py = 0d0
        
! 	center for the potential energy
!      x0 = 0d0
!      y0 = 0d0

!     print out the initial conditions        
      write(*,*) 'Initial Conditions'
      print *,'ax =',ax
      print *,'ay =',ay
      write(*,*) 'Number of trajectories =', Ntraj
      write(*,*) 'Time steps =',kmax,dt
      write(*,*) 'Mass =',am1,am2
      print *,'Initial Momentum',px0,py0
      print *,'friction constant ',fc

!     initial value for action function s

c      aver_s = 0d0    
c      do i=1,Ntraj
c        s(i)=px(i)*(x(i)-qx0)+py(i)*(y(i)-qy0)
c        aver_s=aver_s+w(i)*s(i)
c      enddo



! begin the time step loop
      do 10 k=1,kmax
        
        t=t+dt

c        do 11 i=1,Ntraj
c          px(i) = px(i)+(px(i)*dxpx(i)/m1+py(i)*dxpy(i)/m2)*dt2+
c     +            (fx(i)+qfx(i)-gamma*px(i))*dt2
c          py(i) = py(i)+(py(i)*dypy(i)/m2+px(i)*dypx(i)/m1)*dt2+
c     +            (fy(i)+qfy(i))*dt2
c          x(i)=x(i)+px(i)*dt/m1
c          y(i)=y(i)+py(i)*dt/m2
c11      end do
        call derivs(x,y,Ntraj,fx,fy,pe)
        call fitp(am1,am2,x,y,px,py,rx,ry,w,Ntraj,qp,qfx,qfy,apx,
     +            apy,frx,fry,cr,err)

        do i=1,Ntraj
          x(i) = x(i)+apx(i)/am1*dt
          y(i) = y(i)+apy(i)/am2*dt
!          px(i) = px(i)-(px(i)*dxpx(i)/m1+py(i)*dypx(i)/m2)*dt+
!     +            (fx(i)+qfx(i)-fc*px(i)/m1)*dt
!          py(i) = py(i)-(py(i)*dypy(i)/m2+px(i)*dxpy(i)/m1)*dt+
!     +            (fy(i)+qfy(i)-fc*py(i)/m2)*dt
!          c(i) = c(i)-(px(i)/m1*rx(i)+dxpx(i)/2d0/m1+py(i)*ry(i)/m2+
!     +           dypy(i)/2d0/m2)*dt
          px(i) = px(i)+(fx(i)+qfx(i)-fc*px(i))*dt
          py(i) = py(i)+(fy(i)+qfy(i)-fc*py(i))*dt
!     there might be a problem with the mass symbol
          rx(i) = rx(i)+frx(i)*dt
          ry(i) = ry(i)+fry(i)*dt
        enddo

        ke = 0d0
        do i=1,ntraj
          ke(i) = px(i)**2/(2d0*am1)+py(i)**2/(2d0*am2)
        enddo

! calculate action function
c	do i = 1,Ntraj
c	s(i)=s(i)+(ke(i)-pe(i)-qp(i)-gamma*(s(i)-aver_s))*dt
c	end do

c        aver_s = 0d0    
c        do i=1,Ntraj
c        aver_s=aver_s+w(i)*s(i)
c        enddo

! print out a random trajectory	
        write(104,1000) t,err(1),err(2)
        write(105,1000) t,(cr(i,1),i=1,8),(cr(i,2),i=1,8)
        write(102,1000) t,(x(i),i=1,20),(y(i),i=1,20)

c       norm check
!        ww = 0d0
!        do i=1,ntraj
!          ww = ww+w(i)
!        enddo
!        write(101,1000) t,ww

       
! calculate the total energy, the sum of all the trajectories	
      call aver(ntraj,w,pe,po)
      call aver(ntraj,w,ke,ki)
      call aver(ntraj,w,qp,qu)
      
      tot = po+ki+qu
      write(100,1000) t,ki,po,qu,tot
      

! end of time loop
10	end do

c       final wavefunction
        do i=1,ntraj
          write(103,1000) x(i),px(i),y(i),py(i)
        enddo
      

      write(*,1021) tot
1021  format('total energy ',e14.6)

1000  format(100(e14.7,1x))

      close(100)
      close(101)
      close(102)
      close(103)
!      close(104)
!      close(105)
      

      end program main

c     ---------------------------------       
c     classical force
c     ----------------------------------
      subroutine derivs(x,y,Ntraj,fx,fy,v)
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Ntraj
      real*8,intent(IN)  :: x(Ntraj),y(Ntraj)
      real*8,intent(OUT) :: fx(Ntraj),fy(Ntraj),v(Ntraj)
      
      epi = 0.5d0
      alfa = 0.5d0

      do i=1,Ntraj
        fx(i) = -x(i)-4d0*alfa*x(i)**3/2d0-epi*y(i)
        fy(i) = -y(i)-4d0*alfa*y(i)**3/2d0-epi*x(i)
        v(i) = (x(i)**2+alfa*x(i)**4)/2d0+(y(i)**2+alfa*y(i)**4)/2d0
     +         +epi*x(i)*y(i)
      enddo

!------------2-dim for helium solid---------------------
!      a0 = 0d0
!c      a0 = -2.22384d-2
!      a1 = 6.49557d-6
!      a2 = -4.50827d-6
!      a3 = 1.56397d-6
!      a4 = 3.61582d-8
!      a12 = 6.5047d-6
!
!      do i=1,ntraj
!        fx(i) = -2d0*a1*x(i)-4d0*x(i)**3*a2-6d0*x(i)**5*a3-
!     +           a4*8d0*x(i)**7-a12*y(i)
!        fy(i) = -2d0*a1*y(i)-4d0*y(i)**3*a2-6d0*y(i)**5*a3-
!     +           8d0*a4*y(i)**7-a12*x(i)
!        v(i) = a0+a1*x(i)**2+a2*x(i)**4+a3*x(i)**6+a4*x(i)**8+
!     +         a0+a1*y(i)**2+a2*y(i)**4+a3*y(i)**6+a4*y(i)**8+
!     +         a12*x(i)*y(i)
!        enddo

      return
      end subroutine

!------------------------------------------------
!     average over trajectories
!--------------------------------------------------
      subroutine aver(Ntraj,w,x,y)
        implicit real*8(a-h,o-z)
        real*8 :: x(Ntraj),w(Ntraj)
        y = 0d0

        do i=1,Ntraj 
          y = y+x(i)*w(i) 
        enddo
      
        return
        end subroutine
