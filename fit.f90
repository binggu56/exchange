!     polynomial fitting for {p,r}
!     two steps: 1st global linear fitting; 2nd local fitting with basis up to third order 
      subroutine fit(ntraj,ndim,cp,cr,s1,am,w)

      implicit real*8(a-h, o-z)

      include 'mpif.h'

      integer*4,intent(in)    :: ndim
      real*8, intent(in), dimension(ntraj)  :: w
      real*8, intent(in), dimension(ndim) :: am
!      real*8, intent(in), dimension(ndim,ntraj) :: x,p,r

      real*8, intent(in) :: s1(ndim+1,ndim+1)

!      real*8, intent(out), dimension(ndim,ntraj) :: du,fr,ap
!      real*8, intent(out), dimension(ntraj)  :: u

!      real*8, dimension(ndim,ndim) :: dp,dr
!      real*8, dimension(ndim)      :: ddp,ddr
      real*8, intent(inout) :: cp(ndim+1,ndim),cr(ndim+1,ndim)

      real*8 :: f(ndim+1),c(ndim+1,2*ndim)

      integer info

!---------matrix S1 & c for linear fitting------------------

      nb = ndim+1
      
!      do i = 1,ntraj
!        call basis(ndim,ntraj,i,x,f)
!        do k = 1,ndim
!          do j = 1,ndim+1
!            cp(j,k) = cp(j,k)+f(j)*p(k,i)*w(i) 
!            cr(j,k) = cr(j,k)+f(j)*r(k,i)*w(i)
!          enddo
!        enddo
!      enddo 

      do k=1,ndim
        do j=1,nb
          c(j,k)      = cp(j,k)
          c(j,ndim+k) = cr(j,k)
        enddo
      enddo

! -----matrix S=f*f---------------------------
!      do i=1,ntraj
!        call basis(ndim,ntraj,i,x,f)
!        do k2=1,ndim+1
!          do k1=1,k2
!            s1(k1,k2) = s1(k1,k2)+f(k1)*f(k2)*w(i)
!          enddo
!        enddo
!      enddo
!
!      do k2=1,ndim+1
!        do k1=1,k2
!          s1(k2,k1) = s1(k1,k2)
!        enddo
!      enddo

!---------------solve matrix M*c = f*p-------------
!      time = mpi_wtime()      
      call dposv('U',nb,2*ndim,s1,nb,c,nb,INFO)
      if(info/=0) then
        write(*,*) 'linear fitting of p & r failed.'
        stop
      endif

      do k=1,ndim
        do j=1,nb
          cp(j,k) = c(j,k)
          cr(j,k) = c(j,k+ndim)
        enddo
      enddo
      
!      time = mpi_wtime() -time
!      write(*,*) 'time for linear fitting', time
      return
      end subroutine 


!-------------same for r--------------------------
!      mat = s1
!      call dposv('U',nb,ndim,s1,nb,cr,nb,INFO)
!      if(info/=0) then
!        write(*,*) 'linear fitting of r failed.'
!        stop
!      endif

!-------------output for {ap,ar}--------------------
!      time = mpi_wtime()
!      dp = 0d0
!      dr = 0d0
!      ddp = 0.d0
!      ddr = 0.d0

!------first-order derivative of p,r, second order is 0 for linear fitting---------------
!        do j=1,ndim
!          do k=1,ndim
!            dp(k,j) = cp(k,j)
!            dr(k,j) = cr(k,j)
!          enddo
!        enddo
!      time = mpi_wtime() - time
!      write(*,*) 'time to get approximated p,r',time

!---------------------------------------
!     get approximated p,r; slave
!---------------------------------------
      subroutine aver(ndim,ntrp,ntraj,xp,cp,cr,app,arp)
      
      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: ndim,ntrp,ntraj
      real*8 app(ndim,ntrp),arp(ndim,ntrp),cp(ndim+1,ndim),cr(ndim+1,ndim), &
             xp(ndim,ntrp),f(ndim+1)

      nb = ndim+1

!      cp2 = 0d0
!      cr2 = 0d0
      app = 0d0
      arp = 0d0
      
      do i=1,ntrp
        call basis(ndim,ntrp,i,xp,f)
        do j=1,ndim
          do k=1,nb
            app(j,i) = app(j,i)+cp(k,j)*f(k)
            arp(j,i) = arp(j,i)+cr(k,j)*f(k)
          enddo
        enddo
      enddo

      return
      end subroutine 

! --- collect averages for second-step fitting
!     for every processor 

      subroutine aver_proc(ndim,ntraj_proc,wp,cp2_proc,cr2_proc,s2p,x_proc,p_proc,rp_proc, & 
                           ap_proc,ar_proc) 
      
      implicit real*8(a-h, o-z)

      real*8 :: cp2_proc(4,ndim),cr2_proc(4,ndim),f(4),wp(ntraj_proc)

      real*8 :: ap_proc(ndim,ntraj_proc), ar_proc(ndim,ntraj_proc), & 
                x_proc(ndim,ntraj_proc), p_proc(ndim,ntraj_proc), & 
                rp_proc(ndim,ntraj_proc)

      real*8, intent(out) :: s2p(16,ndim)

      cp2_proc = 0d0
      cr2_proc = 0d0 
      s2p = 0d0
 
      dimloop:do j=1,ndim

        do i=1,ntraj_proc
          f = (/x_proc(j,i),x_proc(j,i)**2,x_proc(j,i)**3,1d0/)
          do k=1,4
            cp2_proc(k,j) = cp2_proc(k,j)+(p_proc(j,i)-ap_proc(j,i))*f(k)*wp(i)
            cr2_proc(k,j) = cr2_proc(k,j)+(rp_proc(j,i)-ar_proc(j,i))*f(k)*wp(i)
          enddo

! ------- would use the symmetry of s2p to reduce the length of
!         s2p(ndim,nb**2) to (ndim,np**2/2), now just use 4**2

          do m=1,4
            do n=1,4
              nn = 4*(m-1)+n
              s2p(nn,j) = s2p(nn,j)+f(m)*f(n)*wp(i) 
            enddo
          enddo

        enddo 

      enddo dimloop
      
      return
      end subroutine

! --- second step fitting, root 

      subroutine fit2(ndim,ntraj,w,cp2,cr2,s2_sum)
      
      implicit real*8(a-h,o-z)

      real*8 :: f2(4),s2_sum(16,ndim),cpr(4,2),w(ntraj),s2(4,4)
             
      real*8, intent(inout) :: cp2(4,ndim),cr2(4,ndim)

      dimloop:do j=1,ndim
        
        do m=1,4
          do n=1,4
            l = 4*(m-1)+n
            s2(m,n) = s2_sum(l,j)
          enddo 
        enddo 

        do i=1,4
          cpr(i,1) = cp2(i,j)
          cpr(i,2) = cr2(i,j)
        enddo 

!        do i=1,ntraj
!          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
!          df2 = (/1d0,2d0*x(j,i),3d0*x(j,i)**2,0d0/)
!        enddo

!        cpr = 0d0
!
!        do i=1,ntraj
!          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
!          do k=1,4
!            cpr(k,1) = cpr(k,1)+(p(j,i)-ap(j,i))*f2(k)*w(i)
!            cpr(k,2) = cpr(k,2)+(r(j,i)-ar(j,i))*f2(k)*w(i)
!          enddo
!        enddo
!          
!        s2 = 0d0`
!        do i=1,ntraj
!          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
!          do m=1,4
!            do n=1,m
!              s2(m,n) = s2(m,n)+f2(m)*f2(n)*w(i) 
!            enddo
!          enddo
!        enddo
!
!        do m=2,4
!          do n=1,m-1
!            s2(n,m) = s2(m,n) 
!          enddo
!        enddo

!       write(*,*) 's2 = ' 
!       do l=1,4
!         write(*,*) (s2(k,l),k=1,4)
!       enddo 
!        
        call dposv('U',4,2,s2,4,cpr,4,info)

        if(info/=0) then
          write(*,*) 8899
8899      format('cubic fitting of r failed')
          stop
        endif

!-------store coefficients---------------------------

        do i=1,4
          cp2(i,j) = cpr(i,1)
          cr2(i,j) = cpr(i,2)
        enddo

      enddo dimloop

!      time = mpi_wtime()-time

!      write(*,*) 'time to do 2nd step fitting',time

      return 
      end subroutine 

! -----------------------------------------------
! --- second step fitting, all nodes
! --------------------------------------------
      subroutine comp(am,wp,ntrajp,ndim,eup,cp,cr,cp2,cr2,xp,app,du,fr)

      use cdat

      implicit real*8(a-h,o-z)

      real*8, dimension(ndim,ndim) :: dp,dr

      real*8, intent(in) :: cp(ndim+1,ndim),cr(ndim+1,ndim),cp2(4,ndim),cr2(4,ndim)

      real*8, intent(in) :: xp(ndim,ntrajp),am(ndim),wp(ntrajp)

      real*8 :: up(ntrajp),ddp(ndim),ddr(ndim),          &
                app(ndim,ntrajp),arp(ndim,ntrajp),du(ndim,ntrajp),fr(ndim,ntrajp)

      real*8 :: f2(4),df2(4),ddf2(4),f(ndim+1)
      
      nb = ndim+1
      app = 0d0
      arp = 0d0

      do i=1,ntrajp
        call basis(ndim,ntrajp,i,xp,f)
        do j=1,ndim
          do k=1,nb
            app(j,i) = app(j,i)+cp(k,j)*f(k)
            arp(j,i) = arp(j,i)+cr(k,j)*f(k)
          enddo
        enddo
      enddo

!------second order derivative only contain diagonal elements--------------------
!-------4-dim array ddp(ndim,ndim,ndim,ntraj) contract to ddp(ndim,ndim,ntraj)----
      do i=1,ntrajp
        f2 = (/xp(j,i),xp(j,i)**2,xp(j,i)**3,1d0/)

        do j=1,ndim
          do k=1,4
            app(j,i) = app(j,i)+cp2(k,j)*f2(k)
            arp(j,i) = arp(j,i)+cr2(k,j)*f2(k)
          enddo
        enddo
      enddo

      up = 0d0
      du = 0d0
      fr = 0d0

!      do i=1,ntrajp
!        do j=1,ndim
!          du(j,i) = 0.d0
!          fr(j,i) = 0.d0
!        enddo
!      enddo

      do j=1,ndim
        do k=1,ndim
          dr(k,j) = cr(k,j)
          dp(k,j) = cp(k,j)
        enddo
      enddo

      traj: do i=1,ntrajp

        do j=1,ndim
            dr(j,j) = cr(j,j)+cr2(1,j)+2d0*cr2(2,j)*xp(j,i)+3d0*cr2(3,j)*xp(j,i)**2
            dp(j,j) = cp(j,j)+cp2(1,j)+2d0*cp2(2,j)*xp(j,i)+3d0*cp2(3,j)*xp(j,i)**2
            ddr(j) = 2d0*cr2(2,j)+6d0*cr2(3,j)*xp(j,i)
            ddp(j) = 2d0*cp2(2,j)+6d0*cp2(3,j)*xp(j,i)
        enddo

        do j=1,ndim
          up(i) = up(i)-(arp(j,i)**2+dr(j,j))/2d0/am(j)
        enddo
        
        do j=1,ndim
          do k=1,ndim
            du(j,i) = du(j,i)-(arp(k,i)*dr(j,k))/am(k)
            fr(j,i) = fr(j,i)-dp(j,k)*arp(k,i)/am(k)
          enddo
        enddo

        do j=1,ndim
          du(j,i) = du(j,i)-ddr(j)/2d0/am(j)
          fr(j,i) = fr(j,i)-ddp(j)/2d0/am(j)
        enddo

      enddo traj

      eup = 0d0
      do i=1,ntrajp
        eup = eup + wp(i)*up(i)
      enddo

      return
      end subroutine
!----------------------------------------
!     linear basis 
!---------------------------------------
      subroutine basis(Ndim,Ntr,i,x,f)

      implicit real*8(a-h,o-z)

      integer*4,intent(in)  :: Ndim,Ntr,i

      real*8      :: x(Ndim,Ntr)
      real*8,intent(OUT) :: f(Ndim+1)

      f = 0d0
!---basis vector f = ((x(i),i=1,Ndim),1) for each trajectory-------------------
      do j=1,Ndim
        f(j)=x(j,i)
      enddo
      f(Ndim+1)=1d0
     
      return
      end subroutine
!--------------------------------------------------------------

