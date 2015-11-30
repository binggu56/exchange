!     polynomial fitting for {p,r}
!     two steps: 1st global linear fitting; 2nd local fitting with basis up to third order 
      subroutine fit(am,x,p,r,w,ndim,ntr,u,du,ap,fr)

      implicit real*8(a-h, o-z)

      integer*4                  :: ntr,ndim
!      integer*4, parameter        :: nb = ndim+1
      real*8, dimension(ntr)      :: w,u
      real*8, dimension(ndim)     :: am
      real*8, dimension(ndim,ntr) :: x,p,r,du,fr,ap,ar
      real*8, dimension(ndim,ndim,ntr) :: ddp,ddr,dp,dr

      real*8 f(ndim+1,ntr),f2(4),df2(4),ddf2(4),s1(ndim+1,ndim+1),s2(4,4),cp(ndim+1,ndim), & 
             mat(ndim+1,ndim+1),cr(ndim+1,ndim),cpr(4,2)
      integer info

!---------matrix S1 & c for linear fitting------------------
      nb = ndim+1
      cp = 0d0
      cr = 0d0
      s1 = 0d0
      err = 0d0

      call basis(Ndim,Ntr,x,f)
      
      print *,'tag1'    

      do j = 1,ndim+1
        do k = 1,ndim
          do i = 1,ntr
            print *,'here'
            cp(j,k) = cp(j,k)+f(j,i)*p(k,i)*w(i) 
            cr(j,k) = cr(j,k)+f(j,i)*r(k,i)*w(i)
          enddo
        enddo
      enddo 

      print *,'tag2'

      do k1=1,nb
        do k2=1,nb
          do i=1,ntr
            s1(k1,k2) = s1(k1,k2)+f(k1,i)*f(k2,i)*w(i)
          enddo
        enddo
      enddo

      print *,'tag3'
!---------------solve matrix M*c = f*p-------------
      mat = s1
      call dposv('U',nb,ndim,mat,nb,cp,nb,INFO)
      if(info/=0) then
        write(*,*) 'linear fitting of p failed.'
        stop
      endif
!-------------same for r--------------------------
      mat = s1
      call dposv('U',nb,ndim,mat,nb,cr,nb,INFO)
      if(info/=0) then
        write(*,*) 'linear fitting of r failed.'
        stop
      endif

      print *,'tag4'
!-------------output for {ap,ar}--------------------
      dp = 0d0
      ap = 0d0
      ar = 0d0
      dr = 0d0
      ddp = 0.d0
      ddr = 0.d0

      do i=1,ntr
        do j=1,ndim
          do k=1,nb
            ap(j,i) = ap(j,i)+cp(k,j)*f(k,i)
            ar(j,i) = ar(j,i)+cr(k,j)*f(k,i)
          enddo
        enddo
      enddo
!------first-order derivative of p,r, second order is 0 for linear fitting---------------
      do i=1,ntr
        do j=1,ndim
          do k=1,ndim
            dp(k,j,i) = cp(k,j)
            dr(k,j,i) = cr(k,j)
          enddo
        enddo
      enddo


!---------------------------------------
!     2nd fitting for each DoF      
!---------------------------------------
      dimloop:do j=1,ndim

        s2 = 0d0
        cpr = 0d0

!        do i=1,ntr
!          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
!          df2 = (/1d0,2d0*x(j,i),3d0*x(j,i)**2,0d0/)
!        enddo

        do k=1,4
          do i=1,ntr
          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
            cpr(k,1) = cpr(k,1)+(p(j,i)-ap(j,i))*f2(k)*w(i)
            cpr(k,2) = cpr(k,2)+(r(j,i)-ar(j,i))*f2(k)*w(i)
          enddo
        enddo
          
          do m=1,4
            do n=1,4
              do i=1,ntr
              f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
              s2(m,n) = s2(m,n)+f2(m)*f2(n)*w(i) 
            enddo
          enddo
        enddo 
        
        call dposv('U',4,2,s2,4,cpr,4,info)
        if(info/=0) then
          write(*,*) 'cubic fitting of r failed.'
          stop
        endif
      

!------second order derivative only contain diagonal elements--------------------
!-------4-dim array ddp(ndim,ndim,ndim,ntraj) contract to ddp(ndim,ndim,ntraj)----
      do i=1,ntr
          f2 = (/x(j,i),x(j,i)**2,x(j,i)**3,1d0/)
          df2 = (/1d0,2d0*x(j,i),3d0*x(j,i)**2,0d0/)
          ddf2 = (/0d0,2d0,6d0*x(j,i),0d0/)

        do k=1,4
          ap(j,i) = ap(j,i)+cpr(k,1)*f2(k)
          dp(j,j,i) = dp(j,j,i)+cpr(k,1)*df2(k)
          ddp(j,j,i) = ddp(j,j,i)+cpr(k,1)*ddf2(k)

          ar(j,i) = ar(j,i)+cpr(k,2)*f2(k)
          dr(j,j,i) = dr(j,j,i)+cpr(k,2)*df2(k)
          ddr(j,j,i) = ddr(j,j,i)+cpr(k,2)*ddf2(k)
        enddo

      enddo

      enddo dimloop


      u = 0.d0
      du = 0.d0
      fr = 0.d0


      do i=1,ntr
        do j=1,ndim
          u(i) = u(i)-(ar(j,i)**2+dr(j,j,i))/2d0/am(j)
        enddo
        
        do j=1,ndim

          do k=1,ndim
            du(j,i) = du(j,i)-(ar(k,i)*dr(j,k,i))/am(k)
            fr(j,i) = fr(j,i)-dp(j,k,i)*ar(k,i)/am(k)
          enddo
      
          du(j,i) = du(j,i)-ddr(j,j,i)/2d0/am(j)
          fr(j,i) = fr(j,i)-ddp(j,j,i)/2d0/am(j)
        enddo

      enddo


      return
      end subroutine
!----------------------------------------
!     linear basis 
!---------------------------------------
      subroutine basis(Ndim,Ntr,x,f)
      implicit real*8(a-h,o-z)
      integer*4  :: Ndim,Ntr
      real*8      :: x(Ndim,Ntr)
      real*8,intent(OUT) :: f(Ndim+1,Ntr)

      f = 0d0

!---basis vector f = ((x(i),i=1,Ndim),1) for each trajectory-------------------
      do i=1,Ntr
        do j=1,Ndim
          f(j,i)=x(j,i)
        enddo
          f(Ndim+1,i)=1d0
       enddo
     
      return
      end subroutine
!--------------------------------------------------------------
