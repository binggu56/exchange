c     QSATS version 1.0 (3 March 2011)

c     file name: vinit.f

c ----------------------------------------------------------------------

c     this subroutine sets up the arrays for linear interpolation of
c     the potential energy function (as a function of the
c     interatomic distance).

c     Morse potential 

c ----------------------------------------------------------------------

      subroutine vinit(rmin, bin)
      
      use cdat

      implicit double precision (a-h, o-z)

      parameter (one=1.0d0)

c --- hartree to kelvin conversion factor

      parameter (hart=315774.65d0)

!      include 'sizes.h'

      include 'qsats.h'

!      r2min=9.0d0 ! Hinde  
      rmin = 0.2d0
      rmax = 10d0   

!      r2min = 4.d0 ! customize 

      bin=(rmax-rmin)/NVBINS   

!      write (6, 6000) sqrt(r2min), bin
!6000  format ('DEFINING potential energy grid parameters'//,
!     +        'minimum R      = ', f10.5, ' bohr'/,
!     +        'R**2 bin size  = ', f10.5, ' bohr**2'/)
!
!      write (6, 6020)
!6020  format ('using HFD-B(He) potential energy curve'/, '[R.A. ',
!     +        'Aziz et al., Mol. Phys. vol. 61, p. 1487 (1987)]'/)

c --- evaluate the potential energy function at grid points in R.

      do i=1, NVBINS
         r = dble(i-1)*bin + rmin
         call ho(r,v0,dv)
         v(1,i) = v0
         v(2,i) = dv 
      enddo

c --- finite difference first-order derivative  

!      do i=1, NVBINS-1
!         v(2, i)=(v(1,i+1)-v(1,i))/bin
!      end do

! --- second order derivative 
        
!      do i=2,NVBINS-1
!         v(3, i)=(v(1, i+1)+v(1, i-1)-2d0*v(1,i))/bin**2
!      end do

c --- debugging output.

      if (idebug.ge.3) then

         vmin=v(1, 1)

         do i=1, NVBINS
            if (v(1, i).lt.vmin) then
               vmin=v(1, i)
               r=dble(i-1)*bin+rmin
            end if
         end do

         vmin=vmin*hart

         write (6, 6100) vmin, r
6100     format ('minimum is ', f12.5, ' K at R = ', f10.5, ' bohr'/)

      end if

      return
      end
! ------------------------
! ----- the following part will be replaced by a simple Morse potential 
      
      subroutine morse(r,v,dv)

      implicit real*8 (a-h, o-z)
      
      d =  0.17429d0
      r0 = 1.4d0 
      a = 1.0435d0
      
      v = d*(exp(-a*(r-r0))-1d0)**2 
      dv = 2d0*d*(-a)*exp(-a*(r-r0))*(exp(-a*(r-r0))-1d0)

      return 
      end subroutine 

      subroutine ho(r,v,dv) 
      implicit real*8 (a-h,o-z) 

      am = 925d0 
      w = 0.0199
      ak = am*w**2 
      r0 = 1.4d0 

      v = ak/2d0*(r-r0)**2
      dv = ak*(r-r0) 
      return 
      end subroutine 

c ----------------------------------------------------------------------

c     this function computes the He-He potential energy as a function
c     of the squared interatomic distance.

c     the potential energy function is described in detail by R.A. Aziz
c     et al., Molecular Physics, volume 61, p. 1487 (1987).

c ----------------------------------------------------------------------

      double precision function hfdbhe(r2)

      implicit real*8 (a-h, o-z)

c --- parameters for the HFD(B) He-He potential energy function

      parameter (astar=1.8443101d5)
      parameter (alstar=10.43329537d0)
      parameter (bestar=-2.27965105d0)
      parameter (d=1.4826d0)
      parameter (c6=1.36745214d0)
      parameter (c8=0.42123807d0)
      parameter (c10=0.17473318d0)

      parameter (rm=5.59926d0)
      parameter (eps=10.948d0)

c --- hartree to kelvin conversion factor

      parameter (hart=315774.65d0)

      r=sqrt(r2)

      x=r/rm

      vstar=astar*exp(-alstar*x+bestar*x**2)

      vd=c6/x**6+c8/x**8+c10/x**10

      if (x.lt.d) vd=vd*exp(-(d/x-1.0d0)**2)

      hfdbhe=(vstar-vd)*eps/hart

      return
      end
