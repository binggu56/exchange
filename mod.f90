! --- modules 

      module cdat 

        implicit real*8(a-h,o-z)

        real*8, public, parameter  :: pi = 4.d0*atan(1.d0)

        complex*16, public, parameter :: im=(0d0,1d0)

        integer*4 :: nb ! number of basis in the linear fitting 
      
      save
      contains 

      double precision function aver_traj(ntraj,w,x)

        integer*4, intent(in) :: Ntraj

        real*8 :: x(Ntraj),w(Ntraj)

        aver_traj = 0d0

        do i=1,Ntraj
          aver_traj = aver_traj + x(i)*w(i)
        enddo

        end function aver_traj

! ----- random number seeds
        subroutine seed(idum,Ndim)

        implicit real*8(a-h,o-z)

        integer*4, intent(in) :: Ndim

        integer*4, intent(out) :: idum(Ndim)

        do i=1,Ndim
          idum(i) = 5 + i
        enddo

        return
        end subroutine

      end module cdat

      module sys 
      
      integer (kind = 4), parameter :: NATOMS = 9, NVBINS = 20000
      
      real*8, parameter :: Req = 1.4d0, RATIO = 1.05
      
      integer*4, parameter :: NIP = 2, NPAIRS = NIP*NATOMS/2


      real (kind = 8) :: v(2, NVBINS), xtal(NATOMS, 3)
      real (kind = 8) :: xlen, ylen, zlen, dxmax, dymax, dzmax

      integer (kind = 4) :: vpvec(3, NPAIRS),ivpair(2, NPAIRS)
      integer (kind = 4) :: ipairs(NIP, NATOMS), npair(NATOMS)
      integer (kind = 4) :: nvpair
      integer (kind = 4) :: iwork(60)


      save 
      end module 
