      subroutine pdf(ndim,ntraj_proc,NBIN_PDF,wp,x_proc,p_gr)

      use cdat, only: pi, bin_pdf, binvrs_pdf

      implicit real*8(a-h, o-z)

      include 'sizes.h'

      include 'qsats.h'

      real*8 :: x(ndim,ntraj_proc), wp(ntraj_proc) 

      real*8, intent(out) :: p_gr(NBIN_PDF)

!      common / bincom/ bin, binvrs, r2min 

      dimension q(NATOM3)

!      bin_pdf = 1d-2 
!      binvrs_pdf = 1d0/bin 

      p_gr = 0d0

      traj: do k=1,ntraj_proc

        do j=1,NATOM3
          q(j) = x_proc(j,k)
        enddo 

        do n=1, nvpair2

          i=ivpair2(1, n)
          j=ivpair2(2, n)

          dx=-((q(3*j-2))+vpvec2(1, n)+(-q(3*i-2)))
          dy=-((q(3*j-1))+vpvec2(2, n)+(-q(3*i-1)))
          dz=-((q(3*j))  +vpvec2(3, n)+(-q(3*i))  )

          r2=dx*dx+dy*dy+dz*dz
          r = dsqrt(r2)

!          ibin=int((r2-r2min)*binvrs)+1
          ibin = int(r*binvrs_pdf)

          if(ibin < NBIN_PDF) p_gr(ibin)=p_gr(ibin)+wp(k)/r2
                 
        end do

      enddo traj 

      return 
      end subroutine
