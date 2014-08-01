        subroutine derivs(m1,m2,k1,k2,x,y,Ntraj,dvdx,dvdy,pe)
        implicit none
        integer*8,intent(IN) :: Ntraj
        real*8,intent(IN) :: x(Ntraj),y(Ntraj)
        real*8,intent(IN) :: m1,m2,k1,k2
        integer*4 :: i
        real*8 :: De,a
        real*8,intent(OUT) :: dvdx(Ntraj),dvdy(Ntraj),pe(Ntraj)
        !common/pot_par/m1,m2
!       morse potential for x coordinate, Eckart potential for y
!       direction  
        De=7d0
        a=0.08d0

        do i=1,Ntraj
        dvdx(i)=4d0*De*(1-exp(-a*x(i)**2))*exp(-a*x(i)**2)*a*x(i)
        dvdy(i)=-43.5968*sinh(1.3624*y(i))/(cosh(1.3624*y(i)))**3
        pe(i)=De*(1-exp(-a*x(i)**2))+16d0/(cosh(1.3624*y(i)))**2
        end do
        return
        end subroutine
