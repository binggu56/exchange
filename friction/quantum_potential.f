        subroutine quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)
        integer*4,intent(IN) :: Ntraj
        real*8, intent(OUT) :: qfx(Ntraj),qfy(Ntraj),qp(Ntraj)
        integer*4 :: i,j
	integer INFO
        real*8,intent(IN) :: m1,m2
        real*8,intent(IN) :: w(Ntraj),x(Ntraj),y(Ntraj)
        real*8 :: f(3),rx(Ntraj),ry(Ntraj)
	real*8 ::  s(3,3),c(3,2)
        !common/pot_par/m1,m2
! define c matrix 
        c(1,1)=-0.5d0
        c(1,2)=0d0
        c(2,1)=0d0
        c(2,2)=-0.5d0
        c(3,1)=0d0
        c(3,2)=0d0

! quantum force
        s=0d0
        do i=1,Ntraj
        f=(/x(i),y(i),1d0/)
        do m=1,3
        do n=1,3
        s(m,n)=w(i)*f(m)*f(n)+s(m,n)
        end do
        end do
        end do

        !do m=1,3
        !print *, s(m,1),s(m,2),s(m,3)
        !enddo
        !STOP
        
! calculate matrix c(t)
        call DPOSV('U',3,2,s,3,c,3,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

! the momentum operator r=cf
        do i=1,Ntraj
        rx(i)=c(1,1)*x(i)+c(2,1)*y(i)+c(3,1)
        ry(i)=c(1,2)*x(i)+c(2,2)*y(i)+c(3,2)

! calculate quantum potential
        qp(i)=-rx(i)**2/(2d0*m1)-ry(i)**2/(2d0*m2)
     &     -c(1,1)/(2d0*m1)-c(2,2)/(2d0*m2)
        
        qfx(i)=rx(i)*c(1,1)/m1+ry(i)*c(1,2)/m2
        qfy(i)=rx(i)*c(2,1)/m1+ry(i)*c(2,2)/m2

        end do

        return
        end subroutine
