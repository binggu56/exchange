        
	program main
	implicit none
	integer :: idum,i
!	integer*4, parameter :: Ntraj=10
!	real*8 :: x(Ntraj),y(Ntraj)
!        x(1)=0
!        y(1)=0
!        do i=1,Ntraj-1
	
!       x(i+1)=x(i)+1d0/Ntraj
!        y(i+1)=y(i)+1d0/Ntraj
	
        real :: gasdev,random(100) 
        idum=200000
        do i=1,100
        random(i)=gasdev(idum)
        write(*,*) i,random(i)        
        end do




	end program main
