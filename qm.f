!     QSATS version 1.0 (3 March 2011)

!     file name: eloc.f

! ----------------------------------------------------------------------

!     this computes the total energy and the expectation value of the
!     potential energy from the snapshots recorded by QSATS.

! ----------------------------------------------------------------------


! --------main program 

      program eloc

      use cdat
      use sys 

      implicit double precision (a-h, o-z)

!      include 'sizes.h'

!      include 'qsats.h'

      include 'mpif.h'
      
      common /bincom/ bin, binvrs, rmin

!------------trajectory initial allocate------------------------

      real*8, dimension(NATOMS)         :: p0,alpha,x0,cf,am 

      real*8, allocatable, dimension(:)   :: w,pe,ke,wp

      real*8, allocatable, dimension(:,:) :: x,p,ap,ar,rp,du,fr

      real*8, allocatable, dimension(:,:) :: ap_proc,du_proc,fr_proc,
     +        x_proc,p_proc,rp_proc,s1,cpp,crp,s1p,ar_proc,s2p,
     +        cp2,cr2, s2_sum,mat,mat_proc 
     
      real*8, allocatable, dimension(:,:) :: cp,cr,cp2_proc,cr2_proc,cpr
      
      complex*16, allocatable, dimension(:,:) :: g
      
      complex*16, allocatable, dimension(:,:) :: cor_proc, cor 

      character*4 :: atom(NATOMS)

      integer     :: myid, ierr, numprocs, root, tag
!------------------------------------------------------

! --- this common block is used to enable interpolation in the potential
!     energy lookup table in the subroutine local below.


!      dimension vtavg(NREPS), vtavg2(NREPS), etavg(NREPS), etavg2(NREPS)

      parameter (half=0.5d0)
      parameter (one=1.0d0)

      tag = 1

! --- initialize mpi environment
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world, myid, ierr)
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      
      root = 0

!---- read variables from IN
      open(10,file='IN', status='old', action='read')

      read(10,*) ntraj
      read(10,*) kmax,dt
      read(10,*) am0
      read(10,*) a0
      read(10,*) cf0
      read(10,*) iread

      close(10)

! --- initialization.

!      call tstamp

      if (myid == 0) then 
      
        open(103, file='en.dat',    action='write')
        open(112, file='traj.dat',  action='write')
        open(110, file='cor1.dat',  action='write')
        open(111, file='cor2.dat',  action='write')
        open(116, file='cor3.dat',  action='write')
        open(114, file='rfit1.dat', action='write')
        open(115, file='rfit2.dat', action='write')

!        write (6, 6001) NREPS, NATOMS, NATOM3, NATOM6, NATOM7,
!     +                NVBINS, RATIO, NIP, NPAIRS
!6001    format ('compile-time parameters:'//,  
!     +       'NREPS  = ', i6/,               
!     +       'NATOMS = ', i6/,               
!     +       'NATOM3 = ', i6/,               
!     +       'NATOM6 = ', i6/,               
!     +       'NATOM7 = ', i6/,               
!     +       'NVBINS = ', i6/,               
!     +       'RATIO  = ', f8.4/,              
!     +       'NIP    = ', i6/,               
!     +       'NPAIRS = ', i6/)    
      
        write (6, 6001) NATOMS,NVBINS, RATIO
6001    format ('compile-time parameters:'//,  
     +       'NATOMS = ', i6/,               
     +       'NVBINS = ', i6/,               
     +       'RATIO  = ', f8.4/)               
      endif   ! master work end


! --- passing input variables to other procs

!      call mpi_bcast(ntraj,1,mpi_integer,root,mpi_comm_world,ierr)
!
!      call mpi_bcast(dt,1,mpi_double_precision,root,
!     +               mpi_comm_world,ierr)
!
!      call mpi_bcast(a0,1,mpi_double_precision,root,
!     +               mpi_comm_world,ierr)
!
!      call mpi_bcast(am0,1,mpi_double_precision,root,
!     +               mpi_comm_world,ierr)
!
!      call mpi_bcast(cf0,1,mpi_double_precision,root,
!     +               mpi_comm_world,ierr)
!
!      call mpi_bcast(kmax,1,mpi_integer,root,mpi_comm_world,ierr)

      ndim = NATOMS ! linear chain model 

      call init_pars(ndim)

!      den = 5.231d-3
!      den = 4.61421d-3 ! Hinde 

      ntraj_proc = ntraj/numprocs

      call vinit(rmin, bin)

      if(myid == root) then

      write (*, 6000) rmin, bin

!6000  format ('DEFINING potential energy grid parameters'//,
!     +        'minimum R      = ', f10.5, ' bohr'/,
!     +        'R**2 bin size  = ', f10.5, ' bohr**2'/)

6000  format ('DEFINING potential energy grid parameters'//,
     +        'minimum R      = ', f10.5, ' bohr'/,
     +        'R bin size  = ', f10.5, ' bohr'/)

!      write (*, 6020)
!6020  format ('using HFD-B(He) potential energy curve'/, '[R.A. ',
!     +        'Aziz et al., Mol. Phys. vol. 61, p. 1487 (1987)]'/)

      endif   ! root print
      
      binvrs=one/bin

! --- read crystal lattice points.

!      ltfile = 'lattice-file-180'

!      open (8, file=ltfile, status='old')

!      read (8, *) nlpts  ! number of atoms in lattice file 
      
!      read (8, *) xlen, ylen, zlen ! --- read the edge lengths of the supercell.

!      if(myid == root) then
      
!        write (6, 6200) ltfile
!6200    format ('READING crystal lattice from ', a16/)

!        if (nlpts.ne.NATOMS) then
!          write (6, *) 'ERROR: number of atoms in lattice file = ', 
!     +                 nlpts
!          write (6, *) 'number of atoms in source code = ', NATOMS
!          stop
!        endif

      xlen = NATOMS*Req 
!      endif

! --- bcast nlpts and xlen,ylen,zlen
!      call MPI_BCAST(nlpts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(xlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(ylen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(zlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! --- compute a distance scaling factor.

!!      den0=dble(NATOMS)/(xlen*ylen*zlen)

! --- scale is a distance scaling factor, computed from the atomic
!     number density specified by the user.

!!      scale=exp(dlog(den/den0)/3.0d0)
      scale = 1d0 

!      xlen=xlen/scale
!      ylen=ylen/scale
!      zlen=zlen/scale


      dxmax=half*xlen
!      dymax=half*ylen
!      dzmax=half*zlen

      do i=1, NATOMS

!         read (8, *) xtal(i, 1), xtal(i, 2), xtal(i, 3)
         xtal(i,1) = half*Req+(i-1)*Req 

!         xtal(i, 1)=xtal(i, 1)/scale
!         xtal(i, 2)=xtal(i, 2)/scale
!         xtal(i, 3)=xtal(i, 3)/scale
      end do

!      close (8)
      
      if (myid == 0) then

!      write (6, 6300) scale
!6300  format ('supercell scaling factor computed from density = ',f12.8/)

      write (6, 6310) xlen
6310  format ('supercell edge lengths [bohr]         = ', f10.5/)

!      write (6, 6320) xtal(NATOMS, 1), xtal(NATOMS, 2),
!     +                xtal(NATOMS, 3)
!6320  format ('final lattice point [bohr]            = ', 3f10.5/)
      
      endif

! --- this variable helps us remember the nearest-neighbor distance.

!      rnnmin=-1.0d0
      rnnmin = Req  


!      do j=2, NATOMS
!
!         dx=xtal(j, 1)-xtal(1, 1)
!!         dy=xtal(j, 2)-xtal(1, 2)
!!         dz=xtal(j, 3)-xtal(1, 3)
!
!! ------ this sequence of if-then-else statements enforces the
!!        minimum image convention.
!
!         if (dx.gt.dxmax) then
!            dx=dx-xlen
!         else if (dx.lt.-dxmax) then
!            dx=dx+xlen
!         end if
!
!         if (dy.gt.dymax) then
!            dy=dy-ylen
!         else if (dy.lt.-dymax) then
!            dy=dy+ylen
!         end if
!
!         if (dz.gt.dzmax) then
!            dz=dz-zlen
!         else if (dz.lt.-dzmax) then
!            dz=dz+zlen
!         end if
!
!         r=sqrt(dx*dx+dy*dy+dz*dz)
!
!         if (r.lt.rnnmin.or.rnnmin.le.0.0d0) rnnmin=r
!
!      end do

      if(myid == 0) then

      write (6, 6330) rnnmin
6330  format ('nearest neighbor (NN) distance [bohr] = ', f10.5/)

!      write (6, 6340) xlen/rnnmin, ylen/rnnmin, zlen/rnnmin
!6340  format ('supercell edge lengths [NN distances] = ', 3f10.5/)

      endif

! --- compute interacting pairs.

      do i=1, NATOMS
         npair(i) = 0
      end do

      nvpair=0

      do i=1, NATOMS
      do j=i+1, NATOMS ! eliminate overcounting 

!         if (j.ne.i) then

            dx=xtal(j, 1)-xtal(i, 1)
!            dy=xtal(j, 2)-xtal(i, 2)
!            dz=xtal(j, 3)-xtal(i, 3)

! --------- this sequence of if-then-else statements enforces the
!           minimum image convention.

            if (dx.gt.dxmax) then
               dx=dx-xlen
            else if (dx.lt.-dxmax) then
               dx=dx+xlen
            end if

!            if (dy.gt.dymax) then
!               dy=dy-ylen
!            else if (dy.lt.-dymax) then
!               dy=dy+ylen
!            end if

!            if (dz.gt.dzmax) then
!               dz=dz-zlen
!            else if (dz.lt.-dzmax) then
!               dz=dz+zlen
!            end if

!            r2=dx*dx+dy*dy+dz*dz
            r2 = dx*dx 
            r=sqrt(r2)

! --------- interacting pairs are those for which r is less than a
!           certain cutoff amount. 

            if (r/rnnmin.lt.RATIO) then

               nvpair = nvpair + 1

               ivpair(1, nvpair)=i
               ivpair(2, nvpair)=j

               vpvec(1, nvpair)= dx
!               vpvec(2, nvpair)= dy
!               vpvec(3, nvpair)= dz

               npair(i) = npair(i)+1

               ipairs(npair(i), i)=nvpair

            end if

!         end if

      end do
      end do
      
      if (myid == 0) then
      write (6, 6400) npair(1), nvpair
6400  format ('atom 1 interacts with ', i3, ' other atoms'//, 
     +        'total number of interacting pairs = ', i6/)
      endif

! --- reduce half the number of pairs 

!      call reduce()

! --- allocate local arrays x,p,rp for root only (only those arrays only
!     exist at root proc) 
      if (myid == 0) then

        allocate(s1(ndim+1,ndim+1),pe(ntraj),ke(ntraj),x(ndim,ntraj), 
     +           p(ndim,ntraj),ap(ndim,ntraj),ar(ndim,ntraj),
     +           rp(ndim,ntraj),w(ntraj),s2_sum(16,ndim),mat(nb,nb), 
     +           cor(3,ndim))
      endif

! --- allocate arrays for all cores 

      allocate(cp2(4,ndim),cr2(4,ndim), mat_proc(nb,nb))

      allocate(wp(ntraj_proc),cpr(ndim+1,2*ndim), s2p(16,ndim))

      allocate(fr_proc(ndim,ntraj_proc),du_proc(ndim,ntraj_proc),
     +         ap_proc(ndim,ntraj_proc),x_proc(ndim,ntraj_proc),
     +         p_proc(ndim,ntraj_proc),rp_proc(ndim,ntraj_proc),
     +         ar_proc(ndim,ntraj_proc))

      allocate(s1p(ndim+1,ndim+1),cpp(ndim+1,ndim),crp(ndim+1,ndim),
     +         cp2_proc(4,ndim),cr2_proc(4,ndim),
     +         cp(ndim+1,ndim),cr(ndim+1,ndim))

      allocate(g(nb,3),cor_proc(3,ndim))
      
! --- initial parameters 

      dt2 = dt/2d0 
      t   = 0d0

      do i=1,Ndim
        am(i) = am0*1836.15d0 
        x0(i) = 0.0d0
        p0(i) = 0d0
        alpha(i) = a0
        cf(i) = cf0
      enddo

      if(myid == root) then

! ----- initial setup for trajectories

        write(*,6005) numprocs, ntraj_proc        
6005    format('Num of cores                 =', i6/,
     +         'Num of trajectories per core =', i6)

! ----- check number of trajectories

        if (mod(ntraj,numprocs) .ne. 0) then
          
          write(*,6004) mod(ntraj,numprocs)
6004      format('Bad number of trajectories.'/, 
     +           'The remainder is',i6/)
          stop
        endif

! --- intial print 

      write(*,1010) Ntraj,Ndim,kmax,dt,cf(1),am(1)
1010  format('Initial Conditions'//,
     +       'Ntraj   = ' , i6/,                          
     +       'Ndim    = ' , i6/,                         
     +       'kmax    = ' , i6/,  'dt      = ' , f8.4/,  
     +       'fric    = ' , f8.4/,'Mass    = ' , e14.6/) 

      write(*,1011) alpha(1),x0(1),p0(1)
1011  format('Initial Wavepacket'//,'alpha0   = ' , f8.4/,  
     +       'center   = ' , f8.4/,'momentum = ' , f8.4/)

! --- intial {x,p,r}, sampling, root

      if(iread == 0) then
      
        call init_traj(ndim,ntraj,alpha,x,x0,p,p0,rp)
       
! --- continue job from last checkpoint
      elseif(iread == 1) then

      open(11, file='init.dat', status='old',action='read')

      do i=1,ntraj
        do j=1,ndim
          read(11,*) x(j,i),p(j,i),rp(j,i)
        enddo
      enddo

      close(11)

      endif

!---- expectation value of x(1,ntraj)--------
      do i=1,ntraj 
        w(i) = 1d0/dble(ntraj)
      enddo 
      
!      g = (0d0,0d0)
!      ic = 64
!      g(3*ic-2,1) = (1d0,0d0) 
!      g(3*ic-1,2) = (1d0,0d0)
!      g(3*ic,  3) = (1d0,0d0) 
      
      av = 0d0
      do i=1,ntraj
        av = av+w(i)*x(1,i)
      enddo
      write(*,6008) av
6008  format('Expectation value of x(1,ntraj) = ', f10.6)

      endif ! root work end 

!      time = mpi_wtime()

      call mpi_barrier(mpi_comm_world, ierr) 

! --- send required info {cf,x,p,r,w} for slave nodes 
      call mpi_scatter(x,ndim*ntraj_proc,mpi_double_precision,
     +                 x_proc,ndim*ntraj_proc,mpi_double_precision,
     +                 root,mpi_comm_world,ierr)

      call mpi_scatter(p,ndim*ntraj_proc,mpi_double_precision,
     +                 p_proc,ndim*ntraj_proc,mpi_double_precision,
     +                 root,mpi_comm_world,ierr)

      call mpi_scatter(rp,ndim*ntraj_proc,mpi_double_precision,
     +                 rp_proc,ndim*ntraj_proc,mpi_double_precision,
     +                 root,mpi_comm_world,ierr)

      call mpi_scatter(w,ntraj_proc,mpi_double_precision,
     +                 wp,ntraj_proc,mpi_double_precision,
     +                 root,mpi_comm_world,ierr)

!       call mpi_bcast(cf,ndim,mpi_double_precision,root,
!     +               mpi_comm_world,ierr)

      if(myid == root) then

!        time = mpi_wtime() - time
        write(*,6678)
6678    format('Scatter work finished.') 
!     +         f12.6/)
        
        write(*,6679) 
6679    format('Now propagate quantum trajectories...')

      endif

! --- compute intial force fields
      call prefit(ntraj_proc,ndim,wp,x_proc,p_proc,rp_proc,
     +              s1p,cpp,crp)
     
      call mpi_reduce(s1p,s1,(ndim+1)*(ndim+1),MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(cpp,cp,(ndim+1)*ndim,MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(crp,cr,(ndim+1)*ndim,MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
      if (myid == root) then
!        write(6,*) 's1',s1(1,1),s1(2,3)
        call fit(ntraj,ndim,cp,cr,s1,am,w)
      endif
      
      call mpi_bcast(cp,(ndim+1)*ndim,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)
      call mpi_bcast(cr,(ndim+1)*ndim,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

      call aver(ndim,ntraj_proc,ntraj,x_proc,cp,cr,ap_proc,ar_proc)

      call aver_proc(ndim,ntraj_proc,wp,cp2_proc,cr2_proc,s2p,
     +                  x_proc,p_proc,rp_proc,ap_proc,ar_proc) 

      call mpi_barrier(mpi_comm_world,ierr)
      
      call mpi_reduce(cp2_proc,cp2,4*ndim,
     +     MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

      call mpi_reduce(cr2_proc,cr2,4*ndim,
     +     MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

      call mpi_reduce(s2p,s2_sum,16*ndim,
     +     MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

      if(myid == 0) then
        call fit2(ndim,ntraj,w,cp2,cr2,s2_sum)
      endif
      
      call mpi_bcast(cp2,ndim*4,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

      call mpi_bcast(cr2,ndim*4,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

      call comp(am,wp,ntraj_proc,ndim,eup,cp,cr,cp2,cr2,
     +          x_proc,ap_proc,du_proc,fr_proc)
! --- end of section to compute initial force fields 
      if(myid == 0) then 
        write(6,*) 'Begin time loop for trjectory propagation...'
      endif 

! --- trajectories propagation      

      time: do kt=1,kmax 

!        kout: do it = 1,kout 

        t=t+dt

        call  traj(myid,dt,ndim,ntraj_proc,cf,am,x_proc,p_proc,
     +                rp_proc,ap_proc,wp,du_proc,fr_proc)

!      if(myid == 0) then 
!        write(6,*) 'inside loop traj'
!      endif 
! ----- compute quantum force x(Ndim,Ntraj),p,r --------

        call prefit(ntraj_proc,ndim,wp,x_proc,p_proc,rp_proc,
     +              s1p,cpp,crp)
     
        call mpi_reduce(s1p,s1,(ndim+1)*(ndim+1),MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call mpi_reduce(cpp,cp,(ndim+1)*ndim,MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      
        call mpi_reduce(crp,cr,(ndim+1)*ndim,MPI_DOUBLE_PRECISION,
     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
!        time = mpi_wtime() - time

! ----- linear fitting for the first step
        if (myid == root) then
!         write(6,*) 'after reduing' 
!           write(*,6688) time
!6688       format('time to collect matrix elements',f12.6/)
!          time = mpi_wtime()
          call fit(ntraj,ndim,cp,cr,s1,am,w)
!      write(6,*) 'first fit' 
!          time = mpi_wtime()-time
!          write(*,6680) time
!6680      format('time for linear fit at root', f12.6/)
        endif
        

     
!        call MPI_BCAST(cpr,(ndim+1)*2*ndim,MPI_DOUBLE_PRECISION,root,
!     +       MPI_COMM_WORLD,ierr)

        call mpi_bcast(cp,(ndim+1)*ndim,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

        call mpi_bcast(cr,(ndim+1)*ndim,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

! ----  get approximate {p,r}, salve 
        
!        time = mpi_wtime()

        call aver(ndim,ntraj_proc,ntraj,x_proc,cp,cr,ap_proc,ar_proc)

! ----- collect approximated {ap,ar}
      
!        call MPI_GATHER(ap_proc,ndim*ntraj_proc,mpi_double_precision,ap,
!     +                  ndim*ntraj_proc,mpi_double_precision,ROOT,
!     +                  MPI_COMM_WORLD,ierr)
!
!        call MPI_GATHER(ar_proc,ndim*ntraj_proc,mpi_double_precision,ar,
!     +                  ndim*ntraj_proc,mpi_double_precision,ROOT,
!     +                  MPI_COMM_WORLD,ierr)
!
! ---- compute averages of f*f, f= (1,x,x^2,x^3)

        call aver_proc(ndim,ntraj_proc,wp,cp2_proc,cr2_proc,s2p,
     +                  x_proc,p_proc,rp_proc,ap_proc,ar_proc) 

        call mpi_barrier(mpi_comm_world,ierr)
        
        call mpi_reduce(cp2_proc,cp2,4*ndim,
     +       MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

        call mpi_reduce(cr2_proc,cr2,4*ndim,
     +       MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

        call mpi_reduce(s2p,s2_sum,16*ndim,
     +       MPI_DOUBLE_PRECISION,mpi_sum,root,MPI_COMM_WORLD,ierr)

! ----  do second fitting, root

        if(myid == root) then

!         time = mpi_wtime()-time
          
!          write(*,*) 's2 on root',(s2_sum(i,1),i=1,16)

!          write(*,6689) time

!6689      format('time to gather approximated p,r',f12.6/)

!          time = mpi_wtime()

! ------- most time consuming part 

          call fit2(ndim,ntraj,w,cp2,cr2,s2_sum)
!          write(6,*) 'test'

!          time = mpi_wtime() - time
!          write(*,6691) time
!6691      format('time to do second fit at root',f12.6/)
        endif

!        call MPI_BARRIER(mpi_comm_world,ierr)

        call mpi_bcast(cp2,ndim*4,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

        call mpi_bcast(cr2,ndim*4,mpi_double_precision,root,
     +                   mpi_comm_world,ierr)

! ----- wait until get info 
        call mpi_barrier(mpi_comm_world,ierr)

        
! ----- send quantum force to slave nodes, compute partial average
        
!        call mpi_scatter(ap,ndim*ntraj_proc,mpi_double_precision,
!     +                   ap_proc,ndim*ntraj_proc,mpi_double_precision,
!     +                   root,mpi_comm_world,ierr)
!
!        call mpi_scatter(du,ndim*ntraj_proc,mpi_double_precision,
!     +                   du_proc,ndim*ntraj_proc,mpi_double_precision,
!     +                   root,mpi_comm_world,ierr)
!
!        call mpi_scatter(fr,ndim*ntraj_proc,mpi_double_precision,
!     +                   fr_proc,ndim*ntraj_proc,mpi_double_precision,
!     +                   root,mpi_comm_world,ierr)

!        call mpi_bcast(du,ndim*ntraj,mpi_double_precision,root,
!     +                   mpi_comm_world,ierr)
!
!        call mpi_bcast(fr,ndim*ntraj,mpi_double_precision,root,
!     +                   mpi_comm_world,ierr)


! ------- propagate trajectory in each proc for one time step
        call comp(am,wp,ntraj_proc,ndim,enu_proc,cp,cr,cp2,cr2,
     +            x_proc,ap_proc,du_proc,fr_proc)

        call prop_p(myid,dt,ndim,ntraj_proc,cf,am,x_proc,p_proc,rp_proc,
     +            wp,du_proc,fr_proc,env_proc,enk_proc)

!        call cor_mat(ntraj_proc,ndim,wp,x_proc,p_proc,mat_proc)

!        call mpi_reduce(mat_proc,mat,nb*nb,MPI_DOUBLE_PRECISION,
!     +                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
!        if(myid == 0) then 
!          call prop_c(s1,mat,am,dt,ndim,ntraj,g) 
!        endif

! ----- transfer coefficients to each processor
         
!        call mpi_bcast(g,3*nb,mpi_double_complex,0,mpi_comm_world,ierr)

! ----- compute local correlation C(AB), here A=B=x
!        call corr(am,dt,ndim,ntraj_proc,wp,x_proc,g,cor_proc) 
        
!        call mpi_reduce(cor_proc,cor,3*ndim,mpi_double_complex,mpi_sum,
!     +       0,mpi_comm_world, ierr)     

! ----- end of correlation section ---------------

! ----- set values to 0 to do mpi_reduce to get the full {x(ndim,ntraj),p,r} matrix
! ----- collect info from other nodes {x,p,r}, then compute quantum potential, by root
!        time = mpi_wtime()
!
!        call MPI_GATHER(x_proc,ndim*ntraj_proc,mpi_double_precision,x,
!     +                  ndim*ntraj_proc,mpi_double_precision,ROOT,
!     +                  MPI_COMM_WORLD,ierr)
!
!        call MPI_GATHER(p_proc,ndim*ntraj_proc,mpi_double_precision,p,
!     +                  ndim*ntraj_proc,mpi_double_precision,ROOT,
!     +                  MPI_COMM_WORLD,ierr)
!
!        call MPI_GATHER(rp_proc,ndim*ntraj_proc,mpi_double_precision,rp,
!     +                  ndim*ntraj_proc,mpi_double_precision,ROOT,
!     +                  MPI_COMM_WORLD,ierr)
!
! ----- root, gather energy components, can be eliminated by mpi_i/o 

        call mpi_reduce(env_proc,env,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                  root,MPI_COMM_WORLD,ierr)

        call mpi_reduce(enu_proc,enu,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                  root,MPI_COMM_WORLD,ierr)

        call mpi_reduce(enk_proc,enk,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                  root,MPI_COMM_WORLD,ierr)


! ----- write data to files 
        if(myid == root) then
          
C          time = mpi_wtime() - time 
C          write(*,6690) time
C6690      format('time to gather {x,p,r}',f12.6/)
          write(103,1000) t,enk,env,enu,(env+enk+enu)
          call flush(103)

          write(112,1000) t,(x_proc(1,i),i=1,20),(p_proc(1,i),i=1,20) 
          write(110,1000) t,(cor(1,j),j=1,ndim)
          write(111,1000) t,(cor(2,j),j=1,ndim)
          write(113,1000) t,(cor(3,j),j=1,ndim)

        endif

      enddo time 

      if (myid == root) then 
      write(6,6691) 
6691  format('End of trajectory propagation, now output information.') 
      endif 


! --- record final data , transfer data to output 
      call MPI_GATHER(x_proc,ndim*ntraj_proc,mpi_double_precision,x,
     +                ndim*ntraj_proc,mpi_double_precision,ROOT,
     +                MPI_COMM_WORLD,ierr)

      call MPI_GATHER(p_proc,ndim*ntraj_proc,mpi_double_precision,p,
     +                ndim*ntraj_proc,mpi_double_precision,ROOT,
     +                MPI_COMM_WORLD,ierr)

      call MPI_GATHER(rp_proc,ndim*ntraj_proc,mpi_double_precision,rp,
     +                ndim*ntraj_proc,mpi_double_precision,ROOT,
     +                MPI_COMM_WORLD,ierr)

      if (myid == root) then

      open(13,file='temp.dat',action='write')
      
        do i=1,Ntraj
          do j=1,Ndim
            write(13,1000) x(j,i),p(j,i),rp(j,i)
          enddo
        enddo

      close(13)


      do j=1,ndim+1 
        write(114,1000) (cr(j,i), i=1,ndim)        
      enddo 
      do j=1,4 
        write(115,1000) (cr2(j,i),i=1,ndim)
      enddo 


      close(11)
      close(112)
      close(110)
      close(111)
      close(113)
      close(114)
      close(115)
      

! --- deallocate arrays
!      deallocate(fr_proc,du_proc,ap_proc,x_proc,p_proc,rp_proc)

      write(6,8888) 
8888  format('*******************************'/,  
     +       '******** Job finished. ********'/,    
     +       '*******************************')

      endif ! root
      
      call mpi_finalize(ierr)

1000  format(2000(e14.7,1x))
      
      stop
      end program 


! ----------------------------------------------------------------------

!     quit is a subroutine used to terminate execution if there is
!     an error.

!     it is needed here because the subroutine that reads the parameters
!     (subroutine input) may call it.

! ----------------------------------------------------------------------

      subroutine quit

        write (6, *) 'termination via subroutine quit'

        stop
      end subroutine

!--------------------------------------
!      end subroutine
! ----------------------------------------------------------------
! --- propagation of trajs for one time step
! --- each processor have part of the full matrix {x,p,r}, 
!     of size (NDIM,NTRAJ_PROC), run independently
! ----------------------------------------------------------------
      subroutine traj(myid,dt,ndim,ntraj_proc,cf,am,x_proc,p_proc,
     +                rp_proc,ap_proc,wp,du_proc,fr_proc)

      use cdat
      use sys

      implicit real*8(a-h, o-z)

      integer*4,intent(in)    :: myid,ntraj_proc,ndim
      real*8,   intent(in)    :: dt
      real*8, intent(in), dimension(ntraj_proc)  :: wp
      real*8, intent(in), dimension(ndim) :: am,cf
      real*8, intent(inout), dimension(ndim,ntraj_proc) :: x_proc,
     +       p_proc,rp_proc

      real*8, dimension(ndim,ntraj_proc) :: du_proc,fr_proc,ap_proc
      real*8, dimension(ndim)       :: dv

      real*8 :: q(ndim)

      dt2 = dt/2d0

! ----- half-step increments of momenta & full step increment of positions
!      do i = myid*ntraj_proc+1,(myid+1)*ntraj_proc
      do i=1,ntraj_proc
        do m=1,ndim
          q(m) = x_proc(m,i)
        enddo
        
        call local(q,vloc,dv)
!        call long_force(q,vlong,dvl)

        do j=1,Ndim
          p_proc(j,i) = p_proc(j,i)+(-dv(j)-
     +                  du_proc(j,i)-cf(j)*p_proc(j,i)/am(j))*dt2

          rp_proc(j,i)=rp_proc(j,i)+fr_proc(j,i)*dt2

          x_proc(j,i)=x_proc(j,i)+p_proc(j,i)*dt/am(j)
        enddo

      enddo

      return
      end subroutine


! --- propagate momentum by half time step 
      subroutine prop_p(myid,dt,ndim,ntraj_proc,cf,am,x_proc,p_proc,
     +                  rp_proc,wp,du_proc,fr_proc,env_proc,enk_proc)
      
      use cdat
      use sys 

      implicit real*8(a-h, o-z)

!      include 'sizes.h'

!      include 'qsats.h'
      
      integer*4,intent(in)    :: myid,ntraj_proc,ndim
      real*8,   intent(in)    :: dt
      real*8, intent(in), dimension(ntraj_proc)  :: wp
      real*8, intent(in), dimension(ndim) :: am,cf
      real*8, intent(inout), dimension(ndim,ntraj_proc) :: p_proc, 
     +                                        x_proc,rp_proc
      real*8, dimension(ndim,ntraj_proc) :: du_proc,fr_proc
      real*8, dimension(ndim)       :: dv,dvl

      real*8, dimension(ntraj_proc)  :: env,enk 

      real*8 :: q(ndim)

      dt2 = dt/2d0

! ----- half-step increments of momenta & full step increment of positions
!      do i = myid*ntraj_proc+1,(myid+1)*ntraj_proc
      do i=1,ntraj_proc
        do m=1,ndim 
          q(m) = x_proc(m,i)
        enddo
        
        call local(q,vloc,dv)
!        call long_force(q,vlong,dvl)

        env(i) = vloc

        do j=1,Ndim
          p_proc(j,i) = p_proc(j,i)+(-dv(j)-
     +                  du_proc(j,i)-cf(j)*p_proc(j,i)/am(j))*dt2
          
          rp_proc(j,i)=rp_proc(j,i)+fr_proc(j,i)*dt2
        enddo

      enddo

! --- update potential, kinetic, and total energy each proc
      enk = 0d0 
      do i=1,ntraj_proc 
        do j=1,ndim
          enk(i) = enk(i) + p_proc(j,i)**2/(2d0*am(j))
        enddo 
      enddo 
      
      enk_proc = aver_traj(ntraj_proc,wp,enk)
      env_proc = aver_traj(ntraj_proc,wp,env)

      return
      end subroutine

!--------------------------------------

!     distribute each processor with part of the work to get S=f*f

!----------------------------------------------------------- 

! --- construct matrix <p*f|df> to compute correlation 

      subroutine cor_mat(ntraj_proc,ndim,w_proc,x_proc,p_proc,mat)
      
      use cdat, only : nb 

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: ntraj_proc,ndim

      real*8, intent(in), dimension(ndim,ntraj_proc) :: x_proc,p_proc
      real*8, intent(out), dimension(ndim+1,ndim+1) :: mat 

      real*8 :: f(ndim+1),w_proc(ntraj_proc),df(ndim,nb)

      mat = 0d0 

      df = 0d0 
      do k=1,ndim 
        df(k,k) = 1d0 
      enddo 

      do i=1,ntraj_proc

        call basis(ndim,ntraj_proc,i,x_proc,f)
        
        do j=1,ndim 
          do k=1,nb 
            mat(k,j) = mat(k,j)+p_proc(j,i)*f(k)*df(j,j)*w_proc(i)
          enddo 
        enddo 

      enddo

      return
      end subroutine

      subroutine prefit(ntraj_proc,ndim,w_proc,x_proc,p_proc,rp_proc,
     +           s1,cp,cr)

      use cdat, only : nb 

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: ntraj_proc,ndim

      real*8, intent(in), dimension(ndim,ntraj_proc) :: x_proc,
     +                                                p_proc,rp_proc

      real*8, intent(out), dimension(ndim+1,ndim+1) :: s1

      real*8 :: f(ndim+1),w_proc(ntraj_proc),cp(ndim+1,ndim),
     +          cr(ndim+1,ndim),df(ndim,nb)

      s1 = 0d0
      cp = 0d0
      cr = 0d0
      df = 0d0

      do k=1,ndim 
        df(k,k) = 1d0 
      enddo 

      do i=1,ntraj_proc

        call basis(ndim,ntraj_proc,i,x_proc,f)
        
        do k2 = 1,ndim+1
          do k1=1,k2
            s1(k1,k2) = s1(k1,k2)+f(k1)*f(k2)*w_proc(i)
          enddo
        enddo

      enddo

      do i=1,ntraj_proc
        
        call basis(ndim,ntraj_proc,i,x_proc,f)
        
        do k=1,ndim
          do j=1,ndim+1
            cp(j,k) = cp(j,k)+f(j)*p_proc(k,i)*w_proc(i)
            cr(j,k) = cr(j,k)+f(j)*rp_proc(k,i)*w_proc(i)
          enddo
        enddo
      enddo

      do k2=1,ndim+1
        do k1=1,k2
          s1(k2,k1) = s1(k1,k2)
        enddo
      enddo


      return
      end subroutine


! --- intialize global values in module cdat 

      subroutine init_pars(ndim)

      use cdat, only : nb 

      implicit real*8(a-h,o-z)

      nb = ndim+1 

      end subroutine
      
! --- initialize  trajectories 
      
      subroutine init_traj(ndim,ntraj,alpha,x,x0,p,p0,rp) 
      
      use cdat, only : seed

      implicit real*8(a-h,o-z) 

      
      real*8, dimension(ndim,ntraj) :: x,p,r,rp
      
      real*8, dimension(ndim) :: x0,p0,alpha 
      
      integer :: idum(ndim) 
      
      call seed(idum,Ndim)
      
      pow = 6d0

! --- initial Lagrangian grid points (evolve with time)
      do i=1,Ntraj
        do j=1,Ndim
1100      x(j,i) = gasdev(idum(j))
          x(j,i) = x(j,i)/dsqrt(4d0*alpha(j))
          if((x(j,i))**2 .gt. pow/2d0/alpha(j)) goto 1100
        enddo
      enddo

! --- initial momentum for QTs and weights

      do i=1,Ntraj
        do j=1,Ndim
            p(j,i) = p0(j)
            rp(j,i) = -2d0*alpha(j)*(x(j,i))
        enddo
      enddo
      
      return 
      end subroutine 
