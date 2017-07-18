!!---------------------------------------------------------------------------
!!
!! A simple ocean model
!! Based on one-layer, linearised shallow water equations
!! Written in Fortran 90
!!
!!
!! AUTHORS: Joakim Kjellsson, Laurent Brodeau, Abubakr Babiker Salih
!!
!!---------------------------------------------------------------------------
!! 
!! HISTORY
!!
!! - sept. 2012, Kjellsson & Salih added MPI support
!!
!! - sept. 2010, Brodeau: added netcdf and namelist support
!!
!! - 2010, original code by Joakim Kjellsson
!!
!!---------------------------------------------------------------------------
!!
PROGRAM SHALLOW
  

  USE netcdf
  USE mpi

  
  IMPLICIT none ! Perhaps the most important line in the code
  
  !!
  !! Parameters set in the namelist
  !!
  REAL ::          &
       & D0    = 1000.,       &        !average depth
       & f0    = 1e-4,        &        !Coriolis constant
       & g     = 9.81,        &        !gravity
       & gamma = 0.1,         &        !Asselin coeffjiient
       & Lx    = 5.*1e+7,     &        !width of domain [m]
       & Ly    = 5.*1e+7,     &
       & Ldx   = 1.*1e+6,     &
       & Ldy   = 1.*1e+6,     &
       & tm    = 60*60*24*100.,&       !length of run [s]
       & alfa  = 0.,          &        !slope of bottom
       & beta  = 1e-11,       &        !slope of Coriolis
       & A     = 0.,          &
       & mu    = 0.,          &
       & tau   = 1e-6        
       
  
  INTEGER :: &
       & imt   = 101, &         !number of global x-points
       & jmt   = 101, &         !number of global y-points
       &  nx, ny, nsub,    &    !number of subdomain points
       & nsubcycles = 100, &    !number of cycles between stored time steps
       & nsponge = 10           !number of sponge points
  
  INTEGER*4                              ::  id_nc, ierr, id_x, id_y,     &
      &                                      id_vx, id_vy, id_t, id_u,    &
      &                                      id_v, id_h, id_f, id_d,      &
      &                                      id_su, id_sv, id_sh,         &
      &                                      ierr_mpi, rank, nranks,      &
      &                                      master = 0, column_type, rp, rm 
  
  CHARACTER(LEN=200) :: file,ncdir
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!
  LOGICAL :: lexist
  !!
  REAL :: dx, dy, dt
  REAL :: pi = 3.1415
  !!
  !! Define matrices
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: u, v, h
  REAL*4, DIMENSION(:,:),   ALLOCATABLE :: spongeu, spongev,spongeh, f, D, wind
  REAL*4, DIMENSION(:),     ALLOCATABLE :: vx, vy, vtime
  !!
  !! Indices for loops
  INTEGER :: nt, i, j, k, l, m, n, ji, im, ip, jj, jp, jm, nc, nm, np, jn, tmp, &
  &          ji1, ji2, jj1, jj2
  REAL    :: du,dv, dh
  REAL    :: value
  
  !! Defining namelist section "nsetup":
  NAMELIST /nsetup/ D0, f0, g, gamma, alfa, beta, A, mu, tau, imt, jmt, Lx,&
                    Ly, Ldx, Ldy, tm, nsubcycles
  !!
  !!
  !!
  !! Declarations are done, starting the program
  !!
  !! First testing if the namelist file is present in the current directory:
  INQUIRE(FILE='namelist', EXIST=lexist )
  IF ( .NOT. lexist ) THEN
     PRINT *, 'ERROR: file "namelist" not found!'; STOP
  END IF
  !!
  !! Opening and reading the namelist:
  OPEN( UNIT=11, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )
  READ(11,nsetup)
  CLOSE(11)
  !!
  !!
  file = 'shallowwater.nc'
  ncdir = '/cfs/ekman/scratch/j/jkjel/hpc_project/'
  
  dx = Lx/REAL(imt-1)  !delta x
  dy = Ly/REAL(jmt-1)  !delta y
  !!
  !!
  !! Time step length is determined by grid box size,
  !! a CFL-number, and phase speed (c^2 = g*D)
  dt = 0.2*min(dx,dy)/sqrt(g*D0)
  !! We add +1 to make sure that the model runs to at least tm
  nt = int((int(tm/dt)+1)/nsubcycles)+1  
  !!
  !!
  !! ====================================================
  !! === Start ===
  !! ====================================================
  print*,' === Shallow water model === '
  print*,' dx = ',dx
  print*,' dy = ',dy
  print*,' dt = ',dt
  !!
  !! Test if solution is potentially unstable
  if(A*dt/dx**2 > 1./8) then
       print*,'      ===WARNING!!!==='
       print*,'Potentially unstable solution'
       print*,'A*dt/dx**2 > 1/8'
       print*,'Reduce diffusion!!!'
       stop
  endif
  if(mu*dt > 1./2) then
       print*,'      ===WARNING!!!==='
       print*,'Potentially unstable solution'
       print*,'mu*dt > 1/2'
       print*,'Reduce friction!!!'
       stop
  endif
  

  
  !!----------------------------------------------------------------------
  !!
  !!  Start parallel 
  !!
  !!----------------------------------------------------------------------
  !!
  CALL MPI_INIT(ierr_mpi) !Initialize
  CALL err_mpi(ierr_mpi)
  
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr_mpi) !Get number of processes
  CALL err_mpi(ierr_mpi)
  
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr_mpi)
  CALL err_mpi(ierr_mpi)
  
  !!
  !!----------------------------------------------------------------------
  
  IF ( nranks == 1 ) THEN
     NX = IMT
  ELSE
     NX = INT(FLOAT(IMT)/FLOAT(nranks-1))
  END IF
  NY = JMT
  nsub = INT(FLOAT(IMT)/FLOAT(NX))
  IF ( nsub /= nranks-1 .AND. nranks /= 1) THEN
     PRINT*,'Number of subdomains ',nsub
     PRINT*,'Number of ranks (excl. master) ',nranks-1
     PRINT*,'NX, IMT ',NX,IMT
     STOP
  END IF
  PRINT*,NX,NY,nranks,IMT,JMT
  !!
  !! Define a vector type that we communicate between submatrices
  !!
  CALL MPI_TYPE_VECTOR(1, NX, 0, MPI_FLOAT, column_type, ierr_mpi)
  CALL err_mpi(ierr_mpi)
  
  CALL MPI_TYPE_COMMIT(column_type, ierr)
  CALL err_mpi(ierr_mpi)
    
  !!-----------------------------------------------------------------------
  
  IF (rank == master) THEN
     CALL system('date')
     PRINT*,' Running shallow water model using MPI '
     PRINT*,' Number of processes: ',nranks
  END IF
  PRINT*,' Current rank: ',rank
  
  !!-----------------------------------------------------------------------
  
  !!
  !! Allocating arrays
  !! Add ghost cells to u,v,h to apply boundary conditions 
  !!
  ALLOCATE ( u(0:nx+1,0:ny,3), v(0:nx+1,0:ny+1,3), h(0:nx+1,0:ny,3), &
       &     vx(0:nx+1), vy(0:ny+1), vtime(nt), &
       &     f(0:nx+1,0:ny), D(0:nx+1,0:ny), wind(0:nx+1,0:ny), &
       &     spongeu(0:nx+1,0:ny), spongev(0:nx+1,0:ny+1), spongeh(0:nx+1, 0:ny) )
  
  !ALLOCATE ( xstart(nranks-1), xend(nranks-1), ystart(nranks-1), yend(nranks-1) )
  
  IF (rank /= master ) THEN
     ji1  = 1 + (rank-1) * NX
     ji2  = rank * NX
     jj1  = 1
     jj2  = NY
  END IF
  IF (nranks == 1) THEN
     ji1 = 1
     jj1 = 1
     ji2 = NX
     jj2 = NY
  END IF  

  !! Building coordinates vectors
  !! Loop over extra spatial points so that the 0 and NX+1 points are equal to 
  !! the boundary points of the adjacent cells. 
  DO ji=1,NX
     vx(ji) = (ji+ji1-1)*dx
  END DO
  DO jj=1,NY
     vy(jj) = (jj+jj1-1)*dy
  END DO
  DO n=1,nt
     vtime(n) = (n - 1)*dt*float(nsubcycles)
  END DO
  !!
  !!
  ! ====================================================
  ! === Reset matrices and variables ===
  ! ====================================================
  !!
  du = 0.
  dv = 0.
  dh = 0.
  u(:,:,:) = 0.
  v(:,:,:) = 0.
  h(:,:,:) = 0.
  D(:,:) = 0.
  f(:,:) = 0.
  wind(:,:) = 0.
  spongeu(:,:) = 0.
  spongev(:,:) = 0.
  spongeh(:,:) = 0.
  !!
  !! ====================================================
  !! === Define topography and beta plane ===
  !! ====================================================
  !!
  DO jj=1,NY
       f(:,jj) = f0 + beta*(vy(jj)-Ly/2.)
       D(:,jj) = D0 + alfa*(vy(jj)-Ly/2.)
  ENDDO
  
  !!
  !! ====================================================
  !! === Wind stress / friction ===
  !! ====================================================
  !!
  !!
  DO jj=1,NY
       wind(:,jj) = tau*sin((vy(jj)-Ly/2)/Ly*pi)
  ENDDO
  !!
  !! ====================================================
  !! === Initial condition ===
  !! ====================================================
  
  !! Initial Gaussian disturbance
  DO ji=1,NX
     DO jj=1,NY
        h(ji,jj,1) = 2. * exp(-( ((vx(ji)-Lx/2)/Ldx)**2 + &
        &                        ((vy(jj)-Ly/2)/Ldy)**2 ))
     END DO
  END DO
  
  !!
  !! Send and receive boundary points to fill ghost cells in initial fields.
  !! After sending f(1,:), we must let all ranks receive and store in f(NX+1,:)
  !! before sending something new.
  !! The master does not send nor receive
  !!
  IF (rank /= master) THEN
     
     rm = rank - 1
     IF (rank == 1) THEN
        rm = nranks
     END IF
     
     rp = rank + 1
     IF (rank == nranks) THEN
        rp = 1
     END IF
     
     CALL MPI_SEND(f(1,:), 1, column_type, rm, rank, MPI_COMM_WORLD, ierr_mpi)
     CALL MPI_RECV(f(NX+1,:), 1, column_type, rp, rp, MPI_COMM_WORLD, ierr_mpi)
     
     CALL MPI_SEND(f(NX,:), 1, column_type, rp, rank, MPI_COMM_WORLD, ierr_mpi)
     CALL MPI_RECV(f(0,:), 1, column_type, rm, rm, MPI_COMM_WORLD, ierr_mpi)
     
  END IF
  
  !!
  !! If we are running only one process, it has to communicate with itself
  !!
  IF (nranks == 1) THEN
     
     f(NX+1,:) = f(1,:)
     D(NX+1,:) = D(1,:)
     wind(NX+1,:) = wind(1,:)
     h(NX+1,:,1) = h(1,:,1)
     
     f(0,:) = f(NX,:)
     D(0,:) = D(NX,:)
     wind(0,:) = wind(NX,:)
     h(0,:,1) = h(NX,:,1)
     
  END IF
  
  !!
  !! Set geostrophic balance
  !!
  DO ji=1,NX
     DO jj=2,NY-2
        ip = ji+1 ; im = ji-1
        jp = jj+1 ; jm = jj-1
        IF (ji == 1) THEN
           im = NX
        ELSE IF (ji == NX) THEN
           ip = 1
        END IF
        u(ji,jj,1) = -g / ( 0.5 * (f(ji,jj) + f(im,jj)) ) * 0.25 / dy *  &
                    ( h(ji,jp,1) + h(im,jp,1) - h(ji,jm,1) - h(im,jm,1) )
        v(ji,jj,1) =  g / ( 0.5 * (f(ji,jj) + f(ji,jm)) ) * 0.25 / dx *  &
                    ( h(ip,jj,1) + h(ip,jm,1) - h(im,jj,1) - h(im,jm,1) )
     END DO
  END DO
  
  !!
  !! If we are running only one process, it has to communicate with itself
  !!
  IF (nranks == 1) THEN
     
     u(NX+1,:,1) = u(1,:,1)
     v(NX+1,:,1) = v(1,:,1)
     
     u(0,:,1) = u(NX,:,1)
     v(0,:,1) = v(NX,:,1)
     
  END IF
  
  !! ====================================================
  !! === Sponge zone ===
  !! ====================================================
  !!  
  
  nsponge = 10.
  IF (NY == nsponge) THEN
     PRINT*,' Warning NY = size of sponge zone'
  END IF  

  DO jj=1,nsponge+1
     value = 0.5 + 0.5 * cos(pi*(jj-1)/FLOAT(nsponge))
     spongeu(:,jj) = value
     spongeu(:,NY-jj) = value
     spongev(:,NY-jj+1) = value
     spongev(:,jj) = value
     spongeh(:,jj) = value
     spongeh(:,NY-jj) = value
  END DO

  !! ====================================================
  !! === Initial netCDF ===
  !! ====================================================
  !!
  IF (rank == master) THEN
     CALL output_netcdf('initialise')
  
  END IF
  
  CALL output_netcdf('first')
  !! 
  ! ====================================================
  ! === Main time loop ===
  ! ====================================================
  
  nc = 1
  nm = 3
  np = 2
  
  jn = 1
  CALL output_netcdf('each')
  
  ! Time loops
  DO jn=1,nt
     
     ! Subcycles
     DO jm=1,nsubcycles
                
        ! Loop over y-points
        DO jj=1,NY-1
           
           ! Loop over x-points
           DO ji=1,NX
              
              du = 0.
              dv = 0.
              dh = 0.
              
              CALL calculate_du(ji-1, ji, ji+1, jj-1, jj, jj+1)
              CALL calculate_dv(ji-1, ji, ji+1, jj-1, jj, jj+1)
              CALL calculate_dh(ji-1, ji, ji+1, jj-1, jj, jj+1)
              
              !!
              !! Initialise by Euler forward
              !!
              IF (jn == 1 .AND. jm == 1) THEN
                 u(ji,jj,np) = (u(ji,jj,nc) + du*dt)*(1.-spongeu(ji,jj))
                 v(ji,jj,np) = (v(ji,jj,nc) + dv*dt)*(1.-spongev(ji,jj))
                 h(ji,jj,np) = (h(ji,jj,nc) + dh*dt)*(1.-spongeh(ji,jj)) 
              
              !!
              !! Leap frog
              !!
              ELSE 
                 u(ji,jj,np) = (u(ji,jj,nm) + du*2.*dt)*(1.-spongeu(ji,jj))
                 v(ji,jj,np) = (v(ji,jj,nm) + dv*2.*dt)*(1.-spongev(ji,jj))
                 h(ji,jj,np) = (h(ji,jj,nm) + dh*2.*dt)*(1.-spongeh(ji,jj))
                 
                 ! Asselin filter
                 u(ji,jj,nc) = u(ji,jj,nc) + gamma * &
                              (u(ji,jj,nm) + u(ji,jj,np) - 2. * u(ji,jj,nc))
                 v(ji,jj,nc) = v(ji,jj,nc) + gamma * &
                              (v(ji,jj,nm) + v(ji,jj,np) - 2. * v(ji,jj,nc))
                 h(ji,jj,nc) = h(ji,jj,nc) + gamma * &
                              (h(ji,jj,nm) + h(ji,jj,np) - 2. * h(ji,jj,nc))
              END IF
              
              
           END DO
           
        END DO
        
        v(:,1,:) = 0.
        
        IF (nranks == 1) THEN
           
           u(NX+1,:,nc) = u(1,:,nc)
           v(NX+1,:,nc) = v(1,:,nc)
           h(NX+1,:,nc) = h(1,:,nc)

           u(0,:,nc)    = u(NX,:,nc)
           v(0,:,nc)    = v(NX,:,nc)
           h(0,:,nc)    = h(NX,:,nc)
           
        END IF
     
 ! END IF
        
 !       !!! === Boundaries ===
 !       !!! First: Western boundary
 ! 
 !       ji = 1 
 ! 
 !       DO jj=2,NY-1
 ! 
 ! 
 !         dv = 0.
 !         dh = 0.
 !   
 !         CALL calculate_dv(ji-1, ji, ji+1, jj-1, jj, jj+1)
 !         CALL calculate_dh(ji-1, ji, ji+1, jj-1, jj, jj+1)
 !         
 !         IF (jn == 1 .AND. jm == 1) THEN
 !            v(ji,jj,np) = v(ji,jj,nc) + dv*dt
 !            h(ji,jj,np) = h(ji,jj,nc) + dh*dt
 !         ELSE
 !            ! Time step
 !            v(ji,jj,np) = v(ji,jj,nm) + dv*2.*dt
 !            h(ji,jj,np) = h(ji,jj,nm) + dh*2.*dt
 !            ! Asselin filter
 !            v(ji,jj,nc) = v(ji,jj,nc) + gamma * &
 !                         (v(ji,jj,nm) + v(ji,jj,np) - 2. * v(ji,jj,nc))
 !            h(ji,jj,nc) = h(ji,jj,nc) + gamma * &
 !                         (h(ji,jj,nm) + h(ji,jj,np) - 2. * h(ji,jj,nc))
 !         END IF
 !   
 !       ENDDO

                
 !       !! Second: Southern boundary
 !       
 !       jj = 1 
 !       
 !       DO ji=1,NX
 !          
 !          
 !          du = 0.
 !          dh = 0.
 !          
 !          CALL calculate_du(ji-1, ji, ji+1, jj-1, jj, jj+1)
 !          CALL calculate_dh(ji-1, ji, ji+1, jj-1, jj, jj+1)
 !          
 !          IF (jn == 1 .AND. jm == 1) THEN
 !             u(ji,jj,np) = u(ji,jj,nc) + du*dt
 !             h(ji,jj,np) = h(ji,jj,nc) + dh*dt
 !          ELSE
 !             ! Time step
 !             u(ji,jj,np) = u(ji,jj,nm) + du*2.*dt
 !             h(ji,jj,np) = h(ji,jj,nm) + dh*2.*dt
 !             ! Asselin filter
 !             u(ji,jj,nc) = u(ji,jj,nc) + gamma * &
 !                          (u(ji,jj,nm) + u(ji,jj,np) - 2. * u(ji,jj,nc))    
 !             h(ji,jj,nc) = h(ji,jj,nc) + gamma * &
 !                          (h(ji,jj,nm) + h(ji,jj,np) - 2. * h(ji,jj,nc))
 !          END IF
 !       END DO
 !       
!        !! Third: South-West corner (u(1,1) = v(1,1) = 0 but not h(1,1))
!        jj = 1 
!        ji = 1 
!     
!        dh = 0.
!  
!        CALL calculate_dh(ji-1, ji, ji+1, jj-1, jj, jj+1)
!        
!        IF (jn == 1 .AND. jm == 1) THEN
!           h(ji,jj,np) = h(ji,jj,nc) + dh*dt
!        ELSE
!           ! Time step
!           h(ji,jj,np) = h(ji,jj,nm) + dh*2.*dt
!           ! Asselin filter
!           h(ji,jj,nc) = h(ji,jj,nc) + gamma * &
!                        (h(ji,jj,nm) + h(ji,jj,np) - 2. * h(ji,jj,nc))
!        END IF
        
        !!
        !! Shift indices so that the new becomes the current, the current 
        !! becomes the old, and the old becomes the one to be calcaulated.
        !!
        tmp = nc
        nc = np
        np = nm
        nm = tmp
        
        !!
        !! Sync after each time step
        !!
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
        CALL err_mpi(ierr_mpi)
        
     END DO
     
     IF (rank == master) THEN
        PRINT*,'Done with n = ',jn,' of ',nt,' timesteps'
     END IF
     CALL output_netcdf('each')
     
  END DO
  
  IF (rank == master) THEN
     PRINT*,' All time steps finished'
  END IF
  
  !!
  !! Sync before closing the file
  !!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
  CALL err_mpi(ierr_mpi)  

  IF (rank == master) THEN
     CALL output_netcdf('finalize')
     PRINT*,' End of netCDF output'
  END IF
  
  !! End parallel stuff
  CALL MPI_FINALIZE(ierr_mpi)
  CALL err_mpi(ierr_mpi)
  
  IF (rank == master) THEN
     PRINT*,' End of MPI'
  END IF
  
  ! === THE END ===

  PRINT*,' Done'

CONTAINS

   SUBROUTINE output_netcdf(str_mode)
      
      CHARACTER(LEN=*) :: str_mode
      
      SELECT CASE(str_mode)
      
      CASE ('initialise') 
         
         ierr = NF90_CREATE(TRIM(ncdir)//TRIM(file), NF90_CLOBBER, id_nc) 
         CALL err(ierr)
         
         ierr = NF90_DEF_DIM(id_nc, 'x', IMT, id_x)
         CALL err(ierr)
         ierr = NF90_DEF_DIM(id_nc, 'y', JMT, id_y)
         CALL err(ierr)
         ierr = NF90_DEF_DIM(id_nc, 'time', NF90_UNLIMITED, id_t)
         CALL err(ierr)
         
         ierr = NF90_DEF_VAR(id_nc, 'vx', NF90_FLOAT, id_x, id_vx)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'vy', NF90_FLOAT, id_y, id_vy)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'u', NF90_FLOAT, [id_x, id_y, id_t], id_u)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'v', NF90_FLOAT, [id_x, id_y, id_t], id_v)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'h', NF90_FLOAT, [id_x, id_y, id_t], id_h)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'f', NF90_FLOAT, [id_x, id_y], id_f)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'D', NF90_FLOAT, [id_x, id_y], id_d)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'sponge_u', NF90_FLOAT, [id_x, id_y], id_su)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'sponge_v', NF90_FLOAT, [id_x, id_y], id_sv)
         CALL err(ierr)
         ierr = NF90_DEF_VAR(id_nc, 'sponge_h', NF90_FLOAT, [id_x, id_y], id_sh)
         CALL err(ierr)
         
         ierr = NF90_PUT_ATT(id_nc, id_vx, 'name', 'x position')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_vy, 'name', 'y position')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_u, 'name', 'zonal velocity')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_v, 'name', 'meridional velocity')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_h, 'name', 'depth perturbation')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_f, 'name', 'Coriolis parameter')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_d, 'name', 'mean depth')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_su, 'name', 'relaxation zone for u')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_sv, 'name', 'relaxation zone for v')
         CALL err(ierr)
         ierr = NF90_PUT_ATT(id_nc, id_sh, 'name', 'relaxation zone for h')
         CALL err(ierr)

         ierr = NF90_ENDDEF(id_nc)
         CALL err(ierr)
         
      CASE ('first')
         
         print*,'ji1 = ',ji1, 'jj1 = ',jj1, NX, NY

         ierr = NF90_PUT_VAR(id_nc, id_vx, vx(1:NX), start=[ji1], count=[NX])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_vy, vy(1:NY), start=[jj1], count=[NY])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_f, f(1:NX,1:NY-1), start=[ji1,jj1], &
              &                                           count=[NX,NY-1])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_d, D(1:NX,1:NY-1), start=[ji1,jj1], &
              &                                           count=[NX,NY-1])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_su, spongeu(1:NX,1:NY-1), start=[ji1,jj1], &
              &                                                  count=[NX,NY-1])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_sv, spongev(1:NX,1:NY), start=[ji1,jj1], &
              &                                                count=[NX,NY])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_sh, spongeh(1:NX,1:NY-1), start=[ji1,jj1], &
              &                                                  count=[NX,NY-1])
         CALL err(ierr)
      
      CASE ('each')
         
         ierr = NF90_PUT_VAR(id_nc, id_u, u(1:NX,1:NY-1,nc), start=[ji1,jj1,jn],   &
              &                                      count=[NX,NY-1,1])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_v, v(1:NX,1:NY,nc), start=[ji1,jj1,jn],   &
              &                                      count=[NX,NY,1])
         CALL err(ierr)
         ierr = NF90_PUT_VAR(id_nc, id_h, h(1:NX,1:NY-1,nc), start=[1,1,jn],       &
              &                                      count=[NX,NY-1,1])
         CALL err(ierr)
         
      CASE ('finalize')
         
         ierr = NF90_CLOSE(id_nc)
         CALL err(ierr)
      
      END SELECT
      
   END SUBROUTINE 
      




!!
!! Subroutine that calculates du
!!
SUBROUTINE calculate_du(im,ic,ip,jm,jc,jp)
   
   INTEGER, INTENT(IN) :: im, ic, ip, jm, jc, jp

   ! Height gradient
   du = du - g/dx*( h(ic,jc,nc)-h(im,jc,nc) )
   
   ! Coriolis
   du = du + (f(ic,jc)+f(im,jc))/8.*&
             ( v(ic,jc,nc)+v(im,jc,nc)+v(ic,jp,nc)+v(im,jp,nc) )
   
   ! Wind forcing
   du = du + wind(ic,jc)
   
   ! Linear friction and diffusion
   ! NOTE: Diffusion uses 2nd order derivative, which does not exist near 
   ! the boundaries. Hence the if-statements below.
   du = du - mu * u(ic,jc,nm)
   du = du + A/dx**2 * (u(ip,jc,nm) - 2. * u(ic,jc,nm) + u(im,jc,nm))
   du = du + A/dy**2 * (u(ic,jp,nm) - 2. * u(ic,jc,nm) + u(ic,jm,nm))
   
END SUBROUTINE

!!
!! Subroutine to calculate dv
!!
SUBROUTINE calculate_dv(im,ic,ip,jm,jc,jp)
   
   INTEGER, INTENT(IN) :: im, ic, ip, jm, jc, jp
   
   dv = dv - g/dy * ( h(ic,jc,nc) - h(ic,jm,nc) )
   
   dv = dv - (f(ic,jc) + f(ic,jm)) / 8. * &
   &         ( u(ic,jc,nc) + u(ic,jm,nc) + u(ip,jc,nc) + u(ip,jm,nc) )
   
   !IF (jn > 1) THEN
   dv = dv - mu * v(ic,jc,nm)
   !   IF (ic > 1 .AND. ic < IMT-1) THEN
   dv = dv + A/dx**2 * ( v(ip,jc,nm) - 2 * v(ic,jc,nm) + v(im,jc,nm) )
   !   END IF
   !   IF (jc > 1 .AND. jc < JMT) THEN
   dv = dv + A/dy**2 * (v(ic,jp,nm) - 2 * v(ic,jc,nm) + v(ic,jm,nm))   
   !   END IF
   !END IF
   
END SUBROUTINE

!!
!! Subroutine to calculate dh
!!
SUBROUTINE calculate_dh(im,ic,ip,jm,jc,jp)
   
   INTEGER, INTENT(IN) :: im, ic, ip, jm, jc, jp
   
   ! Convergence/divergence
   dh = dh - D(ic,jc) * ( (u(ip,jc,nc) - u(ic,jc,nc)) / dx + &
   &                      (v(ic,jp,nc) - v(ic,jc,nc)) / dy )
   !IF (ic > 1 .AND. ic < IMT-1) THEN 
   dh = dh - 0.25/dx * (u(ic,jc,nc)+u(ip,jc,nc)) * (D(ip,jc) - D(im,jc))
   !END IF
   !IF (jc > 1 .AND. jc < JMT-1) THEN 
   dh = dh - 0.25/dy * (v(ic,jc,nc)+v(ic,jp,nc)) * (D(ic,jp) - D(ic,jm))
   !END IF
   
END SUBROUTINE


!!
!! Subroutine for error messaging in netCDF
!!
SUBROUTINE err(ierr2)
   
      INTEGER :: ierr2
   
      IF( ierr2 /= 0) THEN
         PRINT*, NF90_STRERROR(ierr2)
         STOP
      END IF
   
   END SUBROUTINE


  SUBROUTINE err_mpi(ierr_mpi2)
   
     INTEGER :: ierr_mpi2, rc
     CHARACTER*300 :: err_string
     
     IF (ierr_mpi2 /= MPI_SUCCESS) THEN
       !CALL MPI_ERROR_STRING(ierr_mpi2, err_string)
       !PRINT*,err_string
       CALL MPI_ABORT(MPI_COMM_WORLD,rc,ierr_mpi2)
     END IF

  END SUBROUTINE

  !!
END PROGRAM SHALLOW
