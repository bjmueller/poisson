c     ================================================================

      MODULE poisson3d

c     ================================================================
c
c     AUTHOR: Bernhard Mueller & Conrad Chan
c     DATE:   1 June MMXIIX
c
c     This module solves the Poisson equation in 3D spherical
c     polar coordinates using a decomposition into eigenfunctions of the
c     angular component of the discrete Laplacian. The transformation
c     between the real and spectral domains is handled by FFTs for the
c     phi-direction and MVMs for the theta-direction.
c
c     Input:
c
c     (for setup routine)
c     m_in:   Number of radial zones
c     n_in:  Number of zones in theta-direction
c     o_in: Number of zones in r-direction
c     r_if_in: radial coordinates at cell interfaces
c     theta_if, theta: theta coordinates at cell interfaces and centres
c     n_s_in, n_e_in: initial and final index of local domain in
c                     theta-direction
c     o_s_in, o_e_in: initial and final index of local domain in
c                     phi-direction
c
c     (For solver)
c     rho:  Source term in the Poisson equation
c     
c     Output:
c     phi:  Solution field
c
c     To compile, a module mo_mpi is needed, which must provide:
c     - the Fortran MPI module or mpif.f,
c     - a Cartesian communicator cart_comm for domain decomposition
c       in the theta and varphi direction.
c     N_theta, N_varphi, and the number of MPI tasks are assumed
c     to be powers of 2.
c
c     Note that the routine solve_poisson_3d does not provide
c     values for phi in ghost zones outside the local MPI domain.
c
c     This module has been adapted from a working version in the
c     CoCoNuT-FMT code and has not been tested oustide the code.
c     
c     ----------------------------------------------------------------

      IMPLICIT NONE

      SAVE

      INTEGER, PARAMETER :: rk  = C_DOUBLE
      INTEGER, PARAMETER :: ik = C_INT32_T

c     ----------------------------------------------------------------
c     Shared module variables
c     ----------------------------------------------------------------
c     For both the FFT and the poisson solver:
      INTEGER (KIND=ik) :: nr, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc

c     ----------------------------------------------------------------
c     For the FFT
      INTEGER (KIND=ik) :: logn
      REAL (KIND=rk) :: pi

      INTEGER (KIND=ik) :: binary (0:20)
      INTEGER (KIND=ik), ALLOCATABLE :: index (:)


c     For the Poisson solver

      INTEGER :: nth_proc

      REAL (KIND=rk) :: dtheta, dphi
      REAL (KIND=rk) :: lambda_phi

      REAL (KIND=rk), ALLOCATABLE :: r (:)
      REAL (KIND=rk), ALLOCATABLE :: r_if (:)
      REAL (KIND=rk), ALLOCATABLE :: theta_if (:)

      REAL (KIND=rk), ALLOCATABLE :: a (:), b(:), c (:)
      REAL (KIND=rk), ALLOCATABLE :: lap_th (:,:)
      REAL (KIND=rk), ALLOCATABLE :: vol_th (:,:)
      
      REAL (KIND=rk), ALLOCATABLE :: matrix (:,:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: lambda (:,:)

      REAL (KIND=rk), ALLOCATABLE :: da_r (:), dv_r (:)
      REAL (KIND=rk), ALLOCATABLE :: diag0 (:), offdiag0 (:)

      CONTAINS


c     ================================================================

      SUBROUTINE setup_poisson_and_fft (m_in, n_in, o_in,
     &     n_s_in, n_e_in, o_e_in, o_s_in, r_if_in, theta_if, theta)

c     ================================================================

      IMPLICIT NONE

c     ----------------------------------------------------------------
c     Subroutine arguments:
      INTEGER (KIND=ik), INTENT (IN) :: m_in, n_in, o_in,
     &     n_s_in, n_e_in, o_e_in, o_s_in
      REAL (KIND=rk), INTENT (IN) :: r_if_in (0:m_in)
      REAL (KIND=rk), INTENT (IN) :: theta_if (0:n_in)
      REAL (KIND=rk), INTENT (IN) :: theta    (1:n_in)

c     ----------------------------------------------------------------
      REAL (KIND=rk) :: w (n_in/2)
      REAL (KIND=rk) :: vec (n_in/2,n_in/2)
      REAL (KIND=rk) :: work (3*n_in) 

      INTEGER (KIND=ik) :: i, j, jj, k, kk, istat, nth, nth2
      INTEGER (KIND=ik) :: info
c     ----------------------------------------------------------------

      pi = ACOS (-1.0d0)
      PRINT*,'Dimensions',m_in,n_in,o_in,n_s_in,n_e_in,o_s_in,o_e_in
      PRINT*,r_if_in(:)
      PRINT*,theta_if(:)
      nr  = m_in
      n   = n_in
      o   = o_in
      n_s = n_s_in
      n_e = n_e_in
      o_s = o_s_in
      o_e = o_e_in

      n_loc = n_e - n_s + 1
      o_loc = o_e - o_s + 1

c     -----------------------------------------------------------------
c     FFT setup:
c     Determine log_2 (qz), integer arithmetic
      j = o
      logn = 0
      DO i = 1, o
         j = j / 2
         IF (j .EQ. 0) EXIT
         logn = logn + 1
         PRINT *,'log_2 (qz) =',logn
      END DO

      ALLOCATE (index(0:o-1), STAT=istat)

      DO i = 0, o-1
         j = i
         DO k = logn-1, 0, -1
            binary (k) = j / 2 ** k
            j = j - binary (k) * 2 ** k
         END DO

         j=0
         DO k = 0, logn-1
            j = j + binary (logn-1-k) * 2 ** k
         END DO
         index (i) = j
      END DO

c     -----------------------------------------------------------------
c     Setup Poisson solver
      nth = n
      nth2 = n / 2
      nth_proc = n / n_loc

      ALLOCATE (lap_th (2,1:nth/2), STAT = istat)
      ALLOCATE (vol_th (1,1:nth), STAT = istat)
      ALLOCATE (a (1:nth), b (1:nth), c (1:nth), STAT = istat)
      ALLOCATE (matrix (1:nth2,1:nth2,o_s:o_e,0:1) ,STAT = istat)
      ALLOCATE (lambda (1:nth,o_s:o_e), STAT = istat)

      ALLOCATE (diag0 (1:nr), offdiag0 (1:nr-1), STAT = istat)
      ALLOCATE (da_r (0:nr), dv_r (1:nr), STAT = istat)
      ALLOCATE (r_if (0:nr), STAT = istat)

      r_if (0:nr) = r_if_in (0:nr)

      DO i = 0, nr
         da_r (i) = r_if (i) ** 2
      END DO
      DO i = 1, nr
         r (i) = 0.5d0 * (r_if (i) + r_if (i-1))
         dv_r (i) = (r_if (i) ** 3 - r_if (i-1) ** 3) / 3.0d0
      END DO      
      DO i = 1, nr-1
         offdiag0 (i) = - da_r (i) / (r (i+1) - r(i))
         diag0 (i) = - offdiag0 (i)
      END DO
      diag0 (nr) = 0.0d0 
      DO i = 2, nr
         diag0 (i) = diag0 (i) - offdiag0 (i-1)
      END DO

c     equidistant grid in theta assumed at the moment
      dtheta = pi / REAL (n, KIND=rk)
      dphi = 2.0d0 * pi / REAL (o, KIND=rk)
      

c     -----------------------------------------------------------------
c     Compute eigenfunctions
c     -----------------------------------------------------------------

      DO kk = o_s, o_e
         k = index (kk-1)         
c     k runs over frequency indices stored by local task
         
         lambda_phi = 
     &        (2.0d0 *  SIN (0.5d0 * DBLE (k) * dphi) / dphi) ** 2 *
     &        dtheta         


         DO j = 1, nth
            a (j) = SIN (theta_if (j-1)) / dtheta ** 2
            b (j) = - (SIN (theta_if (j-1)) + SIN (theta_if (j))) /
     &           dtheta ** 2 - lambda_phi / SIN (theta (j)) / dtheta 
            c (j) = SIN (theta_if (j)) / dtheta ** 2
         END DO

         lap_th (:,:) = 0.0d0
         vol_th (:,:) = 0.0d0

c     Eigenvectors with EVEN parity
         DO j = 1, nth2
            lap_th (1,j) = b (j)
            IF (j .LE. nth2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth2) lap_th (1,j) = b (j) + c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
            PRINT *,'Poisson-3D (operator matrix)',j,lap_th(:,j),
     &           vol_th(1,j)
         END DO
         
         CALL DSBGV ('V','L',nth2,1,0,lap_th(:,1:nth2),2,
     &        vol_th(1,1:nth2),1,w,vec,nth2,work,info)
         
         DO j = 1, nth
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
         END DO
         DO j = 1, nth2 
            matrix (1:nth2,j,kk,0) = vec(1:nth2,j)
         END DO
         lambda (1:nth2,kk) = DBLE (w (1:nth2))
         IF (info .NE. 0) STOP 'DSBGV failed!'
         
c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth/2
            DO i = 1, nth/2
               PRINT *, k, vec (i,j), 
     &              j, i, COS (theta (i)),'eigenvec (even)'
            END DO
         END DO

c     effective l*(l+1)
         DO j = 1, nth/2
            PRINT *,'l*(l+1), even',j, k, w (j), lambda_phi
         END DO

c     Eigenvectors with ODD parity
         lap_th (:,:) = 0.0d0

         DO j = 1, nth2
            lap_th (1,j) = b (j)
            IF (j .LE. nth2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth2) lap_th (1,j) = b (j) - c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
            PRINT *,'Poisson-3D (operator matrix)',j,lap_th(:,j),
     &           vol_th(1,j)
         END DO
         
         CALL DSBGV ('V','L',nth2,1,0,lap_th(:,1:nth2),2,
     &        vol_th(1,1:nth2),1,w,vec,nth2,work,info)
         
         DO j = 1, nth
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
         END DO
         DO j = 1, nth2 
            matrix (1:nth2,j,kk,1) = vec(1:nth2,j)
         END DO
         lambda (nth2+1:nth,kk) = DBLE (w (1:nth2))
         IF (info .NE. 0) STOP 'DSBGV failed!'
         
c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth2
            DO i = 1, nth2
               PRINT *, k, vec (i,j), 
     &              j, i, COS (theta (i)),'eigenvec (odd)'
            END DO
         END DO

c     effective l*(l+1)
         DO j = 1, nth/2
            PRINT *,'l*(l+1), odd',j, k, w (j), lambda_phi
         END DO

      END DO

      DEALLOCATE (lap_th)

      RETURN
      
      END SUBROUTINE setup_poisson_and_fft

c     ================================================================




c     ================================================================

      SUBROUTINE solve_poisson_3d (rho, phi)

c     ================================================================

      USE mo_mpi

      IMPLICIT NONE

      REAL (KIND=rk), INTENT (IN)  :: rho (nr, n_s:n_e, o_s:o_e)
      REAL (KIND=rk), INTENT (OUT) :: phi (nr, n_s:n_e, o_s:o_e)

      REAL (KIND=rk) :: x (n_s:n_e,2*nr,o_s:o_e)
      REAL (KIND=rk) :: y (n_s:n_e,2*nr,o_s:o_e,0:nth_proc/2-1)

      REAL (KIND=rk) :: src (1:nr,2)
      REAL (KIND=rk) :: scr1
      REAL (KIND=rk) :: diag (1:nr), offdiag (1:nr-1)

      DOUBLE COMPLEX :: fft1 (nr,n_s:n_e,o_s:o_e)

      INTEGER (KIND=ik) :: i, j, jj, k, kk, jk, l, l0, k0, inc, 
     &     ntmp, ndat, nth2, info
      
      REAL (KIND=rk) :: sbuf (n_s:n_e,1:2*nr,o_loc,1:nth_proc/2)
      REAL (KIND=rk) :: rbuf (n_s:n_e,1:2*nr,o_loc,1:nth_proc/2)

      DOUBLE COMPLEX :: sbufp (1:nr,n_s:n_e,o_s:o_e)
      DOUBLE COMPLEX :: rbufp (1:nr,n_s:n_e,o_s:o_e)

      INTEGER (KIND=ik) :: srce, dest, coord (2), ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      INTEGER (KIND=ik), parameter :: tag_ysndl=10000_ik
      INTEGER (KIND=ik), parameter :: tag_ysndr=10001_ik
      INTEGER (KIND=ik), parameter :: tag_zsndl=10002_ik
      INTEGER (KIND=ik), parameter :: tag_zsndr=10003_ik

      REAL (KIND=rk) :: tim1, tim2

      nth2 = n/2

      CALL fft_forward (rho (1:nr,n_s:n_e,o_s:o_e),
     &     fft1 (1:nr,n_s:n_e,o_s:o_e), nr)

c     -----------------------------------------------------------------
c     Projection on vectors of even and odd parity

      IF (n_loc .LT. n) THEN

         l0   = (n_s-1) / n_loc
         k0   = (o_s-1) / o_loc

         coord (1) = nth_proc - 1 - l0
         coord (2) = k0
         CALL MPI_Cart_Rank (cart_comm, coord, 
     &        dest, ierr)

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
         END DO

         ndat = n_loc * nr * o_loc

         IF (l0 .LT. nth_proc/2) THEN !even parity

            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik, 
     &           cart_comm, mpistat, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              0.5d0 * (fft1 (1:nr,j,k) + rbufp (1:nr,jj,k))
            END DO

         ELSE !odd parity
            
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik,
     &           cart_comm, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              0.5d0 * (rbufp (1:nr,jj,k) - fft1 (1:nr,j,k))
            END DO
            
         END IF

      ELSE

         DO jk = 1, n_loc/2 * o_loc
            k = int((jk + n_loc/2 - 1) / (n_loc/2) )
            j = (n_s - 1) + (jk - (k - 1) * (n_loc/2))
            jj = (n_e + 1) - (jk - (k - 1) * (n_loc/2))
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
            fft1 (1:nr,j ,k) =
     &           0.5d0 * (fft1  (1:nr,jj,k) + fft1 (1:nr,j ,k))
            fft1 (1:nr,jj,k) =
     &           0.5d0 * (sbufp (1:nr,j ,k) - fft1 (1:nr,jj,k))
         END DO

      END IF

c     -----------------------------------------------------------------
c     Forward Legendre transform, local part:

      DO k = o_s, o_e

         DO i = 1, nr
            x (n_s:n_e,2*i-1,k) =
     &           DBLE  (fft1(i,n_s:n_e,k)) * vol_th(1,n_s:n_e)
            x (n_s:n_e,2*i  ,k) =
     &           DIMAG (fft1(i,n_s:n_e,k)) * vol_th(1,n_s:n_e)
         END DO

c     WARNING: this assumes that nth_proc is a least 2
         DO l = 0, nth_proc/2-1
            IF (n_e .LE. nth2) THEN
c     even parity
               CALL DGEMM ('T', 'N', n_loc, 2*nr , n_loc, 1.0d0,
     &              matrix (n_s:n_e,n_loc*l+1:n_loc*(l+1),k,0),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            ELSE
c     odd parity
               CALL DGEMM ('T', 'N', n_loc, 2*nr , n_loc, 1.0d0,
     &              matrix (n_s-nth2:n_e-nth2,
     &              n_loc*l+1:n_loc*(l+1),k,1),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            END IF
         END DO

      END DO
      
c     -----------------------------------------------------------------
c     MPI part of matrix-matrix multiplication
      
      inc = 1
      DO WHILE (inc .LE. nth_proc / 4)

         ntmp = nth_proc / (2*inc)
         l0   = (n_s-1) / (n_loc * inc)
         k0   = (o_s-1) / o_loc
         ndat = 2*nr*n_loc*o_loc*ntmp/2
            
         IF (MOD (l0, 2) .EQ. 0) THEN ! send right

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2)

            CALL MPI_Cart_shift (cart_comm, 0_ik,  inc,
     &           srce, dest, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik, 
     &           cart_comm, mpistat, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         ELSE !send left

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2)
            CALL MPI_Cart_shift (cart_comm, 0_ik,  -inc,
     &           srce, dest, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik,
     &           cart_comm, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         END IF

c     Reduction step: after the final iteration,
c     y (n_s:n_e,1:2*nr,o_s:o_e,1) contains one row
c     of the matrix-matrix-product
         DO l = 0, ntmp/2 - 1   !DO NOT PARALLELIZE THIS LOOP
            y (n_s:n_e,1:2*nr,o_s:o_e,l) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l) +
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l+1)
         END DO

         inc = 2 * inc

      END DO
      

c     -----------------------------------------------------------------
c     Now solve the tridiagonal system for each (j,k)

      DO kk = o_s, o_e

         k = index (kk-1) !frequency space index

         DO j = n_s, n_e

            DO i = 1, nr
               diag (i) = diag0 (i) - lambda (j,kk) * dv_r (i) / 
     &              r (i) ** 2
               src (i,1) = - y (j,2*i-1,kk,0) * dv_r (i)
               src (i,2) = - y (j,2*i  ,kk,0) * dv_r (i)
            END DO

            DO i = 1, nr-1
               offdiag (i) = offdiag0 (i)
            END DO

            scr1 = SQRT (1.0d0 - 4.0d0 * lambda (j, kk))
            i = nr
            diag (i) = diag (i) + da_r (i) *
     &           (1.0d0 + scr1) / (2.0d0 * r_if (i)) *
     &           (r (i) / r_if (i)) ** scr1

            CALL DPTSV (nr, 2, diag, offdiag, src, nr, info)

            IF (info .EQ. 0) THEN
                  DO i = 1, nr
                     x (j,2*i-1,kk) = src (i,1)
                     x (j,2*i  ,kk) = src (i,2)
                  END DO
            ELSE
               DO i = 1, nr
                  PRINT *,i,diag(i),offdiag(i),src(i,:),lambda(j,kk)
               END DO
               STOP 'DPTSV failed!'
            END IF

         END DO !j-loop

      END DO !kk- and k-loop


c     -----------------------------------------------------------------
c     Backward Legendre transform, local part:

      DO k = o_s, o_e

         DO l = 0, nth_proc/2-1

            IF (n_e .LE. nth2) THEN
c     even parity
               CALL DGEMM ('N', 'N', n_loc, 2*nr, n_loc, 1.0d0,
     &              matrix (n_loc*l+1:n_loc*(l+1),n_s:n_e,k,0),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            ELSE
c     odd parity
               CALL DGEMM ('N', 'N', n_loc, 2*nr, n_loc, 1.0d0,
     &              matrix (n_loc*l+1:n_loc*(l+1),
     &              n_s-nth2:n_e-nth2,k,1),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            END IF
         END DO

      END DO

c     -----------------------------------------------------------------
c     MPI part of matrix-matrix multiplication (works exactly as before)

c     Pack
      inc = 1
      DO WHILE (inc .LE. nth_proc / 4)

         ntmp = nth_proc / (2*inc)
         l0   = (n_s-1) / (n_loc * inc)
         ndat = 2*nr*n_loc*o_loc*ntmp/2

         IF (MOD (l0, 2) .EQ. 0) THEN ! send right
c     
            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2)

            CALL MPI_Cart_shift (cart_comm, 0_ik,  inc,
     &           srce, dest, ierr)
            CALL MPI_Send (sbuf (:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik, 
     &           cart_comm, mpistat, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         ELSE !send left

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2)
            CALL MPI_Cart_shift (cart_comm, 0_ik,  -inc,
     &           srce, dest, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik,
     &           cart_comm, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         END IF

c     Reduction step: after the final iteration,
c     y (n_s:n_e,1:2*nr,o_s:o_e,1) contains one row
c     of the matrix-matrix-product
         DO l = 0, ntmp/2 - 1   !DO NOT PARALLELIZE THIS LOOP
            y (n_s:n_e,1:2*nr,o_s:o_e,l) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l) +
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l+1) 
         END DO

         inc = 2 * inc

      END DO


c     -----------------------------------------------------------------
c     Pack x into complex array for inverse FFT
      DO k = o_s, o_e
         DO j = n_s, n_e
            DO i = 1, nr
               fft1 (i,j,k) = CMPLX (y (j,2*i-1,k,0), y(j,2*i,k,0),
     &              KIND(0.0d0))
            END DO
         END DO
      END DO

c     -----------------------------------------------------------------
c     Add parts with even and odd parity
      IF (n_loc .LT. n) THEN

         l0   = (n_s-1) / n_loc
         k0   = (o_s-1) / o_loc

         coord (1) = nth_proc - 1 - l0
         coord (2) = k0

         CALL MPI_Cart_Rank (cart_comm, coord, 
     &        dest, ierr)

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
         END DO

         ndat = n_loc * nr * o_loc
         
         IF (l0 .LT. nth_proc/2) THEN !even parity

            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik, 
     &           cart_comm, mpistat, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              fft1 (1:nr,j,k) + rbufp (1:nr,jj,k)
            END DO

         ELSE !odd parity
            
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik,
     &           cart_comm, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              rbufp (1:nr,jj,k) - fft1 (1:nr,j,k)
            END DO
            
         END IF

      ELSE

         DO jk = 1, n_loc/2 * o_loc
            k = int((jk + n_loc/2 - 1) / (n_loc/2) )
            j = (n_s - 1) + (jk - (k - 1) * (n_loc/2))
            jj = (n_e + 1) - (jk - (k - 1) * (n_loc/2))
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
            fft1 (1:nr,j ,k) =
     &           fft1  (1:nr,jj,k) + fft1 (1:nr,j ,k)
            fft1 (1:nr,jj,k) =
     &           sbufp (1:nr,j ,k) - fft1 (1:nr,jj,k)
         END DO

      END IF

c     -----------------------------------------------------------------
c     Inverse FFT
      CALL fft_inverse (fft1 (1:nr,n_s:n_e,o_s:o_e),
     &     phi(1:nr,n_s:n_e,o_s:o_e), nr)

      RETURN

      END SUBROUTINE solve_poisson_3d

c     ================================================================



c     =================================================================

      SUBROUTINE fft_forward (xin, yout, m)

c     =================================================================

      USE mo_mpi

      IMPLICIT NONE

c     -----------------------------------------------------------------
c     Subroutine arguments:

      INTEGER (KIND=ik), INTENT (IN) :: m

      REAL (KIND=rk), INTENT (IN ) :: xin  (1:m,n_s:n_e,o_s:o_e)
      DOUBLE COMPLEX, INTENT (OUT) :: yout (1:m,n_s:n_e,o_s:o_e)

c     -----------------------------------------------------------------
      
      REAL (KIND=rk) :: arg

      DOUBLE COMPLEX :: y  (1:m,n_s:n_e,o_s-1:o_e-1)
      DOUBLE COMPLEX :: y0 (1:m,n_s:n_e,o_s-1:o_e-1)

      DOUBLE COMPLEX :: om, omi

      INTEGER (KIND=ik) :: i, j, ji, jf, k, jk, di, di0, di1, inx, nk
      INTEGER (KIND=ik) :: jj, kk

      INTEGER (KIND=ik) :: src, dest, ndat, ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right

c     -----------------------------------------------------------------

      ndat = m * n_loc * o_loc

c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         y0 (1:m,j,k-1) = xin (1:m,j,k)
      END DO



      DO k = 1, logn
         
         di = 2 ** (logn - k)
         nk = 2 ** (k-1)
         arg = - DBLE (2.0d0*pi/DBLE (2*di))
         om = DCMPLX (COS (arg), SIN (arg))

         IF (2*di .LE. o_loc) THEN !no MPI communication needed
            
            ji = (o_s - 1) / (2*di)
            jf = (o_e - 1) / (2*di)
            
            DO j = ji, jf
               di0 = 2*di*j
               di1 = 2*di*j+di
               DO i = 0, di-1
                  omi = DCMPLX (COS (i * arg), SIN (i * arg))
                  y (:,:,di0+i) =
     &                 y0 (:,:,di0+i) + y0 (:,:,di1+i)
                  y (:,:,di1+i) =
     &                 (y0 (:,:,di0+i) - y0 (:,:,di1+i)) * omi
               END DO
            END DO
            
         ELSE                   !MPI communication needed

            IF (MOD ((o_s-1)/di, 2) .EQ. 0) THEN
               
               send_right = .TRUE.
               CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6000_ik, cart_comm, ierr)
               CALL MPI_Recv (y,  ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6001_ik, cart_comm, mpistat, ierr) 
               
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =
     &                 y0 (:,:,i) + y (:,:,i)
               END DO
               
            ELSE
               
               send_right = .FALSE.
               CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6000_ik, cart_comm, mpistat, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6001_ik, cart_comm, ierr) 
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =
     &                 (y (:,:,i) - y0 (:,:,i)) * omi !om ** i
               END DO
            END IF

         END IF

         y0 = y
         
      END DO
      
      DO jk = 1, n_loc * o_loc
         kk = int((jk + n_loc - 1) / n_loc )
         jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
         kk = kk + o_s - 1
         yout (1:m,jj,kk) = y (1:m,jj,kk-1)
      END DO


      RETURN

      END SUBROUTINE fft_forward
      
c     =================================================================



c     =================================================================

      SUBROUTINE fft_inverse (xin, yout, m)

c     =================================================================

      USE mo_mpi

      IMPLICIT NONE

c     -----------------------------------------------------------------
c     Subroutine arguments:

      INTEGER (KIND=ik), INTENT (IN) :: m 

      DOUBLE COMPLEX, INTENT (IN ) :: xin  (1:m,n_s:n_e,o_s:o_e)
      REAL (KIND=rk), INTENT (OUT) :: yout (1:m,n_s:n_e,o_s:o_e)

c     -----------------------------------------------------------------

      REAL (KIND=rk) :: arg

      DOUBLE COMPLEX :: y  (1:m,n_s:n_e,o_s-1:o_e-1)
      DOUBLE COMPLEX :: y0 (1:m,n_s:n_e,o_s-1:o_e-1)

      DOUBLE COMPLEX :: om, omi

      INTEGER (KIND=ik) :: i, j, ji, jf, k, jk, di, di0, di1, inx, nk
      INTEGER (KIND=ik) :: jj, kk

      INTEGER (KIND=ik) :: src, dest, ndat, ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right
c     -----------------------------------------------------------------

      ndat = m * n_loc * o_loc

c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         y0 (1:m,j,k-1) = xin (1:m,j,k)
      END DO

      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2
         y0 (:,j,k) = y0 (:,j,k) / DBLE (o)
      END DO
      
      DO k = logn, 1, -1
         
         di = 2 ** (logn - k)
         nk = 2 ** (k-1)
         arg = DBLE (2.0d0*pi/DBLE (2*di))
         om = DCMPLX (COS (arg), SIN (arg))         

         IF (2*di .LE. o_loc) THEN !no MPI communication needed

            ji = (o_s - 1) / (2*di)
            jf = (o_e - 1) / (2*di)
            
            DO j = ji, jf
               
               di0 = 2*di*j
               di1 = 2*di*j+di
               DO i = 0, di-1
                  omi = DCMPLX (COS (i * arg), SIN (i * arg))
                  y (:,:,di0+i) = !0.5_rk * 
     &                 (y0 (:,:,di0+i) + y0 (:,:,di1+i) * omi)
                  y (:,:,di1+i) = !0.5_rk *
     &                 (y0 (:,:,di0+i) - y0 (:,:,di1+i) * omi)
               END DO
            END DO
            
         ELSE                   !MPI communication needed

            IF (MOD ((o_s-1)/di, 2) .EQ. 0) THEN
               
               send_right = .TRUE.
               CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6000_ik, cart_comm, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6001_ik, cart_comm, mpistat, ierr)
               
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =   !0.5_rk *
     &                 (y0 (:,:,i) + y (:,:,i) * omi)
               END DO

            ELSE

               send_right = .FALSE.
               CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6000_ik, cart_comm, mpistat, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6001_ik, cart_comm, ierr) 
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =   !0.5_rk *
     &                 (y (:,:,i) - y0 (:,:,i) * omi)
               END DO
            END IF

         END IF
         
         y0 = y

      END DO


c     Unpack array:
      DO jk = 1, n_loc * o_loc
         kk = int((jk + n_loc - 1) / n_loc )
         jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
         kk = kk + o_s - 1

         DO i = 1, m
            yout (i,jj,kk) = REAL (y (i,jj,kk-1), KIND=rk)
         END DO

      END DO


      RETURN

      END SUBROUTINE fft_inverse
      
c     =================================================================


      END MODULE poisson3d

c     ================================================================
