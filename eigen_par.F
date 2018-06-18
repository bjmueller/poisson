c     ================================================================

      PROGRAM poisson3d

c     ================================================================
c
c     AUTHOR: Bernhard Mueller & Conrad Chan
c     VERSION: 1.0
c     DATE:   1 June MMXIIX
c
c     This test program solves the Poisson equation in 3D spherical
c     polar coordinates using a decomposition into eigenfunctions of the
c     angular component of the discrete Laplacian. The transformation
c     between the real and spectral domains is handled by FFTs for the
c     phi-direction and MVMs for the theta-direction.
c
c     Input:
c     nr:   Number of radial zones
c     nth:  Number of zones in theta-direction
c     nphi: Number of zones in r-direction
c     rho:  Source term in the Poisson equation
c     
c     Output:
c     phi:  Solution field
c
c     ----------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'fftw3.f'

      INTEGER, PARAMETER :: nr = 550, nth = 128, nphi = 256
      
      DOUBLE PRECISION :: phi (nr, nth, nphi)
      DOUBLE PRECISION :: rho (nr, nth, nphi)
      DOUBLE PRECISION :: a (nth), b(nth), c (nth)
      DOUBLE PRECISION :: dr, dtheta, dphi
      DOUBLE PRECISION :: r_if (0:nr), r (1:nr)
      DOUBLE PRECISION :: theta_if (0:nth), theta (nth)
      DOUBLE PRECISION :: varphi (1:nphi)
      DOUBLE PRECISION :: pc_pi
      DOUBLE PRECISION :: lambda_phi

      DOUBLE PRECISION :: lap_th (2,nth/2)
      DOUBLE PRECISION :: vol_th (1,nth)
      DOUBLE PRECISION :: w (1:nth)
      DOUBLE PRECISION :: vec (nth/2,nth/2)
      DOUBLE PRECISION :: lap_th2 (2,nth)
      DOUBLE PRECISION :: vec2 (nth,nth)
      DOUBLE PRECISION :: work (3*nth) 
      
      DOUBLE PRECISION :: matrix (nth/2,nth,0:nphi-1)
      DOUBLE PRECISION :: lambda (nth,0:nphi-1)
      DOUBLE PRECISION :: x (nth,2*nr), y (nth,2*nr)

      DOUBLE PRECISION :: da_r (0:nr), dv_r (1:nr)
      DOUBLE PRECISION :: diag (1:nr), offdiag (1:nr-1),
     &     diag0 (1:nr), offdiag0 (1:nr-1), src (1:nr,2)
      DOUBLE PRECISION :: scr1

      DOUBLE COMPLEX :: fft1 (nr,nth,0:nphi/2)
      DOUBLE PRECISION :: yeven, yodd

      INTEGER :: plan_fft, plan_inv, info

      INTEGER :: i, j, k, m, ip1, im1

      DOUBLE PRECISION :: tim1, tim2, tim3

      CALL CPU_TIME (tim1)

      pc_pi = ACOS (-1.0d0)
      dr = 1.0d0
      dphi = 2.0 * pc_pi / DBLE (nphi)

      r_if (0) = 0.0d0
      r_if (1) = 4.0d4
      DO i = 2, nr
         r_if (i) = 1.02d0 * r _if (i-1)
      END DO
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
      DO i = 2, nr
         diag0 (i) = diag0 (i) - offdiag0 (i-1)
      END DO


      DO j = 0, nth
         theta_if (j) = DBLE (j) / DBLE (nth) * pc_pi
      END DO
      DO j = 1, nth
         theta (j) = 0.5d0 * (theta_if (j-1) + theta_if (j))
      END DO
      dtheta = pc_pi / DBLE (nth)

      DO m = 1, nphi
         varphi (m) = (DBLE (m) - 0.5d0) / DBLE (nphi) * 2.0d0 * pc_pi
      END DO
      dphi = 2.0d0 * pc_pi / DBLE (nphi)
      

c     -----------------------------------------------------------------
c     Compute eigenfunctions
c     -----------------------------------------------------------------

      matrix (:,:,:) = 0.0d0

      DO m = 0, nphi/2

         lambda_phi = 
     &        (2.0d0 *  SIN (0.5d0 * DBLE (m) * dphi) / dphi) ** 2 *
     &        dtheta         

         DO j = 1, nth
            a (j) = SIN (theta_if (j-1)) / dtheta ** 2
            b (j) = - (SIN (theta_if (j-1)) + SIN (theta_if (j))) /
     &           dtheta ** 2 - lambda_phi / SIN (theta (j)) /
     &           dtheta ** 1
            c (j) = SIN (theta_if (j)) / dtheta ** 2
         END DO

         lap_th (:,:) = 0.0d0
         vol_th (:,:) = 0.0d0

c     Eigenfunctions with EVEN parity
         DO j = 1, nth/2
            lap_th (1,j) = b (j)
            IF (j .LE. nth/2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth/2) lap_th (1,j) = b (j) + c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta ** 1
         END DO
         
         CALL DSBGV ('V','L',nth/2,1,0,lap_th,2,vol_th(1,1:nth/2),
     &        1,w(1:nth/2),vec,nth/2,work,info)
         
         DO j = 1, nth/2
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta ** 1
            matrix (1:nth/2,j,m) = vec(1:nth/2,j)
         END DO
         lambda (1:nth/2,m) = DBLE (w (1:nth/2))
         IF (info .NE. 0) STOP 'DSBGV failed!'

c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth/2
            DO k = 1, nth/2
               WRITE (20,'(i4,x,e16.7,x,2i4,x,e16.7," even")') 
     &              k, vec (k,j)/SQRT(2.0d0), j, m, COS (theta (k))
            END DO
            WRITE (20,*)
         END DO

c     Eigenfunctions with ODD parity
         DO j = 1, nth/2
            lap_th (1,j) = b (j)
            IF (j .LE. nth/2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth/2) lap_th (1,j) = b (j) - c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta ** 1
         END DO
         
         CALL DSBGV ('V','L',nth/2,1,0,lap_th,2,vol_th(1,1:nth/2),
     &        1,w(nth/2+1:nth),vec,nth/2,work,info)
         
         DO j = 1, nth
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta ** 1
         END DO
         DO j = 1, nth/2
            matrix (1:nth/2,j+nth/2,m) = vec(1:nth/2,j)
         END DO
         lambda (nth/2+1:nth,m) = DBLE (w (nth/2+1:nth))
         IF (info .NE. 0) STOP 'DSBGV failed!'
         
c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth/2
            DO k = 1, nth/2
               WRITE (20,'(i4,x,e16.7,x,2i4,x,e16.7," odd")')
     &              k, vec (k,j)/SQRT(2.0d0), j, m, COS (theta (k))
            END DO
            WRITE (20,*)
         END DO

c     effective l*(l+1)
         DO j = 1, nth
            WRITE (21,*) j, m, w (j), lambda_phi
         END DO

      END DO



c     -----------------------------------------------------------------
c     Solve 3D poisson equation
c     -----------------------------------------------------------------

c     Density field for test
      rho (:,:,:) = 0.0d0
      rho (2*(nr/3),nth/4,nphi/2) = 1.0d0

      CALL CPU_TIME (tim2)


c     Forward FFT
      CALL dfftw_plan_many_dft_r2c (plan_fft, 1, nphi, nr * nth, 
     &     rho (:,:,:), nphi, nr*nth, 1,
     &     fft1(:,:,:), nphi/2+1, nr*nth, 1, FFTW_ESTIMATE)

      CALL dfftw_execute_dft_r2c (plan_fft, rho, fft1)

      DO m = 0, nphi/2

         DO i = 1, nr
            y (:,2*i-1) = DBLE  (fft1(i,:,m)) * vol_th(1,:)
            y (:,2*i  ) = DIMAG (fft1(i,:,m)) * vol_th(1,:)
         END DO

c     Split into components of even and odd parity
         DO i = 1, 2*nr
            DO j = 1, nth/2
               x (j      ,i) = 0.5d0 * (y (j,i) + y (nth+1-j,i))
               x (j+nth/2,i) = 0.5d0 * (y (j,i) - y (nth+1-j,i))
            END DO
         END DO

         CALL DGEMM ('T', 'N', nth/2, 2*nr , nth/2, 1.0d0,
     &        matrix (1:nth/2,1:nth/2,m),
     &        nth/2, x (1:nth/2,:), nth/2, 0.0d0,
     &        y (1:nth/2,:), nth/2)
         CALL DGEMM ('T', 'N', nth/2, 2*nr , nth/2, 1.0d0,
     &        matrix (1:nth/2,nth/2+1:nth,m),
     &        nth/2, x (nth/2+1:nth,:), nth/2, 0.0d0,
     &        y (nth/2+1:nth,:), nth/2)


         DO j = 1, nth
c     Solve tridiagonal system
            DO i = 1, nr
               diag (i) = diag0 (i) - lambda (j,m) * dv_r (i) / 
     &              r (i) ** 2
               src (i,1) = - y (j,2*i-1) * dv_r (i)
               src (i,2) = - y (j,2*i  ) * dv_r (i)
            END DO
            DO i = 1, nr-1
               offdiag (i) = offdiag0 (i)
            END DO
            scr1 = SQRT (1.0d0 - 4.0d0 * lambda (j, m))
            i = nr
            diag (i) = diag (i) + da_r (i) *
     &           (1.0d0 + scr1) / (2.0d0 * r_if (i)) 
     &           *(r(i)/r_if(i))**scr1
            CALL DPTSV (nr, 2, diag, offdiag, src, nr, info)
            IF (info .EQ. 0) THEN
               DO i = 1, nr
                  y (j,2*i-1) = src (i,1)
                  y (j,2*i  ) = src (i,2)
               END DO
            ELSE
               STOP 'DPTSV failed!'
            END IF
         END DO

         CALL DGEMM ('N', 'N', nth/2, 2*nr, nth/2, 1.0d0,
     &        matrix (1:nth/2,1:nth/2,m),
     &        nth/2, y (1:nth/2,:), nth/2, 0.0d0,
     &        x (1:nth/2,:), nth/2)
         CALL DGEMM ('N', 'N', nth/2, 2*nr , nth/2, 1.0d0,
     &        matrix (1:nth/2,nth/2+1:nth,m),
     &        nth/2, y (nth/2+1:nth,:), nth/2, 0.0d0, 
     &        x (nth/2+1:nth,:), nth/2)


c     Combine components of even and odd parity
         DO i = 1, 2*nr
            DO j = 1, nth/2
               y (j      ,i) = x (j,i) + x (j+nth/2,i)
               y (nth+1-j,i) = x (j,i) - x (j+nth/2,i)
            END DO
         END DO

         DO j = 1, nth
            DO i = 1, nr
               fft1 (i,j,m) = CMPLX (y (j,2*i-1), y(j,2*i), KIND(0.0d0))
            END DO
         END DO
      END DO

c     Inverse FFT
      CALL dfftw_plan_many_dft_c2r (plan_inv, 1, nphi, nr * nth, 
     &     fft1(:,:,:), nphi/2+1, nr*nth, 1,
     &     phi (:,:,:), nphi, nr*nth, 1, FFTW_ESTIMATE)

      CALL dfftw_execute_dft_c2r (plan_inv, fft1, phi)

      CALL CPU_TIME (tim3)

      phi = phi / DBLE (nphi) * 4.0d0 * pc_pi

      PRINT '("Setup:   ",f6.2,x,"s")',tim2-tim1
      PRINT '("Poisson: ",f6.2,x,"s")',tim3-tim2

      OPEN (30, FILE = 'phi.dat', FORM = 'UNFORMATTED')
      WRITE (30) phi, r (1:nr), theta (1:nth), varphi (1:nphi)
      CLOSE (30)

      scr1 = 0.0d0
      DO k = 1, nphi
         DO j = 1, nth
            DO i = 1, nr
               scr1 = scr1 +
     &              rho (i,j,k)* dv_r (i) * vol_th (1,j) * dtheta * dphi
            END DO
         END DO
      END DO

      PRINT *
      PRINT '("Total mass: ",e26.17)', scr1

      END PROGRAM poisson3d

c     ================================================================
