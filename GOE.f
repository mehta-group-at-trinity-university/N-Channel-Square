      SUBROUTINE GOE(N,RMAT,ROUTIN)
C
C     ########################################################################
C
C     Valerio Cappellini                                            March 2007
C
C     USES :: rgnf_lux
C     COMPATIBILITY :: Fortran 90
C     VERSION :: 1.0.0
C
C     A random REAL symmetric N x N matrix X is said to belong to the Gaussian
C     Ortogonal Ensemble (GOE) if the diagonal elements x_{ii} and upper diagonal
C     elements x_{ij} are independently chosen with p.d.f.'s
C
C                                    x_{ii}^2
C                      1          - ----------
C     P(x_{ii}) = ------------  e       2
C                  SQRT(2\pi)
C
C      and
C                                    
C                      1         - x_{ij}^2
C     P(x_{ij}) = -----------  e
C                  SQRT(\pi)
C
C     respectively. Equivalently
C
C                   N(N+1)        N
C                - --------    - ---
C     P(X) = \pi      4      2    2   exp[ - Tr X^2 / 2 ]
C
C     In other words, the diagonal entries have distribution N[0,1], where 
C     N[a,b] denotes the "Normal" (or "Gaussian") distribution of mean <a>
C     and variance <b>, while the upper triangular elements have distribution 
C     N[0,1/SQRT(2)]
C         
C
C     X is the output variable <RMAT>.
C
C
C      For N x N matrices taken from GOE, GUE & GSE 
C     (GSE represantable by 2*N x 2*N CPLX matrices)
C      it holds true:
C
C     =========================+===========================+=============+==============
C                              |                           |             |
C     Kind of elements         | # (indep. REAL variables) |  variance   | Prefactor
C                              |                           |             |
C     =========================+===========================+=============+==============
C                              |                           |             |
C     GOE     diagonal entries |             N             |      1      | 1/SQRT(2\pi)
C                              |                           |             |
C     -------------------------+---------------------------+-------------+--------------
C                              |                           |             |
C     GOE off-diagonal entries |          N(N-1)/2         |  1/SQRT(2)  | 1/SQRT(\pi)
C                              |                           |             |
C     -------------------------+---------------------------+-------------+--------------
C                              |                           |             |
C     GUE     diagonal entries |             N             |  1/SQRT(2)  | 1/SQRT(\pi)
C                              |                           |             |
C     -------------------------+---------------------------+-------------+--------------
C                              |                           |             |
C     GUE off-diagonal entries |           N(N-1)          |     1/2     | SQRT(2/\pi)
C                              |                           |             |
C     -------------------------+---------------------------+-------------+--------------
C                              |                           |             |
C     GSE     diagonal entries |             N             |     1/2     | SQRT(2/\pi)
C                              |                           |             |
C     -------------------------+---------------------------+-------------+--------------
C                              |                           |             |
C     GSE off-diagonal entries |         2*N(N-1)          | 1/2*SQRT(2) | 2/SQRT(\pi)
C                              |                           |             |
C     =========================+===========================+=============+==============
C
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     ROUTIN   : Name of the subroutine you are going to use for generating
C                the random numbers (suggested :: rgnf_lux.f)
C
C     IMPORTANT [1] !!! Remember to initialize the Subroutine rgnf_lux.f
C                       by means of "Call Ranset_()"
C     IMPORTANT [2] !!! Remember to put the file "i_seed_rng" into the 
C                       working directory.
C     IMPORTANT [3] !!! Remember to insert main_program: "EXTERNAL rgnf_lux"
C
C     ########################################################################
C
C     $     debug
C     See "A fast normal random number generator" by Joseph Leva for the algorithm used here (NPM-2016)
      
      SAVE  S, T, A, B, R1, R2
      INTEGER N, emme, mu
      DOUBLE PRECISION U(2), RMAT(N,N), S, T, A, B, R1, R2
      DOUBLE PRECISION V, X, Y, Q
      DOUBLE PRECISION twor
      EXTERNAL ROUTIN
      DATA  S, T, A, B / 0.449871d0, -0.386595d0, 0.19600d0, 0.25472d0/
      DATA  R1, R2 / 0.27597d0, 0.27846d0/
C     generate pair of uniform deviates

      
      twor=1.715527769921414d0
      DO emme = 1, N
         DO mu = emme, N
            
 50         CALL ROUTIN(U,2)
C            write(6,*) "pair of uniform deviates returned: ", U(1), U(2)
            V = twor * (U(2) - 0.5d0)
            X = U(1) - S
            Y = ABS(V) - T
            Q = X**2 + Y*(A*Y - B*X)
C     accept P if inside inner ellipse
            IF (Q .LT. R1)  GO TO 100
C     reject P if outside outer ellipse
            IF (Q .GT. R2)  GO TO 50
C     reject P if outside acceptance region
            IF (V**2 .GT. -4.0 *DLOG(U(1)) *U(1)**2)  GO TO 50
C     ratio of P's coordinates is normal deviate
 100        RMAT(emme,mu) = V/U(1)
            IF(emme.ne.mu) THEN 
               RMAT(emme,mu) = RMAT(emme,mu) / SQRT(2.)
               RMAT(mu,emme)=RMAT(emme,mu)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
      
