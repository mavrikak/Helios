c
c     Conjugate gradient subset for program scat3d.
c
c     This is merely a customized version of CCGPAK vers. 1c from 
c     Flatau, Schneider and Evans.  All variables are now in double
c     precision.
c
c
c
      SUBROUTINE cgdoc
C
C                       CCGPAK -- Ver. 1c
C
C
C                       available from:
C                         P. J. Flatau
C              University of California, San Diego
C              Scripps Institution of Oceanography
C                   La Jolla, CA 92093-0221
C                    pflatau@macao.ucsd.edu
C
C                    T. Schneider and F. Evans
C                    Colorado State University
C                 Department of Atmospheric Science
C                   Fort Collins, Colorado 80521
C               schneid@herschel.atmos.colostate.edu
C
C
C  Conjugate gradient package for solving complex matrix equations
C  ----------------------------------------------------------------
C
C     CCGPAK  is designed to solve the system of linear equations
C
C                A*x = y
C
C     where  A  is an  n*n  complex matrix and  y is a given vector.
C     The matrix  A  is not required to be positive definite.
C     The matrix  A  is intended to be large.  It is accessed
C     by means of a subroutine call of the form
C
C                call prod( n,x,y, ... )
C
C     which must return the product  y = A*x  for any given vector  x.
C     It is expected that "prod" routine will be customized to
C     take advantage of structure of A and/or use secondary storage
C     for very large problems. Each routine returns CPU "time" and number
C     "iter" of iterations required to compute an approximate solution to
C     "err" accuracy. "mxiter" iterations are specified on input.
C
C
C Index:
C
C CGSQR ... Complex conjugate gradient-squared algorithm of Sonneveld,
C           for unsymmetric cases.  Encoded from HAV.
C CGSTAB... Complex conjugate gradient squared stabilized
C           encoded from HAV's algorithm Bi-CGSTAB.
C GACG  ... generalized augumented conjugate gradient method for
C           unsymmetric case.
C GBICG ... generalized bi-conjugate gradient.
C GMCG  ... modified conjugate method for unsymmetric complex matrices.
C RCG   ... residual minimized; the residuals are minimized at each
C           iteration, no scaling is introduced.
C SACG  ... augumented conjugate gradient method for symmetric case.
C SBICG ... simplified biconjugate gradient method specialized
C           to symmetric matrices.
C SMCG  ... modified conjugate gradient method for symmetric case.
C SRCG  ... search directions scaled at each iteration and the residuals
C           are minimized.
C XCG   ... the error between the true solution and the approximate
C           solution is minimized at each iteration.
C
C
C On input:
C CA      --- matrix A in Ax=y
C CY      --- known vector y (unchanged on output)
C CX      --- initial guess (can be arbitrary)
C MXN     --- first dimension of matrix A
C N       --- order of problem
C ERR     --- fractional error (e.g. 1.d-6)
C MXITER --- maximum number of iterations
C
C
C On output:
C ITER  ---  number of iterations to get convergence.
C TIME  ---  CPU time (machine dependent, write "timer.f" function).
C CX    ---  solution vector.
C
C
C Externals:
C -- Timer routine.
C -- Random number generator.
C -- "Basic Linear Algebra Package". WARNING this is (PJF)'s
C     implementation of BLAS !
C
C
C References:
C
C HAV= H. A. van der Vorst,
C      "BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR
C      THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS", SIAM J. Sci. Stat.
C      Comput., V.13, #2, pp. 631-644, 1992.
C
C KFE= K. F. Evans, Colorado State University, Dept. Atmos. Sci.
C      For Collins, Colo.
C
C PJF= P. J. Flatau, Scripps Institution of Oceanography,
C      University of California San Diego, La Jolla, California
C
C PPK= M. Petravic, G. Kuo-Petravic, An ILUCG algorithm which
C      minimizes in the Eucleidean norm, J. Comp. Phys., 32, 263--269,
C      1979 (the algorithm there is used in the discrete dipole
C      approximation code by DDSCAT by Draine+Flatau).
C      (it seems to be similiar to SRCG code here).
C
C TKS= T. K. Sarkar, X. Yang, E. Arvas,
C      "A limited survey of various conjugate gradient methods for
C      complex matrix equations arising in electromagnetic wave
C      interactions", Wave motion, 1988, 10, 527---546
C      Dept. Electrical Eng. Syracuse University, Syracuse, NY
C
C TLS= T. L. Schneider, Colorado State University, Dept. Atmos. Sci.
C      For Collins, Colo.
C
      RETURN
      END





      SUBROUTINE cgsqr(ca, cy, cx, mxn, n, err, mxiter, iter, time)
C
C    CGSQR: Complex conjugate gradient-squared algorithm of Sonneveld,
C           for unsymmetric cases.  Encoded from HAV.
C  History: 92/09/03 (TLS)
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      complex*16 alpha, beta, rho, rho0
      real*8 rnorm, ynorm
      INTEGER i
C     .. Local Arrays ..
      complex*16 aw(mxd), cax(mxd), p(mxd), q(mxd), r(mxd), r0(mxd),
     $        u(mxd), v(mxd), w(mxd)
C     .. External Functions ..
      complex*16 cdotpr
      real*8 snorm
      real timer
      EXTERNAL cdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod
C     .. Intrinsic Functions ..
      INTRINSIC dcmplx, dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' cgsqr: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      ynorm = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cax)
      DO 10 i = 1, n
          r(i) = cy(i) - cax(i)
          r0(i) = r(i)
          p(i) = dcmplx(0.d0,0.d0)
          q(i) = dcmplx(0.d0,0.d0)
   10 CONTINUE
      rho0 = dcmplx(1.d0,1.d0)
C
      DO 50 iter = 1, mxiter
          rho = cdotpr(n,r0,'c',r,'u')
          beta = rho/rho0
          DO 20 i = 1, n
              u(i) = r(i) + beta*q(i)
              p(i) = u(i) + beta* (q(i)+beta*p(i))
   20     CONTINUE
          CALL prod(mxn,n,ca,'u',p,'u',v)
          alpha = rho/cdotpr(n,r0,'c',v,'u')
          DO 30 i = 1, n
              q(i) = u(i) - alpha*v(i)
              w(i) = u(i) + q(i)
              cx(i) = cx(i) + alpha*w(i)
   30     CONTINUE
          CALL prod(mxn,n,ca,'u',w,'u',aw)
          DO 40 i = 1, n
              r(i) = r(i) - alpha*aw(i)
   40     CONTINUE
C
          rnorm = snorm(n,r)
c          WRITE (*,FMT=*) 'sqrt(||r||/||y||) =  ', iter,
c     $      dsqrt(rnorm/ynorm)
          IF (dsqrt(rnorm/ynorm) .LT. err) GO TO 60
          rho0 = rho
   50 CONTINUE
C
   60 CONTINUE
      time = timer('get',1)
C
      RETURN
      END





      SUBROUTINE cgstab(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C  CGSTAB:  Complex conjugate gradient squared stabilized
C           encoded from HAV's algorithm Bi-CGSTAB.
C HISTORY:  92/09/08 (TLS)
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      complex*16 alpha, beta, omega, rho, rho0
      real*8 rnorm, ynorm
      INTEGER i
C     .. Local Arrays ..
      complex*16 cax(mxd), p(mxd), r(mxd)
      complex*16 r0(mxd), s(mxd), t(mxd), v(mxd)
C     .. External Functions ..
      complex*16 cdotpr
      real*8 snorm
      EXTERNAL cdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod
C     .. Intrinsic Functions ..
      INTRINSIC dcmplx, dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' cgstab: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      ynorm = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cax)
      DO 10 i = 1, n
          r(i) = cy(i) - cax(i)
          r0(i) = r(i)
          p(i) = dcmplx(0.d0,0.d0)
          v(i) = dcmplx(0.d0,0.d0)
   10 CONTINUE
      rho0 = dcmplx(1.d0,1.d0)
      alpha = dcmplx(1.d0,1.d0)
      omega = dcmplx(1.d0,1.d0)
C
      DO 50 iter = 1, mxiter
          rho = cdotpr(n,r0,'c',r,'u')
          beta = (rho/rho0)* (alpha/omega)
          DO 20 i = 1, n
              p(i) = r(i) + beta* (p(i)-omega*v(i))
   20     CONTINUE
          CALL prod(mxn,n,ca,'u',p,'u',v)
          alpha = rho/cdotpr(n,r0,'c',v,'u')
          DO 30 i = 1, n
              s(i) = r(i) - alpha*v(i)
   30     CONTINUE
          CALL prod(mxn,n,ca,'u',s,'u',t)
          omega = cdotpr(n,t,'c',s,'u')/snorm(n,t)
          DO 40 i = 1, n
              cx(i) = cx(i) + alpha*p(i) + omega*s(i)
              r(i) = s(i) - omega*t(i)
   40     CONTINUE
          rnorm = snorm(n,r)
c          WRITE (*,FMT=*) 'sqrt(||r||/||y||) =  ', iter,
c     $      dsqrt(rnorm/ynorm)
          IF (dsqrt(rnorm/ynorm) .LT. err) GO TO 60
          rho0 = rho
   50 CONTINUE
C
   60 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE gacg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: GACG --- generalized augumented conjugate gradient method
C          for unsymmetric case
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, bk, bw2, q2, sk, sk2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cbq(2*mxd), cbw(2*mxd), cq(2*mxd), cw(2*mxd)
C     .. External Functions ..
      real*8 sdotpr, snorm
      EXTERNAL sdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL czero, prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dconjg, dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' gacg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = 2.*snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cq)
      CALL prod(mxn,n,ca,'c',cx,'u',cq(n+1))
      DO 10 i = 1, n
          cq(i) = cy(i) - cq(i)
          cq(n+i) = dconjg(cy(i)) - cq(n+i)
          cw(i) = cq(i)
          cw(n+i) = cq(n+i)
   10 CONTINUE
      CALL prod(mxn,n,ca,'u',cq(n+1),'u',cbq)
      CALL czero(n,cbq(n+1))
      CALL prod(mxn,n,ca,'c',cq,'u',cbq(n+1))
      sk = 2*sdotpr(n,cbq,'u',cq,'c')
      DO 20 i = 1, 2*n
          cbw(i) = cbq(i)
   20 CONTINUE
      bw2 = snorm(2*n,cbw)
   30 CONTINUE
      ak = sk/bw2
      CALL scincr(n, (ak),cw(n+1),'u',cx)
      CALL scincr(2*n, (-ak),cbw,'u',cq)
      q2 = snorm(2*n,cq)
c      WRITE (*,FMT=*) 'sqrt(q2/ay)= ', iter, dsqrt(q2/ay)
      IF ((dsqrt(q2/ay)) .LT. err) GO TO 50
      CALL prod(mxn,n,ca,'u',cq(n+1),'u',cbq)
      CALL prod(mxn,n,ca,'c',cq,'u',cbq(n+1))
      sk2 = 2.*sdotpr(n,cbq,'u',cq,'c')
      bk = sk2/sk
      DO 40 i = 1, 2*n
          cw(i) = cq(i) + bk*cw(i)
          cbw(i) = cbq(i) + bk*cbw(i)
   40 CONTINUE
      bw2 = snorm(2*n,cbw)
      sk = sk2
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 50
      GO TO 30
   50 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE gbicg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: GBICG --- generalized bi-conjugate gradient
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      complex*16 cak, cbk, csk, csk2
      real*8 ay, ek
      INTEGER i
C     .. Local Arrays ..
      complex*16 cap(mxd), caw(mxd), cp(mxd), cq(mxd), cr(mxd), cw(mxd)
C     .. External Functions ..
      complex*16 cdotpr
      real*8 snorm
      EXTERNAL cdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL cincr, prod
C     .. Intrinsic Functions ..
      INTRINSIC dconjg, dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' gbicg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0.d0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cr)
      csk = (0.d0,0.d0)
      DO 10 i = 1, n
          cr(i) = cy(i) - cr(i)
          cq(i) = dconjg(cr(i))
          cw(i) = cq(i)
          cp(i) = cr(i)
   10 CONTINUE
      csk = cdotpr(n,cr,'u',cr,'u')
   20 CONTINUE
      CALL prod(mxn,n,ca,'c',cw,'u',caw)
      CALL prod(mxn,n,ca,'u',cp,'u',cap)
      cak = cdotpr(n,cap,'u',cw,'c')
      cak = csk/cak
      CALL cincr(n, (cak),cp,'u',cx)
      CALL cincr(n, (-cak),cap,'u',cr)
      CALL cincr(n, (-dconjg(cak)),caw,'u',cq)
      csk2 = cdotpr(n,cr,'u',cq,'c')
      ek = snorm(n,cr)
c      WRITE (*,FMT=*) 'sqrt(ek/ay)= ', iter, dsqrt(ek/ay)
      IF ((dsqrt(ek/ay)) .LT. err) GO TO 40
      cbk = csk2/csk
      DO 30 i = 1, n
          cp(i) = cr(i) + cbk*cp(i)
          cw(i) = cq(i) + dconjg(cbk)*cw(i)
   30 CONTINUE
      csk = csk2
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 40
      GO TO 20
   40 CONTINUE
      time = timer('get',1)
      RETURN
      END




      SUBROUTINE gmcg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: GMCG --- modified conjugate method for unsymmetric complex
C          matrices
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, bk, bk2, q2, w2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cbw(2*mxd), cq(2*mxd), cw(2*mxd), cxh(mxd), cz(mxd)
C     .. External Functions ..
      real*8 sdotpr, snorm
      EXTERNAL sdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dconjg, dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' gmcg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      DO 10 i = 1, n
          cz(i) = cx(i)
          cxh(i) = cx(i)
   10 CONTINUE
      ay = 2.*snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cq)
      CALL prod(mxn,n,ca,'c',cz,'u',cq(n+1))
      DO 20 i = 1, n
          cq(i) = cy(i) - cq(i)
          cq(n+i) = dconjg(cy(i)) - cq(n+i)
          cw(i) = cq(i)
          cw(n+i) = cq(n+i)
   20 CONTINUE
      w2 = snorm(2*n,cw)
   30 CONTINUE
      CALL prod(mxn,n,ca,'u',cw(n+1),'u',cbw)
      CALL prod(mxn,n,ca,'c',cw,'u',cbw(n+1))
      ak = 2.*sdotpr(n,cbw,'u',cw,'c')
c      IF (ak .EQ. 0.d0) THEN
c          PRINT *, ' gmcg.  ak = 0 '
c          STOP
c      END IF
      ak = w2/ak

      CALL scincr(n, (ak),cw(n+1),'u',cx)
      CALL scincr(n, (ak),cw,'u',cz)
      CALL scincr(2*n, (-ak),cbw,'u',cq)
      q2 = snorm(2*n,cq)
      bk = q2/w2
      bk2 = 1. + bk
      DO 40 i = 1, n
          cxh(i) = (cx(i)+bk*cxh(i))/bk2
   40 CONTINUE
      DO 50 i = 1, 2*n
          cw(i) = (cq(i)+bk*cw(i))/bk2
   50 CONTINUE
      w2 = snorm(2*n,cw)
c      WRITE (*,FMT=*) 'sqrt(w2/ay)= ', iter, dsqrt(w2/ay)
      IF ((dsqrt(w2/ay)) .LT. err) GO TO 60
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 60
      GO TO 30
   60 CONTINUE
      DO 70 i = 1, n
          cx(i) = cxh(i)
   70 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE rcg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: RCG --- residual minimized; the residuals are minimized
C          at each iteration, no scaling is introduced.
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, ek, qk, sk, sk2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cp(mxd), cprod(mxd), cr(mxd)
C     .. External Functions ..
      real*8 snorm
      EXTERNAL snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' rcg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cprod)
      DO 10 i = 1, n
          cr(i) = cy(i) - cprod(i)
   10 CONTINUE
      CALL prod(mxn,n,ca,'c',cr,'u',cp)
      sk = snorm(n,cp)
   20 CONTINUE
      CALL prod(mxn,n,ca,'u',cp,'u',cprod)
      ak = snorm(n,cprod)
      ak = sk/ak
      CALL scincr(n, (ak),cp,'u',cx)
      CALL scincr(n, (-ak),cprod,'u',cr)
      ek = snorm(n,cr)
c      WRITE (*,FMT=*) 'sqrt(ek/ay)= ', iter, dsqrt(ek/ay)
      IF ((dsqrt(ek/ay)) .LT. err) GO TO 40
      CALL prod(mxn,n,ca,'c',cr,'u',cprod)
      sk2 = snorm(n,cprod)
      qk = sk2/sk
      DO 30 i = 1, n
          cp(i) = cprod(i) + qk*cp(i)
   30 CONTINUE
      iter = iter + 1
      sk = sk2
      IF (iter .GT. mxiter) GO TO 40
      GO TO 20
   40 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE sacg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: SACG --- augumented conjugate gradient method for symmetric
C          case.
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ap2, ay, bk, r2, sk, sk2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cap(mxd), car(mxd), cp(mxd), cr(mxd)
C     .. External Functions ..
      real*8 sdotpr, snorm
      EXTERNAL sdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL czero, prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' sacg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cr)
      DO 10 i = 1, n
          cr(i) = cy(i) - cr(i)
          cp(i) = cr(i)
   10 CONTINUE
      CALL prod(mxn,n,ca,'u',cr,'c',car)
      sk = sdotpr(n,car,'u',cr,'c')
      DO 20 i = 1, n
          cap(i) = car(i)
   20 CONTINUE
      ap2 = snorm(n,cap)
   30 CONTINUE
      ak = sk/ap2
      CALL scincr(n, (ak),cp,'c',cx)
      CALL scincr(n, (-ak),cap,'u',cr)
      r2 = snorm(n,cr)
      CALL czero(n,car)
c      WRITE (*,FMT=*) 'sqrt(r2/ay)= ', iter, dsqrt(r2/ay)
      IF ((dsqrt(r2/ay)) .LT. err) GO TO 50
      CALL prod(mxn,n,ca,'u',cr,'c',car)
      sk2 = sdotpr(n,car,'u',cr,'c')
      bk = sk2/sk
      DO 40 i = 1, n
          cp(i) = cr(i) + bk*cp(i)
          cap(i) = car(i) + bk*cap(i)
   40 CONTINUE
      ap2 = snorm(n,cap)
      sk = sk2
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 50
      GO TO 30
   50 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE sbicg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: SBICG --- simplified biconjugate gradient method specialized
C          to symmetric matrices
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      complex*16 cak, cbk, csk, csk2
      real*8 ay, ek
      INTEGER i
C     .. Local Arrays ..
      complex*16 cap(mxd), cp(mxd), cr(mxd)
C     .. External Functions ..
      complex*16 cdotpr
      real*8 snorm
      EXTERNAL cdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL cincr, prod
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' sbicg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cr)
      DO 10 i = 1, n
          cr(i) = cy(i) - cr(i)
          cp(i) = cr(i)
   10 CONTINUE
      csk = cdotpr(n,cr,'u',cr,'u')
   20 CONTINUE
      CALL prod(mxn,n,ca,'u',cp,'u',cap)
      cak = cdotpr(n,cap,'u',cp,'u')
      cak = csk/cak
      CALL cincr(n, (cak),cp,'u',cx)
      CALL cincr(n, (-cak),cap,'u',cr)
      csk2 = cdotpr(n,cr,'u',cr,'u')
      ek = snorm(n,cr)
c      WRITE (*,FMT=*) 'sqrt(ek/ay)= ', iter, dsqrt(ek/ay)
      IF ((dsqrt(ek/ay)) .LT. err) GO TO 40
      cbk = csk2/csk
      DO 30 i = 1, n
          cp(i) = cr(i) + cbk*cp(i)
   30 CONTINUE
      csk = csk2
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 40
      GO TO 20
   40 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE smcg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: SMCG --- modified conjugate gradient method for symmetric
C          case.
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, bk, bk2, p2, r2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cap(mxd), cp(mxd), cr(mxd), cxh(mxd)
C     .. External Functions ..
      real*8 sdotpr, snorm
      EXTERNAL sdotpr, snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' smcg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      DO 10 i = 1, n
          cxh(i) = cx(i)
   10 CONTINUE
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cr)
      DO 20 i = 1, n
          cr(i) = cy(i) - cr(i)
          cp(i) = cr(i)
   20 CONTINUE
      p2 = snorm(n,cp)
   30 CONTINUE
      CALL prod(mxn,n,ca,'u',cp,'c',cap)
      ak = sdotpr(n,cap,'u',cp,'c')
      ak = p2/ak
      CALL scincr(n, (ak),cp,'c',cx)
      CALL scincr(n, (-ak),cap,'u',cr)
      r2 = snorm(n,cr)
      bk = r2/p2
      bk2 = 1. + bk
      DO 40 i = 1, n
          cxh(i) = (cx(i)+bk*cxh(i))/bk2
   40 CONTINUE
      DO 50 i = 1, n
          cp(i) = (cr(i)+bk*cp(i))/bk2
   50 CONTINUE
      p2 = snorm(n,cp)
c      WRITE (*,FMT=*) 'sqrt(p2/ay)= ', iter, dsqrt(p2/ay)
      IF ((dsqrt(p2/ay)) .LT. err) GO TO 60
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 60
      GO TO 30
   60 CONTINUE
      DO 70 i = 1, n
          cx(i) = cxh(i)
   70 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE srcg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: SRCG --- search directions scaled at each iteration
C          and the residuals are minimized
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, ek, sk
      INTEGER i
C     .. Local Arrays ..
      complex*16 cp(mxd), cprod(mxd), cr(mxd)
C     .. External Functions ..
      real*8 snorm
      EXTERNAL snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' srcg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cprod)
      DO 10 i = 1, n
          cr(i) = cy(i) - cprod(i)
   10 CONTINUE
      CALL prod(mxn,n,ca,'c',cr,'u',cprod)
      sk = snorm(n,cprod)
      DO 20 i = 1, n
          cp(i) = cprod(i)/sk
   20 CONTINUE
   30 CONTINUE
      CALL prod(mxn,n,ca,'u',cp,'u',cprod)
      ak = snorm(n,cprod)
      ak = 1./ak
      CALL scincr(n, (ak),cp,'u',cx)
      CALL scincr(n, (-ak),cprod,'u',cr)
      ek = snorm(n,cr)
c      WRITE (*,FMT=*) 'sqrt(ek/ay)= ', iter, dsqrt(ek/ay)
      IF ((dsqrt(ek/ay)) .LT. err) GO TO 40
      CALL prod(mxn,n,ca,'c',cr,'u',cprod)
      sk = snorm(n,cprod)
      CALL scincr(n, (1./sk),cprod,'u',cp)
      iter = iter + 1
      IF (iter .GT. mxiter) GO TO 40
      GO TO 30
   40 CONTINUE
      time = timer('get',1)
      RETURN
      END





      SUBROUTINE xcg(ca,cy,cx,mxn,n,err,mxiter,iter,time)
C
C Purpose: XCG --- the error between the true solution and the
C          approximate solution is minimized at each iteration
C
C History: 1988     (TKS); original version
C          92/06/09 (PJF); re-written using "BLAS"-type style
C          92/09/11 (TLS); Removed unnecessary zeroing.
C
C
C     .. Parameters ..
      INTEGER mxd
      PARAMETER (mxd=20000)
C     .. Scalar Arguments ..
      real*8 err
      real time, timer
      INTEGER iter, mxiter, mxn, n
C     .. Array Arguments ..
      complex*16 ca(*), cx(n), cy(n)
C     .. Local Scalars ..
      real*8 ak, ay, qk, sk, sk2
      INTEGER i
C     .. Local Arrays ..
      complex*16 cp(mxd), cprod(mxd), cr(mxd)
C     .. External Functions ..
      real*8 snorm
      EXTERNAL snorm, timer
C     .. External Subroutines ..
      EXTERNAL prod, scincr
C     .. Intrinsic Functions ..
      INTRINSIC dsqrt
c      IF (n .GT. mxd) THEN
c          WRITE (*,FMT='(a,2i4)') ' xcg: correct n or mxd ', n, mxd
c          STOP
c      END IF
      time = timer('set',1)
      iter = 0
      ay = snorm(n,cy)
      CALL prod(mxn,n,ca,'u',cx,'u',cprod)
      DO 10 i = 1, n
          cr(i) = cy(i) - cprod(i)
   10 CONTINUE
      sk = snorm(n,cr)
      CALL prod(mxn,n,ca,'c',cr,'u',cp)
   20 CONTINUE
      CALL prod(mxn,n,ca,'u',cp,'u',cprod)
      ak = snorm(n,cp)
      ak = sk/ak
      CALL scincr(n, (ak),cp,'u',cx)
      CALL scincr(n, (-ak),cprod,'u',cr)
      sk2 = snorm(n,cr)
c      WRITE (*,FMT=*) 'sqrt(sk2/ay)= ', iter, dsqrt(sk2/ay)
      IF ((dsqrt(sk2/ay)) .LT. err) GO TO 40
      CALL prod(mxn,n,ca,'c',cr,'u',cprod)
      qk = sk2/sk
      DO 30 i = 1, n
          cp(i) = cprod(i) + qk*cp(i)
   30 CONTINUE
      iter = iter + 1
      sk = sk2
      IF (iter .GT. mxiter) GO TO 40
      GO TO 20
   40 CONTINUE
      time = timer('get',1)
      RETURN
      END




      SUBROUTINE prod(mxn,n,cmat,chara,cx,charb,cy)
C
C Purpose: calculate y = A x_guess   (cy= cmat cx)
C History: 92/??/?? (PJF)
C          92/09/03 (TLS) Added call to "CZERO".
C
C
C On input:
C
C MXN               --- leading dimension of cmat. To transfer matrix "cmat"
C                       through common block or read it from secondary
C                       storage set "mxn=1" and modify the code.
C N                 --- number of unknowns
C CMAT              --- A(n,n)
C CHARA*1 ('C','U') --- (conjugate/don't conjugate) of matrix elements
C CX(N)             --- current approximation to the solution vector
C CHARB*1 ('C','U') --- (conjugate/don't conjugate) "cx" vector
C
C
C On output:
C
C CY(N)             --- y=Ax
C
C
C     .. Scalar Arguments ..
      INTEGER mxn, n
      CHARACTER chara*1, charb*1
C     .. Array Arguments ..
      complex*16 cmat(mxn,*), cx(n), cy(n)
C     .. Local Scalars ..
      INTEGER i, j
C     .. External Subroutines ..
      EXTERNAL czero
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      CALL czero(n,cy)
      IF (chara .EQ. 'u'.AND.charb .EQ. 'u') THEN
          DO 20 j = 1, n
              DO 10 i = 1, n
                  cy(i) = cy(i) + cmat(i,j)*cx(j)
   10         CONTINUE
   20     CONTINUE
      ELSE IF (chara .EQ. 'u'.AND.charb .EQ. 'c') THEN
          DO 40 j = 1, n
              DO 30 i = 1, n
                  cy(i) = cy(i) + cmat(i,j)*dconjg(cx(j))
   30         CONTINUE
   40     CONTINUE
C CONJUGATE MATRIX "A"
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'u') THEN
          DO 60 i = 1, n
              DO 50 j = 1, n
                  cy(i) = cy(i) + dconjg(cmat(j,i))*cx(j)
   50         CONTINUE
   60     CONTINUE
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'c') THEN
          DO 80 i = 1, n
              DO 70 j = 1, n
                  cy(i) = cy(i) + dconjg(cmat(j,i))*dconjg(cx(j))
   70         CONTINUE
   80     CONTINUE
      END IF
      RETURN
      END





      SUBROUTINE scincr(n,scal,cx,chara,cy)
C
C Purpose: real*8 constant times a vector plus a vector
C          (vectors complex).
C History: 92/06/09 (PJF)
C          92/09/09 (TLS) Added conjugate option.
C
C
C     .. Scalar Arguments ..
      real*8 scal
      INTEGER n
      CHARACTER chara*1
C     .. Array Arguments ..
      complex*16 cx(n), cy(n)
C     .. Local Scalars ..
      INTEGER i
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      IF (chara .EQ. 'c') THEN
          DO 10 i = 1, n
              cy(i) = cy(i) + scal*dconjg(cx(i))
   10     CONTINUE
      ELSE
          DO 20 i = 1, n
              cy(i) = cy(i) + scal*cx(i)
   20     CONTINUE
      END IF
      RETURN
      END





      SUBROUTINE cincr(n,cscal,cx,chara,cy)
C
C Purpose: Complex constant times a vector plus a vector
C          (vectors complex).
C History: 92/06/09 (PJF)
C          92/09/09 (TLS) Added conjugate option.
C
C
C     .. Scalar Arguments ..
      complex*16 cscal
      INTEGER n
      CHARACTER chara*1
C     .. Array Arguments ..
      complex*16 cx(n), cy(n)
C     .. Local Scalars ..
      INTEGER i
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      IF (chara .EQ. 'c') THEN
          DO 10 i = 1, n
              cy(i) = cy(i) + cscal*dconjg(cx(i))
   10     CONTINUE
      ELSE
          DO 20 i = 1, n
              cy(i) = cy(i) + cscal*cx(i)
   20     CONTINUE
      END IF
      RETURN
      END





      SUBROUTINE sincr(n,scal,x,y)
C
C Purpose: real*8 constant times a vector plus a vector
C           (vectors real).
C History: 92/0?/0? (PJF)
C
C
C     .. Scalar Arguments ..
      real*8 scal
      INTEGER n
C     .. Array Arguments ..
      real*8 x(n), y(n)
C     .. Local Scalars ..
      INTEGER i
      DO 10 i = 1, n
          y(i) = y(i) + scal*x(i)
   10 CONTINUE
      RETURN
      END





      SUBROUTINE czero(n,cx)
C
C Purpose: zero vector cx
C History: 92/06/09 (PJF)
C
C
C     .. Scalar Arguments ..
      INTEGER n
C     .. Array Arguments ..
      complex*16 cx(n)
C     .. Local Scalars ..
      INTEGER i
      DO 10 i = 1, n
          cx(i) = (0.d0,0.d0)
   10 CONTINUE
      RETURN
      END





      complex*16 FUNCTION cdotpr(n,cx,chara,cy,charb)
C
C Purpose: a vector times a vector (result is complex)
C History: 92/06/09 (PJF)
C
C
C     .. Scalar Arguments ..
      INTEGER n
      CHARACTER chara*1, charb*1
C     .. Array Arguments ..
      complex*16 cx(n), cy(n)
C     .. Local Scalars ..
      INTEGER i
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      cdotpr = (0.d0,0.d0)
      IF (chara .EQ. 'u'.AND.charb .EQ. 'u') THEN
          DO 10 i = 1, n
              cdotpr = cdotpr + cx(i)*cy(i)
   10     CONTINUE
      ELSE IF (chara .EQ. 'u'.AND.charb .EQ. 'c') THEN
          DO 20 i = 1, n
              cdotpr = cdotpr + cx(i)*dconjg(cy(i))
   20     CONTINUE
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'u') THEN
          DO 30 i = 1, n
              cdotpr = cdotpr + dconjg(cx(i))*cy(i)
   30     CONTINUE
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'c') THEN
          DO 40 i = 1, n
              cdotpr = cdotpr + dconjg(cx(i))*dconjg(cy(i))
   40     CONTINUE
      END IF
      RETURN
      END





      real*8 FUNCTION sdotpr(n,cx,chara,cy,charb)
C
C Purpose: a vector times a vector (result is real)
C History: 92/06/09 (PJF)
C
C
C     .. Scalar Arguments ..
      INTEGER n
      CHARACTER chara*1, charb*1
C     .. Array Arguments ..
      complex*16 cx(n), cy(n)
C     .. Local Scalars ..
      INTEGER i
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      sdotpr = 0.d0
      IF (chara .EQ. 'u'.AND.charb .EQ. 'u') THEN
          DO 10 i = 1, n
              sdotpr = sdotpr + cx(i)*cy(i)
   10     CONTINUE
      ELSE IF (chara .EQ. 'u'.AND.charb .EQ. 'c') THEN
          DO 20 i = 1, n
              sdotpr = sdotpr + cx(i)*dconjg(cy(i))
   20     CONTINUE
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'u') THEN
          DO 30 i = 1, n
              sdotpr = sdotpr + dconjg(cx(i))*cy(i)
   30     CONTINUE
      ELSE IF (chara .EQ. 'c'.AND.charb .EQ. 'c') THEN
          DO 40 i = 1, n
              sdotpr = sdotpr + dconjg(cx(i))*dconjg(cy(i))
   40     CONTINUE
      END IF
      RETURN
      END





      real*8 FUNCTION snorm(n,cx)
C
C Purpose:  norm of a complex vector
C History: 92/06/09 (PJF)
C
C
C     .. Scalar Arguments ..
      INTEGER n
C     .. Array Arguments ..
      complex*16 cx(n)
C     .. Local Scalars ..
      INTEGER i
C     .. Intrinsic Functions ..
      INTRINSIC dconjg
      snorm = 0.d0
      DO 10 i = 1, n
          snorm = snorm + cx(i)*dconjg(cx(i))
   10 CONTINUE
      RETURN
      END




      real FUNCTION timer(func,iwhat)
C
C IBM RS6000 timer. AIX 3.1
C CALL TIMER('SET',IWHAT) --- sets timing for event IWHAT.
C CALL TIMER('GET',IWHAT) --- gets timing for event IWHAT.
C IWHAT --- part of the code. The first call with IWHAT is stored
C           locally in vector IARY with index IWHAT.
C
C
C NOTE:
C Some comments regarding the implementation of "TIMER" function
C on other machines.  If your machine is not listed consult the
C appropriate manuals or return nothing.
C
C    SUN   --- OLDTIM = ETIME(TARRAY) declare real tarray(2)
C    CRAY  --- call second(time1)
C    VMS   --- IAD = 0
C          --- CALL LIB$INIT_TIMER(IAD)   (first call)
C          --- CALL LIB$SHOW_TIMER(IAD,2) (second call)
C    TITAN --- OLDTIM = CPUTIM(0.)
C
C
C
C     .. Scalar Arguments ..
      INTEGER iwhat
      CHARACTER func*3
      real tarray(2)
c      timer=etime(tarray)
      timer=0
      RETURN
      END

