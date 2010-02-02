*DECK DWNNLS
      SUBROUTINE DWNNLS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IWORK, WORK)
C***BEGIN PROLOGUE  DWNNLS
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality constraints and nonnegativity constraints on
C            selected variables.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A
C***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D)
C***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
C             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
C             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem.  Suppose there are given matrices E and A of
C     respective dimensions ME by N and MA by N, and vectors F
C     and B of respective lengths ME and MA.  This subroutine
C     solves the problem
C
C               EX = F, (equations to be exactly satisfied)
C
C               AX = B, (equations to be approximately satisfied,
C                        in the least squares sense)
C
C               subject to components L+1,...,N nonnegative
C
C     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
C
C     The problem is reposed as problem DWNNLS
C
C               (WT*E)X = (WT*F)
C               (   A)    (   B), (least squares)
C               subject to components L+1,...,N nonnegative.
C
C     The subprogram chooses the heavy weight (or penalty parameter) WT.
C
C     The parameters for DWNNLS are
C
C     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     W(*,*),MDW,  The array W(*,*) is double subscripted with first
C     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
C                  discussion let us call M = ME + MA.  Then MDW
C                  must satisfy MDW.GE.M.  The condition MDW.LT.M
C                  is an error.
C
C                  The array W(*,*) contains the matrices and vectors
C
C                       (E  F)
C                       (A  B)
C
C                  in rows and columns 1,...,M and 1,...,N+1
C                  respectively.  Columns 1,...,L correspond to
C                  unconstrained variables X(1),...,X(L).  The
C                  remaining variables are constrained to be
C                  nonnegative. The condition L.LT.0 or L.GT.N is
C                  an error.
C
C     PRGOPT(*)    This double precision array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case LINK=1 and the values KEY and DATA SET
C                  are not referenced. The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=LINK1 (link to first entry of next group)
C               .  PRGOPT(2)=KEY1 (key to the option change)
C               .  PRGOPT(3)=DATA VALUE (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2)=DATA VALUE
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK.GT.NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000 an error
C                  message is printed and the subprogram returns.
C
C                  OPTIONS..
C
C                  KEY=6
C                         Scale the nonzero columns of the
C                  entire data matrix
C                  (E)
C                  (A)
C                  to have length one. The DATA SET for
C                  this option is a single value.  It must
C                  be nonzero if unit length column scaling is
C                  desired.
C
C                  KEY=7
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  with a user-provided diagonal matrix.
C                  The DATA SET for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=8
C                         Change the rank determination tolerance from
C                  the nominal value of SQRT(SRELPR).  This quantity
C                  can be no smaller than SRELPR, The arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The DATA SET for this option
C                  is the new tolerance.
C
C                  KEY=9
C                         Change the blow-up parameter from the
C                  nominal value of SQRT(SRELPR).  The reciprocal of
C                  this parameter is used in rejecting solution
C                  components as too large when a variable is
C                  first brought into the active set.  Too large
C                  means that the proposed component times the
C                  reciprocal of the parameter is not less than
C                  the ratio of the norms of the right-side
C                  vector and the data matrix.
C                  This parameter can be no smaller than SRELPR,
C                  the arithmetic-storage precision.
C
C                  For example, suppose we want to provide
C                  a diagonal matrix to scale the problem
C                  matrix and change the tolerance used for
C                  determining linear dependence of dropped col
C                  vectors.  For these options the dimensions of
C                  PRGOPT(*) must be at least N+6.  The FORTRAN
C                  statements defining these options would
C                  be as follows.
C
C                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
C                  PRGOPT(2)=7 (user-provided scaling key)
C
C                  CALL DCOPY(N,D,1,PRGOPT(3),1) (copy the N
C                  scaling factors from a user array called D(*)
C                  into PRGOPT(3)-PRGOPT(N+2))
C
C                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
C                  PRGOPT(N+4)=8 (linear dependence tolerance key)
C                  PRGOPT(N+5)=... (new value of the tolerance)
C
C                  PRGOPT(N+6)=1 (no more options to change)
C
C
C     IWORK(1),    The amounts of working storage actually allocated
C     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
C                  respectively.  These quantities are compared with
C                  the actual amounts of storage needed for DWNNLS( ).
C                  Insufficient storage allocated for either WORK(*)
C                  or IWORK(*) is considered an error.  This feature
C                  was included in DWNNLS( ) because miscalculating
C                  the storage formulas for WORK(*) and IWORK(*)
C                  might very well lead to subtle and hard-to-find
C                  execution errors.
C
C                  The length of WORK(*) must be at least
C
C                  LW = ME+MA+5*N
C                  This test will not be made if IWORK(1).LE.0.
C
C                  The length of IWORK(*) must be at least
C
C                  LIW = ME+MA+N
C                  This test will not be made if IWORK(2).LE.0.
C
C     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     X(*)         An array dimensioned at least N, which will
C                  contain the N components of the solution vector
C                  on output.
C
C     RNORM        The residual norm of the solution.  The value of
C                  RNORM contains the residual vector length of the
C                  equality constraints and least squares equations.
C
C     MODE         The value of MODE indicates the success or failure
C                  of the subprogram.
C
C                  MODE = 0  Subprogram completed successfully.
C
C                       = 1  Max. number of iterations (equal to
C                            3*(N-L)) exceeded. Nearly all problems
C                            should complete in fewer than this
C                            number of iterations. An approximate
C                            solution and its corresponding residual
C                            vector length are in X(*) and RNORM.
C
C                       = 2  Usage error occurred.  The offending
C                            condition is noted with the error
C                            processing subprogram, XERMSG( ).
C
C     User-designated
C     Working arrays..
C
C     WORK(*)      A double precision working array of length at least
C                  M + 5*N.
C
C     IWORK(*)     An integer-valued working array of length at least
C                  M+N.
C
C***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Report SAND77-0552, Sandia
C                 Laboratories, June 1978.
C               K. H. Haskell and R. J. Hanson, Selected algorithms for
C                 the linearly constrained least squares problem - a
C                 users guide, Report SAND78-1290, Sandia Laboratories,
C                 August 1979.
C               K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Mathematical Programming
C                 21 (1981), pp. 98-118.
C               R. J. Hanson and K. H. Haskell, Two algorithms for the
C                 linearly constrained least squares problem, ACM
C                 Transactions on Mathematical Software, September 1982.
C               C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974.
C***ROUTINES CALLED  DWNLSM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls, change Prologue
C           comments to agree with WNNLS.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DWNNLS
      INTEGER IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW, ME,
     *     MODE, N
      DOUBLE PRECISION  PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
      CHARACTER*8 XERN1
C***FIRST EXECUTABLE STATEMENT  DWNNLS
      MODE = 0
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
C
      IF (IWORK(1).GT.0) THEN
         LW = ME + MA + 5*N
         IF (IWORK(1).LT.LW) THEN
            WRITE (XERN1, '(I8)') LW
            CALL XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (IWORK(2).GT.0) THEN
         LIW = ME + MA + N
         IF (IWORK(2).LT.LIW) THEN
            WRITE (XERN1, '(I8)') LIW
            CALL XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (MDW.LT.ME+MA) THEN
         CALL XERMSG ('SLATEC', 'DWNNLS',
     *      'THE VALUE MDW.LT.ME+MA IS AN ERROR', 1, 1)
         MODE = 2
         RETURN
      ENDIF
C
      IF (L.LT.0 .OR. L.GT.N) THEN
         CALL XERMSG ('SLATEC', 'DWNNLS',
     *      'L.GE.0 .AND. L.LE.N IS REQUIRED', 2, 1)
         MODE = 2
         RETURN
      ENDIF
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
C     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
C     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
C
      L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
C
      CALL DWNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,
     *            IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),
     *            WORK(L4), WORK(L5))
      RETURN
      END
*DECK DWNLSM
      SUBROUTINE DWNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
C***BEGIN PROLOGUE  DWNLSM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DWNNLS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to DWNNLS.
C     The documentation for DWNNLS has complete usage instructions.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine DWNNLS, the following work arrays are used in
C     subroutine DWNLSM  (they are passed through the calling
C     sequence from DWNNLS for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C     Variables of type REAL are DOUBLE PRECISION.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C
C***SEE ALSO  DWNNLS
C***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM,
C                    DROTMG, DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Fixed an error message.  (RWC)
C   900604  DP version created from SP version.  (RWC)
C   900911  Restriction on value of ALAMDA included.  (WRB)
C***END PROLOGUE  DWNLSM
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      DOUBLE PRECISION D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),
     *   W(MDW,*), WD(*), X(*), Z(*)
C
      EXTERNAL D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM, DROTMG,
     *   DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG
      DOUBLE PRECISION D1MACH, DASUM, DNRM2
      INTEGER IDAMAX
C
      DOUBLE PRECISION ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,
     *   DOPE(3), DRELPR, EANORM, FAC, SM, SPARAM(5), T, TAU, WMAX, Z2,
     *   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,
     *   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,
     *   NOPT, NSOLN, NTIMES
      LOGICAL DONE, FEASBL, FIRST, HITCON, POS
C
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DWNLSM
C
C     Initialize variables.
C     DRELPR is the precision for the particular machine
C     being used.  This logic avoids resetting it every entry.
C
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
C
C     Set the nominal tolerance used in the code.
C
      TAU = SQRT(DRELPR)
C
      M = MA + MME
      ME = MME
      MODE = 2
C
C     To process option vector
C
      FAC = 1.D-4
C
C     Set the nominal blow up factor used in the code.
C
      BLOWUP = TAU
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
      CALL DCOPY (N, 1.D0, 0, D, 1)
C
C     Define bound for number of options to change.
C
      NOPT = 1000
C
C     Define bound for positive value of LINK.
C
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMSG ('SLATEC', 'DWNLSM',
     +      'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL XERMSG ('SLATEC', 'DWNLSM',
     +      'IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.',
     +      3, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = DNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               D(J) = T
  110       CONTINUE
         ENDIF
C
         IF (KEY.EQ.7) CALL DCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(DRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(DRELPR,PRGOPT(LAST+2))
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL XERMSG ('SLATEC', 'DWNLSM',
     +         'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL DSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
C
C     Process option vector
C
      DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
C
C     Compute scale factor to apply to equality constraint equations.
C
      DO 130 J = 1,N
         WD(J) = DASUM(M,W(1,J),1)
  130 CONTINUE
C
      IMAX = IDAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = DASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(DRELPR*FAC)
C
C     On machines, such as the VAXes using D floating, with a very
C     limited exponent range for double precision values, the previously
C     computed value of ALAMDA may cause an overflow condition.
C     Therefore, this code further limits the value of ALAMDA.
C
      ALAMDA = MIN(ALAMDA,SQRT(D1MACH(2)))
C
C     Define scaling diagonal matrix for modified Givens usage and
C     classify equation types.
C
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
C
C        When equation I is heavily weighted ITYPE(I)=0,
C        else ITYPE(I)=1.
C
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.D0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
C
C     Set the solution vector X(*) to zero and the column interchange
C     matrix to the identity.
C
      CALL DCOPY (N, 0.D0, 0, X, 1)
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
C
C     Perform initial triangularization in the submatrix
C     corresponding to the unconstrained variables.
C     Set first L components of dual vector to zero because
C     these correspond to the unconstrained variables.
C
      CALL DCOPY (L, 0.D0, 0, WD, 1)
C
C     The arrays IDOPE(*) and DOPE(*) are used to pass
C     information to DWNLIT().  This was done to avoid
C     a long calling sequence or the use of COMMON.
C
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     +            IDOPE, DOPE, DONE)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
C
C     Perform WNNLS algorithm using the following steps.
C
C     Until(DONE)
C        compute search direction and feasible point
C        when (HITCON) add constraints
C        else perform multiplier test and drop a constraint
C        fin
C     Compute-Final-Solution
C
C     To compute search direction and feasible point,
C     solve the triangular system of currently non-active
C     variables and store the solution in Z(*).
C
C     To solve system
C     Copy right hand side into TEMP vector to use overwriting method.
C
  160 IF (DONE) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL DCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL DAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
C
C     Increment iteration counter and check against maximum number
C     of iterations.
C
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         DONE = .TRUE.
      ENDIF
C
C     Check to see if any constraints have become active.
C     If so, calculate an interpolation factor so that all
C     active constraints are removed from the basis.
C
      ALPHA = 2.D0
      HITCON = .FALSE.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.D0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .TRUE.
         ENDIF
  180 CONTINUE
C
C     Compute search direction and feasible point
C
      IF (HITCON) THEN
C
C        To add constraints, use computed ALPHA to interpolate between
C        last feasible solution X(*) and current unconstrained (and
C        infeasible) solution Z(*).
C
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .FALSE.
C
C        Remove column JCON and shift columns JCON+1 through N to the
C        left.  Swap column JCON into the N th position.  This achieves
C        upper Hessenberg form for the nonactive constraints and
C        leaves an upper Hessenberg matrix to retriangularize.
C
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL DCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
C
C        Update permuted index vector to reflect this shift and swap.
C
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
C
C        Similarly permute X(*) vector.
C
         CALL DCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.D0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
C
C        Retriangularize upper Hessenberg matrix after adding
C        constraints.
C
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL DSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL DSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
C
C              Swapped row was formerly a pivot element, so it will
C              be large enough to perform elimination.
C              Zero IP1 to I in column J.
C
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
C
C                 Zero IP1 to I in column J
C
                  IF (W(I+1,J).NE.0.D0) THEN
                     CALL DROTMG (SCALE(I), SCALE(I+1), W(I,J),
     +                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.D0
                     CALL DROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                           SPARAM)
                  ENDIF
               ELSE
                  CALL DSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL DSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.D0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
C
C        See if the remaining coefficients in the solution set are
C        feasible.  They should be because of the way ALPHA was
C        determined.  If any are infeasible, it is due to roundoff
C        error.  Any that are non-positive will be set to zero and
C        removed from the solution set.
C
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.D0) GO TO 250
  240    CONTINUE
         FEASBL = .TRUE.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
C
C        To perform multiplier test and drop a constraint.
C
         CALL DCOPY (NSOLN, Z, 1, X, 1)
         IF (NSOLN.LT.N) CALL DCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
C
C        Reclassify least squares equations as equalities as necessary.
C
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL DSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL DSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
C
C        Form inner product vector WD(*) of dual coefficients.
C
         DO 280 J = NSOLN+1,N
            SM = 0.D0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
C
C        Find J such that WD(J)=WMAX is maximum.  This determines
C        that the incoming column J will reduce the residual vector
C        and be positive.
C
  290    WMAX = 0.D0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.D0) GO TO 330
C
C        Set dual coefficients to zero for incoming column.
C
         WD(IWMAX) = 0.D0
C
C        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set.
C        Perform transformation to retriangularize, and test for near
C        linear dependence.
C
C        Swap column IWMAX into NSOLN-th position to maintain upper
C        Hessenberg form of adjacent columns, and add new column to
C        triangular decomposition.
C
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL DSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.D0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
C
C        Reduce column NSOLN so that the matrix of nonactive constraints
C        variables is triangular.
C
         DO 320 J = M,NIV+1,-1
            JP = J - 1
C
C           When operating near the ME line, test to see if the pivot
C           element is near zero.  If so, use the largest element above
C           it as the pivot.  This is to maintain the sharp interface
C           between weighted and non-weighted rows in all cases.
C
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
C
            IF (W(J,NSOLN).NE.0.D0) THEN
               CALL DROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN),
     +                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.D0
               CALL DROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),
     +                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
C
C        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
C        this is nonpositive or too large.  If this was true or if the
C        pivot term was zero, reject the column as dependent.
C
         IF (W(NIV,NSOLN).NE.0.D0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.D0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
C
C           Try to add row ME+1 as an additional equality constraint.
C           Check size of proposed new solution component.
C           Reject it if it is too large.
C
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.D0) THEN
            ISOL = ME + 1
            IF (POS) THEN
C
C              Swap rows ME+1 and NIV, and scale factors for these rows.
C
               CALL DSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL DSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .FALSE.
         ENDIF
C
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.DONE)) GO TO 290
      ENDIF
      GO TO 160
C
C     Else perform multiplier test and drop a constraint.  To compute
C     final solution.  Solve system, store results in X(*).
C
C     Copy right hand side into TEMP vector to use overwriting method.
C
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL DCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL DAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
C
C     Solve system.
C
      CALL DCOPY (NSOLN, Z, 1, X, 1)
C
C     Apply Householder transformations to X(*) if KRANK.LT.L
C
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL DH12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
C
C     Fill in trailing zeroes for constrained variables not in solution.
C
      IF (NSOLN.LT.N) CALL DCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
C
C     Permute solution vector to natural order.
C
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
C
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL DSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
C
C     Rescale the solution using the column scaling.
C
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
C
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
C
      RNORM = SQRT(RNORM)
      RETURN
      END
*DECK DWNLIT
      SUBROUTINE DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE,
     +   RNORM, IDOPE, DOPE, DONE)
C***BEGIN PROLOGUE  DWNLIT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DWNNLS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to DWNNLS( ).
C     The documentation for DWNNLS( ) has complete usage instructions.
C
C     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
C           B as the (N+1)st col.
C
C     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
C     col interchanges.
C
C***SEE ALSO  DWNNLS
C***ROUTINES CALLED  DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1,
C                    DWNLT2, DWNLT3, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900604  DP version created from SP version. .  (RWC)
C***END PROLOGUE  DWNLIT
      INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
      DOUBLE PRECISION DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
      LOGICAL DONE
C
      EXTERNAL DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1,
     *   DWNLT2, DWNLT3, IDAMAX
      INTEGER IDAMAX
      LOGICAL DWNLT2
C
      DOUBLE PRECISION ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5),
     *   T, TAU
      INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME,
     *   MEND, NIV, NSOLN
      LOGICAL INDEP, RECALC
C
C***FIRST EXECUTABLE STATEMENT  DWNLIT
      ME    = IDOPE(1)
      NSOLN = IDOPE(2)
      L1    = IDOPE(3)
C
      ALSQ   = DOPE(1)
      EANORM = DOPE(2)
      TAU    = DOPE(3)
C
      LB     = MIN(M-1,L)
      RECALC = .TRUE.
      RNORM  = 0.D0
      KRANK  = 0
C
C     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
C     included in the test for column independence.
C
      FACTOR = 1.D0
      LEND = L
      DO 180 I=1,LB
C
C        Set IR to point to the I-th row.
C
         IR = I
         MEND = M
         CALL DWNLT1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                W)
C
C        Update column SS and find pivot column.
C
         CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C        Perform column interchange.
C        Test independence of incoming column.
C
  130    IF (DWNLT2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
C
C           Eliminate I-th column below diagonal using modified Givens
C           transformations applied to (A B).
C
C           When operating near the ME line, use the largest element
C           above it as the pivot.
C
            DO 160 J=M,I+1,-1
               JP = J-1
               IF (J.EQ.ME+1) THEN
                  IMAX = ME
                  AMAX = SCALE(ME)*W(ME,I)**2
                  DO 150 JP=J-1,I,-1
                     T = SCALE(JP)*W(JP,I)**2
                     IF (T.GT.AMAX) THEN
                        IMAX = JP
                        AMAX = T
                     ENDIF
  150             CONTINUE
                  JP = IMAX
               ENDIF
C
               IF (W(J,I).NE.0.D0) THEN
                  CALL DROTMG (SCALE(JP), SCALE(J), W(JP,I), W(J,I),
     +                         SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  160       CONTINUE
         ELSE IF (LEND.GT.I) THEN
C
C           Column I is dependent.  Swap with column LEND.
C           Perform column interchange,
C           and find column in remaining set with largest SS.
C
            CALL DWNLT3 (I, LEND, M, MDW, IPIVOT, H, W)
            LEND = LEND - 1
            IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
            HBAR = H(IMAX)
            GO TO 130
         ELSE
            KRANK = I - 1
            GO TO 190
         ENDIF
  180 CONTINUE
      KRANK = L1
C
  190 IF (KRANK.LT.ME) THEN
         FACTOR = ALSQ
         DO 200 I=KRANK+1,ME
            CALL DCOPY (L, 0.D0, 0, W(I,1), MDW)
  200    CONTINUE
C
C        Determine the rank of the remaining equality constraint
C        equations by eliminating within the block of constrained
C        variables.  Remove any redundant constraints.
C
         RECALC = .TRUE.
         LB = MIN(L+ME-KRANK, N)
         DO 270 I=L+1,LB
            IR = KRANK + I - L
            LEND = N
            MEND = ME
            CALL DWNLT1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H,
     +                   SCALE, W)
C
C           Update col ss and find pivot col
C
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange
C           Eliminate elements in the I-th col.
C
            DO 240 J=ME,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL DROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  240       CONTINUE
C
C           I=column being eliminated.
C           Test independence of incoming column.
C           Remove any redundant or dependent equality constraints.
C
            IF (.NOT.DWNLT2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
               JJ = IR
               DO 260 IR=JJ,ME
                  CALL DCOPY (N, 0.D0, 0, W(IR,1), MDW)
                  RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
                  W(IR,N+1) = 0.D0
                  SCALE(IR) = 1.D0
C
C                 Reclassify the zeroed row as a least squares equation.
C
                  ITYPE(IR) = 1
  260          CONTINUE
C
C              Reduce ME to reflect any discovered dependent equality
C              constraints.
C
               ME = JJ - 1
               GO TO 280
            ENDIF
  270    CONTINUE
      ENDIF
C
C     Try to determine the variables KRANK+1 through L1 from the
C     least squares equations.  Continue the triangularization with
C     pivot element W(ME+1,I).
C
  280 IF (KRANK.LT.L1) THEN
         RECALC = .TRUE.
C
C        Set FACTOR=ALSQ to remove effect of heavy weight from
C        test for column independence.
C
         FACTOR = ALSQ
         DO 350 I=KRANK+1,L1
C
C           Set IR to point to the ME+1-st row.
C
            IR = ME+1
            LEND = L
            MEND = M
            CALL DWNLT1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                   W)
C
C           Update column SS and find pivot column.
C
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange.
C           Eliminate I-th column below the IR-th element.
C
            DO 320 J=M,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL DROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.D0
                  CALL DROTM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  320       CONTINUE
C
C           Test if new pivot element is near zero.
C           If so, the column is dependent.
C           Then check row norm test to be classified as independent.
C
            T = SCALE(IR)*W(IR,I)**2
            INDEP = T .GT. (TAU*EANORM)**2
            IF (INDEP) THEN
               RN = 0.D0
               DO 340 I1=IR,M
                  DO 330 J1=I+1,N
                     RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
               INDEP = T .GT. RN*TAU**2
            ENDIF
C
C           If independent, swap the IR-th and KRANK+1-th rows to
C           maintain the triangular form.  Update the rank indicator
C           KRANK and the equality constraint pointer ME.
C
            IF (.NOT.INDEP) GO TO 360
            CALL DSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
            CALL DSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
C
C           Reclassify the least square equation as an equality
C           constraint and rescale it.
C
            ITYPE(IR) = 0
            T = SQRT(SCALE(KRANK+1))
            CALL DSCAL(N+1, T, W(KRANK+1,1), MDW)
            SCALE(KRANK+1) = ALSQ
            ME = ME+1
            KRANK = KRANK+1
  350    CONTINUE
      ENDIF
C
C     If pseudorank is less than L, apply Householder transformation.
C     from right.
C
  360 IF (KRANK.LT.L) THEN
         DO 370 J=KRANK,1,-1
            CALL DH12 (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1,
     +                J-1)
  370    CONTINUE
      ENDIF
C
      NIV = KRANK + NSOLN - L
      IF (L.EQ.N) DONE = .TRUE.
C
C     End of initial triangularization.
C
      IDOPE(1) = ME
      IDOPE(2) = KRANK
      IDOPE(3) = NIV
      RETURN
      END
*DECK DWNLT1
      SUBROUTINE DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,
     +   SCALE, W)
C***BEGIN PROLOGUE  DWNLT1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To update the column Sum Of Squares and find the pivot column.
C     The column Sum of Squares Vector will be updated at each step.
C     When numerically necessary, these values will be recomputed.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT1
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
C
      EXTERNAL IDAMAX
      INTEGER IDAMAX
C
      INTEGER J, K
C
C***FIRST EXECUTABLE STATEMENT  DWNLT1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
C
C        Update column SS=sum of squares.
C
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
C
C        Test for numerical accuracy.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
C
C     If required, recalculate column SS, using rows IR through MEND.
C
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.D0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
C
C        Find column with largest SS.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END
*DECK DWNLT2
      LOGICAL FUNCTION DWNLT2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
C***BEGIN PROLOGUE  DWNLT2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To test independence of incoming column.
C
C     Test the column IC to determine if it is linearly independent
C     of the columns already in the basis.  In the initial tri. step,
C     we usually want the heavy weight ALAMDA to be included in the
C     test for independence.  In this case, the value of FACTOR will
C     have been set to 1.E0 before this procedure is invoked.
C     In the potentially rank deficient problem, the value of FACTOR
C     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
C     heavy weight from the test for independence.
C
C     Write new column as partitioned vector
C           (A1)  number of components in solution so far = NIV
C           (A2)  M-NIV components
C     And compute  SN = inverse weighted length of A1
C                  RN = inverse weighted length of A2
C     Call the column independent when RN .GT. TAU*SN
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT2
      DOUBLE PRECISION FACTOR, SCALE(*), TAU, WIC(*)
      INTEGER IR, ME, MEND
C
      DOUBLE PRECISION RN, SN, T
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT  DWNLT2
      SN = 0.E0
      RN = 0.E0
      DO 10 J=1,MEND
         T = SCALE(J)
         IF (J.LE.ME) T = T/FACTOR
         T = T*WIC(J)**2
C
         IF (J.LT.IR) THEN
            SN = SN + T
         ELSE
            RN = RN + T
         ENDIF
   10 CONTINUE
      DWNLT2 = RN .GT. SN*TAU**2
      RETURN
      END
*DECK DWNLT3
      SUBROUTINE DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C***BEGIN PROLOGUE  DWNLT3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Perform column interchange.
C     Exchange elements of permuted index vector and perform column
C     interchanges.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  DSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      DOUBLE PRECISION H(*), W(MDW,*)
C
      EXTERNAL DSWAP
C
      DOUBLE PRECISION T
      INTEGER ITEMP
C
C***FIRST EXECUTABLE STATEMENT  DWNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
C
         CALL DSWAP(M, W(1,IMAX), 1, W(1,I), 1)
C
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END
