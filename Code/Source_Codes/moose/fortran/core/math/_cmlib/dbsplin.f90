!===============================================================================
!
!***BEGIN PROLOGUE  BSPDOC
!***DATE WRITTEN   810223   (YYMMDD)
!***REVISION DATE  840425   (YYMMDD)
!***CATEGORY NO.  E,K
!***KEYWORDS  B-SPLINES,DOCUMENTATION,SPLINES
!***AUTHOR  AMOS, D. E., (SNLA)
!***PURPOSE  Comments on B-spline routines
!***DESCRIPTION
!
!     Written by D. E. Amos, May 1980
!
!     References
!         1. Computation with Splines and B-Splines  by D. E.
!            Amos, SAND78-1968, March, 1979.
!         2. Quadrature Subroutines for Splines and B-Splines by
!            D. E. Amos, SAND79-1825, December, 1979.
!         3. A Practical Guide to Splines  by C. de Boor,
!            Applied Math Sci 27, Springer, N.Y., 1978.
!         4. On Calculating with B-Splines  by C. de Boor, J.
!            Approx. Theory,6,50-62(1972)
!         5. Package for Calculating with B-Splines  by C. De
!            Boor, SIAM J. Numer. anal.,14,441-472(1977).
!         6. Constrained Least Squares Curve Fitting to Discrete
!            Data Using B-Splines - A User's Guide  by R. J.
!            Hanson, SAND78-1291, February,1979.
!         7. Monotone Piecewise Cubic Interpolation by F. N. Fritsch
!            and R. E. Carlson, LLNL report UCRL-82453, January, 1979.
!
!     Abstract
!         BSPDOC is a non-executable, B-spline documentary routine.
!         The narrative describes a B-spline and the routines
!         necessary to manipulate B-splines at a fairly high level.
!         The basic package described herein is that of reference
!         5 with names altered to prevent duplication and conflicts
!         with routines from reference 3.  The call lists used here
!         are also different.  Work vectors were added to ensure
!         portability and proper execution in an overlay environ-
!         ment.  These work arrays can be used for other purposes
!         except as noted in BSPVN.  While most of the original
!         routines in reference 5 were restricted to orders 20
!         or less, this restriction was removed from all routines
!         except the quadrature routine BSQAD.  (See the section
!         below on differentiation and integration for details.)
!
!         The subroutines referenced below are single precision
!         routines.  Corresponding double precision versions are also
!         part of the package, and these are referenced by prefixing
!         a D in front of the single precision name.  For example,
!         BVALU and DBVALU are the single and double precision
!         versions for evaluating a B-spline or any of its deriva-
!         tives in the B-representation.
!
!                ****Description of B-Splines****
!
!     A collection of polynomials of fixed degree K-1 defined on a
!     subdivision (X(I),X(I+1)), I=1,...,M-1 of (A,B) with X(1)=A,
!     X(M)=B is called a B-spline of order K.  If the spline has K-2
!     continuous derivatives on (A,B), then the B-spline is simply
!     called a spline of order K.  Each of the M-1 polynomial pieces
!     has K coefficients, making a total of K(M-1) parameters.  This
!     B-spline and its derivatives have M-2 jumps at the subdivision
!     points X(I), I=2,...,M-1.  Continuity requirements at these
!     subdivision points add constraints and reduce the number of free
!     parameters.  If a B-spline is continuous at each of the M-2 sub-
!     division points, there are K(M-1)-(M-2) free parameters; if in
!     addition the B-spline has continuous first derivatives, there
!     are K(M-1)-2(M-2) free parameters, etc., until we get to a
!     spline where we have K(M-1)-(K-1)(M-2) = M+K-2 free parameters.
!     Thus, the principle is that increasing the continuity of
!     derivatives decreases the number of free parameters and
!     conversely.
!
!     The points at which the polynomials are tied together by the
!     continuity conditions are called knots.  If two knots are
!     allowed to come together at some X(I), then we say that we
!     have a knot of multiplicity 2 there, and the knot values are
!     the X(I) value.  If we reverse the procedure of the first
!     paragraph, we find that adding a knot to increase multiplicity
!     increases the number of free parameters and, according to the
!     principle above, we thereby introduce a discontinuity in what
!     was the highest continuous derivative at that knot.  Thus, the
!     number of free parameters is N = NU+K-2 where NU is the sum
!     of multipicities at the X(I) values with X(1) and X(M) of
!     multiplicity 1 (NU = M if all knots are simple, i.e., for a
!     spline, all knots have multiplicity 1.)  Each knot can have a
!     multiplicity of at most K.  A B-spline is commonly written in the
!     B-representation
!
!               Y(X) = sum( A(I)*B(I,X), I=1 , N)
!
!     to show the explicit dependence of the spline on the free
!     parameters or coefficients A(I)=BCOEF(I) and basis functions
!     B(I,X).  These basis functions are themselves special B-splines
!     which are zero except on (at most) K adjoining intervals where
!     each B(I,X) is positive and, in most cases, hat or bell-
!     shaped.  In order for the nonzero part of B(1,X) to be a spline
!     covering (X(1),X(2)), it is necessary to put K-1 knots to the
!     left of A and similarly for B(N,X) to the right of B.  Thus, the
!     total number of knots for this representation is NU+2K-2 = N+K.
!     These knots are carried in an array T(*) dimensioned by at least
!     N+K.  From the construction, A=T(K) and B=T(N+1) and the spline is
!     defined on T(K).LE.X.LE.T(N+1).  The nonzero part of each basis
!     function lies in the  Interval (T(I),T(I+K)).  In many problems
!     where extrapolation beyond A or B is not anticipated, it is common
!     practice to set T(1)=T(2)=...=T(K)=A and T(N+1)=T(N+2)=...=
!     T(N+K)=B.  In summary, since T(K) and T(N+1) as well as
!     interior knots can have multiplicity K, the number of free
!     parameters N = sum of multiplicties - K.  The fact that each
!     B(I,X) function is nonzero over at most K intervals means that
!     for a given X value, there are at most K nonzero terms of the
!     sum.  This leads to banded matrices in linear algebra problems,
!     and references 3 and 6 take advantage of this in con-
!     structing higher level routines to achieve speed and avoid
!     ill-conditioning.
!
!                     ****Basic Routines****
!
!     The basic routines which most casual users will need are those
!     concerned with direct evaluation of splines or B-splines.
!     Since the B-representation, denoted by (T,BCOEF,N,K), is
!     preferred because of numerical stability, the knots T(*), the
!     B-spline coefficients BCOEF(*), the number of coefficients N,
!     and the order K of the polynomial pieces (of degree K-1) are
!     usually given.  While the knot array runs from T(1) to T(N+K),
!     the B-spline is normally defined on the interval T(K).LE.X.LE.
!     T(N+1).  To evaluate the B-spline or any of its derivatives
!     on this interval, one can use
!
!                  Y = BVALU(T,BCOEF,N,K,ID,X,INBV,WORK)
!
!     where ID is an integer for the ID-th derivative, 0.LE.ID.LE.K-1.
!     ID=0 gives the zero-th derivative or B-spline value at X.
!     If X.LT.T(K) or X.GT.T(N+1), whether by mistake or the result
!     of round off accumulation in incrementing X, BVALU gives a
!     diagnostic.  INBV is an initialization parameter which is set
!     to 1 on the first call.  Distinct splines require distinct
!     INBV parameters.  WORK is a scratch vector of length at least
!     3*K.
!
!     When more conventional communication is needed for publication,
!     physical interpretation, etc., the B-spline coefficients can
!     be converted to piecewise polynomial (PP) coefficients.  Thus,
!     the breakpoints (distinct knots) XI(*), the number of
!     polynomial pieces LXI, and the (right) derivatives C(*,J) at
!     each breakpoint XI(J) are needed to define the Taylor
!     expansion to the right of XI(J) on each interval XI(J).LE.
!     X.LT.XI(J+1), J=1,LXI where XI(1)=A and XI(LXI+1)=B.
!     These are obtained from the (T,BCOEF,N,K) representation by
!
!                CALL BSPPP(T,BCOEF,N,K,LDC,C,XI,LXI,WORK)
!
!     where LDC.GE.K is the leading dimension of the matrix C and
!     WORK is a scratch vector of length at least K*(N+3).
!     Then the PP-representation (C,XI,LXI,K) of Y(X), denoted
!     by Y(J,X) on each interval XI(J).LE.X.LT.XI(J+1), is
!
!     Y(J,X) = sum( C(I,J)*((X-XI(J))**(I-1))/factorial(I-1), I=1,K)
!
!     for J=1,...,LXI.  One must view this conversion from the B-
!     to the PP-representation with some skepticism because the
!     conversion may lose significant digits when the B-spline
!     varies in an almost discontinuous fashion.  To evaluate
!     the B-spline or any of its derivatives using the PP-
!     representation, one uses
!
!                Y = PPVAL(LDC,C,XI,LXI,K,ID,X,INPPV)
!
!     where ID and INPPV have the same meaning and useage as ID and
!     INBV in BVALU.
!
!     To determine to what extent the conversion process loses
!     digits, compute the relative error ABS((Y1-Y2)/Y2) over
!     the X interval with Y1 from PPVAL and Y2 from BVALU.  A
!     major reason for considering PPVAL is that evaluation is
!     much faster than that from BVALU.
!
!     Recall that when multiple knots are encountered, jump type
!     discontinuities in the B-spline or its derivatives occur
!     at these knots, and we need to know that BVALU and PPVAL
!     return right limiting values at these knots except at
!     X=B where left limiting values are returned.  These values
!     are used for the Taylor expansions about left end points of
!     breakpoint intervals.  That is, the derivatives C(*,J) are
!     right derivatives.  Note also that a computed X value which,
!     mathematically, would be a knot value may differ from the knot
!     by a round off error.  When this happens in evaluating a dis-
!     continuous B-spline or some discontinuous derivative, the
!     value at the knot and the value at X can be radically
!     different.  In this case, setting X to a T or XI value makes
!     the computation precise.  For left limiting values at knots
!     other than X=B, see the prologues to BVALU and other
!     routines.
!
!                     ****Interpolation****
!
!     BINTK is used to generate B-spline parameters (T,BCOEF,N,K)
!     which will interpolate the data by calls to BVALU.  A similar
!     interpolation can also be done for cubic splines using BINT4
!     or the code in reference 7.  If the PP-representation is given,
!     one can evaluate this representation at an appropriate number of
!     abscissas to create data then use BINTK or BINT4 to generate
!     the B-representation.
!
!               ****Differentiation and Integration****
!
!     Derivatives of B-splines are obtained from BVALU or PPVAL.
!     Integrals are obtained from BSQAD using the B-representation
!     (T,BCOEF,N,K) and PPQAD using the PP-representation (C,XI,LXI,
!     K).  More compicated integrals involving the product of a
!     of a function F and some derivative of a B-spline can be
!     evaluated with BFQAD or PFQAD using the B- or PP- represen-
!     tations respectively.  All quadrature routines, except for PPQAD,
!     are limited in accuracy to 18 digits or working precision,
!     whichever is smaller.  PPQAD is limited to working precision
!     only.  In addition, the order K for BSQAD is limited to 20 or
!     less.  If orders greater than 20 are required, use BFQAD with
!     F(X) = 1.
!
!                      ****Extrapolation****
!
!     Extrapolation outside the interval (A,B) can be accomplished
!     easily by the PP-representation using PPVAL.  However,
!     caution should be exercised, especially when several knots
!     are located at A or B or when the extrapolation is carried
!     significantly beyond A or B.  On the other hand, direct
!     evaluation with BVALU outside A=T(K).LE.X.LE.T(N+1)=B
!     produces an error message, and some manipulation of the knots
!     and coefficients are needed to extrapolate with BVALU.  This
!     process is described in reference 6.
!
!                ****Curve Fitting and Smoothing****
!
!     Unless one has many accurate data points, direct inter-
!     polation is not recommended for summarizing data.  The
!     results are often not in accordance with intuition since the
!     fitted curve tends to oscillate through the set of points.
!     Monotone splines (reference 7) can help curb this undulating
!     tendency but constrained least squares is more likely to give an
!     acceptable fit with fewer parameters.  Subroutine FC, des-
!     cribed in reference 6, is recommended for this purpose.  The
!     output from this fitting process is the B-representation.
!
!              **** Routines in the B-Spline Package ****
!
!                      Single Precision Routines
!
!         The subroutines referenced below are SINGLE PRECISION
!         routines. Corresponding DOUBLE PRECISION versions are also
!         part of the package and these are referenced by prefixing
!         a D in front of the single precision name. For example,
!         BVALU and DBVALU are the SINGLE and DOUBLE PRECISION
!         versions for evaluating a B-spline or any of its deriva-
!         tives in the B-representation.
!
!     BINT4 - interpolates with splines of order 4
!     BINTK - interpolates with splines of order k
!     BSQAD - integrates the B-representation on subintervals
!     PPQAD - integrates the PP-representation
!     BFQAD - integrates the product of a function F and any spline
!             derivative in the B-representation
!     PFQAD - integrates the product of a function F and any spline
!             derivative in the PP-representation
!     BVALU - evaluates the B-representation or a derivative
!     PPVAL - evaluates the PP-representation or a derivative
!     INTRV - gets the largest index of the knot to the left of x
!     BSPPP - converts from B- to PP-representation
!     BSPVD - computes nonzero basis functions and derivatives at x
!     BSPDR - sets up difference array for BSPEV
!     BSPEV - evaluates the B-representation and derivatives
!     BSPVN - called by BSPEV, BSPVD, BSPPP and BINTK for function and
!             derivative evaluations
!                        Auxiliary Routines
!
!       BSGQ8,PPGQ8,BNSLV,BNFAC,XERROR,DBSGQ8,DPPGQ8,DBNSLV,DBNFAC
!
!                    Machine Dependent Routines
!
!                      I1MACH, R1MACH, D1MACH
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  BSPDOC
!
!===============================================================================
MODULE MOOSE_CMLIB_DBSPLIN
  USE MOOSE_CMLIB_XERROR

  CONTAINS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBINTK(X,Y,T,N,K,BCOEF,Q,WORK)
  !***BEGIN PROLOGUE  DBINTK
  !***DATE WRITTEN   800901   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E1A
  !***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
  !         SPLINE
  !***AUTHOR  AMOS, D. E., (SNLA)
  !***PURPOSE  Produces the B-spline coefficients, BCOEF, of the
  !        B-spline of order K with knots T(I), I=1,...,N+K, which
  !        takes on the value Y(I) at X(I), I=1,...,N.
  !***DESCRIPTION
  !
  ! Written by Carl de Boor and modified by D. E. Amos
  !
  ! References
  !
  !     A Practical Guide to Splines by C. de Boor, Applied
  !     Mathematics Series 27, Springer, 1979.
  !
  ! Abstract    **** a double precision routine ****
  !
  !     DBINTK is the SPLINT routine of the reference.
  !
  !     DBINTK produces the B-spline coefficients, BCOEF, of the
  !     B-spline of order K with knots T(I), I=1,...,N+K, which
  !     takes on the value Y(I) at X(I), I=1,...,N.  The spline or
  !     any of its derivatives can be evaluated by calls to DBVALU.
  !
  !     The I-th equation of the linear system A*BCOEF = B for the
  !     coefficients of the interpolant enforces interpolation at
  !     X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
  !     a band matrix with 2K-1 bands if A is invertible.  The matrix
  !     A is generated row by row and stored, diagonal by diagonal,
  !     in the rows of Q, with the main diagonal going into row K.
  !     The banded system is then solved by a call to DBNFAC (which
  !     constructs the triangular factorization for A and stores it
  !     again in Q), followed by a call to DBNSLV (which then
  !     obtains the solution BCOEF by substitution).  DBNFAC does no
  !     pivoting, since the total positivity of the matrix A makes
  !     this unnecessary.  The linear system to be solved is
  !     (theoretically) invertible if and only if
  !             T(I) .LT. X(I) .LT. T(I+K),        for all I.
  !     Equality is permitted on the left for I=1 and on the right
  !     for I=N when K knots are used at X(1) or X(N).  Otherwise,
  !     violation of this condition is certain to lead to an error.
  !
  !     DBINTK calls DBSPVN, DBNFAC, DBNSLV, XERROR
  !
  ! Description of Arguments
  !
  !     Input       X,Y,T are double precision
  !       X       - vector of length N containing data point abscissa
  !                 in strictly increasing order.
  !       Y       - corresponding vector of length N containing data
  !                 point ordinates.
  !       T       - knot vector of length N+K
  !                 Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
  !                 .GE. X(N), this leaves only N-K knots (not nec-
  !                 essarily X(I) values) interior to (X(1),X(N))
  !       N       - number of data points, N .GE. K
  !       K       - order of the spline, K .GE. 1
  !
  !     Output      BCOEF,Q,WORK are double precision
  !       BCOEF   - a vector of length N containing the B-spline
  !                 coefficients
  !       Q       - a work vector of length (2*K-1)*N, containing
  !                 the triangular factorization of the coefficient
  !                 matrix of the linear system being solved.  The
  !                 coefficients for the interpolant of an
  !                 additional data set (X(I),yY(I)), I=1,...,N
  !                 with the same abscissa can be obtained by loading
  !                 YY into BCOEF and then executing
  !                     CALL DBNSLV(Q,2K-1,N,K-1,K-1,BCOEF)
  !       WORK    - work vector of length 2*K
  !
  ! Error Conditions
  !     Improper input is a fatal error
  !     Singular system of equations is a fatal error
  !***REFERENCES  C. DE BOOR, *A PRACTICAL GUIDE TO SPLINES*, APPLIED
  !             MATHEMATICS SERIES 27, SPRINGER, 1979.
  !           D.E. AMOS, *COMPUTATION WITH SPLINES AND B-SPLINES*,
  !             SAND78-1968,SANDIA LABORATORIES,MARCH,1979.
  !***ROUTINES CALLED  DBNFAC,DBNSLV,DBSPVN,XERROR
  !***END PROLOGUE  DBINTK
  !
  !
  INTEGER :: IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT, &
        LENQ, NP1
  DOUBLE PRECISION :: BCOEF(N), Y(N), Q(*), T(*), X(N), XI, WORK(*)
  ! DIMENSION Q(2*K-1,N), T(N+K)
  !***FIRST EXECUTABLE STATEMENT  DBINTK
  IF(K.LT.1) GO TO 100
  IF(N.LT.K) GO TO 105
  JJ = N - 1
  IF(JJ.EQ.0) GO TO 6
  DO I=1,JJ
  IF(X(I).GE.X(I+1)) GO TO 110
  END DO
    6   CONTINUE
  NP1 = N + 1
  KM1 = K - 1
  KPKM2 = 2*KM1
  LEFT = K
             ! ZERO OUT ALL ENTRIES OF Q
  LENQ = N*(K+KM1)
  DO I=1,LENQ
    Q(I) = 0.0D0
  END DO
  !
  !  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
  DO I=1,N
    XI = X(I)
    ILP1MX = MIN0(I+K,NP1)
     ! *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
     !         T(LEFT) .LE. X(I) .LT. T(LEFT+1)
     ! MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
    LEFT = MAX0(LEFT,I)
    IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
    LEFT = LEFT + 1
    IF (LEFT.LT.ILP1MX) GO TO 20
    LEFT = LEFT - 1
    IF (XI.GT.T(LEFT+1)) GO TO 80
     ! *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
     ! A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
     ! LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
     ! ARE RETURNED, IN  BCOEF (USED FOR TEMP.STORAGE HERE), BY THE
     ! FOLLOWING
   30   CALL DBSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
     ! WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
     ! A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
     ! A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
     ! AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
     ! DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
     ! ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
     ! DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
     ! ENTRY
     !     I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
     !            =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
     ! OF  Q .
    JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
    DO J=1,K
      JJ = JJ + KPKM2
      Q(JJ) = BCOEF(J)
    END DO
  END DO
  !
  ! ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
  CALL DBNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
  GO TO (60, 90), IFLAG
  ! *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60   DO 70 I=1,N
    BCOEF(I) = Y(I)
   70   CONTINUE
  CALL DBNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
  RETURN
  !
  !
  80   CONTINUE
  CALL XERROR( ' DBINTK,  SOME ABSCISSA WAS NOT IN THE SUPPORT OF TH&
        &E CORRESPONDING BASIS FUNCTION AND THE SYSTEM IS SINGULAR.',109,2, &
          1)
  RETURN
  90     CONTINUE
  CALL XERROR( ' DBINTK,  THE SYSTEM OF SOLVER DETECTS A SINGULAR SY&
          &STEM ALTHOUGH THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATIS&
          &FIED.',123,8,1)
  RETURN
  100     CONTINUE
  CALL XERROR( ' DBINTK,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
  RETURN
  105     CONTINUE
  CALL XERROR( ' DBINTK,  N DOES NOT SATISFY N.GE.K', 35, 2, 1)
  RETURN
  110     CONTINUE
  CALL XERROR( ' DBINTK,  X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', 57, 2, 1)
  RETURN
  END SUBROUTINE DBINTK
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBNFAC(W,NROWW,NROW,NBANDL,NBANDU,IFLAG)
  !***BEGIN PROLOGUE  DBNFAC
  !***REFER TO  DBINT4,DBINTK
  !
  !  DBNFAC is the BANFAC routine from
  !    * A Practical Guide to Splines *  by C. de Boor
  !
  !  DBNFAC is a double precision routine
  !
  !  Returns in  W  the LU-factorization (without pivoting) of the banded
  !  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
  !  onals in the work array  W .
  !
  ! *****  I N P U T  ****** W is double precision
  !  W.....Work array of size  (NROWW,NROW)  containing the interesting
  !    part of a banded matrix  A , with the diagonals or bands of  A
  !    stored in the rows of  W , while columns of  A  correspond to
  !    columns of  W . This is the storage mode used in  LINPACK  and
  !    results in efficient innermost loops.
  !       Explicitly,  A  has  NBANDL  bands below the diagonal
  !                        +     1     (main) diagonal
  !                        +   NBANDU  bands above the diagonal
  !    and thus, with    MIDDLE = NBANDU + 1,
  !      A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
  !                                          J=1,...,NROW .
  !    For example, the interesting entries of A (1,2)-banded matrix
  !    of order  9  would appear in the first  1+1+2 = 4  rows of  W
  !    as follows.
  !                      13 24 35 46 57 68 79
  !                   12 23 34 45 56 67 78 89
  !                11 22 33 44 55 66 77 88 99
  !                21 32 43 54 65 76 87 98
  !
  !    All other entries of  W  not identified in this way with an en-
  !    try of  A  are never referenced .
  !  NROWW.....Row dimension of the work array  W .
  !    must be  .GE.  NBANDL + 1 + NBANDU  .
  !  NBANDL.....Number of bands of  A  below the main diagonal
  !  NBANDU.....Number of bands of  A  above the main diagonal .
  !
  ! *****  O U T P U T  ****** W is double precision
  !  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
  ! If  IFLAG = 1, then
  !  W.....contains the LU-factorization of  A  into a unit lower triangu-
  !    lar matrix  L  and an upper triangular matrix  U (both banded)
  !    and stored in customary fashion over the corresponding entries
  !    of  A . This makes it possible to solve any particular linear
  !    system  A*X = B  for  X  by a
  !          CALL DBNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
  !    with the solution X  contained in  B  on return .
  ! If  IFLAG = 2, then
  !    one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
  !    one of the potential pivots was found to be zero indicating
  !    that  A  does not have an LU-factorization. This implies that
  !    A  is singular in case it is totally positive .
  !
  ! *****  M E T H O D  ******
  ! Gauss elimination  W I T H O U T  pivoting is used. The routine is
  !  intended for use with matrices  A  which do not require row inter-
  !  changes during factorization, especially for the  T O T A L L Y
  !  P O S I T I V E  matrices which occur in spline calculations.
  ! The routine should NOT be used for an arbitrary banded matrix.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  DBNFAC
  !
  INTEGER :: IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K, &
        KMAX, MIDDLE, MIDMK, NROWM1
  DOUBLE PRECISION :: W(NROWW,NROW), FACTOR, PIVOT
  !
  !***FIRST EXECUTABLE STATEMENT  DBNFAC
  IFLAG = 1
  MIDDLE = NBANDU + 1
                      ! W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
  NROWM1 = NROW - 1
  IF (NROWM1 < 0) THEN
     GO TO 120
  ELSEIF (NROWM1 == 0) THEN
     GO TO 110
  ENDIF

   10   IF (NBANDL.GT.0) GO TO 30
             ! A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
  DO I=1,NROWM1
    IF (W(MIDDLE,I).EQ.0.0D0) GO TO 120
  END DO
  GO TO 110
   30   IF (NBANDU.GT.0) GO TO 60
           ! A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
           !    DIVIDE EACH COLUMN BY ITS DIAGONAL .
  DO I=1,NROWM1
    PIVOT = W(MIDDLE,I)
    IF (PIVOT.EQ.0.0D0) GO TO 120
    JMAX = MIN0(NBANDL,NROW-I)
    DO J=1,JMAX
      W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
    END DO
  END DO
  RETURN
  !
  !    A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60   DO 100 I=1,NROWM1
                               ! W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
    PIVOT = W(MIDDLE,I)
    IF (PIVOT.EQ.0.0D0) GO TO 120
              ! JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
              !     BELOW THE DIAGONAL .
    JMAX = MIN0(NBANDL,NROW-I)
           ! DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
    DO J=1,JMAX
      W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
    END DO
              ! KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
              !     THE RIGHT OF THE DIAGONAL .
    KMAX = MIN0(NBANDU,NROW-I)
               ! SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
               ! (BELOW ROW  I ) .
    DO K=1,KMAX
      IPK = I + K
      MIDMK = MIDDLE - K
      FACTOR = W(MIDMK,IPK)
      DO J=1,JMAX
        W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
      END DO
    END DO
  100   CONTINUE
                                    ! CHECK THE LAST DIAGONAL ENTRY .
  110   IF (W(MIDDLE,NROW).NE.0.0D0) RETURN
  120   IFLAG = 2
  RETURN
  END SUBROUTINE DBNFAC
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBNSLV(W,NROWW,NROW,NBANDL,NBANDU,B)
  !***BEGIN PROLOGUE  DBNSLV
  !***REFER TO  DBINT4,DBINTK
  !
  !  DBNSLV is the BANSLV routine from
  !    * A Practical Guide to Splines *  by C. de Boor
  !
  !  DBNSLV is a double precision routine
  !
  !  Companion routine to  DBNFAC . It returns the solution  X  of the
  !  linear system  A*X = B  in place of  B , given the LU-factorization
  !  for  A  in the work array  W from DBNFAC.
  !
  ! *****  I N P U T  ****** W,B are DOUBLE PRECISION
  !  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
  !    banded matrix  A  of order  NROW  as constructed in  DBNFAC .
  !    For details, see  DBNFAC .
  !  B.....Right side of the system to be solved .
  !
  ! *****  O U T P U T  ****** B is DOUBLE PRECISION
  !  B.....Contains the solution  X , of order  NROW .
  !
  ! *****  M E T H O D  ******
  ! (With  A = L*U, as stored in  W,) the unit lower triangular system
  !  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
  !  upper triangular system  U*X = Y  is solved for  X  . The calcul-
  !  ations are so arranged that the innermost loops stay within columns.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  DBNSLV
  !
  INTEGER :: NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
  DOUBLE PRECISION :: W(NROWW,NROW), B(NROW)
  !***FIRST EXECUTABLE STATEMENT  DBNSLV
  MIDDLE = NBANDU + 1
  IF (NROW.EQ.1) GO TO 80
  NROWM1 = NROW - 1
  IF (NBANDL.EQ.0) GO TO 30
                              ! FORWARD PASS
         ! FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
         ! OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
  DO I=1,NROWM1
    JMAX = MIN0(NBANDL,NROW-I)
    DO J=1,JMAX
      B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
    END DO
  END DO
                              ! BACKWARD PASS
         ! FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
         ! ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
         ! OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30   IF (NBANDU.GT.0) GO TO 50
                             ! A  IS LOWER TRIANGULAR .
  DO I=1,NROW
    B(I) = B(I)/W(1,I)
  END DO
  RETURN
   50   I = NROW
   60   B(I) = B(I)/W(MIDDLE,I)
  JMAX = MIN0(NBANDU,I-1)
  DO J=1,JMAX
    B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
  END DO
  I = I - 1
  IF (I.GT.1) GO TO 60
   80   B(1) = B(1)/W(MIDDLE,1)
  RETURN
  END SUBROUTINE DBNSLV
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBSPVD(T,K,NDERIV,X,ILEFT,LDVNIK,VNIKX,WORK)
  !***BEGIN PROLOGUE  DBSPVD
  !***DATE WRITTEN   800901   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E3,K6
  !***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
  !         SPLINE
  !***AUTHOR  AMOS, D. E., (SNLA)
  !***PURPOSE  Calculates the value and all derivatives of order less than
  !        NDERIV of all basis functions which do not vanish at X.
  !***DESCRIPTION
  !
  ! Written by Carl de Boor and modified by D. E. Amos
  !
  ! Reference
  !     SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
  !
  ! Abstract    **** a double precision routine ****
  !
  !     DBSPVD is the BSPLVD routine of the reference.
  !
  !     DBSPVD calculates the value and all derivatives of order
  !     less than NDERIV of all basis functions which do not
  !     (possibly) vanish at X.  ILEFT is input such that
  !     T(ILEFT) .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X,
  !     ILO,ILEFT,MFLAG) will produce the proper ILEFT.  The output of
  !     DBSPVD is a matrix VNIKX(I,J) of dimension at least (K,NDERIV)
  !     whose columns contain the K nonzero basis functions and
  !     their NDERIV-1 right derivatives at X, I=1,K, J=1,NDERIV.
  !     These basis functions have indices ILEFT-K+I, I=1,K,
  !     K .LE. ILEFT .LE. N.  The nonzero part of the I-th basis
  !     function lies in (T(I),T(I+K)), I=1,N).
  !
  !     If X=T(ILEFT+1) then VNIKX contains left limiting values
  !     (left derivatives) at T(ILEFT+1).  In particular, ILEFT = N
  !     produces left limiting values at the right end point
  !     X=T(N+1).  To obtain left limiting values at T(I), I=K+1,N+1,
  !         set X= next lower distinct knot, call INTRV to get ILEFT,
  !     set X=T(I), and then call DBSPVD.
  !
  !     DBSPVD calls DBSPVN
  !
  ! Description of Arguments
  !     Input      T,X are double precision
  !      T       - knot vector of length N+K, where
  !                N = number of B-spline basis functions
  !                N = sum of knot multiplicities-K
  !      K       - order of the B-spline, K .GE. 1
  !      NDERIV  - number of derivatives = NDERIV-1,
  !                1 .LE. NDERIV .LE. K
  !      X       - argument of basis functions,
  !                T(K) .LE. X .LE. T(N+1)
  !      ILEFT   - largest integer such that
  !                T(ILEFT) .LE. X .LT.  T(ILEFT+1)
  !      LDVNIK  - leading dimension of matrix VNIKX
  !
  !     Output     VNIKX,WORK are double precision
  !      VNIKX   - matrix of dimension at least (K,NDERIV) contain-
  !                ing the nonzero basis functions at X and their
  !                derivatives columnwise.
  !      WORK    - a work vector of length (K+1)*(K+2)/2
  !
  ! Error Conditions
  !     Improper input is a fatal error
  !***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
  !             SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
  !             JUNE 1977, PP. 441-472.
  !***ROUTINES CALLED  DBSPVN,XERROR
  !***END PROLOGUE  DBSPVD
  !
  !
  INTEGER :: I,IDERIV,ILEFT,IPKMD,J,JJ,JLOW,JM,JP1MID,K,KMD, KP1, L, &
        LDUMMY, M, MHIGH, NDERIV
  DOUBLE PRECISION :: FACTOR, FKMD, T, V, VNIKX, WORK, X
  ! DIMENSION T(ILEFT+K), WORK((K+1)*(K+2)/2)
  ! A(I,J) = WORK(I+J*(J+1)/2),  I=1,J+1  J=1,K-1
  ! A(I,K) = W0RK(I+K*(K-1)/2)  I=1.K
  ! WORK(1) AND WORK((K+1)*(K+2)/2) ARE NOT USED.
  DIMENSION T(*), VNIKX(LDVNIK,NDERIV), WORK(*)
  !***FIRST EXECUTABLE STATEMENT  DBSPVD
  IF(K.LT.1) GO TO 200
  IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 205
  IF(LDVNIK.LT.K) GO TO 210
  IDERIV = NDERIV
  KP1 = K + 1
  JJ = KP1 - IDERIV
  CALL DBSPVN(T, JJ, K, 1, X, ILEFT, VNIKX, WORK, IWORK)
  IF (IDERIV.EQ.1) GO TO 100
  MHIGH = IDERIV
  DO M=2,MHIGH
    JP1MID = 1
    DO J=IDERIV,K
      VNIKX(J,IDERIV) = VNIKX(JP1MID,1)
      JP1MID = JP1MID + 1
    END DO
    IDERIV = IDERIV - 1
    JJ = KP1 - IDERIV
    CALL DBSPVN(T, JJ, K, 2, X, ILEFT, VNIKX, WORK, IWORK)
  END DO
  !
  JM = KP1*(KP1+1)/2
  DO L = 1,JM
    WORK(L) = 0.0D0
  END DO
  ! A(I,I) = WORK(I*(I+3)/2) = 1.0       I = 1,K
  L = 2
  J = 0
  DO I = 1,K
    J = J + L
    WORK(J) = 1.0D0
    L = L + 1
  END DO
  KMD = K
  DO M=2,MHIGH
    KMD = KMD - 1
    FKMD = FLOAT(KMD)
    I = ILEFT
    J = K
    JJ = J*(J+1)/2
    JM = JJ - J
    DO LDUMMY=1,KMD
      IPKMD = I + KMD
      FACTOR = FKMD/(T(IPKMD)-T(I))
      DO L=1,J
        WORK(L+JJ) = (WORK(L+JJ)-WORK(L+JM))*FACTOR
      END DO
      I = I - 1
      J = J - 1
      JJ = JM
      JM = JM - J
    END DO
  !
    DO I=1,K
      V = 0.0D0
      JLOW = MAX0(I,M)
      JJ = JLOW*(JLOW+1)/2
      DO J=JLOW,K
        V = WORK(I+JJ)*VNIKX(J,M) + V
        JJ = JJ + J + 1
      END DO
      VNIKX(I,M) = V
    END DO
  END DO
  100   RETURN
  !
  !
  200   CONTINUE
  CALL XERROR( ' DBSPVD,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
  RETURN
  205   CONTINUE
  CALL XERROR( ' DBSPVD,  NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K', &
        50, 2, 1)
  RETURN
  210   CONTINUE
  CALL XERROR( ' DBSPVD,  LDVNIK DOES NOT SATISFY LDVNIK.GE.K',45, &
        2, 1)
  RETURN
  END SUBROUTINE DBSPVD
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBSPVN(T,JHIGH,K,INDEX,X,ILEFT,VNIKX,WORK,IWORK)
  !***BEGIN PROLOGUE  DBSPVN
  !***DATE WRITTEN   800901   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E3,K6
  !***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
  !         SPLINE
  !***AUTHOR  AMOS, D. E., (SNLA)
  !***PURPOSE  Calculates the value of all (possibly) nonzero basis
  !        functions at X.
  !***DESCRIPTION
  !
  ! Written by Carl de Boor and modified by D. E. Amos
  !
  ! Reference
  !     SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
  !
  ! Abstract    **** a double precision routine ****
  !     DBSPVN is the BSPLVN routine of the reference.
  !
  !     DBSPVN calculates the value of all (possibly) nonzero basis
  !     functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K)
  !     .LE. X .LE. T(N+1) and J=IWORK is set inside the routine on
  !     the first call when INDEX=1.  ILEFT is such that T(ILEFT) .LE.
  !     X .LT. T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG)
  !     produces the proper ILEFT.  DBSPVN calculates using the basic
  !     algorithm needed in DBSPVD.  If only basis functions are
  !     desired, setting JHIGH=K and INDEX=1 can be faster than
  !     calling DBSPVD, but extra coding is required for derivatives
  !     (INDEX=2) and DBSPVD is set up for this purpose.
  !
  !     Left limiting values are set up as described in DBSPVD.
  !
  ! Description of Arguments
  !
  !     Input      T,X are double precision
  !      T       - knot vector of length N+K, where
  !                N = number of B-spline basis functions
  !                N = sum of knot multiplicities-K
  !      JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
  !      K       - highest possible order
  !      INDEX   - INDEX = 1 gives basis functions of order JHIGH
  !                      = 2 denotes previous entry with work, IWORK
  !                          values saved for subsequent calls to
  !                          DBSPVN.
  !      X       - argument of basis functions,
  !                T(K) .LE. X .LE. T(N+1)
  !      ILEFT   - largest integer such that
  !                T(ILEFT) .LE. X .LT.  T(ILEFT+1)
  !
  !     Output     VNIKX, WORK are double precision
  !      VNIKX   - vector of length K for spline values.
  !      WORK    - a work vector of length 2*K
  !      IWORK   - a work parameter.  Both WORK and IWORK contain
  !                information necessary to continue for INDEX = 2.
  !                When INDEX = 1 exclusively, these are scratch
  !                variables and can be used for other purposes.
  !
  ! Error Conditions
  !     Improper input is a fatal error.
  !***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
  !             SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
  !             JUNE 1977, PP. 441-472.
  !***ROUTINES CALLED  XERROR
  !***END PROLOGUE  DBSPVN
  !
  !
  INTEGER :: ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
  DOUBLE PRECISION :: T, VM, VMPREV, VNIKX, WORK, X
  ! DIMENSION T(ILEFT+JHIGH)
  DIMENSION T(*), VNIKX(K), WORK(*)
  ! CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
  ! WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
  !***FIRST EXECUTABLE STATEMENT  DBSPVN
  IF(K.LT.1) GO TO 90
  IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
  IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
  IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
  GO TO (10, 20), INDEX
   10   IWORK = 1
  VNIKX(1) = 1.0D0
  IF (IWORK.GE.JHIGH) GO TO 40
  !
   20   IPJ = ILEFT + IWORK
  WORK(IWORK) = T(IPJ) - X
  IMJP1 = ILEFT - IWORK + 1
  WORK(K+IWORK) = X - T(IMJP1)
  VMPREV = 0.0D0
  JP1 = IWORK + 1
  DO L=1,IWORK
    JP1ML = JP1 - L
    VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
    VNIKX(L) = VM*WORK(L) + VMPREV
    VMPREV = VM*WORK(K+JP1ML)
  END DO
  VNIKX(JP1) = VMPREV
  IWORK = JP1
  IF (IWORK.LT.JHIGH) GO TO 20
  !
   40   RETURN
  !
  !
   90   CONTINUE
  CALL XERROR( ' DBSPVN,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
  RETURN
  100   CONTINUE
  CALL XERROR( ' DBSPVN,  JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K', &
        48, 2, 1)
  RETURN
  105   CONTINUE
  CALL XERROR( ' DBSPVN,  INDEX IS NOT 1 OR 2',29,2,1)
  RETURN
  110   CONTINUE
  CALL XERROR( ' DBSPVN,  X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEF&
        &T+1)', 56, 2, 1)
  RETURN
  END SUBROUTINE DBSPVN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION DBVALU(T,A,N,K,IDERIV,X,INBV,WORK)
  !***BEGIN PROLOGUE  DBVALU
  !***DATE WRITTEN   800901   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E3,K6
  !***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
  !         SPLINE
  !***AUTHOR  AMOS, D. E., (SNLA)
  !***PURPOSE  Evaluates the B-representation of a B-spline at X for the
  !        function value or any of its derivatives.
  !***DESCRIPTION
  !
  ! Written by Carl de Boor and modified by D. E. Amos
  !
  ! Reference
  !     SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
  !
  ! Abstract   **** a double precision routine ****
  !     DBVALU is the BVALUE function of the reference.
  !
  !     DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
  !     at X for the function value on IDERIV=0 or any of its
  !     derivatives on IDERIV=1,2,...,K-1.  Right limiting values
  !     (right derivatives) are returned except at the right end
  !     point X=T(N+1) where left limiting values are computed.  The
  !     spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
  !     a fatal error message when X is outside of this interval.
  !
  !     To compute left derivatives or left limiting values at a
  !     knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
  !
  !     DBVALU calls DINTRV
  !
  ! Description of Arguments
  !
  !     Input      T,A,X are double precision
  !      T       - knot vector of length N+K
  !      A       - B-spline coefficient vector of length N
  !      N       - number of B-spline coefficients
  !                N = sum of knot multiplicities-K
  !      K       - order of the B-spline, K .GE. 1
  !      IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
  !                IDERIV = 0 returns the B-spline value
  !      X       - argument, T(K) .LE. X .LE. T(N+1)
  !      INBV    - an initialization parameter which must be set
  !                to 1 the first time DBVALU is called.
  !
  !     Output     WORK,DBVALU are double precision
  !      INBV    - INBV contains information for efficient process-
  !                ing after the initial call and INBV must not
  !                be changed by the user.  Distinct splines require
  !                distinct INBV parameters.
  !      WORK    - work vector of length 3*K.
  !      DBVALU  - value of the IDERIV-th derivative at X
  !
  ! Error Conditions
  !     An improper input is a fatal error
  !***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
  !             SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
  !             JUNE 1977, PP. 441-472.
  !***ROUTINES CALLED  DINTRV,XERROR
  !***END PROLOGUE  DBVALU
  !
  !
  INTEGER :: I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ, &
        IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
  DOUBLE PRECISION :: A, FKMJ, T, WORK, X
  DIMENSION T(*), A(N), WORK(*)
  !***FIRST EXECUTABLE STATEMENT  DBVALU
  DBVALU = 0.0D0
  IF(K.LT.1) GO TO 102
  IF(N.LT.K) GO TO 101
  IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
  KMIDER = K - IDERIV
  !
  ! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
  ! (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
  KM1 = K - 1
  CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
  IF (X.LT.T(K)) GO TO 120
  IF (MFLAG.EQ.0) GO TO 20
  IF (X.GT.T(I)) GO TO 130
   10   IF (I.EQ.K) GO TO 140
  I = I - 1
  IF (X.EQ.T(I)) GO TO 10
  !
  ! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
  ! WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
  !
   20   IMK = I - K
  DO J=1,K
    IMKPJ = IMK + J
    WORK(J) = A(IMKPJ)
  END DO
  IF (IDERIV.EQ.0) GO TO 60
  DO J=1,IDERIV
    KMJ = K - J
    FKMJ = DBLE(FLOAT(KMJ))
    DO JJ=1,KMJ
      IHI = I + JJ
      IHMKMJ = IHI - KMJ
      WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
    END DO
  END DO
  !
  ! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
  ! GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60   IF (IDERIV.EQ.KM1) GO TO 100
  IP1 = I + 1
  KPK = K + K
  J1 = K + 1
  J2 = KPK + 1
  DO J=1,KMIDER
    IPJ = I + J
    WORK(J1) = T(IPJ) - X
    IP1MJ = IP1 - J
    WORK(J2) = X - T(IP1MJ)
    J1 = J1 + 1
    J2 = J2 + 1
  END DO
  IDERP1 = IDERIV + 1
  DO J=IDERP1,KM1
    KMJ = K - J
    ILO = KMJ
    DO JJ=1,KMJ
      WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ) &
            *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
      ILO = ILO - 1
    END DO
  END DO
  100   DBVALU = WORK(1)
  RETURN
  !
  !
  101   CONTINUE
  CALL XERROR( ' DBVALU,  N DOES NOT SATISFY N.GE.K',35,2,1)
  RETURN
  102   CONTINUE
  CALL XERROR( ' DBVALU,  K DOES NOT SATISFY K.GE.1',35,2,1)
  RETURN
  110   CONTINUE
  CALL XERROR( ' DBVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', &
        50, 2, 1)
  RETURN
  120   CONTINUE
  CALL XERROR( ' DBVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)', &
        48, 2, 1)
  RETURN
  130   CONTINUE
  CALL XERROR( ' DBVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)', &
        47, 2, 1)
  RETURN
  140   CONTINUE
  CALL XERROR( ' DBVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED A&
        &T T(K)',    58, 2, 1)
  RETURN
  END FUNCTION DBVALU
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DINTRV(XT,LXT,X,ILO,ILEFT,MFLAG)
  !***BEGIN PROLOGUE  DINTRV
  !***DATE WRITTEN   800901   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  E3,K6
  !***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
          ! SPLINE
  !***AUTHOR  AMOS, D. E., (SNLA)
  !***PURPOSE  Computes the largest integer ILEFT in 1.LE.ILEFT.LE.LXT
         ! such that XT(ILEFT).LE.X where XT(*) is a subdivision of
         ! the X interval.
  !***DESCRIPTION
  !
  ! Written by Carl de Boor and modified by D. E. Amos
  !
  ! Reference
  !     SIAM J.  Numerical Analysis, 14, No. 3, June 1977, pp.441-472.
  !
  ! Abstract    **** a double precision routine ****
  !     DINTRV is the INTERV routine of the reference.
  !
  !     DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
  !     LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
  !     the X interval.  Precisely,
  !
  !                  X .LT. XT(1)                1         -1
  !     if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
  !       XT(LXT) .LE. X                         LXT        1,
  !
  !     That is, when multiplicities are present in the break point
  !     to the left of X, the largest index is taken for ILEFT.
  !
  ! Description of Arguments
  !
  !     Input      XT,X are double precision
  !      XT      - XT is a knot or break point vector of length LXT
  !      LXT     - length of the XT vector
  !      X       - argument
  !      ILO     - an initialization parameter which must be set
  !                to 1 the first time the spline array XT is
  !                processed by DINTRV.
  !
  !     Output
  !      ILO     - ILO contains information for efficient process-
  !                ing after the initial call and ILO must not be
  !                changed by the user.  Distinct splines require
  !                distinct ILO parameters.
  !      ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
  !      MFLAG   - signals when X lies out of bounds
  !
  ! Error Conditions
  !     None
  !***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
  !             SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
  !             JUNE 1977, PP. 441-472.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  DINTRV
  !
  !
  INTEGER :: IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
  DOUBLE PRECISION :: X, XT
  DIMENSION XT(LXT)
  !***FIRST EXECUTABLE STATEMENT  DINTRV
  IHI = ILO + 1
  IF (IHI.LT.LXT) GO TO 10
  IF (X.GE.XT(LXT)) GO TO 110
  IF (LXT.LE.1) GO TO 90
  ILO = LXT - 1
  IHI = LXT
  !
   10   IF (X.GE.XT(IHI)) GO TO 40
  IF (X.GE.XT(ILO)) GO TO 100
  !
  ! *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
  ISTEP = 1
   20   IHI = ILO
  ILO = IHI - ISTEP
  IF (ILO.LE.1) GO TO 30
  IF (X.GE.XT(ILO)) GO TO 70
  ISTEP = ISTEP*2
  GO TO 20
   30   ILO = 1
  IF (X.LT.XT(1)) GO TO 90
  GO TO 70
  ! *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40   ISTEP = 1
   50   ILO = IHI
  IHI = ILO + ISTEP
  IF (IHI.GE.LXT) GO TO 60
  IF (X.LT.XT(IHI)) GO TO 70
  ISTEP = ISTEP*2
  GO TO 50
   60   IF (X.GE.XT(LXT)) GO TO 110
  IHI = LXT
  !
  ! *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70   MIDDLE = (ILO+IHI)/2
  IF (MIDDLE.EQ.ILO) GO TO 100
  ! NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
  IF (X.LT.XT(MIDDLE)) GO TO 80
  ILO = MIDDLE
  GO TO 70
   80   IHI = MIDDLE
  GO TO 70
  ! *** SET OUTPUT AND RETURN
   90   MFLAG = -1
  ILEFT = 1
  RETURN
  100   MFLAG = 0
  ILEFT = ILO
  RETURN
  110   MFLAG = 1
  ILEFT = LXT
  RETURN
  END SUBROUTINE DINTRV
  !-----------------------------------------------------------------------------

END MODULE MOOSE_CMLIB_DBSPLIN
