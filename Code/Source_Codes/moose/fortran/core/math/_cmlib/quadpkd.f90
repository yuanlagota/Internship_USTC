!===============================================================================
!
!               Q U A D P K D   P A C K A G E
!
!
!  QUADPKD is a package of Fortran programs for computing definite
!  integrals of various forms.  These programs are the double
!  precision analogues of programs in the QUADPKS package. Please
!  refer to documentation of QUADPKS for an overview of its
!  capabilities.
!
!  The double precision subprogram names are the same as their
!  single precision counterparts except that they are preceded by
!  a "D".  Thus the single precision program QAG is DQAG in
!  double precision.
!
!
!***BEGIN PROLOGUE  QPDOC
!***DATE WRITTEN   810401   (YYMMDD)
!***REVISION DATE  840417   (YYMMDD)
!***CATEGORY NO.  H2
!***KEYWORDS  GUIDELINES FOR SELECTION,QUADPACK,SURVEY OF INTEGRATORS
!***AUTHOR  PIESSENS, ROBERT(APPL. MATH. AND PROGR. DIV.- K.U.LEUVEN)
!           DE DONKER, ELISE(APPL. MATH. AND PROGR. DIV.- K.U.LEUVEN
!           KAHANER,DAVID(NATIONAL BUREAU OF STANDARDS)
!***PURPOSE  QUADPACK documentation routine.
!***DESCRIPTION
!
! 1. Introduction
!    ------------
!    QUADPACK is a FORTRAN subroutine package for the numerical
!    computation of definite one-dimensional integrals. It originated
!    from a joint project of R. Piessens and E. de Doncker (Appl.
!    Math. and Progr. Div.- K.U.Leuven, Belgium), C. Ueberhuber (Inst.
!    Fuer Math.- Techn.U.Wien, Austria), and D. Kahaner (Nation. Bur.
!    of Standards- Washington D.C., U.S.A.).
!    Documentation routine QPDOC describes the package in the form it
!    was released from A.M.P.D.- Leuven, for adherence to the SLATEC
!    library in May 1981. Apart from a survey of the integrators, some
!    guidelines will be given in order to help the QUADPACK user with
!    selecting an appropriate routine or a combination of several
!    routines for handling his problem.
!
!    In the LONG DESCRIPTION of QPDOC it is demonstrated how to call
!    the integrators, by means of small example calling programs.
!
!    For precise guidelines involving the use of each routine in
!    particular, we refer to the extensive introductory comments
!    within each routine.
!
! 2. Survey
!    ------
!    The following list gives an overview of the QUADPACK integrators.
!    The routine names for the DOUBLE PRECISION versions are preceded
!    by the letter D.
!
!    - QNG  : Is a simple non-adaptive automatic integrator, based on
!             a sequence of rules with increasing degree of algebraic
!             precision (Patterson, 1968).
!
!    - QAG  : Is a simple globally adaptive integrator using the
!             strategy of Aind (Piessens, 1973). It is possible to
!             choose between 6 pairs of Gauss-Kronrod quadrature
!             formulae for the rule evaluation component. The pairs
!             of high degree of precision are suitable for handling
!             integration difficulties due to a strongly oscillating
!             integrand.
!
!    - QAGS : Is an integrator based on globally adaptive interval
!             subdivision in connection with extrapolation (de Doncker,
!             1978) by the Epsilon algorithm (Wynn, 1956).
!
!    - QAGP : Serves the same purposes as QAGS, but also allows
!             for eventual user-supplied information, i.e. the
!             abscissae of internal singularities, discontinuities
!             and other difficulties of the integrand function.
!             The algorithm is a modification of that in QAGS.
!
!    - QAGI : Handles integration over infinite intervals. The
!             infinite range is mapped onto a finite interval and
!             then the same strategy as in QAGS is applied.
!
!    - QAWO : Is a routine for the integration of COS(OMEGA*X)*F(X)
!             or SIN(OMEGA*X)*F(X) over a finite interval (A,B).
!             OMEGA is is specified by the user
!             The rule evaluation component is based on the
!             modified Clenshaw-Curtis technique.
!             An adaptive subdivision scheme is used connected with
!             an extrapolation procedure, which is a modification
!             of that in QAGS and provides the possibility to deal
!             even with singularities in F.
!
!    - QAWF : Calculates the Fourier cosine or Fourier sine
!             transform of F(X), for user-supplied interval (A,
!             INFINITY), OMEGA, and F. The procedure of QAWO is
!             used on successive finite intervals, and convergence
!             acceleration by means of the Epsilon algorithm (Wynn,
!             1956) is applied to the series of the integral
!             contributions.
!
!    - QAWS : Integrates W(X)*F(X) over (A,B) with A.LT.B finite,
!             and   W(X) = ((X-A)**ALFA)*((B-X)**BETA)*V(X)
!             where V(X) = 1 or LOG(X-A) or LOG(B-X)
!                            or LOG(X-A)*LOG(B-X)
!             and   ALFA.GT.(-1), BETA.GT.(-1).
!             The user specifies A, B, ALFA, BETA and the type of
!             the function V.
!             A globally adaptive subdivision strategy is applied,
!             with modified Clenshaw-Curtis integration on the
!             subintervals which contain A or B.
!
!    - QAWC : Computes the Cauchy Principal Value of F(X)/(X-C)
!             over a finite interval (A,B) and for
!             user-determined C.
!             The strategy is globally adaptive, and modified
!             Clenshaw-Curtis integration is used on the subranges
!             which contain the point X = C.
!
!  Each of the routines above also has a "more detailed" version
!    with a name ending in E, as QAGE.  These provide more
!    information and control than the easier versions.
!
!
!   The preceeding routines are all automatic.  That is, the user
!      inputs his problem and an error tolerance.  The routine
!      attempts to perform the integration to within the requested
!      absolute or relative error.
!   There are, in addition, a number of non-automatic integrators.
!      These are most useful when the problem is such that the
!      user knows that a fixed rule will provide the accuracy
!      required.  Typically they return an error estimate but make
!      no attempt to satisfy any particular input error request.
!
!      QK15
!      QK21
!      QK31
!      QK41
!      QK51
!      QK61
!           Estimate the integral on [a,b] using 15, 21,..., 61
!           point rule and return an error estimate.
!      QK15I 15 point rule for (semi)infinite interval.
!      QK15W 15 point rule for special singular weight functions.
!      QC25C 25 point rule for Cauchy Principal Values
!      QC25F 25 point rule for sin/cos integrand.
!      QMOMO Integrates k-th degree Chebychev polynomial times
!            function with various explicit singularities.
!
! 3. Guidelines for the use of QUADPACK
!    ----------------------------------
!    Here it is not our purpose to investigate the question when
!    automatic quadrature should be used. We shall rather attempt
!    to help the user who already made the decision to use QUADPACK,
!    with selecting an appropriate routine or a combination of
!    several routines for handling his problem.
!
!    For both quadrature over finite and over infinite intervals,
!    one of the first questions to be answered by the user is
!    related to the amount of computer time he wants to spend,
!    versus his -own- time which would be needed, for example, for
!    manual subdivision of the interval or other analytic
!    manipulations.
!
!    (1) The user may not care about computer time, or not be
!        willing to do any analysis of the problem. especially when
!        only one or a few integrals must be calculated, this attitude
!        can be perfectly reasonable. In this case it is clear that
!        either the most sophisticated of the routines for finite
!        intervals, QAGS, must be used, or its analogue for infinite
!        intervals, GAGI. These routines are able to cope with
!        rather difficult, even with improper integrals.
!        This way of proceeding may be expensive. But the integrator
!        is supposed to give you an answer in return, with additional
!        information in the case of a failure, through its error
!        estimate and flag. Yet it must be stressed that the programs
!        cannot be totally reliable.
!        ------
!
!    (2) The user may want to examine the integrand function.
!        If bad local difficulties occur, such as a discontinuity, a
!        singularity, derivative singularity or high peak at one or
!        more points within the interval, the first advice is to
!        split up the interval at these points. The integrand must
!        then be examinated over each of the subintervals separately,
!        so that a suitable integrator can be selected for each of
!        them. If this yields problems involving relative accuracies
!        to be imposed on -finite- subintervals, one can make use of
!        QAGP, which must be provided with the positions of the local
!        difficulties. However, if strong singularities are present
!        and a high accuracy is requested, application of QAGS on the
!        subintervals may yield a better result.
!
!        For quadrature over finite intervals we thus dispose of QAGS
!        and
!        - QNG for well-behaved integrands,
!        - QAG for functions with an oscillating behaviour of a non
!          specific type,
!        - QAWO for functions, eventually singular, containing a
!          factor COS(OMEGA*X) or SIN(OMEGA*X) where OMEGA is known,
!        - QAWS for integrands with Algebraico-Logarithmic end point
!          singularities of known type,
!        - QAWC for Cauchy Principal Values.
!
!        Remark
!        ------
!        On return, the work arrays in the argument lists of the
!        adaptive integrators contain information about the interval
!        subdivision process and hence about the integrand behaviour:
!        the end points of the subintervals, the local integral
!        contributions and error estimates, and eventually other
!        characteristics. For this reason, and because of its simple
!        globally adaptive nature, the routine QAG in particular is
!        well-suited for integrand examination. Difficult spots can
!        be located by investigating the error estimates on the
!        subintervals.
!
!        For infinite intervals we provide only one general-purpose
!        routine, QAGI. It is based on the QAGS algorithm applied
!        after a transformation of the original interval into (0,1).
!        Yet it may eventuate that another type of transformation is
!        more appropriate, or one might prefer to break up the
!        original interval and use QAGI only on the infinite part
!        and so on. These kinds of actions suggest a combined use of
!        different QUADPACK integrators. Note that, when the only
!        difficulty is an integrand singularity at the finite
!        integration limit, it will in general not be necessary to
!        break up the interval, as QAGI deals with several types of
!        singularity at the boundary point of the integration range.
!        It also handles slowly convergent improper integrals, on
!        the condition that the integrand does not oscillate over
!        the entire infinite interval. If it does we would advise
!        to sum succeeding positive and negative contributions to
!        the integral -e.g. integrate between the zeros- with one
!        or more of the finite-range integrators, and apply
!        convergence acceleration eventually by means of QUADPACK
!        subroutine QELG which implements the Epsilon algorithm.
!        Such quadrature problems include the Fourier transform as
!        a special case. Yet for the latter we have an automatic
!        integrator available, QAWF.
!***LONG DESCRIPTION
!
! 4. Example Programs
!    ----------------
! 4.1. Calling Program for QNG
!      -----------------------
!
!            REAL A,ABSERR,B,F,EPSABS,EPSREL,RESULT
!            INTEGER IER,NEVAL
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            CALL QNG(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = EXP(X)/(X*X+0.1E+01)
!            RETURN
!            END
!
! 4.2. Calling Program for QAG
!      -----------------------
!
!            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
!            INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
!            DIMENSION IWORK(100),WORK(400)
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            KEY = 6
!            LIMIT = 100
!            LENW = LIMIT*4
!            CALL QAG(F,A,B,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,
!           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 2.0E0/(2.0E0+SIN(31.41592653589793E0*X))
!            RETURN
!            END
!
! 4.3. Calling Program for QAGS
!      ------------------------
!
!            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
!            INTEGER IER,IWORK,LAST,LENW,LIMIT,NEVAL
!            DIMENSION IWORK(100),WORK(400)
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            LIMIT = 100
!            LENW = LIMIT*4
!            CALL QAGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER,
!           *  LIMIT,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 0.0E0
!            IF(X.GT.0.0E0) F = 1.0E0/SQRT(X)
!            RETURN
!            END
!
! 4.4. Calling Program for QAGP
!      ------------------------
!
!            REAL A,ABSERR,B,EPSABS,EPSREL,F,POINTS,RESULT,WORK
!            INTEGER IER,IWORK,LAST,LENIW,LENW,LIMIT,NEVAL,NPTS2
!            DIMENSION IWORK(204),POINTS(4),WORK(404)
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            NPTS2 = 4
!            POINTS(1) = 1.0E0/7.0E0
!            POINTS(2) = 2.0E0/3.0E0
!            LIMIT = 100
!            LENIW = LIMIT*2+NPTS2
!            LENW = LIMIT*4+NPTS2
!            CALL QAGP(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,RESULT,ABSERR,
!           *  NEVAL,IER,LENIW,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 0.0E+00
!            IF(X.NE.1.0E0/7.0E0.AND.X.NE.2.0E0/3.0E0) F =
!           *  ABS(X-1.0E0/7.0E0)**(-0.25E0)*
!           *  ABS(X-2.0E0/3.0E0)**(-0.55E0)
!            RETURN
!            END
!
! 4.5. Calling Program for QAGI
!      ------------------------
!
!            REAL ABSERR,BOUN,EPSABS,EPSREL,F,RESULT,WORK
!            INTEGER IER,INF,IWORK,LAST,LENW,LIMIT,NEVAL
!            DIMENSION IWORK(100),WORK(400)
!            EXTERNAL F
!            BOUN = 0.0E0
!            INF = 1
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            LIMIT = 100
!            LENW = LIMIT*4
!            CALL QAGI(F,BOUN,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
!           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 0.0E0
!            IF(X.GT.0.0E0) F = SQRT(X)*ALOG(X)/
!           *                   ((X+1.0E0)*(X+2.0E0))
!            RETURN
!            END
!
! 4.6. Calling Program for QAWO
!      ------------------------
!
!            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,OMEGA,WORK
!            INTEGER IER,INTEGR,IWORK,LAST,LENIW,LENW,LIMIT,MAXP1,NEVAL
!            DIMENSION IWORK(200),WORK(925)
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            OMEGA = 10.0E0
!            INTEGR = 1
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            LIMIT = 100
!            LENIW = LIMIT*2
!            MAXP1 = 21
!            LENW = LIMIT*4+MAXP1*25
!            CALL QAWO(F,A,B,OMEGA,INTEGR,EPSABS,EPSREL,RESULT,ABSERR,
!           *  NEVAL,IER,LENIW,MAXP1,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 0.0E0
!            IF(X.GT.0.0E0) F = EXP(-X)*ALOG(X)
!            RETURN
!            END
!
! 4.7. Calling Program for QAWF
!      ------------------------
!
!            REAL A,ABSERR,EPSABS,F,RESULT,OMEGA,WORK
!            INTEGER IER,INTEGR,IWORK,LAST,LENIW,LENW,LIMIT,LIMLST,
!           *  LST,MAXP1,NEVAL
!            DIMENSION IWORK(250),WORK(1025)
!            EXTERNAL F
!            A = 0.0E0
!            OMEGA = 8.0E0
!            INTEGR = 2
!            EPSABS = 1.0E-3
!            LIMLST = 50
!            LIMIT = 100
!            LENIW = LIMIT*2+LIMLST
!            MAXP1 = 21
!            LENW = LENIW*2+MAXP1*25
!            CALL QAWF(F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
!           *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            IF(X.GT.0.0E0) F = SIN(50.0E0*X)/(X*SQRT(X))
!            RETURN
!            END
!
! 4.8. Calling Program for QAWS
!      ------------------------
!
!            REAL A,ABSERR,ALFA,B,BETA,EPSABS,EPSREL,F,RESULT,WORK
!            INTEGER IER,INTEGR,IWORK,LAST,LENW,LIMIT,NEVAL
!            DIMENSION IWORK(100),WORK(400)
!            EXTERNAL F
!            A = 0.0E0
!            B = 1.0E0
!            ALFA = -0.5E0
!            BETA = -0.5E0
!            INTEGR = 1
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            LIMIT = 100
!            LENW = LIMIT*4
!            CALL QAWS(F,A,B,ALFA,BETA,INTEGR,EPSABS,EPSREL,RESULT,
!           *  ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = SIN(10.0E0*X)
!            RETURN
!            END
!
! 4.9. Calling Program for QAWC
!      ------------------------
!
!            REAL A,ABSERR,B,C,EPSABS,EPSREL,F,RESULT,WORK
!            INTEGER IER,IWORK,LAST,LENW,LIMIT,NEVAL
!            DIMENSION IWORK(100),WORK(400)
!            EXTERNAL F
!            A = -1.0E0
!            B = 1.0E0
!            C = 0.5E0
!            EPSABS = 0.0E0
!            EPSREL = 1.0E-3
!            LIMIT = 100
!            LENW = LIMIT*4
!            CALL QAWC(F,A,B,C,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
!           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
!      C  INCLUDE WRITE STATEMENTS
!            STOP
!            END
!      C
!            REAL FUNCTION F(X)
!            REAL X
!            F = 1.0E0/(X*X+1.0E-4)
!            RETURN
!            END
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  QPDOC
!
!===============================================================================
MODULE MOOSE_CMLIB_QUADPKD
  USE MOOSE_CMLIB_MACHCON
  USE MOOSE_CMLIB_XERROR

  CONTAINS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQAGE(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR, &
        NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
  !***BEGIN PROLOGUE  DQAGE
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DABS to generic MAX/ABS
  !***CATEGORY NO.  H2A1A1
  !***KEYWORDS  AUTOMATIC INTEGRATOR,GAUSS-KRONROD,GENERAL-PURPOSE,
          ! GLOBALLY ADAPTIVE,INTEGRAND EXAMINATOR
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  The routine calculates an approximation result to a given
        !  definite integral   I = Integral of F over (A,B),
        !  hopefully satisfying following claim for accuracy
        !  ABS(I-RESLT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***DESCRIPTION
  !
  !    Computation of a definite integral
  !    Standard fortran subroutine
  !    Double precision version
  !
  !    PARAMETERS
  !     ON ENTRY
  !        F      - Double precision
  !                 Function subprogram defining the integrand
  !                 function F(X). The actual name for F needs to be
  !                 declared E X T E R N A L in the driver program.
  !
  !        A      - Double precision
  !                 Lower limit of integration
  !
  !        B      - Double precision
  !                 Upper limit of integration
  !
  !        EPSABS - Double precision
  !                 Absolute accuracy requested
  !        EPSREL - Double precision
  !                 Relative accuracy requested
  !                 If  EPSABS.LE.0
  !                 and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                 the routine will end with IER = 6.
  !
  !        KEY    - Integer
  !                 Key for choice of local integration rule
  !                 A Gauss-Kronrod pair is used with
  !                      7 - 15 points if KEY.LT.2,
  !                     10 - 21 points if KEY = 2,
  !                     15 - 31 points if KEY = 3,
  !                     20 - 41 points if KEY = 4,
  !                     25 - 51 points if KEY = 5,
  !                     30 - 61 points if KEY.GT.5.
  !
  !        LIMIT  - Integer
  !                 Gives an upperbound on the number of subintervals
  !                 in the partition of (A,B), LIMIT.GE.1.
  !
  !     ON RETURN
  !        RESULT - Double precision
  !                 Approximation to the integral
  !
  !        ABSERR - Double precision
  !                 Estimate of the modulus of the absolute error,
  !                 which should equal or exceed ABS(I-RESULT)
  !
  !        NEVAL  - Integer
  !                 Number of integrand evaluations
  !
  !        IER    - Integer
  !                 IER = 0 Normal and reliable termination of the
  !                         routine. It is assumed that the requested
  !                         accuracy has been achieved.
  !                 IER.GT.0 Abnormal termination of the routine
  !                         The estimates for result and error are
  !                         less reliable. It is assumed that the
  !                         requested accuracy has not been achieved.
  !        ERROR MESSAGES
  !                 IER = 1 Maximum number of subdivisions allowed
  !                         has been achieved. One can allow more
  !                         subdivisions by increasing the value
  !                         of LIMIT.
  !                         However, if this yields no improvement it
  !                         is rather advised to analyze the integrand
  !                         in order to determine the integration
  !                         difficulties. If the position of a local
  !                         difficulty can be determined(e.g.
  !                         SINGULARITY, DISCONTINUITY within the
  !                         interval) one will probably gain from
  !                         splitting up the interval at this point
  !                         and calling the integrator on the
  !                         subranges. If possible, an appropriate
  !                         special-purpose integrator should be used
  !                         which is designed for handling the type of
  !                         difficulty involved.
  !                     = 2 The occurrence of roundoff error is
  !                         detected, which prevents the requested
  !                         tolerance from being achieved.
  !                     = 3 Extremely bad integrand behaviour occurs
  !                         at some points of the integration
  !                         interval.
  !                     = 6 The input is invalid, because
  !                         (EPSABS.LE.0 and
  !                          EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                         RESULT, ABSERR, NEVAL, LAST, RLIST(1) ,
  !                         ELIST(1) and IORD(1) are set to zero.
  !                         ALIST(1) and BLIST(1) are set to A and B
  !                         respectively.
  !
  !        ALIST   - Double precision
  !                  Vector of dimension at least LIMIT, the first
  !                   LAST  elements of which are the left
  !                  end points of the subintervals in the partition
  !                  of the given integration range (A,B)
  !
  !        BLIST   - Double precision
  !                  Vector of dimension at least LIMIT, the first
  !                   LAST  elements of which are the right
  !                  end points of the subintervals in the partition
  !                  of the given integration range (A,B)
  !
  !        RLIST   - Double precision
  !                  Vector of dimension at least LIMIT, the first
  !                   LAST  elements of which are the
  !                  integral approximations on the subintervals
  !
  !        ELIST   - Double precision
  !                  Vector of dimension at least LIMIT, the first
  !                   LAST  elements of which are the moduli of the
  !                  absolute error estimates on the subintervals
  !
  !        IORD    - Integer
  !                  Vector of dimension at least LIMIT, the first K
  !                  elements of which are pointers to the
  !                  error estimates over the subintervals,
  !                  such that ELIST(IORD(1)), ...,
  !                  ELIST(IORD(K)) form a decreasing sequence,
  !                  with K = LAST if LAST.LE.(LIMIT/2+2), and
  !                  K = LIMIT+1-LAST otherwise
  !
  !        LAST    - Integer
  !                  Number of subintervals actually produced in the
  !                  subdivision process
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH,DQK15,DQK21,DQK31,DQK41,DQK51,DQK61,DQPSRT
  !***END PROLOGUE  DQAGE
  !
  DOUBLE PRECISION :: A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B, &
        BLIST,B1,B2,DEFABS,DEFAB1,DEFAB2,ELIST,EPMACH, &
        EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F, &
        RESABS,RESULT,RLIST,UFLOW
  INTEGER :: IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL, &
        NRMAX
  !
  DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT), &
        RLIST(LIMIT)
  !
  EXTERNAL F
  !
  !        LIST OF MAJOR VARIABLES
  !        -----------------------
  !
  !       ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
  !                   CONSIDERED UP TO NOW
  !       BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
  !                   CONSIDERED UP TO NOW
  !       RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
  !                  (ALIST(I),BLIST(I))
  !       ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
  !       MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
  !                   ERROR ESTIMATE
  !       ERRMAX    - ELIST(MAXERR)
  !       AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
  !       ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
  !       ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
  !                   ABS(RESULT))
  !       *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
  !       *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
  !       LAST      - INDEX FOR SUBDIVISION
  !
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH  IS THE LARGEST RELATIVE SPACING.
  !       UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQAGE
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  !       TEST ON VALIDITY OF PARAMETERS
  !       ------------------------------
  !
  IER = 0
  NEVAL = 0
  LAST = 0
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  ALIST(1) = A
  BLIST(1) = B
  RLIST(1) = 0.0D+00
  ELIST(1) = 0.0D+00
  IORD(1) = 0
  IF(EPSABS.LE.0.0D+00.AND. &
        EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)) IER = 6
  IF(IER.EQ.6) GO TO 999
  !
  !       FIRST APPROXIMATION TO THE INTEGRAL
  !       -----------------------------------
  !
  KEYF = KEY
  IF(KEY.LE.0) KEYF = 1
  IF(KEY.GE.7) KEYF = 6
  NEVAL = 0
  IF(KEYF.EQ.1) CALL DQK15(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  IF(KEYF.EQ.2) CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  IF(KEYF.EQ.3) CALL DQK31(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  IF(KEYF.EQ.4) CALL DQK41(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  IF(KEYF.EQ.5) CALL DQK51(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  IF(KEYF.EQ.6) CALL DQK61(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  LAST = 1
  RLIST(1) = RESULT
  ELIST(1) = ABSERR
  IORD(1) = 1
  !
  !       TEST ON ACCURACY.
  !
  ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
  IF(ABSERR.LE.0.5D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
  IF(LIMIT.EQ.1) IER = 1
  IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS) &
        .OR.ABSERR.EQ.0.0D+00) GO TO 60
  !
  !       INITIALIZATION
  !       --------------
  !
  !
  ERRMAX = ABSERR
  MAXERR = 1
  AREA = RESULT
  ERRSUM = ABSERR
  NRMAX = 1
  IROFF1 = 0
  IROFF2 = 0
  !
  !       MAIN DO-LOOP
  !       ------------
  !
  DO LAST = 2,LIMIT
  !
  !       BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
  !
    A1 = ALIST(MAXERR)
    B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
    A2 = B1
    B2 = BLIST(MAXERR)
    IF(KEYF.EQ.1) CALL DQK15(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.2) CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.3) CALL DQK31(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.4) CALL DQK41(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.5) CALL DQK51(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.6) CALL DQK61(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    IF(KEYF.EQ.1) CALL DQK15(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    IF(KEYF.EQ.2) CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    IF(KEYF.EQ.3) CALL DQK31(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    IF(KEYF.EQ.4) CALL DQK41(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    IF(KEYF.EQ.5) CALL DQK51(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    IF(KEYF.EQ.6) CALL DQK61(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
  !
  !       IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
  !       AND ERROR AND TEST FOR ACCURACY.
  !
    NEVAL = NEVAL+1
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
    IF(ABS(RLIST(MAXERR)-AREA12).LE.0.1D-04*ABS(AREA12) &
          .AND.ERRO12.GE.0.99D+00*ERRMAX) IROFF1 = IROFF1+1
    IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
    IF(ERRSUM.LE.ERRBND) GO TO 8
  !
  !       TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
  !
    IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
  !
  !       SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
  !       EQUALS LIMIT.
  !
    IF(LAST.EQ.LIMIT) IER = 1
  !
  !       SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
  !       AT A POINT OF THE INTEGRATION RANGE.
  !
    IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03* &
          EPMACH)*(ABS(A2)+0.1D+04*UFLOW)) IER = 3
  !
  !       APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
  !
    8   IF(ERROR2.GT.ERROR1) GO TO 10
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    GO TO 20
   10   ALIST(MAXERR) = A2
    ALIST(LAST) = A1
    BLIST(LAST) = B1
    RLIST(MAXERR) = AREA2
    RLIST(LAST) = AREA1
    ELIST(MAXERR) = ERROR2
    ELIST(LAST) = ERROR1
  !
  !       CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
  !       IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
  !       WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
  !
   20   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
  ! ***JUMP OUT OF DO-LOOP
    IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
  END DO
  !
  !       COMPUTE FINAL RESULT.
  !       ---------------------
  !
   40   RESULT = 0.0D+00
  DO K=1,LAST
    RESULT = RESULT+RLIST(K)
  END DO
  ABSERR = ERRSUM
   60   IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
  IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
  999   RETURN
  END SUBROUTINE DQAGE
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQAG(F,A,B,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER, &
        LIMIT,LENW,LAST,IWORK,WORK)
  !***BEGIN PROLOGUE  DQAG
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***CATEGORY NO.  H2A1A1
  !***KEYWORDS  AUTOMATIC INTEGRATOR,GAUSS-KRONROD,GENERAL-PURPOSE,
          ! GLOBALLY ADAPTIVE,INTEGRAND EXAMINATOR
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  The routine calculates an approximation result to a given
        !  definite integral I = integral of F over (A,B),
        !  hopefully satisfying following claim for accuracy
        !  ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***DESCRIPTION
  !
  !    Computation of a definite integral
  !    Standard fortran subroutine
  !    Double precision version
  !
  !        F      - Double precision
  !                 Function subprogam defining the integrand
  !                 Function F(X). The actual name for F needs to be
  !                 Declared E X T E R N A L in the driver program.
  !
  !        A      - Double precision
  !                 Lower limit of integration
  !
  !        B      - Double precision
  !                 Upper limit of integration
  !
  !        EPSABS - Double precision
  !                 Absolute accoracy requested
  !        EPSREL - Double precision
  !                 Relative accuracy requested
  !                 If  EPSABS.LE.0
  !                 And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                 The routine will end with IER = 6.
  !
  !        KEY    - Integer
  !                 Key for choice of local integration rule
  !                 A GAUSS-KRONROD PAIR is used with
  !                   7 - 15 POINTS If KEY.LT.2,
  !                  10 - 21 POINTS If KEY = 2,
  !                  15 - 31 POINTS If KEY = 3,
  !                  20 - 41 POINTS If KEY = 4,
  !                  25 - 51 POINTS If KEY = 5,
  !                  30 - 61 POINTS If KEY.GT.5.
  !
  !     ON RETURN
  !        RESULT - Double precision
  !                 Approximation to the integral
  !
  !        ABSERR - Double precision
  !                 Estimate of the modulus of the absolute error,
  !                 Which should EQUAL or EXCEED ABS(I-RESULT)
  !
  !        NEVAL  - Integer
  !                 Number of integrand evaluations
  !
  !        IER    - Integer
  !                 IER = 0 Normal and reliable termination of the
  !                         routine. It is assumed that the requested
  !                         accuracy has been achieved.
  !                 IER.GT.0 Abnormal termination of the routine
  !                         The estimates for RESULT and ERROR are
  !                         Less reliable. It is assumed that the
  !                         requested accuracy has not been achieved.
  !                  ERROR MESSAGES
  !                 IER = 1 Maximum number of subdivisions allowed
  !                         has been achieved. One can allow more
  !                         subdivisions by increasing the value of
  !                         LIMIT (and taking the according dimension
  !                         adjustments into account). HOWEVER, If
  !                         this yield no improvement it is advised
  !                         to analyze the integrand in order to
  !                         determine the integration difficulaties.
  !                         If the position of a local difficulty can
  !                         be determined (I.E.SINGULARITY,
  !                         DISCONTINUITY WITHIN THE INTERVAL) One
  !                         will probably gain from splitting up the
  !                         interval at this point and calling the
  !                         INTEGRATOR on the SUBRANGES. If possible,
  !                         AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
  !                         should be used which is designed for
  !                         handling the type of difficulty involved.
  !                     = 2 The occurrence of roundoff error is
  !                         detected, which prevents the requested
  !                         tolerance from being achieved.
  !                     = 3 Extremely bad integrand behaviour occurs
  !                         at some points of the integration
  !                         interval.
  !                     = 6 The input is invalid, because
  !                         (EPSABS.LE.0 AND
  !                          EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                         OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
  !                         RESULT, ABSERR, NEVAL, LAST are set
  !                         to zero.
  !                         EXCEPT when LENW is invalid, IWORK(1),
  !                         WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
  !                         set to zero, WORK(1) is set to A and
  !                         WORK(LIMIT+1) to B.
  !
  !     DIMENSIONING PARAMETERS
  !        LIMIT - Integer
  !                Dimensioning parameter for IWORK
  !                Limit determines the maximum number of subintervals
  !                in the partition of the given integration interval
  !                (A,B), LIMIT.GE.1.
  !                If LIMIT.LT.1, the routine will end with IER = 6.
  !
  !        LENW  - Integer
  !                Dimensioning parameter for work
  !                LENW must be at least LIMIT*4.
  !                IF LENW.LT.LIMIT*4, the routine will end with
  !                IER = 6.
  !
  !        LAST  - Integer
  !                On return, LAST equals the number of subintervals
  !                produced in the subdiviosion process, which
  !                determines the number of significant elements
  !                actually in the WORK ARRAYS.
  !
  !     WORK ARRAYS
  !        IWORK - Integer
  !                Vector of dimension at least limit, the first K
  !                elements of which contain pointers to the error
  !                estimates over the subintervals, such that
  !                WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K))
  !                form a decreasing sequence with K = LAST If
  !                LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
  !
  !        WORK  - Double precision
  !                Vector of dimension at least LENW
  !                on return
  !                WORK(1), ..., WORK(LAST) contain the left end
  !                points of the subintervals in the partition of
  !                 (A,B),
  !                WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
  !                 right end points,
  !                WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                 the integral approximations over the subintervals,
  !                WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
  !                 the error estimates.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DQAGE,XERROR
  !***END PROLOGUE  DQAG
  DOUBLE PRECISION :: A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
  INTEGER :: IER,IWORK,KEY,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
  !
  DIMENSION IWORK(LIMIT),WORK(LENW)
  !
  EXTERNAL F
  !
  !     CHECK VALIDITY OF LENW.
  !
  !***FIRST EXECUTABLE STATEMENT  DQAG
  IER = 6
  NEVAL = 0
  LAST = 0
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
  !
  !     PREPARE CALL FOR DQAGE.
  !
  L1 = LIMIT+1
  L2 = LIMIT+L1
  L3 = LIMIT+L2
  !
  CALL DQAGE(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR,NEVAL, &
        IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
  !
  !     CALL ERROR HANDLER IF NECESSARY.
  !
  LVL = 0
10   IF(IER.EQ.6) LVL = 1
  IF(IER.NE.0) CALL XERROR( 'ABNORMAL RETURN FROM DQAG ',26,IER,LVL)
  RETURN
  END SUBROUTINE DQAG
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK15(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK15
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  15-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B), with error
        !                 estimate
        !             J = integral of ABS(F) over (A,B)
  !***DESCRIPTION
  !
  !       Integration rules
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS
  !        ON ENTRY
  !          F      - Double precision
  !                   Function subprogram defining the integrand
  !                   FUNCTION F(X). The actual name for F needs to be
  !                   Declared E X T E R N A L in the calling program.
  !
  !          A      - Double precision
  !                   Lower limit of integration
  !
  !          B      - Double precision
  !                   Upper limit of integration
  !
  !        ON RETURN
  !          RESULT - Double precision
  !                   Approximation to the integral I
  !                   Result is computed by applying the 15-POINT
  !                   KRONROD RULE (RESK) obtained by optimal addition
  !                   of abscissae to the7-POINT GAUSS RULE(RESG).
  !
  !          ABSERR - Double precision
  !                   Estimate of the modulus of the absolute error,
  !                   which should not exceed ABS(I-RESULT)
  !
  !          RESABS - Double precision
  !                   Approximation to the integral J
  !
  !          RESASC - Double precision
  !                   Approximation to the integral of ABS(F-I/(B-A))
  !                   over (A,B)
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK15
  !
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !       BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !       CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
  !                XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
  !                GAUSS RULE
  !                XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                ADDED TO THE 7-POINT GAUSS RULE
  !
  !       WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
  !
  !       WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.129484966168869693270611432679082D0 /
  DATA WG  (  2) / 0.279705391489276667901467771423780D0 /
  DATA WG  (  3) / 0.381830050505118944950369775488975D0 /
  DATA WG  (  4) / 0.417959183673469387755102040816327D0 /
  !
  DATA XGK (  1) / 0.991455371120812639206854697526329D0 /
  DATA XGK (  2) / 0.949107912342758524526189684047851D0 /
  DATA XGK (  3) / 0.864864423359769072789712788640926D0 /
  DATA XGK (  4) / 0.741531185599394439863864773280788D0 /
  DATA XGK (  5) / 0.586087235467691130294144838258730D0 /
  DATA XGK (  6) / 0.405845151377397166906606412076961D0 /
  DATA XGK (  7) / 0.207784955007898467600689403773245D0 /
  DATA XGK (  8) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.022935322010529224963732008058970D0 /
  DATA WGK (  2) / 0.063092092629978553290700663189204D0 /
  DATA WGK (  3) / 0.104790010322250183839876322541518D0 /
  DATA WGK (  4) / 0.140653259715525918745189590510238D0 /
  DATA WGK (  5) / 0.169004726639267902826583426598550D0 /
  DATA WGK (  6) / 0.190350578064785409913256402421014D0 /
  DATA WGK (  7) / 0.204432940075298892414161999234649D0 /
  DATA WGK (  8) / 0.209482141084727828012999174891714D0 /
  !
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       ABSC   - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
  !       RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
  !       RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                I.E. TO I/(B-A)
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK15
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
  !       THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  FC = F(CENTR)
  RESG = FC*WG(4)
  RESK = FC*WGK(8)
  RESABS = ABS(RESK)
  DO J=1,3
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J = 1,4
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(8)*ABS(FC-RESKH)
  DO J=1,7
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK15
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK21(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK21
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  21-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B), with error
        !                 estimate
        !             J = Integral of ABS(F) over (A,B)
  !***DESCRIPTION
  !
  !       Integration rules
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS
  !        ON ENTRY
  !          F      - Double precision
  !                   Function subprogram defining the integrand
  !                   FUNCTION F(X). The actual name for F needs to be
  !                   Declared E X T E R N A L in the driver program.
  !
  !          A      - Double precision
  !                   Lower limit of integration
  !
  !          B      - Double precision
  !                   Upper limit of integration
  !
  !        ON RETURN
  !          RESULT - Double precision
  !                   Approximation to the integral I
  !                   RESULT is computed by applying the 21-POINT
  !                   KRONROD RULE (RESK) obtained by optimal addition
  !                   of abscissae to the 10-POINT GAUSS RULE (RESG).
  !
  !          ABSERR - Double precision
  !                   Estimate of the modulus of the absolute error,
  !                   which should not exceed ABS(I-RESULT)
  !
  !          RESABS - Double precision
  !                   Approximation to the integral J
  !
  !          RESASC - Double precision
  !                   Approximation to the integral of ABS(F-I/(B-A))
  !                   over (A,B)
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK21
  !
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !       BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !       CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
  !                XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
  !                GAUSS RULE
  !                XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                ADDED TO THE 10-POINT GAUSS RULE
  !
  !       WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
  !
  !       WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.066671344308688137593568809893332D0 /
  DATA WG  (  2) / 0.149451349150580593145776339657697D0 /
  DATA WG  (  3) / 0.219086362515982043995534934228163D0 /
  DATA WG  (  4) / 0.269266719309996355091226921569469D0 /
  DATA WG  (  5) / 0.295524224714752870173892994651338D0 /
  !
  DATA XGK (  1) / 0.995657163025808080735527280689003D0 /
  DATA XGK (  2) / 0.973906528517171720077964012084452D0 /
  DATA XGK (  3) / 0.930157491355708226001207180059508D0 /
  DATA XGK (  4) / 0.865063366688984510732096688423493D0 /
  DATA XGK (  5) / 0.780817726586416897063717578345042D0 /
  DATA XGK (  6) / 0.679409568299024406234327365114874D0 /
  DATA XGK (  7) / 0.562757134668604683339000099272694D0 /
  DATA XGK (  8) / 0.433395394129247190799265943165784D0 /
  DATA XGK (  9) / 0.294392862701460198131126603103866D0 /
  DATA XGK ( 10) / 0.148874338981631210884826001129720D0 /
  DATA XGK ( 11) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.011694638867371874278064396062192D0 /
  DATA WGK (  2) / 0.032558162307964727478818972459390D0 /
  DATA WGK (  3) / 0.054755896574351996031381300244580D0 /
  DATA WGK (  4) / 0.075039674810919952767043140916190D0 /
  DATA WGK (  5) / 0.093125454583697605535065465083366D0 /
  DATA WGK (  6) / 0.109387158802297641899210590325805D0 /
  DATA WGK (  7) / 0.123491976262065851077958109831074D0 /
  DATA WGK (  8) / 0.134709217311473325928054001771707D0 /
  DATA WGK (  9) / 0.142775938577060080797094273138717D0 /
  DATA WGK ( 10) / 0.147739104901338491374841515972068D0 /
  DATA WGK ( 11) / 0.149445554002916905664936468389821D0 /
  !
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       ABSC   - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
  !       RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
  !       RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                I.E. TO I/(B-A)
  !
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK21
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
  !       THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  RESG = 0.0D+00
  FC = F(CENTR)
  RESK = WGK(11)*FC
  RESABS = ABS(RESK)
  DO J=1,5
    JTW = 2*J
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J = 1,5
    JTWM1 = 2*J-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(11)*ABS(FC-RESKH)
  DO J=1,10
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK21
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK31(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK31
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  31-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B) with error
        !                 estimate
        !             J = Integral of ABS(F) over (A,B)
  !***DESCRIPTION
  !
  !       Integration rules
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS
  !        ON ENTRY
  !          F      - Double precision
  !                   Function subprogram defining the integrand
  !                   FUNCTION F(X). The actual name for F needs to be
  !                   Declared E X T E R N A L in the calling program.
  !
  !          A      - Double precision
  !                   Lower limit of integration
  !
  !          B      - Double precision
  !                   Upper limit of integration
  !
  !        ON RETURN
  !          RESULT - Double precision
  !                   Approximation to the integral I
  !                   RESULT is computed by applying the 31-POINT
  !                   GAUSS-KRONROD RULE (RESK), obtained by optimal
  !                   addition of abscissae to the 15-POINT GAUSS
  !                   RULE (RESG).
  !
  !          ABSERR - Double precison
  !                   Estimate of the modulus of the modulus,
  !                   which should not exceed ABS(I-RESULT)
  !
  !          RESABS - Double precision
  !                   Approximation to the integral J
  !
  !          RESASC - Double precision
  !                   Approximation to the integral of ABS(F-I/(B-A))
  !                   over (A,B)
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK31
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !       BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !       CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
  !                XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
  !                GAUSS RULE
  !                XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                ADDED TO THE 15-POINT GAUSS RULE
  !
  !       WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
  !
  !       WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.030753241996117268354628393577204D0 /
  DATA WG  (  2) / 0.070366047488108124709267416450667D0 /
  DATA WG  (  3) / 0.107159220467171935011869546685869D0 /
  DATA WG  (  4) / 0.139570677926154314447804794511028D0 /
  DATA WG  (  5) / 0.166269205816993933553200860481209D0 /
  DATA WG  (  6) / 0.186161000015562211026800561866423D0 /
  DATA WG  (  7) / 0.198431485327111576456118326443839D0 /
  DATA WG  (  8) / 0.202578241925561272880620199967519D0 /
  !
  DATA XGK (  1) / 0.998002298693397060285172840152271D0 /
  DATA XGK (  2) / 0.987992518020485428489565718586613D0 /
  DATA XGK (  3) / 0.967739075679139134257347978784337D0 /
  DATA XGK (  4) / 0.937273392400705904307758947710209D0 /
  DATA XGK (  5) / 0.897264532344081900882509656454496D0 /
  DATA XGK (  6) / 0.848206583410427216200648320774217D0 /
  DATA XGK (  7) / 0.790418501442465932967649294817947D0 /
  DATA XGK (  8) / 0.724417731360170047416186054613938D0 /
  DATA XGK (  9) / 0.650996741297416970533735895313275D0 /
  DATA XGK ( 10) / 0.570972172608538847537226737253911D0 /
  DATA XGK ( 11) / 0.485081863640239680693655740232351D0 /
  DATA XGK ( 12) / 0.394151347077563369897207370981045D0 /
  DATA XGK ( 13) / 0.299180007153168812166780024266389D0 /
  DATA XGK ( 14) / 0.201194093997434522300628303394596D0 /
  DATA XGK ( 15) / 0.101142066918717499027074231447392D0 /
  DATA XGK ( 16) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.005377479872923348987792051430128D0 /
  DATA WGK (  2) / 0.015007947329316122538374763075807D0 /
  DATA WGK (  3) / 0.025460847326715320186874001019653D0 /
  DATA WGK (  4) / 0.035346360791375846222037948478360D0 /
  DATA WGK (  5) / 0.044589751324764876608227299373280D0 /
  DATA WGK (  6) / 0.053481524690928087265343147239430D0 /
  DATA WGK (  7) / 0.062009567800670640285139230960803D0 /
  DATA WGK (  8) / 0.069854121318728258709520077099147D0 /
  DATA WGK (  9) / 0.076849680757720378894432777482659D0 /
  DATA WGK ( 10) / 0.083080502823133021038289247286104D0 /
  DATA WGK ( 11) / 0.088564443056211770647275443693774D0 /
  DATA WGK ( 12) / 0.093126598170825321225486872747346D0 /
  DATA WGK ( 13) / 0.096642726983623678505179907627589D0 /
  DATA WGK ( 14) / 0.099173598721791959332393173484603D0 /
  DATA WGK ( 15) / 0.100769845523875595044946662617570D0 /
  DATA WGK ( 16) / 0.101330007014791549017374792767493D0 /
  !
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       ABSC   - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
  !       RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
  !       RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                I.E. TO I/(B-A)
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !***FIRST EXECUTABLE STATEMENT  DQK31
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
  !       THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  FC = F(CENTR)
  RESG = WG(8)*FC
  RESK = WGK(16)*FC
  RESABS = ABS(RESK)
  DO J=1,7
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J = 1,8
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(16)*ABS(FC-RESKH)
  DO J=1,15
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK31
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK41(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK41
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  41-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B), with error
        !                 estimate
        !             J = Integral of ABS(F) over (A,B)
  !***DESCRIPTION
  !
  !       Integration rules
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS
  !        ON ENTRY
  !          F      - Double precision
  !                   Function subprogram defining the integrand
  !                   FUNCTION F(X). The actual name for F needs to be
  !                   declared E X T E R N A L in the calling program.
  !
  !          A      - Double precision
  !                   Lower limit of integration
  !
  !          B      - Double precision
  !                   Upper limit of integration
  !
  !        ON RETURN
  !          RESULT - Double precision
  !                   Approximation to the integral I
  !                   RESULT is computed by applying the 41-POINT
  !                   GAUSS-KRONROD RULE (RESK) obtained by optimal
  !                   addition of abscissae to the 20-POINT GAUSS
  !                   RULE (RESG).
  !
  !          ABSERR - Double precision
  !                   Estimate of the modulus of the absolute error,
  !                   which should not exceed ABS(I-RESULT)
  !
  !          RESABS - Double precision
  !                   Approximation to the integral J
  !
  !          RESASC - Double precision
  !                   Approximation to the integal of ABS(F-I/(B-A))
  !                   over (A,B)
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK41
  !
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(20),FV2(20),XGK(21),WGK(21),WG(10)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !       BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !       CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
  !                XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
  !                GAUSS RULE
  !                XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                ADDED TO THE 20-POINT GAUSS RULE
  !
  !       WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
  !
  !       WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.017614007139152118311861962351853D0 /
  DATA WG  (  2) / 0.040601429800386941331039952274932D0 /
  DATA WG  (  3) / 0.062672048334109063569506535187042D0 /
  DATA WG  (  4) / 0.083276741576704748724758143222046D0 /
  DATA WG  (  5) / 0.101930119817240435036750135480350D0 /
  DATA WG  (  6) / 0.118194531961518417312377377711382D0 /
  DATA WG  (  7) / 0.131688638449176626898494499748163D0 /
  DATA WG  (  8) / 0.142096109318382051329298325067165D0 /
  DATA WG  (  9) / 0.149172986472603746787828737001969D0 /
  DATA WG  ( 10) / 0.152753387130725850698084331955098D0 /
  !
  DATA XGK (  1) / 0.998859031588277663838315576545863D0 /
  DATA XGK (  2) / 0.993128599185094924786122388471320D0 /
  DATA XGK (  3) / 0.981507877450250259193342994720217D0 /
  DATA XGK (  4) / 0.963971927277913791267666131197277D0 /
  DATA XGK (  5) / 0.940822633831754753519982722212443D0 /
  DATA XGK (  6) / 0.912234428251325905867752441203298D0 /
  DATA XGK (  7) / 0.878276811252281976077442995113078D0 /
  DATA XGK (  8) / 0.839116971822218823394529061701521D0 /
  DATA XGK (  9) / 0.795041428837551198350638833272788D0 /
  DATA XGK ( 10) / 0.746331906460150792614305070355642D0 /
  DATA XGK ( 11) / 0.693237656334751384805490711845932D0 /
  DATA XGK ( 12) / 0.636053680726515025452836696226286D0 /
  DATA XGK ( 13) / 0.575140446819710315342946036586425D0 /
  DATA XGK ( 14) / 0.510867001950827098004364050955251D0 /
  DATA XGK ( 15) / 0.443593175238725103199992213492640D0 /
  DATA XGK ( 16) / 0.373706088715419560672548177024927D0 /
  DATA XGK ( 17) / 0.301627868114913004320555356858592D0 /
  DATA XGK ( 18) / 0.227785851141645078080496195368575D0 /
  DATA XGK ( 19) / 0.152605465240922675505220241022678D0 /
  DATA XGK ( 20) / 0.076526521133497333754640409398838D0 /
  DATA XGK ( 21) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.003073583718520531501218293246031D0 /
  DATA WGK (  2) / 0.008600269855642942198661787950102D0 /
  DATA WGK (  3) / 0.014626169256971252983787960308868D0 /
  DATA WGK (  4) / 0.020388373461266523598010231432755D0 /
  DATA WGK (  5) / 0.025882133604951158834505067096153D0 /
  DATA WGK (  6) / 0.031287306777032798958543119323801D0 /
  DATA WGK (  7) / 0.036600169758200798030557240707211D0 /
  DATA WGK (  8) / 0.041668873327973686263788305936895D0 /
  DATA WGK (  9) / 0.046434821867497674720231880926108D0 /
  DATA WGK ( 10) / 0.050944573923728691932707670050345D0 /
  DATA WGK ( 11) / 0.055195105348285994744832372419777D0 /
  DATA WGK ( 12) / 0.059111400880639572374967220648594D0 /
  DATA WGK ( 13) / 0.062653237554781168025870122174255D0 /
  DATA WGK ( 14) / 0.065834597133618422111563556969398D0 /
  DATA WGK ( 15) / 0.068648672928521619345623411885368D0 /
  DATA WGK ( 16) / 0.071054423553444068305790361723210D0 /
  DATA WGK ( 17) / 0.073030690332786667495189417658913D0 /
  DATA WGK ( 18) / 0.074582875400499188986581418362488D0 /
  DATA WGK ( 19) / 0.075704497684556674659542775376617D0 /
  DATA WGK ( 20) / 0.076377867672080736705502835038061D0 /
  DATA WGK ( 21) / 0.076600711917999656445049901530102D0 /
  !
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       ABSC   - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
  !       RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
  !       RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
  !                TO I/(B-A)
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK41
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
  !       THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  RESG = 0.0D+00
  FC = F(CENTR)
  RESK = WGK(21)*FC
  RESABS = ABS(RESK)
  DO J=1,10
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J = 1,10
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(21)*ABS(FC-RESKH)
  DO J=1,20
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK41
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK51(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK51
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
               ! 960627   Missing WGK(26) restored (RFB).
               ! 000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  51-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B) with error
        !                 estimate
        !             J = Integral of ABS(F) over (A,B)
  !***DESCRIPTION
  !
  !       Integration rules
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS
  !        ON ENTRY
  !          F      - Double precision
  !                   Function subroutine defining the integrand
  !                   function F(X). The actual name for F needs to be
  !                   declared E X T E R N A L in the calling program.
  !
  !          A      - Double precision
  !                   Lower limit of integration
  !
  !          B      - Double precision
  !                   Upper limit of integration
  !
  !        ON RETURN
  !          RESULT - Double precision
  !                   Approximation to the integral I
  !                   RESULT is computed by applying the 51-point
  !                   Kronrod rule (RESK) obtained by optimal addition
  !                   of abscissae to the 25-point Gauss rule (RESG).
  !
  !          ABSERR - Double precision
  !                   Estimate of the modulus of the absolute error,
  !                   which should not exceed ABS(I-RESULT)
  !
  !          RESABS - Double precision
  !                   Approximation to the integral J
  !
  !          RESASC - Double precision
  !                   Approximation to the integral of ABS(F-I/(B-A))
  !                   over (A,B)
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK51
  !
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(25),FV2(25),XGK(26),WGK(26),WG(13)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
  !       BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
  !       CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
  !                XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
  !                GAUSS RULE
  !                XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
  !                ADDED TO THE 25-POINT GAUSS RULE
  !
  !       WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
  !
  !       WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.011393798501026287947902964113235D0 /
  DATA WG  (  2) / 0.026354986615032137261901815295299D0 /
  DATA WG  (  3) / 0.040939156701306312655623487711646D0 /
  DATA WG  (  4) / 0.054904695975835191925936891540473D0 /
  DATA WG  (  5) / 0.068038333812356917207187185656708D0 /
  DATA WG  (  6) / 0.080140700335001018013234959669111D0 /
  DATA WG  (  7) / 0.091028261982963649811497220702892D0 /
  DATA WG  (  8) / 0.100535949067050644202206890392686D0 /
  DATA WG  (  9) / 0.108519624474263653116093957050117D0 /
  DATA WG  ( 10) / 0.114858259145711648339325545869556D0 /
  DATA WG  ( 11) / 0.119455763535784772228178126512901D0 /
  DATA WG  ( 12) / 0.122242442990310041688959518945852D0 /
  DATA WG  ( 13) / 0.123176053726715451203902873079050D0 /
  !
  DATA XGK (  1) / 0.999262104992609834193457486540341D0 /
  DATA XGK (  2) / 0.995556969790498097908784946893902D0 /
  DATA XGK (  3) / 0.988035794534077247637331014577406D0 /
  DATA XGK (  4) / 0.976663921459517511498315386479594D0 /
  DATA XGK (  5) / 0.961614986425842512418130033660167D0 /
  DATA XGK (  6) / 0.942974571228974339414011169658471D0 /
  DATA XGK (  7) / 0.920747115281701561746346084546331D0 /
  DATA XGK (  8) / 0.894991997878275368851042006782805D0 /
  DATA XGK (  9) / 0.865847065293275595448996969588340D0 /
  DATA XGK ( 10) / 0.833442628760834001421021108693570D0 /
  DATA XGK ( 11) / 0.797873797998500059410410904994307D0 /
  DATA XGK ( 12) / 0.759259263037357630577282865204361D0 /
  DATA XGK ( 13) / 0.717766406813084388186654079773298D0 /
  DATA XGK ( 14) / 0.673566368473468364485120633247622D0 /
  DATA XGK ( 15) / 0.626810099010317412788122681624518D0 /
  DATA XGK ( 16) / 0.577662930241222967723689841612654D0 /
  DATA XGK ( 17) / 0.526325284334719182599623778158010D0 /
  DATA XGK ( 18) / 0.473002731445714960522182115009192D0 /
  DATA XGK ( 19) / 0.417885382193037748851814394594572D0 /
  DATA XGK ( 20) / 0.361172305809387837735821730127641D0 /
  DATA XGK ( 21) / 0.303089538931107830167478909980339D0 /
  DATA XGK ( 22) / 0.243866883720988432045190362797452D0 /
  DATA XGK ( 23) / 0.183718939421048892015969888759528D0 /
  DATA XGK ( 24) / 0.122864692610710396387359818808037D0 /
  DATA XGK ( 25) / 0.061544483005685078886546392366797D0 /
  DATA XGK ( 26) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.001987383892330315926507851882843D0 /
  DATA WGK (  2) / 0.005561932135356713758040236901066D0 /
  DATA WGK (  3) / 0.009473973386174151607207710523655D0 /
  DATA WGK (  4) / 0.013236229195571674813656405846976D0 /
  DATA WGK (  5) / 0.016847817709128298231516667536336D0 /
  DATA WGK (  6) / 0.020435371145882835456568292235939D0 /
  DATA WGK (  7) / 0.024009945606953216220092489164881D0 /
  DATA WGK (  8) / 0.027475317587851737802948455517811D0 /
  DATA WGK (  9) / 0.030792300167387488891109020215229D0 /
  DATA WGK ( 10) / 0.034002130274329337836748795229551D0 /
  DATA WGK ( 11) / 0.037116271483415543560330625367620D0 /
  DATA WGK ( 12) / 0.040083825504032382074839284467076D0 /
  DATA WGK ( 13) / 0.042872845020170049476895792439495D0 /
  DATA WGK ( 14) / 0.045502913049921788909870584752660D0 /
  DATA WGK ( 15) / 0.047982537138836713906392255756915D0 /
  DATA WGK ( 16) / 0.050277679080715671963325259433440D0 /
  DATA WGK ( 17) / 0.052362885806407475864366712137873D0 /
  DATA WGK ( 18) / 0.054251129888545490144543370459876D0 /
  DATA WGK ( 19) / 0.055950811220412317308240686382747D0 /
  DATA WGK ( 20) / 0.057437116361567832853582693939506D0 /
  DATA WGK ( 21) / 0.058689680022394207961974175856788D0 /
  DATA WGK ( 22) / 0.059720340324174059979099291932562D0 /
  DATA WGK ( 23) / 0.060539455376045862945360267517565D0 /
  DATA WGK ( 24) / 0.061128509717053048305859030416293D0 /
  DATA WGK ( 25) / 0.061471189871425316661544131965264D0 /
  DATA WGK ( 26) / 0.061580818067832935078759824240055D0 /
  !
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       ABSC   - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
  !       RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
  !       RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
  !                I.E. TO I/(B-A)
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK51
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
  !       THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  FC = F(CENTR)
  RESG = WG(13)*FC
  RESK = WGK(26)*FC
  RESABS = ABS(RESK)
  DO J=1,12
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J = 1,13
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(26)*ABS(FC-RESKH)
  DO J=1,25
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK51
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQK61(F,A,B,RESULT,ABSERR,RESABS,RESASC)
  !***BEGIN PROLOGUE  DQK61
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  830518   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A2
  !***KEYWORDS  61-POINT GAUSS-KRONROD RULES
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  To compute I = Integral of F over (A,B) with error
        !                 estimate
        !             J = Integral of DABS(F) over (A,B)
  !***DESCRIPTION
  !
  !    Integration rule
  !    Standard fortran subroutine
  !    Double precision version
  !
  !
  !    PARAMETERS
  !     ON ENTRY
  !       F      - Double precision
  !                Function subprogram defining the integrand
  !                function F(X). The actual name for F needs to be
  !                declared E X T E R N A L in the calling program.
  !
  !       A      - Double precision
  !                Lower limit of integration
  !
  !       B      - Double precision
  !                Upper limit of integration
  !
  !     ON RETURN
  !       RESULT - Double precision
  !                Approximation to the integral I
  !                RESULT is computed by applying the 61-point
  !                Kronrod rule (RESK) obtained by optimal addition of
  !                abscissae to the 30-point Gauss rule (RESG).
  !
  !       ABSERR - Double precision
  !                Estimate of the modulus of the absolute error,
  !                which should equal or exceed DABS(I-RESULT)
  !
  !       RESABS - Double precision
  !                Approximation to the integral J
  !
  !       RESASC - Double precision
  !                Approximation to the integral of DABS(F-I/(B-A))
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***END PROLOGUE  DQK61
  !
  DOUBLE PRECISION :: A,DABSC,ABSERR,B,CENTR,ABS,DHLGTH, &
        EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
        RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER :: J,JTW,JTWM1
  EXTERNAL F
  !
  DIMENSION FV1(30),FV2(30),XGK(31),WGK(31),WG(15)
  !
  !       THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
  !       INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
  !       ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
  !
  !       XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
  !               XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
  !               GAUSS RULE
  !               XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
  !               TO THE 30-POINT GAUSS RULE
  !
  !       WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
  !
  !       WG    - WEIGTHS OF THE 30-POINT GAUSS RULE
  !
  !
  ! GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
  ! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
  ! BELL LABS, NOV. 1981.
  !
  DATA WG  (  1) / 0.007968192496166605615465883474674D0 /
  DATA WG  (  2) / 0.018466468311090959142302131912047D0 /
  DATA WG  (  3) / 0.028784707883323369349719179611292D0 /
  DATA WG  (  4) / 0.038799192569627049596801936446348D0 /
  DATA WG  (  5) / 0.048402672830594052902938140422808D0 /
  DATA WG  (  6) / 0.057493156217619066481721689402056D0 /
  DATA WG  (  7) / 0.065974229882180495128128515115962D0 /
  DATA WG  (  8) / 0.073755974737705206268243850022191D0 /
  DATA WG  (  9) / 0.080755895229420215354694938460530D0 /
  DATA WG  ( 10) / 0.086899787201082979802387530715126D0 /
  DATA WG  ( 11) / 0.092122522237786128717632707087619D0 /
  DATA WG  ( 12) / 0.096368737174644259639468626351810D0 /
  DATA WG  ( 13) / 0.099593420586795267062780282103569D0 /
  DATA WG  ( 14) / 0.101762389748405504596428952168554D0 /
  DATA WG  ( 15) / 0.102852652893558840341285636705415D0 /
  !
  DATA XGK (  1) / 0.999484410050490637571325895705811D0 /
  DATA XGK (  2) / 0.996893484074649540271630050918695D0 /
  DATA XGK (  3) / 0.991630996870404594858628366109486D0 /
  DATA XGK (  4) / 0.983668123279747209970032581605663D0 /
  DATA XGK (  5) / 0.973116322501126268374693868423707D0 /
  DATA XGK (  6) / 0.960021864968307512216871025581798D0 /
  DATA XGK (  7) / 0.944374444748559979415831324037439D0 /
  DATA XGK (  8) / 0.926200047429274325879324277080474D0 /
  DATA XGK (  9) / 0.905573307699907798546522558925958D0 /
  DATA XGK ( 10) / 0.882560535792052681543116462530226D0 /
  DATA XGK ( 11) / 0.857205233546061098958658510658944D0 /
  DATA XGK ( 12) / 0.829565762382768397442898119732502D0 /
  DATA XGK ( 13) / 0.799727835821839083013668942322683D0 /
  DATA XGK ( 14) / 0.767777432104826194917977340974503D0 /
  DATA XGK ( 15) / 0.733790062453226804726171131369528D0 /
  DATA XGK ( 16) / 0.697850494793315796932292388026640D0 /
  DATA XGK ( 17) / 0.660061064126626961370053668149271D0 /
  DATA XGK ( 18) / 0.620526182989242861140477556431189D0 /
  DATA XGK ( 19) / 0.579345235826361691756024932172540D0 /
  DATA XGK ( 20) / 0.536624148142019899264169793311073D0 /
  DATA XGK ( 21) / 0.492480467861778574993693061207709D0 /
  DATA XGK ( 22) / 0.447033769538089176780609900322854D0 /
  DATA XGK ( 23) / 0.400401254830394392535476211542661D0 /
  DATA XGK ( 24) / 0.352704725530878113471037207089374D0 /
  DATA XGK ( 25) / 0.304073202273625077372677107199257D0 /
  DATA XGK ( 26) / 0.254636926167889846439805129817805D0 /
  DATA XGK ( 27) / 0.204525116682309891438957671002025D0 /
  DATA XGK ( 28) / 0.153869913608583546963794672743256D0 /
  DATA XGK ( 29) / 0.102806937966737030147096751318001D0 /
  DATA XGK ( 30) / 0.051471842555317695833025213166723D0 /
  DATA XGK ( 31) / 0.000000000000000000000000000000000D0 /
  !
  DATA WGK (  1) / 0.001389013698677007624551591226760D0 /
  DATA WGK (  2) / 0.003890461127099884051267201844516D0 /
  DATA WGK (  3) / 0.006630703915931292173319826369750D0 /
  DATA WGK (  4) / 0.009273279659517763428441146892024D0 /
  DATA WGK (  5) / 0.011823015253496341742232898853251D0 /
  DATA WGK (  6) / 0.014369729507045804812451432443580D0 /
  DATA WGK (  7) / 0.016920889189053272627572289420322D0 /
  DATA WGK (  8) / 0.019414141193942381173408951050128D0 /
  DATA WGK (  9) / 0.021828035821609192297167485738339D0 /
  DATA WGK ( 10) / 0.024191162078080601365686370725232D0 /
  DATA WGK ( 11) / 0.026509954882333101610601709335075D0 /
  DATA WGK ( 12) / 0.028754048765041292843978785354334D0 /
  DATA WGK ( 13) / 0.030907257562387762472884252943092D0 /
  DATA WGK ( 14) / 0.032981447057483726031814191016854D0 /
  DATA WGK ( 15) / 0.034979338028060024137499670731468D0 /
  DATA WGK ( 16) / 0.036882364651821229223911065617136D0 /
  DATA WGK ( 17) / 0.038678945624727592950348651532281D0 /
  DATA WGK ( 18) / 0.040374538951535959111995279752468D0 /
  DATA WGK ( 19) / 0.041969810215164246147147541285970D0 /
  DATA WGK ( 20) / 0.043452539701356069316831728117073D0 /
  DATA WGK ( 21) / 0.044814800133162663192355551616723D0 /
  DATA WGK ( 22) / 0.046059238271006988116271735559374D0 /
  DATA WGK ( 23) / 0.047185546569299153945261478181099D0 /
  DATA WGK ( 24) / 0.048185861757087129140779492298305D0 /
  DATA WGK ( 25) / 0.049055434555029778887528165367238D0 /
  DATA WGK ( 26) / 0.049795683427074206357811569379942D0 /
  DATA WGK ( 27) / 0.050405921402782346840893085653585D0 /
  DATA WGK ( 28) / 0.050881795898749606492297473049805D0 /
  DATA WGK ( 29) / 0.051221547849258772170656282604944D0 /
  DATA WGK ( 30) / 0.051426128537459025933862879215781D0 /
  DATA WGK ( 31) / 0.051494729429451567558340433647099D0 /
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTERVAL
  !       DABSC  - ABSCISSA
  !       FVAL*  - FUNCTION VALUE
  !       RESG   - RESULT OF THE 30-POINT GAUSS RULE
  !       RESK   - RESULT OF THE 61-POINT KRONROD RULE
  !       RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
  !                OVER (A,B), I.E. TO I/(B-A)
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQK61
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  CENTR = 0.5D+00*(B+A)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  !
  !       COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
  !       INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
  !
  RESG = 0.0D+00
  FC = F(CENTR)
  RESK = WGK(31)*FC
  RESABS = ABS(RESK)
  DO J=1,15
    JTW = J*2
    DABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-DABSC)
    FVAL2 = F(CENTR+DABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  DO J=1,15
    JTWM1 = J*2-1
    DABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-DABSC)
    FVAL2 = F(CENTR+DABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  END DO
  RESKH = RESK*0.5D+00
  RESASC = WGK(31)*ABS(FC-RESKH)
  DO J=1,30
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
  END DO
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
        ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  RETURN
  END SUBROUTINE DQK61
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQNG(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER)
  !***BEGIN PROLOGUE  DQNG
  !***DATE WRITTEN   800101   (YYMMDD)
  !***REVISION DATE  810101   (YYMMDD)
  !***REVISION HISTORY (YYMMDD)
  !   000601   Changed DMAX1/DMIN1/DABS to generic MAX/MIN/ABS
  !***CATEGORY NO.  H2A1A1
  !***KEYWORDS  AUTOMATIC INTEGRATOR,GAUSS-KRONROD(PATTERSON),
          ! NON-ADAPTIVE,SMOOTH INTEGRAND
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
        ! KAHANER, DAVID, NBS - MODIFIED (2/82)
  !***PURPOSE  The routine calculates an approximation result to a
        !  given definite integral I = integral of F over (A,B),
        !  hopefully satisfying following claim for accuracy
        !  ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***DESCRIPTION
  !
  ! NON-ADAPTIVE INTEGRATION
  ! STANDARD FORTRAN SUBROUTINE
  ! DOUBLE PRECISION VERSION
  !
  !       F      - Double precision
  !                Function subprogram defining the integrand function
  !                F(X). The actual name for F needs to be declared
  !                E X T E R N A L in the driver program.
  !
  !       A      - Double precision
  !                Lower limit of integration
  !
  !       B      - Double precision
  !                Upper limit of integration
  !
  !       EPSABS - Double precision
  !                Absolute accuracy requested
  !       EPSREL - Double precision
  !                Relative accuracy requested
  !                If  EPSABS.LE.0
  !                And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                The routine will end with IER = 6.
  !
  !     ON RETURN
  !       RESULT - Double precision
  !                Approximation to the integral I
  !                Result is obtained by applying the 21-POINT
  !                GAUSS-KRONROD RULE (RES21) obtained by optimal
  !                addition of abscissae to the 10-POINT GAUSS RULE
  !                (RES10), or by applying the 43-POINT RULE (RES43)
  !                obtained by optimal addition of abscissae to the
  !                21-POINT GAUSS-KRONROD RULE, or by applying the
  !                87-POINT RULE (RES87) obtained by optimal addition
  !                of abscissae to the 43-POINT RULE.
  !
  !       ABSERR - Double precision
  !                Estimate of the modulus of the absolute error,
  !                which should EQUAL or EXCEED ABS(I-RESULT)
  !
  !       NEVAL  - Integer
  !                Number of integrand evaluations
  !
  !       IER    - IER = 0 normal and reliable termination of the
  !                        routine. It is assumed that the requested
  !                        accuracy has been achieved.
  !                IER.GT.0 Abnormal termination of the routine. It is
  !                            assumed that the requested accuracy has
  !                        not been achieved.
  !       ERROR MESSAGES
  !                IER = 1 The maximum number of steps has been
  !                        executed. The integral is probably too
  !                        difficult to be calculated by DQNG.
  !                    = 6 The input is invalid, because
  !                        EPSABS.LE.0 AND
  !                        EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
  !                        RESULT, ABSERR and NEVAL are set to zero.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH,XERROR
  !***END PROLOGUE  DQNG
  !
  DOUBLE PRECISION :: A,ABSC,ABSERR,B,CENTR,DHLGTH, &
        EPMACH,EPSABS,EPSREL,F,FCENTR,FVAL,FVAL1,FVAL2,FV1,FV2, &
        FV3,FV4,HLGTH,RESULT,RES10,RES21,RES43,RES87,RESABS,RESASC, &
        RESKH,SAVFUN,UFLOW,W10,W21A,W21B,W43A,W43B,W87A,W87B,X1,X2,X3,X4
  INTEGER :: IER,IPX,K,L,NEVAL
  EXTERNAL F
  !
  DIMENSION FV1(5),FV2(5),FV3(5),FV4(5),X1(5),X2(5),X3(11),X4(22), &
        W10(5),W21A(5),W21B(6),W43A(10),W43B(12),W87A(21),W87B(23), &
        SAVFUN(21)
  !
  !       THE FOLLOWING DATA STATEMENTS CONTAIN THE
  !       ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED.
  !
  !       X1      ABSCISSAE COMMON TO THE 10-, 21-, 43- AND 87-
  !               POINT RULE
  !       X2      ABSCISSAE COMMON TO THE 21-, 43- AND 87-POINT RULE
  !       X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT RULE
  !       X4      ABSCISSAE OF THE 87-POINT RULE
  !       W10     WEIGHTS OF THE 10-POINT FORMULA
  !       W21A    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X1
  !       W21B    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X2
  !       W43A    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X1, X3
  !           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X3
  !       W87A    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X1,
  !               X2, X3
  !       W87B    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X4
  !
  !
  ! GAUSS-KRONROD-PATTERSON QUADRATURE COEFFICIENTS FOR USE IN
  ! QUADPACK ROUTINE QNG.  THESE COEFFICIENTS WERE CALCULATED WITH
  ! 101 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, BELL LABS, NOV 1981.
  !
  DATA X1    (  1) / 0.973906528517171720077964012084452D0 /
  DATA X1    (  2) / 0.865063366688984510732096688423493D0 /
  DATA X1    (  3) / 0.679409568299024406234327365114874D0 /
  DATA X1    (  4) / 0.433395394129247190799265943165784D0 /
  DATA X1    (  5) / 0.148874338981631210884826001129720D0 /
  DATA W10   (  1) / 0.066671344308688137593568809893332D0 /
  DATA W10   (  2) / 0.149451349150580593145776339657697D0 /
  DATA W10   (  3) / 0.219086362515982043995534934228163D0 /
  DATA W10   (  4) / 0.269266719309996355091226921569469D0 /
  DATA W10   (  5) / 0.295524224714752870173892994651338D0 /
  !
  DATA X2    (  1) / 0.995657163025808080735527280689003D0 /
  DATA X2    (  2) / 0.930157491355708226001207180059508D0 /
  DATA X2    (  3) / 0.780817726586416897063717578345042D0 /
  DATA X2    (  4) / 0.562757134668604683339000099272694D0 /
  DATA X2    (  5) / 0.294392862701460198131126603103866D0 /
  DATA W21A  (  1) / 0.032558162307964727478818972459390D0 /
  DATA W21A  (  2) / 0.075039674810919952767043140916190D0 /
  DATA W21A  (  3) / 0.109387158802297641899210590325805D0 /
  DATA W21A  (  4) / 0.134709217311473325928054001771707D0 /
  DATA W21A  (  5) / 0.147739104901338491374841515972068D0 /
  DATA W21B  (  1) / 0.011694638867371874278064396062192D0 /
  DATA W21B  (  2) / 0.054755896574351996031381300244580D0 /
  DATA W21B  (  3) / 0.093125454583697605535065465083366D0 /
  DATA W21B  (  4) / 0.123491976262065851077958109831074D0 /
  DATA W21B  (  5) / 0.142775938577060080797094273138717D0 /
  DATA W21B  (  6) / 0.149445554002916905664936468389821D0 /
  !
  DATA X3    (  1) / 0.999333360901932081394099323919911D0 /
  DATA X3    (  2) / 0.987433402908088869795961478381209D0 /
  DATA X3    (  3) / 0.954807934814266299257919200290473D0 /
  DATA X3    (  4) / 0.900148695748328293625099494069092D0 /
  DATA X3    (  5) / 0.825198314983114150847066732588520D0 /
  DATA X3    (  6) / 0.732148388989304982612354848755461D0 /
  DATA X3    (  7) / 0.622847970537725238641159120344323D0 /
  DATA X3    (  8) / 0.499479574071056499952214885499755D0 /
  DATA X3    (  9) / 0.364901661346580768043989548502644D0 /
  DATA X3    ( 10) / 0.222254919776601296498260928066212D0 /
  DATA X3    ( 11) / 0.074650617461383322043914435796506D0 /
  DATA W43A  (  1) / 0.016296734289666564924281974617663D0 /
  DATA W43A  (  2) / 0.037522876120869501461613795898115D0 /
  DATA W43A  (  3) / 0.054694902058255442147212685465005D0 /
  DATA W43A  (  4) / 0.067355414609478086075553166302174D0 /
  DATA W43A  (  5) / 0.073870199632393953432140695251367D0 /
  DATA W43A  (  6) / 0.005768556059769796184184327908655D0 /
  DATA W43A  (  7) / 0.027371890593248842081276069289151D0 /
  DATA W43A  (  8) / 0.046560826910428830743339154433824D0 /
  DATA W43A  (  9) / 0.061744995201442564496240336030883D0 /
  DATA W43A  ( 10) / 0.071387267268693397768559114425516D0 /
  DATA W43B  (  1) / 0.001844477640212414100389106552965D0 /
  DATA W43B  (  2) / 0.010798689585891651740465406741293D0 /
  DATA W43B  (  3) / 0.021895363867795428102523123075149D0 /
  DATA W43B  (  4) / 0.032597463975345689443882222526137D0 /
  DATA W43B  (  5) / 0.042163137935191811847627924327955D0 /
  DATA W43B  (  6) / 0.050741939600184577780189020092084D0 /
  DATA W43B  (  7) / 0.058379395542619248375475369330206D0 /
  DATA W43B  (  8) / 0.064746404951445885544689259517511D0 /
  DATA W43B  (  9) / 0.069566197912356484528633315038405D0 /
  DATA W43B  ( 10) / 0.072824441471833208150939535192842D0 /
  DATA W43B  ( 11) / 0.074507751014175118273571813842889D0 /
  DATA W43B  ( 12) / 0.074722147517403005594425168280423D0 /
  !
  DATA X4    (  1) / 0.999902977262729234490529830591582D0 /
  DATA X4    (  2) / 0.997989895986678745427496322365960D0 /
  DATA X4    (  3) / 0.992175497860687222808523352251425D0 /
  DATA X4    (  4) / 0.981358163572712773571916941623894D0 /
  DATA X4    (  5) / 0.965057623858384619128284110607926D0 /
  DATA X4    (  6) / 0.943167613133670596816416634507426D0 /
  DATA X4    (  7) / 0.915806414685507209591826430720050D0 /
  DATA X4    (  8) / 0.883221657771316501372117548744163D0 /
  DATA X4    (  9) / 0.845710748462415666605902011504855D0 /
  DATA X4    ( 10) / 0.803557658035230982788739474980964D0 /
  DATA X4    ( 11) / 0.757005730685495558328942793432020D0 /
  DATA X4    ( 12) / 0.706273209787321819824094274740840D0 /
  DATA X4    ( 13) / 0.651589466501177922534422205016736D0 /
  DATA X4    ( 14) / 0.593223374057961088875273770349144D0 /
  DATA X4    ( 15) / 0.531493605970831932285268948562671D0 /
  DATA X4    ( 16) / 0.466763623042022844871966781659270D0 /
  DATA X4    ( 17) / 0.399424847859218804732101665817923D0 /
  DATA X4    ( 18) / 0.329874877106188288265053371824597D0 /
  DATA X4    ( 19) / 0.258503559202161551802280975429025D0 /
  DATA X4    ( 20) / 0.185695396568346652015917141167606D0 /
  DATA X4    ( 21) / 0.111842213179907468172398359241362D0 /
  DATA X4    ( 22) / 0.037352123394619870814998165437704D0 /
  DATA W87A  (  1) / 0.008148377384149172900002878448190D0 /
  DATA W87A  (  2) / 0.018761438201562822243935059003794D0 /
  DATA W87A  (  3) / 0.027347451050052286161582829741283D0 /
  DATA W87A  (  4) / 0.033677707311637930046581056957588D0 /
  DATA W87A  (  5) / 0.036935099820427907614589586742499D0 /
  DATA W87A  (  6) / 0.002884872430211530501334156248695D0 /
  DATA W87A  (  7) / 0.013685946022712701888950035273128D0 /
  DATA W87A  (  8) / 0.023280413502888311123409291030404D0 /
  DATA W87A  (  9) / 0.030872497611713358675466394126442D0 /
  DATA W87A  ( 10) / 0.035693633639418770719351355457044D0 /
  DATA W87A  ( 11) / 0.000915283345202241360843392549948D0 /
  DATA W87A  ( 12) / 0.005399280219300471367738743391053D0 /
  DATA W87A  ( 13) / 0.010947679601118931134327826856808D0 /
  DATA W87A  ( 14) / 0.016298731696787335262665703223280D0 /
  DATA W87A  ( 15) / 0.021081568889203835112433060188190D0 /
  DATA W87A  ( 16) / 0.025370969769253827243467999831710D0 /
  DATA W87A  ( 17) / 0.029189697756475752501446154084920D0 /
  DATA W87A  ( 18) / 0.032373202467202789685788194889595D0 /
  DATA W87A  ( 19) / 0.034783098950365142750781997949596D0 /
  DATA W87A  ( 20) / 0.036412220731351787562801163687577D0 /
  DATA W87A  ( 21) / 0.037253875503047708539592001191226D0 /
  DATA W87B  (  1) / 0.000274145563762072350016527092881D0 /
  DATA W87B  (  2) / 0.001807124155057942948341311753254D0 /
  DATA W87B  (  3) / 0.004096869282759164864458070683480D0 /
  DATA W87B  (  4) / 0.006758290051847378699816577897424D0 /
  DATA W87B  (  5) / 0.009549957672201646536053581325377D0 /
  DATA W87B  (  6) / 0.012329447652244853694626639963780D0 /
  DATA W87B  (  7) / 0.015010447346388952376697286041943D0 /
  DATA W87B  (  8) / 0.017548967986243191099665352925900D0 /
  DATA W87B  (  9) / 0.019938037786440888202278192730714D0 /
  DATA W87B  ( 10) / 0.022194935961012286796332102959499D0 /
  DATA W87B  ( 11) / 0.024339147126000805470360647041454D0 /
  DATA W87B  ( 12) / 0.026374505414839207241503786552615D0 /
  DATA W87B  ( 13) / 0.028286910788771200659968002987960D0 /
  DATA W87B  ( 14) / 0.030052581128092695322521110347341D0 /
  DATA W87B  ( 15) / 0.031646751371439929404586051078883D0 /
  DATA W87B  ( 16) / 0.033050413419978503290785944862689D0 /
  DATA W87B  ( 17) / 0.034255099704226061787082821046821D0 /
  DATA W87B  ( 18) / 0.035262412660156681033782717998428D0 /
  DATA W87B  ( 19) / 0.036076989622888701185500318003895D0 /
  DATA W87B  ( 20) / 0.036698604498456094498018047441094D0 /
  DATA W87B  ( 21) / 0.037120549269832576114119958413599D0 /
  DATA W87B  ( 22) / 0.037334228751935040321235449094698D0 /
  DATA W87B  ( 23) / 0.037361073762679023410321241766599D0 /
  !
  !       LIST OF MAJOR VARIABLES
  !       -----------------------
  !
  !       CENTR  - MID POINT OF THE INTEGRATION INTERVAL
  !       HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
  !       FCENTR - FUNCTION VALUE AT MID POINT
  !       ABSC   - ABSCISSA
  !       FVAL   - FUNCTION VALUE
  !       SAVFUN - ARRAY OF FUNCTION VALUES WHICH HAVE ALREADY BEEN
  !                COMPUTED
  !       RES10  - 10-POINT GAUSS RESULT
  !       RES21  - 21-POINT KRONROD RESULT
  !       RES43  - 43-POINT RESULT
  !       RES87  - 87-POINT RESULT
  !       RESABS - APPROXIMATION TO THE INTEGRAL OF ABS(F)
  !       RESASC - APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
  !
  !       MACHINE DEPENDENT CONSTANTS
  !       ---------------------------
  !
  !       EPMACH IS THE LARGEST RELATIVE SPACING.
  !       UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQNG
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  !
  !       TEST ON VALIDITY OF PARAMETERS
  !       ------------------------------
  !
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  NEVAL = 0
  IER = 6
  IF(EPSABS.LE.0.0D+00.AND.EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)) &
        GO TO 80
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  CENTR = 0.5D+00*(B+A)
  FCENTR = F(CENTR)
  NEVAL = 21
  IER = 1
  !
  !      COMPUTE THE INTEGRAL USING THE 10- AND 21-POINT FORMULA.
  !
  DO L = 1,3
  GO TO (5,25,45),L
    5   RES10 = 0.0D+00
  RES21 = W21B(6)*FCENTR
  RESABS = W21B(6)*ABS(FCENTR)
  DO K=1,5
    ABSC = HLGTH*X1(K)
    FVAL1 = F(CENTR+ABSC)
    FVAL2 = F(CENTR-ABSC)
    FVAL = FVAL1+FVAL2
    RES10 = RES10+W10(K)*FVAL
    RES21 = RES21+W21A(K)*FVAL
    RESABS = RESABS+W21A(K)*(ABS(FVAL1)+ABS(FVAL2))
    SAVFUN(K) = FVAL
    FV1(K) = FVAL1
    FV2(K) = FVAL2
  END DO
  IPX = 5
  DO K=1,5
    IPX = IPX+1
    ABSC = HLGTH*X2(K)
    FVAL1 = F(CENTR+ABSC)
    FVAL2 = F(CENTR-ABSC)
    FVAL = FVAL1+FVAL2
    RES21 = RES21+W21B(K)*FVAL
    RESABS = RESABS+W21B(K)*(ABS(FVAL1)+ABS(FVAL2))
    SAVFUN(IPX) = FVAL
    FV3(K) = FVAL1
    FV4(K) = FVAL2
  END DO
  !
  !      TEST FOR CONVERGENCE.
  !
  RESULT = RES21*HLGTH
  RESABS = RESABS*DHLGTH
  RESKH = 0.5D+00*RES21
  RESASC = W21B(6)*ABS(FCENTR-RESKH)
  DO K = 1,5
    RESASC = RESASC+W21A(K)*(ABS(FV1(K)-RESKH)+ABS(FV2(K)-RESKH)) &
          +W21B(K)*(ABS(FV3(K)-RESKH)+ABS(FV4(K)-RESKH))
  END DO
  ABSERR = ABS((RES21-RES10)*HLGTH)
  RESASC = RESASC*DHLGTH
  GO TO 65
  !
  !      COMPUTE THE INTEGRAL USING THE 43-POINT FORMULA.
  !
   25   RES43 = W43B(12)*FCENTR
  NEVAL = 43
  DO K=1,10
    RES43 = RES43+SAVFUN(K)*W43A(K)
  END DO
  DO K=1,11
    IPX = IPX+1
    ABSC = HLGTH*X3(K)
    FVAL = F(ABSC+CENTR)+F(CENTR-ABSC)
    RES43 = RES43+FVAL*W43B(K)
    SAVFUN(IPX) = FVAL
  END DO
  !
  !      TEST FOR CONVERGENCE.
  !
  RESULT = RES43*HLGTH
  ABSERR = ABS((RES43-RES21)*HLGTH)
  GO TO 65
  !
  !      COMPUTE THE INTEGRAL USING THE 87-POINT FORMULA.
  !
   45   RES87 = W87B(23)*FCENTR
  NEVAL = 87
  DO K=1,21
    RES87 = RES87+SAVFUN(K)*W87A(K)
  END DO
  DO K=1,22
    ABSC = HLGTH*X4(K)
    RES87 = RES87+W87B(K)*(F(ABSC+CENTR)+F(CENTR-ABSC))
  END DO
  RESULT = RES87*HLGTH
  ABSERR = ABS((RES87-RES43)*HLGTH)
   65   IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) &
              ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  IF (RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
        ((EPMACH*0.5D+02)*RESABS,ABSERR)
  IF (ABSERR.LE.MAX(EPSABS,EPSREL*ABS(RESULT))) IER = 0
  ! ***JUMP OUT OF DO-LOOP
  IF (IER.EQ.0) GO TO 999
  END DO
   80   CALL XERROR( 'ABNORMAL RETURN FROM DQNG ',26,IER,0)
  999   RETURN
  END SUBROUTINE DQNG
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DQPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
  !***BEGIN PROLOGUE  DQPSRT
  !***REFER TO  DQAGE,DQAGIE,DQAGPE,DQAWSE
  !***ROUTINES CALLED  (NONE)
  !***REVISION DATE  810101   (YYMMDD)
  !***KEYWORDS  SEQUENTIAL SORTING
  !***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
          ! K. U. LEUVEN
        ! DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
        !   K. U. LEUVEN
  !***PURPOSE  This routine maintains the descending ordering in the
        !  list of the local error estimated resulting from the
        !  interval subdivision process. At each call two error
        !  estimates are inserted using the sequential search
        !  method, top-down for the largest error estimate and
        !  bottom-up for the smallest error estimate.
  !***DESCRIPTION
  !
  !       Ordering routine
  !       Standard fortran subroutine
  !       Double precision version
  !
  !       PARAMETERS (MEANING AT OUTPUT)
  !          LIMIT  - Integer
  !                   Maximum number of error estimates the list
  !                   can contain
  !
  !          LAST   - Integer
  !                   Number of error estimates currently in the list
  !
  !          MAXERR - Integer
  !                   Maxerr points to the NRMAX-th largest error
  !                   estimate currently in the list
  !
  !          ERMAX  - Double precision
  !                   NRMAX-th largest error estimate
  !                   ERMAX = ELIST(MAXERR)
  !
  !          ELIST  - Double precision
  !                   Vector of dimension LAST containing
  !                   the error estimates
  !
  !          IORD   - Integer
  !                   Vector of dimension LAST, the first K elements
  !                   of which contain pointers to the error
  !                   estimates, such that
  !                   ELIST(IORD(1)),...,  ELIST(IORD(K))
  !                   form a decreasing sequence, with
  !                   K = LAST if LAST.LE.(LIMIT/2+2), and
  !                   K = LIMIT+1-LAST otherwise
  !
  !          NRMAX  - Integer
  !                   MAXERR = IORD(NRMAX)
  !***END PROLOGUE  DQPSRT
  !
  DOUBLE PRECISION :: ELIST,ERMAX,ERRMAX,ERRMIN
  INTEGER :: I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR, &
        NRMAX
  DIMENSION ELIST(LAST),IORD(LAST)
  !
  !       CHECK WHETHER THE LIST CONTAINS MORE THAN
  !       TWO ERROR ESTIMATES.
  !
  !***FIRST EXECUTABLE STATEMENT  DQPSRT
  IF(LAST.GT.2) GO TO 10
  IORD(1) = 1
  IORD(2) = 2
  GO TO 90
  !
  !       THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
  !       DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
  !       ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
  !       START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
  !
   10   ERRMAX = ELIST(MAXERR)
  IF(NRMAX.EQ.1) GO TO 30
  IDO = NRMAX-1
  DO I = 1,IDO
    ISUCC = IORD(NRMAX-1)
  ! ***JUMP OUT OF DO-LOOP
    IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
    IORD(NRMAX) = ISUCC
    NRMAX = NRMAX-1
  END DO
  !
  !       COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
  !       IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
  !       SUBDIVISIONS STILL ALLOWED.
  !
   30   JUPBN = LAST
  IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
  ERRMIN = ELIST(LAST)
  !
  !       INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
  !       STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
  !
  JBND = JUPBN-1
  IBEG = NRMAX+1
  IF(IBEG.GT.JBND) GO TO 50
  DO I=IBEG,JBND
    ISUCC = IORD(I)
  ! ***JUMP OUT OF DO-LOOP
    IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
    IORD(I-1) = ISUCC
  END DO
   50   IORD(JBND) = MAXERR
  IORD(JUPBN) = LAST
  GO TO 90
  !
  !       INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
  !
   60   IORD(I-1) = MAXERR
  K = JBND
  DO J=I,JBND
    ISUCC = IORD(K)
  ! ***JUMP OUT OF DO-LOOP
    IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
    IORD(K+1) = ISUCC
    K = K-1
  END DO
  IORD(I) = LAST
  GO TO 90
   80   IORD(K+1) = LAST
  !
  !       SET MAXERR AND ERMAX.
  !
   90   MAXERR = IORD(NRMAX)
  ERMAX = ELIST(MAXERR)
  RETURN
  END SUBROUTINE DQPSRT
  !-----------------------------------------------------------------------------

END MODULE MOOSE_CMLIB_QUADPKD
