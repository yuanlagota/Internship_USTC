!===============================================================================
!
! P A C K A G E      DTENSBS
!
! (Version 1982 )
!
!
! Subprograms for interpolation of two and three  dimensional  gridded
! data using tensor products of B-spline basis functions.  This  is  a
! double precision version of the package TENSBS.
!
! By two dimensional gridded data we mean data of the form
!
!          (x(i), y(j), f(x(i),y(j)))  i=1,..,nx, j=1,..,ny.
!
! The subprograms in this package  determine  a  piecewise  polynomial
! function S(x,y) such that
!
!         S(x(i),y(j)) = f(x(i),y(j))  i=1,..,nx, j=1,..,ny.
!
! The function S takes the form
!
!                         nx   ny
!             S(x,y)  =  SUM  SUM  a   U (x) V (y)
!                        i=1  j=1   ij  i     j
!
! where U(i) and V(j) are fixed one-dimensional  piecewise  polynomial
! functions (the B-spline basis functions of the reference). The  user
! specifies the order (degree+1) of the polynomial pieces that  define
! the function S in each direction.  The  resulting  interpolant  will
! have continuous derivatives of up to order-2 in each direction.  For
! example, if the user specifies order 4 in x and order 3 in  y,  then
! the functions U(i) will be piecewise  cubic  polynomials  while  the
! functions  V(j)  will  be  piecewise   quadratics.   The   resulting
! interpolating function will have continuous first and second partial
! derivatives  with  respect  to  x  and  continuous   first   partial
! derivative with respect to y. (Lower continuity can be  obtained  by
! using the option for user-specified "knots" -- see the reference.)
!
! The subroutines in this package are
!
!
! DB2INK.........computes   parameters   that   define   a   piecewise
!                polynomial function that interpolates a given set  of
!                two-dimensional gridded data.
!
! DB2VAL.........evaluates the interpolating  function  determined  by
!                DB2INK or one of its derivatives.
!
! DB3INK.........computes   parameters   that   define   a   piecewise
!                polynomial function that interpolates a given set  of
!                three-dimensional gridded data.
!
! DB3VAL.........evaluates the interpolating  function  determined  by
!                DB3INK or one of its derivatives.
!
!
! Reference
! Carl de Boor, A Practical Guide  to  Splines,  Springer-Verlag,  New
!      York, 1978.
!
!===============================================================================
MODULE MOOSE_CMLIB_DTENSBS
  USE MOOSE_CMLIB_XERROR
  USE MOOSE_CMLIB_DBSPLIN

  CONTAINS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DB2INK(X,NX,Y,NY,FCN,LDF,KX,KY,TX,TY,BCOEF,WORK,IFLAG)
  !***BEGIN PROLOGUE  DB2INK
  !***DATE WRITTEN   25 MAY 1982
  !***REVISION DATE  25 MAY 1982
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E1A
  !***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
  !         PIECEWISE POLYNOMIALS
  !***AUTHOR  BOISVERT, RONALD, NBS
  !         SCIENTIFIC COMPUTING DIVISION
  !         NATIONAL BUREAU OF STANDARDS
  !         WASHINGTON, DC 20234
  !***PURPOSE  DOUBLE PRECISION VERSION OF B2INK.
  !        DB2INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
  !        INTERPOLATES TWO-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
  !        THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
  !        (OPTIONALLY) THE KNOT SEQUENCE.
  !***DESCRIPTION
  !
  !   DB2INK determines the parameters of a  function  that  interpolates
  !   the two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX
  !   and j=1,..,NY. The interpolating function and its  derivatives  may
  !   subsequently be evaluated by the function DB2VAL.
  !
  !   The interpolating  function  is  a  piecewise  polynomial  function
  !   represented as a tensor product of one-dimensional  B-splines.  The
  !   form of this function is
  !
  !                      NX   NY
  !          S(x,y)  =  SUM  SUM  a   U (x) V (y)
  !                     i=1  j=1   ij  i     j
  !
  !   where the functions U(i)  and  V(j)  are  one-dimensional  B-spline
  !   basis functions. The coefficients a(i,j) are chosen so that
  !
  !     S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY
  !
  !   Note that  for  each  fixed  value  of  y  S(x,y)  is  a  piecewise
  !   polynomial function of x alone, and for each fixed value of x  S(x,
  !   y) is a piecewise polynomial function of y alone. In one  dimension
  !   a piecewise polynomial may  be  created  by  partitioning  a  given
  !   interval into subintervals and defining a distinct polynomial piece
  !   on each one. The points where adjacent subintervals meet are called
  !   knots. Each of the functions U(i) and V(j)  above  is  a  piecewise
  !   polynomial.
  !
  !   Users of DB2INK choose  the  order  (degree+1)  of  the  polynomial
  !   pieces used to define the piecewise polynomial in each of the x and
  !   y directions (KX and KY). Users also  may  define  their  own  knot
  !   sequence in x and y separately (TX and TY).  If  IFLAG=0,  however,
  !   DB2INK will choose sequences of knots that result  in  a  piecewise
  !   polynomial interpolant with KX-2 continuous partial derivatives  in
  !   x and KY-2 continuous partial derivatives in y. (KX knots are taken
  !   near each endpoint in the x direction,  not-a-knot  end  conditions
  !   are used, and the remaining knots are placed at data points  if  KX
  !   is even or at midpoints between data points if KX  is  odd.  The  y
  !   direction is treated similarly.)
  !
  !   After a call to DB2INK, all information  necessary  to  define  the
  !   interpolating function are contained in the parameters NX, NY,  KX,
  !   KY, TX, TY, and BCOEF. These quantities should not be altered until
  !   after the last call of the evaluation routine DB2VAL.
  !
  !
  !   I N P U T
  !   ---------
  !
  !   X       Double precision 1D array (size NX)
  !       Array of x abcissae. Must be strictly increasing.
  !
  !   NX      Integer scalar (.GE. 3)
  !       Number of x abcissae.
  !
  !   Y       Double precision 1D array (size NY)
  !       Array of y abcissae. Must be strictly increasing.
  !
  !   NY      Integer scalar (.GE. 3)
  !       Number of y abcissae.
  !
  !   FCN     Double precision 2D array (size LDF by NY)
  !       Array of function values to interpolate. FCN(I,J) should
  !       contain the function value at the point (X(I),Y(J))
  !
  !   LDF     Integer scalar (.GE. NX)
  !       The actual leading dimension of FCN used in the calling
  !       calling program.
  !
  !   KX      Integer scalar (.GE. 2, .LT. NX)
  !       The order of spline pieces in x.
  !       (Order = polynomial degree + 1)
  !
  !   KY      Integer scalar (.GE. 2, .LT. NY)
  !       The order of spline pieces in y.
  !       (Order = polynomial degree + 1)
  !
  !
  !   I N P U T   O R   O U T P U T
  !   -----------------------------
  !
  !   TX      Double precision 1D array (size NX+KX)
  !       The knots in the x direction for the spline interpolant.
  !       If IFLAG=0 these are chosen by DB2INK.
  !       If IFLAG=1 these are specified by the user.
  !                  (Must be non-decreasing.)
  !
  !   TY      Double precision 1D array (size NY+KY)
  !       The knots in the y direction for the spline interpolant.
  !       If IFLAG=0 these are chosen by DB2INK.
  !       If IFLAG=1 these are specified by the user.
  !                  (Must be non-decreasing.)
  !
  !
  !   O U T P U T
  !   -----------
  !
  !   BCOEF   Double precision 2D array (size NX by NY)
  !       Array of coefficients of the B-spline interpolant.
  !       This may be the same array as FCN.
  !
  !
  !   M I S C E L L A N E O U S
  !   -------------------------
  !
  !   WORK    Double precision 1D array (size NX*NY + max( 2*KX*(NX+1),
  !                                         2*KY*(NY+1) ))
  !       Array of working storage.
  !
  !   IFLAG   Integer scalar.
  !       On input:  0 == knot sequence chosen by DB2INK
  !                  1 == knot sequence chosen by user.
  !       On output: 1 == successful execution
  !                  2 == IFLAG out of range
  !                  3 == NX out of range
  !                  4 == KX out of range
  !                  5 == X not strictly increasing
  !                      6 == TX not non-decreasing
  !                  7 == NY out of range
  !                  8 == KY out of range
  !                  9 == Y not strictly increasing
  !                     10 == TY not non-decreasing
  !
  !***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
  !             SPRINGER-VERLAG, NEW YORK, 1978.
  !           CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
  !             PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !             VOL. 5 (1979), PP. 173-182.
  !***ROUTINES CALLED  DBTPCF,DBKNOT
  !***END PROLOGUE  DB2INK
  !
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        NX, NY, LDF, KX, KY, IFLAG
  DOUBLE PRECISION :: &
        X(NX), Y(NY), FCN(LDF,NY), TX(*), TY(*), BCOEF(NX,NY), &
        WORK(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        I, IW, NPK
  !
  !  -----------------------
  !  CHECK VALIDITY OF INPUT
  !  -----------------------
  !
  !***FIRST EXECUTABLE STATEMENT
  IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
  IF (NX .LT. 3)  GO TO 930
  IF (NY .LT. 3)  GO TO 970
  IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
  IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
  DO I=2,NX
     IF (X(I) .LE. X(I-1))  GO TO 950
  END DO
  DO I=2,NY
     IF (Y(I) .LE. Y(I-1))  GO TO 990
  END DO
  IF (IFLAG .EQ. 0)  GO TO 50
     NPK = NX + KX
     DO I=2,NPK
        IF (TX(I) .LT. TX(I-1))  GO TO 960
     END DO
     NPK = NY + KY
     DO I=2,NPK
        IF (TY(I) .LT. TY(I-1))  GO TO 1000
     END DO
   50   CONTINUE
  !
  !  ------------
  !  CHOOSE KNOTS
  !  ------------
  !
  IF (IFLAG .NE. 0)  GO TO 100
     CALL DBKNOT(X,NX,KX,TX)
     CALL DBKNOT(Y,NY,KY,TY)
  100   CONTINUE
  !
  !  -------------------------------
  !  CONSTRUCT B-SPLINE COEFFICIENTS
  !  -------------------------------
  !
  IFLAG = 1
  IW = NX*NY + 1
  CALL DBTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,WORK(IW))
  CALL DBTPCF(Y,NY,WORK,NY,NX,TY,KY,BCOEF,WORK(IW))
  GO TO 9999
  !
  !  -----
  !  EXITS
  !  -----
  !
  920   CONTINUE
  CALL XERRWV('DB2INK -  IFLAG=I1 IS OUT OF RANGE.', &
        36,2,1,1,IFLAG,I2,0,R1,R2)
  IFLAG = 2
  GO TO 9999
  !
  930   CONTINUE
  IFLAG = 3
  CALL XERRWV('DB2INK -  NX=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,NX,I2,0,R1,R2)
  GO TO 9999
  !
  940   CONTINUE
  IFLAG = 4
  CALL XERRWV('DB2INK -  KX=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,KX,I2,0,R1,R2)
  GO TO 9999
  !
  950   CONTINUE
  IFLAG = 5
  CALL XERRWV('DB2INK -  X ARRAY MUST BE STRICTLY INCREASING.', &
        46,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
  960   CONTINUE
  IFLAG = 6
  CALL XERRWV('DB2INK -  TX ARRAY MUST BE NON-DECREASING.', &
        42,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
  970   CONTINUE
  IFLAG = 7
  CALL XERRWV('DB2INK -  NY=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,NY,I2,0,R1,R2)
  GO TO 9999
  !
  980   CONTINUE
  IFLAG = 8
  CALL XERRWV('DB2INK -  KY=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,KY,I2,0,R1,R2)
  GO TO 9999
  !
  990   CONTINUE
  IFLAG = 9
  CALL XERRWV('DB2INK -  Y ARRAY MUST BE STRICTLY INCREASING.', &
        46,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 1000   CONTINUE
  IFLAG = 10
  CALL XERRWV('DB2INK -  TY ARRAY MUST BE NON-DECREASING.', &
        42,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 9999   CONTINUE
  RETURN
  END SUBROUTINE DB2INK
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION DB2VAL(XVAL,YVAL,IDX,IDY,TX,TY,NX,NY, &
        KX,KY,ILOY,INBVX,BCOEF,WORK)
  !***BEGIN PROLOGUE  DB2VAL
  !***DATE WRITTEN   25 MAY 1982
  !***REVISION DATE  25 MAY 1982
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E1A
  !***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
  !         PIECEWISE POLYNOMIALS
  !***AUTHOR  BOISVERT, RONALD, NBS
  !         SCIENTIFIC COMPUTING DIVISION
  !         NATIONAL BUREAU OF STANDARDS
  !         WASHINGTON, DC 20234
  !***PURPOSE  DB2VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
  !        FUNCTION CONSTRUCTED BY THE ROUTINE DB2INK OR ONE OF ITS
  !        PARTIAL DERIVATIVES.
  !        DOUBLE PRECISION VERSION OF B2VAL.
  !***DESCRIPTION
  !
  !   DB2VAL  evaluates   the   tensor   product   piecewise   polynomial
  !   interpolant constructed  by  the  routine  DB2INK  or  one  of  its
  !   derivatives at the point (XVAL,YVAL). To evaluate  the  interpolant
  !   itself, set IDX=IDY=0, to evaluate the first partial  with  respect
  !   to x, set IDX=1,IDY=0, and so on.
  !
  !   DB2VAL returns 0.0E0 if (XVAL,YVAL) is out of range. That is, if
  !        XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
  !        YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY)
  !   If the knots TX  and  TY  were  chosen  by  DB2INK,  then  this  is
  !   equivalent to
  !        XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
  !        YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY
  !   where EPSX = 0.1*(X(NX)-X(NX-1)) and EPSY = 0.1*(Y(NY)-Y(NY-1)).
  !
  !   The input quantities TX, TY, NX, NY, KX, KY, and  BCOEF  should  be
  !   unchanged since the last call of DB2INK.
  !
  !
  !   I N P U T
  !   ---------
  !
  !   XVAL    Double precision scalar
  !       X coordinate of evaluation point.
  !
  !   YVAL    Double precision scalar
  !       Y coordinate of evaluation point.
  !
  !   IDX     Integer scalar
  !       X derivative of piecewise polynomial to evaluate.
  !
  !   IDY     Integer scalar
  !       Y derivative of piecewise polynomial to evaluate.
  !
  !   TX      Double precision 1D array (size NX+KX)
  !       Sequence of knots defining the piecewise polynomial in
  !       the x direction.  (Same as in last call to DB2INK.)
  !
  !   TY      Double precision 1D array (size NY+KY)
  !       Sequence of knots defining the piecewise polynomial in
  !       the y direction.  (Same as in last call to DB2INK.)
  !
  !   NX      Integer scalar
  !       The number of interpolation points in x.
  !       (Same as in last call to DB2INK.)
  !
  !   NY      Integer scalar
  !       The number of interpolation points in y.
  !       (Same as in last call to DB2INK.)
  !
  !   KX      Integer scalar
  !       Order of polynomial pieces in x.
  !       (Same as in last call to DB2INK.)
  !
  !   KY      Integer scalar
  !       Order of polynomial pieces in y.
  !       (Same as in last call to DB2INK.)
  !
  !   BCOEF   Double precision 2D array (size NX by NY)
  !       The B-spline coefficients computed by DB2INK.
  !
  !   WORK    Double precision 1D array (size 3*max(KX,KY) + KY)
  !       A working storage array.
  !
  !***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
  !             SPRINGER-VERLAG, NEW YORK, 1978.
  !***ROUTINES CALLED  DINTRV,DBVALU
  !***END PROLOGUE  DB2VAL
  !
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !
  !   MODIFICATION
  !   ------------
  !
  !   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
  !
  !   R.F. BOISVERT, NIST
  !   22 FEB 00
  !
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        IDX, IDY, NX, NY, KX, KY, ILOY, INBVX
  DOUBLE PRECISION :: &
        XVAL, YVAL, TX(*), TY(*), BCOEF(NX,NY), WORK(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        INBV, K, LEFTY, MFLAG, KCOL, IW
  !
  ! DATA ILOY /1/,  INBVX /1/
  ! SAVE ILOY    ,  INBVX
  !
  !
  !***FIRST EXECUTABLE STATEMENT
  DB2VAL = 0.0D0
  !  NEXT STATEMENT - RFB MOD
  IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR. &
        YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY))      GO TO 100
  CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
  IF (MFLAG .NE. 0)  GO TO 100
     IW = KY + 1
     KCOL = LEFTY - KY
     DO K=1,KY
        KCOL = KCOL + 1
        WORK(K) = DBVALU(TX,BCOEF(1,KCOL),NX,KX,IDX,XVAL,INBVX, &
              WORK(IW))
     END DO
     INBV = 1
     KCOL = LEFTY - KY + 1
     DB2VAL = DBVALU(TY(KCOL),WORK,KY,KY,IDY,YVAL,INBV,WORK(IW))
  100   CONTINUE
  RETURN
  END FUNCTION DB2VAL
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DB3INK(X,NX,Y,NY,Z,NZ,FCN,LDF1,LDF2,KX,KY,KZ,TX,TY,TZ, &
        BCOEF,WORK,IFLAG)
  !***BEGIN PROLOGUE  DB3INK
  !***DATE WRITTEN   25 MAY 1982
  !***REVISION DATE  25 MAY 1982
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E1A
  !***KEYWORDS  INTERPOLATION, THREE-DIMENSIONS, GRIDDED DATA, SPLINES,
  !         PIECEWISE POLYNOMIALS
  !***AUTHOR  BOISVERT, RONALD, NBS
  !         SCIENTIFIC COMPUTING DIVISION
  !         NATIONAL BUREAU OF STANDARDS
  !         WASHINGTON, DC 20234
  !***PURPOSE  DOUBLE PRECISION VERSION OF DB3INK
  !        DB3INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
  !        INTERPOLATES THREE-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
  !        THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
  !        (OPTIONALLY) THE KNOT SEQUENCE.
  !***DESCRIPTION
  !
  !   DB3INK determines the parameters of a  function  that  interpolates
  !   the three-dimensional gridded data (X(i),Y(j),Z(k),FCN(i,j,k))  for
  !   i=1,..,NX, j=1,..,NY, and k=1,..,NZ. The interpolating function and
  !   its derivatives may  subsequently  be  evaluated  by  the  function
  !   DB3VAL.
  !
  !   The interpolating  function  is  a  piecewise  polynomial  function
  !   represented as a tensor product of one-dimensional  B-splines.  The
  !   form of this function is
  !
  !                  NX   NY   NZ
  !    S(x,y,z)  =  SUM  SUM  SUM  a   U (x) V (y) W (z)
  !                 i=1  j=1  k=1   ij  i     j     k
  !
  !   where the functions U(i), V(j), and  W(k)  are  one-dimensional  B-
  !   spline basis functions. The coefficients a(i,j) are chosen so that
  !
  !   S(X(i),Y(j),Z(k)) = FCN(i,j,k)  for i=1,..,NX, j=1,..,NY, k=1,..,NZ
  !
  !   Note that for fixed values of y  and  z  S(x,y,z)  is  a  piecewise
  !   polynomial function of x alone, for fixed values of x and z  S(x,y,
  !   z) is a piecewise polynomial function of y  alone,  and  for  fixed
  !   values of x and y S(x,y,z)  is  a  function  of  z  alone.  In  one
  !   dimension a piecewise polynomial may be created by  partitioning  a
  !   given interval into subintervals and defining a distinct polynomial
  !   piece on each one. The points where adjacent subintervals meet  are
  !   called knots. Each of the functions U(i), V(j), and W(k) above is a
  !   piecewise polynomial.
  !
  !   Users of DB3INK choose  the  order  (degree+1)  of  the  polynomial
  !   pieces used to define the piecewise polynomial in each of the x, y,
  !   and z directions (KX, KY, and KZ). Users also may define their  own
  !   knot sequence in x, y, and z separately (TX, TY, and TZ). If IFLAG=
  !   0, however, DB3INK will choose sequences of knots that result in  a
  !   piecewise  polynomial  interpolant  with  KX-2  continuous  partial
  !   derivatives in x, KY-2 continuous partial derivatives in y, and KZ-
  !   2 continuous partial derivatives in z. (KX  knots  are  taken  near
  !   each endpoint in x, not-a-knot end conditions  are  used,  and  the
  !   remaining knots are placed at data points  if  KX  is  even  or  at
  !   midpoints between data points if KX is odd. The y and z  directions
  !   are treated similarly.)
  !
  !   After a call to DB3INK, all information  necessary  to  define  the
  !   interpolating function are contained in the parameters NX, NY,  NZ,
  !   KX, KY, KZ, TX, TY, TZ, and BCOEF. These quantities should  not  be
  !   altered until after the last call of the evaluation routine DB3VAL.
  !
  !
  !   I N P U T
  !   ---------
  !
  !   X       Double precision 1D array (size NX)
  !       Array of x abcissae. Must be strictly increasing.
  !
  !   NX      Integer scalar (.GE. 3)
  !       Number of x abcissae.
  !
  !   Y       Double precision 1D array (size NY)
  !       Array of y abcissae. Must be strictly increasing.
  !
  !   NY      Integer scalar (.GE. 3)
  !       Number of y abcissae.
  !
  !   Z       Double precision 1D array (size NZ)
  !       Array of z abcissae. Must be strictly increasing.
  !
  !   NZ      Integer scalar (.GE. 3)
  !       Number of z abcissae.
  !
  !   FCN     Double precision 3D array (size LDF1 by LDF2 by NY)
  !       Array of function values to interpolate. FCN(I,J,K) should
  !       contain the function value at the point (X(I),Y(J),Z(K))
  !
  !   LDF1    Integer scalar (.GE. NX)
  !       The actual first dimension of FCN used in the
  !       calling program.
  !
  !   LDF2    Integer scalar (.GE. NY)
  !       The actual second dimension of FCN used in the calling
  !       program.
  !
  !   KX      Integer scalar (.GE. 2, .LT. NX)
  !       The order of spline pieces in x.
  !       (Order = polynomial degree + 1)
  !
  !   KY      Integer scalar (.GE. 2, .LT. NY)
  !       The order of spline pieces in y.
  !       (Order = polynomial degree + 1)
  !
  !   KZ      Integer scalar (.GE. 2, .LT. NZ)
  !       The order of spline pieces in z.
  !       (Order = polynomial degree + 1)
  !
  !
  !   I N P U T   O R   O U T P U T
  !   -----------------------------
  !
  !   TX      Double precision 1D array (size NX+KX)
  !       The knots in the x direction for the spline interpolant.
  !       If IFLAG=0 these are chosen by DB3INK.
  !       If IFLAG=1 these are specified by the user.
  !                  (Must be non-decreasing.)
  !
  !   TY      Double precision 1D array (size NY+KY)
  !       The knots in the y direction for the spline interpolant.
  !       If IFLAG=0 these are chosen by DB3INK.
  !       If IFLAG=1 these are specified by the user.
  !                  (Must be non-decreasing.)
  !
  !   TZ      Double precision 1D array (size NZ+KZ)
  !       The knots in the z direction for the spline interpolant.
  !       If IFLAG=0 these are chosen by DB3INK.
  !       If IFLAG=1 these are specified by the user.
  !                  (Must be non-decreasing.)
  !
  !
  !   O U T P U T
  !   -----------
  !
  !   BCOEF   Double precision 3D array (size NX by NY by NZ)
  !       Array of coefficients of the B-spline interpolant.
  !       This may be the same array as FCN.
  !
  !
  !   M I S C E L L A N E O U S
  !   -------------------------
  !
  !   WORK    Double precision 1D array (size NX*NY*NZ + max( 2*KX*(NX+1),
  !                         2*KY*(NY+1), 2*KZ*(NZ+1) )
  !       Array of working storage.
  !
  !   IFLAG   Integer scalar.
  !       On input:  0 == knot sequence chosen by B2INK
  !                  1 == knot sequence chosen by user.
  !       On output: 1 == successful execution
  !                  2 == IFLAG out of range
  !                  3 == NX out of range
  !                  4 == KX out of range
  !                  5 == X not strictly increasing
  !                      6 == TX not non-decreasing
  !                  7 == NY out of range
  !                  8 == KY out of range
  !                  9 == Y not strictly increasing
  !                     10 == TY not non-decreasing
  !                 11 == NZ out of range
  !                 12 == KZ out of range
  !                 13 == Z not strictly increasing
  !                     14 == TY not non-decreasing
  !
  !***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
  !             SPRINGER-VERLAG, NEW YORK, 1978.
  !           CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
  !             PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !             VOL. 5 (1979), PP. 173-182.
  !***ROUTINES CALLED  DBTPCF,DBKNOT
  !***END PROLOGUE  DB3INK
  !
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        NX, NY, NZ, LDF1, LDF2, KX, KY, KZ, IFLAG
  DOUBLE PRECISION :: &
        X(NX), Y(NY), Z(NZ), FCN(LDF1,LDF2,NZ), TX(*), TY(*), TZ(*), &
        BCOEF(NX,NY,NZ), WORK(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        I, J, LOC, IW, NPK
  !
  !  -----------------------
  !  CHECK VALIDITY OF INPUT
  !  -----------------------
  !
  !***FIRST EXECUTABLE STATEMENT
  IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
  IF (NX .LT. 3)  GO TO 930
  IF (NY .LT. 3)  GO TO 970
  IF (NZ .LT. 3)  GO TO 1010
  IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
  IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
  IF ((KZ .LT. 2) .OR. (KZ .GE. NZ))  GO TO 1020
  DO I=2,NX
     IF (X(I) .LE. X(I-1))  GO TO 950
  END DO
  DO I=2,NY
     IF (Y(I) .LE. Y(I-1))  GO TO 990
  END DO
  DO I=2,NZ
     IF (Z(I) .LE. Z(I-1))  GO TO 1030
  END DO
  IF (IFLAG .EQ. 0)  GO TO 70
     NPK = NX + KX
     DO I=2,NPK
        IF (TX(I) .LT. TX(I-1))  GO TO 960
     END DO
     NPK = NY + KY
     DO I=2,NPK
        IF (TY(I) .LT. TY(I-1))  GO TO 1000
     END DO
     NPK = NZ + KZ
     DO I=2,NPK
        IF (TZ(I) .LT. TZ(I-1))  GO TO 1040
     END DO
   70   CONTINUE
  !
  !  ------------
  !  CHOOSE KNOTS
  !  ------------
  !
  IF (IFLAG .NE. 0)  GO TO 100
     CALL DBKNOT(X,NX,KX,TX)
     CALL DBKNOT(Y,NY,KY,TY)
     CALL DBKNOT(Z,NZ,KZ,TZ)
  100   CONTINUE
  !
  !  -------------------------------
  !  CONSTRUCT B-SPLINE COEFFICIENTS
  !  -------------------------------
  !
  IFLAG = 1
  IW = NX*NY*NZ + 1
  !
  ! COPY FCN TO WORK IN PACKED FOR DBTPCF
  LOC = 0
  DO K=1,NZ
     DO J=1,NY
        DO I=1,NX
           LOC = LOC + 1
           WORK(LOC) = FCN(I,J,K)
        END DO
     END DO
  END DO
  !
  CALL DBTPCF(X,NX,WORK,NX,NY*NZ,TX,KX,BCOEF,WORK(IW))
  CALL DBTPCF(Y,NY,BCOEF,NY,NX*NZ,TY,KY,WORK,WORK(IW))
  CALL DBTPCF(Z,NZ,WORK,NZ,NX*NY,TZ,KZ,BCOEF,WORK(IW))
  GO TO 9999
  !
  !  -----
  !  EXITS
  !  -----
  !
  920   CONTINUE
  CALL XERRWV('DB3INK -  IFLAG=I1 IS OUT OF RANGE.', &
        35,2,1,1,IFLAG,I2,0,R1,R2)
  IFLAG = 2
  GO TO 9999
  !
  930   CONTINUE
  IFLAG = 3
  CALL XERRWV('DB3INK -  NX=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,NX,I2,0,R1,R2)
  GO TO 9999
  !
  940   CONTINUE
  IFLAG = 4
  CALL XERRWV('DB3INK -  KX=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,KX,I2,0,R1,R2)
  GO TO 9999
  !
  950   CONTINUE
  IFLAG = 5
  CALL XERRWV('DB3INK -  X ARRAY MUST BE STRICTLY INCREASING.', &
        46,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
  960   CONTINUE
  IFLAG = 6
  CALL XERRWV('DB3INK -  TX ARRAY MUST BE NON-DECREASING.', &
        42,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
  970   CONTINUE
  IFLAG = 7
  CALL XERRWV('DB3INK -  NY=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,NY,I2,0,R1,R2)
  GO TO 9999
  !
  980   CONTINUE
  IFLAG = 8
  CALL XERRWV('DB3INK -  KY=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,KY,I2,0,R1,R2)
  GO TO 9999
  !
  990   CONTINUE
  IFLAG = 9
  CALL XERRWV('DB3INK -  Y ARRAY MUST BE STRICTLY INCREASING.', &
        46,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 1000   CONTINUE
  IFLAG = 10
  CALL XERRWV('DB3INK -  TY ARRAY MUST BE NON-DECREASING.', &
        42,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 1010   CONTINUE
  IFLAG = 11
  CALL XERRWV('DB3INK -  NZ=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,NZ,I2,0,R1,R2)
  GO TO 9999
  !
 1020   CONTINUE
  IFLAG = 12
  CALL XERRWV('DB3INK -  KZ=I1 IS OUT OF RANGE.', &
        32,IFLAG,1,1,KZ,I2,0,R1,R2)
  GO TO 9999
  !
 1030   CONTINUE
  IFLAG = 13
  CALL XERRWV('DB3INK -  Z ARRAY MUST BE STRICTLY INCREASING.', &
        46,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 1040   CONTINUE
  IFLAG = 14
  CALL XERRWV('DB3INK -  TZ ARRAY MUST BE NON-DECREASING.', &
        42,IFLAG,1,0,I1,I2,0,R1,R2)
  GO TO 9999
  !
 9999   CONTINUE
  RETURN
  END SUBROUTINE DB3INK
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION DB3VAL(XVAL,YVAL,ZVAL,IDX,IDY,IDZ, &
        TX,TY,TZ,NX,NY,NZ,KX,KY,KZ,ILOY,ILOZ,INBVX,BCOEF,WORK)
  !***BEGIN PROLOGUE  DB3VAL
  !***DATE WRITTEN   25 MAY 1982
  !***REVISION DATE  25 MAY 1982
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  E1A
  !***KEYWORDS  INTERPOLATION, THREE-DIMENSIONS, GRIDDED DATA, SPLINES,
  !         PIECEWISE POLYNOMIALS
  !***AUTHOR  BOISVERT, RONALD, NBS
  !         SCIENTIFIC COMPUTING DIVISION
  !         NATIONAL BUREAU OF STANDARDS
  !         WASHINGTON, DC 20234
  !***PURPOSE  DB3VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
  !        FUNCTION CONSTRUCTED BY THE ROUTINE B3INK OR ONE OF ITS
  !        PARTIAL DERIVATIVES.
  !        DOUBLE PRECISION VERSION OF B3VAL.
  !***DESCRIPTION
  !
  !   DB3VAL  evaluates   the   tensor   product   piecewise   polynomial
  !   interpolant constructed  by  the  routine  DB3INK  or  one  of  its
  !   derivatives  at  the  point  (XVAL,YVAL,ZVAL).  To   evaluate   the
  !   interpolant  itself,  set  IDX=IDY=IDZ=0,  to  evaluate  the  first
  !   partial with respect to x, set IDX=1,IDY=IDZ=0, and so on.
  !
  !   DB3VAL returns 0.0D0 if (XVAL,YVAL,ZVAL) is out of range. That is,
  !        XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
  !        YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
  !        ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)
  !   If the knots TX, TY, and TZ were chosen by  DB3INK,  then  this  is
  !   equivalent to
  !        XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
  !        YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY .OR.
  !        ZVAL.LT.Z(1) .OR. ZVAL.GT.Z(NZ)+EPSZ
  !   where EPSX = 0.1*(X(NX)-X(NX-1)), EPSY =  0.1*(Y(NY)-Y(NY-1)),  and
  !   EPSZ = 0.1*(Z(NZ)-Z(NZ-1)).
  !
  !   The input quantities TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and  BCOEF
  !   should remain unchanged since the last call of DB3INK.
  !
  !
  !   I N P U T
  !   ---------
  !
  !   XVAL    Double precision scalar
  !       X coordinate of evaluation point.
  !
  !   YVAL    Double precision scalar
  !       Y coordinate of evaluation point.
  !
  !   ZVAL    Double precision scalar
  !       Z coordinate of evaluation point.
  !
  !   IDX     Integer scalar
  !       X derivative of piecewise polynomial to evaluate.
  !
  !   IDY     Integer scalar
  !       Y derivative of piecewise polynomial to evaluate.
  !
  !   IDZ     Integer scalar
  !       Z derivative of piecewise polynomial to evaluate.
  !
  !   TX      Double precision 1D array (size NX+KX)
  !       Sequence of knots defining the piecewise polynomial in
  !       the x direction.  (Same as in last call to DB3INK.)
  !
  !   TY      Double precision 1D array (size NY+KY)
  !       Sequence of knots defining the piecewise polynomial in
  !       the y direction.  (Same as in last call to DB3INK.)
  !
  !   TZ      Double precision 1D array (size NZ+KZ)
  !       Sequence of knots defining the piecewise polynomial in
  !       the z direction.  (Same as in last call to DB3INK.)
  !
  !   NX      Integer scalar
  !       The number of interpolation points in x.
  !       (Same as in last call to DB3INK.)
  !
  !   NY      Integer scalar
  !       The number of interpolation points in y.
  !       (Same as in last call to DB3INK.)
  !
  !   NZ      Integer scalar
  !       The number of interpolation points in z.
  !       (Same as in last call to DB3INK.)
  !
  !   KX      Integer scalar
  !       Order of polynomial pieces in x.
  !       (Same as in last call to DB3INK.)
  !
  !   KY      Integer scalar
  !       Order of polynomial pieces in y.
  !       (Same as in last call to DB3INK.)
  !
  !   KZ      Integer scalar
  !       Order of polynomial pieces in z.
  !       (Same as in last call to DB3INK.)
  !
  !   BCOEF   Double precision 2D array (size NX by NY by NZ)
  !       The B-spline coefficients computed by DB3INK.
  !
  !   WORK    Double precision 1D array (size KY*KZ+3*max(KX,KY,KZ)+KZ)
  !       A working storage array.
  !
  !***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
  !             SPRINGER-VERLAG, NEW YORK, 1978.
  !***ROUTINES CALLED  DINTRV,DBVALU
  !***END PROLOGUE  DB3VAL
  !
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !
  !   MODIFICATION
  !   ------------
  !
  !   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
  !
  !   R.F. BOISVERT, NIST
  !   22 FEB 00
  !
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        IDX, IDY, IDZ, NX, NY, NZ, KX, KY, KZ, ILOY, ILOZ, INBVX
  DOUBLE PRECISION :: &
        XVAL, YVAL, ZVAL, TX(*), TY(*), TZ(*), BCOEF(NX,NY,NZ), &
        WORK(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        INBV1, INBV2, LEFTY, LEFTZ, MFLAG, &
        KCOLY, KCOLZ, IZ, IZM1, IW, I, J, K
  !
  ! DATA ILOY /1/,  ILOZ /1/,  INBVX /1/
  ! SAVE ILOY    ,  ILOZ    ,  INBVX
  !
  !
  !***FIRST EXECUTABLE STATEMENT
  DB3VAL = 0.0D0
  !  NEXT STATEMENT - RFB MOD
  IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR. &
        YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR. &
        ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)) GO TO 100
  CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
  IF (MFLAG .NE. 0)  GO TO 100
  CALL DINTRV(TZ,NZ+KZ,ZVAL,ILOZ,LEFTZ,MFLAG)
  IF (MFLAG .NE. 0)  GO TO 100
     IZ = 1 + KY*KZ
     IW = IZ + KZ
     KCOLZ = LEFTZ - KZ
     I = 0
     DO K=1,KZ
        KCOLZ = KCOLZ + 1
        KCOLY = LEFTY - KY
        DO J=1,KY
           I = I + 1
           KCOLY = KCOLY + 1
           WORK(I) = DBVALU(TX,BCOEF(1,KCOLY,KCOLZ),NX,KX,IDX,XVAL, &
                 INBVX,WORK(IW))
        END DO
     END DO
     INBV1 = 1
     IZM1 = IZ - 1
     KCOLY = LEFTY - KY + 1
     DO K=1,KZ
        I = (K-1)*KY + 1
        J = IZM1 + K
        WORK(J) = DBVALU(TY(KCOLY),WORK(I),KY,KY,IDY,YVAL, &
              INBV1,WORK(IW))
     END DO
     INBV2 = 1
     KCOLZ = LEFTZ - KZ + 1
     DB3VAL = DBVALU(TZ(KCOLZ),WORK(IZ),KZ,KZ,IDZ,ZVAL,INBV2, &
           WORK(IW))
  100   CONTINUE
  RETURN
  END FUNCTION DB3VAL
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBKNOT(X,N,K,T)
  !***BEGIN PROLOGUE  DBKNOT
  !***REFER TO  DB2INK,DB3INK
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  DBKNOT
  !
  !  --------------------------------------------------------------------
  !  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
  !  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
  !  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
  !  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
  !  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
  !  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
  !  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
  !  DOUBLE PRECISION VERSION OF BKNOT.
  !  --------------------------------------------------------------------
  !
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        N, K
  DOUBLE PRECISION :: &
        X(N), T(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        I, J, IPJ, NPJ, IP1
  DOUBLE PRECISION :: &
        RNOT
  !
  !
  !  ----------------------------
  !  PUT K KNOTS AT EACH ENDPOINT
  !  ----------------------------
  !
  ! (SHIFT RIGHT ENPOINTS SLIGHTLY -- SEE PG 350 OF REFERENCE)
  RNOT = X(N) + 0.10D0*( X(N)-X(N-1) )
  DO J=1,K
     T(J) = X(1)
     NPJ = N + J
     T(NPJ) = RNOT
  END DO
  !
  !  --------------------------
  !  DISTRIBUTE REMAINING KNOTS
  !  --------------------------
  !
  IF (MOD(K,2) .EQ. 1)  GO TO 150
  !
  ! CASE OF EVEN K --  KNOTS AT DATA POINTS
  !
  I = (K/2) - K
  JSTRT = K+1
  DO J=JSTRT,N
     IPJ = I + J
     T(J) = X(IPJ)
  END DO
  GO TO 200
  !
  ! CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
  !
  150   CONTINUE
  I = (K-1)/2 - K
  IP1 = I + 1
  JSTRT = K + 1
  DO J=JSTRT,N
     IPJ = I + J
     T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  END DO
  200   CONTINUE
  !
  RETURN
  END SUBROUTINE DBKNOT
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE DBTPCF(X,N,FCN,LDF,NF,T,K,BCOEF,WORK)
  !***BEGIN PROLOGUE  DBTPCF
  !***REFER TO  DB2INK,DB3INK
  !***ROUTINES CALLED  DBINTK,DBNSLV
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  DBTPCF
  !
  !  -----------------------------------------------------------------
  !  DBTPCF COMPUTES B-SPLINE INTERPOLATION COEFFICIENTS FOR NF SETS
  !  OF DATA STORED IN THE COLUMNS OF THE ARRAY FCN. THE B-SPLINE
  !  COEFFICIENTS ARE STORED IN THE ROWS OF BCOEF HOWEVER.
  !  EACH INTERPOLATION IS BASED ON THE N ABCISSA STORED IN THE
  !  ARRAY X, AND THE N+K KNOTS STORED IN THE ARRAY T. THE ORDER
  !  OF EACH INTERPOLATION IS K. THE WORK ARRAY MUST BE OF LENGTH
  !  AT LEAST 2*K*(N+1).
  !  DOUBLE PRECISION VERSION OF BTPCF.
  !  -----------------------------------------------------------------
  !
  !  ------------
  !  DECLARATIONS
  !  ------------
  !
  !  PARAMETERS
  !
  INTEGER :: &
        N, LDF, K
  DOUBLE PRECISION :: &
        X(N), FCN(LDF,NF), T(*), BCOEF(NF,N), WORK(*)
  !
  !  LOCAL VARIABLES
  !
  INTEGER :: &
        I, J, K1, K2, IQ, IW
  !
  !  ---------------------------------------------
  !  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
  !  ---------------------------------------------
  !
  !***FIRST EXECUTABLE STATEMENT
  IF (NF .LE. 0)  GO TO 500
  K1 = K - 1
  K2 = K1 + K
  IQ = 1 + N
  IW = IQ + K2*N+1
  !
  !  -----------------------------
  !  COMPUTE B-SPLINE COEFFICIENTS
  !  -----------------------------
  !
  !
  !   FIRST DATA SET
  !
  CALL DBINTK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
  DO I=1,N
     BCOEF(1,I) = WORK(I)
  END DO
  !
  !  ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
  !
  IF (NF .EQ. 1)  GO TO 500
  DO J=2,NF
     DO I=1,N
        WORK(I) = FCN(I,J)
     END DO
     CALL DBNSLV(WORK(IQ),K2,N,K1,K1,WORK)
     DO I=1,N
        BCOEF(J,I) = WORK(I)
     END DO
  END DO
  !
  !  ----
  !  EXIT
  !  ----
  !
  500   CONTINUE
  RETURN
  END SUBROUTINE DBTPCF
  !-----------------------------------------------------------------------------

END MODULE MOOSE_CMLIB_DTENSBS
