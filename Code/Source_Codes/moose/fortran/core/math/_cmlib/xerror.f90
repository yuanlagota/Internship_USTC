!===============================================================================
!
!   P A C K A G E      XERROR
!
!   (Version 1982 )
!
!
!   A package of Fortran subprograms  for  the  processing  of  error
!   messages.  It  is  used  by  many  of  the   prewritten   library
!   subroutines in CMLIB.
!
!       A typical user of these library programs may  never  need  to
!   call any of the subroutines or functions in  the  XERROR  package
!   itself. Rather, these routines  will  make  themselves  known  by
!   their actions. These actions include printouts, dumps, tables and
!   tracebacks. The subroutine user may on occasion  want  to  change
!   the volume of this output, alter the file it goes to,  etc.  This
!   requires callling one or more modules in the XERROR package.  The
!   effect of these calls may be immediate or they  may  be  delayed.
!   Delayed  effects  are  caused  by  setting  flags  and  typically
!   consists  of  increasing/decreasing   frequency   of   subsequent
!   printout.
!
!       From a user's point of view there are only  a  few  important
!   things to know about XERROR. When a library  routine  detects  an
!   error situation XERROR will  be  called.  This  call  includes  a
!   message which describes the difficulty, an  error  number  and  a
!   seriousness  level.  The  XERROR  package  (usually)  prints  the
!   message, records the error number in a table  for  later  summary
!   and  then  may  take  additional  action   depending   upon   the
!   seriousness level. By convention, the error message  contains  in
!   its first few characters the name  of  the  routine  causing  the
!   error.
!
!       The seriousness level can be FATAL, RECOVERABLE, or  WARNING.
!   A FATAL error causes a program abort and in this situation XERROR
!   never returns to the subroutine which called it.  Thus  the  math
!   library subroutine writer who calls XERROR  with  a  FATAL  error
!   does not need to deal with what to  do  next.  For  example,  not
!   dimensioning an array large enough is a typical  FATAL  error.  A
!   RECOVERABLE error is something that the user  might  be  able  to
!   correct. For example, a data value out of range. By declaring  an
!   error RECOVERABLE, the library subroutine writer must  anticipate
!   that XERROR will return and must code for  the  possibility  that
!   the computation should be able to proceed. Actually,  XERROR  may
!   or may not return depending upon the user's choice.  (More  about
!   this shortly.) A WARNING error is just like a  RECOVERABLE  error
!   except that only unusual action by the user will  prevent  XERROR
!   from returning to the library routine. In  practice  most  errors
!   are declared either FATAL or WARNING.
!
!       Two things happen when XERROR receives a message. The message
!   is written out to a file, a traceback of where the call was  made
!   from is initiated and ultimately an abort may occur. The  default
!   output for error messages is the user's terminal. This is  almost
!   always the correct place for these messages  but  sometimes  they
!   can be annoying because they interrupt the normal flow of output.
!   This is especially  true  of  WARNING  messages  which  generally
!   require no action on the user's  part.  A  user  can  change  the
!   destination of these messages by inserting a call to
!
!                            CALL XSETUN(N)
!
!   This  diverts  all  messages  to  logical  unit  N  and   remains
!   operational until changed again or until the program ends.
!
!       The traceback process can also be  annoying.  Its  effect  is
!   system dependent but commonly it will drop the user into  a  dump
!   routine. If you do not want a traceback call
!
!                            CALL XSETF(-1)
!
!   Again, this remains operational until reset to +1  or  until  the
!   program ends.
!
!       These calls are only  two  examples  of  a  great  many  user
!   callable routines in the XERROR package, but they  are  the  most
!   useful. The individual  routines  describe  their  function  more
!   exactly. They are
!
!      CALL XERROR to process an error message
!      CALL XERRWV to process an error message with numeric values
!      CALL XSETF  to set the control variable, KONTRL (default=2)
!      CALL XGETF  to get the current value of KONTRL
!      CALL XERMAX to set limit on the number of times to print
!                  a message
!      CALL XERDMP to print error summary and clear tables
!      NER= NUMXER to get most recent message number
!      CALL XERCLR to clear current message number
!      CALL XSETUA to set up to 5 output unit numbers
!      CALL XGETUA to get current output unit numbers
!      CALL XSETUN to set one output unit for all messages
!      CALL XGETUN to get the one output unit
!      CALL XERABT to terminate and print a traceback
!      CALL XERCTL to perform special processing
!
!       The XERROR package is described in the  report  "XERROR,  The
!   SLATEC Error-Handling Package", by  R.E. Jones and D.K.  Kahaner,
!   Sandia  National  Laboratories,  Albuquerque,  New  Mexico  87185
!   (Sandia Report SAND82-0800 UC-32).
!
!===============================================================================
MODULE MOOSE_CMLIB_XERROR
  USE MOOSE_CMLIB_MACHCON

  CONTAINS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE FDUMP
  !***BEGIN PROLOGUE  FDUMP
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  Z
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Symbolic dump (should be locally written).
  !***DESCRIPTION
     ! ***Note*** Machine Dependent Routine
     ! FDUMP is intended to be replaced by a locally written
     ! version which produces a symbolic dump.  Failing this,
     ! it should be replaced by a version which prints the
     ! subprogram nesting list.  Note that this dump must be
     ! printed on each of up to five files, as indicated by the
     ! XGETUA routine.  See XSETUA and XGETUA for details.
  !
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  ! Latest revision ---  23 May 1979
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  FDUMP
  !***FIRST EXECUTABLE STATEMENT  FDUMP
  RETURN
  END SUBROUTINE FDUMP
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
  !***BEGIN PROLOGUE  J4SAVE
  !***REFER TO  XERROR
  ! Abstract
  !    J4SAVE saves and recalls several global variables needed
  !    by the library error handling routines.
  !
  ! Description of Parameters
  !  --Input--
  !    IWHICH - Index of item desired.
  !            = 1 Refers to current error number.
  !            = 2 Refers to current error control flag.
  !             = 3 Refers to current unit number to which error
  !                messages are to be sent.  (0 means use standard.)
  !             = 4 Refers to the maximum number of times any
  !                 message is to be printed (as set by XERMAX).
  !             = 5 Refers to the total number of units to which
  !                 each error message is to be written.
  !             = 6 Refers to the 2nd unit for error messages
  !             = 7 Refers to the 3rd unit for error messages
  !             = 8 Refers to the 4th unit for error messages
  !             = 9 Refers to the 5th unit for error messages
  !    IVALUE - The value to be set for the IWHICH-th parameter,
  !             if ISET is .TRUE. .
  !    ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
  !             given the value, IVALUE.  If ISET=.FALSE., the
  !             IWHICH-th parameter will be unchanged, and IVALUE
  !             is a dummy parameter.
  !  --Output--
  !    The (old) value of the IWHICH-th parameter will be returned
  !    in the function value, J4SAVE.
  !
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !    Adapted from Bell Laboratories PORT Library Error Handler
  ! Latest revision ---  23 MAY 1979
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  J4SAVE
  LOGICAL :: ISET
  INTEGER :: IPARAM(9)
  SAVE IPARAM
  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
  !***FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = IPARAM(IWHICH)
  IF (ISET) IPARAM(IWHICH) = IVALUE
  RETURN
  END FUNCTION J4SAVE
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERABT(MESSG,NMESSG)
  !***BEGIN PROLOGUE  XERABT
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  R3C
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Aborts program execution and prints error message.
  !***DESCRIPTION
  ! Abstract
  !    ***Note*** machine dependent routine
  !    XERABT aborts the execution of the program.
  !    The error message causing the abort is given in the calling
  !    sequence, in case one needs it for printing on a dayfile,
  !    for example.
  !
  ! Description of Parameters
  !    MESSG and NMESSG are as in XERROR, except that NMESSG may
  !    be zero, in which case no message is being supplied.
  !
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  ! Latest revision ---  19 MAR 1980
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  XERABT
  CHARACTER(len=*) :: MESSG
  !***FIRST EXECUTABLE STATEMENT  XERABT
  STOP
  END SUBROUTINE XERABT
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
  !***BEGIN PROLOGUE  XERCTL
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  R3C
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Allows user control over handling of individual errors.
  !***DESCRIPTION
  ! Abstract
  !    Allows user control over handling of individual errors.
  !    Just after each message is recorded, but before it is
  !    processed any further (i.e., before it is printed or
  !    a decision to abort is made), a call is made to XERCTL.
  !    If the user has provided his own version of XERCTL, he
  !    can then override the value of KONTROL used in processing
  !    this message by redefining its value.
  !    KONTRL may be set to any value from -2 to 2.
  !    The meanings for KONTRL are the same as in XSETF, except
  !    that the value of KONTRL changes only for this message.
  !    If KONTRL is set to a value outside the range from -2 to 2,
  !    it will be moved back into that range.
  !
  ! Description of Parameters
  !
  !  --Input--
  !    MESSG1 - the first word (only) of the error message.
  !    NMESSG - same as in the call to XERROR or XERRWV.
  !    NERR   - same as in the call to XERROR or XERRWV.
  !    LEVEL  - same as in the call to XERROR or XERRWV.
  !    KONTRL - the current value of the control flag as set
  !             by a call to XSETF.
  !
  !  --Output--
  !    KONTRL - the new value of KONTRL.  If KONTRL is not
  !             defined, it will remain at its original value.
  !             This changed value of control affects only
  !             the current occurrence of the current message.
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  XERCTL
  CHARACTER(len=20) :: MESSG1
  !***FIRST EXECUTABLE STATEMENT  XERCTL
  RETURN
  END SUBROUTINE XERCTL
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERPRT(MESSG,NMESSG)
  !***BEGIN PROLOGUE  XERPRT
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  Z
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Prints error messages.
  !***DESCRIPTION
  ! Abstract
  !    Print the Hollerith message in MESSG, of length NMESSG,
  !    on each file indicated by XGETUA.
  ! Latest revision ---  19 MAR 1980
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
  !***END PROLOGUE  XERPRT
  INTEGER :: LUN(5)
  CHARACTER(len=*) :: MESSG
  ! OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
  !***FIRST EXECUTABLE STATEMENT  XERPRT
  CALL XGETUA(LUN,NUNIT)
  LENMES = LEN(MESSG)
  DO KUNIT=1,NUNIT
     IUNIT = LUN(KUNIT)
     IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
     DO ICHAR=1,LENMES,72
        LAST = MIN0(ICHAR+71 , LENMES)
        WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
     END DO
  END DO
  RETURN
  END SUBROUTINE XERPRT
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
  !***BEGIN PROLOGUE  XERROR
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  R3C
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Processes an error (diagnostic) message.
  !***DESCRIPTION
  ! Abstract
  !    XERROR processes a diagnostic message, in a manner
  !    determined by the value of LEVEL and the current value
  !    of the library error control flag, KONTRL.
  !    (See subroutine XSETF for details.)
  !
  ! Description of Parameters
  !  --Input--
  !    MESSG - the Hollerith message to be processed, containing
  !            no more than 72 characters.
  !    NMESSG- the actual number of characters in MESSG.
  !    NERR  - the error number associated with this message.
  !            NERR must not be zero.
  !    LEVEL - error category.
  !            =2 means this is an unconditionally fatal error.
  !            =1 means this is a recoverable error.  (I.e., it is
  !               non-fatal if XSETF has been appropriately called.)
  !            =0 means this is a warning message only.
  !            =-1 means this is a warning message which is to be
  !               printed at most once, regardless of how many
  !               times this call is executed.
  !
  ! Examples
  !    CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
  !    CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
  !                43,2,1)
  !    CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
  !    1ULLY COLLAPSED.',65,3,0)
  !    CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
  !
  ! Latest revision ---  19 MAR 1980
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  XERRWV
  !***END PROLOGUE  XERROR
  CHARACTER(len=*) :: MESSG
  !***FIRST EXECUTABLE STATEMENT  XERROR
  CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
  RETURN
  END SUBROUTINE XERROR
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
  !***BEGIN PROLOGUE  XERRWV
  !***DATE WRITTEN   800319   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  R3C
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Processes error message allowing 2 integer and two real
         ! values to be included in the message.
  !***DESCRIPTION
  ! Abstract
  !    XERRWV processes a diagnostic message, in a manner
  !    determined by the value of LEVEL and the current value
  !    of the library error control flag, KONTRL.
  !    (See subroutine XSETF for details.)
  !    In addition, up to two integer values and two real
  !    values may be printed along with the message.
  !
  ! Description of Parameters
  !  --Input--
  !    MESSG - the Hollerith message to be processed.
  !    NMESSG- the actual number of characters in MESSG.
  !    NERR  - the error number associated with this message.
  !            NERR must not be zero.
  !    LEVEL - error category.
  !            =2 means this is an unconditionally fatal error.
  !            =1 means this is a recoverable error.  (I.e., it is
  !               non-fatal if XSETF has been appropriately called.)
  !            =0 means this is a warning message only.
  !            =-1 means this is a warning message which is to be
  !               printed at most once, regardless of how many
  !               times this call is executed.
  !    NI    - number of integer values to be printed. (0 to 2)
  !    I1    - first integer value.
  !    I2    - second integer value.
  !    NR    - number of real values to be printed. (0 to 2)
  !    R1    - first real value.
  !    R2    - second real value.
  !
  ! Examples
  !    CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
  !    1   1,NUM,0,0,0.,0.)
  !    CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
  !    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
  !
  ! Latest revision ---  19 MAR 1980
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
  !                XGETUA
  !***END PROLOGUE  XERRWV
  CHARACTER(len=*) :: MESSG
  CHARACTER(len=20) :: LFIRST
  CHARACTER(len=37) :: FORM
  DIMENSION LUN(5)
  ! GET FLAGS
  !***FIRST EXECUTABLE STATEMENT  XERRWV
  LKNTRL = J4SAVE(2,0,.FALSE.)
  MAXMES = J4SAVE(4,0,.FALSE.)
  ! CHECK FOR VALID INPUT
  IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND. &
        (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
     IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
     CALL XERPRT('XERROR -- INVALID INPUT',23)
     IF (LKNTRL.GT.0) CALL FDUMP
     IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.', &
           29)
     IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
     CALL XERABT('XERROR -- INVALID INPUT',23)
     RETURN
   10   CONTINUE
  ! RECORD MESSAGE
  JUNK = J4SAVE(1,NERR,.TRUE.)
  CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
  ! LET USER OVERRIDE
  LFIRST = MESSG
  LMESSG = NMESSG
  LERR = NERR
  LLEVEL = LEVEL
  CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
  ! RESET TO ORIGINAL VALUES
  LMESSG = NMESSG
  LERR = NERR
  LLEVEL = LEVEL
  LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
  MKNTRL = IABS(LKNTRL)
  ! DECIDE WHETHER TO PRINT MESSAGE
  IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
  IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES))) &
        .OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES)) &
        .OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1)) &
        .OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
     IF (LKNTRL.LE.0) GO TO 20
        CALL XERPRT(' ',1)
        ! INTRODUCTION
        IF (LLEVEL.EQ.(-1)) CALL XERPRT &
              ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
        IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
        IF (LLEVEL.EQ.1) CALL XERPRT &
              ('RECOVERABLE ERROR IN...',23)
        IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20   CONTINUE
     ! MESSAGE
     CALL XERPRT(MESSG,LMESSG)
     CALL XGETUA(LUN,NUNIT)
     ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
     ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
     DO KUNIT=1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
        DO I=1,MIN(NI,2)
           WRITE (FORM,21) I,ISIZEI
   21       FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
           IF (I.EQ.1) WRITE (IUNIT,FORM) I1
           IF (I.EQ.2) WRITE (IUNIT,FORM) I2
        END DO
        DO I=1,MIN(NR,2)
           WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23       FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E', &
                  I2,'.',I2,')')
           IF (I.EQ.1) WRITE (IUNIT,FORM) R1
           IF (I.EQ.2) WRITE (IUNIT,FORM) R2
        END DO
        IF (LKNTRL.LE.0) GO TO 40
           ! ERROR NUMBER
           WRITE (IUNIT,30) LERR
   30       FORMAT (15H ERROR NUMBER =,I10)
   40    CONTINUE
     END DO
     ! TRACE-BACK
     IF (LKNTRL.GT.0) CALL FDUMP
  100   CONTINUE
  IFATAL = 0
  IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2))) &
        IFATAL = 1
  ! QUIT HERE IF MESSAGE IS NOT FATAL
  IF (IFATAL.LE.0) RETURN
  IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
     ! PRINT REASON FOR ABORT
     IF (LLEVEL.EQ.1) CALL XERPRT &
           ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
     IF (LLEVEL.EQ.2) CALL XERPRT &
           ('JOB ABORT DUE TO FATAL ERROR.',29)
     ! PRINT ERROR SUMMARY
     CALL XERSAV(' ',-1,0,0,KDUMMY)
  120   CONTINUE
  ! ABORT
  IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
  CALL XERABT(MESSG,LMESSG)
  RETURN
  END SUBROUTINE XERRWV
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
  !***BEGIN PROLOGUE  XERSAV
  !***DATE WRITTEN   800319   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  Z
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Records that an error occurred.
  !***DESCRIPTION
  ! Abstract
  !    Record that this error occurred.
  !
  ! Description of Parameters
  ! --Input--
  !   MESSG, NMESSG, NERR, LEVEL are as in XERROR,
  !   except that when NMESSG=0 the tables will be
  !   dumped and cleared, and when NMESSG is less than zero the
  !   tables will be dumped and not cleared.
  ! --Output--
  !   ICOUNT will be the number of times this message has
  !   been seen, or zero if the table has overflowed and
  !   does not contain this message specifically.
  !   When NMESSG=0, ICOUNT will not be altered.
  !
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  ! Latest revision ---  19 Mar 1980
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
  !***END PROLOGUE  XERSAV
  INTEGER :: LUN(5)
  CHARACTER(len=*) :: MESSG
  CHARACTER(len=20) :: MESTAB(10),MES
  DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
  SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
  ! NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
  ! ERROR TABLE INITIALLY
  DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5), &
        KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10) &
        /0,0,0,0,0,0,0,0,0,0/
  DATA KOUNTX/0/
  !***FIRST EXECUTABLE STATEMENT  XERSAV
  IF (NMESSG.GT.0) GO TO 80
  ! DUMP THE TABLE
     IF (KOUNT(1).EQ.0) RETURN
     ! PRINT TO EACH UNIT
     CALL XGETUA(LUN,NUNIT)
     DO KUNIT=1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
        ! PRINT TABLE HEADER
        WRITE (IUNIT,10)
   10    FORMAT (32H0          ERROR MESSAGE SUMMARY/ &
             51  H MESSAGE START             NERR     LEVEL     COUNT)
        ! PRINT BODY OF TABLE
        DO I=1,10
           IF (KOUNT(I).EQ.0) GO TO 30
           WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15       FORMAT (1X,A20,3I10)
        END DO
   30    CONTINUE
        ! PRINT NUMBER OF OTHER ERRORS
        IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40    FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
        WRITE (IUNIT,50)
   50    FORMAT (1X)
     END DO
     IF (NMESSG.LT.0) RETURN
     ! CLEAR THE ERROR TABLES
     DO I=1,10
       KOUNT(I) = 0
     END DO
     KOUNTX = 0
     RETURN
   80   CONTINUE
  ! PROCESS A MESSAGE...
  ! SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
  ! OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
  MES = MESSG
  DO I=1,10
     II = I
     IF (KOUNT(I).EQ.0) GO TO 110
     IF (MES.NE.MESTAB(I)) GO TO 90
     IF (NERR.NE.NERTAB(I)) GO TO 90
     IF (LEVEL.NE.LEVTAB(I)) GO TO 90
     GO TO 100
   90 CONTINUE
  END DO
  ! THREE POSSIBLE CASES...
  ! TABLE IS FULL
     KOUNTX = KOUNTX+1
     ICOUNT = 1
     RETURN
  ! MESSAGE FOUND IN TABLE
  100   KOUNT(II) = KOUNT(II) + 1
     ICOUNT = KOUNT(II)
     RETURN
  ! EMPTY SLOT FOUND FOR NEW MESSAGE
  110   MESTAB(II) = MES
     NERTAB(II) = NERR
     LEVTAB(II) = LEVEL
     KOUNT(II)  = 1
     ICOUNT = 1
     RETURN
  END SUBROUTINE XERSAV
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE XGETUA(IUNITA,N)
  !***BEGIN PROLOGUE  XGETUA
  !***DATE WRITTEN   790801   (YYMMDD)
  !***REVISION DATE  820801   (YYMMDD)
  !***CATEGORY NO.  R3C
  !***KEYWORDS  ERROR,XERROR PACKAGE
  !***AUTHOR  JONES, R. E., (SNLA)
  !***PURPOSE  Returns unit number(s) to which error messages are being
         ! sent.
  !***DESCRIPTION
  ! Abstract
  !    XGETUA may be called to determine the unit number or numbers
  !    to which error messages are being sent.
  !    These unit numbers may have been set by a call to XSETUN,
  !    or a call to XSETUA, or may be a default value.
  !
  ! Description of Parameters
  !  --Output--
  !    IUNIT - an array of one to five unit numbers, depending
  !            on the value of N.  A value of zero refers to the
  !            default unit, as defined by the I1MACH machine
  !            constant routine.  Only IUNIT(1),...,IUNIT(N) are
  !            defined by XGETUA.  The values of IUNIT(N+1),...,
  !            IUNIT(5) are not defined (for N .LT. 5) or altered
  !            in any way by XGETUA.
  !    N     - the number of units to which copies of the
  !            error messages are being sent.  N will be in the
  !            range from 1 to 5.
  !
  ! Latest revision ---  19 MAR 1980
  ! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
  !             HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
  !             1982.
  !***ROUTINES CALLED  J4SAVE
  !***END PROLOGUE  XGETUA
  DIMENSION IUNITA(5)
  !***FIRST EXECUTABLE STATEMENT  XGETUA
  N = J4SAVE(5,0,.FALSE.)
  DO I=1,N
     INDEX = I+4
     IF (I.EQ.1) INDEX = 3
     IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
  END DO
  RETURN
  END SUBROUTINE XGETUA
  !-----------------------------------------------------------------------------

END MODULE MOOSE_CMLIB_XERROR
