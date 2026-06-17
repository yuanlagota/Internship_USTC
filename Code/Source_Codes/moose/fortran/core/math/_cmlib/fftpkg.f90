!=============================================================================== 
!
!***BEGIN PROLOGUE  FFTDOC
!***DATE WRITTEN   780201   (YYMMDD)
!***REVISION DATE  830701   (YYMMDD)
!***CATEGORY NO.  J1
!***KEYWORDS  DOCUMENTATION,FAST FOURIER TRANSFORM,FFT
!***AUTHOR  SWATZTRAUBER, PAUL, (NCAR)
!***PURPOSE  Documentation for FFT package.
!***DESCRIPTION
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                       Version 3  June 1979
!
!          A Package of Fortran Subprograms for The Fast Fourier
!           Transform of Periodic and Other Symmetric Sequences
!                              By
!                       Paul N Swarztrauber
!
!       National Center For Atmospheric Research  Boulder,Colorado 8030
!        which is sponsored by the National Science Foundation
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     This package consists of programs which perform Fast Fourier
!     Transforms for both complex and real periodic sequences and
!     certain other symmetric sequences that are listed below.
!
!     1.   RFFTI     Initialize RFFTF and RFFTB
!     2.   RFFTF     Forward transform of a real periodic sequence
!     3.   RFFTB     Backward transform of a real coefficient array
!
!     4.   CFFTI     Initialize CFFTF and CFFTB
!     5.   CFFTF     Forward transform of a complex periodic sequence
!     6.   CFFTB     Unnormalized inverse of CFFTF
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  FFTDOC
!
!=============================================================================== 
MODULE MOOSE_CMLIB_FFTPKG
  IMPLICIT DOUBLE PRECISION (A-H, O-Z)

  CONTAINS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTB1(N,C,CH,WA,WFAC)
  !***BEGIN PROLOGUE  CFFTB1
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  PASSB,PASSB2,PASSB3,PASSB4,PASSB5
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  CFFTB1
  DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WFAC(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTB1
  NF = INT(WFAC(2))
  NA = 0
  L1 = 1
  IW = 1
  DO K1=1,NF
     IP = INT(WFAC(K1+2))
     L2 = IP*L1
     IDO = N/L2
     IDOT = IDO+IDO
     IDL1 = IDOT*L1
     IF (IP .NE. 4) GO TO 103
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IF (NA .NE. 0) GO TO 101
     CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     GO TO 102
  101   CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102   NA = 1-NA
     GO TO 115
  103   IF (IP .NE. 2) GO TO 106
     IF (NA .NE. 0) GO TO 104
     CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
     GO TO 105
  104   CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105   NA = 1-NA
     GO TO 115
  106   IF (IP .NE. 3) GO TO 109
     IX2 = IW+IDOT
     IF (NA .NE. 0) GO TO 107
     CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
     GO TO 108
  107   CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108   NA = 1-NA
     GO TO 115
  109   IF (IP .NE. 5) GO TO 112
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IX4 = IX3+IDOT
     IF (NA .NE. 0) GO TO 110
     CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 111
  110   CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111   NA = 1-NA
     GO TO 115
  112   IF (NA .NE. 0) GO TO 113
     CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     GO TO 114
  113   CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114   IF (NAC .NE. 0) NA = 1-NA
  115   L1 = L2
     IW = IW+(IP-1)*IDOT
  END DO
  IF (NA .EQ. 0) RETURN
  N2 = N+N
  DO I=1,N2
     C(I) = CH(I)
  END DO
  RETURN
  END SUBROUTINE CFFTB1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTB(N,C,WSAVE)
  !***BEGIN PROLOGUE  CFFTB
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  830401   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  J1A2
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Unnormalized inverse of CFFTF.
  !***DESCRIPTION
  !
  !  Subroutine CFFTB computes the backward complex discrete Fourier
  !  transform (the Fourier synthesis).  Equivalently, CFFTB computes
  !  a complex periodic sequence from its Fourier coefficients.
  !  The transform is defined below at output parameter C.
  !
  !  A call of CFFTF followed by a call of CFFTB will multiply the
  !  sequence by N.
  !
  !  The array WSAVE which is used by subroutine CFFTB must be
  !  initialized by calling subroutine CFFTI(N,WSAVE).
  !
  !  Input Parameters
  !
  !
  !  N      the length of the complex sequence C.  The method is
  !     more efficient when N is the product of small primes.
  !
  !  C      a complex array of length N which contains the sequence
  !
  !  WSAVE   a real work array which must be dimensioned at least 4*N+15
  !      in the program that calls CFFTB.  The WSAVE array must be
  !      initialized by calling subroutine CFFTI(N,WSAVE), and a
  !      different WSAVE array must be used for each different
  !      value of N.  This initialization does not have to be
  !      repeated so long as N remains unchanged.  Thus subsequent
  !      transforms can be obtained faster than the first.
  !      The same WSAVE array can be used by CFFTF and CFFTB.
  !
  !  Output Parameters
  !
  !  C      For J=1,...,N
  !
  !         C(J)=the sum from K=1,...,N of
  !
  !               C(K)*EXP(I*J*K*2*PI/N)
  !
  !                     where I=SQRT(-1)
  !
  !  WSAVE   contains initialization calculations which must not be
  !      destroyed between calls of subroutine CFFTF or CFFTB
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CFFTB1
  !***END PROLOGUE  CFFTB
  DIMENSION       C(*)       ,WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTB
  IF (N .EQ. 1) RETURN
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
  RETURN
  END SUBROUTINE CFFTB
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTF1(N,C,CH,WA,WFAC)
  !***BEGIN PROLOGUE  CFFTF1
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  PASSF,PASSF2,PASSF3,PASSF4,PASSF5
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  CFFTF1
  DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WFAC(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTF1
  NF = INT(WFAC(2))
  NA = 0
  L1 = 1
  IW = 1
  DO K1=1,NF
     IP = INT(WFAC(K1+2))
     L2 = IP*L1
     IDO = N/L2
     IDOT = IDO+IDO
     IDL1 = IDOT*L1
     IF (IP .NE. 4) GO TO 103
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IF (NA .NE. 0) GO TO 101
     CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     GO TO 102
  101   CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102   NA = 1-NA
     GO TO 115
  103   IF (IP .NE. 2) GO TO 106
     IF (NA .NE. 0) GO TO 104
     CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
     GO TO 105
  104   CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105   NA = 1-NA
     GO TO 115
  106   IF (IP .NE. 3) GO TO 109
     IX2 = IW+IDOT
     IF (NA .NE. 0) GO TO 107
     CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
     GO TO 108
  107   CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108   NA = 1-NA
     GO TO 115
  109   IF (IP .NE. 5) GO TO 112
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IX4 = IX3+IDOT
     IF (NA .NE. 0) GO TO 110
     CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 111
  110   CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111   NA = 1-NA
     GO TO 115
  112   IF (NA .NE. 0) GO TO 113
     CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     GO TO 114
  113   CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114   IF (NAC .NE. 0) NA = 1-NA
  115   L1 = L2
     IW = IW+(IP-1)*IDOT
  END DO
  IF (NA .EQ. 0) RETURN
  N2 = N+N
  DO I=1,N2
     C(I) = CH(I)
  END DO
  RETURN
  END SUBROUTINE CFFTF1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTF(N,C,WSAVE)
  !***BEGIN PROLOGUE  CFFTF
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  800626   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  J1A2
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Forward transform of a complex, periodic sequence.
  !***DESCRIPTION
  !
  !  Subroutine CFFTF computes the forward complex discrete Fourier
  !  transform (the Fourier analysis).  Equivalently, CFFTF computes
  !  the Fourier coefficients of a complex periodic sequence.
  !  The transform is defined below at output parameter C.
  !
  !  The transform is not normalized.  To obtain a normalized transform
  !  the output must be divided by N.  Otherwise a call of CFFTF
  !  followed by a call of CFFTB will multiply the sequence by N.
  !
  !  The array WSAVE which is used by subroutine CFFTF must be
  !  initialized by calling subroutine CFFTI(N,WSAVE).
  !
  !  Input Parameters
  !
  !
  !  N      the length of the complex sequence C.  The method is
  !     more efficient when N is the product of small primes.
  !
  !  C      a complex array of length N which contains the sequence
  !
  !  WSAVE   a real work array which must be dimensioned at least 4*N+15
  !      in the program that calls CFFTF.  The WSAVE array must be
  !      initialized by calling subroutine CFFTI(N,WSAVE), and a
  !      different WSAVE array must be used for each different
  !      value of N.  This initialization does not have to be
  !      repeated so long as N remains unchanged.  Thus subsequent
  !      transforms can be obtained faster than the first.
  !      The same WSAVE array can be used by CFFTF and CFFTB.
  !
  !  Output Parameters
  !
  !  C      for J=1,...,N
  !
  !         C(J)=the sum from K=1,...,N of
  !
  !               C(K)*EXP(-I*J*K*2*PI/N)
  !
  !                     where I=SQRT(-1)
  !
  !  WSAVE   contains initialization calculations which must not be
  !      destroyed between calls of subroutine CFFTF or CFFTB
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CFFTF1
  !***END PROLOGUE  CFFTF
  DIMENSION       C(*)       ,WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTF
  IF (N .EQ. 1) RETURN
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
  RETURN
  END SUBROUTINE CFFTF
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTI1(N,WA,WFAC)
  !***BEGIN PROLOGUE  CFFTI1
  !***REFER TO  CFFTI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  CFFTI1
  DIMENSION       WA(*)      ,WFAC(*)    ,NTRYH(4)
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
  !***FIRST EXECUTABLE STATEMENT  CFFTI1
  NL = N
  NF = 0
  J = 0
  101   J = J+1
  IF (J-4) 102,102,103
  102   NTRY = NTRYH(J)
  GO TO 104
  103   NTRY = NTRY+2
  104   NQ = NL/NTRY
  NR = NL-NTRY*NQ
  IF (NR) 101,105,101
  105   NF = NF+1
  WFAC(NF+2) = NTRY
  NL = NQ
  IF (NTRY .NE. 2) GO TO 107
  IF (NF .EQ. 1) GO TO 107
  DO I=2,NF
     IB = NF-I+2
     WFAC(IB+2) = WFAC(IB+1)
  END DO
  WFAC(3) = 2
  107   IF (NL .NE. 1) GO TO 104
  WFAC(1) = N
  WFAC(2) = NF
  TPI = 6.28318530717959D0
  ARGH = TPI/FLOAT(N)
  I = 2
  L1 = 1
  DO K1=1,NF
     IP = INT(WFAC(K1+2))
     LD = 0
     L2 = L1*IP
     IDO = N/L2
     IDOT = IDO+IDO+2
     IPM = IP-1
     DO J=1,IPM
        I1 = I
        WA(I-1) = 1.D0
        WA(I) = 0.D0
        LD = LD+L1
        FI = 0.
        ARGLD = FLOAT(LD)*ARGH
        DO II=4,IDOT,2
           I = I+2
           FI = FI+1.D0
           ARG = FI*ARGLD
           WA(I-1) = COS(ARG)
           WA(I) = SIN(ARG)
        END DO
        IF (IP .LE. 5) GO TO 109
        WA(I1-1) = WA(I-1)
        WA(I1) = WA(I)
  109  CONTINUE
     END DO
     L1 = L2
  END DO
  RETURN
  END SUBROUTINE CFFTI1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE CFFTI(N,WSAVE)
  !***BEGIN PROLOGUE  CFFTI
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  830401   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  J1A2
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Initialize for CFFTF and CFFTB.
  !***DESCRIPTION
  !
  !  Subroutine CFFTI initializes the array WSAVE which is used in
  !  both CFFTF and CFFTB.  The prime factorization of N together with
  !  a tabulation of the trigonometric functions are computed and
  !  stored in WSAVE.
  !
  !  Input Parameter
  !
  !  N       the length of the sequence to be transformed
  !
  !  Output Parameter
  !
  !  WSAVE   a work array which must be dimensioned at least 4*N+15.
  !      The same work array can be used for both CFFTF and CFFTB
  !      as long as N remains unchanged.  Different WSAVE arrays
  !      are required for different values of N.  The contents of
  !      WSAVE must not be changed between calls of CFFTF or CFFTB.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CFFTI1
  !***END PROLOGUE  CFFTI
  DIMENSION       WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTI
  IF (N .EQ. 1) RETURN
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
  RETURN
  END SUBROUTINE CFFTI
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSB2(IDO,L1,CC,CH,WA1)
  !***BEGIN PROLOGUE  PASSB2
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSB2
  DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
        WA1(*)
  !***FIRST EXECUTABLE STATEMENT  PASSB2
  IF (IDO .GT. 2) GO TO 102
  DO K=1,L1
     CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
     CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
     CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
     CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSB2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSB3(IDO,L1,CC,CH,WA1,WA2)
  !***BEGIN PROLOGUE  PASSB3
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSB3
  DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
        WA1(*)     ,WA2(*)
  DATA TAUR,TAUI /-.5,.866025403784439/
  !***FIRST EXECUTABLE STATEMENT  PASSB3
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TR2 = CC(1,2,K)+CC(1,3,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     TI2 = CC(2,2,K)+CC(2,3,K)
     CI2 = CC(2,1,K)+TAUR*TI2
     CH(2,K,1) = CC(2,1,K)+TI2
     CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
     CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
     CH(2,K,2) = CI2+CR3
     CH(2,K,3) = CI2-CR3
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSB3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSB4(IDO,L1,CC,CH,WA1,WA2,WA3)
  !***BEGIN PROLOGUE  PASSB4
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSB4
  DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)
  !***FIRST EXECUTABLE STATEMENT  PASSB4
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TI1 = CC(2,1,K)-CC(2,3,K)
     TI2 = CC(2,1,K)+CC(2,3,K)
     TR4 = CC(2,4,K)-CC(2,2,K)
     TI3 = CC(2,2,K)+CC(2,4,K)
     TR1 = CC(1,1,K)-CC(1,3,K)
     TR2 = CC(1,1,K)+CC(1,3,K)
     TI4 = CC(1,2,K)-CC(1,4,K)
     TR3 = CC(1,2,K)+CC(1,4,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,3) = TR2-TR3
     CH(2,K,1) = TI2+TI3
     CH(2,K,3) = TI2-TI3
     CH(1,K,2) = TR1+TR4
     CH(1,K,4) = TR1-TR4
     CH(2,K,2) = TI1+TI4
     CH(2,K,4) = TI1-TI4
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,4,K)-CC(I,2,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,4,K)-CC(I,2,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSB4
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  !***BEGIN PROLOGUE  PASSB5
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSB5
  DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
        -.809016994374947,.587785252292473/
  !***FIRST EXECUTABLE STATEMENT  PASSB5
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TI5 = CC(2,2,K)-CC(2,5,K)
     TI2 = CC(2,2,K)+CC(2,5,K)
     TI4 = CC(2,3,K)-CC(2,4,K)
     TI3 = CC(2,3,K)+CC(2,4,K)
     TR5 = CC(1,2,K)-CC(1,5,K)
     TR2 = CC(1,2,K)+CC(1,5,K)
     TR4 = CC(1,3,K)-CC(1,4,K)
     TR3 = CC(1,3,K)+CC(1,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CH(2,K,1) = CC(2,1,K)+TI2+TI3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
     CR5 = TI11*TR5+TI12*TR4
     CI5 = TI11*TI5+TI12*TI4
     CR4 = TI12*TR5-TI11*TR4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,5) = CR2+CI5
     CH(2,K,2) = CI2+CR5
     CH(2,K,3) = CI3+CR4
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(2,K,4) = CI3-CR4
     CH(2,K,5) = CI2-CR5
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSB5
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  !***BEGIN PROLOGUE  PASSB
  !***REFER TO  CFFTB
  !***ROUTINES CALLED  (NONE)
  !***END PROLOGUE  PASSB
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
        C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
        CH2(IDL1,IP)
  !***FIRST EXECUTABLE STATEMENT  PASSB
  IDOT = IDO/2
  NT = IP*IDL1
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
  !
  IF (IDO .LT. L1) GO TO 106
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=1,IDO
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
     END DO
  END DO
  DO K=1,L1
  !DIR$ IVDEP
     DO I=1,IDO
        CH(I,K,1) = CC(I,1,K)
     END DO
  END DO
  GO TO 112
  106   DO 109 J=2,IPPH
     JC = IPP2-J
     DO I=1,IDO
  !DIR$ IVDEP
        DO K=1,L1
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
     END DO
  109   CONTINUE
  DO I=1,IDO
  !DIR$ IVDEP
     DO K=1,L1
        CH(I,K,1) = CC(I,1,K)
     END DO
  END DO
  112   IDL = 2-IDO
  INC = 0
  DO L=2,IPPH
     LC = IPP2-L
     IDL = IDL+IDO
  !DIR$ IVDEP
     DO IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
        C2(IK,LC) = WA(IDL)*CH2(IK,IP)
     END DO
     IDLJ = IDL
     INC = INC+IDO
     DO J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
        WAR = WA(IDLJ-1)
        WAI = WA(IDLJ)
  !DIR$ IVDEP
        DO IK=1,IDL1
           C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
           C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
        END DO
     END DO
  END DO
  DO J=2,IPPH
  !DIR$ IVDEP
     DO IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
     END DO
  END DO
  DO J=2,IPPH
     JC = IPP2-J
  !DIR$ IVDEP
     DO IK=2,IDL1,2
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
     END DO
  END DO
  NAC = 1
  IF (IDO .EQ. 2) RETURN
  NAC = 0
  DO IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  END DO
  DO J=2,IP
  !DIR$ IVDEP
     DO K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
     END DO
  END DO
  IF (IDOT .GT. L1) GO TO 127
  IDIJ = 0
  DO J=2,IP
     IDIJ = IDIJ+2
     DO I=4,IDO,2
        IDIJ = IDIJ+2
  !DIR$ IVDEP
        DO K=1,L1
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  RETURN
  127   IDJ = 2-IDO
  DO J=2,IP
     IDJ = IDJ+IDO
     DO K=1,L1
        IDIJ = IDJ
  !DIR$ IVDEP
        DO I=4,IDO,2
           IDIJ = IDIJ+2
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  RETURN
  END SUBROUTINE PASSB
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSF2(IDO,L1,CC,CH,WA1)
  !***BEGIN PROLOGUE  PASSF2
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSF2
  DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
        WA1(*)
  !***FIRST EXECUTABLE STATEMENT  PASSF2
  IF (IDO .GT. 2) GO TO 102
  DO K=1,L1
     CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
     CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
     CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
     CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
  DO K=1,L1
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSF2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSF3(IDO,L1,CC,CH,WA1,WA2)
  !***BEGIN PROLOGUE  PASSF3
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSF3
  DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
        WA1(*)     ,WA2(*)
  DATA TAUR,TAUI /-.5,-.866025403784439/
  !***FIRST EXECUTABLE STATEMENT  PASSF3
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TR2 = CC(1,2,K)+CC(1,3,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     TI2 = CC(2,2,K)+CC(2,3,K)
     CI2 = CC(2,1,K)+TAUR*TI2
     CH(2,K,1) = CC(2,1,K)+TI2
     CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
     CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
     CH(2,K,2) = CI2+CR3
     CH(2,K,3) = CI2-CR3
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSF3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSF4(IDO,L1,CC,CH,WA1,WA2,WA3)
  !***BEGIN PROLOGUE  PASSF4
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSF4
  DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)
  !***FIRST EXECUTABLE STATEMENT  PASSF4
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TI1 = CC(2,1,K)-CC(2,3,K)
     TI2 = CC(2,1,K)+CC(2,3,K)
     TR4 = CC(2,2,K)-CC(2,4,K)
     TI3 = CC(2,2,K)+CC(2,4,K)
     TR1 = CC(1,1,K)-CC(1,3,K)
     TR2 = CC(1,1,K)+CC(1,3,K)
     TI4 = CC(1,4,K)-CC(1,2,K)
     TR3 = CC(1,2,K)+CC(1,4,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,3) = TR2-TR3
     CH(2,K,1) = TI2+TI3
     CH(2,K,3) = TI2-TI3
     CH(1,K,2) = TR1+TR4
     CH(1,K,4) = TR1-TR4
     CH(2,K,2) = TI1+TI4
     CH(2,K,4) = TI1-TI4
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,2,K)-CC(I,4,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,4,K)-CC(I-1,2,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,2,K)-CC(I,4,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,4,K)-CC(I-1,2,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSF4
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  !***BEGIN PROLOGUE  PASSF5
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSF5
  DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,-.951056516295154, &
        -.809016994374947,-.587785252292473/
  !***FIRST EXECUTABLE STATEMENT  PASSF5
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1
     TI5 = CC(2,2,K)-CC(2,5,K)
     TI2 = CC(2,2,K)+CC(2,5,K)
     TI4 = CC(2,3,K)-CC(2,4,K)
     TI3 = CC(2,3,K)+CC(2,4,K)
     TR5 = CC(1,2,K)-CC(1,5,K)
     TR2 = CC(1,2,K)+CC(1,5,K)
     TR4 = CC(1,3,K)-CC(1,4,K)
     TR3 = CC(1,3,K)+CC(1,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CH(2,K,1) = CC(2,1,K)+TI2+TI3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
     CR5 = TI11*TR5+TI12*TR4
     CI5 = TI11*TI5+TI12*TI4
     CR4 = TI12*TR5-TI11*TR4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,5) = CR2+CI5
     CH(2,K,2) = CI2+CR5
     CH(2,K,3) = CI3+CR4
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(2,K,4) = CI3-CR4
     CH(2,K,5) = CI2-CR5
  END DO
  RETURN
  102   IF(IDO/2.LT.L1) GO TO 105
  DO K=1,L1
  !DIR$ IVDEP
     DO I=2,IDO,2
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
     END DO
  END DO
  RETURN
  105   DO 107 I=2,IDO,2
  !DIR$ IVDEP
     DO K=1,L1
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
     END DO
  107   CONTINUE
  RETURN
  END SUBROUTINE PASSF5
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE PASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  !***BEGIN PROLOGUE  PASSF
  !***REFER TO  CFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  PASSF
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
        C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP), &
        CH2(IDL1,IP)
  !***FIRST EXECUTABLE STATEMENT  PASSF
  IDOT = IDO/2
  NT = IP*IDL1
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
  !
  IF (IDO .LT. L1) GO TO 106
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=1,IDO
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
     END DO
  END DO
  DO K=1,L1
  !DIR$ IVDEP
     DO I=1,IDO
        CH(I,K,1) = CC(I,1,K)
     END DO
  END DO
  GO TO 112
  106   DO 109 J=2,IPPH
     JC = IPP2-J
     DO I=1,IDO
  !DIR$ IVDEP
        DO K=1,L1
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
        END DO
     END DO
  109   CONTINUE
  DO I=1,IDO
  !DIR$ IVDEP
     DO K=1,L1
        CH(I,K,1) = CC(I,1,K)
     END DO
  END DO
  112   IDL = 2-IDO
  INC = 0
  DO L=2,IPPH
     LC = IPP2-L
     IDL = IDL+IDO
  !DIR$ IVDEP
     DO IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
        C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
     END DO
     IDLJ = IDL
     INC = INC+IDO
     DO J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
        WAR = WA(IDLJ-1)
        WAI = WA(IDLJ)
  !DIR$ IVDEP
        DO IK=1,IDL1
           C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
           C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
        END DO
     END DO
  END DO
  DO J=2,IPPH
  !DIR$ IVDEP
     DO IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
     END DO
  END DO
  DO J=2,IPPH
     JC = IPP2-J
  !DIR$ IVDEP
     DO IK=2,IDL1,2
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
     END DO
  END DO
  NAC = 1
  IF (IDO .EQ. 2) RETURN
  NAC = 0
  !DIR$ IVDEP
  DO IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  END DO
  DO J=2,IP
  !DIR$ IVDEP
     DO K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
     END DO
  END DO
  IF (IDOT .GT. L1) GO TO 127
  IDIJ = 0
  DO J=2,IP
     IDIJ = IDIJ+2
     DO I=4,IDO,2
        IDIJ = IDIJ+2
  !DIR$ IVDEP
        DO K=1,L1
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  RETURN
  127   IDJ = 2-IDO
  DO J=2,IP
     IDJ = IDJ+IDO
     DO K=1,L1
        IDIJ = IDJ
  !DIR$ IVDEP
        DO I=4,IDO,2
           IDIJ = IDIJ+2
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  RETURN
  END SUBROUTINE PASSF
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADB2(IDO,L1,CC,CH,WA1)
  !***BEGIN PROLOGUE  RADB2
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADB2
  DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
        WA1(*)
  !***FIRST EXECUTABLE STATEMENT  RADB2
  DO K=1,L1
     CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
     CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  END DO
  IF (IDO-2) 107,105,102
  102   IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 108
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
        TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
        CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
        TI2 = CC(I,1,K)+CC(IC,2,K)
        CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
        CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
     END DO
  END DO
  GO TO 111
  108   DO 110 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
        TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
        CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
        TI2 = CC(I,1,K)+CC(IC,2,K)
        CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
        CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
     END DO
  110   CONTINUE
  111   IF (MOD(IDO,2) .EQ. 1) RETURN
  105   DO 106 K=1,L1
     CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
     CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106   CONTINUE
  107   RETURN
  END SUBROUTINE RADB2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADB3(IDO,L1,CC,CH,WA1,WA2)
  !***BEGIN PROLOGUE  RADB3
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADB3
  DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
        WA1(*)     ,WA2(*)
  DATA TAUR,TAUI /-.5,.866025403784439/
  !***FIRST EXECUTABLE STATEMENT  RADB3
  DO K=1,L1
     TR2 = CC(IDO,2,K)+CC(IDO,2,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
  END DO
  IF (IDO .EQ. 1) RETURN
  IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 104
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,3,K)-CC(IC,2,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
        CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
     END DO
  END DO
  RETURN
  104   DO 106 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,3,K)-CC(IC,2,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
        CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
     END DO
  106   CONTINUE
  RETURN
  END SUBROUTINE RADB3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADB4(IDO,L1,CC,CH,WA1,WA2,WA3)
  !***BEGIN PROLOGUE  RADB4
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADB4
  DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)
  DATA SQRT2 /1.414213562373095/
  !***FIRST EXECUTABLE STATEMENT  RADB4
  DO K=1,L1
     TR1 = CC(1,1,K)-CC(IDO,4,K)
     TR2 = CC(1,1,K)+CC(IDO,4,K)
     TR3 = CC(IDO,2,K)+CC(IDO,2,K)
     TR4 = CC(1,3,K)+CC(1,3,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,2) = TR1-TR4
     CH(1,K,3) = TR2-TR3
     CH(1,K,4) = TR1+TR4
  END DO
  IF (IDO-2) 107,105,102
  102   IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 108
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        TI1 = CC(I,1,K)+CC(IC,4,K)
        TI2 = CC(I,1,K)-CC(IC,4,K)
        TI3 = CC(I,3,K)-CC(IC,2,K)
        TR4 = CC(I,3,K)+CC(IC,2,K)
        TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
        TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
        TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1-TR4
        CR4 = TR1+TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
        CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
        CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
        CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
        CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
        CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
     END DO
  END DO
  GO TO 111
  108   DO 110 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        TI1 = CC(I,1,K)+CC(IC,4,K)
        TI2 = CC(I,1,K)-CC(IC,4,K)
        TI3 = CC(I,3,K)-CC(IC,2,K)
        TR4 = CC(I,3,K)+CC(IC,2,K)
        TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
        TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
        TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1-TR4
        CR4 = TR1+TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
        CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
        CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
        CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
        CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
        CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
     END DO
  110   CONTINUE
  111   IF (MOD(IDO,2) .EQ. 1) RETURN
  105   DO 106 K=1,L1
     TI1 = CC(1,2,K)+CC(1,4,K)
     TI2 = CC(1,4,K)-CC(1,2,K)
     TR1 = CC(IDO,1,K)-CC(IDO,3,K)
     TR2 = CC(IDO,1,K)+CC(IDO,3,K)
     CH(IDO,K,1) = TR2+TR2
     CH(IDO,K,2) = SQRT2*(TR1-TI1)
     CH(IDO,K,3) = TI2+TI2
     CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106   CONTINUE
  107   RETURN
  END SUBROUTINE RADB4
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  !***BEGIN PROLOGUE  RADB5
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADB5
  DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
        -.809016994374947,.587785252292473/
  !***FIRST EXECUTABLE STATEMENT  RADB5
  DO K=1,L1
     TI5 = CC(1,3,K)+CC(1,3,K)
     TI4 = CC(1,5,K)+CC(1,5,K)
     TR2 = CC(IDO,2,K)+CC(IDO,2,K)
     TR3 = CC(IDO,4,K)+CC(IDO,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI5 = TI11*TI5+TI12*TI4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(1,K,5) = CR2+CI5
  END DO
  IF (IDO .EQ. 1) RETURN
  IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 104
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        TI5 = CC(I,3,K)+CC(IC,2,K)
        TI2 = CC(I,3,K)-CC(IC,2,K)
        TI4 = CC(I,5,K)+CC(IC,4,K)
        TI3 = CC(I,5,K)-CC(IC,4,K)
        TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
        TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
        CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
        CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
        CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
        CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
     END DO
  END DO
  RETURN
  104   DO 106 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        TI5 = CC(I,3,K)+CC(IC,2,K)
        TI2 = CC(I,3,K)-CC(IC,2,K)
        TI4 = CC(I,5,K)+CC(IC,4,K)
        TI3 = CC(I,5,K)-CC(IC,4,K)
        TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
        TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
        CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
        CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
        CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
        CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
     END DO
  106   CONTINUE
  RETURN
  END SUBROUTINE RADB5
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADBG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  !***BEGIN PROLOGUE  RADBG
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADBG
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
        C1(IDO,L1,IP)          ,C2(IDL1,IP), &
        CH2(IDL1,IP)           ,WA(*)
  DATA TPI/6.28318530717959/
  !***FIRST EXECUTABLE STATEMENT  RADBG
  ARG = TPI/FLOAT(IP)
  DCP = COS(ARG)
  DSP = SIN(ARG)
  IDP2 = IDO+2
  NBD = (IDO-1)/2
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IF (IDO .LT. L1) GO TO 103
  DO K=1,L1
     DO I=1,IDO
        CH(I,K,1) = CC(I,1,K)
     END DO
  END DO
  GO TO 106
  103   DO 105 I=1,IDO
     DO K=1,L1
        CH(I,K,1) = CC(I,1,K)
     END DO
  105   CONTINUE
  106   DO 108 J=2,IPPH
     JC = IPP2-J
     J2 = J+J
     DO K=1,L1
        CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
        CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
     END DO
  108   CONTINUE
  IF (IDO .EQ. 1) GO TO 116
  IF (NBD .LT. L1) GO TO 112
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=3,IDO,2
           IC = IDP2-I
           CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
           CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
           CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
           CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
        END DO
     END DO
  END DO
  GO TO 116
  112   DO 115 J=2,IPPH
     JC = IPP2-J
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        DO K=1,L1
           CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
           CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
           CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
           CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
        END DO
     END DO
  115   CONTINUE
  116   AR1 = 1.
  AI1 = 0.
  DO L=2,IPPH
     LC = IPP2-L
     AR1H = DCP*AR1-DSP*AI1
     AI1 = DCP*AI1+DSP*AR1
     AR1 = AR1H
     DO IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
        C2(IK,LC) = AI1*CH2(IK,IP)
     END DO
     DC2 = AR1
     DS2 = AI1
     AR2 = AR1
     AI2 = AI1
     DO J=3,IPPH
        JC = IPP2-J
        AR2H = DC2*AR2-DS2*AI2
        AI2 = DC2*AI2+DS2*AR2
        AR2 = AR2H
        DO IK=1,IDL1
           C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
           C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
        END DO
     END DO
  END DO
  DO J=2,IPPH
     DO IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
     END DO
  END DO
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
        CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
        CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
     END DO
  END DO
  IF (IDO .EQ. 1) GO TO 132
  IF (NBD .LT. L1) GO TO 128
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=3,IDO,2
           CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
           CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
           CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
           CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
        END DO
     END DO
  END DO
  GO TO 132
  128   DO 131 J=2,IPPH
     JC = IPP2-J
     DO I=3,IDO,2
        DO K=1,L1
           CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
           CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
           CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
           CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
        END DO
     END DO
  131   CONTINUE
  132   CONTINUE
  IF (IDO .EQ. 1) RETURN
  DO IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  END DO
  DO J=2,IP
     DO K=1,L1
        C1(1,K,J) = CH(1,K,J)
     END DO
  END DO
  IF (NBD .GT. L1) GO TO 139
  IS = -IDO
  DO J=2,IP
     IS = IS+IDO
     IDIJ = IS
     DO I=3,IDO,2
        IDIJ = IDIJ+2
        DO K=1,L1
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  GO TO 143
  139   IS = -IDO
  DO J=2,IP
     IS = IS+IDO
     DO K=1,L1
        IDIJ = IS
  !DIR$ IVDEP
        DO I=3,IDO,2
           IDIJ = IDIJ+2
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        END DO
     END DO
  END DO
  143   RETURN
  END SUBROUTINE RADBG
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADF2(IDO,L1,CC,CH,WA1)
  !***BEGIN PROLOGUE  RADF2
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADF2
  DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           , &
        WA1(*)
  !***FIRST EXECUTABLE STATEMENT  RADF2
  DO K=1,L1
     CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
     CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  END DO
  IF (IDO-2) 107,105,102
  102   IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 108
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CH(I,1,K) = CC(I,K,1)+TI2
        CH(IC,2,K) = TI2-CC(I,K,1)
        CH(I-1,1,K) = CC(I-1,K,1)+TR2
        CH(IC-1,2,K) = CC(I-1,K,1)-TR2
     END DO
  END DO
  GO TO 111
  108   DO 110 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CH(I,1,K) = CC(I,K,1)+TI2
        CH(IC,2,K) = TI2-CC(I,K,1)
        CH(I-1,1,K) = CC(I-1,K,1)+TR2
        CH(IC-1,2,K) = CC(I-1,K,1)-TR2
     END DO
  110   CONTINUE
  111   IF (MOD(IDO,2) .EQ. 1) RETURN
  105   DO 106 K=1,L1
     CH(1,2,K) = -CC(IDO,K,2)
     CH(IDO,1,K) = CC(IDO,K,1)
  106   CONTINUE
  107   RETURN
  END SUBROUTINE RADF2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADF3(IDO,L1,CC,CH,WA1,WA2)
  !***BEGIN PROLOGUE  RADF3
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADF3
  DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           , &
        WA1(*)     ,WA2(*)
  DATA TAUR,TAUI /-.5,.866025403784439/
  !***FIRST EXECUTABLE STATEMENT  RADF3
  DO K=1,L1
     CR2 = CC(1,K,2)+CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2
     CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
     CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  END DO
  IF (IDO .EQ. 1) RETURN
  IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 104
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
     END DO
  END DO
  RETURN
  104   DO 106 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
     END DO
  106   CONTINUE
  RETURN
  END SUBROUTINE RADF3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADF4(IDO,L1,CC,CH,WA1,WA2,WA3)
  !***BEGIN PROLOGUE  RADF4
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADF4
  DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)
  DATA HSQT2 /.7071067811865475/
  !***FIRST EXECUTABLE STATEMENT  RADF4
  DO K=1,L1
     TR1 = CC(1,K,2)+CC(1,K,4)
     TR2 = CC(1,K,1)+CC(1,K,3)
     CH(1,1,K) = TR1+TR2
     CH(IDO,4,K) = TR2-TR1
     CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
     CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  END DO
  IF (IDO-2) 107,105,102
  102   IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 111
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        TR1 = CR2+CR4
        TR4 = CR4-CR2
        TI1 = CI2+CI4
        TI4 = CI2-CI4
        TI2 = CC(I,K,1)+CI3
        TI3 = CC(I,K,1)-CI3
        TR2 = CC(I-1,K,1)+CR3
        TR3 = CC(I-1,K,1)-CR3
        CH(I-1,1,K) = TR1+TR2
        CH(IC-1,4,K) = TR2-TR1
        CH(I,1,K) = TI1+TI2
        CH(IC,4,K) = TI1-TI2
        CH(I-1,3,K) = TI4+TR3
        CH(IC-1,2,K) = TR3-TI4
        CH(I,3,K) = TR4+TI3
        CH(IC,2,K) = TR4-TI3
     END DO
  END DO
  GO TO 110
  111   DO 109 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        TR1 = CR2+CR4
        TR4 = CR4-CR2
        TI1 = CI2+CI4
        TI4 = CI2-CI4
        TI2 = CC(I,K,1)+CI3
        TI3 = CC(I,K,1)-CI3
        TR2 = CC(I-1,K,1)+CR3
        TR3 = CC(I-1,K,1)-CR3
        CH(I-1,1,K) = TR1+TR2
        CH(IC-1,4,K) = TR2-TR1
        CH(I,1,K) = TI1+TI2
        CH(IC,4,K) = TI1-TI2
        CH(I-1,3,K) = TI4+TR3
        CH(IC-1,2,K) = TR3-TI4
        CH(I,3,K) = TR4+TI3
        CH(IC,2,K) = TR4-TI3
     END DO
  109   CONTINUE
  110   IF (MOD(IDO,2) .EQ. 1) RETURN
  105   DO 106 K=1,L1
     TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
     TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
     CH(IDO,1,K) = TR1+CC(IDO,K,1)
     CH(IDO,3,K) = CC(IDO,K,1)-TR1
     CH(1,2,K) = TI1-CC(IDO,K,3)
     CH(1,4,K) = TI1+CC(IDO,K,3)
  106   CONTINUE
  107   RETURN
  END SUBROUTINE RADF4
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  !***BEGIN PROLOGUE  RADF5
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADF5
  DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           , &
        WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
        -.809016994374947,.587785252292473/
  !***FIRST EXECUTABLE STATEMENT  RADF5
  DO K=1,L1
     CR2 = CC(1,K,5)+CC(1,K,2)
     CI5 = CC(1,K,5)-CC(1,K,2)
     CR3 = CC(1,K,4)+CC(1,K,3)
     CI4 = CC(1,K,4)-CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2+CR3
     CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
     CH(1,3,K) = TI11*CI5+TI12*CI4
     CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
     CH(1,5,K) = TI12*CI5-TI11*CI4
  END DO
  IF (IDO .EQ. 1) RETURN
  IDP2 = IDO+2
  IF((IDO-1)/2.LT.L1) GO TO 104
  DO K=1,L1
  !DIR$ IVDEP
     DO I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
        DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
        CR2 = DR2+DR5
        CI5 = DR5-DR2
        CR5 = DI2-DI5
        CI2 = DI2+DI5
        CR3 = DR3+DR4
        CI4 = DR4-DR3
        CR4 = DI3-DI4
        CI3 = DI3+DI4
        CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
        CH(I,1,K) = CC(I,K,1)+CI2+CI3
        TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
        TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
        TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
        TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
        TR5 = TI11*CR5+TI12*CR4
        TI5 = TI11*CI5+TI12*CI4
        TR4 = TI12*CR5-TI11*CR4
        TI4 = TI12*CI5-TI11*CI4
        CH(I-1,3,K) = TR2+TR5
        CH(IC-1,2,K) = TR2-TR5
        CH(I,3,K) = TI2+TI5
        CH(IC,2,K) = TI5-TI2
        CH(I-1,5,K) = TR3+TR4
        CH(IC-1,4,K) = TR3-TR4
        CH(I,5,K) = TI3+TI4
        CH(IC,4,K) = TI4-TI3
     END DO
  END DO
  RETURN
  104   DO 106 I=3,IDO,2
     IC = IDP2-I
  !DIR$ IVDEP
     DO K=1,L1
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
        DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
        CR2 = DR2+DR5
        CI5 = DR5-DR2
        CR5 = DI2-DI5
        CI2 = DI2+DI5
        CR3 = DR3+DR4
        CI4 = DR4-DR3
        CR4 = DI3-DI4
        CI3 = DI3+DI4
        CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
        CH(I,1,K) = CC(I,K,1)+CI2+CI3
        TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
        TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
        TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
        TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
        TR5 = TI11*CR5+TI12*CR4
        TI5 = TI11*CI5+TI12*CI4
        TR4 = TI12*CR5-TI11*CR4
        TI4 = TI12*CI5-TI11*CI4
        CH(I-1,3,K) = TR2+TR5
        CH(IC-1,2,K) = TR2-TR5
        CH(I,3,K) = TI2+TI5
        CH(IC,2,K) = TI5-TI2
        CH(I-1,5,K) = TR3+TR4
        CH(IC-1,4,K) = TR3-TR4
        CH(I,5,K) = TI3+TI4
        CH(IC,4,K) = TI4-TI3
     END DO
  106   CONTINUE
  RETURN
  END SUBROUTINE RADF5
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RADFG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  !***BEGIN PROLOGUE  RADFG
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RADFG
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
        C1(IDO,L1,IP)          ,C2(IDL1,IP), &
        CH2(IDL1,IP)           ,WA(*)
  DATA TPI/6.28318530717959/
  !***FIRST EXECUTABLE STATEMENT  RADFG
  ARG = TPI/FLOAT(IP)
  DCP = COS(ARG)
  DSP = SIN(ARG)
  IPPH = (IP+1)/2
  IPP2 = IP+2
  IDP2 = IDO+2
  NBD = (IDO-1)/2
  IF (IDO .EQ. 1) GO TO 119
  DO IK=1,IDL1
     CH2(IK,1) = C2(IK,1)
  END DO
  DO J=2,IP
     DO K=1,L1
        CH(1,K,J) = C1(1,K,J)
     END DO
  END DO
  IF (NBD .GT. L1) GO TO 107
  IS = -IDO
  DO J=2,IP
     IS = IS+IDO
     IDIJ = IS
     DO I=3,IDO,2
        IDIJ = IDIJ+2
        DO K=1,L1
           CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
           CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
        END DO
     END DO
  END DO
  GO TO 111
  107   IS = -IDO
  DO J=2,IP
     IS = IS+IDO
     DO K=1,L1
        IDIJ = IS
  !DIR$ IVDEP
        DO I=3,IDO,2
           IDIJ = IDIJ+2
           CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
           CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
        END DO
     END DO
  END DO
  111   IF (NBD .LT. L1) GO TO 115
  DO J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=3,IDO,2
           C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
           C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
           C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
           C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
        END DO
     END DO
  END DO
  GO TO 121
  115   DO 118 J=2,IPPH
     JC = IPP2-J
     DO I=3,IDO,2
        DO K=1,L1
           C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
           C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
           C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
           C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
        END DO
     END DO
  118   CONTINUE
  GO TO 121
  119   DO 120 IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  120   CONTINUE
  121   DO 123 J=2,IPPH
     JC = IPP2-J
     DO K=1,L1
        C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
        C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
     END DO
  123   CONTINUE
  !
  AR1 = 1.
  AI1 = 0.
  DO L=2,IPPH
     LC = IPP2-L
     AR1H = DCP*AR1-DSP*AI1
     AI1 = DCP*AI1+DSP*AR1
     AR1 = AR1H
     DO IK=1,IDL1
        CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
        CH2(IK,LC) = AI1*C2(IK,IP)
     END DO
     DC2 = AR1
     DS2 = AI1
     AR2 = AR1
     AI2 = AI1
     DO J=3,IPPH
        JC = IPP2-J
        AR2H = DC2*AR2-DS2*AI2
        AI2 = DC2*AI2+DS2*AR2
        AR2 = AR2H
        DO IK=1,IDL1
           CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
           CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
        END DO
     END DO
  END DO
  DO J=2,IPPH
     DO IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+C2(IK,J)
     END DO
  END DO
  !
  IF (IDO .LT. L1) GO TO 132
  DO K=1,L1
     DO I=1,IDO
        CC(I,1,K) = CH(I,K,1)
     END DO
  END DO
  GO TO 135
  132   DO 134 I=1,IDO
     DO K=1,L1
        CC(I,1,K) = CH(I,K,1)
     END DO
  134   CONTINUE
  135   DO 137 J=2,IPPH
     JC = IPP2-J
     J2 = J+J
     DO K=1,L1
        CC(IDO,J2-2,K) = CH(1,K,J)
        CC(1,J2-1,K) = CH(1,K,JC)
     END DO
  137   CONTINUE
  IF (IDO .EQ. 1) RETURN
  IF (NBD .LT. L1) GO TO 141
  DO J=2,IPPH
     JC = IPP2-J
     J2 = J+J
     DO K=1,L1
  !DIR$ IVDEP
        DO I=3,IDO,2
           IC = IDP2-I
           CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
           CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
           CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
           CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
        END DO
     END DO
  END DO
  RETURN
  141   DO 144 J=2,IPPH
     JC = IPP2-J
     J2 = J+J
     DO I=3,IDO,2
        IC = IDP2-I
        DO K=1,L1
           CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
           CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
           CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
           CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
        END DO
     END DO
  144   CONTINUE
  RETURN
  END SUBROUTINE RADFG
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTB1(N,C,CH,WA,WFAC)
  !***BEGIN PROLOGUE  RFFTB1
  !***REFER TO  RFFTB
  !***ROUTINES CALLED  RADB2,RADB3,RADB4,RADB5,RADBG
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RFFTB1
  DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WFAC(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTB1
  NF = INT(WFAC(2))
  NA = 0
  L1 = 1
  IW = 1
  DO K1=1,NF
     IP = INT(WFAC(K1+2))
     L2 = IP*L1
     IDO = N/L2
     IDL1 = IDO*L1
     IF (IP .NE. 4) GO TO 103
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IF (NA .NE. 0) GO TO 101
     CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     GO TO 102
  101   CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102   NA = 1-NA
     GO TO 115
  103   IF (IP .NE. 2) GO TO 106
     IF (NA .NE. 0) GO TO 104
     CALL RADB2 (IDO,L1,C,CH,WA(IW))
     GO TO 105
  104   CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105   NA = 1-NA
     GO TO 115
  106   IF (IP .NE. 3) GO TO 109
     IX2 = IW+IDO
     IF (NA .NE. 0) GO TO 107
     CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
     GO TO 108
  107   CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108   NA = 1-NA
     GO TO 115
  109   IF (IP .NE. 5) GO TO 112
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IX4 = IX3+IDO
     IF (NA .NE. 0) GO TO 110
     CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 111
  110   CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111   NA = 1-NA
     GO TO 115
  112   IF (NA .NE. 0) GO TO 113
     CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     GO TO 114
  113   CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114   IF (IDO .EQ. 1) NA = 1-NA
  115   L1 = L2
     IW = IW+(IP-1)*IDO
  END DO
  IF (NA .EQ. 0) RETURN
  DO I=1,N
     C(I) = CH(I)
  END DO
  RETURN
  END SUBROUTINE RFFTB1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTB(N,R,WSAVE)
  !***BEGIN PROLOGUE  RFFTB
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  830401   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  J1A1
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Backward transform of a real coefficient array.
  !***DESCRIPTION
  !
  !  Subroutine RFFTB computes the real perodic sequence from its
  !  Fourier coefficients (Fourier synthesis).  The transform is defined
  !  below at output parameter R.
  !
  !  Input Parameters
  !
  !  N       the length of the array R to be transformed.  The method
  !      is most efficient when N is a product of small primes.
  !      N may change so long as different work arrays are provided.
  !
  !  R       a real array of length N which contains the sequence
  !      to be transformed
  !
  !  WSAVE   a work array which must be dimensioned at least 2*N+15
  !      in the program that calls RFFTB.  The WSAVE array must be
  !      initialized by calling subroutine RFFTI(N,WSAVE), and a
  !      different WSAVE array must be used for each different
  !      value of N.  This initialization does not have to be
  !      repeated so long as N remains unchanged.  Thus subsequent
  !      transforms can be obtained faster than the first.
  !      The same WSAVE array can be used by RFFTF and RFFTB.
  !
  !
  !  Output Parameters
  !
  !  R       For N even and For I = 1,...,N
  !
  !           R(I) = R(1)+(-1)**(I-1)*R(N)
  !
  !                plus the sum from K=2 to K=N/2 of
  !
  !                 2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !      For N odd and For I = 1,...,N
  !
  !           R(I) = R(1) plus the sum from K=2 to K=(N+1)/2 of
  !
  !                2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !               -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !   *****  Note:
  !           This transform is unnormalized since a call of RFFTF
  !           followed by a call of RFFTB will multiply the input
  !           sequence by N.
  !
  !  WSAVE   contains results which must not be destroyed between
  !      calls of RFFTB or RFFTF.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  RFFTB1
  !***END PROLOGUE  RFFTB
  DIMENSION       R(*)       ,WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTB
  IF (N .EQ. 1) RETURN
  CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
  RETURN
  END SUBROUTINE RFFTB
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTF1(N,C,CH,WA,WFAC)
  !***BEGIN PROLOGUE  RFFTF1
  !***REFER TO  RFFTF
  !***ROUTINES CALLED  RADF2,RADF3,RADF4,RADF5,RADFG
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RFFTF1
  DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,WFAC(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTF1
  NF = INT(WFAC(2))
  NA = 1
  L2 = N
  IW = N
  DO K1=1,NF
     KH = NF-K1
     IP = INT(WFAC(KH+3))
     L1 = L2/IP
     IDO = N/L2
     IDL1 = IDO*L1
     IW = IW-(IP-1)*IDO
     NA = 1-NA
     IF (IP .NE. 4) GO TO 102
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IF (NA .NE. 0) GO TO 101
     CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     GO TO 110
  101   CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
     GO TO 110
  102   IF (IP .NE. 2) GO TO 104
     IF (NA .NE. 0) GO TO 103
     CALL RADF2 (IDO,L1,C,CH,WA(IW))
     GO TO 110
  103   CALL RADF2 (IDO,L1,CH,C,WA(IW))
     GO TO 110
  104   IF (IP .NE. 3) GO TO 106
     IX2 = IW+IDO
     IF (NA .NE. 0) GO TO 105
     CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
     GO TO 110
  105   CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
     GO TO 110
  106   IF (IP .NE. 5) GO TO 108
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IX4 = IX3+IDO
     IF (NA .NE. 0) GO TO 107
     CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 110
  107   CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 110
  108   IF (IDO .EQ. 1) NA = 1-NA
     IF (NA .NE. 0) GO TO 109
     CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     NA = 1
     GO TO 110
  109   CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
     NA = 0
  110   L2 = L1
  END DO
  IF (NA .EQ. 1) RETURN
  DO I=1,N
     C(I) = CH(I)
  END DO
  RETURN
  END SUBROUTINE RFFTF1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTF(N,R,WSAVE)
  !***BEGIN PROLOGUE  RFFTF
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  830401   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***CATEGORY NO.  J1A1
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Forward transform of a real, periodic sequence.
  !***DESCRIPTION
  !
  !  Subroutine RFFTF computes the Fourier coefficients of a real
  !  perodic sequence (Fourier analysis).  The transform is defined
  !  below at output parameter R.
  !
  !  Input Parameters
  !
  !  N       the length of the array R to be transformed.  The method
  !      is most efficient when N is a product of small primes.
  !      N may change so long as different work arrays are provided
  !
  !  R       a real array of length N which contains the sequence
  !      to be transformed
  !
  !  WSAVE   a work array which must be dimensioned at least 2*N+15
  !      in the program that calls RFFTF.  The WSAVE array must be
  !      initialized by calling subroutine RFFTI(N,WSAVE), and a
  !      different WSAVE array must be used for each different
  !      value of N.  This initialization does not have to be
  !      repeated so long as N remains unchanged.  Thus subsequent
  !      transforms can be obtained faster than the first.
  !      the same WSAVE array can be used by RFFTF and RFFTB.
  !
  !
  !  Output Parameters
  !
  !  R       R(1) = the sum from I=1 to I=N of R(I)
  !
  !      If N is even set L = N/2; if N is odd set L = (N+1)/2
  !
  !        then for K = 2,...,L
  !
  !           R(2*K-2) = the sum from I = 1 to I = N of
  !
  !                R(I)*COS((K-1)*(I-1)*2*PI/N)
  !
  !           R(2*K-1) = the sum from I = 1 to I = N of
  !
  !               -R(I)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !      If N is even
  !
  !           R(N) = the sum from I = 1 to I = N of
  !
  !                (-1)**(I-1)*R(I)
  !
  !   *****  Note:
  !           This transform is unnormalized since a call of RFFTF
  !           followed by a call of RFFTB will multiply the input
  !           sequence by N.
  !
  !  WSAVE   contains results which must not be destroyed between
  !      calls of RFFTF or RFFTB.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  RFFTF1
  !***END PROLOGUE  RFFTF
  DIMENSION       R(*)       ,WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTF
  IF (N .EQ. 1) RETURN
  CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
  RETURN
  END SUBROUTINE RFFTF
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTI1(N,WA,WFAC)
  !***BEGIN PROLOGUE  RFFTI1
  !***REFER TO  RFFTI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !
  !***END PROLOGUE  RFFTI1
  DIMENSION       WA(*)      ,WFAC(*)    ,NTRYH(4)
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
  !***FIRST EXECUTABLE STATEMENT  RFFTI1
  NL = N
  NF = 0
  J = 0
  101   J = J+1
  IF (J-4) 102,102,103
  102   NTRY = NTRYH(J)
  GO TO 104
  103   NTRY = NTRY+2
  104   NQ = NL/NTRY
  NR = NL-NTRY*NQ
  IF (NR) 101,105,101
  105   NF = NF+1
  WFAC(NF+2) = NTRY
  NL = NQ
  IF (NTRY .NE. 2) GO TO 107
  IF (NF .EQ. 1) GO TO 107
  DO I=2,NF
     IB = NF-I+2
     WFAC(IB+2) = WFAC(IB+1)
  END DO
  WFAC(3) = 2
  107   IF (NL .NE. 1) GO TO 104
  WFAC(1) = N
  WFAC(2) = NF
  TPI = 6.28318530717959D0
  ARGH = TPI/FLOAT(N)
  IS = 0
  NFM1 = NF-1
  L1 = 1
  IF (NFM1 .EQ. 0) RETURN
  DO K1=1,NFM1
     IP = INT(WFAC(K1+2))
     LD = 0
     L2 = L1*IP
     IDO = N/L2
     IPM = IP-1
     DO J=1,IPM
        LD = LD+L1
        I = IS
        ARGLD = FLOAT(LD)*ARGH
        FI = 0.
        DO II=3,IDO,2
           I = I+2
           FI = FI+1.D0
           ARG = FI*ARGLD
           WA(I-1) = COS(ARG)
           WA(I) = SIN(ARG)
        END DO
        IS = IS+IDO
     END DO
     L1 = L2
  END DO
  RETURN
  END SUBROUTINE RFFTI1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE RFFTI(N,WSAVE)
  !***BEGIN PROLOGUE  RFFTI
  !***DATE WRITTEN   790601   (YYMMDD)
  !***REVISION DATE  830401   (YYMMDD)
  !***REVISION HISTORY  (YYMMDD)
  !   000330  Modified array declarations.  (JEC)
  !***CATEGORY NO.  J1A1
  !***KEYWORDS  FOURIER TRANSFORM
  !***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
  !***PURPOSE  Initialize for RFFTF and RFFTB.
  !***DESCRIPTION
  !
  !  Subroutine RFFTI initializes the array WSAVE which is used in
  !  both RFFTF and RFFTB.  The prime factorization of N together with
  !  a tabulation of the trigonometric functions are computed and
  !  stored in WSAVE.
  !
  !  Input Parameter
  !
  !  N       the length of the sequence to be transformed.
  !
  !  Output Parameter
  !
  !  WSAVE   a work array which must be dimensioned at least 2*N+15.
  !      The same work array can be used for both RFFTF and RFFTB
  !      as long as N remains unchanged.  Different WSAVE arrays
  !      are required for different values of N.  The contents of
  !      WSAVE must not be changed between calls of RFFTF or RFFTB.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  RFFTI1
  !***END PROLOGUE  RFFTI
  DIMENSION       WSAVE(*)
  !***FIRST EXECUTABLE STATEMENT  RFFTI
  IF (N .EQ. 1) RETURN
  CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
  RETURN
  END SUBROUTINE RFFTI
  !-----------------------------------------------------------------------------

END MODULE MOOSE_CMLIB_FFTPKG
