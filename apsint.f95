!*******************************************************************************
!
!     SUBROUTINE TO CALCULATE THE INTEGRAL THAT APPEARS IN THE ASYMPTOTIC PERIOD
!     SPACING FORMULA
!     --------------------------------------------------------------------------
!
!     AUTHOR NAME   : KULDEEP VERMA
!     EMAIL ADDRESS : kuldeep@phys.au.dk, kuldeepv89@gmail.com
!
!*******************************************************************************
!
!     nmesh : (input) Number of mesh points in the model.
!     rx : (input) Radial coordinate profile.
!     NX : (input) Brunt Vaisala frequency profile.
!     psint : (output) Period spacing integral.
!     repsint : (output) Relative accuracy achieved on the integral.
!     nerr : (output) Error parameter.
!
!-------------------------------------------------------------------------------



      SUBROUTINE APSINT(nmesh,rx,NX,psint,repsint,nerr)
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: nbsp_order = 2
      INTEGER :: i, nmesh, nfun, nmax, nerr
      INTEGER :: num_of_knots, calc_flag
      INTEGER :: inc(nmesh)

      REAL*8, PARAMETER :: reps = 1.d-8
      REAL*8, PARAMETER :: aeps = 0.d0
      REAL*8 :: rx(nmesh), NX(nmesh), N_div_r(nmesh)
      REAL*8 :: aa(nmesh,3*nbsp_order), coeff(nmesh)
      REAL*8 :: knots(nmesh+2*nbsp_order+1), wk(5*nmesh)
      REAL*8 :: psint, repsint, rx0, rxn
!f2py intent(in)    :: rx, NX
!f2py intent(out)   :: psint, repsint, nerr
!f2py depend(nmesh) :: rx


      !Compute N/r
      DO i = 1, nmesh
        N_div_r(i) = NX(i)/rx(i)
      ENDDO


      !Compute coefficients for B-spline interpolation
      calc_flag = 0
      CALL BSPINT(nmesh,rx,N_div_r,nbsp_order,aa,nmesh,coeff,&
           knots,num_of_knots,calc_flag,inc,wk,nerr)
      !WRITE (*,'(3I10)') nerr, nmesh, num_of_knots
      IF (nerr .GT. 200) THEN
        nerr = -1
        RETURN
      ENDIF


      !Compute the integral
      rx0  = rx(1)
      rxn  = rx(nmesh)
      nmax = 1000000
      CALL ADPINT(psint,rx0,rxn,reps,aeps,repsint,nerr,nfun,nmax)
      !WRITE (*,'(3I10,2E16.8)') nerr, nfun, nmax, psint, repsint
      repsint = repsint/psint
      nerr = 0
      IF (repsint .GT. 1.d-2) nerr = -2
      IF (psint .LT. 0.d0) nerr = -3



      CONTAINS



!-------------------------------------------------------------------------------
!     To integrate a function over finite interval using adaptive control of
!     step size
!-------------------------------------------------------------------------------
!
!	RINT : (output) Calculated value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies specified accuracy was not achieved on
!			at least one subinterval
!		IER=32 implies that this failure occurred more than IFMAX (=5) times
!		IER=325 implies that subroutine failed to attain required
!			accuracy using NMAX function evaluations
!		In all cases DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by subroutine
!	NMAX : (input/output) Maximum number of function evaluations to be tried
!		If NMAX.LE.0 it is set to MAXPT (=100000)
!
!		FUNCTION F(X) must be supplied by the user.
!
!	Required routines : KRONRD (or GAUS16), F

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,IER,NPT,NMAX)
      IMPLICIT NONE

      LOGICAL :: Q

      INTEGER, PARAMETER :: IPMAX=100, IFMAX=5, MAXPT=100000
      INTEGER :: IER, NPT, NMAX, IFAIL, IU, NP

      REAL*8 :: RINT, FINT, XL, XU, REPS, AEPS, AEPSL, DIF, DIF0
      REAL*8 :: RL, RM, RU
      REAL*8 :: XU1(IPMAX)

      IER=0
      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      NPT=0
      RL=XL
      RU=XU
      IU=0

!	To evaluate the integral over [RL,RU]
1000  CALL KRONRD(FINT,RL,RU,DIF0,NP)
!1000  CALL GAUS16(FINT,RL,RU,DIF0,NP)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
!	Q=.TRUE. if the interval cannot be divided further
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU

      IF(DIF0.LT.MAX(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
!	Accept the value of FINT if adequate convergence or if the interval
!	cannot be subdivided further
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.MAX(ABS(RINT)*REPS,AEPSL)) THEN
!	Integration fails to converge on this subinterval. Go to the next subinterval
          IER=31
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
!	If failure is frequent then adjust the convergence criterion.
            IER=32
            AEPSL=DIF*0.5
          ENDIF
        ENDIF

!	If all subintervals are exhausted then return
        IF(IU.LE.0) RETURN

!	otherwise try next subinterval
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE

!	Subdivide the current interval and try again
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
!	If the number of function evaluations has exceeded the limit then
!	try a last call to estimate the integral over the remaining interval
      IER=325
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP)
!      CALL GAUS16(FINT,RL,RU,DIF0,NP)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END SUBROUTINE ADPINT



!-------------------------------------------------------------------------------
!     To integrate a function over a finite interval using Gauss-Kronrod
!     formula For use with ADPINT
!-------------------------------------------------------------------------------
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	N : (output) Number of function evaluations used by subroutine
!	F : (input) Name of the function routine to calculate the integrand
!
!	FUNCTION F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE KRONRD(RI,A,B,DIF,N)
      IMPLICIT NONE

      INTEGER :: N, K

      REAL*8 :: RI, A, B, DIF, AT, BT, F1, F2, FBT, R1
      REAL*8 :: W7(4), A7(4), WK7(4), WK15(4), AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,&
               0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,&
               0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,&
               0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15/0.02293532201052922496D0, 0.10479001032225018383D0,&
               0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15/0.99145537112081263920D0, 0.86486442335976907278D0,&
               0.58608723546769113029D0, 0.20778495500789846760D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
!	7-point Gauss-Legendre formula
        R1=R1+W7(K)*(F1+F2)
!	15-point Kronrod formula
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END SUBROUTINE KRONRD



!-------------------------------------------------------------------------------
!     To calculate coefficients for B-spline interpolation
!-------------------------------------------------------------------------------
!
!     N : (input) Number of entries in the table
!     X : (input) Array of length N containing the abscissas
!     F : (input) Array of length N containing the function values
!     	F(I) is the tabulated function value at X(I).
!     K : (input) Order of B-spline required. K=4 gives cubic B-splines
!     A : (input/output) Real array of length LA*3K containing the
!     	triangular decomposition of equation matrix in band form
!     	For IFLG=2, this array must be supplied, for other values
!     	of IFLG it is calculated by the subroutine
!     LA : (input) The first dimension of A as specified in calling program
!     	LA.GE.N
!     C : (output) Coefficients of expansion, which will be calculated
!     	provided IFLG.NE.1
!     XF : (input/output) Real array of size NO, containing
!     	the knots used for B-spline calculations.
!     	The knots must be distinct and in ascending order.
!     	For IFLG=2, this array must be supplied, for other values
!     	of IFLG it is calculated by the subroutine
!     NO : (input/output) Number of knots for B-splines
!     	For IFLG=2, this number must be supplied, for other values
!     	of IFLG it is calculated by the subroutine
!     IFLG : (input/output) Integer specifying the type of calculation required
!     	IFLG=0 The matrix will be calculated and solved for coefficients
!     	IFLG=1 The matrix will be calculated and triangular decomposition
!     		is obtained, but coefficients are not calculated
!     	IFLG=2 The triangular decomposition of matrix is assumed
!     		to be available in A and coefficients C are calculated
!     	IFLG=-1 same as 0, except that no pivoting will be used
!     INC : (input/output) Integer array containing information about
!     	pivoting during solution of system of linear equations
!     	For IFLG=2, this array must be supplied, for other values
!     	of IFLG it is calculated by the subroutine
!     WK : Scratch array of length 3*N+K+7
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	IER=204 implies N<K or K<2
!     	other values may be set by BSPLIN or GAUBND
!
!     Required routines : BSPLIN, GAUBND

      SUBROUTINE BSPINT(N,X,F,K,A,LA,C,XF,NO,IFLG,INC,WK,IER)
      IMPLICIT NONE

      INTEGER :: N, K, LA, NO, IFLG, IER, I, J, KB, KL, KU
      INTEGER :: IDET, LEFT, NDB, NDERIV, NUM
      INTEGER :: INC(N)

      REAL*8 :: XB, DET
      REAL*8 :: X(N), A(LA,3*K), WK(3*N+K+7), C(N), F(N), XF(NO)

      IF(N.LT.K.OR.K.LT.2) THEN
        IER=204
        RETURN
      ENDIF

      IF(IFLG.LE.1) THEN
      !set up the knots for B-splines by dropping points near the ends
        XF(1)=X(1)
        KL=(K-2)/2
        KU=(K-1)/2
        DO 2000 I=2+KL,N-1-KU
          XF(I-KL)=X(I)
2000    CONTINUE
        XF(N-KL-KU)=X(N)
        NO=N-KL-KU
        NDB=N+2
        NDERIV=0

      !Set up the equation matrix for calculating coefficients of expansion
      !The matrix is in band form A_{i,j} is stored in A(I,J-I+K)
        DO 2500 I=1,N
          XB=X(I)
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB+2),LEFT,IER,&
               WK(2*NDB+2))
          IF(IER.GT.100) RETURN
          DO 2200 J=MAX(1,I-K+1),MIN(N,I+K-1)
            A(I,J-I+K)=WK(J)
2200      CONTINUE
2500    CONTINUE
      ENDIF

      !Solve the system of equations for a band matrix
      NUM=1
      KB=K-1
      DO 3000 I=1,N
        C(I)=F(I)
3000  CONTINUE
      CALL GAUBND(N,KB,NUM,A,C,DET,IDET,INC,LA,IER,IFLG,WK)
      END SUBROUTINE BSPINT



!-------------------------------------------------------------------------------
!     To calculate the B-spline basis functions at a specified point
!-------------------------------------------------------------------------------
!
!     X : (input) Real array of length NX containing the knots.
!     	The knots must be distinct and in ascending order.
!     NX : (input) Number of knots
!     K : (input) Order of B-spline, 0< K <KMAX+1
!     	K=4 gives cubic B-splines
!     XB : (input) The point at which B-spline basis functions are to be evaluated
!     NDERIV : (input) Number of derivatives required
!     	NDERIV.LE.0 only B-splines are calculated
!     	NDERIV=1 first derivative is also calculated
!     	NDERIV>1 first and second derivatives are also calculated
!     B : (output) Array of length NX+K-2 containing the value of
!     	B-spline basis functions
!     DB : (output) Array of length NX+K-2 containing the value of
!     	the first derivative of B-spline basis functions (if NDERIV>0)
!     DDB : (output) Array of length NX+K-2 containing the value of
!     	the second derivative of B-spline basis functions (if NDERIV>1)
!     LEFT : (output) XB is located between X(LEFT) and X(LEFT+1)
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	IER=26 implies XB > X(NX)
!     	IER=27 implies XB < X(1)
!     	IER=203 implies NX<2, K<1 or K>KMAX
!     WK : Real array of length NX+2K+1 used as scratch space
!
!     Required routines : None

      SUBROUTINE BSPLIN(X,NX,K,XB,NDERIV,B,DB,DDB,LEFT,IER,WK)
      IMPLICIT NONE

      INTEGER, PARAMETER :: KMAX = 20
      INTEGER :: NX, K, NDERIV, LEFT, IER
      INTEGER :: I, J , IGH, LOW, LX, MID, NIGH

      REAL*8 :: XB, P1, P2, T1, T2, T3
      REAL*8 :: X(NX), B(NX+K-2), DR(KMAX), DL(KMAX), DB(NX+K-2)
      REAL*8 :: DDB(NX+K-2), WK(-K:NX+K)

      SAVE
      DATA LOW/0/

      IF(NX.LE.1.OR.K.LT.1.OR.K.GT.KMAX) THEN
        IER=203
        RETURN
      ENDIF

      IER=0

      IF(LOW.LT.1.OR.LOW.GE.NX) THEN
      !If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=NX
      ELSE
        IGH=LOW+1
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.&
        (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
      !Extend the range
        IF(XB.GT.X(LOW)) THEN
      !Extend the range on higher side
          IF(IGH.GE.NX) THEN
            IER=26
            LOW=NX-1
          ELSE
            NIGH=MIN(NX,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
      !Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=27
          ELSE
            NIGH=LOW
            LOW=MAX(1,LOW-2*(IGH-LOW))
            IGH=NIGH
            GO TO 1000
          ENDIF
        ENDIF
      ELSE

      !Once the point is bracketed between two tabular points locate it by bisection
1500    IF(IGH-LOW.GT.1.AND.XB.NE.X(LOW)) THEN
          MID=(LOW+IGH)/2
          IF(XB.LE.X(MID).EQV.XB.LE.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

      !Evaluate the B-spline basis functions

      !Define the extra knots on either side of table
      !Note that the program assumes knots from -K+2 to NX+K-1
      !and the B-splines B_{i,k}, i ranges from 1 to NX+K-2
      !The knots are stored in scratch array WK.
      DO 1700 I=1,NX
        WK(I)=X(I)
1700  CONTINUE
      DO 1800 I=1,K
        WK(1-I)=X(1)
        WK(NX+I)=X(NX)
1800  CONTINUE

      DO 1900 I=1,NX+K-2
        B(I)=0.0
        DB(I)=0.0
        DDB(I)=0.0
1900  CONTINUE
      LEFT=LOW
      LX=LOW-1
      J=1
      B(LX+1)=1.

      !The recurrence relation for B-splines
      DO 3000 J=1,K-1
        DR(J)=WK(LOW+J)-XB
        DL(J)=XB-WK(LOW+1-J)
        T1=0.0
        DO 2000 I=1,J
          T2=B(LX+I)/(DR(I)+DL(J+1-I))
          B(LX+I)=T1+T2*DR(I)
          T1=T2*DL(J+1-I)
2000    CONTINUE
        B(LX+J+1)=T1

      !Calculate the first derivative using recurrence relations
        IF(J.EQ.K-2.AND.NDERIV.GT.0) THEN
          T1=0.0
          DO 2200 I=1,J+1
            T2=B(LX+I)/(WK(LOW+I)-WK(LOW+I+1-K))
            DB(LX+I)=(K-1)*(T1-T2)
            T1=T2
2200      CONTINUE
          DB(LX+J+2)=(K-1)*T1
        ENDIF

      !Calculate the second derivative using recurrence relations
        IF(J.EQ.K-3.AND.NDERIV.GT.1) THEN
          T2=0.0
          P1=0.0
          DO 2400 I=1,J+1
            T3=B(LX+I)/(WK(LOW+I)-WK(LOW+I+2-K))
            P2=(T2-T3)/(WK(LOW+I)-WK(LOW+I-K+1))
            DDB(LX+I)=(K-2)*(K-1)*(P1-P2)
            T2=T3
            P1=P2
2400      CONTINUE
          P2=T2/(WK(LOW+J+2)-WK(LOW+J+3-K))
          DDB(LX+J+2)=(K-2)*(K-1)*(P1-P2)
          DDB(LX+J+3)=(K-2)*(K-1)*P2
        ENDIF
3000  CONTINUE

      !For K=2 the first derivative has to be calculated outside the loop
      IF(K.EQ.2.AND.NDERIV.GT.0) THEN
        T2=1/(WK(LOW+1)-WK(LOW+2-K))
        DB(LX+1)=-T2
        DB(LX+2)=T2
      ENDIF

      !For K=3 the second derivative has to be calculated outside the loop
      IF(K.EQ.3.AND.NDERIV.GT.1) THEN
        T3=1./(WK(LOW+1)-WK(LOW+3-K))
        P2= -T3/(WK(LOW+1)-WK(LOW-K+2))
        DDB(LX+1)=-2.*P2
        P1=P2
        P2=T3/(WK(LOW+2)-WK(LOW+3-K))
        DDB(LX+2)=2.*(P1-P2)
        DDB(LX+3)=2.*P2
      ENDIF

      END SUBROUTINE BSPLIN



!-------------------------------------------------------------------------------
!     Solution of a system of linear equations using gaussian elimination for
!     a band matrix
!-------------------------------------------------------------------------------
!
!     N : (input) Number of equations to be solved
!     KB : (input) Bandwidth of matrix A(I,J)=0 if ABS(I-J)>KB
!     NUM : (input) Number of different sets (each with N equations) of
!     	equations to be solved
!     A : (input/output) The matrix of coefficient of size LJ*(3*KB+1)
!     	A(I,J-I+KB+1) is the coefficient of x_j in Ith equation
!     	at output it will contain the triangular decomposition
!     X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!     	X(I,J) is the Ith element of Jth right hand side
!     	at output it will contain the solutions
!     DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!     INC : (output) Integer array of length N containing information about
!     	interchanges performed during elimination
!     LJ : (input) First dimension of arrays A and X in calling program
!     IER : (output) Error flag, IER=0 signifies successful execution
!     	IER=104 implies (N.LE.0 or N.GT.LJ or KB.GT.N)
!     	IER=124 implies some pivot turned out to be zero and hence
!          		matrix must be nearly singular
!     IFLG : (input) Integer variable used as a flag to specify the type
!     	of computation required
!     		If IFLG=-1, both elimination and solution are calculated
!     			without pivoting and IFLG is set to 2
!     	If IFLG=0, both elimination and solution are computed
!     			with partial pivoting and IFLG is set to 2
!     		If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
!     		If IFLG.GE.2 only solution is calculated, the triangular
!     			decomposition should have been calculated earlier
!     WK : Real array of length 3*KB+1 used as scratch space
!
!     Required routines : None

      SUBROUTINE GAUBND(N,KB,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
      IMPLICIT NONE

      INTEGER :: N, KB, NUM, LJ, IDET, IER, IFLG
      INTEGER :: I, J, K, KB1, KM, L, L1
      INTEGER :: INC(N)

      REAL*8 :: DET, R1, T1
      REAL*8 :: A(LJ,3*KB+1), X(LJ,NUM), WK(3*KB+1)

      IF(N.LE.0.OR.N.GT.LJ.OR.KB.GT.N) THEN
        IER=104
        RETURN
      ENDIF

      KB1=KB+1
      IER=124
      IF(IFLG.LE.1) THEN
      !Perform elimination
        DO 2000 I=1,N
          DO 2000 J=2*KB+2,3*KB+1
            A(I,J)=0.0
2000    CONTINUE

        DET=1.0
        IDET=0
        DO 2600 K=1,N-1
      !Find the maximum element in the Kth column
          R1=0.0
          KM=K
          IF(IFLG.GE.0) THEN
            DO 2200 L=K,MIN(N,K+KB)
              IF(ABS(A(L,K-L+KB1)).GT.R1) THEN
                R1=ABS(A(L,K-L+KB1))
                KM=L
              ENDIF
2200        CONTINUE
          ENDIF

          INC(K)=KM
          IF(KM.NE.K) THEN
      !Interchange the rows if needed
            DO 2300 L=K,MIN(N,2*KB+K)
              WK(L-K+1)=A(K,L-K+KB1)
2300        CONTINUE
            DO 2400 L=K,MIN(N,2*KB+K)
              A(K,L-K+KB1)=A(KM,L-KM+KB1)
2400        A(KM,L-KM+KB1)=WK(L-K+1)
            DET=-DET
          ENDIF

          DET=DET*A(K,KB1)
          IF(A(K,KB1).EQ.0.0) RETURN
      !To check for singular or nearly singular matrices replace this
      !statement by, where REPS is approximately \hcross*Max(A(I,J))
      !    IF(ABS(A(K,KB1)).LT.REPS) RETURN
          IF(DET.NE.0.0) THEN

      !Scale the value of the determinant DET
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125D0
              IDET=IDET+5
              GO TO 2350
            ENDIF

2370        IF(ABS(DET).LT.0.03125D0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF

          DO 2500 L=K+1,MIN(N,K+KB)
            A(L,K-L+KB1)=A(L,K-L+KB1)/A(K,KB1)
            DO 2500 L1=K+1,MIN(N,2*KB+K)
2500      A(L,L1-L+KB1)=A(L,L1-L+KB1)-A(L,K-L+KB1)*A(K,L1-K+KB1)
2600    CONTINUE
        DET=DET*A(N,KB1)
        INC(N)=N
      !If pivot is zero then return, IER has been set to 124
        IF(A(N,KB1).EQ.0.0) RETURN
      !To check for singular or nearly singular matrices replace this
      !statement by, where REPS is approximately \hcross*Max(A(I,J))
      !    IF(ABS(A(N,kB1)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
      !Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
      !Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,MIN(N,K+KB)
3000    X(L,J)=X(L,J)-A(L,K-L+KB1)*X(K,J)

      !back-substitution
        X(N,J)=X(N,J)/A(N,KB1)
        DO 3300 K=N-1,1,-1
          DO 3200 L=MIN(N,K+2*KB),K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L-K+KB1)
3300    X(K,J)=X(K,J)/A(K,KB1)
5000  CONTINUE
      END SUBROUTINE GAUBND



!-------------------------------------------------------------------------------
!     TO CALCULATE FUNCTION VALUE USING B-SPLINE EXPANSION.
!-------------------------------------------------------------------------------
!
!     N : (input) Number of knots to define B-splines
!     X : (input) Real array of length N+2K+1 containing the knots.
!     	The knots must be distinct and in ascending order.
!     K : (input) Order of B-splines, K=4 for cubic B-splines
!     NDERIV : (input) Number of derivatives required
!     	For NDERIV.LE.0 only function value is calculated
!     	For NDERIV=1 first derivative is also calculated
!     	For NDERIV>1 both first and second derivatives are calculated
!     WT : (input) Coefficients of B-spline expansion
!     X0 : (input) The point at which expansion has to be evaluated
!     DF : (output) First derivative of function at X0
!     DDF : (output) Second derivative of function at X0
!     WK : Scratch array of length 4N+5K+2
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	Nonzero values of IER may be set by BSPLIN which is called
!
!     BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
!     where \phi_i(x) are B-spline basis functions on knots X
!
!     Required routines : BSPLIN

      FUNCTION BSPEVL(N,X,K,NDERIV,WT,X0,DF,DDF,WK,IER)
      IMPLICIT NONE

      INTEGER :: N, K, NDERIV, IER, I, LEFT, N1, N2, NK

      REAL*8 :: X0, DF, DDF, F, BSPEVL
      REAL*8 :: X(N+K), WT(N+K-2), WK(4*N+5*K+2)

      BSPEVL=0.0
      NK=(N+K)
      CALL BSPLIN(X,N,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LEFT,IER,WK(3*NK))
      IF(IER.GT.100) RETURN

      F=0.0
      DF=0.0
      DDF=0.0
      N1=N+K-1
      N2=2*(N+K)-1
      DO 2000 I=LEFT,LEFT+K-1
        F=F+WT(I)*WK(I)
        DF=DF+WT(I)*WK(N1+I)
        DDF=DDF+WT(I)*WK(N2+I)
2000  CONTINUE
      BSPEVL=F
      END FUNCTION BSPEVL



!-------------------------------------------------------------------------------
!     To compute the integrand in the asymptotic period spacing integral at an
!     arbitrary radial coordinate
!-------------------------------------------------------------------------------
!
!     x : (input) radial coordinate
!
!	Required routines : None

      FUNCTION F(x)
      IMPLICIT NONE

      INTEGER :: nerr
      REAL*8 :: x, dx, ddx, F

      F = BSPEVL(num_of_knots,knots,nbsp_order,0,coeff,x,dx,ddx,wk,nerr)
      !WRITE (1,'(I10,2E16.8)') nerr, x, F

      END FUNCTION F

      END SUBROUTINE APSINT
