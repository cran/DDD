
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Example how to apply Fortran code with variable-length parameters
! ---------------------------------------------------------------------
! ------------------------------------------------------------------------

!==========================================================================
! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      INTEGER  :: kk
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)  
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE initmod (steadyparms)
      USE dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 2  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)
      kk = INT(parms(2) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P(3 * (N + 2 + 2 * kk)))

      initialised = .FALSE.
       
      END SUBROUTINE initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE runmod (neq, t, Conc, dConc, yout, ip)
      USE dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec(N + 2 + 2 * kk),muvec(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FF2, FF3

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL fill1d(P, 3 * (N + 2 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

 !  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] 
    !  * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]

      V(1) = 0
      DO I = 2, N + 1 
        V(I) = Conc(I - 1)
      ENDDO
      V(N + 2) = 0
      DO I = 1, N + 2 + 2 * kk
       lavec(I) = P(I)
       muvec(I) = P(I + N + 2 + 2 * kk)
       nn(I)    = P(I + 2 * (N + 2 + 2 * kk))
      ENDDO

      DO I = 2, N + 1 
        FF1 = lavec(I + kk - 1) * nn(I + 2 * kk - 1) * V(I - 1)
        FF2 = muvec(I + kk + 1) * nn(I + 1) * V(I + 1)
        FF3 = (lavec(I + kk) + muvec(I + kk)) * nn(I + kk) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
  
      END SUBROUTINE runmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE runmodbw (neq, t, Conc, dConc, yout, ip)
      USE dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec(N + 2 + 2 * kk),muvec(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FF2, FF3

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL fill1d(P, 3 * (N + 2 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

! dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
! dG = x[1 + (kk == 0)]


      V(1) = 0
      DO I = 2, N 
        V(I) = Conc(I - 1)
      ENDDO
      V(N + 1) = 0
      V(N + 2) = 0
      DO I = 1, N + 2 + 2 * kk
       lavec(I) = P(I)
       muvec(I) = P(I + N + 2 + 2 * kk)
       nn(I)    = P(I + 2 * (N + 2 + 2 * kk))
      ENDDO

      DO I = 2, N - 1
        FF1 = lavec(I + kk) * nn(I + 2 * kk) * V(I + 1)
        FF2 = muvec(I + kk) * nn(I) * V(I - 1)
        FF3 = (lavec(I + kk) + muvec(I + kk)) * nn(I + kk) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
      I = N
      FF1 = lavec(I + kk) * nn(I + 2 * kk) * V(I + 1)
      FF2 = muvec(I + kk) * nn(I) * V(I - 1)
      FF3 = (0 + muvec(I + kk)) * nn(I + kk) * V(I)
      dConc(I - 1) = FF1 + FF2 - FF3
      I = N + 1
      IF (kk .eq. 0) THEN
         dConc(I - 1) = Conc(2)
      ELSE
         dConc(I - 1) = Conc(1)
      ENDIF      

      END SUBROUTINE runmodbw
