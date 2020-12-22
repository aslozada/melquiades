!----------------------------------------------------------------------
!         MELQUIADES:   Metropolis Monte Carlo Program
!----------------------------------------------------------------------
!bop
!
!   Module: m_random
!
!   !Description: This module contains the random number generator RANLUX,
!see Computer Physics communications 79 100. This version of RANLUX code is 
!a translate from FORTRAN77 to Fortran 90 by L. Meissner. In MELQUIADES these
!ranlux is called from a function m_rand()
!\\
!\\
!   !Interface:
!
    module m_random
    use m_kind
!
!   !Public member functions:
    public :: m_rand
!
!   !Revision history:
!   16Aug 2015 Asdrubal Lozada
!
!eop
!----------------------------------------------------------------------
!    implicit none

    integer, parameter :: NSeeds = 25,&
                        & MaxLev = 4,&
                        & LxDflt = 4 ! Luxury level caotic 24 bit
    real, parameter :: TwoP12 = 4096.0
    integer, parameter :: IGiga = 1000000000,&
!                        & JSDFlt = 314159265,&
                        & ITwo24 = 2 ** 24,&
                        & ICons = 2147483563
    integer :: JSDFlt
    integer :: II
    integer, parameter :: Next(NSeeds - 1) = (/ NSeeds - 1, (II, II = 1, NSeeds - 2) /) 
    integer :: I24 = 24,&
             & J24 = 10,&
             & In24 = 0,&
             & Kount = 0,&
             & LuxLev = LxDflt,&
             & MKount = 0 
    integer, dimension(0: MaxLev) :: NDSkip = (/ 0, 24, 73, 199, 365 /)
    integer, save :: NSkip, InSeed
    real :: Carry = 0.0 
    real, save :: Seeds(NSeeds - 1), TwoM24, TwoM12
    logical, save :: NotYet = .true.
    real :: Uni

    private :: RCarry

    contains 

    subroutine RanLux (RVec)
    real, intent(out) :: RVec(:)
    integer :: ISeeds(NSeeds - 1), I, IVec, JSeed, K, LEnv, LP


    JSDFlt = m_seed

    LEnv = SIZE (RVec)
    if( NotYet ) then
      NotYet = .false.
      JSeed = JSDFlt
      InSeed = JSeed
      LuxLev = LxDflt
      NSkip = NDSkip(LuxLev)
      LP = NSkip + NSeeds - 1
      In24 = 0
      Kount = 0
      MKount = 0
      TwoM24 = 1.0

      do I = 1, NSeeds - 1
        TwoM24 = TwoM24 * 0.5
        K = JSeed / 53668
        JSeed = 40014 * (JSeed - K * 53668) - K * 12211
        if ( JSeed < 0) JSeed = JSeed + ICons
        ISeeds(I) = mod(JSeed, ITwo24)
!---------------------------------------------
! Bulid new seed: Exported to box simulation
        m_seed = ISeeds(I)
!---------------------------------------------
      end do
      TwoM12 = TwoM24 * 4096.0
      Seeds = REAL (ISeeds) * TwoM24
      I24 = NSeeds - 1
      J24 = 10
      Carry = MERGE (TwoM24, 0.0, Seeds(NSeeds - 1) == 0.0)
    end if

    do IVec = 1, LEnv
      RVec(IVec) = RCarry (1)
      In24 = In24 + 1
      if ( In24 == NSeeds - 1 ) then
        In24 = 0
        Kount = Kount + NSkip
        Uni = RCarry (NSkip)
      end if
   end do

   where (RVec < TwoM12) RVec = RVec + TwoM24 * Seeds(J24)
       where (Rvec == 0.0) RVec = TwoM24 * TwoM24
    Kount = Kount + LEnv
    if( Kount >= IGiga ) then
      MKount = MKount + 1
      Kount = Kount - IGiga
   end if
  return
  end subroutine RanLux

  SUBROUTINE RLuxIn (ISDext)
    INTEGER, INTENT(in) :: ISDext(:)
    INTEGER :: I, ISD
! start subroutine RLuxIn
    IF (SIZE(ISDext) /= NSeeds) THEN
      RETURN
    END IF
    ! The following IF block added by Phillip Helbig, based on conversation with Fred James;
    ! an equivalent correction has been published by James.
    IF (NotYet) THEN
      NotYet = .FALSE.
    END IF
    TwoM24 = 1.0
    DO I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
    END DO
    TwoM12 = TwoM24 * 4096.0
    Seeds = REAL (ISDext(: NSeeds - 1)) * TwoM24
    Carry = 0.0
    IF (ISDext(NSeeds) < 0) Carry = TwoM24
    ISD = ABS (ISDext(NSeeds))
    I24 = MOD (ISD, 100)
    ISD = ISD / 100
    J24 = MOD (ISD, 100)
    ISD = ISD / 100
    In24 = MOD (ISD, 100)
    ISD = ISD / 100
    LuxLev = ISD
    IF (LuxLev <= MaxLev) THEN
      NSkip = NDSkip(LuxLev)
    ELSE IF (LuxLev >= NSeeds - 1) THEN
      NSkip = LuxLev - NSeeds + 1
    ELSE
      NSkip = NDSkip(MaxLev)
      LuxLev = MaxLev
    END IF
    InSeed = - 1
    RETURN
  END SUBROUTINE RLuxIn

! Ouput Seeds as integers
  SUBROUTINE RLuxUt (ISDext)
    INTEGER, INTENT(out) :: ISDext(:)
! start subroutine RLuxUt
    IF (SIZE(ISDext) /= NSeeds) THEN
      ISDext = 0
      RETURN
    END IF
    ISDext(: NSeeds - 1) = INT (Seeds * TwoP12 * TwoP12)
    ISDext(NSeeds) = MERGE (-ISDext(NSeeds), I24 + 100 * J24 + 10000 * In24 + 1000000 * LuxLev, Carry > 0.0)
    RETURN
  END SUBROUTINE RLuxUt

! Output the "convenient" restart point
  SUBROUTINE RLuxAt (LOut, InOut, K1, K2)
    INTEGER, INTENT(out) :: LOut, InOut, K1, K2
! start subroutine RLuxAt
    LOut = LuxLev
    InOut = InSeed
    K1 = Kount
    K2 = MKount
    RETURN
  END SUBROUTINE RLuxAt

! Initialize from one or three integers
  SUBROUTINE RLuxGo (Lux, Int, K1, K2)
    INTEGER, INTENT(in) :: Lux, Int, K1, K2
    INTEGER :: ISeeds(NSeeds - 1), ILx, I, IOuter, IZip, IZip2, JSeed, K
! start subroutine RLuxGo
    IF (Lux < 0) THEN
      LuxLev = LxDflt
    ELSE IF (Lux <= MaxLev) THEN
      LuxLev = Lux
    ELSE IF (Lux < NSeeds - 1 .OR. Lux > 2000) THEN
      LuxLev = MaxLev
    ELSE
      LuxLev = Lux
      DO ILx = 0, MaxLev
        IF (Lux == NDSkip(ILx) + NSeeds - 1) THEN
          LuxLev = ILx
        END IF
      END DO
    END IF
    IF (LuxLev <= MaxLev) THEN
      NSkip = NDSkip(LuxLev)
    ELSE
      NSkip = LuxLev - 24
    END IF
    In24 = 0
    IF (Int < 0) THEN
    ELSE IF (Int > 0) THEN
      JSeed = Int
    ELSE
      JSeed = JSDFlt
    END IF
    InSeed = JSeed
    NotYet = .FALSE.
    TwoM24 = 1.0
    DO I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
      K = JSeed / 53668
      JSeed = 40014 * (JSeed - K * 53668) - K * 12211
      IF (JSeed < 0) JSeed = JSeed + ICons
      ISeeds(I) = MOD (JSeed, ITwo24)
    END DO
    TwoM12 = TwoM24 * 4096.0
    Seeds = REAL (ISeeds) * TwoM24
    I24 = NSeeds - 1
    J24 = 10
    Carry = MERGE (TwoM24, 0.0, Seeds(NSeeds - 1) == 0.0)

    ! If restarting at a break point, skip K1 + IGIGA * K2
    ! Note that this is the number of numbers delivered to the user PLUS the number skipped (if Luxury > 0) .
    Kount = ABS (K1)
    MKount = ABS (K2)
    IF (Kount + MKount /= 0) THEN
      DO IOuter = 1, MKount + 1
        Uni = RCarry (MERGE (Kount, IGiga, IOuter == MKount + 1))
      END DO
      ! Get the right value of IN24 by direct calculation
      In24 = MOD (Kount, NSkip + NSeeds - 1)
      IF (MKount > 0) THEN
        IZip = MOD (IGiga, NSkip + NSeeds - 1)
        IZip2 = MKount * IZip + In24
        In24 = MOD (IZip2, NSkip + NSeeds - 1)
      END IF
      ! Now IN24 had better be between zero and 23 inclusive
      IF ((In24 < 1) .OR. (In24 >= NSeeds - 1)) THEN
        In24 = 0
      END IF
    END IF
    RETURN
  END SUBROUTINE RLuxGo

  FUNCTION RCarry (N) RESULT (Uni)  ! Private (in module); generates a sequence of N uniform random numbers; returns the last one.
    REAL :: Uni
    INTEGER, INTENT(in) :: N
    INTEGER :: Many
! start function RCarry
    DO Many = 1, N
    ! The Generator proper: "Subtract-with-borrow", as proposed by Marsaglia and Zaman, Florida State University, March, 1989
      Uni = Seeds(J24) - Seeds(I24) - Carry
      IF (Uni < 0.0) THEN
        Uni = Uni + 1.0
        Carry = TwoM24
      ELSE
        Carry = 0.0
      END IF
      Seeds(I24) = Uni
      I24 = Next(I24)
      J24 = Next(J24)
    END DO
    RETURN
  END FUNCTION RCarry

! Function to call ranlux in MELQUIADES

  function m_rand()
  
!  implicit none
  real :: m_rand
  real, allocatable :: RVec(:)

  allocate (Rvec(100))
  
  call RanLux(RVec)

  m_rand = RVec(1)  
  end function m_rand

end module m_random
