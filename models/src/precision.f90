  !*****************************************************************************
  !** This file is bo be included in each module of the mhd shock code.       **
  !** It contains parameters used for the precision of real and integer data  **
  !** Do : real(KIND=DP) :: x                                                 **
  !**      integer(KIND=LONG) :: y                                            **
  !** The variable minus_infinity is the smallest negative number available   **
  !** The variable plus_infinity is the greatest positive number available    **
  !*****************************************************************************

  INTEGER, PARAMETER :: DP=selected_real_kind(P=15) ! real kind
  INTEGER, PARAMETER :: LONG=KIND(1)                ! integer kind
  REAL(KIND=DP),PARAMETER :: minus_infinity = -HUGE(1._DP) ! smallest < 0. number
  REAL(KIND=DP),PARAMETER :: plus_infinity  =  HUGE(1._DP) ! greatest > 0. number
