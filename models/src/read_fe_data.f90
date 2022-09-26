!-----------------------------------------------------------------------

  MODULE MODULE_READ_FE_DATA

! Subroutine called in INITIALIZE to read in data
! Data used in LINE_EXCIT, subroutine COLFEPL
!
! Module to read in data required for subroutine FE_LEVEL_POPULATIONS
! Reads in atomic data: 
!                       effective collision strengths (gamfepl)
!                       radiative decay rates         (aijfepl)
!                       upper level                   (iupfepl)
!                       lower level                   (jlofepl)
!----------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY      : GET_FILE_NUMBER 

    IMPLICIT NONE
    INCLUDE "precision.f90"

    INTEGER (KIND=LONG),               PUBLIC :: gam_n,aij_n
    INTEGER (KIND=LONG),              PRIVATE :: i,j, up, down
    REAL (KIND=DP),                   PRIVATE :: gam, A
    REAL (KIND=DP),DIMENSION(35,35),   PUBLIC :: gamfepl
    REAL (KIND=DP),DIMENSION(256),     PUBLIC :: aijfepl
    INTEGER, DIMENSION(256),           PUBLIC :: iupfepl, jlofepl

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


  CONTAINS
    
  SUBROUTINE READ_FE_DATA
  
! Initialise
  gamfepl  = 0E0_DP
  aijfepl  = 0E0_DP
  iupfepl  = 0E0_DP
  jlofepl  = 0E0_DP

! Effective collisions strengths (Fe-electron collisions) for 1000 K
  gam_n=GET_FILE_NUMBER()
  open (gam_n, file='input/gamma_eFe', status='old')        
      
        do i=1,595
          read (gam_n,*) down, up, gam
	  gamfepl(down, up) = gam
        end do

	close (gam_n)

! Radiative decay rates, upper and lower levels and transition energies
   
  aij_n=GET_FILE_NUMBER()           
  open (aij_n, file='input/Aval_Fe', status='old') 

        do i=1,256
          read (aij_n,*) up, down, A
          aijfepl(i) = A
          iupfepl(i) = up
          jlofepl(i) = down
       end do

	close (aij_n)

        END SUBROUTINE READ_FE_DATA
  END MODULE MODULE_READ_FE_DATA
