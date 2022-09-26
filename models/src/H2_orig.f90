MODULE MODULE_H2
  !*****************************************************************************
  !** The module 'MODULE_H2' contains variables and subroutines related to    **
  !** the H2 molecule, except op_H2 which is in MODULE_PHYS_VAR     **
  !**     * levels, ortho/para, reaction coefficients                         **
  !**     * Aij, quadrupole lines                                             **
  !** It contains also variables and subroutines needed to compute H2 cooling **
  !** and ortho/para ratio.                                                   **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-----------------------------------------
  ! energy levels
  ! variables defined in READ_H2_LEVELS
  !-----------------------------------------

  INTEGER(KIND=LONG)            :: Vmin_H2, Vmax_H2 ! min. and max. values for V (vibration)
  INTEGER(KIND=LONG)            :: Jmin_H2, Jmax_H2 ! min. and max. values for J (rotation)

  LOGICAL :: op_LTE ! .TRUE. if o/p H2 has been initialized to LTE value

  !-----------------------
  ! data type of one level
  !-----------------------
  TYPE TYPE_H2_LEVEL
     INTEGER(KIND=LONG) :: V, J    ! numbers of vibration and rotation
     REAL(KIND=DP)      :: Weight  ! statistical weight : 2J+1 (para) or 3(2J+1) (ortho)
     REAL(KIND=DP)      :: Energy  ! energy (K)
     REAL(KIND=DP)      :: Density ! density of the level (cm-3)
     REAL(KIND=DP)      :: Dens_old ! density at last call to DRIVE
     REAL(KIND=DP)      :: Col_dens ! column density of the level (cm-2)
     REAL(KIND=DP)      :: Form_gr  ! fraction of H2 formed on grains on level V J
  END TYPE TYPE_H2_LEVEL

  ! initialized in READ_H2_LEVELS
  TYPE(TYPE_H2_LEVEL),DIMENSION(:), ALLOCATABLE  :: H2_lev ! vector containing the H2 levels
  INTEGER(KIND=LONG), DIMENSION(:,:), SAVE, ALLOCATABLE :: index_VJ_H2 ! index (1..NH2_lev) of 1 level (V,J)

  INTEGER(KIND=LONG), PRIVATE, PARAMETER :: Num_undefined=-1 ! useful to initialize index_VJ_H2

  ! density of para and ortho-H2 (cm-3)
  REAL(KIND=DP) :: Dens_paraH2, Dens_orthoH2

  !--------------------------
  ! data type of one H2 line
  !--------------------------
  TYPE TYPE_H2_LINE
     CHARACTER(LEN=9)   :: name           ! name of the line
     INTEGER(KIND=LONG) :: Nup, Nlow      ! index of the upper and the lower levels
     REAL(KIND=DP)      :: Aij            ! Aij (s-1)
     REAL(KIND=DP)      :: DeltaE         ! energy of the line : Eup-Elow (K)
     REAL(KIND=DP)      :: emiss     ! emiss of the line (erg/cm3/s)
     REAL(KIND=DP)      :: emiss_old ! emissivity at the last call to DRIVE
     REAL(KIND=DP)      :: intensity      ! intensity integrated along the shock (erg/cm2/s/sr)
  END TYPE TYPE_H2_LINE

  !--- vector containing the H2 lines (dimension NH2_lines) ---
  INTEGER(KIND=LONG) :: NH2_lines      ! number of H2 lines
  TYPE (TYPE_H2_LINE), DIMENSION(:), SAVE, ALLOCATABLE :: H2_lines ! lines

  !---------------------------------------------
  ! collision rates
  ! read in READ_H2_RATES, used in EVOLUTION_H2
  !---------------------------------------------

  ! their dimension used are (4,NH2_lev,NH2_lev), all values outside are rejected
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_H_H2  ! collisions H-H2
  LOGICAL, PRIVATE, DIMENSION(:,:),         ALLOCATABLE :: mask_H_H2  ! collisions H-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_He_H2 ! collisions He-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_H2_H2 ! collisions H2-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(25,0:16,0:36) :: r_raw_GR_H2 ! collisions GR-H2 (raw rates)
  REAL(KIND=DP), PRIVATE, DIMENSION(25,0:16,0:36) :: d2r_raw_GR_H2 ! collisions GR-H2 (raw rates)
  REAL(KIND=DP), PRIVATE, DIMENSION(0:16,0:36) :: vin_GR_H2 ! corresponding velocity
  REAL(KIND=DP), PRIVATE, DIMENSION(0:16,0:36) :: pgr0_GR_H2 ! corresponding rate

  !--------------------------------------------
  ! evolution terms for radiative transitions
  ! read in READ_H2_LINES
  !--------------------------------------------

  ! Aij_H2(Nup,Nlow) : Aij (s-1) of the line Nup -> Nlow (N is the index of the level)
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Aij_H2
  ! sum of the Aij (s-1) of the lines starting from the level Nup
  REAL(KIND=DP),DIMENSION(:), ALLOCATABLE :: sum_Aij_H2


  !-------------------------------------------------
  ! Data from Gerlich for ortho-para transfer by H+
  ! with David Flower modification for ortho/para
  ! alternation (3/1).
  !-------------------------------------------------
  REAL(KIND=DP), PRIVATE, PARAMETER, DIMENSION(29) :: Gerlich = &
       (/2.024D-10, 9.575D-10, 2.729D-10, 7.910D-10 &
       , 2.246D-10, 6.223D-10, 1.828D-10, 4.975D-10, 1.484D-10 &
       , 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10 &
       , 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10 &
       , 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10 &
       , 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10 /)

  !--------------------------------------------------
  ! source terms for H2 level population (cm-3.s-1)
  ! and H2 internal energy (erg.cm-3.s-1)
  !--------------------------------------------------
  REAL(KIND=DP), public, parameter         :: H2_dissoc = 4.4781  ! H2 dissociation energy (in eV)
  REAL(KIND=DP), public                    :: H2_int_E            ! H2 internal energy at formation on grains (in K)
  REAL(KIND=DP)                            :: H2_energy=0.0_DP    ! calculated in COMPUTE_H2
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YN_rovib_H2         ! change in H2 levels number density (cm-3.s-1)

  !--------------------------------------------------------------
  ! source terms of level selective chemial reactions involving H2
  ! 22 juin 2001 - JLB
  ! One type only yet: collisional dissociation of H2
  ! Allocated in : INITIALIZE_ROVIB_H2 (H2.f90)
  ! Used in : DIFFUN & CHEMISTRY
  !---------------------------------------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_ch_H2        ! total selective reaction rate
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_ne_H2        ! total selective reaction rate (neutrals only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_io_H2        ! total selective reaction rate (ions only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_el_H2        ! total selective reaction rate (electrons only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_rx_H2        ! partial selective reaction rate
  REAL(KIND=DP) :: Sel_tot_H2       ! Collisional dissociation of H2 (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tne_H2       ! Collisional dissociation of H2 by neutrals (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tio_H2       ! Collisional dissociation of H2 by ions (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tel_H2       ! Collisional dissociation of H2 by electrons (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: For_gr_H2        ! Formation of H2 on grains (cm-3 s-1)

  !-----------------------------------
  ! cooling rate (erg.cm-3.s-1)
  ! due to H2 lines
  ! calculated in COMPUTE_OP_H2
  !-----------------------------------
  REAL(KIND=DP) :: cooling_H2 ! cooling rate of the neutral fluid due to H2 lines (erg/cm3/s)

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE READ_H2_LEVELS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 levels : V, J, weight, energy.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    H2_lev, index_VJ_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR, ONLY: NH2_lev, NH2_lines_out
    IMPLICIT NONE

    CHARACTER(len=1) :: charact
    INTEGER(KIND=LONG) :: i, ii, V, J
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_lev='input/H2_levels_Evj.in'
    INTEGER                     :: file_H2_lev

    ! initialization
    ALLOCATE (H2_lev(NH2_lev))
    H2_lev(:)%V=0
    H2_lev(:)%J=0
    H2_lev(:)%weight=Zero
    H2_lev(:)%Energy=Zero
    H2_lev(:)%density=Zero
    H2_lev(:)%Dens_old=Zero
    H2_lev(:)%Col_dens=Zero
    H2_lev(:)%Form_gr=Zero

    ! file opening
    file_H2_lev = GET_FILE_NUMBER()
    OPEN(file_H2_lev,file=name_file_H2_lev,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! comments
    DO i=1,6
       READ(file_H2_lev,'(A1)')charact
    END DO

    DO i=1, NH2_lev
       READ(file_H2_lev,*) &
            ii, &
            H2_lev(i)%V, &
            H2_lev(i)%J, &
            H2_lev(i)%weight, &
            H2_lev(i)%Energy
       ! check if all lines have been correctly read, we must have i=ii
       IF (i /= ii) STOP "*** WARNING, error de read des niveaux de H2"
    END DO

    ! min. and max. values of V and J for these levels
    Vmin_H2=MINVAL(H2_lev(:)%V)
    Vmax_H2=MAXVAL(H2_lev(:)%V)
    Jmin_H2=MINVAL(H2_lev(:)%J)
    Jmax_H2=MAXVAL(H2_lev(:)%J)

    ! allocation and filling of the table index_VJ_H2
    ALLOCATE(index_VJ_H2(Vmin_H2:Vmax_H2, Jmin_H2:Jmax_H2))
    index_VJ_H2(:,:)=num_undefined
    DO i=1, NH2_lev
       index_VJ_H2(H2_lev(i)%V,H2_lev(i)%J)=i
    ENDDO

    ! file closure
    CLOSE(file_H2_lev)

  END SUBROUTINE READ_H2_LEVELS


  SUBROUTINE INITIALIZE_ROVIB_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    initialize population of the ro-vibrational levels of H2
    !    2 possibilities :
    !          * op_H2=op_H2 read in READ_PARAMETERS
    !          * op_H2=op LTE at temperature Tn (if op_H2 > op_H2_LTE)
    !    note : op_H2 is (re)-calculated
    ! subroutine/function needed :
    !    COMPUTE_OP_H2
    ! input variables :
    ! ouput variables :
    ! results :
    !    H2_lev, op_H2
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : Tn, op_H2, NH2_lev, iforH2
    USE MODULE_CHEMICAL_SPECIES, ONLY : speci, ind_H2
    USE MODULE_CONSTANTS, ONLY : Zero, kB, EVerg

    IMPLICIT NONE
    CHARACTER(len=7) :: op_choice
    INTEGER(KIND=LONG) :: i
    REAL(KIND=DP) :: Zortho, Zpara, weight_ortho, weight_para
    REAL(KIND=DP),PARAMETER :: op_H2_LTE=999._DP
    REAL(KIND=DP),PARAMETER :: population_limit=1.D-20 ! lower limit for H2 population
    REAL(KIND=DP) :: int_enrg, T_ini, tform, tf_min, tf_max, res

    ALLOCATE (YN_rovib_H2(NH2_lev))
    ALLOCATE (Sel_ch_H2(NH2_lev))
    ALLOCATE (Sel_ne_H2(NH2_lev))
    ALLOCATE (Sel_io_H2(NH2_lev))
    ALLOCATE (Sel_el_H2(NH2_lev))
    ALLOCATE (Sel_rx_H2(NH2_lev))

    ! which o/p to take ?
    IF (op_H2 > op_H2_LTE) THEN
       op_LTE=.TRUE.
    ELSE
       op_LTE=.FALSE.
    END IF

    ! computes the populations, given the o/p choice

    IF (OP_LTE) THEN
       !--------------------------------------------------------
       !--- the o/p has it's LTE value at the temperature Tn ---
       !--------------------------------------------------------
       H2_lev(:)%density=H2_lev(:)%weight * &
            EXP(-(H2_lev(:)%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)

    ELSE
       !-----------------------------------------------------
       !--- the o/p has the value read in READ_PARAMETERS ---
       !-----------------------------------------------------
       ! partition functions and weights for ortho-H2 and para-H2
       Zortho=Zero
       T_ini = Tn
       DO i=index_VJ_H2(0,1),NH2_lev,2 ! first ortho level : (V=0,J=1)
          Zortho = Zortho + H2_lev(i)%weight * &
               EXP(-(H2_lev(i)%energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)
       END DO
       weight_ortho=op_H2/(1._DP+op_H2)/Zortho

       Zpara=Zero
       DO i=index_VJ_H2(0,0),NH2_lev,2 ! first para level : (V=0,J=0)
          Zpara = Zpara + H2_lev(i)%weight * &
               EXP(-(H2_lev(i)%energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)
       END DO
       weight_para=1._DP/(1._DP+op_H2)/Zpara

       ! populations for ortho levels
       WHERE (MOD(H2_lev(:)%J,2) > 0)
          H2_lev(:)%density=weight_ortho * &
               H2_lev(:)%weight * &
               EXP(-(H2_lev(:)%energy-H2_lev(index_VJ_H2(0,0))%energy)/Tn)
       END WHERE

       ! populations for para levels
       WHERE (MOD(H2_lev(:)%J,2) == 0)
          H2_lev(:)%density=weight_para * &
               H2_lev(:)%weight * &
               EXP(-(H2_lev(:)%energy-H2_lev(index_VJ_H2(0,0))%energy)/Tn)
       END WHERE

    END IF

    ! avoids too small numbers => set lower limit to the populations
    WHERE (H2_lev(:)%density < population_limit)
       H2_lev(:)%density=population_limit
    END WHERE

!   print *, H2_lev%density

!  Formation of H2 on grains :
!   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)
!   1: Proportional to Boltzman distribution at 17249 K
!   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
!   3: v = 6, J = 0,1
!   4: fraction = relative populations at t, initialised as H2_lev%density
!                 and changed during integration

    IF (iforH2 == 0) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      tform = int_enrg
      tf_min = 0.0_dp
      tf_max = 3.0_DP * int_enrg
      res = 1.0_dp

      do while (abs(res) > 1.0e-3_DP)
        res = SUM(H2_lev%weight * (int_enrg-H2_lev%Energy) &
            * EXP(-H2_lev%Energy/tform))
        if (res < 0.0_dp) then
          tf_max = tform
          tform = 0.5_DP * (tf_min + tf_max)
        else
          tf_min = tform
          tform = 0.5_DP * (tf_min + tf_max)
        endif
      end do

      int_enrg = tform

      H2_lev%Form_gr = H2_lev%weight * &
           EXP(-(H2_lev%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ELSE IF (iforH2 == 1) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      H2_lev%Form_gr = H2_lev%weight * &
           EXP(-(H2_lev%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ELSE IF (iforH2 == 2) then

      H2_lev(index_VJ_H2(14,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(14,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 3) then

      H2_lev(index_VJ_H2(6,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(6,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 4) then

      H2_lev%Form_gr = H2_lev%density
      H2_lev%Form_gr = H2_lev%Form_gr / SUM(H2_lev%Form_gr)

    ENDIF

    H2_int_E = SUM(H2_lev%Form_gr * H2_lev%Energy)

    ! normalization at H2 density
    H2_lev(:)%density=H2_lev(:)%density &
         *speci(ind_H2)%density/(SUM(DBLE(H2_lev(:)%density)))

    ! op_H2 is (re)-calculated for checking
    CALL COMPUTE_OP_H2
    IF (op_LTE) WRITE(*,'("--- LTE chosen for ortho:para-H2, &
         &o/p = ",ES7.1," ---")')op_H2

  END SUBROUTINE INITIALIZE_ROVIB_H2



  SUBROUTINE READ_H2_RATES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read collision rates for H-H2, He-H2 and H2-H2.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    rate_H_H2, rate_He_H2, rate_H2_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR, ONLY : NH2_lev, H_H2_flag
    USE NUM_REC, ONLY: spline
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: Nrate, i, LevU, LevL, Vup, Jup, Vlow, Jlow
    CHARACTER(LEN=25)           :: name_file_H_H2
    INTEGER                     :: file_H_H2
    CHARACTER(LEN=*), PARAMETER :: name_file_He_H2='input/coeff_He_H2.in'
    INTEGER                     :: file_He_H2
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_H2='input/coeff_H2_H2.in'
    INTEGER                     :: file_H2_H2
    CHARACTER(LEN=*), PARAMETER :: name_file_GR_H2='input/coeff_GR_H2.in'
    INTEGER                     :: file_GR_H2
    integer :: ii, jj, kk
    integer :: iju, ivu
    REAL(KIND=DP) :: yp1, ypn
    REAL(KIND=DP) :: toto1, toto2, toto3, toto4
    REAL(KIND=DP), DIMENSION(1:25) :: vin
    data yp1, ypn / 1.0d30, 1.0d30 /
    data vin / 2.d0,4.d0,6.d0,8.d0,&
               10.d0,12.d0,14.d0,16.d0,18.d0,&
               20.d0,22.d0,24.d0,26.d0,28.d0,&
               30.d0,32.d0,34.d0,36.d0,38.d0,&
               40.d0,42.d0,44.d0,46.d0,48.d0,50.d0/


    ! initialization
    ALLOCATE (rate_H_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (mask_H_H2(NH2_lev,NH2_lev))
    ALLOCATE (rate_He_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (rate_H2_H2(4,NH2_lev,NH2_lev))
    rate_H_H2(1,:,:)=-50._DP
    rate_H_H2(2:4,:,:)=Zero
    mask_H_H2(:,:)= .FALSE.
    rate_He_H2(1,:,:)=-50._DP
    rate_He_H2(2:4,:,:)=Zero
    rate_H2_H2(1,:,:)=-50._DP
    rate_H2_H2(2:4,:,:)=Zero
    r_raw_GR_H2(:,:,:)=Zero
    d2r_raw_GR_H2(:,:,:)=Zero
    vin_GR_H2(:,:)=0

    !-----------------------
    ! rate H-H2
    !-----------------------

    ! We read either DRF file or MM file or BOTH
    ! If BOTH, then first MM then DRF so that Quantum rates override Semiclassical ones

    if (H_H2_flag == "MM" .OR. H_H2_flag == "BOTH") then
      name_file_H_H2 = 'input/coeff_H_H2_MM.in'
      file_H_H2 = GET_FILE_NUMBER()
      OPEN(file_H_H2,file=name_file_H_H2,status='OLD',&
           access='SEQUENTIAL',form='FORMATTED',action='READ')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    if (H_H2_flag == "DRF" .OR. H_H2_flag == "BOTH") then
      name_file_H_H2 = 'input/coeff_H_H2_DRF.in'
      file_H_H2 = GET_FILE_NUMBER()
      OPEN(file_H_H2,file=name_file_H_H2,status='OLD',&
           access='SEQUENTIAL',form='FORMATTED',action='READ')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
           mask_H_H2(LevU,LevL) = .TRUE.
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    !-------------
    ! rate He-H2
    !-------------
    ! file opening
    file_He_H2 = GET_FILE_NUMBER()
    OPEN(file_He_H2,file=name_file_He_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! commentss
    DO i=1,4
       READ(file_He_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_He_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_He_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_He_H2(1,LevU,LevL) = toto1
         rate_He_H2(2,LevU,LevL) = toto2
         rate_He_H2(3,LevU,LevL) = toto3
         rate_He_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_He_H2)

    !-------------
    ! rate H2-H2
    !-------------
    ! file opening
    file_H2_H2 = GET_FILE_NUMBER()
    OPEN(file_H2_H2,file=name_file_H2_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    DO i=1,4
       READ(file_H2_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_H2_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_H2_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_H2_H2(1,LevU,LevL) = toto1
         rate_H2_H2(2,LevU,LevL) = toto2
         rate_H2_H2(3,LevU,LevL) = toto3
         rate_H2_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_H2_H2)

    ! Raw rates - file opening Para, then Ortho
    file_GR_H2 = GET_FILE_NUMBER()
    OPEN(file_GR_H2,file="input/PH2GR.DAT",status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    DO ii=1,25
       READ(file_GR_H2,*)
       DO kk=0,36,2
         READ(file_GR_H2,*) iju, (r_raw_GR_H2(ii,jj,kk),jj=0,16)
       END DO
    END DO
    ! file closure
    CLOSE(file_GR_H2)

    file_GR_H2 = GET_FILE_NUMBER()
    OPEN(file_GR_H2,file="input/OH2GR.DAT",status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    DO ii=1,25
       READ(file_GR_H2,*)
       DO kk=1,35,2
         READ(file_GR_H2,*) iju, (r_raw_GR_H2(ii,jj,kk),jj=0,16)
       END DO
    END DO
    ! file closure
    CLOSE(file_GR_H2)

    DO ii=3,NH2_lev
      iju = H2_lev(ii)%J
      ivu = H2_lev(ii)%V
      DO jj=1,25
        if (r_raw_GR_H2(jj,ivu,iju) > 0.0_DP) then
          vin_GR_H2(ivu,iju) = vin(jj)
          pgr0_GR_H2(ivu,iju) = r_raw_GR_H2(jj,ivu,iju)
          exit
        endif
      END DO
      call spline (vin(jj:),r_raw_GR_H2(jj:,ivu,iju),yp1,ypn,d2r_raw_GR_H2(jj:,ivu,iju))
    END DO

  END SUBROUTINE READ_H2_RATES


  SUBROUTINE READ_H2_LINES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 lines (levels, energies, Aij).
    ! subroutine/function needed :
    !    NAME_H2_LINE
    ! input variables :
    ! ouput variables :
    ! results :
    !    Aij_H2, sum_Aij_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR, ONLY : NH2_lev, NH2_lines_out
    IMPLICIT NONE
    REAL(KIND=DP) :: Aij_S, Aij_Q, Aij_O
    INTEGER(KIND=LONG) :: i, Vup, Vlow, Jup, Jlow, Nup, Nlow, line_number
    CHARACTER(LEN=1) :: line_type
    CHARACTER(LEN=*), PARAMETER :: name_file_Aij_H2='input/Aij_H2.in'
    CHARACTER(LEN=*), PARAMETER :: format_Aij_H2='(3I4,3D15.6)'
    INTEGER                     :: file_Aij_H2

    ! initialization
    ALLOCATE (Aij_H2(NH2_lev,NH2_lev), sum_Aij_H2(NH2_lev))
    Aij_H2(:,:)=Zero
    sum_Aij_H2(:)=Zero

    ! file opening
    file_Aij_H2 = GET_FILE_NUMBER()
    OPEN(file_Aij_H2,file=name_file_Aij_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    DO i=1,5
       READ(file_Aij_H2,*)
    END DO

    ! read line per line. There is max. 3 Aij per level (lines S, Q, O)
    ! the file is in increasing Vup order
    Vup=0
    DO WHILE (Vup <= Vmax_H2)
       READ(file_Aij_H2,format_Aij_H2) &
            Vup, Vlow, Jup, Aij_S, Aij_Q, Aij_O
       ! chock if the level is in the model
       IF ((Vup <= Vmax_H2) .AND. (Vlow <= Vmax_H2) .AND. (Jup <=Jmax_H2)) THEN
          Nup=index_VJ_H2(Vup,Jup) ! index of the upper level
          IF ((Nup > 0) .AND. (Nup <= NH2_lev)) THEN
             ! S line
             Jlow=Jup-2
             IF ((Jlow >=0) .AND. (Aij_S /= Zero)) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_S
             ENDIF
             ! Q line
             Jlow=Jup
             IF (Aij_Q /= Zero) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_Q
             ENDIF
             ! O line
             Jlow=Jup+2
             IF ((Jlow <=JMAX_H2) .AND. (Aij_O /= Zero)) THEN
                Nlow=index_VJ_H2(Vlow,Jlow)
                IF ((Nlow > 0) .AND. (Nlow <= NH2_lev)) Aij_H2(Nup,Nlow)=Aij_O
             ENDIF
          ENDIF
       ENDIF
    END DO

    ! file closure
    CLOSE(file_Aij_H2)

    ! computes the sum of the Aij of the lines starting from the level Nup
    ! starts at (V=0,J=2)
    DO Nup = index_VJ_H2(0,2), NH2_lev
       sum_Aij_H2(Nup)=SUM(Aij_H2(Nup,1:Nup-1))
    END DO
    ! counts the number of transitions in the model (Aij_H2 > 0)
    NH2_lines=COUNT(mask=(Aij_H2 > Zero))
    NH2_lines_out = min(NH2_lines,NH2_lines_out)
    ! allocation and initialization of H2_lines
    ALLOCATE(H2_lines(NH2_lines))
    H2_lines(:)%name=''
    H2_lines(:)%Nup=0
    H2_lines(:)%Nlow=0
    H2_lines(:)%emiss=Zero
    H2_lines(:)%emiss_old=Zero
    H2_lines(:)%intensity=Zero
    H2_lines(:)%DeltaE=Zero
    H2_lines(:)%Aij=Zero

    ! vector H2_lines is filled in order of increasing Eup
    line_number=0
    DO Nup=index_VJ_H2(0,2),NH2_lev
       DO Nlow=1, Nup-1
          IF (Aij_H2(Nup,Nlow) > Zero) THEN ! tests if line exists
             line_number=line_number+1
             Vup=H2_lev(Nup)%V
             Jup=H2_lev(Nup)%J
             Vlow=H2_lev(Nlow)%V
             Jlow=H2_lev(Nlow)%J
             ! finds the type of the line (S, Q, O)
             SELECT CASE(Jup-Jlow)
             CASE (2)
                line_type='S'
             CASE (0)
                line_type='Q'
             CASE (-2)
                line_type='O'
             CASE DEFAULT
             END SELECT
             ! fills H2_lines
             H2_lines(line_number)%name=NAME_H2_LINE(Vup=Vup,Vlow=Vlow,Jlow=Jlow,&
                  line_type=line_type)
             H2_lines(line_number)%Nup=Nup
             H2_lines(line_number)%Nlow=Nlow
             H2_lines(line_number)%Aij=Aij_H2(Nup,Nlow)
             H2_lines(line_number)%DeltaE=H2_lev(Nup)%energy-H2_lev(Nlow)%energy
          ENDIF
       END DO
    END DO

  END SUBROUTINE READ_H2_LINES


  SUBROUTINE COMPUTE_OP_H2
    !---------------------------------------------------------------------------
    ! called by :
    !     INIT_ROVIB_H2
    !     DIFFUN
    ! purpose :
    !     Computes ortho:para H2 ratio from H2_lev%density. Densities of
    !     ortho-H2 and para-H2 are also calculated.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     op_H2, Dens_orthoH2, Dens_paraH2
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : op_H2
    IMPLICIT NONE

    Dens_paraH2  = SUM(H2_lev(:)%density,mask=(MOD(H2_lev(:)%J,2)==0))
    Dens_orthoH2 = SUM(H2_lev(:)%density,mask=(MOD(H2_lev(:)%J,2)>0))
    op_H2 = Dens_orthoH2/Dens_paraH2

  END SUBROUTINE COMPUTE_OP_H2


  FUNCTION NAME_H2_LINE(Vup,Vlow,Jlow,line_type) RESULT(name)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_H2_LINES
    ! purpose :
    !    computes the name of one line of the form 1-0S(1)
    ! subroutine/function needed :
    ! input variables :
    !    * Vup, Jup  -> (V,J) of the upper level  (two digits : < 100)
    !    * Jlow      -> J of th lower level       (two digits : < 100)
    !    * line_type -> 'S', 'O', 'Q'
    ! output variables :
    ! results :
    !    name -> name of the line, without blanks
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: Vup, Vlow, Jlow
    CHARACTER(LEN=1), INTENT(in)   :: line_type ! 'S', 'O' or 'Q'
    CHARACTER(LEN=2) :: str_Vup, str_Vlow,str_Jlow
    CHARACTER(LEN=9) :: name

    WRITE(str_Vup,'(I2)')Vup
    WRITE(str_Vlow,'(I2)')Vlow
    WRITE(str_Jlow,'(I2)')Jlow

    ! return the name without blanks
    name = TRIM(ADJUSTL(str_Vup)) // '-' // TRIM(ADJUSTL(str_Vlow))  // line_type // &
         '(' // TRIM(ADJUSTL(str_Jlow)) // ')'

  END FUNCTION NAME_H2_LINE



  SUBROUTINE EVOLUTION_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    COMPUTE_H2
    ! purpose :
    !    Calculate the source term (cm-3.s-1) for rovibrational levels of H2,
    !    using Aij_H2 (radiation, s-1) read in READ_H2_LINE and Cij_H2
    !    (collisions, s-1) calculated in this subroutine.
    !    The result, YN_rovib_H2, is used in COMPUTE_H2 and in DIFFUN.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : Tn, ABS_DeltaV, NH2_lev, distance, H_H2_flag
    USE MODULE_GRAINS, ONLY : Rsquare_grain, Dens_grain
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_H, Dens_H2, Dens_Hplus, Dens_He
    USE MODULE_CONSTANTS, ONLY: pi, Zero
    USE NUM_REC, ONLY: splint
    USE MODULE_DEBUG_JLB

    IMPLICIT NONE
    !--- probability of excitation by collision (s-1) ---
    REAL(KIND=DP), DIMENSION(NH2_lev,NH2_lev) :: Cij_H2
    !--- sum of Cij_H2 (s-1) starting from a given level ---
    REAL(KIND=DP), DIMENSION(NH2_lev) :: sum_Cij_H2
    !--- sum of the radiative and collisional terms ---
    REAL(KIND=DP), DIMENSION(NH2_lev,NH2_lev) :: Aij_plus_Cij_H2
    INTEGER(KIND=LONG) :: Nup, Jup, Vup ! upper level
    INTEGER(KIND=LONG) :: Nlow, Jlow, Vlow ! lower level
    INTEGER(KIND=LONG) :: ABS_deltaJ, Level1, Level2, i
    REAL(KIND=DP) :: deexcit_to_excit, excit_rate, deexcit_rate
    REAL(KIND=DP) :: decal, Tdd, Tdd2, gam, gam1, gam2, delta_E
    REAL(KIND=DP) :: coeff_H_NR, coeff_H_R, coeff_H, coeff_He, coeff_H2, coeff_Hplus
    REAL(KIND=DP) :: dv_kms, Vth, prob_excit, excit_grain, Vin_H2
    REAL(KIND=DP), DIMENSION(1:25) :: vin

    data vin / 2.d0,4.d0,6.d0,8.d0,&
               10.d0,12.d0,14.d0,16.d0,18.d0,&
               20.d0,22.d0,24.d0,26.d0,28.d0,&
               30.d0,32.d0,34.d0,36.d0,38.d0,&
               40.d0,42.d0,44.d0,46.d0,48.d0,50.d0/

    !--- initialization ---
    Cij_H2          = Zero ! collisions (2 dimensions)
    sum_Cij_H2      = Zero ! collisions (1 dimension)
    Aij_plus_Cij_H2 = Zero ! collisions + radiation (2 dimensions)
    YN_rovib_H2     = Zero ! source terms for H2-level populations (cm-3.s-1)

    !---------------------------------------------------------
    ! A. Calculate the probability of excitation by collisions
    !    Cij_H2 (s-1)
    !---------------------------------------------------------
    DO Nup = 2, NH2_lev ! loop on upper levels
       Jup = H2_lev(Nup)%J
       Vup = H2_lev(Nup)%V
       DO Nlow = 1, Nup-1 ! loop on lower levels
          Jlow = H2_lev(Nlow)%J
          Vlow = H2_lev(Nlow)%V
          ABS_deltaJ = ABS(Jup - Jlow)
          delta_E = H2_lev(Nup)%energy - H2_lev(Nlow)%energy

          !--- deexcit_to_excit is the Boltzmann factor allowing to ---
          !--- get the excitation rate from the de-excitation rate. ---
!         deexcit_to_excit = EXP(-MIN(delta_E/Tn, 180._DP)) * &
!              H2_lev(Nup)%weight/ H2_lev(Nlow)%weight
          deexcit_to_excit = EXP(-delta_E/Tn) * &
               H2_lev(Nup)%weight/ H2_lev(Nlow)%weight

          !--------------------------------------------------------------------------
          ! 1) fit of H2 + H collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          if (H_H2_flag == "DRF") then
            decal = 1.0_DP
          else
            IF (MOD(ABS_deltaJ,2) == 0) THEN
              decal = 1.0_DP
            ELSE
              decal = 0.3_DP
            ENDIF
          endif

!  JLB - 5 V 2003 - pb avec MM a TRES haute temperature
!         Tdd = decal + Tn * 1.0e-3_dp
          Tdd = decal + min(Tn * 1.0e-3_dp, 5.0e1_dp)
          Tdd2 = Tdd*Tdd

          ! 1.1) non-reactive collisions (always possible)
          !      Warning: reactive collision are mixed in Martin & Mandy,
          !               and are thus included in (badly-named) coeff_H_NR

          gam = rate_H_H2(1,Nup,Nlow)           &
               + rate_H_H2(2,Nup,Nlow) / Tdd     &
               + rate_H_H2(3,Nup,Nlow) / Tdd2    &
               + rate_H_H2(4,Nup,Nlow) * Tdd
          coeff_H_NR = 10.0_DP**gam

          coeff_H_R = 0.0_DP
          ! 1.2) reactive collisions are added (David Flower's rule) if H_H2_flag = DRF
          !      and for Delta J = 2 only if H_H2_flag =  BOTH (complement quantum rates)
          IF (H_H2_flag == "DRF" .AND. Vup /= Vlow) THEN
             IF (MOD(ABS_deltaJ,2) == 0) THEN
                coeff_H_R = 10.0_DP**gam * &
                     EXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
             ELSE
                IF (Jlow > 0) THEN
                   Level1 = index_VJ_H2(Vlow,Jlow-1)
                ELSE
                   Level1 = index_VJ_H2(Vlow,Jlow+1)
                ENDIF
                Level2 = index_VJ_H2(Vlow,Jlow+1)
                IF (Level2 == num_undefined) THEN
                   Level2 = index_VJ_H2(Vlow,Jlow-1)
                ENDIF

                gam1 = rate_H_H2(1,Nup,Level1)         &
                     + rate_H_H2(2,Nup,Level1) / Tdd   &
                     + rate_H_H2(3,Nup,Level1) / Tdd2  &
                     + rate_H_H2(4,Nup,Level1) * Tdd
                gam2 = rate_H_H2(1,Nup,Level2)         &
                     + rate_H_H2(2,Nup,Level2) / Tdd   &
                     + rate_H_H2(3,Nup,Level2) / Tdd2  &
                     + rate_H_H2(4,Nup,Level2) * Tdd
                coeff_H_R = (10.0_DP**gam1 + 10.0_DP**gam2) * 0.5_DP * &
                     EXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
                IF (MOD(Jup,2) == 1) THEN
                   coeff_H_R = coeff_H_R / 3.0_DP
                ENDIF

             ENDIF
          ENDIF
!         IF (H_H2_flag == "BOTH" .AND. Vup /= Vlow  .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow) == .TRUE.) THEN
          IF (H_H2_flag == "BOTH" .AND. Vup /= Vlow  .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow)) THEN
             coeff_H_R = 10.0_DP**gam * &
                     EXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
          ENDIF

          ! 1.3) ortho-para transfer with H from Schofield (only if MM are not included)
          !    + David Flower prescription (only for Delta v=0)
          IF (H_H2_flag == "DRF" &
!       .OR. (H_H2_flag == "BOTH" .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow) == .TRUE.)) THEN
        .OR. (H_H2_flag == "BOTH" .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow))) THEN
            IF (Vup == Vlow .AND. ABS_deltaJ < 3) THEN
               IF (MOD(Jup,2) == 1 .AND. MOD(ABS_deltaJ,2) == 1) THEN
                  coeff_H_R = 8.0D-11 * EXP(-3900.0_DP/Tn) / 3.0_DP
               ELSE
                  coeff_H_R = 8.0D-11 * EXP(-3900.0_DP/Tn)
               ENDIF
            ENDIF
          ENDIF

          ! coeff_H_R is 0.0 if H_H2_flag = MM

          coeff_H = coeff_H_NR + coeff_H_R

          !--------------------------------------------------------------------------
          ! 2) fit of H2 + He collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_He_H2(1,Nup,Nlow)          &
               + rate_He_H2(2,Nup,Nlow) / Tdd   &
               + rate_He_H2(3,Nup,Nlow) / Tdd2  &
               + rate_He_H2(4,Nup,Nlow) * Tdd
          coeff_He = 10.0_DP**gam

          !--------------------------------------------------------------------------
          ! 3) fit of H2 + para-H2 collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_H2_H2(1,Nup,Nlow)          &
               + rate_H2_H2(2,Nup,Nlow) / Tdd  &
               + rate_H2_H2(3,Nup,Nlow) / Tdd2 &
               + rate_H2_H2(4,Nup,Nlow) * Tdd
          coeff_H2 = 10.0_DP**gam

          !------------------------------------------------
          ! 4) ortho-para transfer with ions (Gerlich data)
          !    only in v=0
          !------------------------------------------------
          IF (Vup == 0 .AND. ABS_deltaJ == 1) THEN
             coeff_Hplus = Gerlich(Jup)
          ELSE
             coeff_Hplus = 0.0_DP
          ENDIF

          !---------------------------
          ! 5) Calculate Cij_H2 (s-1)
          !---------------------------
          deexcit_rate= coeff_H2 * Dens_H2 &
               + coeff_He * Dens_He        &
               + coeff_Hplus * Dens_Hplus  &
               + coeff_H * Dens_H
          excit_rate  = deexcit_rate * deexcit_to_excit

          !----------------------------------
          ! 6) Calculate excitation by grains
          !----------------------------------
! Excitation probabilites are fitted by a spline, see f77 , for
! DeltaV > VIN (first impact velocity at which excitation
! probability becomes non zero), VTH = VIN-2km/s.
! Exponential form at DeltaV < VIN, scaled by PRG0 to pass through point at VIN
! Raw excitation rates are in r_raw_GR_H2
! Second derivatives, as computed by SPLINE are in d2r_raw_GR_H2

          dv_kms = ABS_DeltaV  * 1.0D-5
          Vin_H2 = vin_GR_H2(Vup,Jup)
          Vth = Vin_H2 - 2.0_DP

          excit_grain = 0.0_DP
          prob_excit = 0.0_DP
!         if ((Vlow==0) .and. ((Jlow==0) .or. (Jlow==1)) .and. (mod(Jup-Jlow,2)==0)) then
          if ((Vlow==0) .and. (Jlow <= 7) .and. (mod(Jup-Jlow,2)==0)) then
            if (dv_kms >= Vin_H2) then
               prob_excit = splint(vin,r_raw_GR_H2(:,Vup,Jup),d2r_raw_GR_H2(:,Vup,Jup),dv_kms)
               prob_excit = max (prob_excit,0.0d0)

            else if (dv_kms > 1.0e-20) then
!
! 8.02*dv_kms**2 = 0.1*Ts where 3/2kTs = 1/2mH2 dv**2
! extrapolation to velocity difference lower than the first point in the excitation grid
! at Vin = vin_GR_H2(Vup,Jup)

               prob_excit = EXP(delta_E/8.02*(1.0_DP/vin_GR_H2(Vup,Jup)**2 - 1.0_DP/dv_kms**2))
               prob_excit = prob_excit*pgr0_GR_H2(Vup,Jup)

            else

               prob_excit = 0.0_DP

            end if
            excit_grain = prob_excit * Dens_grain * pi * Rsquare_grain * ABS_DeltaV

          else
            excit_grain = 0.0_DP
          end if

          Cij_H2(Nlow,Nup) = Cij_H2(Nlow,Nup) + excit_rate + excit_grain
          Cij_H2(Nup,Nlow) = Cij_H2(Nup,Nlow) + deexcit_rate

       END DO ! end of loop on Nlow
    END DO ! end of loop on Nup

    !--- calculate sum_Cij_H2 ---
    DO Nup=1,NH2_lev
       sum_Cij_H2(Nup) = SUM(Cij_H2(Nup,1:NH2_lev))
    END DO

    !---------------------------------------------------------------
    ! B. Add (de)excitation probabilities from collisions (Cij_H2)
    !    and from radiation (Aij_H2)
    !---------------------------------------------------------------
    Aij_plus_Cij_H2 = Aij_H2 + Cij_H2 ! two dimensions
    DO i=1, NH2_lev ! diagonal terms
       Aij_plus_Cij_H2(i,i)=Aij_plus_Cij_H2(i,i) - sum_Aij_H2(i) - sum_Cij_H2(i)
    END DO

    !----------------------------------------------------------
    ! C. Calculate evolution terms for populations in cm-3.s-1
    !    = matrix multiplication of population density and
    !    (de)excitation probabilities.
    !----------------------------------------------------------
    YN_rovib_H2(1:NH2_lev) = &
         MATMUL(H2_lev(1:NH2_lev)%density, &
         Aij_plus_Cij_H2(1:NH2_lev,1:NH2_lev))


  END SUBROUTINE EVOLUTION_H2



  SUBROUTINE COMPUTE_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    SOURCE
    ! purpose :
    !    Calculate the source terms for H2 levels (cm-3.s-1), H2 internal energy
    !    (erg.cm-3.s-1), the H2 lines (erg.cm-3.s-1) and the corresponding
    !    cooling rate (erg.cm-3.s-1).
    ! subroutine/function needed :
    !    EVOLUTION_H2
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2, H2_Energy, H2_lines%emiss, cooling_H2
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : kB
    USE MODULE_PHYS_VAR, ONLY : NH2_lev
    USE MODULE_DEBUG_JLB
    IMPLICIT NONE

    !----------------------------------------------------------
    ! In EVOLUTION_H2, the source terms YN_rovib_H2 (cm-3.s-1)
    ! is calculated, with treatment of radiation and collision
    !----------------------------------------------------------
    CALL EVOLUTION_H2

    !-----------------------------------------
    ! Calculate the emissivity (erg.cm-3.s-1)
    ! of each H2 quadrupolar line.
    ! result in H2_lines%emiss
    !-----------------------------------------
    ! first in K.cm-3.s-1
    H2_lines(1:NH2_lines)%emiss = &
         H2_lines(1:NH2_lines)%DeltaE * H2_lines(1:NH2_lines)%Aij * &
         H2_lev(H2_lines(1:NH2_lines)%Nup)%density
    ! change into erg.cm-3.s-1
    H2_lines%emiss=H2_lines%emiss*kB

    !-------------------------------------------------------------
    ! cooling rate cooling_H2 (erg.cm-3.s-1) = sum of emissivities
    !-------------------------------------------------------------
    cooling_H2=SUM(DBLE(H2_lines(1:NH2_lines)%emiss))

    !---------------------------------------------------------------
    ! source term for internal energy of H2 H2_energy (erg.cm-3.s-1)
    !---------------------------------------------------------------
    H2_energy = kB * SUM(YN_rovib_H2(1:NH2_lev) * H2_lev(1:NH2_lev)%energy)

  END SUBROUTINE COMPUTE_H2

END MODULE MODULE_H2
