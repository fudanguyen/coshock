MODULE MODULE_DUST_TREATMENT
   !*****************************************************************************
   !** The module 'MODULE_DUST' contains variables and subroutines             **
   !** related to the dust properties :                                        **
   !**     * the type RF_TYPE                                                  **
   !**     * the parameters for the radiation field grid                       **
   !**     * RF_INIT (radiation field initialization subroutine)               **
   !**     * RF_CREATE (create an RF_TYPE variable)                            **
   !**     * RF_FREE (free an RF_TYPE variable)                                **
   !**     * RF_LOCATE (locate a value in an RF_TYPE variable)                 **
   !**     * RF_EXTEND (extend an RF_TYPE variable to a new value)             **
   !**     * RF_REVERS (reorder indexes of an rf_type variable)                **
   !**     * auxiliary subroutines (integral, ...)                             **
   !** Most of the routines & structures are taken from the Meudon PDR code    **
   !*****************************************************************************
 
   IMPLICIT NONE
   INCLUDE "precision.f90"

   ! --------------------------------------------------------------------------------------
   ! Draine's files give log grids for Q_abs, Q_scat, g.
   ! Change the following PARAMETERs if you change file.
   ! --------------------------------------------------------------------------------------
   CHARACTER(LEN=*), PARAMETER                     :: name_file_qcoef = 'Mix_Qcoef_BG_10_2013.dat'
   CHARACTER(LEN=*), PARAMETER                     :: name_file_hcapa = 'Mix_hcapa_BG_10_2013.dat'
   CHARACTER(LEN=*), PARAMETER                     :: name_file_dielc = 'Draine1985-Dielec-func.dat'

   ! --------------------------------------------------------------------------------------
   ! Grains data base. Read in READ_DRAINE_FILES
   ! --------------------------------------------------------------------------------------
   INTEGER                                         :: np_r                 ! Number of radius pts READ
   INTEGER                                         :: np_w                 ! Number of wlgth pts READ
   REAL (KIND=dp)                                  :: r_min, r_max
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: radgr
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: wavgrd
   INTEGER                                         :: jwlmin_gr, jwlmax_gr ! range of wl used in wavgrd

   ! --------------------------------------------------------------------------------------
   ! Optical properties
   ! --------------------------------------------------------------------------------------
   REAL (KIND=dp),   DIMENSION (:,:), ALLOCATABLE  :: Qabs_Gra, Qsca_Gra, g_Gra
   REAL (KIND=dp),   DIMENSION (:,:), ALLOCATABLE  :: Qabs_Sil, Qsca_Sil, g_Sil
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: Q_abs, Q_sca, Q_emi
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: kappa_dst

   ! --------------------------------------------------------------------------------------
   ! Heat capacities: from a mix of dustem's data over log grids of dust temperature
   ! --------------------------------------------------------------------------------------
   INTEGER, PUBLIC                                 :: np_t                 ! Number of temperature pts READ
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: tempgrd
   REAL (KIND=dp),   DIMENSION (:,:), ALLOCATABLE  :: Cheat_Gra,Cheat_Sil
   REAL (KIND=dp),   DIMENSION (:),   ALLOCATABLE  :: Cheat

   ! --------------------------------------------------------------------------------------
   ! Positions of the photometric bands
   ! --------------------------------------------------------------------------------------
   REAL (KIND=dp),   PARAMETER :: B_band_wl = 4450.0_dp
   REAL (KIND=dp),   PARAMETER :: V_band_wl = 5510.0_dp
   INTEGER                     :: ilocb, ilocv

   PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



SUBROUTINE DEALLOCATE_DUST
   !---------------------------------------------------------------------------
   ! called by :
   !     MHD_VODE
   ! purpose :
   !     Deallocate all the tables allocated to store the dust properties
   ! subroutine/function needed :
   ! input variables :
   ! output variables :
   ! results :
   !---------------------------------------------------------------------------
   IMPLICIT none

   DEALLOCATE ( radgr )
   DEALLOCATE ( wavgrd )
   DEALLOCATE ( Qabs_Gra )
   DEALLOCATE ( Qsca_Gra )
   DEALLOCATE ( g_Gra )
   DEALLOCATE ( Qabs_Sil )
   DEALLOCATE ( Qsca_Sil )
   DEALLOCATE ( g_Sil )
   DEALLOCATE ( Q_abs )
   DEALLOCATE ( Q_sca )
   DEALLOCATE ( Q_emi )
   DEALLOCATE ( kappa_dst )
   DEALLOCATE ( tempgrd )
   DEALLOCATE ( Cheat_Gra )
   DEALLOCATE ( Cheat_Sil )
   DEALLOCATE ( Cheat )

END SUBROUTINE DEALLOCATE_DUST



SUBROUTINE READ_DRAINE_FILES
   !---------------------------------------------------------------------------
   ! called by :
   !     INITIALIZE
   ! purpose :
   !     Read the Draine files containing the grain absorption coefficients
   !     and heat capacities depending on the size and composition of dust
   !        - allocation of tables
   !        - unit conversions
   ! subroutine/function needed :
   ! input variables :
   ! output variables :
   ! results :
   !     radgr, wavgrd, tempgrd
   !     Qabs_Gra, Qsca_Gra, g_Gra
   !     Qabs_Sil, Qsca_Sil, g_Sil
   !     Cheat_Gra,Cheat_Sil
   !---------------------------------------------------------------------------
   USE MODULE_TECHCONFIG, ONLY : fichier, data_dir, grains_dir
   USE MODULE_CONSTANTS

   IMPLICIT NONE

   INTEGER                                    :: i
   INTEGER                                    :: j
   REAL (KIND=dp)                             :: toto
   REAL (KIND=dp)                             :: dum
   REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: radgr_hc

   !======================================================================================
   !             READING  GRAPHITE and SILICATE IN FILE:  Mix_Qcoef_BG_10_2013
   !======================================================================================
   fichier = TRIM(data_dir)//TRIM(grains_dir)//TRIM(name_file_qcoef)

   ! -------------------------------------------------------------------------------------
   ! Read and test the nb of points (in dust radius and wavelength) contained in datafile
   ! Allocate all the related tables
   ! -------------------------------------------------------------------------------------
   OPEN (iwrtmp, file=fichier, status='old')
   READ (iwrtmp,'(3/,i4,2e10.3)') np_r, r_min, r_max
   READ (iwrtmp,*)                np_w
   ALLOCATE ( radgr(1:np_r) )
   ALLOCATE ( wavgrd(0:np_w+1) )
   ALLOCATE ( Qabs_Gra(1:np_r,0:np_w+1) )
   ALLOCATE ( Qsca_Gra(1:np_r,0:np_w+1) )
   ALLOCATE ( g_Gra(1:np_r,0:np_w+1) )
   ALLOCATE ( Qabs_Sil(1:np_r,0:np_w+1) )
   ALLOCATE ( Qsca_Sil(1:np_r,0:np_w+1) )
   ALLOCATE ( g_Sil(1:np_r,0:np_w+1) )
   ALLOCATE ( Q_abs(0:np_w+1) )
   ALLOCATE ( Q_sca(0:np_w+1) )
   ALLOCATE ( Q_emi(0:np_w+1) )
   ALLOCATE ( kappa_dst(0:np_w+1) )

   ! -------------------------------------------------------------------------------------
   ! Read dust properties in the datafile (graphite and silicate)
   ! -------------------------------------------------------------------------------------
   DO i = 1, np_r
      READ (iwrtmp,*) radgr(i)
      READ (iwrtmp,*)
      DO j = np_w, 1, -1
         READ (iwrtmp,*) wavgrd(j), Qabs_Gra(i,j), Qsca_Gra(i,j), g_Gra(i,j), &
                       & toto,      Qabs_Sil(i,j), Qsca_Sil(i,j), g_Sil(i,j)
         IF (toto /= wavgrd(j)) THEN
            PRINT *, " pb lecture!"
         ENDIF
      ENDDO
   ENDDO
   CLOSE(iwrtmp)
   radgr(1:np_r)    = radgr(1:np_r) * 1.0e-4_dp ! Convert from micrometer to cm
   wavgrd(0:np_w+1) = wavgrd(0:np_w+1) * 1.0e4_dp ! Convert from micrometer to Angstrom

   ! -------------------------------------------------------------------------------------
   ! Locate the B and V band in the wavelength grid
   ! -------------------------------------------------------------------------------------
   dum = ABS(B_band_wl - wavgrd(1)) / B_band_wl
   DO i = 2, np_w
      IF( ABS(B_band_wl - wavgrd(i)) / B_band_wl < dum ) THEN
         ilocb = i
         dum = ABS(B_band_wl - wavgrd(i)) / B_band_wl
      ENDIF
   ENDDO
   dum = ABS(V_band_wl - wavgrd(1)) / V_band_wl
   DO i = 2, np_w
      IF( ABS(V_band_wl - wavgrd(i)) / V_band_wl < dum ) THEN
         ilocv = i
         dum = ABS(V_band_wl - wavgrd(i)) / V_band_wl
      ENDIF
   ENDDO
   WRITE(*,*) B_band_wl, wavgrd(ilocb)
   WRITE(*,*) V_band_wl, wavgrd(ilocv)

   !======================================================================================
   !             READING HEAT CAPACITIES IN FILE: Mix_hcapa_BG_10_2013
   !======================================================================================
   fichier = TRIM(data_dir)//TRIM(grains_dir)//TRIM(name_file_hcapa)

   ! -------------------------------------------------------------------------------------
   ! Read and test the nb of points (in temperature) contained in datafile
   ! Allocate all the related tables
   ! -------------------------------------------------------------------------------------
   OPEN (iwrtmp, file=fichier, status='old')
   READ (iwrtmp,*) np_t
   ALLOCATE ( radgr_hc(1:np_r) )
   ALLOCATE ( tempgrd(1:np_t) )
   ALLOCATE ( Cheat_Gra(1:np_r,1:np_t) )
   ALLOCATE ( Cheat_Sil(1:np_r,1:np_t) )
   ALLOCATE ( Cheat(1:np_t) )

   ! -------------------------------------------------------------------------------------
   ! Read heat capacities in the datafile  (graphite and silicate)
   ! -------------------------------------------------------------------------------------
   DO i = 1, np_r
      READ (iwrtmp,*) radgr_hc(i)
      radgr_hc(i) = radgr_hc(i) * 1.0e-4_dp ! Convert from micrometer to cm
      IF(ABS(radgr_hc(i) - radgr(i))/radgr(i) > 1.0e-10_dp) THEN
         PRINT *, "grains radii in Mix_hcapa_BG_10_2013.dat differ from those in Mix_Qcoef_BG_10_2013.dat"
         STOP
      ENDIF
      READ (iwrtmp,*)
      DO j = 1, np_t
         READ (iwrtmp,*) tempgrd(j), Cheat_Gra(i,j), Cheat_Sil(i,j)
      ENDDO
      READ(iwrtmp,*)
      ! multiply by the grain volume
      Cheat_Gra(i,1:np_t) = Cheat_Gra(i,1:np_t) * 4.0_dp/3.0_dp * pi * radgr(i)**3.0_dp
      Cheat_Sil(i,1:np_t) = Cheat_Sil(i,1:np_t) * 4.0_dp/3.0_dp * pi * radgr(i)**3.0_dp
   ENDDO
   CLOSE(iwrtmp)

   DEALLOCATE ( radgr_hc )

END SUBROUTINE READ_DRAINE_FILES



REAL(KIND=dp)  FUNCTION Q_inter(Qx,x,a,i0)
   !---------------------------------------------------------------------------
   ! called by :
   !      COMPUTE_QCOEF
   ! purpose :
   !     Interpolate grain coefficients between two sizes
   ! subroutine/function needed :
   ! input variables :
   !     i0, i0+1 -> lower and upper sizes indices in input files
   !     x        -> grains sizes in input files
   !     a        -> grain size
   !     Q        -> Q coef in input files
   ! output variables :
   ! results :
   !---------------------------------------------------------------------------

   IMPLICIT NONE

   REAL(KIND=dp), DIMENSION(:), INTENT (IN) :: Qx
   REAL(KIND=dp), DIMENSION(:), INTENT (IN) :: x
   REAL(KIND=dp),               INTENT (IN) :: a
   INTEGER,                     INTENT (IN) :: i0

   Q_inter = Qx(i0) + (Qx(i0+1)-Qx(i0)) * ((a-x(i0)) / (x(i0+1)-x(i0)))

END FUNCTION Q_inter



SUBROUTINE COMPUTE_QCOEF
   !---------------------------------------------------------------------------
   ! called by :
   !     INITIALIZE
   ! purpose :
   !     Interpolate the absorption and scattering coefficients assuming a
   !     grain of size r_grm (root of the mean square of mantle radius)
   !     We adopt graphite absorption coefficient as a test
   ! subroutine/function needed :
   ! input variables :
   ! output variables :
   ! results :
   !      Q_abs, Q_sca
   !      kappa_dst
   !      conv_coldens (depending on F_invAv)
   !      inv_Av_fac   (depending on F_invAv)
   !---------------------------------------------------------------------------
   USE MODULE_GRAINS,    ONLY : r_grm, rsq_grm, dens_grc, rv, cdunit
   USE MODULE_CONSTANTS, ONLY : pi
   USE MODULE_PHYS_VAR,  ONLY : nH, F_invAv, inv_Av_fac, conv_coldens

   IMPLICIT none

   INTEGER       :: i0
   INTEGER       :: i

   ! ---------------------------------------
   ! Find the closest size of grain in 
   ! Draine's file smallest than r_grm
   ! ---------------------------------------
   i0 = 1
   DO
      IF (radgr(i0+1) > r_grm .AND. radgr(i0) <= r_grm) THEN
         EXIT
      ENDIF
      i0 = i0 + 1
   ENDDO

   DO i = 1, np_w
      Q_abs(i) = Q_inter(Qabs_Gra(1:np_r,i),radgr,r_grm,i0)
      Q_sca(i) = Q_inter(Qsca_Gra(1:np_r,i),radgr,r_grm,i0)
   ENDDO

   ! ---------------------------------------
   ! Take only the absorption coefficients 
   ! to compute kappa_dust (approximation)
   ! ---------------------------------------
   kappa_dst(1:np_w) = pi * rsq_grm * Q_abs(1:np_w) * dens_grc

   ! ---------------------------------------
   ! Compute the rv and cdunit obtained with 
   ! grain properties
   ! ---------------------------------------
   rv     = kappa_dst(ilocv) / (kappa_dst(ilocb)-kappa_dst(ilocv))
   cdunit = 1.0_dp / ( 2.5_dp * LOG10(exp(1.0)) * (kappa_dst(ilocb)-kappa_dst(ilocv)) ) * nH

   ! ---------------------------------------
   ! Scale the grain properties coefficients
   ! to reproduce the NH/Av chosen in input
   ! - constant scaling for all lambdas
   ! - recompute cdunit
   ! Or not
   ! - recompute conv_coldens and inv_Av_fac
   !   based on grain coefficients
   ! ---------------------------------------
   IF(F_invAv == 1) THEN
      kappa_dst(1:np_w) = kappa_dst(1:np_w) * inv_Av_fac / (rv / cdunit)
      Q_abs(1:np_w)     = Q_abs(1:np_w) * inv_Av_fac / (rv / cdunit)
      cdunit = 1.0 / ( 2.5 * LOG10(exp(1.0)) * (kappa_dst(ilocb)-kappa_dst(ilocv)) ) * nH
   ELSE
      conv_coldens = cdunit / rv
      inv_Av_fac   = 1.0_dp / conv_coldens
   ENDIF

END SUBROUTINE COMPUTE_QCOEF


END MODULE MODULE_DUST_TREATMENT

