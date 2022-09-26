MODULE MODULE_DRIVE

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-------------------------------------------------------
  ! values saved before call to DRIVE
  !-------------------------------------------------------
  REAL (KIND=DP)            :: Vn_old          = 0.0_DP
  REAL (KIND=DP)            :: Vi_old          = 0.0_DP
  REAL (KIND=DP)            :: dVn_old         = 0.0_DP
  REAL (KIND=DP)            :: dVi_old         = 0.0_DP
  REAL (KIND=DP)            :: old_Tn          = 0.0_DP
  REAL (KIND=DP)            :: Tinit           = 0.0_DP ! Detects whether Tn decreases.
  LOGICAL                   :: increase        = .false.
  LOGICAL                   :: maximum_reached = .false.

  !-------------------------------------------------------
  ! initial time scales
  !-------------------------------------------------------
  REAL (KIND=DP)            :: time_scale
  REAL (KIND=DP)            :: tiny = 1.0e-40_DP
  INTEGER                   :: i_time

  !-------------------------------------------------------
  ! Computation timescales
  !-------------------------------------------------------
  REAL (KIND=DP)            :: time1
  REAL (KIND=DP)            :: time2

  !-------------------------------------------------------
  ! Variables for CJ-type shocks
  !-------------------------------------------------------
  ! find jump conditions
  LOGICAL                   :: cj_active   = .false.
  LOGICAL                   :: cs_active   = .false.
  REAL (KIND=DP)            :: t_min_jump
  REAL (KIND=DP)            :: t_max_jump
  REAL (KIND=DP)            :: t_jump
  INTEGER                   :: i_min_jump
  INTEGER                   :: i_max_jump
  INTEGER                   :: i_jump
  ! adjust trajectory
  LOGICAL                   :: adjust_traj  = .false.
  LOGICAL                   :: adjust_start = .false.
  LOGICAL                   :: test_sonic   = .false.
  LOGICAL                   :: sonic_point  = .false.
  LOGICAL                   :: found_high   = .false.
  LOGICAL                   :: found_down   = .false.
  REAL (KIND=DP), PARAMETER :: eps_traj = 1.e-2_dp
  REAL (KIND=DP)            :: t_limit
  INTEGER                   :: i_limit
  INTEGER                   :: i_adjst                ! QUESTION - REMOVE i_adjst ? -> ONLY DEAL WITH i_limit
  INTEGER                   :: i_down
  INTEGER                   :: nadj = 0
  REAL (KIND=DP)            :: alpha
  REAL (KIND=DP)            :: angle_traj = 0.0_dp
  REAL (KIND=DP)            :: dum
  REAL (KIND=DP)            :: dum1
  REAL (KIND=DP)            :: dum2
  REAL (KIND=DP)            :: coef1
  REAL (KIND=DP)            :: coef2
  REAL (KIND=DP)            :: Vn_high
  REAL (KIND=DP)            :: Vn_down
  REAL (KIND=DP)            :: Vi_high
  REAL (KIND=DP)            :: Vi_down
  REAL (KIND=DP)            :: step_cross = 0.0_dp
  REAL (KIND=DP)            :: sonic_jump = 0.0_dp
  ! recoupling of fluids
  LOGICAL                   :: fluid_decoup = .false.
  LOGICAL                   :: fluid_recoup = .false.
  REAL (KIND=DP)            :: decoup_dist  = 0.0_dp
  REAL (KIND=DP)            :: decoup_strg  = 0.0_dp
  REAL (KIND=DP)            :: corr_flux_dum= 0.0_dp

  INTEGER                   :: countersave

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS


  SUBROUTINE RESET_INTEGRATION(is, again)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    reset integration for the next slab
     ! subroutine/function needed :
     ! input variables :
     !    slab - integer of current slabl
     !    again - boolean
     ! ouput variables :
     !    again - boolean : false if the shock code must stop
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_GRAINS
     USE MODULE_H2
     USE MODULE_CHEMICAL_SPECIES
     USE MODULE_CHEM_REACT
     USE MODULE_DUST_TREATMENT
     USE MODULE_RADIATION
     USE MODULE_EVOLUTION
     USE MODULE_VAR_VODE

     IMPLICIT none

     INTEGER,        INTENT(in)    :: is
     LOGICAL,        INTENT(inout) :: again
     INTEGER                       :: i

     REAL (KIND=dp)                :: dumy
     REAL (KIND=dp)                :: bt_s_bd
     REAL (KIND=dp)                :: cdinf

     LOGICAL,        SAVE          :: first=.true.

     ! -------------------------------------------------
     ! set up the origin of distances
     ! -------------------------------------------------
     IF (first) THEN
        first    = .false.
        distance = posgrd(is-1)
     ENDIF

     ! -------------------------------------------------
     ! set the type of shock to run
     ! -------------------------------------------------
     IF ( posgrd(is-1) < 0 ) THEN
        shock_type = "S1"
     ELSE
        shock_type = shock_type_lock
     ENDIF
     IF      ( shock_type(1:1) == "S" .OR.&
               shock_type(1:1) == "P" ) THEN
        Nfluids    = 1
        viscosity  = .FALSE.
     ELSE IF ( shock_type(1:1) == "J" ) THEN
        Nfluids    = 1
        viscosity  = viscosity_lock
     ELSE IF ( shock_type(1:1) == "C" ) THEN
        Nfluids    = 3
        viscosity  = .FALSE.
     ENDIF

     ! -------------------------------------------------
     ! IF F_COUP_RAD = 2
     !    Set the radiation field to its value at the
     !    current slab and reset the FGK coef and
     !    photodestruction rates with this radiation
     ! -------------------------------------------------
     IF (F_COUP_RAD == 2) THEN
        CALL SET_RAD_GRID(is)
        CALL FGKCOEF

        DO i = 1, Ncross
           phdest_rate(i)    = SECT_INT(phdest(i))
           phheating_rate(i) = HEATING_PHOTOCHEM_SECT_INT(phdest(i))
          ! -- debug --
          ! WRITE(*,*) speci(phdest(i)%ispe)%name, phheating_rate(i), phdest_rate(i)
        ENDDO
     ENDIF

     ! -------------------------------------------------
     ! Reset dvode arguments and call diffun
     ! -------------------------------------------------
     again         = .TRUE.
     distance_old  = 0.0_dp                  ! distance old
     T0_V          = 0.0_dp                  ! distance : value of the integration variable
     Tout_V        = posgrd(is)-distance     ! step to reach

     Hnext         = 0.0d0
     H0_V          = XLL * 1.d-5             ! step length - not really useful
     MF_V          = 22                      ! chosen numerical procedure
     Itask_V       = 1                       ! Which Task DVODE should do
     Istate_V      = 1                       ! How was it done (2 on normal return, 1 to initialise dvode)
     count_istat4  = 0

     CALL diffun(d_v_var,T0_v,v_lvariab,v_dvariab)
     dVn = v_dvariab(iv_Vn)
     dVi = v_dvariab(iv_Vi)

     ! -------------------------------------------------
     ! Info on chemical rates
     ! -------------------------------------------------
     ! PRINT *,'Species name,   density,    Chemical rate (1/s):'
     ! DO i = 1, Nspec
     !    PRINT *, speci(i)%name, speci(i)%density, YN(i)
     ! ENDDO
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Print info on initial length scales
     !-------------------------------------------------------
     PRINT *
     time_scale =  1.0_DP / MAXVAL( abs(v_dvariab(1:d_v_var) ) + tiny )
     i_time = MAXLOC( abs(v_dvariab), dim = 1 )
     PRINT *,' Length scale estimates AT START:'
     PRINT "(' lenght_scale=',1pe13.3,' cm')", time_scale
     PRINT "(' Max_var number i=',i4)", i_time
     IF (i_time.ge.bv_H2_lev.and.i_time.le.ev_H2_lev) THEN
        PRINT *,' this corresponds to H2 level number ',i_time-bv_H2_lev+1
     ENDIF
     IF (i_time.ge.bv_speci.and.i_time.le.ev_speci) THEN
        PRINT *,' this corresponds to species ',speci(i_time-bv_speci+1)%name
     ENDIF
     PRINT "(' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
          exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
     PRINT *
     rwork_V(5) = MIN(1e-4*time_scale,1e-4_dp) ! H0_V           ! step size on first step
     rwork_V(6) = Tout_V                       ! 0 !1.0D17      ! HMAX
     rwork_V(7) = MIN(1e-8*time_scale,1e-8_dp) ! H0_V * 1.0D-10 ! HMIN

  END SUBROUTINE RESET_INTEGRATION



  SUBROUTINE TEST_CONTINUE(again, message)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if the shock code must stop
     ! subroutine/function needed :
     ! input variables :
     !    again - boolean
     ! ouput variables :
     !    again - boolean : false if the shock code must stop
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_PHYS_VAR,          ONLY : Nstep_max, counter, stop_code, timeN,&
                                          coldens_h2, N_H2_0, coldens_co, N_CO_0,&
                                          shock_type, distance, Av, inv_Av_fac, Tn
     USE MODULE_VAR_VODE,          ONLY : duration_max, length_max

     IMPLICIT none

     LOGICAL,           INTENT(inout) :: again
     CHARACTER(LEN=80), INTENT(out)   :: message

     !----------------------------------------------------
     ! STOP integration if :
     !----------------------------------------------------

     ! (1) Maximal number of steps has been reached
     IF (counter == Nstep_max) THEN
        again = .FALSE.
        WRITE (message, '("Maximal number of steps has been reached.")' )
        stop_code = 1
     END IF

     ! (2) Maximal evolution time has been reached
     IF (timeN >= duration_max) THEN
        again = .FALSE.
        WRITE (message,'("Maximal evolution time has been reached. NH2=",1pe13.3," +",1pe13.3)') coldens_h2,N_H2_0
        stop_code = 2
     END IF

     ! (3) Maximal shock length has been reached
     IF (distance >= length_max) THEN
        again = .FALSE.
        PRINT *,'length_max',length_max
        PRINT *,'distance',distance
        WRITE (message,'("Maximal shock length has been reached. NH2=",1pe13.3,A1,1pe13.3)') coldens_h2,'+',N_H2_0
        stop_code = 3
     END IF

     ! (4) drift velocity < DeltaVmin
     ! IF (ABS_DeltaV < DeltaVmin * 1.0d-1) THEN
     ! ! IF ((Vn - Vi) < 10.0_DP) THEN
     !    again = .FALSE.
     !    message_end = "Equilibrium has been reached."
     !    stop_code = 4
     ! END IF

     ! (5) we reach a sonic point
     ! IF (Vn < (Vsound + 1.0D+1)) THEN
     !    again = .FALSE.
     !    message_end = "Sonic Point!"
     !    stop_code = 5
     ! END IF

     ! (6) Temperature is below 50.0 K
     ! IF (counter > 1000 .AND. Tn < 50.0_DP) THEN
     !    again = .FALSE.
     !    message_end = "temperature lower than 50 K"
     !    stop_code = 6
     ! END IF

     ! (6) Temperature decreases below initial T
     IF (old_Tn ==0 ) THEN ! detect first step
        Tinit  = Tn
        old_Tn = Tn
     ENDIF
     IF (Tn > old_Tn) THEN
        increase=.true.
     ENDIF
     IF (Tn < old_Tn .and. increase) THEN
        maximum_reached=.true.
        print *,'Reach maximum !'
        increase = .false.
     ENDIF
     IF ( (.false.) .and. (shock_type(1:1) /= 'S') .and. (shock_type(1:1) /= 'P') ) THEN
        again = .FALSE.
        WRITE (message, '("temperature decreases below 0.9*initial T.")')
        PRINT *, "T=", Tinit, 'Tn=', Tn
        stop_code = 6
     ELSE
        old_Tn = Tn
     ENDIF
 
     ! (7) Numerical instabilities prevent conservations of ions
     !!!!!!! IF (abs((DensityI-SUM(v_variab(bv_ion:ev_ion)))/DensityI) > 1.0e-2_DP) THEN
     !!!!!!!    again = .FALSE.
     !!!!!!!    message_end = "Ions are NOT conserved"
     !!!!!!!    stop_code = 7
     !!!!!!! ENDIF

     ! (8) Av is too small and we are integrating "backward"
     IF (Av<0_DP.and.inv_Av_fac<0_DP) THEN
        again = .false.
        WRITE (message, '("Av<0, that s a bit weird...")')
        stop_code = 8
     ENDIF

     ! (9) NH2 is too small as a result of backward integration
     IF ( (inv_Av_fac<0_DP ).and.(coldens_h2+N_H2_0<1e10_DP) ) THEN
        again = .false.
        PRINT *, 'NH2+N_H2_0=', coldens_h2+N_H2_0
        WRITE (message,'("NH2 residual < 1e10, we will soon reach NH2=0. Relative age=",1pe13.3)') timeN/duration_max
        stop_code = 9
     ENDIF

     ! (10) NCO is too small as a result of backward integration
     IF ( (inv_Av_fac<0_DP ).and.(coldens_co+N_CO_0<1e05_DP) ) THEN
        again = .false.
        PRINT *, 'NCO+N_CO_0=', coldens_co+N_CO_0
        WRITE (message,'("NCO residual < 1e05, we will soon reach NCO=0. Relative age=",1pe13.3)') timeN/duration_max
        stop_code = 10
     ENDIF

  END SUBROUTINE TEST_CONTINUE


  SUBROUTINE TEST_COUPLING(restart)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if ions and neutrals decouple or recouple or if a non-stationary
     !    CJ type shock needs to be run
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     !    restart - boolean : send code to restart the integration
     !                         with a J-type shock
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_RADIATION
     USE MODULE_VAR_VODE
     USE MODULE_H2,            ONLY : H2_lev
     USE MODULE_CO,            ONLY : CO_lev
     USE MODULE_EVOLUTION
     USE MODULE_ENERGETICS,    ONLY : mag_flux_corr
     USE MODULE_PROFIL_TABLES
     USE MODULE_SWITCH_CJ

     IMPLICIT none

     LOGICAL, INTENT(out) :: restart
     INTEGER              :: i

     restart = .false.

     !====================================================================================
     ! Test fluid decoupling strength and fluid recoupling hypothesis
     !====================================================================================
     IF ( shock_type == 'C' ) THEN
        IF     ( .NOT.fluid_decoup ) THEN
           CALL ESTIMATE_DECOUPLING( counter-1, eps_traj, decoup_strg )
           IF ( ABS(decoup_strg) / Vs_cm > eps_traj ) THEN
              fluid_decoup = .true.
              decoup_dist  = distance
           ENDIF
        ELSE IF ( ( ( distance > decoup_dist * (1.0_dp + eps_traj) ).AND.( .NOT.cj_active ).AND.( Vn > Vsound ) ).OR.&
                  ( ( distance > sonic_jump  * (1.0_dp + eps_traj) ).AND.( sonic_point    ) ) ) THEN
           CALL ESTIMATE_DECOUPLING( counter-1, eps_traj, decoup_strg )
           IF ( ABS(decoup_strg) / Vs_cm < eps_traj ) THEN
              fluid_recoup = .true.
           ENDIF
        ENDIF

        IF ( ( timeI >= timeJ .OR. fluid_recoup ).AND.&
             ( .NOT.cj_active .OR. sonic_point  ) ) THEN
           ! -------------------------------------
           ! set up a J-type shock
           ! -------------------------------------
           fluid_recoup = .true.

           ! prevent from switching back to C shocks
           shock_type_lock = "J"

           ! -------------------------------------
           ! reinitialize physical variables
           ! -------------------------------------
           CALL REINIT_PHYS_VARIABLES(traj_main(counter-1))
           dVn_old = dVn
           dVi_old = dVi
           IF (F_AV /= 0) THEN
              DO i = 1, nangle
                 irf(i,1:nwlg) = irf_old(i,1:nwlg)
              ENDDO
           ENDIF
           DO i = 1, NH2_lev
              H2_lev(i)%cd_l = H2_lev(i)%cd_l_old
           ENDDO
           DO i = 1, NCO_lev
              CO_lev(i)%cd_l = CO_lev(i)%cd_l_old
           ENDDO

           ! -------------------------------------
           ! set neutral velocity to ion velocity
           ! assuming neutrals recouple with ions
           ! and not the opposite
           ! -------------------------------------
           ! CALL FALSE_JUMP(Vi)
           ! Vn = Vi
           ! v_variab(iv_Vn) = Vn
           ! v_lvariab(iv_Vn) = LOG(Vn)
           ! -------------------------------------
           ! set ion velocity to neutral velocity
           ! assuming ions recouple with neutrals
           ! and not the opposite
           ! -------------------------------------
           corr_flux_dum = 1.0_DP / Vi
           Vi = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum

           ! -------------------------------------
           ! restart integration
           ! -------------------------------------
           counter = counter - 1
           restart = .true.
        ENDIF
     ENDIF

  END SUBROUTINE TEST_COUPLING


  SUBROUTINE TEST_STAT_CJ(restart)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if ions and neutrals decouple or recouple or if a non-stationary
     !    CJ type shock needs to be run
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     !    restart - boolean : send code to restart the integration
     !                         with a J-type shock
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_RADIATION
     USE MODULE_VAR_VODE
     USE MODULE_EVOLUTION
     USE MODULE_ENERGETICS,    ONLY : mag_flux_corr
     USE MODULE_PROFIL_TABLES
     USE MODULE_SWITCH_CJ

     IMPLICIT none

     LOGICAL, INTENT(out) :: restart
     INTEGER              :: i

     restart = .false.

     !====================================================================================
     ! Switch to CJ-type or C*-type
     !   The procedure to reconstruct the trajectory of these shocks is used if
     !     - dVn changes sign or
     !     - dVi changes sign and we are already in a CJ or C* shock, or
     !     - Vn smaller than Vi and we are already in a CJ or C* shock, or
     !   Note : the procedure is not built for integrating Av
     !====================================================================================
     IF ( ( shock_type == 'C'                  ).AND.&
          ( .NOT.sonic_point                   ).AND.&
          ( ( Vn - Vsound ) / Vs_cm < eps_traj ).AND.&
          ( ( dVn*dVn_old < 0.0_dp ).OR.( dVi*dVi_old < 0.0_dp ) ) ) THEN

        IF ( F_AV == 1 ) THEN
           WRITE(*,*) "This procedure doesn't work if we integrate Av"
           WRITE(*,*) "-> code stops"
           STOP
        ENDIF

        ! -------------------------------------------
        ! Find minimum and maximum times to jump
        ! -------------------------------------------
        IF ( .NOT.cj_active .AND. .NOT.cs_active ) THEN
           ! ----------------------------------------
           ! save iteration number of original trajec
           ! ----------------------------------------
           counter_main = counter - 1
           ! ----------------------------------------
           ! Find the potential jump conditions
           ! ----------------------------------------
           CALL FIND_JUMP_RANGE(t_min_jump, t_max_jump, i_min_jump, i_max_jump, cs_active)
           cj_active = .true.
           ! ----------------------------------------
           ! Shortcut jump if C* type shock
           ! Useful to debug
           ! ----------------------------------------
           ! IF ( cs_active ) THEN
           !    i_min_jump = i_max_jump - 1
           ! ENDIF
        ELSE
           ! ----------------------------------------
           ! Dichotomie - converge on jumping time
           ! ----------------------------------------
           IF( dVn*dVn_old < 0.0_dp ) THEN
              t_max_jump = t_jump
              IF( i_max_jump - i_min_jump > 1 ) THEN
                 i_max_jump = i_jump
                 cs_active  = .false.
              ENDIF
              ! -------------------------------------
              ! Save upper trajectory
              ! -------------------------------------
              traj_high(1:i_limit)           = traj_main(1:i_limit)
              traj_high(i_limit+1:counter-1) = traj_curr(i_limit+1:counter-1)
              found_high = .true.
              alpha = alpha - 1.0_dp / 2.0_dp**(nadj+1)
              counter_high = counter - 1
           ELSE
              t_min_jump = t_jump
              IF( i_max_jump - i_min_jump > 1 ) THEN
                 i_min_jump = i_jump
              ENDIF
              ! -------------------------------------
              ! Save lower trajectory
              ! -------------------------------------
              traj_down(1:i_limit)           = traj_main(1:i_limit)
              traj_down(i_limit+1:counter-1) = traj_curr(i_limit+1:counter-1)
              found_down = .true.
              alpha = alpha + 1.0_dp / 2.0_dp**(nadj+1)
              counter_down = counter - 1
           ENDIF
        ENDIF

        ! -------------------------------------------
        ! Test - recouple fluids as soon as possible
        ! if no possible jump conditions
        ! SHOULD BE REMOVED - NOT PHYSICAL
        ! -------------------------------------------
        IF (i_min_jump >= i_max_jump) THEN
           ! ----------------------------------------
           ! set up a J-type shock
           ! ----------------------------------------
           fluid_recoup = .true.

           ! prevent from switching back to C shock
           shock_type_lock = "J"

           !--- Find the jump index where ion and neutral recouple
           CALL FIND_RECOUPLING_JUMP(i_jump)
           ! i_jump = i_min_jump
           
           !--- reinitialize physical variables
           CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
           t_jump = distance
           dVn_old = dVn
           dVi_old = dVi
           i_limit = i_jump
           t_limit = t_jump
           !--- Apply Rankine-Hugoniot relations for an adiabatic jump
           !--- hydrodynamical shock - on neutral fluid only
           WRITE(*,*) "pre jump = ", Vn, Vi, Tn, Ti
           CALL RH_JUMP_HYDRO
           WRITE(*,*) "postjump = ", Vn, Vi, Tn, Ti
           READ(*,*) 
           ! ----------------------------------------
           ! set neutral velocity to ion velocity
           ! assuming neutrals recouple with ions
           ! and not the opposite
           ! ----------------------------------------
           ! Vn = Vi
           ! v_variab(iv_Vn) = Vn
           ! v_lvariab(iv_Vn) = LOG(Vn)
           ! -------------------------------------
           ! set ion velocity to neutral velocity
           ! assuming ions recouple with neutrals
           ! and not the opposite
           ! -------------------------------------
           corr_flux_dum = 1.0_DP / Vi
           Vi = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
        ENDIF

        ! -------------------------------------------
        ! Reinitialize physical variables at i_jump
        ! or interpolate between i_jump and i_jump+1
        ! -------------------------------------------
        IF( .NOT.adjust_traj .AND. .NOT.adjust_start .AND. .NOT.fluid_recoup ) THEN
           IF ( i_max_jump - i_min_jump > 1 ) THEN
              i_jump = ( i_min_jump + i_max_jump ) / 2
              !--- reinitialize physical variables
              CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
              t_jump = distance
              dVn_old = dVn
              dVi_old = dVi
           ELSE
              t_jump = ( t_min_jump + t_max_jump ) / 2.0_dp
              i_jump  = i_min_jump
              CALL INTERP_JUMP_CONDITIONS(t_jump,i_jump)
           ENDIF
           i_limit = i_jump
           t_limit = t_jump
           !--- Apply Rankine-Hugoniot relations for an adiabatic jump
           !--- hydrodynamical shock - on neutral fluid only
           CALL RH_JUMP_HYDRO
           i_adjst = i_limit
        ENDIF

        ! -------------------------------------------
        ! Treatment of C* shocks
        ! -------------------------------------------
        IF ( ( ( ( i_max_jump - i_min_jump == 1 ).AND.( cs_active ) ).OR.( adjust_start ) ).AND.( .NOT.adjust_traj ) ) THEN
           ! ----------------------------------------
           ! scan backward the main trajectory and 
           ! select a point where to apply small 
           ! variations of Vn and Vi
           ! ----------------------------------------
           counter_main_old = counter_main
           IF ( ( .NOT. adjust_start ).OR.&
                ( adjust_start .AND. (.NOT.found_down .OR. .NOT.found_high) .AND. (alpha > 0.9_dp .OR. alpha < 0.1_dp) ) ) THEN
              IF( (counter_main - i_jump) / 5 /= 0 ) THEN
                 counter_main = counter_main - (counter_main - i_jump) / 5
              ELSE
                 counter_main = counter_main - 1
              ENDIF
           ENDIF
           i_limit = counter_main
           i_adjst = i_limit
           
           ! ----------------------------------------
           ! stop the code if we reach i_jump, the
           ! minimal index where the perturbation
           ! method can work
           ! ----------------------------------------
           IF (i_limit == i_jump) THEN
              WRITE(*,*) "Problem in starting C* trajectory"
              STOP
           ENDIF

           ! ----------------------------------------
           ! reinitialize all quantities at i_limit
           ! ----------------------------------------
           CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
           corr_flux_dum = 1.0_DP / Vi
           t_limit = distance
           dVn_old = dVn
           dVi_old = dVi

           ! ----------------------------------------
           ! Compute the maximum perturbations on 
           ! ion and neutral velocities
           ! ----------------------------------------
           IF ( counter_main_old /= counter_main ) THEN
              ! -------------------------------------
              ! reinitialize identification of up and
              ! down trajectories
              ! -------------------------------------
              found_high = .false.
              found_down = .false.
              ! -------------------------------------
              ! switch to an adjusted trajectory
              ! -------------------------------------
              adjust_start = .true.
              nadj =  0
              angle_traj = 0.0_dp
              CALL COMPUTE_ORTHO_TRAJ(i_limit,eps_traj,coef1,coef2)
              
              dum = 1.0_dp
              DO
                 IF( dum < 1e-7_dp ) THEN
                    CALL ROTATE_ORTHO_TRAJ(coef1,coef2, -angle_traj)
                    IF ( angle_traj < 0.0_dp ) THEN
                       angle_traj = - angle_traj
                    ELSE
                       angle_traj = - angle_traj - pi / 12.0_dp
                    ENDIF
                    IF ( ABS(angle_traj) > pi / 2.0_dp ) THEN
                       WRITE(*,*) "Problem in starting C* trajectory (trajectory rotation)"
                       STOP
                    ENDIF
                    CALL ROTATE_ORTHO_TRAJ(coef1,coef2, angle_traj)
                    dum = 1.0_dp
                 ENDIF
                 Vn_high = traj_main(i_limit)%Vn + dum * 1.5_dp * eps_traj * coef1
                 Vi_high = traj_main(i_limit)%Vi - dum * 1.5_dp * eps_traj * coef2
                 Vn = Vn_high
                 Vi = Vi_high
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO
              DO
                 IF( dum < 1e-7_dp ) THEN
                    WRITE(*,*) "Problem in starting C* trajectory (variation amplitude)"
                    STOP
                 ENDIF
                 Vn_down = traj_main(i_limit)%Vn - dum * 1.5_dp * eps_traj * coef1
                 Vi_down = traj_main(i_limit)%Vi + dum * 1.5_dp * eps_traj * coef2
                 Vn = Vn_down
                 Vi = Vi_down
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO

              Vn_high = traj_main(i_limit)%Vn + dum * 1.5_dp * eps_traj * coef1
              Vi_high = traj_main(i_limit)%Vi - dum * 1.5_dp * eps_traj * coef2
              Vn_down = traj_main(i_limit)%Vn - dum * 1.5_dp * eps_traj * coef1
              Vi_down = traj_main(i_limit)%Vi + dum * 1.5_dp * eps_traj * coef2
           ENDIF

           ! ----------------------------------------
           ! Apply perturbations on velocities
           ! ----------------------------------------
           IF (nadj == 0) THEN
              alpha = 1.0_dp / 2.0_dp**(nadj+1)
              nadj = nadj + 1
           ELSE
              nadj = nadj + 1
           ENDIF
           Vn = alpha * Vn_high + (1.0_dp - alpha) * Vn_down
           Vi = alpha * Vi_high + (1.0_dp - alpha) * Vi_down
           v_variab(iv_Vn) = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vn) = LOG(Vn)
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
        ENDIF

        ! -------------------------------------------
        ! Adjust trajectory after the jump to find
        ! the only solution : either the second sonic
        ! point or the downstream point
        ! -------------------------------------------
        IF ( ( ( i_max_jump - i_min_jump == 1 ).OR.( adjust_traj ) ).AND.( .NOT. sonic_point ) ) THEN
           i = i_limit + 1
           IF (found_high .AND. found_down) THEN
              DO
                 CALL FIND_DISTANCE(traj_high(i)%distance, Nstep_max, traj_down, counter_down, i_down)
                 WRITE(*,*) traj_high(i-1)%distance,      traj_high(i)%distance,      traj_high(i+1)%distance,&
                            traj_down(i_down-1)%distance, traj_down(i_down)%distance, traj_down(i_down+1)%distance
                 i_down = MAX(i_down, i_limit+1)
                 ! i_down = MAX(i_down, i_adjst+1)
                 dum = traj_down(i_down+1)%distance - traj_down(i_down)%distance
                 IF ( dum /= 0.0_dp .AND. i_down /= counter_down) THEN
                    dum1 = (traj_high(i)%distance        - traj_down(i_down)%distance) / dum &
                         * (traj_down(i_down+1)%Vn - traj_down(i_down)%Vn) + traj_down(i_down)%Vn
                    dum2 = (traj_high(i)%distance        - traj_down(i_down)%distance) / dum &
                         * (traj_down(i_down+1)%Vi - traj_down(i_down)%Vi) + traj_down(i_down)%Vi
                 ELSE
                    dum1 = traj_down(i_down)%Vn
                    dum2 = traj_down(i_down)%Vi
                 ENDIF
                 IF ( ( ABS( traj_high(i)%Vn - dum1 ) / traj_high(i)%Vn > eps_traj ).OR.&
                      ( ABS( traj_high(i)%Vi - dum2 ) / traj_high(i)%Vi > eps_traj ).OR.&
                      ( i == counter_high ) ) THEN
                    EXIT
                 ENDIF
                 i = i + 1
              ENDDO
           ENDIF
           i = i - 1

           IF ( i > i_limit ) THEN
              ! -------------------------------------
              ! switch to an adjusted trajectory
              ! -------------------------------------
              adjust_traj = .true.
              ! -------------------------------------
              ! average up and down trajectory and
              ! copy in main trajectory
              ! -------------------------------------
              CALL AVERAGE_TRAJEC ( traj_high, traj_down, traj_main, Nstep_max, counter_down, i_limit+1, i, i_adjst)
              counter_main = i
              nadj  = 0
           ENDIF

           IF ( adjust_traj .AND. (nadj == 0 .OR. nadj == 8) ) THEN
              ! -------------------------------------
              ! reinitialize identification of up and
              ! down trajectories
              ! reinitialize or change angle
              ! reinitialize alpha
              ! reinitialize nadj
              ! Careful - order is important
              ! -------------------------------------
              found_high = .false.
              found_down = .false.
              IF ( nadj == 0 ) THEN
                 angle_traj = 0.0_dp
              ENDIF
              IF ( nadj == 8 ) THEN
                 IF ( angle_traj < 0.0_dp ) THEN
                    angle_traj = - angle_traj
                 ELSE
                    angle_traj = - angle_traj - pi / 12.0_dp
                 ENDIF
              ENDIF
              IF ( ABS(angle_traj) > pi / 2.0_dp ) THEN
                 WRITE(*,*) "Problem in rotating trajectory"
                 STOP
              ENDIF
              nadj  = 0
              alpha = 1.0_dp / 2.0_dp**(nadj+1)

              ! -------------------------------------
              ! Compute the maximum perturbations on 
              ! ion and neutral velocities
              ! -------------------------------------
              CALL REINIT_PHYS_VARIABLES(traj_main(i))
              dVn_old = dVn
              dVi_old = dVi
              CALL COMPUTE_ORTHO_TRAJ(i,eps_traj,coef1,coef2)
              CALL ROTATE_ORTHO_TRAJ(coef1,coef2,angle_traj)

              dum = 1.0_dp
              DO
                 Vn_high = traj_main(i)%Vn + dum * 1.5_dp * eps_traj * coef1
                 Vi_high = traj_main(i)%Vi - dum * 1.5_dp * eps_traj * coef2

                 Vn = Vn_high
                 Vi = Vi_high
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO
              DO
                 Vn_down = traj_main(i)%Vn - dum * 1.5_dp * eps_traj * coef1
                 Vi_down = traj_main(i)%Vi + dum * 1.5_dp * eps_traj * coef2

                 Vn = Vn_down
                 Vi = Vi_down
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO

              Vn_high = traj_main(i)%Vn + dum * 1.5_dp * eps_traj * coef1
              Vi_high = traj_main(i)%Vi - dum * 1.5_dp * eps_traj * coef2
              Vn_down = traj_main(i)%Vn - dum * 1.5_dp * eps_traj * coef1
              Vi_down = traj_main(i)%Vi + dum * 1.5_dp * eps_traj * coef2
           ENDIF

           IF ( adjust_traj ) THEN
              nadj = nadj + 1
              CALL REINIT_PHYS_VARIABLES(traj_main(i))
              corr_flux_dum = 1.0_DP / Vi
              dVn_old = dVn
              dVi_old = dVi
              t_limit = traj_main(i)%distance
              ! -------------------------------------
              ! slightly change ion and neutral
              ! velocities between Vn_high & Vn_down
              ! and between Vi_high & Vi_down
              ! -------------------------------------
              Vn = alpha * Vn_high + (1.0_dp - alpha) * Vn_down
              Vi = alpha * Vi_high + (1.0_dp - alpha) * Vi_down
              v_variab(iv_Vn) = Vn
              v_variab(iv_Vi) = Vi
              v_lvariab(iv_Vn) = LOG(Vn)
              v_lvariab(iv_Vi) = LOG(Vi)
              corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
              mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
           ENDIF
           i_limit = i
        ENDIF

        ! -------------------------------------------
        ! Criteria to end the trajectory - either
        ! crossing the sonic point or recoupling the
        ! fluid
        ! -------------------------------------------
        IF ( adjust_traj ) THEN
           IF ( ( ABS(traj_main(i_limit)%Vn - traj_main(i_limit)%Vsound) / Vs_cm < eps_traj ).AND.&
                ( traj_main(i_limit)%Vn < traj_main(i_limit)%Vsound ) ) THEN
              ! -------------------------------------
              ! test if we can cross the sonic point
              ! or not
              ! -------------------------------------
              CALL TEST_CROSS_SONIC(i_limit, eps_traj, test_sonic)
              IF ( test_sonic ) THEN
                 ! ----------------------------------
                 ! reinitialize variables
                 ! ----------------------------------
                 sonic_point = .true.
                 CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
                 dVn_old = dVn
                 dVi_old = dVi
                 ! ----------------------------------
                 ! extrapolate all variables at
                 ! d_extrapo = d + 2*(d_cross - d)
                 ! ----------------------------------
                 step_cross = 1.0_dp
                 CALL CROSS_SONIC_POINT(i_limit, eps_traj, t_limit, step_cross)
                 distance  = t_limit
                 dist_step = distance - traj_main(i_limit)%distance
                 sonic_jump = distance
                 ! ----------------------------------
                 ! recouple the fluids is necessary
                 ! DOES NOT WORK - recoupling just
                 ! after the jump is numerically
                 ! unstable - dont know why yet
                 ! ----------------------------------
                 ! IF( Vn < Vi .OR. ABS(Vn-Vi) / Vs_cm < eps_traj ) THEN
                 !    ! -------------------------------------
                 !    ! set neutral velocity to ion velocity
                 !    ! assuming neutrals recouple with ions
                 !    ! and not the opposite
                 !    ! -------------------------------------
                 !    Vn = Vi
                 !    v_variab(iv_Vn) = Vn
                 !    v_lvariab(iv_Vn) = LOG(Vn)
                 !    ! -------------------------------------
                 !    ! set up a J-type shock
                 !    ! -------------------------------------
                 !    ! prevent from switching back to C shocks
                 !    shock_type_lock = "J"
                 ! ENDIF
                 DO WHILE ( Vn < Vi )
                    step_cross = step_cross / 2.0_dp
                    CALL CROSS_SONIC_POINT(i_limit, eps_traj, t_limit, step_cross)
                    distance  = t_limit
                    dist_step = distance - traj_main(i_limit)%distance
                    sonic_jump = distance
                    IF( step_cross < eps_traj ) THEN
                       WRITE(*,*) "Problem in jumping sonic point"
                       STOP
                    ENDIF
                 ENDDO
              ENDIF
           ELSE IF ( ABS(traj_main(i_limit)%Vn - traj_main(i_limit)%Vi) / Vs_cm < eps_traj ) THEN
              ! -------------------------------------
              ! set up a J-type shock
              ! -------------------------------------
              fluid_recoup = .true.
              ! prevent from switching back to C shocks
              shock_type_lock = "J"
              ! -------------------------------------
              ! reinitialize variables
              ! -------------------------------------
              CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
              dVn_old = dVn
              dVi_old = dVi
              t_limit = traj_main(i_limit)%distance
              distance = t_limit
              ! -------------------------------------
              ! set neutral velocity to ion velocity
              ! assuming neutrals recouple with ions
              ! and not the opposite
              ! -------------------------------------
              ! CALL FALSE_JUMP(Vi)
              ! Vn = Vi
              ! v_variab(iv_Vn) = Vn
              ! v_lvariab(iv_Vn) = LOG(Vn)
              ! -------------------------------------
              ! set ion velocity to neutral velocity
              ! assuming ions recouple with neutrals
              ! and not the opposite
              ! -------------------------------------
              corr_flux_dum = 1.0_DP / Vi
              Vi = Vn
              v_variab(iv_Vi) = Vi
              v_lvariab(iv_Vi) = LOG(Vi)
              corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
              mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
          ENDIF
        ENDIF

        ! ----------------------------------------
        ! restart integration
        ! ----------------------------------------
        counter = i_limit
        restart = .true.
     ENDIF

  END SUBROUTINE TEST_STAT_CJ


END MODULE MODULE_DRIVE
