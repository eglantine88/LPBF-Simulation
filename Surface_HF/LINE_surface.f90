
! Laser Heat Flux
SUBROUTINE DFLUX(FLUX, SOL, KSTEP, KINC, TIME, NOEL, NPT, COORDS, JLTYP, TEMP, PRESS, SNAME)
  INCLUDE 'ABA_PARAM.INC'
  INCLUDE 'CONSTANT_DFLUX.INC'

  ! Arguments
  DOUBLE PRECISION :: FLUX(2), SOL(*), TIME(2), COORDS(3), TEMP, PRESS
  INTEGER :: KSTEP, KINC, NOEL, NPT, JLTYP
  CHARACTER(LEN=80) :: SNAME

  ! Constants
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589D0
  DOUBLE PRECISION, PARAMETER :: R = 8.314D0           ! J/mol-K

  ! Variables
  DOUBLE PRECISION :: X0, Z0, X, Y, Z, T, TP, XC, YC, ZC, DIST2, R02, Q_laser ! Q_evap, PSAT
  INTEGER :: STEP, MODUL, LAYER

  X = COORDS(1)
  Y = COORDS(2)
  Z = COORDS(3)
  T = TIME(1)

  ! EVAPORATION
  TP = SOL(1)
  PSAT = 101325 * 10**(6.1127 - 18868 / TP) ! To change depending on the material
  Q_evap = - 0.82 * DHV * PSAT / SQRT(2.0 * M * PI * R * TP)

  ! LASER
  
  ! Step info
  STEP = KSTEP-2
  MODUL = INT(MOD(STEP,(2*NB_LASER))/2) ! n° of laser path
  LAYER = INT(STEP/(2*NB_LASER)) ! n° of current layer

  ! Laser Origin at the very beginning
  X0 = - H_SPACING
  Z0 = H_SPACING / 2

  ! Coordinates of the moving laser center
  XC = X0 + LASER_SPEED * T
  YC = LAYER_THICKNESS*(LAYER+1)
  ZC = Z0 + MODUL*H_SPACING

  Q_laser = 0.0

  DIST2 = (X - XC)**2 + (Z - ZC)**2
  R02 = LASER_RADIUS**2
  Q_laser = 2 * LASER_POWER / (PI * R02) * EXP(-2 * DIST2 / R02)
  FLUX(1) = Q_laser * ABSORB + Q_evap

  RETURN
END SUBROUTINE DFLUX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Used to adjust the density wheither it is powder or solid/Liquid/Mushy -> we suppose that mushy zone has the same density as the solid
SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
    INCLUDE 'ABA_PARAM.INC'

    DOUBLE PRECISION :: FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),T(3,3),TIME(2)
    DOUBLE PRECISION :: ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
    CHARACTER(LEN=80) :: CMNAME,ORNAME
    CHARACTER(LEN=3) :: FLGRAY(15)
    
    IF (STATEV(1) .LT. 1.0) THEN 
        FIELD(1) = 0   ! Powder, Mushy
    ELSE                       
        FIELD(1) = 1 ! Solid, Mushy, Liquid
    ENDIF

END SUBROUTINE USDFLD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize state (0=POWDER, 1=SOLID, 2=LIQUID) (not useful)
SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
  INCLUDE 'ABA_PARAM.INC' 

  DOUBLE PRECISION :: STATEV(NSTATV),COORDS(NCRDS)
  STATEV(1) = 0 ! Powder
  IF (COORDS(2) .LT. 0.0D0) THEN ! Support (y<0) = Solid
    STATEV(1) = 2 ! Solid
  END IF

  RETURN
END SUBROUTINE SDVINI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the conductivity and internal energy depending on the state (powder, solid, liquid)
SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
  INCLUDE 'ABA_PARAM.INC'
  INCLUDE 'CONSTANT_UMATHT.INC'

  ! Arguments
  DOUBLE PRECISION :: DUDG(NTGRD), FLUX(NTGRD), DFDT(NTGRD), DFDG(NTGRD,NTGRD), STATEV(NSTATV), DTEMDX(NTGRD)
  DOUBLE PRECISION :: TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DTEMP, DTIME, TEMP, U
  CHARACTER(LEN=80) :: CMNAME

  ! Variables
  DOUBLE PRECISION :: TEMP_dT, LATENT_FRAC, LATENT_FRAC_dT, LATENT_FRAC_INT, dU, SLOPE, TEMP_INT, K, DKDT, TL0, TL1, DT, K0, K1, DK, C0, C1, DC, F
  DOUBLE PRECISION :: LATENT_HEAT
  INTEGER :: LOW, HIGH, MID, LEN

  ! Size of table
  LEN = SIZE(TABLE_K, 2)

  ! Material tables
  LATENT_HEAT = PROPS(1) ! PROPS are the constants defined in 'UserMaterial' 

  ! Temperature after the current increment
  TEMP_dT = TEMP + DTEMP
  
  ! Material Properties
  IF (TEMP_dT .LE. TABLE_K(2,1)) THEN
      K = MERGE((1.0-POROSITY), 1.0 , STATEV(1) == 0) * TABLE_K(1,1) ! (1-Porosity) * Conductivity if powder
      DUDT = TABLE_C(1,1)  ! Specific heat
      DKDT = 0.0d0
  ELSEIF (TEMP_dT .GE. TABLE_K(2,LEN)) THEN
      K = MERGE((1.0-POROSITY), 1.0 , STATEV(1) == 0) * TABLE_K(1,LEN) ! (1-Porosity) * Conductivity if powder
      DUDT = TABLE_C(1,LEN)
      DKDT = 0.d0
    ! If temperature within table range :
  ELSEIF (TEMP_dT .GT. TABLE_K(2,1) .AND. TEMP_dT .LT. TABLE_K(2,LEN)) THEN
      LOW = 1
      HIGH = LEN
      ! Binary search of temperature interval TL0 ≤ TEMP_dT < TL1
      DO WHILE (HIGH - LOW > 1)
        MID = (LOW + HIGH) / 2
        IF (TEMP_dT < TABLE_K(2,MID)) THEN
          HIGH = MID
        ELSE
          LOW = MID
        END IF
      END DO

      ! Compute conductivity and Specific Heat
      TL0 = TABLE_K(2,LOW)
      TL1 = TABLE_K(2,HIGH)
      DT = TL1 - TL0

      K0 =  MERGE(TABLE_K_P(1,LOW), TABLE_K(1,LOW) , STATEV(1) == 0) 
      K1 =  MERGE(TABLE_K_P(1,HIGH), TABLE_K(1,HIGH) , STATEV(1) == 0) 
      DK = K1 - K0 
      DKDT = DK / DT
      K = DKDT * (TEMP_dT - TL0) + K0

      C0 = TABLE_C(1,LOW)
      C1 = TABLE_C(1,HIGH)
      DC = C1 - C0
      DUDT = (DC / DT) * (TEMP_dT - TL0) + C0
  END IF

  ! Change of state (0=POWDER, 2=SOLID, 1=LIQUID)
  IF (TEMP_dT .GE. LIQUIDUS_TEMP) THEN
      STATEV(1) = 1                      ! Solid/Powder -> Liquid
  ELSE IF (TEMP_dT .LT. SOLIDUS_TEMP) THEN
      IF (STATEV(1) .EQ. 1) THEN         ! Liquid -> Solid
          STATEV(1) = 2
      END IF
  END IF

  ! Variation of internal energy due to sensible heat
  dU = DUDT * DTEMP

  ! Latent heat effects
  IF (TEMP .LT. LIQUIDUS_TEMP .AND. TEMP .GT. SOLIDUS_TEMP) THEN
    LATENT_FRAC = LATENT_HEAT * (TEMP - SOLIDUS_TEMP) / (LIQUIDUS_TEMP - SOLIDUS_TEMP)
  ELSEIF (TEMP .GE. LIQUIDUS_TEMP) THEN
    LATENT_FRAC = LATENT_HEAT
  END IF

  IF (TEMP_dT .LT. LIQUIDUS_TEMP .AND. TEMP_dT .GT. SOLIDUS_TEMP) THEN
    LATENT_FRAC_dT = LATENT_HEAT * (TEMP_dT - SOLIDUS_TEMP) / (LIQUIDUS_TEMP - SOLIDUS_TEMP)
    SLOPE = LATENT_HEAT / (LIQUIDUS_TEMP - SOLIDUS_TEMP)
  ELSEIF (TEMP_dT .GE. LIQUIDUS_TEMP) THEN
    LATENT_FRAC_dT = LATENT_HEAT
    SLOPE = 0.0D0
  END IF

  ! If phase change occurs, add latent heat contribution
  IF (LATENT_FRAC .NE. LATENT_FRAC_dT) THEN
    dU = dU + LATENT_FRAC_dT - LATENT_FRAC
    DUDT = DUDT + SLOPE

    ! Smoothing energy slope to prevent discontinuity
    IF (SLOPE .EQ. 0.0D0) THEN
      TEMP_INT = TEMP_dT - 0.25 * DTEMP
      IF (TEMP_INT .LT. LIQUIDUS_TEMP .AND. TEMP_INT .GT. SOLIDUS_TEMP) THEN
        LATENT_FRAC_INT = LATENT_HEAT * (TEMP_INT - SOLIDUS_TEMP) / (LIQUIDUS_TEMP - SOLIDUS_TEMP)
        SLOPE = LATENT_HEAT / (LIQUIDUS_TEMP - SOLIDUS_TEMP)
      ELSEIF (TEMP_INT .GE. LIQUIDUS_TEMP) THEN
        LATENT_FRAC_INT = LATENT_HEAT
        SLOPE = 0.0D0
      END IF

      IF (LATENT_FRAC_INT .NE. LATENT_FRAC_dT) THEN
        DUDT = DUDT + SLOPE
      END IF
    END IF
  END IF

  ! Update internal energy
  U = U + dU

  ! Heat Flux (-k * dT/dx)
  DO I = 1, NTGRD
    FLUX(I) = -K * DTEMDX(I)
    DFDG(I, I) = -K
    DFDT(I) = -DKDT * DTEMDX(I)
  END DO

  RETURN

END SUBROUTINE UMATHT