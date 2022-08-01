!---!------------------------------------------------------------------!------!
!---!------------------------------------------------------------------!------!
! ideal J2 plasticity in FORTRAN                                              !
! Denis Novokshanov                                                           !
! 01.08.2022                                                                  !
! source code to the document                                                 !
! Ideal J2 plasticity in 9 steps - A3 - Version 2022.06.03.pdf                !
!---!------------------------------------------------------------------!------!
      PROGRAM J2
      WRITE(*,*) "START!"
!----------------------------------------------------------------------!
! INPUT
!----------------------------------------------------------------------!
      EMOD  = 210000.0 ! in MPa
      ENU   = 0.3      !
      SIGY0 = 500.0    ! in MPa
!----------------------------------------------------------------------!
      TEEL1  =   0.01
      TEEL2  =  -0.004
      TEEL3  =  -0.004
      TEEL4  =  0.0
      TEEL5  =  0.0
      TEEL6  =  0.0
!----------------------------------------------------------------------!
! CALCULATION
!----------------------------------------------------------------------!
! Step 1
!----------------------------------------------------------------------!
      ELAM = (EMOD*ENU)/((1.0+ENU)*(1.0-2.0*ENU))
      EMU  = (EMOD)/(2.0*(1.0+ENU))
!----------------------------------------------------------------------!
      TST1 = (ELAM+2.0*EMU)*TEEL1 + ELAM*TEEL2 + ELAM*TEEL3
      TST2 = ELAM*TEEL1 + (ELAM+2.0*EMU)*TEEL2 + ELAM*TEEL3
      TST3 = ELAM*TEEL1 + ELAM*TEEL2 + (ELAM+2.0*EMU)*TEEL3
      TST4 = EMU*2.0*TEEL4
      TST5 = EMU*2.0*TEEL5
      TST6 = EMU*2.0*TEEL6
!----------------------------------------------------------------------!
      WRITE(*,*) "Step 1"
      WRITE(*,*) "Trial Stress Tensor TST1 =", TST1
      WRITE(*,*) "Trial Stress Tensor TST2 =", TST2
      WRITE(*,*) "Trial Stress Tensor TST3 =", TST3
      WRITE(*,*) "Trial Stress Tensor TST4 =", TST4
      WRITE(*,*) "Trial Stress Tensor TST5 =", TST5
      WRITE(*,*) "Trial Stress Tensor TST6 =", TST6
!----------------------------------------------------------------------!
! Step 2  
!----------------------------------------------------------------------!
      MST = 1.0/3.0*(TST1 + TST2 + TST3)
!
      TDST1 = TST1 - MST
      TDST2 = TST2 - MST
      TDST3 = TST3 - MST
      TDST4 = TST4
      TDST5 = TST5
      TDST6 = TST6
!----------------------------------------------------------------------!
      WRITE(*,*) "Step 2"
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST1 =", TDST1
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST2 =", TDST2
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST3 =", TDST3
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST4 =", TDST4
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST5 =", TDST5
      WRITE(*,*) "Trial Deviatoric Stress Tensor TDST6 =", TDST6
!----------------------------------------------------------------------!
! Step 3
!----------------------------------------------------------------------!
       SEFNM = SQRT(TDST1**2.0 + TDST2**2.0 + TDST3**2.0 + 2.0*TDST4**2.0 + 2.0*TDST5**2.0 + 2.0*TDST6**2.0)
!
       SIGV = SQRT(3.0/2.0)*SEFNM
!
       WRITE(*,*) "Step 3"
       WRITE(*,*) "Sigma_V =", SIGV
!----------------------------------------------------------------------!
! Step 4
!----------------------------------------------------------------------!
       IF (SIGV .LT. SIGY0) THEN
!----------------------------------------------------------------------!
! Step 5
!----------------------------------------------------------------------!
      WRITE(*,*) "Step 5 - elastic step"
      EEL1 = TEEL1
      EEL2 = TEEL2
      EEL3 = TEEL3
      EEL4 = TEEL4
      EEL5 = TEEL5
      EEL6 = TEEL6
       END IF
!
       IF (SIGV .GT. SIGY0) THEN
!----------------------------------------------------------------------!
! Step 6 - flow direction
!----------------------------------------------------------------------!
      WRITE(*,*) "Step 6 - plastic step"
      FDIR1 = TDST1/SEFNM
      FDIR2 = TDST2/SEFNM
      FDIR3 = TDST3/SEFNM
      FDIR4 = TDST4/SEFNM
      FDIR5 = TDST5/SEFNM
      FDIR6 = TDST6/SEFNM
!----------------------------------------------------------------------!
      WRITE(*,*) "Flow Direction FDIR1", FDIR1
      WRITE(*,*) "Flow Direction FDIR2", FDIR2
      WRITE(*,*) "Flow Direction FDIR3", FDIR3
      WRITE(*,*) "Flow Direction FDIR4", FDIR4
      WRITE(*,*) "Flow Direction FDIR5", FDIR5
      WRITE(*,*) "Flow Direction FDIR6", FDIR6
!----------------------------------------------------------------------!
! Step 7 - plastic strain multiplyer 
!----------------------------------------------------------------------!
      GAMMA = (SEFNM - SQRT(2.0/3.0)*SIGY0)/(2.0*EMU)
      WRITE(*,*)"Step 7"
      WRITE(*,*)"Plastic Strain Multiplyer GAMMA =", GAMMA
      WRITE(*,*)"Eqvivalent Plastic Strain VEPS =", SQRT(2.0/3.0)*GAMMA 
!----------------------------------------------------------------------!
! Step 8 - plastic strain multiplyer 
!----------------------------------------------------------------------!
      EEL1 = TEEL1 - FDIR1*GAMMA
      EEL2 = TEEL2 - FDIR2*GAMMA
      EEL3 = TEEL3 - FDIR3*GAMMA
      EEL4 = TEEL4 - FDIR4*GAMMA
      EEL5 = TEEL5 - FDIR5*GAMMA
      EEL6 = TEEL6 - FDIR6*GAMMA
!
      EPL1 = FDIR1*GAMMA
      EPL2 = FDIR2*GAMMA
      EPL3 = FDIR3*GAMMA
      EPL4 = FDIR4*GAMMA
      EPL5 = FDIR5*GAMMA
      EPL6 = FDIR6*GAMMA
!
      END IF
!----------------------------------------------------------------------!
! Step 9 - calculation of updated stress 
!----------------------------------------------------------------------!
      ST1 = (ELAM+2.0*EMU)*EEL1 + ELAM*EEL2 + ELAM*EEL3
      ST2 = ELAM*EEL1 + (ELAM+2.0*EMU)*EEL2 + ELAM*EEL3
      ST3 = ELAM*EEL1 + ELAM*EEL2 + (ELAM+2.0*EMU)*EEL3
      ST4 = EMU*2.0*EEL4
      ST5 = EMU*2.0*EEL5
      ST6 = EMU*2.0*EEL6
!----------------------------------------------------------------------!
! Step 10 - check of calculation
! von Mises stress at the end of increment (SIGV2) should be smaller 
! than yield stress of the material (SIGY0)
!----------------------------------------------------------------------!
       SIGM = 1.0/3.0*(ST1 + ST2 + ST3)
!
       DST1 = ST1 - SIGM
       DST2 = ST2 - SIGM
       DST3 = ST3 - SIGM
       DST4 = ST4
       DST5 = ST5
       DST6 = ST6
!
       SEFNM = SQRT(DST1**2.0 + DST2**2.0 + DST3**2.0 + 2.0*DST4**2.0 + 2.0*DST5**2.0 + 2.0*DST6**2.0)
!
       SIGV2 = SQRT(3.0/2.0)*SEFNM
!
       WRITE(*,*) "Step 10 - check"
       IF (SIGV2 .LE. (SIGY0+1.0)) THEN
       WRITE(*,*) "Check is OK"
       ELSE
       WRITE(*,*) "ERROR !!!"
       END IF
!----------------------------------------------------------------------!
       WRITE(*,*) 
       WRITE(*,*) "SIGV2 =", SIGV2
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
      WRITE(*,*) "Output **********************************************"  
      WRITE(*,*) "Index,        eps_el,      eps_pl,             stress"
      WRITE(*,*)  " 1  ",  EEL1, "**" , EPL1, "**", ST1
      WRITE(*,*)  " 2  ",  EEL2, "**" , EPL2, "**", ST2
      WRITE(*,*)  " 3  ",  EEL3, "**" , EPL3, "**", ST3
      WRITE(*,*)  " 4  ",  EEL4, "**" , EPL4, "**", ST4
      WRITE(*,*)  " 5  ",  EEL5, "**" , EPL5, "**", ST5
      WRITE(*,*)  " 6  ",  EEL6, "**" , EPL6, "**", ST6
!----------------------------------------------------------------------!
      WRITE(*,*) "END!"
      END PROGRAM