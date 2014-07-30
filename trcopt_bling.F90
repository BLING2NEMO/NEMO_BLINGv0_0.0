MODULE trcopt_bling

#if defined key_bling

   !!======================================================================
   !!                     ***  MODULE trcopt_bling  ***
   !! TOP :  BLINGv0 Compute the light availability in the water column
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   !! Available light calculation in BLINGv0
   !!-----------------------------------------------------------------------
   !! There are multiple types of light.
   !!   IRR_INST is the instantaneous irradiance field.
   !!   IRR_MIX  is the same, but with the irr_inst averaged throughout the  
   !! mixed layer (turbocline*) to account for mixing directly below the boundary 
   !! layer. This quantity is intended to represent the light to which phytoplankton 
   !! subject to turbulent transport in the mixed-layer would be exposed.
   !!   IRR_MEM  is a temporally smoothed field carried between timesteps, to 
   !! represent photoadaptation.
   !!-----------------------------------------------------------------------
   !! *The turbocline depth is the depth at which the
   !! vertical eddy diffusivity coefficient (resulting from the vertical physics
   !! alone, not the isopycnal part, see trazdf.F) fall below a given value
   !! defined locally (avt_c here taken equal to 5 cm/s2). Check subroutine
   !! OPA_SRC/ZDF/zdfmxl.F90 for further explanation.
   !!-----------------------------------------------------------------------

   USE oce_trc
   USE trc
   USE prtctl_trc      ! Print control for debbuging
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_opt_bling ! called by trc_sms_bling

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_opt_bling (kt, irr_inst, irr_mix)

      INTEGER , INTENT(in) :: kt
      REAL(wp), POINTER, DIMENSION(:,:,:), INTENT(out) :: irr_mix, irr_inst

      INTEGER  :: ji, jj, jk, kblt          
      REAL(wp) :: zpig, zkr, zkb, zparr, zparb, zcoef, parfc, sumz_irrad_ML, sumz_hblt

      CALL wrk_alloc( jpi, jpj, jpk, irr_inst, irr_mix )

      IF( nn_timing == 1 )  CALL timing_start('trc_opt_bling')

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_opt_bling : BLING optic-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~ '
      ENDIF

      ! light attenuation by water and phytoplankton
      zcoef=12*redf/rcchl/rpig
      parfc=0.43d0/2.d0
      DO ji=1, jpi
         DO jj=1, jpj

            !One-wavelength biopt
            zparr=qsr(ji,jj)*0.43d0
            irr_inst(ji,jj,1)=zparr
            irr_mix (ji,jj,1)=zparr

            kblt=1
            !Irradiance at the surface
            sumz_irrad_ML=irr_mix(ji,jj,1)*fse3t(ji,jj,1)
            sumz_hblt    =fse3t(ji,jj,1)

            DO jk=2, jpk

               zkr   = xkr0 + xkrp * chl_bling(ji,jj,jk-1)
               zparr = zparr * EXP(-zkr*fse3t(ji,jj,jk-1))
               
               irr_inst(ji,jj,jk)=zparr
               irr_mix (ji,jj,jk)=irr_inst(ji,jj,jk)

               ! irradiance sum within the MLD
               IF (sumz_hblt < hmld(ji,jj)) THEN
                  kblt=kblt+1
                  sumz_irrad_ML=sumz_irrad_ML+irr_mix(ji,jj,jk)*fse3t(ji,jj,jk)
                  sumz_hblt    =sumz_hblt+fse3t(ji,jj,jk)
               ENDIF

            END DO

            ! irradiance mean average within the MLD
            irr_mix(ji,jj,1:kblt)=sumz_irrad_ML / MAX(1.0e-6,sumz_hblt)

         END DO
      END DO

      ! Initialize memory irradiance
      IF (kt==nittrc000) irr_mem(:,:,:)=irr_mix(:,:,:)

      ! Phytoplankton photoadaptation timescale
      ! Forward time-stepping for memory irradiance
      DO jk=1, jpk
         DO jj=1, jpj
            DO ji=1, jpi
               irr_mem(ji,jj,jk)=irr_mem(ji,jj,jk)                                                    &
                                 + (irr_mix (ji,jj,jk)-irr_mem(ji,jj,jk))*MIN(1.d0,gam_irr_mem*rfact)*tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL wrk_dealloc( jpi, jpj, jpk, irr_inst, irr_mix)

      IF( nn_timing == 1 )  CALL timing_stop ('trc_opt_bling')

   END SUBROUTINE trc_opt_bling

#else
   !!======================================================================
   !!  Dummy module :                                      No BLINGv0 model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_opt_bling   ! Empty routine
      WRITE(*,*) 'trc_opt_bling: You should have not entered this routine, error?'
   END SUBROUTINE trc_opt_bling 

#endif

END MODULE trcopt_bling
