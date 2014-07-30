MODULE trcsed_bling

#if defined key_bling

   !!======================================================================
   !!                     ***  MODULE trcsed_bling  ***
   !! TOP :  Phosphorous, iron, and oxygen mass exchange with bottom sediments
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!======================================================================

   USE oce_trc
   USE par_trc
   USE trc

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_sed_bling

# include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_sed_bling 

      INTEGER  :: ji, jj, ikt

      REAL(wp) :: zrfact, f_oxy, fe_2_p_sed
      REAL(wp) :: bfpo4   , bfoxy   , bffed

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sed_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      !---------------------------------------------------------------------
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor. For P, the bottom flux  
      ! puts the sinking flux reaching the bottom cell into the water column 
      ! through diffusion. For iron, the sinking flux disappears into the 
      ! sediments if bottom waters are oxic (assumed adsorbed as oxides),
      ! while an efflux of dissolved iron occurs dependent on the supply of
      ! reducing organic matter (scaled by the org-P sedimentation rate).
      ! If bottom waters are anoxic, the sinking flux of Fe is returned to
      ! the water column. Note this is not appropriate for very long runs
      ! with an anoxic ocean (iron will keep accumulating forever).
      ! For oxygen, the consumption of oxidant required to respire  
      ! the settling flux of organic matter (in support of the
      ! PO4 bottom flux) diffuses from the bottom water into the sediment.
      !---------------------------------------------------------------------

      fe_2_p_sed=1.e-4 * 106.0

      sumblpo4=0.d0
      sumblfed=0.d0
      sumbloxy=0.d0

      DO ji=1, jpi
         DO jj=1, jpj

            ! mbkt is a matrix containing the vertical index of the
            ! bottom layer at each horizontal point
            ikt    = mbkt(ji,jj)
            zrfact = rfact/fse3t(ji,jj,ikt)

            ! Phosphate
            bfpo4 =-       fpop_b(ji,jj)*zrfact

            ! Oxygen
            bfoxy = o2_2_p*fpop_b(ji,jj)*zrfact

            ! Iron
            f_oxy = trn(ji,jj,ikt,jpOxy_bling)
            IF (f_oxy>oxy_min) THEN
               bffed=- fe_2_p_sed*fpop_b(ji,jj)*zrfact 
            ELSE
               bffed=-(fe_2_p_sed*fpop_b(ji,jj)+fpofe_b(ji,jj))*zrfact
            ENDIF

            ! Total mass exchange w/ sediments
            sumblpo4=sumblpo4+bfpo4*cvol(ji,jj,ikt)
            sumbloxy=sumbloxy+bfoxy*cvol(ji,jj,ikt)
            sumblfed=sumblfed+bffed*cvol(ji,jj,ikt)

            ! Add mass exchange w/ sediments to trends
            tra(ji,jj,ikt,jpPO4_bling) = tra(ji,jj,ikt,jpPO4_bling) + bfpo4
            tra(ji,jj,ikt,jpFed_bling) = tra(ji,jj,ikt,jpFed_bling) + bffed
            tra(ji,jj,ikt,jpOxy_bling) = tra(ji,jj,ikt,jpOxy_bling) + bfoxy

         ENDDO
      ENDDO

   END SUBROUTINE trc_sed_bling

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sed_bling 
      WRITE(*,*) 'trc_sed_bling: You should not have seen this print',  kt
   END SUBROUTINE trc_sed_bling

#endif

END MODULE trcsed_bling
