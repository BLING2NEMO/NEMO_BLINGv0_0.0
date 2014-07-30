MODULE trcwri_bling

#if defined key_top && key_bling && defined key_iomput

   USE trc
   USE iom
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_bling

CONTAINS

   SUBROUTINE trc_wri_bling

      INTEGER           :: jn
      CHARACTER(len=20) :: cltra

      DO jn=jp_bling0, jp_bling1
         cltra = TRIM( ctrcnm(jn) )  
         CALL iom_put( cltra, trn(:,:,:,jn)*1.0e+6 )
      END DO

      ! diagnostic tracers output
      IF ( ln_diatrc ) THEN
         CALL iom_put( "CHL_bling" , chl_bling(:,:,:) * tmask(:,:,:) )
         CALL iom_put( "dum_bling" ,      dum1(:,:,:) * tmask(:,:,:) )
      ENDIF

   END SUBROUTINE trc_wri_bling

#else

   PUBLIC trc_wri_bling
CONTAINS
   SUBROUTINE trc_wri_bling
   END SUBROUTINE trc_wri_bling

#endif

END MODULE trcwri_bling
