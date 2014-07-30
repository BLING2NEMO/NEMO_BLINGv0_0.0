MODULE trcfe_bling

#if defined key_bling

   USE oce_trc
   USE trc
   USE prtctl_trc      ! Print control for debbuging
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_feup_bling ! called by trc_sms_bling

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_feup_bling (kt)

      INTEGER, INTENT(in) :: kt
      INTEGER  :: ji, jj, jk
      REAL(wp), POINTER, DIMENSION(:,:,:) :: fe2p

      CALL wrk_alloc( jpi, jpj, jpk, fe2p )

      IF( nn_timing == 1 )  CALL timing_start('trc_feup_bling')

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_feup_bling : Iron uptake'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~ '
      ENDIF

      CALL wrk_dealloc( jpi, jpj, jpk, fe2p)

      IF( nn_timing == 1 )  CALL timing_stop ('trc_feup_bling')

   END SUBROUTINE trc_feup_bling

#else

   SUBROUTINE trc_feup_bling (kt)  ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_fe_bling: You should have not entered this routine, error?'
   END SUBROUTINE trc_feup_bling 

#endif

END MODULE trcfe_bling
