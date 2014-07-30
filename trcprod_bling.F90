MODULE trcpod_bling

#if defined key_bling

   USE oce_trc
   USE trc
   USE prtctl_trc      ! Print control for debbuging
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_pod_bling ! called by trc_sms_bling

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_pod_bling (kt)

      INTEGER :: kt

      IF( nn_timing == 1 )  CALL timing_start('trc_pod_bling')

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_pod_bling : BLING production'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~ '
      ENDIF


      IF( nn_timing == 1 )  CALL timing_stop ('trc_pod_bling')

   END SUBROUTINE trc_pod_bling

#else

   SUBROUTINE trc_pod_bling (kt)  ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_pod_bling: You should have not entered this routine, error?'
   END SUBROUTINE trc_pod_bling 

#endif

END MODULE trcpod_bling
