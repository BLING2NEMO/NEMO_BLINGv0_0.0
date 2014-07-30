MODULE trcini_bling
   !!======================================================================
   !!                         ***  MODULE trcini_bling  ***
   !! TOP :   initialisation of the BLINGv0 tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_bling
   !!----------------------------------------------------------------------
   !!   'key_bling'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_bling   : BLINGv0 model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcsms_bling
   !: PISCES indices for PO4, DOC, O2, F to initialize BLING fields
   USE par_pisces , ONLY : jppo4, jpdoc, jpoxy, jpfer 

   ! BLINGv0 modules
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_bling   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_bling.F90 -1   $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_bling
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_bling  ***  
      !!
      !! ** Purpose :   initialization for BLINGv0 model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------

      !                       ! Allocate BLINGv0 arrays

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_bling: initialisation of BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! Allocate BLINGv0 arrays
      IF( bling_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_bling: unable to allocate BLINGv0 arrays' )
      
      rfact   = rdttrc(1)                          

      IF(lwp) WRITE(numout,*) '    Passive Tracer time step    rfact  = ', rfact, ' rdt = ', rdttra(1)

      IF( .NOT. ln_rsttr ) THEN
        ! initialize prognostic tracers
        trn(:,:,:,jpPO4_bling) = trn(:,:,:,jppo4) * tmask(:,:,:) / 106.d0
        trn(:,:,:,jpDOP_bling) = trn(:,:,:,jpdoc) * tmask(:,:,:) / 106.d0
        trn(:,:,:,jpFed_bling) = trn(:,:,:,jpfer) * tmask(:,:,:)
        trn(:,:,:,jpOxy_bling) = trn(:,:,:,jpoxy) * tmask(:,:,:)

!        trn(:,:,:,jpPO4_bling) = 1.d-6 * tmask(:,:,:) 
!        trn(:,:,:,jpDOP_bling) = 2.d-6 * tmask(:,:,:) 
!        trn(:,:,:,jpPO4_bling) = trn(:,:,:,jppo4) * tmask(:,:,:)
         
        ! initialize diagnostic tracers
        chl_bling(:,:,:) = 0.d0
        dum1     (:,:,:) = 0.d0
      ENDIF

      !
   END SUBROUTINE trc_ini_bling

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_bling             ! Empty routine
   END SUBROUTINE trc_ini_bling
#endif

   !!======================================================================
END MODULE trcini_bling
