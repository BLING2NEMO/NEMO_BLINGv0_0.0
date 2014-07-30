MODULE trcrst_bling
   !!======================================================================
   !!                       ***  MODULE trcrst_bling  ***
   !! TOP :   create, write, read the restart files of BLINGv0 tracer
   !!======================================================================
   !! History :   1.0  !  2010-01 (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_bling
   !!----------------------------------------------------------------------
   !!   'key_bling'                                               CFC tracers
   !!----------------------------------------------------------------------
   !!   trc_rst_read_bling   : read  restart file
   !!   trc_rst_wri_bling    : write restart file
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   PUBLIC  trc_rst_read_bling   ! called by trcini.F90 module
   PUBLIC  trc_rst_wri_bling   ! called by trcini.F90 module

CONTAINS
   
   SUBROUTINE trc_rst_read_bling( knum ) 
     INTEGER, INTENT(in)  :: knum
     WRITE(*,*) 'trc_rst_read_bling: No specific variables to read on unit', knum
   END SUBROUTINE trc_rst_read_bling

   SUBROUTINE trc_rst_wri_bling( kt, kirst, knum )
     INTEGER, INTENT(in)  :: kt, kirst, knum
     WRITE(*,*) 'trc_rst_wri_bling: No specific variables to write on unit' ,knum, ' at time ', kt, kirst
   END SUBROUTINE trc_rst_wri_bling

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read_bling( knum )
     INTEGER, INTENT(in)  :: knum
     WRITE(*,*) 'trc_rst_read_bling: You should not have seen this print! error?', knum
   END SUBROUTINE trc_rst_read_bling

   SUBROUTINE trc_rst_wri_bling( kt, kirst, knum )
     INTEGER, INTENT(in)  :: kt, kirst, knum
     WRITE(*,*) 'trc_rst_wri_bling: You should not have seen this print! error?', kt, kirst, knum
   END SUBROUTINE trc_rst_wri_bling
#endif

   !!======================================================================
END MODULE trcrst_bling
