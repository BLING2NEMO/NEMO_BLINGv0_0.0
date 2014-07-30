MODULE trcsms_bling

#if defined key_bling
   !!======================================================================
   !!                     ***  MODULE trcsms_bling  ***
   !! TOP :  BLINGv0 model main routine. Computes the sources and sinks
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE iom

   ! BLINGv0 specific modules
   USE trcopt_bling
   USE trcsed_bling
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_bling       ! called by trcsms.F90 module

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_sms_bling( kt )
      !!------------------------------------------------------------------------
      !!                     ***  trc_sms_bling  ***
      !!
      !! ** History : default leap-frog scheme(MY_TRC) changed to FORWARD scheme
      !!------------------------------------------------------------------------
      
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      INTEGER  :: ji, jj, jk, jn

      REAL(wp) :: f_po4, f_dop, f_fer, f_oxy
      REAL(wp) :: ztra
      ! Irradiance k
      REAL(wp) :: expkT, po4_up, fe2p_up, def_fe, pc_m, thetamax_fe, alpha_chl, irrk
      ! Production
      REAL(wp) :: pc_tot, mu, biomass_p_ts, theta, chl_dia, mulamb0expkT
      ! Phosphorous
      REAL(wp) :: phi_sm, phi_lg, kappa_remin, phi_dop 
      REAL(wp) :: gamma_dop, gamma_pop
      REAL(wp) :: wsink0_z, wsink0, wsink_acc, koxy, o2_2_p, remin_min !param
      REAL(wp) :: jp_uptake, frac_pop,jdop,jp_recycle
      REAL(wp) :: zzz, wsink, oxy_up, zremin, fpop, fpopkm1, jpremin
      REAL(wp) :: sum_phosp, ave_phosp

      ! Iron
      REAL(wp) :: jfe_uptake, jpofe, jfe_recycle, kfe_eq_lig
      REAL(wp) :: fpofe, fpofekm1, jferemin, jfe_ads_inorg, jfe_ads_org
      REAL(wp) :: dum5, dum2, dum3

      REAL(wp), POINTER, DIMENSION(:)     :: j_pop, feprime, dum4
      REAL(wp), POINTER, DIMENSION(:,:,:) :: irr_inst, irr_mix, j_po4, j_dop, j_fed, j_oxy
      REAL(wp), POINTER, DIMENSION(:,:,:) :: xnegtr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: wrk1, wrk2, wrk3, wrk4
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_bling')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL wrk_alloc( jpk, j_pop, feprime, dum4 )
      CALL wrk_alloc( jpi, jpj, jpk, xnegtr )
      CALL wrk_alloc( jpi, jpj, jpk, irr_inst, irr_mix, j_po4, j_dop, j_fed, j_oxy )
      CALL wrk_alloc( jpi, jpj, jpk, wrk1, wrk2, wrk3, wrk4 )

      CALL trc_opt_bling (kt, irr_inst, irr_mix)  ! optical model (red and blue wavelengths)

      !print*, 'bling...', rfact
      !jk=20
      !write(*,'(I3,3F14.7)'), kt, trb(60,60,jk,jpPO4_bling)*1.d6, &
      !                            trn(60,60,jk,jpPO4_bling)*1.d6, &
      !                            tra(60,60,jk,jpPO4_bling)*1.d6

      DO jj=1, jpj
          DO ji=1, jpi

              dum4(:) = 0.d0

              DO jk=1, jpk

               ! ----------------------------------------------------------
               ! negative trophic variables DO not contribute to the fluxes
               ! ----------------------------------------------------------

               f_po4 = MAX( 0.e0, trn(ji,jj,jk,jpPO4_bling) )
               f_dop = MAX( 0.e0, trn(ji,jj,jk,jpDOP_bling) )
               f_fer = MAX( 0.e0, trn(ji,jj,jk,jpFed_bling) )
               f_oxy = trn(ji,jj,jk,jpOxy_bling)

               ! ----------------------------------------------------------
               ! TEMPERATURE DEPENDENCE
               ! NB: The temperature effect of Eppley (1972) is used instead 
               !     of that in Geider et al (1997) for both simplicity and 
               !     to incorporate combined effects on uptake, incorporation
               !     into organic matter and photorespiration.  Values of PCmax
               !     are normalized to 0C rather than 20C in Geider et al.(1997)
               ! ----------------------------------------------------------

               expkT=EXP(kappa_eppley*tsn(ji,jj,jk,jp_tem))
 
               ! ----------------------------------------------------------
               ! Phytoplankton are assumed to grow according to the general properties 
               ! described in Geider (1997). This formulation gives a biomass-specific 
               ! growthrate as a function of light, nutrient limitation, and 
               ! temperature. We modify this relationship slightly here, as described 
               ! below, and also use the assumption of steady state growth vs. loss to 
               ! derive a simple relationship between growth rate, biomass and uptake.
               ! ----------------------------------------------------------
               ! First, we calculate the limitation terms for PO4 and Fe, and the 
               ! Fe-limited Chl:C maximum.
               ! The light-saturated maximal photosynthesis rate term (pc_m) is simply 
               ! the product of a prescribed maximal photosynthesis rate (pc_0), the 
               ! Eppley temperature dependence, and a Liebig limitation (the minimum
               ! of Michaelis-Menton PO4-limitation, or iron-limitation).
               ! The iron limitation term has a lower limit of def_fe_min 
               ! and is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
               ! so that it approaches 1 as fed approaches infinity. Thus, 
               ! it's of comparable magnitude to the PO4 limitation term.
               !
               ! Fe limitation acts in two additional mechanisms:
               ! 1. By reducing the maximum achievable Chl:C ratio 
               ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
               ! approach a prescribed minimum Chl:C (thetamin) under extreme
               ! Fe-limitation.
               ! 2. By reducing alpha (the initial slope of the P-I curve) under Fe-
               ! limitation.
               ! ----------------------------------------------------------

               ! Iron uptake

               fe2p_up = fe2p_max * f_fer / (kfe + f_fer)
               def_fe  = def_fe_min + (1.d0-def_fe_min)*fe2p_up/(kfe2p_up+fe2p_up)*(kfe2p_up+fe2p_max)/fe2p_max 

               ! Phosphate uptake
               po4_up = f_po4 /( kpo4 + f_po4 )

               ! Maximum production 
               pc_m = pc_0 * expkT * MIN(po4_up,def_fe) 

               ! Iron limitation on photosyntesis machinery
               thetamax_fe=thetamax_lo + (thetamax_hi - thetamax_lo)*def_fe
               alpha_chl  =alpha_min   + (alpha_max   - alpha_min  )*def_fe

               !-----------------------------------------------------------------------
               ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
               ! calculated. This requires a prescribed quantum yield, alpha.
               ! The iron deficiency term is included here as a multiplier of the 
               ! thetamax_fe to represent the importance of Fe in forming chlorophyll
               ! accessory antennae, which do not affect the Chl:C but still affect the
               ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
               ! 2004).
               !-----------------------------------------------------------------------

               ! I_k
               irrk = pc_m / (epsln + alpha_chl*thetamax_fe) + irr_mem(ji,jj,jk)*0.5d0

               !-----------------------------------------------------------------------
               ! Now we can calculate the carbon-specific photosynthesis rate, pc_tot.
               !-----------------------------------------------------------------------

               pc_tot = pc_m*(1.d0-EXP(-irr_mix(ji,jj,jk)/(irrk+epsln)))

               !-----------------------------------------------------------------------
               ! Next, we account for the maintenance effort that phytoplankton must 
               ! exert in order to combat decay. This is prescribed as a fraction of the
               ! light-saturated photosynthesis rate, resp_frac. The result of this is 
               ! to set a level of energy availability below which net growth (and 
               ! therefore nutrient uptake) is zero, given by resp_frac * pc_m.
               !-----------------------------------------------------------------------

               ! Net total production
               mu = MAX (0.d0,pc_tot-resp_frac*pc_m)

               !-----------------------------------------------------------------------
               ! We now must convert this net carbon-specific growth rate to nutrient 
               ! uptake rates, the quantities we are interested in. Since we have no 
               ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
               ! to calculate an implicit biomass from the uptake rate through the  
               ! application of a simple idealized grazing law. This has the effect of 
               ! reducing uptake in low growth-rate regimes and increasing uptake in 
               ! high growth-rate regimes - essentially a non-linear amplification of 
               ! the growth rate variability. The result is:
               !-----------------------------------------------------------------------

               ! Biomass
               mulamb0expkT = mu/(lambda0*expkT)
               biomass_p_ts = p_star*mulamb0expkT*(1.d0+(mulamb0expkT)**2)

               IF (kt==nittrc000) biomass_p(ji,jj,jk)=biomass_p_ts

               biomass_p(ji,jj,jk) =   biomass_p(ji,jj,jk) &
                                    + (biomass_p_ts-biomass_p(ji,jj,jk))*MIN(1.d0,gam_biomass*rdt)*tmask(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! We can now use the diagnostic biomass to calculate the chlorophyll
               ! concentration:
               !-----------------------------------------------------------------------

               ! Chl:C ration
               theta   = thetamax_fe / (1.d0 + (thetamax_fe*alpha_chl*irr_mem(ji,jj,jk))/(2.d0*pc_m+epsln) )

               ! Chl biomass
               chl_dia = biomass_p(ji,jj,jk) * c2p * 12.011e+6 * theta * tmask(ji,jj,jk) 
               chl_bling(ji,jj,jk) = MAX(chl_min, chl_dia)

               !--------------------------------------------------
               ! PHOSPHORUS CYCLE
               !--------------------------------------------------
               ! The uptake of nutrients is assumed to contribute to the growth of
               ! phytoplankton, which subsequently die and are consumed by heterotrophs.
               ! This can involve the transfer of nutrient elements between many
               ! organic pools, both particulate and dissolved, with complex histories.
               ! We take a simple approach here, partitioning the total uptake into two
               ! fractions - sinking and non-sinking - as a function of temperature, 
               ! following Dunne et al. (2005). 
               ! Then, the non-sinking fraction is further subdivided, such that the 
               ! majority is recycled instantaneously to the inorganic nutrient pool,
               ! representing the fast turnover of labile dissolved organic matter via
               ! the microbial loop, and the remainder is converted to semi-labile
               ! dissolved organic matter. Iron and phosphorus are treated identically 
               ! for the first step, but all iron is recycled instantaneously in the
               ! second step (i.e. there is no dissolved organic iron pool).
               !-----------------------------------------------------------------------

               ! Phosphorous uptake flux
               jp_uptake=biomass_p(ji,jj,jk)*mu

               ! 
               frac_pop=(phi_sm+phi_lg*(mulamb0expkT)**2) / (1+(mulamb0expkT)**2) * EXP(kappa_remin*tsn(ji,jj,jk,jp_tem))

               !
               j_pop(jk)=frac_pop*jp_uptake

               !-----------------------------------------------------------------------
               ! Then the remainder is divided between instantaneously recycled and
               ! long-lived dissolved organic matter,
               !-----------------------------------------------------------------------
               !
               jdop=phi_dop*(jp_uptake-j_pop(jk))

               ! 
               jp_recycle=jp_uptake-j_pop(jk)-jdop

               !---------------------------------------------------------------------
               ! IRON
               !---------------------------------------------------------------------
               ! Iron is then taken up as a function of PO4 uptake and iron limitation,
               ! with a maximum Fe:P uptake ratio of fe2p_max:
               !-----------------------------------------------------------------------

               ! Iron uptake flux
               jfe_uptake =jp_uptake*fe2p_up
               jpofe      =frac_pop*jfe_uptake
               jfe_recycle=jfe_uptake-jpofe

               !-----------------------------------------------------------------------
               ! Calculate free and inorganically associated iron concentration for
               ! scavenging.
               ! We assume that there is a 
               ! spectrum of iron ligands present in seawater, with varying binding
               ! strengths and whose composition varies with light and iron 
               ! concentrations. For example, photodissocation of ligand complexes 
               ! occurs under bright light, weakening the binding strength 
               ! (e.g. Barbeau et al., Nature 2001), while at very low iron 
               ! concentrations (order kfe_eq_lig_femin), siderophores are thought
               ! to be produced as a response to extreme
               ! iron stress.
               ! In anoxic waters, iron should be reduced, and therefore mostly 
               ! immune to scavenging. Easiest way to do this is to skip the feprime
               ! calculation if oxygen is less than 0.
               !-----------------------------------------------------------------------

               if (f_oxy > oxy_min ) then
                 dum5       = irr_inst(ji,jj,jk)**2/(kfe_eq_lig_irr**2+irr_inst(ji,jj,jk)**2)
                 dum2       = max(  0.d0, min( 1.d0, 1.2d0*(f_fer-kfe_eq_lig_femin)/(epsln+f_fer) )  )
                 kfe_eq_lig = kfe_eq_lig_max -(kfe_eq_lig_max-kfe_eq_lig_min)*dum5*dum2
                 feprime(jk)= 1.d0 + kfe_eq_lig*(felig_bkg - f_fer)
                 feprime(jk)= (-feprime(jk) + sqrt(feprime(jk)**2 + 4.d0*kfe_eq_lig*f_fer) )/(2.d0*kfe_eq_lig)
               else
                 feprime(jk) = 0.d0
               endif

               jfe_ads_inorg = min( 0.5d0/rdt, kfe_inorg*sqrt(feprime(jk)) )*feprime(jk)

               dum4(jk) = jpofe + jfe_ads_inorg

               !--------------------------------------------------
               ! COMPUTE TRENDS w/o remineralization processes
               !--------------------------------------------------
               !
               j_po4(ji,jj,jk) =   jp_recycle + gamma_dop*f_dop - jp_uptake
               j_dop(ji,jj,jk) = - gamma_dop*f_dop + phi_dop*(jp_uptake-j_pop(jk))
               j_fed(ji,jj,jk) =   jfe_recycle-jfe_uptake-jfe_ads_inorg

               ! checks
               ! wrk1(ji,jj,jk)=mulamb0expkT 
               ! wrk2(ji,jj,jk)=irrk
               ! wrk3(ji,jj,jk)=f_po4+f_dop
               ! wrk4(ji,jj,jk)=0.d0
            ENDDO

            !-----------------------------------------------------------------------
            ! SINKING AND REMINERALIZATION
            !-----------------------------------------------------------------------
            ! Calculate the remineralization lengthscale matrix, zremin, a function 
            ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
            ! then  increases linearly with depth.
            ! The remineralization rate is a function of oxygen concentrations,
            ! following a Holling type 2 dependence, decreasing to a minimum value
            ! of remin_min. This is ad hoc, following work by Bianchi, Sarmiento,
            ! Galbraith and Kwon (unpublished).
            !-----------------------------------------------------------------------
            ! In general, the flux at the bottom of a grid cell should equal
            ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
            ! where Ft is the flux at the top, and prod*dz is the integrated 
            ! production of new sinking particles within the layer.
            ! Since Ft=0 in the first layer,
            !---------------------------------------------------------------------

            ! k=1: surface layer
            zzz=fse3t(ji,jj,1)
            IF (zzz .lt. wsink0_z) THEN
               wsink=wsink0
            ELSE
               wsink=wsink0+wsink_acc*(zzz-wsink0_z)
            ENDIF

            ! Phosphorous
            f_po4   = MAX( 0.e0, trn(ji,jj,1,jpPO4_bling) )
            oxy_up =(f_po4*o2_2_p)**2 / (koxy**2 + (f_po4*o2_2_p)**2)
            zremin =(gamma_pop*oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)
         
            fpop   = j_pop(1)*fse3t(ji,jj,1)/(1.d0+fse3t(ji,jj,1)*zremin) 
            jpremin=(j_pop(1)*fse3t(ji,jj,1)-fpop)/(epsln+fse3t(ji,jj,1))

            !-----------------------------------------------------------------------
            ! Now, calculate the Fe adsorption using this fpop:
            ! The absolute first order rate constant is calculated from the 
            ! concentration of organic particles, after Parekh et al. (2005). Never
            !  allowed to be greater than 1/2dt for numerical stability.
            !-----------------------------------------------------------------------

            !Iron
            dum3        = (fpop*c2p*12.011d0/(epsln+wsink))**(0.58)
            jfe_ads_org = min (0.5d0/rdt, kfe_org*dum3)*feprime(1)
            dum4(1)     = dum4(1)+jfe_ads_org
            fpofe       = dum4(1)*fse3t(ji,jj,1)/(1.d0+fse3t(ji,jj,1)*zremin)
            jferemin    = dum4(1)*( 1.d0-1.d0/(1.d0+fse3t(ji,jj,1)*zremin) )

            ! Add remineralization terms to trends
            j_po4(ji,jj,1)=j_po4(ji,jj,1)+(1.d0-phi_dop)*jpremin
            j_dop(ji,jj,1)=j_dop(ji,jj,1)+       phi_dop*jpremin
            j_fed(ji,jj,1)=j_fed(ji,jj,1)+              jferemin-jfe_ads_org

            !checks
            !wrk1(ji,jj,1)= j_po4(ji,jj,1)+j_dop(ji,jj,1)
            !wrk2(ji,jj,1)= -(j_pop(1)-jpremin)
            !dum1(ji,jj,1)=dum1(ji,jj,1)+wrk1(ji,jj,1)
            !if (ji==130 .AND. jj==110) write(*,*) 1, fse3t(ji,jj,1), zzz, wsink

            ! k=2:NK: rest of the water column
            DO jk=2, jpk

               fpopkm1 =fpop
               fpofekm1=fpofe
    
               zzz=zzz+fse3t(ji,jj,jk)
               IF (zzz .lt. wsink0_z) THEN
                  wsink=wsink0
               ELSE
                  wsink=wsink0+wsink_acc*(zzz-wsink0_z)
               ENDIF

               ! Phosphorous
               f_po4 = MAX( 0.e0, trn(ji,jj,jk,jpPO4_bling) )
               oxy_up =(f_po4*o2_2_p)**2 / (koxy**2 + (f_po4*o2_2_p)**2)
               zremin =(gamma_pop*oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)

               fpop   =(fpopkm1+j_pop(jk)*fse3t(ji,jj,jk))/(1.d0+fse3t(ji,jj,jk)*zremin) 
               jpremin=(fpopkm1+j_pop(jk)*fse3t(ji,jj,jk)-fpop)/(epsln+fse3t(ji,jj,jk))

               ! Iron
               dum3        = (fpop*c2p*12.011d0/(epsln+wsink))**(0.58)
               jfe_ads_org = min (0.5d0/rdt, kfe_org*dum3)*feprime(jk)
               dum4(jk)    = dum4(jk)+jfe_ads_org
               fpofe       = (fpofekm1 + dum4(jk)*fse3t(ji,jj,jk))/(1.d0+fse3t(ji,jj,jk)*zremin)
               jferemin    = ( 1.d0-1.d0/(1.d0+fse3t(ji,jj,jk)*zremin) ) * (fpofekm1/fse3t(ji,jj,jk)+dum4(jk))

               ! Save fPOP and fPOFe at the bottom grid cell to compute bottom fluxes
               IF (jk==mbkt(ji,jj)) THEN
                  fpop_b (ji,jj)=fpop
                  fpofe_b(ji,jj)=fpofe
               ENDIF
           
               ! Add remineralization terms to trends
               j_po4(ji,jj,jk)=j_po4(ji,jj,jk)+(1.d0-phi_dop)*jpremin
               j_dop(ji,jj,jk)=j_dop(ji,jj,jk)+      phi_dop *jpremin
               j_fed(ji,jj,jk)=j_fed(ji,jj,jk)+              jferemin-jfe_ads_org

               !checks
               !wrk1(ji,jj,jk)= j_po4(ji,jj,jk)+j_dop(ji,jj,jk)
               !wrk2(ji,jj,jk)= -(j_pop(jk)-jpremin)
               !dum1(ji,jj,jk)=dum1(ji,jj,jk)+wrk1(ji,jj,jk)
               !wrk3(ji,jj,jk)=0.d0
               !wrk4(ji,jj,jk)=0.d0
               !if (ji==130 .AND. jj==110) write(*,*) jk, fse3t(ji,jj,jk), zzz, wsink

            ENDDO

            !-----------------------------------------------------------------------
            ! OXYGEN
            !-----------------------------------------------------------------------
            ! Assuming constant P:O ratio.
            ! Optional prevention of negative oxygen (does not conserve ocean 
            ! redox potential) or alternatively it can be allowed to go negative, 
            ! keeping track of an implicit nitrate deficit 
            ! plus sulfate reduction.
            !-----------------------------------------------------------------------

            DO jk=1, jpk
               IF ( (ln_prev_o2lt0) .and. (f_oxy<oxy_min) ) then
                   j_oxy(ji,jj,jk)=0.d0
               ELSE
                   j_oxy(ji,jj,jk)=-oxy2p*j_po4(ji,jj,jk)
               ENDIF
            ENDDO

         ENDDO
      ENDDO

      !checks
      !dum1(:,:,:)=1.d0*tmask(:,:,:)-1.d0
      !write(*,*), 'BLING...', kt, MINVAL(j_po4), MAXVAL(j_po4), MINVAL(j_dop), MAXVAL(j_dop)
      !j_po4(:,:,:)=2.d-6/(rfact*14.)
      !j_dop(:,:,:)=3.d-6/(rfact*14.)

      tra(:,:,:,jpPO4_bling) = tra(:,:,:,jpPO4_bling) + j_po4(:,:,:)*rfact
      tra(:,:,:,jpDOP_bling) = tra(:,:,:,jpDOP_bling) + j_dop(:,:,:)*rfact
      tra(:,:,:,jpFed_bling) = tra(:,:,:,jpFed_bling) + j_fed(:,:,:)*rfact
      tra(:,:,:,jpOxy_bling) = tra(:,:,:,jpOxy_bling) + j_oxy(:,:,:)*rfact

      ! add bottom fluxes
      CALL trc_sed_bling

      ! total mass of phosphate
      sum_phosp = glob_sum ( ( trn(:,:,:,jpPO4_bling) + trn(:,:,:,jpDOP_bling) ) *cvol(:,:,:)  )
      ave_phosp = sum_phosp / glob_sum(cvol(:,:,:))

      !test if concentrations fall below 0
      xnegtr(:,:,:) = 1.e0
      DO jn = jp_bling0, jp_bling1
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ( trn(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN 
                     ztra             = ABS(  ( trn(ji,jj,jk,jn) - rtrn ) &
                                            / ( tra(ji,jj,jk,jn) + rtrn ) )
                     xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                  ENDIF
              END DO
            END DO
         END DO
      END DO

      ! Prgonostic tracer fields
      DO jn = jp_bling0, jp_bling1
         trn(:,:,:,jn) = trn(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
      END DO

      tra(:,:,:,jp_bling0:jp_bling1) = 0.e0

      ! Copy new arrays to trb (tracer fields before) to avoid the Asselin filter
      DO jn=jp_bling0, jp_bling1
         trb(:,:,:,jn)=trn(:,:,:,jn)
      ENDDO

      !checks
      !jk=20
      !write(*,'(I3,3F14.7)'), kt, trb(60,60,jk,jpPO4_bling)*1.d6, &
      !                            trn(60,60,jk,jpPO4_bling)*1.d6, &
      !                            tra(60,60,jk,jpPO4_bling)*1.d6
      !IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
      !   WRITE(charout, FMT="('bio')")
      !   CALL prt_ctl_trc_info(charout)
      !   CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      !ENDIF
      !write(*,'(I3,4F14.7)'), kt,sumblpo4, sumbloxy, sumblfed, sum_phosp
      !write(*,*), kt,sumblpo4, sumbloxy, sumblfed, sum_phosp
      wrk4(:,:,10)=e1t(:,:)*e2t(:,:)*tmask_i(:,:)
      wrk4(:,:, 9)=mbkt(:,:)
      wrk4(:,:, 8)=tmask_i(:,:)
      wrk4(100,10,11)=sumblpo4
      wrk4(100,12,11)=sumbloxy
      wrk4(100,13,11)=sumblfed
      wrk4(100,14,11)=sum_phosp
      wrk4(100,15,11)=ave_phosp
      wrk1(:,:,:)=cvol(:,:,:)
      wrk2(:,:,:)=fse3t(:,:,:)

      IF( lk_iomput ) THEN
            CALL iom_put( "wrk1", wrk1(:,:,:)  )  
            CALL iom_put( "wrk2", wrk2(:,:,:)  )  
            CALL iom_put( "wrk3", wrk3(:,:,:)  )  
            CALL iom_put( "wrk4", wrk4(:,:,:)  )  
      ENDIF

      CALL wrk_dealloc( jpk, j_pop, feprime, dum4 )
      CALL wrk_dealloc( jpi, jpj, jpk, xnegtr )
      CALL wrk_dealloc( jpi, jpj, jpk, irr_inst, irr_mix, j_po4, j_dop, j_fed, j_oxy )
      CALL wrk_dealloc( jpi, jpj, jpk, wrk1, wrk2, wrk3, wrk4 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_bling')
      !
   END SUBROUTINE trc_sms_bling

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                       No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_bling             ! Empty routine
      WRITE(*,*) 'trc_sms_bling: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_bling
#endif

   !!======================================================================
END MODULE trcsms_bling
