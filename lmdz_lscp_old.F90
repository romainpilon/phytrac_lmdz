! $Id: lmdz_lscp_old.F90 4674 2023-09-05 10:13:50Z fhourdin $
!
!
MODULE lmdz_lscp_old
CONTAINS
SUBROUTINE fisrtilp(klon,klev,dtime,paprs,pplay,t,q,ptconv,ratqs, &
     d_t, d_q, d_ql, d_qi, rneb,rneblsvol,radliq, rain, snow,          &
     pfrac_impa, pfrac_nucl, pfrac_1nucl,               &
     frac_impa, frac_nucl, beta,                        &
     prfl, psfl, rhcl, zqta, fraca,                     &
     ztv, zpspsk, ztla, zthl, iflag_cld_th,             &
     iflag_ice_thermo,                                  &
     cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv)


  !
  USE icefrac_lsc_mod ! compute ice fraction (JBM 3/14)
  USE lmdz_cloudth, only : cloudth, cloudth_v3, cloudth_v6

  USE lmdz_lscp_ini, ONLY: prt_level, lunout
  USE lmdz_lscp_ini, ONLY : fl_cor_ebil 
  USE lmdz_lscp_ini, ONLY: iflag_t_glace,t_glace_min, t_glace_max, exposant_glace
  USE lmdz_lscp_ini, ONLY : seuil_neb, rain_int_min, iflag_evap_prec, iflag_oldbug_fisrtilp,a_tr_sca
  USE lmdz_lscp_ini, ONLY: iflag_cloudth_vert, iflag_rain_incloud_vol 
  USE lmdz_lscp_ini, ONLY: coef_eva, coef_eva_i, ffallv_lsc, ffallv_con
  USE lmdz_lscp_ini, ONLY: cld_tau_lsc, cld_tau_con, cld_lc_lsc, cld_lc_con
  USE lmdz_lscp_ini, ONLY: reevap_ice, iflag_bergeron, iflag_fisrtilp_qsat, iflag_pdf 



  IMPLICIT none
  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS)
  ! Date: le 20 mars 1995
  ! Objet: condensation et precipitation stratiforme.
  !        schema de nuage
  ! Fusion de fisrt (physique sursaturation, P. LeVan K. Laval)
  !             et ilp (il pleut, L. Li)
  ! Principales parties:
  ! P0> Thermalisation des precipitations venant de la couche du dessus
  ! P1> Evaporation de la precipitation (qui vient du niveau k+1)
  ! P2> Formation du nuage (en k) 
  ! P2.A.0> Calcul des grandeurs nuageuses une pdf en creneau
  ! P2.A.1> Avec les nouvelles PDFs, calcul des grandeurs nuageuses pour 
  ! les valeurs de T et Q initiales
  ! P2.A.2> Prise en compte du couplage entre eau condensee et T.
  ! P2.A.3> Calcul des valeures finales associees a la formation des nuages
  ! P2.B> Nuage "tout ou rien"
  ! P2.C> Prise en compte de la Chaleur latente apres formation nuage
  ! P3> Formation de la precipitation (en k)
  !======================================================================
  ! JLD: 
  ! * Routine probablement fausse (au moins incoherente) si thermcep = .false.
  ! * fl_cor_ebil doit etre > 0 ; 
  !   fl_cor_ebil= 0 pour reproduire anciens bugs
  !======================================================================
  include "YOMCST.h"
  include "YOETHF.h"
  include "FCTTRE.h"
  !
  ! Principaux inputs:
  !
  REAL, INTENT(IN)                              :: dtime  ! intervalle du temps (s)
  INTEGER, INTENT(IN)                           :: klon, klev
  REAL, DIMENSION(klon,klev+1),    INTENT(IN)   :: paprs  ! pression a inter-couche
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: pplay  ! pression au milieu de couche
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: t      ! temperature (K)
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: q      ! humidite specifique (kg/kg)
  LOGICAL, DIMENSION(klon,klev),   INTENT(IN)   :: ptconv ! points ou le schema de conv. prof. est actif
  INTEGER,                         INTENT(IN)   :: iflag_cld_th
  INTEGER,                         INTENT(IN)   :: iflag_ice_thermo
  !
  ! Inputs lies aux thermiques
  ! 
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: ztv
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: zqta, fraca 
  REAL, DIMENSION(klon,klev),      INTENT(IN)   :: zpspsk, ztla
  REAL, DIMENSION(klon,klev),      INTENT(INOUT)   :: zthl
  !
  !  Input/output
  REAL, DIMENSION(klon,klev),      INTENT(INOUT):: ratqs  ! determine la largeur de distribution de vapeur
  !
  ! Principaux outputs:
  !
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: d_t  ! incrementation de la temperature (K)
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: d_q  ! incrementation de la vapeur d'eau
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: d_ql ! incrementation de l'eau liquide
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: d_qi ! incrementation de l'eau glace
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: rneb, rneblsvol ! fraction nuageuse
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: beta ! taux de conversion de l'eau cond
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: radliq ! eau liquide utilisee dans rayonnements
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: rhcl ! humidite relative en ciel clair
  REAL, DIMENSION(klon),           INTENT(OUT)  :: rain
  REAL, DIMENSION(klon),           INTENT(OUT)  :: snow
  REAL, DIMENSION(klon,klev+1),    INTENT(OUT)  :: prfl
  REAL, DIMENSION(klon,klev+1),    INTENT(OUT)  :: psfl 

  !AA
  ! Coeffients de fraction lessivee : pour OFF-LINE
  !
  REAL, DIMENSION(klon,klev),      INTENT(INOUT)  :: pfrac_nucl
  REAL, DIMENSION(klon,klev),      INTENT(INOUT)  :: pfrac_1nucl
  REAL, DIMENSION(klon,klev),      INTENT(INOUT)  :: pfrac_impa
  !
  ! Fraction d'aerosols lessivee par impaction et par nucleation
  ! POur ON-LINE
  !
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: frac_impa
  REAL, DIMENSION(klon,klev),      INTENT(OUT)  :: frac_nucl
  REAL, DIMENSION(klon,klev), INTENT(OUT) :: cloudth_sth,cloudth_senv
  REAL, DIMENSION(klon,klev), INTENT(OUT) :: cloudth_sigmath,cloudth_sigmaenv
  !AA
  ! --------------------------------------------------------------------------------
  !
  ! Options du programme:
  !

  REAL :: smallestreal

  INTEGER, PARAMETER :: ninter=5 ! sous-intervals pour la precipitation
  LOGICAL, PARAMETER :: cpartiel=.TRUE. ! condensation partielle
  REAL, PARAMETER :: t_coup=234.0
  REAL, PARAMETER :: DDT0=.01
  REAL, PARAMETER :: ztfondue=278.15
  ! --------------------------------------------------------------------------------
  !
  ! Variables locales:
  !
  INTEGER :: i, k, n, kk
  REAL :: qsl, qsi
  REAL :: zct      ,zcl
  INTEGER :: ncoreczq  
  REAL, DIMENSION(klon,klev) :: ctot,ctot_vol
  REAL, DIMENSION(klon) :: zqs, zdqs, zdqsdT_raw, Tbef,qlbef,DT
  REAL :: zdelta, zcor, zcvm5  
  REAL ::num,denom

  LOGICAL, DIMENSION(klon) :: lognormale,convergence
  LOGICAL :: ice_thermo
  INTEGER, DIMENSION(klon) :: n_i
  INTEGER :: iter
  REAL :: cste

  REAL, DIMENSION(klon) :: zpdf_sig,zpdf_k,zpdf_delta, Zpdf_a,zpdf_b,zpdf_e1,zpdf_e2, qcloud
  REAL :: erf   
 
  REAL :: zqev, zqevt, zqev0,zqevi, zqevti, zdelq 
  REAL, DIMENSION(klon) :: zrfl(klon), zrfln(klon), zrflclr(klon), zrflcld(klon), d_zrfl_clr_cld(klon), d_zifl_clr_cld(klon), d_zrfl_cld_clr(klon), d_zifl_cld_clr(klon) 

  REAL, DIMENSION(klon) :: zifl, zifln, ziflclr, ziflcld, zoliq, zcond, zq, zqn, zoliqp, zoliqi, zt
! JBM (3/14) nexpo is replaced by exposant_glace
! REAL nexpo ! exponentiel pour glace/eau
! INTEGER, PARAMETER :: nexpo=6
  INTEGER :: exposant_glace_old
  REAL :: t_glace_min_old, ztot
  REAL, DIMENSION(klon) ::  zdz,zrho , zrhol, zfice,zneb,znebprecip
  REAL :: zchau      ,zfroi      
  REAL, DIMENSION(klon) :: znebprecipclr, znebprecipcld, tot_zneb, tot_znebn, d_tot_zneb, d_znebprecip_clr_cld, d_znebprecip_cld_clr, dzfice
  REAL :: zmelt, zpluie, zice
  REAL :: zsolid
!!!!
!  Variables pour Bergeron
  REAL :: zcp, coef1, DeltaT, Deltaq, Deltaqprecl
  REAL, DIMENSION(klon) :: zqpreci, zqprecl
! Variable pour conservation enegie des precipitations
  REAL, DIMENSION(klon) :: zmqc
  !
  LOGICAL, SAVE :: appel1er=.TRUE.
  !$OMP THREADPRIVATE(appel1er)
  !
! iflag_oldbug_fisrtilp=0 enleve le BUG par JYG : tglace_min -> tglace_max
! iflag_oldbug_fisrtilp=1 ajoute le BUG
  !---------------------------------------------------------------
  !
  ! Fonctions en ligne:
  !
  REAL ::  fallvs,fallvc, zzz ! Vitesse de chute pour cristaux de glace
                     ! (Heymsfield & Donner, 1990)
  fallvc (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_con
  fallvs (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_lsc
  !
  !---------------------------------------------------------------
  !AA Variables traceurs:
  !AA  Provisoire !!! Parametres alpha du lessivage
  !AA  A priori on a 4 scavenging # possibles
  !
  ! Variables intermediaires
  !
  REAL :: zalpha_tr, zfrac_lessi
  REAL, DIMENSION(klon) :: zprec_cond
  !AA
! RomP >>> 15 nov 2012
! RomP <<<
  REAL, DIMENSION(klon) :: zmair
  REAL :: zcpair, zcpeau
  !     Pour la conversion eau-neige
  REAL, DIMENSION(klon) :: zlh_solid
  REAL :: zm_solid


  !
  !ym
!CR: pour iflag_ice_thermo=2, on active que la convection
!  ice_thermo = iflag_ice_thermo .GE. 1

  
  znebprecip(:)=0.

!<LTP
  smallestreal=1.e-9
  znebprecipclr(:)=0.
  znebprecipcld(:)=0.
!>LTP

  ice_thermo = (iflag_ice_thermo .EQ. 1).OR.(iflag_ice_thermo .GE. 3)
  zdelq=0.0
  ctot_vol(1:klon,1:klev)=0.0
  rneblsvol(1:klon,1:klev)=0.0

  if (prt_level>9)write(lunout,*)'NUAGES4 A. JAM'
  IF (appel1er) THEN
     WRITE(lunout,*) 'fisrtilp, ninter:', ninter
     WRITE(lunout,*) 'fisrtilp, cpartiel:', cpartiel
     WRITE(lunout,*) 'FISRTILP VERSION LUDO'
     
     IF (ABS(dtime/REAL(ninter)-360.0).GT.0.001) THEN
        WRITE(lunout,*) 'fisrtilp: Ce n est pas prevu, voir Z.X.Li', dtime
        WRITE(lunout,*) 'Je prefere un sous-intervalle de 6 minutes'
        !         CALL abort
     ENDIF
     appel1er = .FALSE.
     !
     !cdir collapse
     DO k = 1, klev
        DO i = 1, klon
           pfrac_nucl(i,k)=1.
           pfrac_1nucl(i,k)=1.
           pfrac_impa(i,k)=1.
           beta(i,k)=0.  !RomP initialisation
        ENDDO
     ENDDO

  ENDIF          !  test sur appel1er
  !
  !MAf Initialisation a 0 de zoliq
  !      DO i = 1, klon
  !         zoliq(i)=0.
  !      ENDDO 
  ! Determiner les nuages froids par leur temperature
  !  nexpo regle la raideur de la transition eau liquide / eau glace.
  !
!CR: on est oblige de definir des valeurs fisrt car les valeurs de newmicro ne sont pas les memes par defaut
  IF (iflag_t_glace.EQ.0) THEN
!   ztglace = RTT - 15.0
    t_glace_min_old = RTT - 15.0
    !AJ<
    IF (ice_thermo) THEN
!     nexpo = 2
      exposant_glace_old = 2
    ELSE
!     nexpo = 6
      exposant_glace_old = 6
    ENDIF
    
  ENDIF
  
!!  RLVTT = 2.501e6 ! pas de redefinition des constantes physiques (jyg)
!!  RLSTT = 2.834e6 ! pas de redefinition des constantes physiques (jyg)
!>AJ
  !cc      nexpo = 1
  !
  ! Initialiser les sorties:
  !
  !cdir collapse
  DO k = 1, klev+1
     DO i = 1, klon
        prfl(i,k) = 0.0
        psfl(i,k) = 0.0
     ENDDO
  ENDDO

  !cdir collapse

  DO k = 1, klev
     DO i = 1, klon
        d_t(i,k) = 0.0
        d_q(i,k) = 0.0
        d_ql(i,k) = 0.0
        d_qi(i,k) = 0.0
        rneb(i,k) = 0.0
        radliq(i,k) = 0.0
        frac_nucl(i,k) = 1. 
        frac_impa(i,k) = 1. 
     ENDDO
  ENDDO
  DO i = 1, klon
     rain(i) = 0.0
     snow(i) = 0.0
     zoliq(i)=0.
     !     ENDDO
     !
     ! Initialiser le flux de precipitation a zero
     !
     !     DO i = 1, klon
     zrfl(i) = 0.0
     zifl(i) = 0.0
!<LTP
     zrflclr(i) = 0.0
     ziflclr(i) = 0.0
     zrflcld(i) = 0.0
     ziflcld(i) = 0.0
     tot_zneb(i) = 0.0
     tot_znebn(i) = 0.0
     d_tot_zneb(i) = 0.0
!>LTP

     zneb(i) = seuil_neb
  ENDDO
  !
  !
  !AA Pour plus de securite 

  zalpha_tr   = 0.
  zfrac_lessi = 0.

  !AA==================================================================
  !
  ncoreczq=0
  ! BOUCLE VERTICALE (DU HAUT VERS LE BAS)
  !
  DO k = klev, 1, -1
     !
     !AA===============================================================
     !
     ! Initialisation temperature et vapeur
     DO i = 1, klon
        zt(i)=t(i,k)
        zq(i)=q(i,k)
     ENDDO
     !
     ! ----------------------------------------------------------------
     ! P0> Thermalisation des precipitations venant de la couche du dessus
     ! ----------------------------------------------------------------
     ! Calculer la varition de temp. de l'air du a la chaleur sensible
     ! transporter par la pluie. On thermalise la pluie avec l'air de la couche.
     ! Cette quantite de pluie qui est thermalisee, et devra continue a l'etre lors
     ! des differentes transformations thermodynamiques. Cette masse d'eau doit
     ! donc etre ajoute a l'humidite de la couche lorsque l'on calcule la variation
     ! de l'enthalpie  de la couche avec la temperature
     ! Variables calculees ou modifiees:
     !   -  zt: temperature de la cocuhe
     !   - zmqc: masse de precip qui doit etre thermalisee
     !
     IF(k.LE.klev-1) THEN         
        DO i = 1, klon
           !IM
           zmair(i)=(paprs(i,k)-paprs(i,k+1))/RG
           ! il n'y a pas encore d'eau liquide ni glace dans la maiille, donc zq suffit
           zcpair=RCPD*(1.0+RVTMP2*zq(i))
           zcpeau=RCPD*RVTMP2
         if (fl_cor_ebil .GT. 0) then
           ! zmqc: masse de precip qui doit etre thermalisee avec l'air de la couche atm 
           ! pour s'assurer que la precip arrivant au sol aura bien la temperature de la 
           ! derniere couche
           zmqc(i) = (zrfl(i)+zifl(i))*dtime/zmair(i)
           ! t(i,k+1)+d_t(i,k+1): nouvelle temp de la couche au dessus
           zt(i) = ( (t(i,k+1)+d_t(i,k+1))*zmqc(i)*zcpeau + zcpair*zt(i) ) &
                 / (zcpair + zmqc(i)*zcpeau)
         else ! si on maintient les anciennes erreurs
           zt(i) = ( (t(i,k+1)+d_t(i,k+1))*zrfl(i)*dtime*zcpeau &
                + zmair(i)*zcpair*zt(i) ) &
                / (zmair(i)*zcpair + zrfl(i)*dtime*zcpeau)
         end if
        ENDDO
     ELSE  ! IF(k.LE.klev-1) 
        DO i = 1, klon 
           zmair(i)=(paprs(i,k)-paprs(i,k+1))/RG 
           zmqc(i) = 0.
        ENDDO
     ENDIF ! end IF(k.LE.klev-1)
     !
     ! ----------------------------------------------------------------
     ! P1> Calcul de l'evaporation de la precipitation
     ! ----------------------------------------------------------------
     ! On evapore une partie des precipitations venant de la maille du dessus.
     ! On calcule l'evaporation et la sublimation des precipitations, jusqu'a
     ! ce que la fraction de cette couche qui est sous le nuage soit saturee.
     ! Variables calculees ou modifiees:
     !   - zrfl et zifl: flux de precip liquide et glace
     !   - zq, zt: humidite et temperature de la cocuhe
     !   - zmqc: masse de precip qui doit etre thermalisee
     !
     IF (iflag_evap_prec>=1) THEN
        DO i = 1, klon
!          S'il y a des precipitations
           IF (zrfl(i)+zifl(i).GT.0.) THEN
              ! Calcul du qsat
              IF (thermcep) THEN
                 zdelta=MAX(0.,SIGN(1.,RTT-zt(i)))
                 zqs(i)= R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
                 zqs(i)=MIN(0.5,zqs(i))
                 zcor=1./(1.-RETV*zqs(i))
                 zqs(i)=zqs(i)*zcor
              ELSE
                 IF (zt(i) .LT. t_coup) THEN
                    zqs(i) = qsats(zt(i)) / pplay(i,k)
                 ELSE
                    zqs(i) = qsatl(zt(i)) / pplay(i,k)
                 ENDIF
              ENDIF
           ENDIF ! (zrfl(i)+zifl(i).GT.0.)
        ENDDO
!AJ<

       IF (.NOT. ice_thermo) THEN
        DO i = 1, klon
!          S'il y a des precipitations
           IF (zrfl(i)+zifl(i).GT.0.) THEN
                ! Evap max pour ne pas saturer la fraction sous le nuage
                ! Evap max jusqu'à atteindre la saturation dans la partie
                ! de la maille qui est sous le nuage de la couche du dessus
                !!! On ne tient compte de cette fraction que sous une seule
                !!! couche sous le nuage
                zqev = MAX (0.0, (zqs(i)-zq(i))*zneb(i) )
             ! Ajout de la prise en compte des precip a thermiser
             ! avec petite reecriture 
             if  (fl_cor_ebil .GT. 0) then ! nouveau
                ! Calcul de l'evaporation du flux de precip herite
                !   d'au-dessus
                zqevt = coef_eva * (1.0-zq(i)/zqs(i)) * SQRT(zrfl(i)) &
                     * zmair(i)/pplay(i,k)*zt(i)*RD
                zqevt = MAX(0.0,MIN(zqevt,zrfl(i))) * dtime/zmair(i)

                ! Seuil pour ne pas saturer la fraction sous le nuage
                zqev = MIN (zqev, zqevt)
                ! Nouveau flux de precip
                zrfln(i) = zrfl(i) - zqev*zmair(i)/dtime
                ! Aucun flux liquide pour T < t_coup, on reevapore tout.
                IF (zt(i) .LT. t_coup.and.reevap_ice) THEN
                  zrfln(i)=0.
                  zqev = (zrfl(i)-zrfln(i))/zmair(i)*dtime
                END IF
                ! Nouvelle vapeur 
                zq(i) = zq(i) + zqev
                zmqc(i) = zmqc(i)-zqev
                ! Nouvelle temperature (chaleur latente)
                zt(i) = zt(i) - zqev &
                     * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
!!JLD debut de partie a supprimer a terme
            else ! if  (fl_cor_ebil .GT. 0)
                ! Calcul de l'evaporation du flux de precip herite
                !   d'au-dessus
                zqevt = coef_eva * (1.0-zq(i)/zqs(i)) * SQRT(zrfl(i)) &
                     * (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
                zqevt = MAX(0.0,MIN(zqevt,zrfl(i))) &
                     * RG*dtime/(paprs(i,k)-paprs(i,k+1))
                ! Seuil pour ne pas saturer la fraction sous le nuage
                zqev = MIN (zqev, zqevt)
                ! Nouveau flux de precip
                zrfln(i) = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1)) &
                     /RG/dtime
                ! Aucun flux liquide pour T < t_coup
                IF (zt(i) .LT. t_coup.and.reevap_ice) zrfln(i)=0.
                ! Nouvelle vapeur 
                zq(i) = zq(i) - (zrfln(i)-zrfl(i)) &
                     * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
                ! Nouvelle temperature (chaleur latente)
                zt(i) = zt(i) + (zrfln(i)-zrfl(i)) &
                     * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                     * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
              end if ! if  (fl_cor_ebil .GT. 0)
!!JLD fin de partie a supprimer a terme
                zrfl(i) = zrfln(i)
                zifl(i) = 0.
           ENDIF ! (zrfl(i)+zifl(i).GT.0.)
        ENDDO
!
       ELSE ! (.NOT. ice_thermo) 
!      ================================
!      Avec thermodynamique de la glace
!      ================================
        DO i = 1, klon


!AJ<
!        S'il y a des precipitations
         IF (zrfl(i)+zifl(i).GT.0.) THEN

        !LTP<
        !On ne tient compte que du flux de précipitation en ciel clair dans le calcul de l'évaporation.
                IF (iflag_evap_prec==4) THEN
                        zrfl(i) = zrflclr(i)
                        zifl(i) = ziflclr(i)
                ENDIF
        
        !>LTP

         IF (iflag_evap_prec==1) THEN
            znebprecip(i)=zneb(i)
         ELSE
            znebprecip(i)=MAX(zneb(i),znebprecip(i))
         ENDIF
         
         IF (iflag_evap_prec==4) THEN
        ! Evap max pour ne pas saturer toute la maille
         zqev0 = MAX (0.0, zqs(i)-zq(i))
         ELSE 
        ! Evap max pour ne pas saturer la fraction sous le nuage
         zqev0 = MAX (0.0, (zqs(i)-zq(i))*znebprecip(i) )
         ENDIF

         !JAM
         ! On differencie qsat pour l'eau et la glace
         ! Si zdelta=1. --> glace
         ! Si zdelta=0. --> eau liquide
        
         ! Calcul du qsat par rapport a l'eau liquide
         qsl= R2ES*FOEEW(zt(i),0.)/pplay(i,k)
         qsl= MIN(0.5,qsl)
         zcor= 1./(1.-RETV*qsl)
         qsl= qsl*zcor
         
         ! Calcul de l'evaporation du flux de precip venant du dessus
         ! Formulation en racine du flux de precip
         ! (Klemp & Wilhelmson, 1978; Sundqvist, 1988)
         IF (iflag_evap_prec==3) THEN
         zqevt = znebprecip(i)*coef_eva*(1.0-zq(i)/qsl) &
              *SQRT(zrfl(i)/max(1.e-4,znebprecip(i))) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
!<LTP
         ELSE IF (iflag_evap_prec==4) THEN
         zqevt = znebprecipclr(i)*coef_eva*(1.0-zq(i)/qsl) &
              *SQRT(zrfl(i)/max(1.e-8,znebprecipclr(i))) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
!>LTP
        ELSE 
         zqevt = 1.*coef_eva*(1.0-zq(i)/qsl)*SQRT(zrfl(i)) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
         ENDIF


         zqevt = MAX(0.0,MIN(zqevt,zrfl(i))) &
              *RG*dtime/(paprs(i,k)-paprs(i,k+1)) 
         
         ! Calcul du qsat par rapport a la glace
         qsi= R2ES*FOEEW(zt(i),1.)/pplay(i,k)
         qsi= MIN(0.5,qsi)
         zcor= 1./(1.-RETV*qsi)
         qsi= qsi*zcor

         ! Calcul de la sublimation du flux de precip solide herite
         !   d'au-dessus
         IF (iflag_evap_prec==3) THEN
         zqevti = znebprecip(i)*coef_eva*(1.0-zq(i)/qsi) &
              *SQRT(zifl(i)/max(1.e-4,znebprecip(i))) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
!<LTP
         ELSE IF (iflag_evap_prec==4) THEN
         zqevti = znebprecipclr(i)*coef_eva*(1.0-zq(i)/qsi) &
              *SQRT(zifl(i)/max(1.e-8,znebprecipclr(i))) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
!>LTP
         ELSE
         zqevti = 1.*coef_eva*(1.0-zq(i)/qsi)*SQRT(zifl(i)) &
              *(paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
         ENDIF
         zqevti = MAX(0.0,MIN(zqevti,zifl(i))) &
              *RG*dtime/(paprs(i,k)-paprs(i,k+1))   

       
        !JAM
        ! Limitation de l'evaporation. On s'assure qu'on ne sature pas 
        ! la fraction de la couche sous le nuage sinon on repartit zqev0 
        ! en conservant la proportion liquide / glace
     
         IF (zqevt+zqevti.GT.zqev0) THEN
            zqev=zqev0*zqevt/(zqevt+zqevti)
            zqevi=zqev0*zqevti/(zqevt+zqevti)
         ELSE
!JLD je ne comprends pas les lignes ci-dessous. On repartit les precips
!       liquides et solides meme si on ne sature pas la couche. 
!       A mon avis, le test est inutile, et il faudrait tout remplacer par:
!            zqev=zqevt
!            zqevi=zqevti
             IF (zqevt+zqevti.GT.0.) THEN
                zqev=MIN(zqev0*zqevt/(zqevt+zqevti),zqevt)
                zqevi=MIN(zqev0*zqevti/(zqevt+zqevti),zqevti)
             ELSE
             zqev=0.
             zqevi=0.
             ENDIF
         ENDIF

         ! Nouveaux flux de precip liquide et solide 
         zrfln(i) = Max(0.,zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1)) &
                                 /RG/dtime)
         zifln(i) = Max(0.,zifl(i) - zqevi*(paprs(i,k)-paprs(i,k+1)) &
                                 /RG/dtime)

         ! Mise a jour de la vapeur, temperature et flux de precip
         zq(i) = zq(i) - (zrfln(i)+zifln(i)-zrfl(i)-zifl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
       if (fl_cor_ebil .GT. 0) then ! avec correction thermalisation des precips
         zmqc(i) = zmqc(i) + (zrfln(i)+zifln(i)-zrfl(i)-zifl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
         zt(i) = zt(i) + (zrfln(i)-zrfl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                  * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i))) &
                  + (zifln(i)-zifl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                  * RLSTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
       else ! sans correction thermalisation des precips
         zt(i) = zt(i) + (zrfln(i)-zrfl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                  * RLVTT/RCPD/(1.0+RVTMP2*zq(i)) &
                  + (zifln(i)-zifl(i)) &
                  * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                  * RLSTT/RCPD/(1.0+RVTMP2*zq(i))
       end if
        ! Nouvelles vaeleurs des precips liquides et solides 
         zrfl(i) = zrfln(i)
         zifl(i) = zifln(i)

!<LTP
        IF (iflag_evap_prec==4) THEN
                zrflclr(i) = zrfl(i)
                ziflclr(i) = zifl(i)    
                IF(zrflclr(i) + ziflclr(i) .LE. 0) THEN
                        znebprecipclr(i) = 0.
                ENDIF   
                zrfl(i) = zrflclr(i) + zrflcld(i)
                zifl(i) = ziflclr(i) + ziflcld(i) 
        ENDIF
!>LTP        


!CR ATTENTION: deplacement de la fonte de la glace
!jyg : Bug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! jyg
!!!        zmelt = ((zt(i)-273.15)/(ztfondue-273.15))**2  !!!!!!!!! jyg
!jyg : Bug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! jyg
           zmelt = ((zt(i)-273.15)/(ztfondue-273.15))             ! jyg
           ! fraction de la precip solide qui est fondue
           zmelt = MIN(MAX(zmelt,0.),1.)
           ! Fusion de la glace
!<LTP
           IF (iflag_evap_prec==4) THEN
                   zrflclr(i)=zrflclr(i)+zmelt*ziflclr(i)
                   zrflcld(i)=zrflcld(i)+zmelt*ziflcld(i)
                   zrfl(i)=zrflclr(i)+zrflcld(i)
!>LTP        
           ELSE
                   zrfl(i)=zrfl(i)+zmelt*zifl(i)
           ENDIF
           if (fl_cor_ebil .LE. 0) then
             ! the following line should not be here. Indeed, if zifl is modified
             ! now, zifl(i)*zmelt is no more the amount of ice that has melt
             ! and therefore the change in temperature computed below is wrong
             zifl(i)=zifl(i)*(1.-zmelt)
           end if
           ! Chaleur latente de fusion
        if (fl_cor_ebil .GT. 0) then ! avec correction thermalisation des precips
           zt(i)=zt(i)-zifl(i)*zmelt*(RG*dtime)/(paprs(i,k)-paprs(i,k+1)) &
                      *RLMLT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
        else ! sans correction thermalisation des precips
           zt(i)=zt(i)-zifl(i)*zmelt*(RG*dtime)/(paprs(i,k)-paprs(i,k+1)) &
                      *RLMLT/RCPD/(1.0+RVTMP2*zq(i))
        end if
           if (fl_cor_ebil .GT. 0) then ! correction bug, deplacement ligne precedente
!<LTP
             IF (iflag_evap_prec==4) THEN
                   ziflclr(i)=ziflclr(i)*(1.-zmelt)
                   ziflcld(i)=ziflcld(i)*(1.-zmelt)
                   zifl(i)=ziflclr(i)+ziflcld(i)
!>LTP        
             ELSE
                   zifl(i)=zifl(i)*(1.-zmelt)
             ENDIF
           end if

           ELSE
              ! Si on n'a plus de pluies, on reinitialise a 0 la farcion
              ! sous nuageuse utilisee pour la pluie.
              znebprecip(i)=0.
           ENDIF ! (zrfl(i)+zifl(i).GT.0.)
        ENDDO

      ENDIF ! (.NOT. ice_thermo) 
     
     ! ----------------------------------------------------------------
     ! Fin evaporation de la precipitation
     ! ----------------------------------------------------------------
     ENDIF ! (iflag_evap_prec>=1)
     !
     ! Calculer Qs et L/Cp*dQs/dT:
     !
     IF (thermcep) THEN
        DO i = 1, klon
           zdelta = MAX(0.,SIGN(1.,RTT-zt(i)))
           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
       if  (fl_cor_ebil .GT. 0) then ! nouveau
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
       else    
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i))
       endif
           zqs(i) = R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
           zqs(i) = MIN(0.5,zqs(i))
           zcor = 1./(1.-RETV*zqs(i))
           zqs(i) = zqs(i)*zcor
           zdqs(i) = FOEDE(zt(i),zdelta,zcvm5,zqs(i),zcor)
           zdqsdT_raw(i) = zdqs(i)*  &
         &         RCPD*(1.0+RVTMP2*zq(i)) / (RLVTT*(1.-zdelta) + RLSTT*zdelta)
        ENDDO
     ELSE
        DO i = 1, klon
           IF (zt(i).LT.t_coup) THEN
              zqs(i) = qsats(zt(i))/pplay(i,k)
              zdqs(i) = dqsats(zt(i),zqs(i))
           ELSE
              zqs(i) = qsatl(zt(i))/pplay(i,k)
              zdqs(i) = dqsatl(zt(i),zqs(i))
           ENDIF
        ENDDO
     ENDIF
     !
     ! Determiner la condensation partielle et calculer la quantite
     ! de l'eau condensee:
     !
!verification de la valeur de iflag_fisrtilp_qsat pour iflag_ice_thermo=1
!       if ((iflag_ice_thermo.eq.1).and.(iflag_fisrtilp_qsat.ne.0)) then
!         write(*,*) " iflag_ice_thermo==1 requires iflag_fisrtilp_qsat==0", &
!        " but iflag_fisrtilp_qsat=",iflag_fisrtilp_qsat, ". Might as well stop here."
!         stop
!       endif

     ! ----------------------------------------------------------------
     ! P2> Formation du nuage
     ! ----------------------------------------------------------------
     ! Variables calculees:
     !   rneb  : fraction nuageuse
     !   zcond : eau condensee moyenne dans la maille.
     !   rhcl: humidite relative ciel-clair
     !   zt : temperature de la maille
     ! ----------------------------------------------------------------
     !
     IF (cpartiel) THEN
        ! -------------------------
        ! P2.A> Nuage fractionnaire
        ! -------------------------
        !
        !   Calcul de l'eau condensee et de la fraction nuageuse et de l'eau
        !   nuageuse a partir des PDF de Sandrine Bony.
        !   rneb  : fraction nuageuse
        !   zqn   : eau totale dans le nuage
        !   zcond : eau condensee moyenne dans la maille.
        !  on prend en compte le réchauffement qui diminue la partie
        ! condensee
        !
        !   Version avec les raqts

        ! ----------------------------------------------------------------
        ! P2.A.0> Calcul des grandeurs nuageuses une pdf en creneau
        ! ----------------------------------------------------------------
        if (iflag_pdf.eq.0) then

           ! version creneau de (Li, 1998)
           do i=1,klon
              zdelq = min(ratqs(i,k),0.99) * zq(i)
              rneb(i,k) = (zq(i)+zdelq-zqs(i)) / (2.0*zdelq)
              zqn(i) = (zq(i)+zdelq+zqs(i))/2.0
           enddo

        else !  if (iflag_pdf.eq.0)
           ! ----------------------------------------------------------------
           ! P2.A.1> Avec les nouvelles PDFs, calcul des grandeurs nuageuses pour 
           ! les valeurs de T et Q initiales
           ! ----------------------------------------------------------------
           do i=1,klon
              if(zq(i).lt.1.e-15) then
                 ncoreczq=ncoreczq+1
                 zq(i)=1.e-15
              endif
           enddo

           if (iflag_cld_th>=5) then

              if (iflag_cloudth_vert<=2) then
               call cloudth(klon,klev,k,ztv, &
                   zq,zqta,fraca, &
                   qcloud,ctot,zpspsk,paprs,pplay,ztla,zthl, &
                   ratqs,zqs,t, &
                   cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv)

              elseif (iflag_cloudth_vert>=3 .and. iflag_cloudth_vert<=5) then
               call cloudth_v3(klon,klev,k,ztv, &
                   zq,zqta,fraca, &
                   qcloud,ctot,ctot_vol,zpspsk,paprs,pplay,ztla,zthl, &
                   ratqs,zqs,t, &
                   cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv)
              !----------------------------------
              !Version these Jean Jouhaud, Decembre 2018 
              !----------------------------------              
              elseif (iflag_cloudth_vert==6) then
               call cloudth_v6(klon,klev,k,ztv, &
                   zq,zqta,fraca, &
                   qcloud,ctot,ctot_vol,zpspsk,paprs,pplay,ztla,zthl, &
                   ratqs,zqs,t, &
                   cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv)

              endif
              do i=1,klon
                 rneb(i,k)=ctot(i,k)
                 rneblsvol(i,k)=ctot_vol(i,k)
                 zqn(i)=qcloud(i)
              enddo

           endif

           if (iflag_cld_th <= 4) then
              lognormale = .true.
           elseif (iflag_cld_th >= 6) then
              ! lognormale en l'absence des thermiques
              lognormale = fraca(:,k) < 1e-10
           else
              ! Dans le cas iflag_cld_th=5, on prend systématiquement la
              ! bi-gaussienne
              lognormale = .false.
           end if

!CR: variation de qsat avec T en presence de glace ou non
!initialisations
           do i=1,klon
              DT(i) = 0.
              n_i(i)=0
              Tbef(i)=zt(i)
              qlbef(i)=0.
           enddo

        ! ----------------------------------------------------------------
        ! P2.A.2> Prise en compte du couplage entre eau condensee et T.
        ! Calcul des grandeurs nuageuses en tenant compte de l'effet de
        ! la condensation sur T, et donc sur qsat et sur les grandeurs nuageuses
        ! qui en dependent. Ce changement de temperature est provisoire, et
        ! la valeur definitive sera calcule plus tard.
        ! Variables calculees: 
        !   rneb : nebulosite
        !   zcond: eau condensee en moyenne dans la maille
        ! note JLD: si on n'a pas de pdf lognormale, ce qui se passe ne me semble
        ! pas clair, il n'y a probablement pas de prise en compte de l'effet de
        ! T sur qsat
        ! ----------------------------------------------------------------

!Boucle iterative: ATTENTION, l'option -1 n'est plus activable ici
           if (iflag_fisrtilp_qsat.ge.0) then
             ! Iteration pour condensation avec variation de qsat(T) 
             ! -----------------------------------------------------
             do iter=1,iflag_fisrtilp_qsat+1
                
               do i=1,klon
!                 do while ((abs(DT(i)).gt.DDT0).or.(n_i(i).eq.0))
                 ! !! convergence = .true. tant que l'on n'a pas converge !!
                 ! ------------------------------
                 convergence(i)=abs(DT(i)).gt.DDT0
                 if ((convergence(i).or.(n_i(i).eq.0)).and.lognormale(i)) then
                 ! si on n'a pas converge
                 !
                 ! P2.A.2.1> Calcul de la fraction nuageuse et de la quantite d'eau condensee
                 ! ---------------------------------------------------------------
                 ! Variables calculees: 
                 ! rneb : nebulosite
                 ! zqn : eau condensee, dans le nuage (in cloud water content)
                 ! zcond: eau condensee en moyenne dans la maille
                 ! rhcl: humidite relative ciel-clair
                 !
                 Tbef(i)=Tbef(i)+DT(i) ! nouvelle temperature
                 if (.not.ice_thermo) then   
                    zdelta = MAX(0.,SIGN(1.,RTT-Tbef(i)))
                 else
                    if (iflag_t_glace.eq.0) then
                    zdelta = MAX(0.,SIGN(1.,t_glace_min_old-Tbef(i)))
                    else if (iflag_t_glace.ge.1) then
                       if (iflag_oldbug_fisrtilp.EQ.0) then
                          zdelta = MAX(0.,SIGN(1.,t_glace_max-Tbef(i)))
                       else
                          !avec bug : zdelta = MAX(0.,SIGN(1.,t_glace_min-Tbef(i)))
                          zdelta = MAX(0.,SIGN(1.,t_glace_min-Tbef(i)))
                       endif
                    endif
                 endif
                 ! Calcul de rneb, qzn et zcond pour les PDF lognormales
                 zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
               if (fl_cor_ebil .GT. 0) then
                 zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
               else
                 zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i))
               end if
                 zqs(i) = R2ES*FOEEW(Tbef(i),zdelta)/pplay(i,k)
                 zqs(i) = MIN(0.5,zqs(i))
                 zcor = 1./(1.-RETV*zqs(i))
                 zqs(i) = zqs(i)*zcor
                 zdqs(i) = FOEDE(Tbef(i),zdelta,zcvm5,zqs(i),zcor)
                 zpdf_sig(i)=ratqs(i,k)*zq(i)
                 zpdf_k(i)=-sqrt(log(1.+(zpdf_sig(i)/zq(i))**2))
                 zpdf_delta(i)=log(zq(i)/zqs(i))
                 zpdf_a(i)=zpdf_delta(i)/(zpdf_k(i)*sqrt(2.))
                 zpdf_b(i)=zpdf_k(i)/(2.*sqrt(2.))
                 zpdf_e1(i)=zpdf_a(i)-zpdf_b(i)
                 zpdf_e1(i)=sign(min(abs(zpdf_e1(i)),5.),zpdf_e1(i))
                 zpdf_e1(i)=1.-erf(zpdf_e1(i))
                 zpdf_e2(i)=zpdf_a(i)+zpdf_b(i)
                 zpdf_e2(i)=sign(min(abs(zpdf_e2(i)),5.),zpdf_e2(i))
                 zpdf_e2(i)=1.-erf(zpdf_e2(i))
             
                 if (zpdf_e1(i).lt.1.e-10) then
                    rneb(i,k)=0.
                    zqn(i)=zqs(i)
                 else
                    rneb(i,k)=0.5*zpdf_e1(i)
                    zqn(i)=zq(i)*zpdf_e2(i)/zpdf_e1(i)
                 endif

                 ! If vertical heterogeneity, change fraction by volume as well
                 if (iflag_cloudth_vert>=3) then
                   ctot_vol(i,k)=rneb(i,k)
                   rneblsvol(i,k)=ctot_vol(i,k)
                 endif

                 endif !convergence

               enddo ! boucle en i

                 ! P2.A.2.2> Calcul APPROCHE de la variation de temperature DT 
                 !         due a la condensation.
                 ! ---------------------------------------------------------------
                 ! Variables calculees: 
                 ! DT : variation de temperature due a la condensation 

                 if (.not. ice_thermo) then
                 ! --------------------------
                 do i=1,klon
                 if ((convergence(i).or.(n_i(i).eq.0)).and.lognormale(i)) then

                 qlbef(i)=max(0.,zqn(i)-zqs(i))
               if (fl_cor_ebil .GT. 0) then
                 num=-Tbef(i)+zt(i)+rneb(i,k)*RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))*qlbef(i)
               else
                 num=-Tbef(i)+zt(i)+rneb(i,k)*RLVTT/RCPD/(1.0+RVTMP2*zq(i))*qlbef(i)
               end if
                 denom=1.+rneb(i,k)*zdqs(i)
                 DT(i)=num/denom
                 n_i(i)=n_i(i)+1
                 endif
                 enddo

                 else ! if (.not. ice_thermo)
                 ! --------------------------
                 if (iflag_t_glace.ge.1) then
                 CALL icefrac_lsc(klon,zt(:),pplay(:,k)/paprs(:,1),zfice(:))
                 endif

                 do i=1,klon
                 if ((convergence(i).or.(n_i(i).eq.0)).and.lognormale(i)) then
                 
                 if (iflag_t_glace.eq.0) then
                    zfice(i) = 1.0 - (Tbef(i)-t_glace_min_old) / (RTT-t_glace_min_old)
                    zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
                    zfice(i) = zfice(i)**exposant_glace_old
                    dzfice(i)= exposant_glace_old * zfice(i)**(exposant_glace_old-1) &
          &                     / (t_glace_min_old - RTT)
                 endif
                 
                 if (iflag_t_glace.ge.1.and.zfice(i)>0.) then
                 dzfice(i)= exposant_glace * zfice(i)**(exposant_glace-1) &
          &                    / (t_glace_min - t_glace_max)
                 endif
               
                 if ((zfice(i).eq.0).or.(zfice(i).eq.1)) then
                    dzfice(i)=0.
                 endif

                 if (zfice(i).lt.1) then
                    cste=RLVTT
                 else
                    cste=RLSTT
                 endif

                 qlbef(i)=max(0.,zqn(i)-zqs(i))
               if (fl_cor_ebil .GT. 0) then
                 num = -Tbef(i)+zt(i)+rneb(i,k)*((1-zfice(i))*RLVTT &
           &          +zfice(i)*RLSTT)/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))*qlbef(i)
                 denom = 1.+rneb(i,k)*((1-zfice(i))*RLVTT+zfice(i)*RLSTT)/cste*zdqs(i) &
                         -(RLSTT-RLVTT)/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))*rneb(i,k)    &
           &               *qlbef(i)*dzfice(i)
               else
                 num = -Tbef(i)+zt(i)+rneb(i,k)*((1-zfice(i))*RLVTT &
           &         +zfice(i)*RLSTT)/RCPD/(1.0+RVTMP2*zq(i))*qlbef(i)
                 denom = 1.+rneb(i,k)*((1-zfice(i))*RLVTT+zfice(i)*RLSTT)/cste*zdqs(i) &
                         -(RLSTT-RLVTT)/RCPD/(1.0+RVTMP2*zq(i))*rneb(i,k)*qlbef(i)*dzfice(i)
               end if
                 DT(i)=num/denom
                 n_i(i)=n_i(i)+1

                 endif ! fin convergence
                 enddo ! fin boucle i

                 endif !ice_thermo

             enddo ! iter=1,iflag_fisrtilp_qsat+1
             ! Fin d'iteration pour condensation avec variation de qsat(T) 
             ! -----------------------------------------------------------
           endif !  if (iflag_fisrtilp_qsat.ge.0)
     ! ----------------------------------------------------------------
     ! Fin de P2.A.2> la prise en compte du couplage entre eau condensee et T
     ! ----------------------------------------------------------------

        endif ! iflag_pdf

!        if (iflag_fisrtilp_qsat.eq.-1) then 
!------------------------------------------
!CR: ATTENTION option fausse mais a existe:
! pour la re-activer, prendre iflag_fisrtilp_qsat=0 et
! activer les lignes suivantes: 
       IF (1.eq.0) THEN 
       DO i=1,klon
           IF (rneb(i,k) .LE. 0.0) THEN
              zqn(i) = 0.0
              rneb(i,k) = 0.0
              zcond(i) = 0.0
              rhcl(i,k)=zq(i)/zqs(i)
           ELSE IF (rneb(i,k) .GE. 1.0) THEN
              zqn(i) = zq(i)
              rneb(i,k) = 1.0                 
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))/(1+zdqs(i))
              rhcl(i,k)=1.0
           ELSE
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)/(1+zdqs(i))
              rhcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
           ENDIF
       ENDDO
       ENDIF
!------------------------------------------

!        ELSE
        ! ----------------------------------------------------------------
        ! P2.A.3> Calcul des valeures finales associees a la formation des nuages
        ! Variables calculees: 
        !   rneb : nebulosite
        !   zcond: eau condensee en moyenne dans la maille
        !   zq : eau vapeur dans la maille
        !   zt : temperature de la maille
        !   rhcl: humidite relative ciel-clair
        ! ----------------------------------------------------------------
        !
        ! Bornage de l'eau in-cloud (zqn) et de la fraction nuageuse (rneb) 
        ! Calcule de l'eau condensee moyenne dans la maille (zcond),
        ! et de l'humidite relative ciel-clair (rhcl)
        DO i=1,klon
           IF (rneb(i,k) .LE. 0.0) THEN
              zqn(i) = 0.0
              rneb(i,k) = 0.0
              zcond(i) = 0.0
              rhcl(i,k)=zq(i)/zqs(i)
           ELSE IF (rneb(i,k) .GE. 1.0) THEN
              zqn(i) = zq(i)
              rneb(i,k) = 1.0
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))
              rhcl(i,k)=1.0
           ELSE
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)
              rhcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
           ENDIF
        ENDDO

       
        ! If vertical heterogeneity, change fraction by volume as well
        if (iflag_cloudth_vert>=3) then
          ctot_vol(1:klon,k)=min(max(ctot_vol(1:klon,k),0.),1.)
          rneblsvol(1:klon,k)=ctot_vol(1:klon,k)
        endif

!        ENDIF

     ELSE ! de IF (cpartiel) 
        ! -------------------------
        ! P2.B> Nuage "tout ou rien"
        ! -------------------------
        ! note JLD: attention, rhcl non calcule. Ca peut avoir des consequences?
        DO i = 1, klon
           IF (zq(i).GT.zqs(i)) THEN
              rneb(i,k) = 1.0
           ELSE
              rneb(i,k) = 0.0
           ENDIF
           zcond(i) = MAX(0.0,zq(i)-zqs(i))/(1.+zdqs(i))
        ENDDO
     ENDIF ! de IF (cpartiel) 
     !
     ! Mise a jour vapeur d'eau
     ! -------------------------
     DO i = 1, klon
        zq(i) = zq(i) - zcond(i)
        !         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD
     ENDDO
!AJ<
     ! ------------------------------------
     ! P2.C> Prise en compte de la Chaleur latente apres formation nuage
     ! -------------------------------------
     ! Variable calcule:
     !   zt : temperature de la maille
     !
     IF (.NOT. ice_thermo) THEN
        if (iflag_fisrtilp_qsat.lt.1) then
           DO i = 1, klon
             zt(i) = zt(i) + zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
           ENDDO
        else if (iflag_fisrtilp_qsat.gt.0) then
           DO i= 1, klon
    if (fl_cor_ebil .GT. 0) then
             zt(i) = zt(i) + zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)+zcond(i)))
    else
             zt(i) = zt(i) + zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zcond(i)))
    end if
           ENDDO
        endif
     ELSE
         if (iflag_t_glace.ge.1) then
            CALL icefrac_lsc(klon,zt(:),pplay(:,k)/paprs(:,1),zfice(:))
         endif
         if (iflag_fisrtilp_qsat.lt.1) then
           DO i = 1, klon
! JBM: icefrac_lsc is now a function contained in icefrac_lsc_mod
!              zfice(i) = icefrac_lsc(zt(i), t_glace_min, &
!                                     t_glace_max, exposant_glace)
              if (iflag_t_glace.eq.0) then
                    zfice(i) = 1.0 - (zt(i)-t_glace_min_old) / (RTT-t_glace_min_old)
                    zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
                    zfice(i) = zfice(i)**exposant_glace_old
              endif
              zt(i) = zt(i) + (1.-zfice(i))*zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*zq(i)) &
                       +zfice(i)*zcond(i) * RLSTT/RCPD/(1.0+RVTMP2*zq(i))
           ENDDO
         else
           DO i=1, klon
! JBM: icefrac_lsc is now a function contained in icefrac_lsc_mod
!              zfice(i) = icefrac_lsc(zt(i), t_glace_min, &
!                                     t_glace_max, exposant_glace)
              if (iflag_t_glace.eq.0) then
                    zfice(i) = 1.0 - (zt(i)-t_glace_min_old) / (RTT-t_glace_min_old)
                    zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
                    zfice(i) = zfice(i)**exposant_glace_old
              endif
        if (fl_cor_ebil .GT. 0) then
              zt(i) = zt(i) + (1.-zfice(i))*zcond(i) & 
           &             * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)+zcond(i))) &
                      +zfice(i)*zcond(i) * RLSTT/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)+zcond(i)))
        else 
              zt(i) = zt(i) + (1.-zfice(i))*zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*(zq(i)+zcond(i))) &
                      +zfice(i)*zcond(i) * RLSTT/RCPD/(1.0+RVTMP2*(zq(i)+zcond(i)))
        end if
           ENDDO
         endif
!         print*,zt(i),zrfl(i),zifl(i),'temp1'
       ENDIF
!>AJ

     ! ----------------------------------------------------------------
     ! P3> Formation des precipitations
     ! ----------------------------------------------------------------
     !
     ! Partager l'eau condensee en precipitation et eau liquide nuageuse
     !

!<LTP

IF (iflag_evap_prec==4) THEN
        !Partitionnement des precipitations venant du dessus en précipitations nuageuses
        !et précipitations ciel clair

        !0) Calculate tot_zneb, la fraction nuageuse totale au-dessus du nuage
        !en supposant un recouvrement maximum aléatoire (voir Jakob and Klein, 2000)
        
        DO i=1, klon
                tot_znebn(i) = 1 - (1-tot_zneb(i))*(1 - max(rneb(i,k),zneb(i))) &
                        /(1-min(zneb(i),1-smallestreal))
                d_tot_zneb(i) = tot_znebn(i) - tot_zneb(i)
                tot_zneb(i) = tot_znebn(i)


                !1) Cloudy to clear air
                d_znebprecip_cld_clr(i) = znebprecipcld(i) - min(rneb(i,k),znebprecipcld(i)) 
                IF (znebprecipcld(i) .GT. 0) THEN
                        d_zrfl_cld_clr(i) = d_znebprecip_cld_clr(i)/znebprecipcld(i)*zrflcld(i) 
                        d_zifl_cld_clr(i) = d_znebprecip_cld_clr(i)/znebprecipcld(i)*ziflcld(i) 
                ELSE 
                        d_zrfl_cld_clr(i) = 0. 
                        d_zifl_cld_clr(i) = 0. 
                ENDIF

                !2) Clear to cloudy air
                d_znebprecip_clr_cld(i) = max(0., min(znebprecipclr(i), rneb(i,k) &
                        - d_tot_zneb(i) - zneb(i)))
                IF (znebprecipclr(i) .GT. 0) THEN
                        d_zrfl_clr_cld(i) = d_znebprecip_clr_cld(i)/znebprecipclr(i)*zrflclr(i) 
                        d_zifl_clr_cld(i) = d_znebprecip_clr_cld(i)/znebprecipclr(i)*ziflclr(i)
                ELSE
                        d_zrfl_clr_cld(i) = 0. 
                        d_zifl_clr_cld(i) = 0. 
                ENDIF 

                !Update variables
                znebprecipcld(i) = znebprecipcld(i) + d_znebprecip_clr_cld(i) - d_znebprecip_cld_clr(i)  
                znebprecipclr(i) = znebprecipclr(i) + d_znebprecip_cld_clr(i) - d_znebprecip_clr_cld(i) 
                zrflcld(i) = zrflcld(i) + d_zrfl_clr_cld(i) - d_zrfl_cld_clr(i) 
                ziflcld(i) = ziflcld(i) + d_zifl_clr_cld(i) - d_zifl_cld_clr(i)
                zrflclr(i) = zrflclr(i) + d_zrfl_cld_clr(i) - d_zrfl_clr_cld(i)
                ziflclr(i) = ziflclr(i) + d_zifl_cld_clr(i) - d_zifl_clr_cld(i)

        ENDDO
ENDIF

!>LTP



     ! Initialisation de zoliq (eau condensee moyenne dans la maille)
     DO i = 1, klon
        IF (rneb(i,k).GT.0.0) THEN
           zoliq(i) = zcond(i)
           zrho(i) = pplay(i,k) / zt(i) / RD
           zdz(i) = (paprs(i,k)-paprs(i,k+1)) / (zrho(i)*RG)
        ENDIF
     ENDDO
!AJ<
     IF (.NOT. ice_thermo) THEN
       IF (iflag_t_glace.EQ.0) THEN
         DO i = 1, klon
            IF (rneb(i,k).GT.0.0) THEN
               zfice(i) = 1.0 - (zt(i)-t_glace_min_old) / (273.13-t_glace_min_old)
               zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
               zfice(i) = zfice(i)**exposant_glace_old
!              zfice(i) = zfice(i)**nexpo
         !!      zfice(i)=0.
            ENDIF
         ENDDO
       ELSE ! of IF (iflag_t_glace.EQ.0)
         CALL icefrac_lsc(klon,zt(:),pplay(:,k)/paprs(:,1),zfice(:))
!         DO i = 1, klon
!            IF (rneb(i,k).GT.0.0) THEN
! JBM: icefrac_lsc is now a function contained in icefrac_lsc_mod
!              zfice(i) = icefrac_lsc(zt(i), t_glace_min, &
!                                     t_glace_max, exposant_glace)
!            ENDIF
!         ENDDO
       ENDIF
     ENDIF

     ! Calcul de radliq (eau condensee pour le rayonnement)
     ! Iteration pour realiser une moyenne de l'eau nuageuse lors de la precip
     ! Remarque: ce n'est donc pas l'eau restante en fin de precip mais une
     ! eau moyenne restante dans le nuage sur la duree du pas de temps qui est
     ! transmise au rayonnement;
     ! ----------------------------------------------------------------
     DO i = 1, klon
        IF (rneb(i,k).GT.0.0) THEN
           zneb(i) = MAX(rneb(i,k), seuil_neb)
     !      zt(i) = zt(i)+zcond(i)*zfice(i)*RLMLT/RCPD/(1.0+RVTMP2*zq(i))  
!           print*,zt(i),'fractionglace'
!>AJ
           radliq(i,k) = zoliq(i)/REAL(ninter+1)
        ENDIF
     ENDDO
     !
     DO n = 1, ninter
        DO i = 1, klon
           IF (rneb(i,k).GT.0.0) THEN
              zrhol(i) = zrho(i) * zoliq(i) / zneb(i)
              ! Initialization of zpluie and zice:
              zpluie=0
              zice=0
              IF (zneb(i).EQ.seuil_neb) THEN
                 ztot = 0.0
              ELSE
                 !  quantite d'eau a eliminer: zchau (Sundqvist, 1978)
                 !  meme chose pour la glace: zfroi (Zender & Kiehl, 1997)
                 if (ptconv(i,k)) then
                    zcl   =cld_lc_con
                    zct   =1./cld_tau_con
                    zfroi    = dtime/REAL(ninter)/zdz(i)*zoliq(i) &
                         *fallvc(zrhol(i)) * zfice(i)
                 else
                    zcl   =cld_lc_lsc
                    zct   =1./cld_tau_lsc
                    zfroi    = dtime/REAL(ninter)/zdz(i)*zoliq(i) &
                         *fallvs(zrhol(i)) * zfice(i)
                 endif

                 ! si l'heterogeneite verticale est active, on utilise 
                 ! la fraction volumique "vraie" plutot que la fraction
                 ! surfacique modifiee, qui est plus grande et reduit
                 ! sinon l'eau in-cloud de facon artificielle
                 if ((iflag_cloudth_vert>=3).AND.(iflag_rain_incloud_vol==1)) then
                    zchau    = zct   *dtime/REAL(ninter) * zoliq(i) &
                         *(1.0-EXP(-(zoliq(i)/ctot_vol(i,k)/zcl   )**2)) *(1.-zfice(i))
                 else
                    zchau    = zct   *dtime/REAL(ninter) * zoliq(i) &
                         *(1.0-EXP(-(zoliq(i)/zneb(i)/zcl   )**2)) *(1.-zfice(i))
                 endif
!AJ<
                 IF (.NOT. ice_thermo) THEN
                   ztot    = zchau + zfroi
                 ELSE
                   zpluie = MIN(MAX(zchau,0.0),zoliq(i)*(1.-zfice(i)))
                   zice = MIN(MAX(zfroi,0.0),zoliq(i)*zfice(i)) 
                   ztot    = zpluie    + zice
                 ENDIF
!>AJ
                 ztot    = MAX(ztot   ,0.0)
              ENDIF
              ztot    = MIN(ztot,zoliq(i))
!AJ<
     !         zoliqp = MAX(zoliq(i)*(1.-zfice(i))-1.*zpluie   , 0.0)
     !         zoliqi = MAX(zoliq(i)*zfice(i)-1.*zice   , 0.0)
!JLD : les 2 variables zoliqp et zoliqi crorresponent a des pseudo precip
!      temporaires et ne doivent pas etre calcule (alors qu'elles le sont
!      si iflag_bergeron <> 2
!      A SUPPRIMER A TERME
              zoliqp(i) = MAX(zoliq(i)*(1.-zfice(i))-1.*zpluie  , 0.0)
              zoliqi(i) = MAX(zoliq(i)*zfice(i)-1.*zice  , 0.0)
              zoliq(i) = MAX(zoliq(i)-ztot   , 0.0)
!>AJ
              radliq(i,k) = radliq(i,k) + zoliq(i)/REAL(ninter+1)
           ENDIF
        ENDDO  ! i = 1,klon
     ENDDO     ! n = 1,ninter

     ! ----------------------------------------------------------------
     !
     IF (.NOT. ice_thermo) THEN
       DO i = 1, klon
         IF (rneb(i,k).GT.0.0) THEN
           d_ql(i,k) = zoliq(i)
           zrfl(i) = zrfl(i)+ MAX(zcond(i)-zoliq(i),0.0) &
                * (paprs(i,k)-paprs(i,k+1))/(RG*dtime)
         ENDIF
       ENDDO
     ELSE
!
!CR&JYG<
! On prend en compte l'effet Bergeron dans les flux de precipitation :
! Si T < 0 C, alors les precipitations liquides sont converties en glace, ce qui
! provoque un accroissement de temperature DeltaT. L'effet de DeltaT sur le condensat
! et les precipitations est grossierement pris en compte en linearisant les equations 
! et en approximant le processus de precipitation liquide par un processus a seuil.
! On fait l'hypothese que le condensat nuageux n'est pas modifié dans cette opération.
! Le condensat precipitant liquide est supprime (dans la limite DeltaT<273-T).
! Le condensat precipitant solide est augmente.
! La vapeur d'eau est augmentee.
!
       IF ((iflag_bergeron .EQ. 2)) THEN
         DO i = 1, klon
           IF (rneb(i,k) .GT. 0.0) THEN
             zqpreci(i)=(zcond(i)-zoliq(i))*zfice(i)
             zqprecl(i)=(zcond(i)-zoliq(i))*(1.-zfice(i))
           if (fl_cor_ebil .GT. 0) then
             zcp=RCPD*(1.0+RVTMP2*(zq(i)+zmqc(i)+zcond(i)))
             coef1 = rneb(i,k)*RLSTT/zcp*zdqsdT_raw(i)
!            Calcul de DT si toute les precips liquides congelent
             DeltaT = RLMLT*zqprecl(i) / (zcp*(1.+coef1))
!            On ne veut pas que T devienne superieur a la temp. de congelation.
!            donc que Delta > RTT-zt(i
             DeltaT = max( min( RTT-zt(i), DeltaT) , 0. )
             zt(i)      = zt(i)      + DeltaT
!            Eau vaporisee du fait de l'augmentation de T
             Deltaq = rneb(i,k)*zdqsdT_raw(i)*DeltaT
!            on reajoute cette eau vaporise a la vapeur et on l'enleve des precips
             zq(i) = zq(i) +  Deltaq
!            Les 3 max si dessous prtotegent uniquement des erreurs d'arrondies
             zcond(i)   = max( zcond(i)- Deltaq, 0. )
!            precip liquide qui congele ou qui s'evapore
             Deltaqprecl = -zcp/RLMLT*(1.+coef1)*DeltaT
             zqprecl(i) = max( zqprecl(i) + Deltaqprecl, 0. )
!            bilan eau glacee
             zqpreci(i) = max (zqpreci(i) - Deltaqprecl - Deltaq, 0.)
           else ! if (fl_cor_ebil .GT. 0)
!            ancien calcul
             zcp=RCPD*(1.0+RVTMP2*(zq(i)+zcond(i)))
             coef1 = RLMLT*zdqs(i)/RLVTT
             DeltaT = max( min( RTT-zt(i), RLMLT*zqprecl(i)/zcp/(1.+coef1) ) , 0.)
             zqpreci(i) = zqpreci(i) + zcp/RLMLT*DeltaT
             zqprecl(i) = max( zqprecl(i) - zcp/RLMLT*(1.+coef1)*DeltaT, 0. )
             zcond(i)   = max( zcond(i)   - zcp/RLVTT*zdqs(i)*DeltaT, 0. )
             zq(i)      = zq(i)      + zcp/RLVTT*zdqs(i)*DeltaT
             zt(i)      = zt(i)      + DeltaT
           end if ! if (fl_cor_ebil .GT. 0)
           ENDIF  ! rneb(i,k) .GT. 0.0
         ENDDO
         DO i = 1, klon
           IF (rneb(i,k).GT.0.0) THEN
             d_ql(i,k) = (1-zfice(i))*zoliq(i)
             d_qi(i,k) = zfice(i)*zoliq(i)
!<LTP
             IF (iflag_evap_prec == 4) THEN
                zrflcld(i) = zrflcld(i)+zqprecl(i) &
                 *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)
                ziflcld(i) = ziflcld(i)+ zqpreci(i) &
                      *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)
                znebprecipcld(i) = rneb(i,k)
                zrfl(i) = zrflcld(i) + zrflclr(i)
                zifl(i) = ziflcld(i) + ziflclr(i)        
!>LTP
             ELSE
                zrfl(i) = zrfl(i)+ zqprecl(i) &
                 *(paprs(i,k)-paprs(i,k+1))/(RG*dtime) 
                zifl(i) = zifl(i)+ zqpreci(i) &
                      *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)  
             
             ENDIF !iflag_evap_prec==4

           ENDIF                      
         ENDDO
!!
       ELSE  ! iflag_bergeron
!>CR&JYG
!!
       DO i = 1, klon
         IF (rneb(i,k).GT.0.0) THEN
!CR on prend en compte la phase glace
!JLD inutile car on ne passe jamais ici si .not.ice_thermo
!           if (.not.ice_thermo) then
!           d_ql(i,k) = zoliq(i)
!           d_qi(i,k) = 0.
!           else
           d_ql(i,k) = (1-zfice(i))*zoliq(i)
           d_qi(i,k) = zfice(i)*zoliq(i)
!           endif
!<LTP
             IF (iflag_evap_prec == 4) THEN
                zrflcld(i) = zrflcld(i)+ MAX(zcond(i)*(1.-zfice(i))-zoliqp(i),0.0) &
                       *(paprs(i,k)-paprs(i,k+1))/(RG*dtime) 
                ziflcld(i) = ziflcld(i)+ MAX(zcond(i)*zfice(i)-zoliqi(i),0.0) &
                        *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)  
                znebprecipcld(i) = rneb(i,k)
                zrfl(i) = zrflcld(i) + zrflclr(i)
                zifl(i) = ziflcld(i) + ziflclr(i)        
!>LTP
             ELSE
!AJ<
                   zrfl(i) = zrfl(i)+ MAX(zcond(i)*(1.-zfice(i))-zoliqp(i),0.0) &
                       *(paprs(i,k)-paprs(i,k+1))/(RG*dtime) 
                        zifl(i) = zifl(i)+ MAX(zcond(i)*zfice(i)-zoliqi(i),0.0) &
                        *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)  
     !      zrfl(i) = zrfl(i)+  zpluie                         &
     !          *(paprs(i,k)-paprs(i,k+1))/(RG*dtime) 
     !      zifl(i) = zifl(i)+  zice                    &
     !               *(paprs(i,k)-paprs(i,k+1))/(RG*dtime)                                    
             ENDIF !iflag_evap_prec == 4             

!CR : on prend en compte l'effet Bergeron dans les flux de precipitation
           IF ((iflag_bergeron .EQ. 1) .AND. (zt(i) .LT. 273.15)) THEN
!<LTP
                IF (iflag_evap_prec == 4) THEN
                     zsolid = zrfl(i)
                     ziflclr(i) = ziflclr(i) +zrflclr(i)
                     ziflcld(i) = ziflcld(i) +zrflcld(i)
                     zifl(i) = ziflclr(i)+ziflcld(i)
                     zrflcld(i)=0.
                     zrflclr(i)=0.   
                     zrfl(i) = zrflclr(i)+zrflcld(i)
!>LTP
                ELSE
                     zsolid = zrfl(i)
                     zifl(i) = zifl(i)+zrfl(i)
                     zrfl(i) = 0.
                 ENDIF!iflag_evap_prec==4

           if (fl_cor_ebil .GT. 0) then
              zt(i)=zt(i)+zsolid*(RG*dtime)/(paprs(i,k)-paprs(i,k+1)) &
                      *(RLSTT-RLVTT)/RCPD/(1.0+RVTMP2*(zq(i)+zmqc(i)))
           else
              zt(i)=zt(i)+zsolid*(RG*dtime)/(paprs(i,k)-paprs(i,k+1)) &
                      *(RLSTT-RLVTT)/RCPD/(1.0+RVTMP2*zq(i))
           end if
           ENDIF  ! (iflag_bergeron .EQ. 1) .AND. (zt(i) .LT. 273.15)
!RC   

         ENDIF ! rneb(i,k).GT.0.0
       ENDDO

       ENDIF  ! iflag_bergeron .EQ. 2
     ENDIF  ! .NOT. ice_thermo

!CR: la fonte est faite au debut
!      IF (ice_thermo) THEN
!       DO i = 1, klon
!           zmelt = ((zt(i)-273.15)/(ztfondue-273.15))**2
!           zmelt = MIN(MAX(zmelt,0.),1.)
!           zrfl(i)=zrfl(i)+zmelt*zifl(i)
!           zifl(i)=zifl(i)*(1.-zmelt)
!           print*,zt(i),'octavio1'
!           zt(i)=zt(i)-zifl(i)*zmelt*(RG*dtime)/(paprs(i,k)-paprs(i,k+1)) &
!                      *RLMLT/RCPD/(1.0+RVTMP2*zq(i))
!           print*,zt(i),zrfl(i),zifl(i),zmelt,'octavio2'
!       ENDDO
!     ENDIF


!<LTP

!Limitation de la fraction surfacique couverte par les précipitations lorsque l'intensité locale du flux de précipitation descend en
!dessous de rain_int_min
    IF (iflag_evap_prec==4) THEN
        DO i=1, klon
            IF (zrflclr(i) + ziflclr(i) .GT. 0 ) THEN
               znebprecipclr(i) = min(znebprecipclr(i), max(zrflclr(i) &
                    /(znebprecipclr(i)*rain_int_min), &
                    ziflclr(i)/(znebprecipclr(i)*rain_int_min)))
            ELSE
                znebprecipclr(i)=0.
            ENDIF

            IF (zrflcld(i) + ziflcld(i) .GT. 0 ) THEN
               znebprecipcld(i) = min(znebprecipcld(i), &
                    max(zrflcld(i)/(znebprecipcld(i)*rain_int_min), &
                    ziflcld(i)/(znebprecipcld(i)*rain_int_min)))
            ELSE
                znebprecipcld(i)=0.
            ENDIF
       ENDDO
    ENDIf

!>LTP




       
     IF (.NOT. ice_thermo) THEN
       DO i = 1, klon
         IF (zt(i).LT.RTT) THEN
           psfl(i,k)=zrfl(i)
         ELSE
           prfl(i,k)=zrfl(i)
         ENDIF
       ENDDO
     ELSE
     ! JAM*************************************************
     ! Revoir partie ci-dessous: a quoi servent psfl et prfl?
     ! *****************************************************

       DO i = 1, klon
     !   IF (zt(i).LT.RTT) THEN
           psfl(i,k)=zifl(i) 
     !   ELSE
           prfl(i,k)=zrfl(i)
     !   ENDIF
!>AJ
       ENDDO
     ENDIF
     ! ----------------------------------------------------------------
     ! Fin de formation des precipitations
     ! ----------------------------------------------------------------
     !
     ! Calculer les tendances de q et de t:
     !
     DO i = 1, klon
        d_q(i,k) = zq(i) - q(i,k)
        d_t(i,k) = zt(i) - t(i,k)
     ENDDO
     !
     !AA--------------- Calcul du lessivage stratiforme  -------------

     DO i = 1,klon
        !
        if(zcond(i).gt.zoliq(i)+1.e-10) then
         beta(i,k) = (zcond(i)-zoliq(i))/zcond(i)/dtime
        else
         beta(i,k) = 0.
        endif
        zprec_cond(i) = MAX(zcond(i)-zoliq(i),0.0) &
             * (paprs(i,k)-paprs(i,k+1))/RG
        IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
           !AA lessivage nucleation LMD5 dans la couche elle-meme
          IF (iflag_t_glace.EQ.0) THEN
           if (t(i,k) .GE. t_glace_min_old) THEN
              zalpha_tr = a_tr_sca(3)
           else
              zalpha_tr = a_tr_sca(4)
           endif
          ELSE ! of IF (iflag_t_glace.EQ.0)
           if (t(i,k) .GE. t_glace_min) THEN
              zalpha_tr = a_tr_sca(3)
           else
              zalpha_tr = a_tr_sca(4)
           endif
          ENDIF
           zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
           pfrac_nucl(i,k)=pfrac_nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
           frac_nucl(i,k)= 1.-zneb(i)*zfrac_lessi 
           !
           ! nucleation avec un facteur -1 au lieu de -0.5
           zfrac_lessi = 1. - EXP(-zprec_cond(i)/zneb(i))
           pfrac_1nucl(i,k)=pfrac_1nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
        ENDIF
        !
     ENDDO      ! boucle sur i
     !
     !AA Lessivage par impaction dans les couches en-dessous
     DO kk = k-1, 1, -1
        DO i = 1, klon
           IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
             IF (iflag_t_glace.EQ.0) THEN
              if (t(i,kk) .GE. t_glace_min_old) THEN
                 zalpha_tr = a_tr_sca(1)
              else
                 zalpha_tr = a_tr_sca(2)
              endif
             ELSE ! of IF (iflag_t_glace.EQ.0)
              if (t(i,kk) .GE. t_glace_min) THEN
                 zalpha_tr = a_tr_sca(1)
              else
                 zalpha_tr = a_tr_sca(2)
              endif
             ENDIF
              zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
              pfrac_impa(i,kk)=pfrac_impa(i,kk)*(1.-zneb(i)*zfrac_lessi)
              frac_impa(i,kk)= 1.-zneb(i)*zfrac_lessi
           ENDIF
        ENDDO
     ENDDO
     !
     !AA===============================================================
     !                     FIN DE LA BOUCLE VERTICALE  
  end DO
  !
  !AA==================================================================
  !
  ! Pluie ou neige au sol selon la temperature de la 1ere couche
  !
!CR: si la thermo de la glace est active, on calcule zifl directement
  IF (.NOT.ice_thermo) THEN
  DO i = 1, klon
     IF ((t(i,1)+d_t(i,1)) .LT. RTT) THEN
!AJ<
!        snow(i) = zrfl(i)
        snow(i) = zrfl(i)+zifl(i)
!>AJ
        zlh_solid(i) = RLSTT-RLVTT
     ELSE
        rain(i) = zrfl(i)
        zlh_solid(i) = 0.
     ENDIF
  ENDDO

  ELSE
     DO i = 1, klon
        snow(i) = zifl(i)
        rain(i) = zrfl(i)
     ENDDO
   
   ENDIF
  !
  ! For energy conservation : when snow is present, the solification
  ! latent heat is considered.
!CR: si thermo de la glace, neige deja prise en compte
  IF (.not.ice_thermo) THEN
  DO k = 1, klev
     DO i = 1, klon
        zcpair=RCPD*(1.0+RVTMP2*(q(i,k)+d_q(i,k)))
        zmair(i)=(paprs(i,k)-paprs(i,k+1))/RG
        zm_solid = (prfl(i,k)-prfl(i,k+1)+psfl(i,k)-psfl(i,k+1))*dtime
        d_t(i,k) = d_t(i,k) + zlh_solid(i) *zm_solid / (zcpair*zmair(i))
     END DO
  END DO
  ENDIF
  !

  if (ncoreczq>0) then
     WRITE(lunout,*)'WARNING : ZQ dans fisrtilp ',ncoreczq,' val < 1.e-15.'
  endif

RETURN
END SUBROUTINE fisrtilp
END MODULE lmdz_lscp_old
