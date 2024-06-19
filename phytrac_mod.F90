!$Id: phytrac_mod.F90 4601 2023-06-30 22:07:30Z dcugnet $
MODULE phytrac_mod
!=================================================================================
! Interface between the LMDZ physical package and tracer computation.
! Chemistry modules (INCA, Reprobus or the more specific traclmdz routine)
! are called from phytrac.
!
!======================================================================
! Auteur(s) FH
! Objet: Moniteur general des tendances traceurs
!
! iflag_vdf_trac : Options for activating transport by vertical diffusion :
!     1. notmal
!     0. emission is injected in the first layer only, without diffusion
!    -1  no emission & no diffusion
! Modification 2013/07/22 : transformed into a module to pass tendencies to
!     physics outputs. Additional keys for controling activation of sub processes.
! Modification R. Pilon 10 octobre 2012 large scale scavenging incloud_scav + bc_scav
! Modification R. Pilon 01 janvier 2012 transport+scavenging KE scheme : cvltr
!=================================================================================

!
! Tracer tendencies, for outputs
!-------------------------------
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cl  ! Td couche limite/traceur
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_dec                            !RomP
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_cv  ! Td convection/traceur
! RomP >>>
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_insc
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_bcscav
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_evapls
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_ls
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_trsp
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sscav
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_sat
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_uscav
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qPr,qDi ! concentration tra dans pluie,air descente insaturee
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qPa,qMel
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: qTrdi,dtrcvMA ! conc traceur descente air insaturee et td convective MA
! RomP <<<
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_th  ! Td thermique
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_lessi_impa ! Td du lessivage par impaction
  REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: d_tr_lessi_nucl ! Td du lessivage par nucleation
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: qPrls      !jyg: concentration tra dans pluie LS a la surf.
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: d_tr_dry ! Td depot sec/traceur (1st layer),ALLOCATABLE,SAVE  jyg
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: flux_tr_dry ! depot sec/traceur (surface),ALLOCATABLE,SAVE    jyg

!$OMP THREADPRIVATE(qPa,qMel,qTrdi,dtrcvMA,d_tr_th,d_tr_lessi_impa,d_tr_lessi_nucl)
!$OMP THREADPRIVATE(d_tr_trsp,d_tr_sscav,d_tr_sat,d_tr_uscav,qPr,qDi)
!$OMP THREADPRIVATE(d_tr_insc,d_tr_bcscav,d_tr_evapls,d_tr_ls,qPrls)
!$OMP THREADPRIVATE(d_tr_cl,d_tr_dry,flux_tr_dry,d_tr_dec,d_tr_cv)

CONTAINS

  SUBROUTINE phytrac_init()

    USE dimphy
    USE infotrac_phy, ONLY: nbtr, type_trac
    USE tracco2i_mod, ONLY: tracco2i_init
    
    IMPLICIT NONE

    ALLOCATE(d_tr_cl(klon,klev,nbtr),d_tr_dry(klon,nbtr))
    ALLOCATE(flux_tr_dry(klon,nbtr),d_tr_dec(klon,klev,nbtr),d_tr_cv(klon,klev,nbtr))
    ALLOCATE(d_tr_insc(klon,klev,nbtr),d_tr_bcscav(klon,klev,nbtr))
    ALLOCATE(d_tr_evapls(klon,klev,nbtr),d_tr_ls(klon,klev,nbtr))
    ALLOCATE(qPrls(klon,nbtr),d_tr_trsp(klon,klev,nbtr))
    ALLOCATE(d_tr_sscav(klon,klev,nbtr),d_tr_sat(klon,klev,nbtr))
    ALLOCATE(d_tr_uscav(klon,klev,nbtr),qPr(klon,klev,nbtr),qDi(klon,klev,nbtr))
    ALLOCATE(qPa(klon,klev,nbtr),qMel(klon,klev,nbtr))
    ALLOCATE(qTrdi(klon,klev,nbtr),dtrcvMA(klon,klev,nbtr))
    ALLOCATE(d_tr_th(klon,klev,nbtr))
    ALLOCATE(d_tr_lessi_impa(klon,klev,nbtr),d_tr_lessi_nucl(klon,klev,nbtr))

    !===============================================================================
    !    -- Do specific treatment according to chemestry model or local LMDZ tracers
    !      
    !===============================================================================
    !   -- CO2 interactif --
    IF (ANY(type_trac == ['co2i','inco'])) CALL tracco2i_init()

       !   -- type_trac == 'co2i' ! PC 
       !   -- CO2 interactif --
       !   -- source is updated with FF and BB emissions 
       !   -- and net fluxes from ocean and orchidee 
       !   -- sign convention : positive into the atmosphere

  END SUBROUTINE phytrac_init

  SUBROUTINE phytrac(                                 &
       nstep,     julien,   gmtime,   debutphy,       &
       lafin,     pdtphys,  u, v,     t_seri,         &
       paprs,     pplay,    pmfu,     pmfd,           &
       pen_u,     pde_u,    pen_d,    pde_d,          &
       cdragh,    coefh,    fm_therm, entr_therm,     &
       yu1,       yv1,      ftsol,    pctsrf,         &
       ustar,     u10m,      v10m,                    &
       wstar,     ale_bl,      ale_wake,              &
       xlat,      xlon,                               &
       frac_impa,frac_nucl,beta_fisrt,beta_v1,        &
       presnivs,  pphis,    pphi,     albsol,         &
       sh,        ch, rh,   cldfra,   rneb,           &
       diafra,    cldliq,   itop_con, ibas_con,       &
       pmflxr,    pmflxs,   prfl,     psfl,           &
       da,        phi,      mp,       upwd,           &
       phi2,      d1a,      dam,      sij, wght_cvfd, &   ! RomP +RL
       wdtrainA,  wdtrainM, sigd,     clw, elij,      &   ! RomP 
       evap,      ep,       epmlmMm,  eplaMm,         &   ! RomP
       dnwd,      aerosol_couple,     flxmass_w,      &
       tau_aero,  piz_aero,  cg_aero, ccm,            &
       rfname,                                        &
       d_tr_dyn,                                      &   ! RomP
       tr_seri, init_source)         
    !
    !======================================================================
    ! Auteur(s) FH
    ! Objet: Moniteur general des tendances traceurs
    ! Modification R. Pilon 01 janvier 2012 transport+scavenging KE scheme : cvltr
    ! Modification R. Pilon 10 octobre 2012 large scale scavenging incloud_scav + bc_scav
    !======================================================================

    USE ioipsl
    USE phys_cal_mod, only : hour
    USE dimphy
    USE infotrac_phy, ONLY: nbtr, nqCO2, type_trac, conv_flg, pbl_flg
    USE strings_mod,  ONLY: int2str
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
    USE iophy
    USE traclmdz_mod
    USE tracinca_mod
    USE tracreprobus_mod
    USE indice_sol_mod
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    USE print_control_mod, ONLY: lunout
    USE aero_mod, ONLY : naero_grp
    USE lmdz_thermcell_dq, ONLY : thermcell_dq

    USE tracco2i_mod

#ifdef CPP_StratAer
    USE traccoag_mod
    USE phys_local_var_mod, ONLY: mdw
    USE phys_local_var_mod, ONLY: budg_dep_dry_ocs,   budg_dep_wet_ocs
    USE phys_local_var_mod, ONLY: budg_dep_dry_so2,   budg_dep_wet_so2
    USE phys_local_var_mod, ONLY: budg_dep_dry_h2so4, budg_dep_wet_h2so4
    USE phys_local_var_mod, ONLY: budg_dep_dry_part,  budg_dep_wet_part
    USE infotrac_phy, ONLY: nbtr_sulgas, id_OCS_strat, id_SO2_strat, id_H2SO4_strat
    USE strataer_nuc_mod, ONLY : tracstrataer_init
    USE aerophys
#endif

    IMPLICIT NONE

    INCLUDE "YOMCST.h"
    INCLUDE "clesphys.h"
    !==========================================================================
    !                   -- ARGUMENT DESCRIPTION --
    !==========================================================================

    ! Input arguments
    !----------------
    !Configuration grille,temps:
    INTEGER,INTENT(IN) :: nstep      ! Appel physique
    INTEGER,INTENT(IN) :: julien     ! Jour julien
    REAL,INTENT(IN)    :: gmtime     ! Heure courante
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    LOGICAL,INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique
    LOGICAL,INTENT(IN) :: lafin      ! le flag de la fin de la physique

    REAL,DIMENSION(klon),INTENT(IN) :: xlat    ! latitudes pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN) :: xlon    ! longitudes pour chaque point 
    !
    !Physique: 
    !--------
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: u       ! variable not used 
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: v       ! variable not used 
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rh      ! humidite relative
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: ch      ! eau liquide (+ glace si le traceur existe) 
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis
    REAL,DIMENSION(klev),INTENT(IN)        :: presnivs 
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldliq  ! eau condensee totale
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldfra  ! fraction nuageuse (tous les nuages)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: diafra  ! fraction nuageuse (convection ou stratus artificiels)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rneb    ! fraction nuageuse (grande echelle)
    !
    REAL                                   :: ql_incl ! contenu en eau liquide nuageuse dans le nuage ! ql_incl=oliq/rneb
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: beta_fisrt ! taux de conversion de l'eau cond (de fisrtilp)
    REAL,DIMENSION(klon,klev),INTENT(out)  :: beta_v1    ! -- (originale version)

    !
    INTEGER,DIMENSION(klon),INTENT(IN)     :: itop_con
    INTEGER,DIMENSION(klon),INTENT(IN)     :: ibas_con
    REAL,DIMENSION(klon),INTENT(IN)        :: albsol  ! albedo surface
    !
    !Dynamique
    !--------
    REAL,DIMENSION(klon,klev,nbtr),INTENT(IN)    :: d_tr_dyn
    !
    !Convection:
    !----------
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfu  ! flux de masse dans le panache montant
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfd  ! flux de masse dans le panache descendant
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_u ! flux entraine dans le panache montant
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_u ! flux detraine dans le panache montant
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_d ! flux entraine dans le panache descendant
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_d ! flux detraine dans le panache descendant

    !...Tiedke     
    REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: pmflxr, pmflxs ! Flux precipitant de pluie, neige aux interfaces [convection]
    REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: prfl, psfl ! Flux precipitant de pluie, neige aux interfaces [large-scale]

    LOGICAL,INTENT(IN)                       :: aerosol_couple
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: flxmass_w
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: tau_aero
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: piz_aero
    REAL,DIMENSION(klon,klev,naero_grp,2),INTENT(IN) :: cg_aero
    CHARACTER(len=4),DIMENSION(naero_grp),INTENT(IN) :: rfname 
    REAL,DIMENSION(klon,klev,2),INTENT(IN)   :: ccm 
    !... K.Emanuel
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: da
    REAL,DIMENSION(klon,klev,klev),INTENT(IN):: phi
    ! RomP >>>
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: d1a,dam
    REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi2
    !
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainA
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainM
    REAL,DIMENSION(klon),INTENT(IN)           :: sigd
    ! ---- RomP flux entraine, detraine et precipitant kerry Emanuel
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: evap
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: ep
    REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: sij
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: wght_cvfd          !RL
    REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: elij
    REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: epmlmMm
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: eplaMm
    REAL,DIMENSION(klon,klev),INTENT(IN)      :: clw
    ! RomP <<<
    !
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: mp
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: upwd      ! saturated updraft mass flux
    REAL,DIMENSION(klon,klev),INTENT(IN)     :: dnwd      ! saturated downdraft mass flux
    !
    !Thermiques:
    !----------
    REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: fm_therm
    REAL,DIMENSION(klon,klev),INTENT(INOUT)     :: entr_therm
    !
    !Couche limite:
    !--------------
    !
    REAL,DIMENSION(:),INTENT(IN)   :: cdragh          ! (klon) coeff drag pour T et Q
    REAL,DIMENSION(:,:),INTENT(IN) :: coefh           ! (klon,klev) coeff melange CL (m**2/s)
    REAL,DIMENSION(:),INTENT(IN)   :: ustar,u10m,v10m ! (klon) u* & vent a 10m (m/s)
    REAL,DIMENSION(:),INTENT(IN)   :: wstar,ale_bl,ale_wake ! (klon) w* and Avail. Lifting Ener.
    REAL,DIMENSION(:),INTENT(IN)   :: yu1             ! (klon) vents au premier niveau
    REAL,DIMENSION(:),INTENT(IN)   :: yv1             ! (klon) vents au premier niveau
    !
    !Lessivage:
    !----------
    !
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ccntrAA
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ccntrENV
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: coefcoli
    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: flag_cvltr
!$OMP THREADPRIVATE(ccntrAA,ccntrENV,coefcoli,flag_cvltr)
    REAL, DIMENSION(klon,klev) :: ccntrAA_3d
    REAL, DIMENSION(klon,klev) :: ccntrENV_3d
    REAL, DIMENSION(klon,klev) :: coefcoli_3d
    !
    ! pour le ON-LINE
    !
    REAL,DIMENSION(klon,klev),INTENT(IN) :: frac_impa ! fraction d'aerosols non impactes
    REAL,DIMENSION(klon,klev),INTENT(IN) :: frac_nucl ! fraction d'aerosols non nuclees

    ! Arguments necessaires pour les sources et puits de traceur:
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol  ! Temperature du sol (surf)(Kelvin)
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf ! Pourcentage de sol (nature du sol)

#ifdef CPP_StratAer
    REAL,DIMENSION(klon)           :: v_dep_dry !dry deposition velocity of stratospheric sulfate in m/s
#endif
    ! Output argument
    !----------------
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT) :: tr_seri ! Concentration Traceur [U/KgA]
    REAL,DIMENSION(klon,klev)                    :: sourceBE 
    REAL,DIMENSION(klon,nbtr), INTENT(IN) :: init_source

    !=======================================================================================
    !                        -- LOCAL VARIABLES --
    !=======================================================================================

    INTEGER :: i, k, it
    INTEGER :: nsplit

    !Sources et Reservoirs de traceurs (ex:Radon):
    !--------------------------------------------
    !
    REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: source  ! a voir lorsque le flux de surface est prescrit
!$OMP THREADPRIVATE(source)

    !
    !Entrees/Sorties:
    !---------------
    INTEGER                   :: iiq, ierr
    INTEGER                   :: nhori, nvert
    REAL                      :: zsto, zout, zjulian
    INTEGER,SAVE              :: nid_tra     ! pointe vers le fichier histrac.nc         
!$OMP THREADPRIVATE(nid_tra)
    REAL,DIMENSION(klon)      :: zx_tmp_fi2d ! variable temporaire grille physique
    INTEGER                   :: itau_w      ! pas de temps ecriture = nstep + itau_phy
    LOGICAL,PARAMETER         :: ok_sync=.TRUE.
    !
    ! Nature du traceur
    !------------------
    LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: aerosol  ! aerosol(it) = true  => aerosol => lessivage
!$OMP THREADPRIVATE(aerosol)                        ! aerosol(it) = false => gaz
    REAL,DIMENSION(klon,klev)             :: delp     ! epaisseur de couche (Pa)
    !
    ! Tendances de traceurs (Td) et flux de traceurs:
    !------------------------
    REAL,DIMENSION(klon,klev)      :: d_tr     ! Td dans l'atmosphere
    REAL,DIMENSION(klon,klev)      :: Mint
    REAL,DIMENSION(klon,klev,nbtr) :: zmfd1a
    REAL,DIMENSION(klon,klev,nbtr) :: zmfdam
    REAL,DIMENSION(klon,klev,nbtr) :: zmfphi2
    ! Physique
    !---------- 
    REAL,DIMENSION(klon,klev,nbtr) :: flestottr ! flux de lessivage dans chaque couche 
    REAL,DIMENSION(klon,klev)      :: zmasse    ! densite atmospherique Kg/m2
    REAL,DIMENSION(klon,klev)      :: ztra_th
    !PhH
    REAL,DIMENSION(klon,klev)      :: zrho
    REAL,DIMENSION(klon,klev)      :: zdz
    REAL                           :: evaplsc,dx,beta ! variable pour lessivage Genthon
    REAL,DIMENSION(klon)           :: his_dh          ! ---
    ! in-cloud scav variables
    REAL           :: ql_incloud_ref     ! ref value of in-cloud condensed water content

    !Controles:
    !---------
    INTEGER,SAVE :: iflag_vdf_trac,iflag_con_trac,iflag_the_trac
    INTEGER,SAVE  :: iflag_con_trac_omp, iflag_vdf_trac_omp,iflag_the_trac_omp
!$OMP THREADPRIVATE(iflag_vdf_trac,iflag_con_trac,iflag_the_trac)

    LOGICAL,SAVE :: lessivage
!$OMP THREADPRIVATE(lessivage)

    !RomP >>>
    INTEGER,SAVE  :: iflag_lscav_omp,iflag_lscav
    REAL, SAVE ::   ccntrAA_in,ccntrAA_omp
    REAL, SAVE ::   ccntrENV_in,ccntrENV_omp
    REAL, SAVE ::   coefcoli_in,coefcoli_omp

    LOGICAL,SAVE  :: convscav_omp,convscav
!$OMP THREADPRIVATE(iflag_lscav)
!$OMP THREADPRIVATE(ccntrAA_in,ccntrENV_in,coefcoli_in)
!$OMP THREADPRIVATE(convscav)
    !RomP <<<
    !######################################################################
    !                    -- INITIALIZATION --
    !######################################################################

    DO k=1,klev
       DO i=1,klon
          sourceBE(i,k)=0.
          Mint(i,k)=0.
          zrho(i,k)=0.
          zdz(i,k)=0.
       ENDDO
    ENDDO

    DO it=1, nbtr
       DO k=1,klev
          DO i=1,klon
             d_tr_insc(i,k,it)=0.
             d_tr_bcscav(i,k,it)=0.
             d_tr_evapls(i,k,it)=0.
             d_tr_ls(i,k,it)=0.
             d_tr_cv(i,k,it)=0.
             d_tr_cl(i,k,it)=0.
             d_tr_trsp(i,k,it)=0.
             d_tr_sscav(i,k,it)=0.
             d_tr_sat(i,k,it)=0.
             d_tr_uscav(i,k,it)=0.
             d_tr_lessi_impa(i,k,it)=0.
             d_tr_lessi_nucl(i,k,it)=0.
             qDi(i,k,it)=0.
             qPr(i,k,it)=0.
             qPa(i,k,it)=0.
             qMel(i,k,it)=0.
             qTrdi(i,k,it)=0.
             dtrcvMA(i,k,it)=0.
             zmfd1a(i,k,it)=0.
             zmfdam(i,k,it)=0.
             zmfphi2(i,k,it)=0.
          ENDDO
       ENDDO
    ENDDO

    DO it=1, nbtr
       DO i=1,klon
          d_tr_dry(i,it)=0.
          flux_tr_dry(i,it)=0.
       ENDDO
    ENDDO

    DO k = 1, klev
       DO i = 1, klon
          delp(i,k) = paprs(i,k)-paprs(i,k+1)
       ENDDO
    ENDDO

    IF (debutphy) THEN
       !!jyg
!$OMP BARRIER
       ecrit_tra=86400. ! frequence de stokage en dur
       ! obsolete car remplace par des ecritures dans phys_output_write
       !RomP >>>
       !
       !Config Key  = convscav
       !Config Desc = Convective scavenging switch: 0=off, 1=on.
       !Config Def  = .FALSE.
       !Config Help = 
       !
!$OMP MASTER
       convscav_omp=.FALSE.
       call getin('convscav', convscav_omp)
       iflag_vdf_trac_omp=1
       call getin('iflag_vdf_trac', iflag_vdf_trac_omp)
       iflag_con_trac_omp=1
       call getin('iflag_con_trac', iflag_con_trac_omp)
       iflag_the_trac_omp=1
       call getin('iflag_the_trac', iflag_the_trac_omp)
!$OMP END MASTER
!$OMP BARRIER
       convscav=convscav_omp
       iflag_vdf_trac=iflag_vdf_trac_omp
       iflag_con_trac=iflag_con_trac_omp
       iflag_the_trac=iflag_the_trac_omp
       write(lunout,*) 'phytrac passage dans routine conv avec lessivage', convscav
       !
       !Config Key  = iflag_lscav
       !Config Desc = Large scale scavenging parametrization: 0=none, 1=old(Genthon92),
       !              2=1+PHeinrich, 3=Reddy_Boucher2004, 4=3+RPilon.
       !Config Def  = 1
       !Config Help = 
       !
!$OMP MASTER
       iflag_lscav_omp=1
       call getin('iflag_lscav', iflag_lscav_omp)
       ccntrAA_omp=1
       ccntrENV_omp=1.
       coefcoli_omp=0.001
       call getin('ccntrAA', ccntrAA_omp)
       call getin('ccntrENV', ccntrENV_omp)
       call getin('coefcoli', coefcoli_omp)
!$OMP END MASTER
!$OMP BARRIER
       iflag_lscav=iflag_lscav_omp
       ccntrAA_in=ccntrAA_omp
       ccntrENV_in=ccntrENV_omp
       coefcoli_in=coefcoli_omp
       !
       SELECT CASE(iflag_lscav)
       CASE(0)
          WRITE(lunout,*)  'Large scale scavenging: none'
       CASE(1)
          WRITE(lunout,*)  'Large scale scavenging: C. Genthon, Tellus(1992), 44B, 371-389'
       CASE(2)
          WRITE(lunout,*)  'Large scale scavenging: C. Genthon, modified P. Heinrich'
       CASE(3)
          WRITE(lunout,*)  'Large scale scavenging: M. Shekkar Reddy and O. Boucher, JGR(2004), 109, D14202'
       CASE(4)
          WRITE(lunout,*)  'Large scale scavenging: Reddy and Boucher, modified R. Pilon'
       END SELECT
       !RomP <<<
       WRITE(*,*) 'FIRST TIME IN PHYTRAC : pdtphys(sec) = ',pdtphys,'ecrit_tra (sec) = ',ecrit_tra
       ALLOCATE( source(klon,nbtr), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('phytrac', 'pb in allocation 1',1)

       ALLOCATE( aerosol(nbtr), stat=ierr)
       IF (ierr /= 0) CALL abort_physic('phytrac', 'pb in allocation 2',1)


       ! Initialize module for specific tracers
       IF (type_trac == 'inca') THEN
          source(:,:)=init_source(:,:)
          CALL tracinca_init(aerosol,lessivage)
       ELSE IF (type_trac == 'repr') THEN
          source(:,:)=0.
       ELSE IF (type_trac == 'co2i') THEN
          source(:,:)=0.
          lessivage  = .FALSE.
          aerosol(:) = .FALSE.
          pbl_flg(:) = 1
          iflag_the_trac= 1
          iflag_vdf_trac= 1
          iflag_con_trac= 1
       ELSE IF (type_trac == 'inco') THEN
          source(:,1:nqCO2) = 0.                          ! from CO2i ModThL
          source(:,nqCO2+1:nbtr)=init_source(:,:)         ! from INCA ModThL
          aerosol(1:nqCO2) = .FALSE.                      ! from CO2i ModThL
          CALL tracinca_init(aerosol(nqCO2+1:nbtr),lessivage)     ! from INCA ModThL
          pbl_flg(1:nqCO2) = 1                            ! From CO2i ModThL
          iflag_the_trac = 1                              ! From CO2i
          iflag_vdf_trac = 1                              ! From CO2i
          iflag_con_trac = 1                              ! From CO2i
#ifdef CPP_StratAer
       ELSE IF (type_trac == 'coag') THEN
          source(:,:)=0.
          CALL tracstrataer_init(aerosol,lessivage)       ! init aerosols and lessivage param
#endif
       ELSE IF (type_trac == 'lmdz') THEN
          CALL traclmdz_init(pctsrf,xlat,xlon,ftsol,tr_seri,t_seri,pplay,sh,pdtphys,aerosol,lessivage)
       ENDIF

       !
       !--initialising coefficients for scavenging in the case of NP
       !
       ALLOCATE(flag_cvltr(nbtr))
       IF (iflag_con.EQ.3) THEN
          !
          ALLOCATE(ccntrAA(nbtr))
          ALLOCATE(ccntrENV(nbtr))
          ALLOCATE(coefcoli(nbtr))
          !
          DO it=1, nbtr
             IF (type_trac == 'repr') THEN
                 flag_cvltr(it)=.FALSE.
             ELSE IF (type_trac == 'inca') THEN
!                IF ((it.EQ.id_Rn222) .OR. ((it.GE.id_SO2) .AND. (it.LE.id_NH3)) ) THEN 
!                   !--gas-phase species
!                   flag_cvltr(it)=.FALSE.
!
!                ELSEIF ( (it.GE.id_CIDUSTM) .AND. (it.LE.id_AIN) ) THEN
!                   !--insoluble aerosol species
!                   flag_cvltr(it)=.TRUE.
!                   ccntrAA(it)=0.7
!                   ccntrENV(it)=0.7
!                   coefcoli(it)=0.001
!                ELSEIF ( (it.EQ.id_Pb210) .OR. ((it.GE.id_CSSSM) .AND. (it.LE.id_SSN))) THEN
!                   !--soluble aerosol species
!                   flag_cvltr(it)=.TRUE.
!                   ccntrAA(it)=0.9
!                   ccntrENV(it)=0.9 
!                   coefcoli(it)=0.001
!                ELSE 
!                   WRITE(lunout,*) 'pb it=', it
!                   CALL abort_physic('phytrac','pb it scavenging',1)
!                ENDIF
                !--test OB 
                !--for now we do not scavenge in cvltr
                flag_cvltr(it)=.FALSE.
             ELSE IF (type_trac == 'co2i') THEN
                !--co2 tracers are not scavenged
                flag_cvltr(it)=.FALSE.
             ELSE IF (type_trac == 'inco') THEN     ! Add ThL
                flag_cvltr(it)=.FALSE.
#ifdef CPP_StratAer
             ELSE IF (type_trac == 'coag') THEN
                IF (convscav.and.aerosol(it)) THEN
                   flag_cvltr(it)=.TRUE.
                   ccntrAA(it) =ccntrAA_in    
                   ccntrENV(it)=ccntrENV_in
                   coefcoli(it)=coefcoli_in
                ELSE
                   flag_cvltr(it)=.FALSE.
                ENDIF
#endif
             ELSE IF (type_trac == 'lmdz') THEN
                IF (convscav.and.aerosol(it)) THEN
                   flag_cvltr(it)=.TRUE.
                   ccntrAA(it) =ccntrAA_in    !--a modifier par JYG a lire depuis fichier 
                   ccntrENV(it)=ccntrENV_in
                   coefcoli(it)=coefcoli_in
                ELSE 
                   flag_cvltr(it)=.FALSE.
                ENDIF
             ENDIF
          ENDDO
          !
       ELSE ! iflag_con .ne. 3 
          flag_cvltr(:) = .FALSE. 
       ENDIF
       !
       ! print out all tracer flags 
       !
       WRITE(lunout,*) 'print out all tracer flags'
       WRITE(lunout,*) 'type_trac      =', type_trac
       WRITE(lunout,*) 'config_inca    =', config_inca
       WRITE(lunout,*) 'iflag_con_trac =', iflag_con_trac
       WRITE(lunout,*) 'iflag_con      =', iflag_con
       WRITE(lunout,*) 'convscav       =', convscav
       WRITE(lunout,*) 'iflag_lscav    =', iflag_lscav
       WRITE(lunout,*) 'aerosol        =', aerosol
       WRITE(lunout,*) 'iflag_the_trac =', iflag_the_trac
       WRITE(lunout,*) 'iflag_thermals =', iflag_thermals
       WRITE(lunout,*) 'iflag_vdf_trac =', iflag_vdf_trac
       WRITE(lunout,*) 'pbl_flg        =', pbl_flg
       WRITE(lunout,*) 'lessivage      =', lessivage
       write(lunout,*)  'flag_cvltr    = ', flag_cvltr

       IF (lessivage .AND. ANY(type_trac == ['inca','inco'])) &
          CALL abort_physic('phytrac', 'lessivage=T config_inca=inca impossible',1)
       !
    ENDIF ! of IF (debutphy)
    !############################################ END INITIALIZATION #######

    DO k=1,klev
       DO i=1,klon
          zmasse(i,k)=(paprs(i,k)-paprs(i,k+1))/rg
       ENDDO
    ENDDO
    !
    IF (id_be .GT. 0) THEN
       DO k=1,klev
          DO i=1,klon
             sourceBE(i,k)=srcbe(i,k)       !RomP  -> pour sortie histrac
          ENDDO
       ENDDO
    ENDIF

    !===============================================================================
    !    -- Do specific treatment according to chemestry model or local LMDZ tracers
    !      
    !===============================================================================
    IF (type_trac == 'inca') THEN
       !    -- CHIMIE INCA  config_inca = aero or chem --
       ! Appel fait en fin de phytrac pour avoir les emissions modifiees par 
       ! la couche limite et la convection avant le calcul de la chimie 

    ELSE IF (type_trac == 'repr') THEN
       !   -- CHIMIE REPROBUS --
       CALL tracreprobus(pdtphys, gmtime, debutphy, julien, &
            presnivs, xlat, xlon, pphis, pphi, &
            t_seri, pplay, paprs, sh , &
            tr_seri)

    ELSE IF (type_trac == 'co2i') THEN
       !   -- CO2 interactif --
       !   -- source is updated with FF and BB emissions 
       !   -- and net fluxes from ocean and orchidee 
       !   -- sign convention : positive into the atmosphere

       CALL tracco2i(pdtphys, debutphy, &
            xlat, xlon, pphis, pphi, &
            t_seri, pplay, paprs, tr_seri, source)
    ELSE IF (type_trac == 'inco') THEN      ! Add ThL
       CALL tracco2i(pdtphys, debutphy, &
            xlat, xlon, pphis, pphi, &
            t_seri, pplay, paprs, tr_seri, source)

#ifdef CPP_StratAer
    ELSE IF (type_trac == 'coag') THEN
       !   --STRATOSPHERIC AER IN THE STRAT --
       CALL traccoag(pdtphys, gmtime, debutphy, julien, &
            presnivs, xlat, xlon, pphis, pphi, &
            t_seri, pplay, paprs, sh, rh , &
            tr_seri)
#endif
    ELSE IF (type_trac == 'lmdz') THEN
       !    -- Traitement des traceurs avec traclmdz
       CALL traclmdz(nstep, julien, gmtime, pdtphys, t_seri, paprs, pplay, &
            cdragh,  coefh, yu1, yv1, ftsol, pctsrf, xlat, xlon,iflag_vdf_trac>=0,sh, &
            rh, pphi, ustar, wstar, ale_bl, ale_wake,  u10m, v10m, &
            tr_seri, source, d_tr_cl,d_tr_dec, zmasse)               !RomP
    ENDIF
    !======================================================================
    !       -- Calcul de l'effet de la convection --
    !======================================================================

    IF (iflag_con_trac==1) THEN

       DO it=1, nbtr
          IF ( conv_flg(it) == 0 ) CYCLE
          IF (iflag_con.LT.2) THEN
             !--pas de transport convectif
             d_tr_cv(:,:,it)=0.

          ELSE IF (iflag_con.EQ.2) THEN
             !--ancien transport convectif de Tiedtke

             CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
                  pplay, paprs, tr_seri(:,:,it), d_tr_cv(:,:,it))
          ELSE   
             !--nouveau transport convectif de Emanuel

             IF (flag_cvltr(it)) THEN
                !--nouveau transport convectif de Emanuel avec lessivage convectif
                !
                !
                ccntrAA_3d(:,:) =ccntrAA(it)
                ccntrENV_3d(:,:)=ccntrENV(it)
                coefcoli_3d(:,:)=coefcoli(it)

                !--beware this interface is a bit weird because it is called for each tracer 
                !--with the full array tr_seri even if only item it is processed

                CALL cvltr_scav(pdtphys, da, phi,phi2,d1a,dam, mp,ep,         &
                     sigd,sij,wght_cvfd,clw,elij,epmlmMm,eplaMm,              &     
                     pmflxr,pmflxs,evap,t_seri,wdtrainA,wdtrainM,             &   
                     paprs,it,tr_seri,upwd,dnwd,itop_con,ibas_con,            &
                     ccntrAA_3d,ccntrENV_3d,coefcoli_3d,                      &
                     d_tr_cv,d_tr_trsp,d_tr_sscav,d_tr_sat,d_tr_uscav,qDi,qPr,&
                     qPa,qMel,qTrdi,dtrcvMA,Mint,                             &
                     zmfd1a,zmfphi2,zmfdam)


             ELSE  !---flag_cvltr(it).EQ.FALSE
                !--nouveau transport convectif de Emanuel mais pas de lessivage convectif

                !--beware this interface is a bit weird because it is called for each tracer 
                !--with the full array tr_seri even if only item it is processed
                !
                CALL cvltr_noscav(it,pdtphys, da, phi,mp,wght_cvfd,paprs,pplay, &  !jyg
                     tr_seri,upwd,dnwd,d_tr_cv)                                !jyg

             ENDIF

          ENDIF !--iflag

          !--on ajoute les tendances

          DO k = 1, klev
             DO i = 1, klon        
                tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_cv(i,k,it)
             ENDDO
          ENDDO

          CALL minmaxqfi(tr_seri(:,:,it),0.,1.e33,'convection it = '//TRIM(int2str(it)))

       ENDDO ! nbtr

#ifdef CPP_StratAer
       IF (type_trac=='coag') THEN
         ! initialize wet deposition flux of sulfur
         budg_dep_wet_ocs(:)=0.0
         budg_dep_wet_so2(:)=0.0
         budg_dep_wet_h2so4(:)=0.0
         budg_dep_wet_part(:)=0.0
         ! compute wet deposition flux of sulfur (sum over gases and particles)
         ! and convert to kg(S)/m2/s
         DO i = 1, klon
         DO k = 1, klev
         DO it = 1, nbtr
         !do not include SO2 because most of it comes trom the troposphere
           IF (it==id_OCS_strat) THEN
             budg_dep_wet_ocs(i)=budg_dep_wet_ocs(i)+d_tr_cv(i,k,it)*(mSatom/mOCSmol) &
                            & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
           ELSEIF (it==id_SO2_strat) THEN
             budg_dep_wet_so2(i)=budg_dep_wet_so2(i)+d_tr_cv(i,k,it)*(mSatom/mSO2mol) &
                            & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
           ELSEIF (it==id_H2SO4_strat) THEN
             budg_dep_wet_h2so4(i)=budg_dep_wet_h2so4(i)+d_tr_cv(i,k,it)*(mSatom/mH2SO4mol) &
                            & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
           ELSEIF (it.GT.nbtr_sulgas) THEN
             budg_dep_wet_part(i)=budg_dep_wet_part(i)+d_tr_cv(i,k,it)*(mSatom/mH2SO4mol)  &
                            & *dens_aer_dry*4./3.*RPI*(mdw(it-nbtr_sulgas)/2.)**3 &
                            & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
           ENDIF
         ENDDO
         ENDDO
         ENDDO
       ENDIF
#endif

    ENDIF ! convection

    !======================================================================
    !    -- Calcul de l'effet des thermiques --
    !======================================================================

    DO it=1,nbtr
       DO k=1,klev
          DO i=1,klon
             d_tr_th(i,k,it)=0.
             tr_seri(i,k,it)=MAX(tr_seri(i,k,it),0.)
! the next safeguard causes some problem for stratospheric aerosol tracers (particle number)
! and there is little justification for it so it is commented out (4 December 2017) by OB
! if reinstated please keep the ifndef CPP_StratAer
!#ifndef CPP_StratAer
!             tr_seri(i,k,it)=MIN(tr_seri(i,k,it),1.e10)
!#endif
          ENDDO
       ENDDO
    ENDDO

    IF (iflag_thermals.GT.0.AND.iflag_the_trac>0) THEN   

       DO it=1, nbtr

          CALL thermcell_dq(klon,klev,1,pdtphys,fm_therm,entr_therm, &
               zmasse,tr_seri(1:klon,1:klev,it),        &
               d_tr_th(1:klon,1:klev,it),ztra_th,0 )

          DO k=1,klev
             DO i=1,klon
                d_tr_th(i,k,it)=pdtphys*d_tr_th(i,k,it)
                tr_seri(i,k,it)=MAX(tr_seri(i,k,it)+d_tr_th(i,k,it),0.)
             ENDDO
          ENDDO

       ENDDO ! it

    ENDIF ! Thermiques

    !======================================================================
    !     -- Calcul de l'effet de la couche limite --
    !======================================================================

    IF (iflag_vdf_trac==1) THEN

       !  Injection during BL mixing
       !
#ifdef CPP_StratAer
       IF (type_trac=='coag') THEN

         ! initialize dry deposition flux of sulfur
         budg_dep_dry_ocs(:)=0.0
         budg_dep_dry_so2(:)=0.0
         budg_dep_dry_h2so4(:)=0.0
         budg_dep_dry_part(:)=0.0

         ! compute dry deposition velocity as function of surface type (numbers
         ! from IPSL note 23, 2002)
         v_dep_dry(:) =  pctsrf(:,is_ter) * 2.5e-3 &
                     & + pctsrf(:,is_oce) * 0.5e-3 &
                     & + pctsrf(:,is_lic) * 2.5e-3 &
                     & + pctsrf(:,is_sic) * 2.5e-3

         ! compute surface dry deposition flux
         zrho(:,1)=pplay(:,1)/t_seri(:,1)/RD

         DO it=1, nbtr
          source(:,it) = - v_dep_dry(:) * tr_seri(:,1,it) * zrho(:,1)
         ENDDO

       ENDIF
#endif

       DO it=1, nbtr
          !
          IF (pbl_flg(it) /= 0) THEN
             !
             CALL cltrac(pdtphys, coefh,t_seri,       &
                  tr_seri(:,:,it), source(:,it),      &
                  paprs, pplay, delp,                 &
                  d_tr_cl(:,:,it),d_tr_dry(:,it),flux_tr_dry(:,it))
             !
             tr_seri(:,:,it)=tr_seri(:,:,it)+d_tr_cl(:,:,it)
             !
#ifdef CPP_StratAer
             IF (type_trac=='coag') THEN
               ! compute dry deposition flux of sulfur (sum over gases and particles)
               IF (it==id_OCS_strat) THEN
                 budg_dep_dry_ocs(:)=budg_dep_dry_ocs(:)-source(:,it)*(mSatom/mOCSmol)
               ELSEIF (it==id_SO2_strat) THEN
                 budg_dep_dry_so2(:)=budg_dep_dry_so2(:)-source(:,it)*(mSatom/mSO2mol)
               ELSEIF (it==id_H2SO4_strat) THEN
                 budg_dep_dry_h2so4(:)=budg_dep_dry_h2so4(:)-source(:,it)*(mSatom/mH2SO4mol)
               ELSEIF (it.GT.nbtr_sulgas) THEN
                 budg_dep_dry_part(:)=budg_dep_dry_part(:)-source(:,it)*(mSatom/mH2SO4mol)*dens_aer_dry &
                                & *4./3.*RPI*(mdw(it-nbtr_sulgas)/2.)**3
               ENDIF
             ENDIF
#endif
             !
          ENDIF
          !
       ENDDO
       !
    ELSE IF (iflag_vdf_trac==0) THEN
       !
       !   Injection of source in the first model layer
       !
       DO it=1,nbtr
          d_tr_cl(:,1,it)=source(:,it)*RG/delp(:,1)*pdtphys
          tr_seri(:,1,it)=tr_seri(:,1,it)+d_tr_cl(:,1,it)
       ENDDO
       d_tr_cl(:,2:klev,1:nbtr)=0.
       !
    ELSE IF (iflag_vdf_trac==-1) THEN
       !
       ! Nothing happens
       d_tr_cl=0.
       !
    ELSE
       !
       CALL abort_physic('iflag_vdf_trac', 'cas non prevu',1)
       !
    ENDIF ! couche limite

    !======================================================================
    !   Calcul de l'effet de la precipitation grande echelle
    !   POUR INCA le lessivage est fait directement dans INCA
    !======================================================================

    IF (lessivage) THEN

       ql_incloud_ref = 10.e-4
       ql_incloud_ref =  5.e-4


       ! calcul du contenu en eau liquide au sein du nuage
       ql_incl = ql_incloud_ref
       ! choix du lessivage
       !
       IF (iflag_lscav .EQ. 3 .OR. iflag_lscav .EQ. 4) THEN
          ! ********  Olivier Boucher version (3) possibly with modified ql_incl (4)
          !
          DO it = 1, nbtr

             IF (aerosol(it)) THEN
             !  incloud scavenging and removal by large scale rain ! orig : ql_incl was replaced by 0.5e-3 kg/kg
             ! the value 0.5e-3 kg/kg is from Giorgi and Chameides (1986), JGR
             ! Liu (2001) proposed to use 1.5e-3 kg/kg

             CALL lsc_scav(pdtphys,it,iflag_lscav,aerosol,ql_incl,prfl,psfl,rneb,beta_fisrt,  &
                           beta_v1,pplay,paprs,t_seri,tr_seri,d_tr_insc,d_tr_bcscav,d_tr_evapls,qPrls)

             !large scale scavenging tendency
             DO k = 1, klev
                DO i = 1, klon
                   d_tr_ls(i,k,it)=d_tr_insc(i,k,it)+d_tr_bcscav(i,k,it)+d_tr_evapls(i,k,it)
                   tr_seri(i,k,it)=tr_seri(i,k,it)+d_tr_ls(i,k,it)
                ENDDO
             ENDDO
             CALL minmaxqfi(tr_seri(:,:,it),0.,1.e33,'lsc scav it = '//TRIM(int2str(it)))
             ENDIF

          ENDDO  !tr

#ifdef CPP_StratAer
         IF (type_trac=='coag') THEN
           ! compute wet deposition flux of sulfur (sum over gases and
           ! particles) and convert to kg(S)/m2/s
           ! adding contribution of d_tr_ls to d_tr_cv (above)
           DO i = 1, klon
           DO k = 1, klev
           DO it = 1, nbtr
             IF (it==id_OCS_strat) THEN
               budg_dep_wet_ocs(i)=budg_dep_wet_ocs(i)+d_tr_ls(i,k,it)*(mSatom/mOCSmol) &
                              & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
             ELSEIF (it==id_SO2_strat) THEN
               budg_dep_wet_so2(i)=budg_dep_wet_so2(i)+d_tr_ls(i,k,it)*(mSatom/mSO2mol) &
                              & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
             ELSEIF (it==id_H2SO4_strat) THEN
               budg_dep_wet_h2so4(i)=budg_dep_wet_h2so4(i)+d_tr_ls(i,k,it)*(mSatom/mH2SO4mol) &
                              & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
             ELSEIF (it.GT.nbtr_sulgas) THEN
               budg_dep_wet_part(i)=budg_dep_wet_part(i)+d_tr_ls(i,k,it)*(mSatom/mH2SO4mol)  &
                              & *dens_aer_dry*4./3.*RPI*(mdw(it-nbtr_sulgas)/2.)**3 &
                              & *(paprs(i,k)-paprs(i,k+1))/RG/pdtphys
             ENDIF
           ENDDO
           ENDDO
           ENDDO
         ENDIF
#endif

       ELSE IF (iflag_lscav .EQ. 2) THEN ! frac_impa, frac_nucl
          ! *********   modified  old version

          d_tr_lessi_nucl(:,:,:) = 0. 
          d_tr_lessi_impa(:,:,:) = 0.
          flestottr(:,:,:) = 0. 
          ! Tendance des aerosols nuclees et impactes 
          DO it = 1, nbtr
             IF (aerosol(it)) THEN
                his_dh(:)=0.
                DO k = 1, klev
                   DO i = 1, klon
                      !PhH
                      zrho(i,k)=pplay(i,k)/t_seri(i,k)/RD
                      zdz(i,k)=(paprs(i,k)-paprs(i,k+1))/zrho(i,k)/RG
                      !
                   ENDDO
                ENDDO

                DO k=klev-1, 1, -1
                   DO i=1, klon
                      !             d_tr_ls(i,k,it)=tr_seri(i,k,it)*(frac_impa(i,k)*frac_nucl(i,k)-1.)
                      dx=d_tr_ls(i,k,it)
                      his_dh(i)=his_dh(i)-dx*zrho(i,k)*zdz(i,k)/pdtphys !  kg/m2/s
                      evaplsc = prfl(i,k) - prfl(i,k+1) + psfl(i,k) - psfl(i,k+1)
                      ! Evaporation Partielle -> Liberation Partielle 0.5*evap
                      IF ( evaplsc .LT.0..and.abs(prfl(i,k+1)+psfl(i,k+1)).gt.1.e-10) THEN
                         evaplsc = (-evaplsc)/(prfl(i,k+1)+psfl(i,k+1))
                         ! evaplsc est donc positif, his_dh(i) est positif
                         !--------------  
                         d_tr_evapls(i,k,it)=0.5*evaplsc*(d_tr_lessi_nucl(i,k+1,it) &
                              +d_tr_lessi_impa(i,k+1,it))
                         !-------------   d_tr_evapls(i,k,it)=-0.5*evaplsc*(d_tr_lsc(i,k+1,it))
                         beta=0.5*evaplsc
                         if ((prfl(i,k)+psfl(i,k)).lt.1.e-10) THEN
                            beta=1.0*evaplsc
                         endif
                         dx=beta*his_dh(i)/zrho(i,k)/zdz(i,k)*pdtphys
                         his_dh(i)=(1.-beta)*his_dh(i)   ! tracer from
                         d_tr_evapls(i,k,it)=dx
                      ENDIF
                      d_tr_ls(i,k,it)=tr_seri(i,k,it)*(frac_impa(i,k)*frac_nucl(i,k)-1.) &
                           +d_tr_evapls(i,k,it)

                      !--------------
                      d_tr_lessi_nucl(i,k,it) = d_tr_lessi_nucl(i,k,it) +    &
                           ( 1 - frac_nucl(i,k) )*tr_seri(i,k,it)
                      d_tr_lessi_impa(i,k,it) = d_tr_lessi_impa(i,k,it) +    &
                           ( 1 - frac_impa(i,k) )*tr_seri(i,k,it)
                      !
                      ! Flux lessivage total 
                      flestottr(i,k,it) = flestottr(i,k,it) -   &
                           ( d_tr_lessi_nucl(i,k,it)   +        &
                           d_tr_lessi_impa(i,k,it) ) *          &
                           ( paprs(i,k)-paprs(i,k+1) ) /        &
                           (RG * pdtphys)
                      !! Mise a jour des traceurs due a l'impaction,nucleation 
                      !                 tr_seri(i,k,it)=tr_seri(i,k,it)*frac_impa(i,k)*frac_nucl(i,k)
                      !!  calcul de la tendance liee au lessivage stratiforme
                      !                 d_tr_ls(i,k,it)=tr_seri(i,k,it)*&
                      !                                (1.-1./(frac_impa(i,k)*frac_nucl(i,k)))
                      !--------------
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          ! *********   end modified old version

       ELSE IF (iflag_lscav .EQ. 1) THEN ! frac_impa, frac_nucl
          ! *********    old version

          d_tr_lessi_nucl(:,:,:) = 0. 
          d_tr_lessi_impa(:,:,:) = 0.
          flestottr(:,:,:) = 0. 
          !=========================
          ! LESSIVAGE LARGE SCALE : 
          !=========================

          ! Tendance des aerosols nuclees et impactes 
          ! -----------------------------------------
          DO it = 1, nbtr
             IF (aerosol(it)) THEN
                DO k = 1, klev
                   DO i = 1, klon
                      d_tr_lessi_nucl(i,k,it) = d_tr_lessi_nucl(i,k,it) +    &
                           ( 1 - frac_nucl(i,k) )*tr_seri(i,k,it)
                      d_tr_lessi_impa(i,k,it) = d_tr_lessi_impa(i,k,it) +    &
                           ( 1 - frac_impa(i,k) )*tr_seri(i,k,it)

                      !
                      ! Flux lessivage total 
                      ! ------------------------------------------------------------
                      flestottr(i,k,it) = flestottr(i,k,it) -   &
                           ( d_tr_lessi_nucl(i,k,it)   +        &
                           d_tr_lessi_impa(i,k,it) ) *          &
                           ( paprs(i,k)-paprs(i,k+1) ) /        &
                           (RG * pdtphys)
                      !
                      ! Mise a jour des traceurs due a l'impaction,nucleation 
                      ! ----------------------------------------------------------------------
                      tr_seri(i,k,it)=tr_seri(i,k,it)*frac_impa(i,k)*frac_nucl(i,k)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

          ! *********   end old version
       ENDIF  !  iflag_lscav . EQ. 1, 2, 3 or 4
       !
    ENDIF !  lessivage


    !    -- CHIMIE INCA  config_inca = aero or chem --
    IF (ANY(type_trac == ['inca','inco'])) THEN  ! ModThL

       CALL tracinca(&
            nstep,    julien,   gmtime,         lafin,     &
            pdtphys,  t_seri,   paprs,          pplay,     &
            pmfu,     upwd,     ftsol,  pctsrf, pphis,     &
            pphi,     albsol,   sh,    ch,     rh,        &
            cldfra,   rneb,     diafra,         cldliq,    &
            itop_con, ibas_con, pmflxr,         pmflxs,    &
            prfl,     psfl,     aerosol_couple, flxmass_w, &
            tau_aero, piz_aero, cg_aero,        ccm,       &
            rfname,                                        &
            tr_seri(:,:,1+nqCO2:nbtr),  source(:,1+nqCO2:nbtr))  ! ModThL  
    ENDIF

  END SUBROUTINE phytrac

END MODULE
