!
! $Id: phys_state_var_mod.F90 4677 2023-09-07 11:07:27Z idelkadi $
!
      MODULE phys_state_var_mod
! Variables sauvegardees pour le startphy.nc
!======================================================================
!
!
!======================================================================
! Declaration des variables
      USE dimphy
      USE netcdf, only: nf90_fill_real
      INTEGER, PARAMETER :: nlevSTD=17
      INTEGER, PARAMETER :: nlevSTD8=8
      INTEGER, PARAMETER :: nlevSTD3=3
      INTEGER, PARAMETER :: nout=3
      INTEGER, PARAMETER :: napisccp=1
      INTEGER, SAVE :: radpas  ! radiation is called every "radpas" step
      INTEGER, SAVE :: cvpas   ! convection is called every "cvpas" step
      INTEGER, SAVE :: cvpas_0 = 1 ! reference value for cvpas
      INTEGER, SAVE :: wkpas   ! wake scheme is called every "wkpas" step
      REAL, PARAMETER :: missing_val_nf90=nf90_fill_real
!$OMP THREADPRIVATE(radpas)
!$OMP THREADPRIVATE(cvpas)
!$OMP THREADPRIVATE(cvpas_0)
!$OMP THREADPRIVATE(wkpas)
      REAL, SAVE :: phys_tstep=0, solaire_etat0
!$OMP THREADPRIVATE(phys_tstep, solaire_etat0)

      REAL, ALLOCATABLE, SAVE :: pctsrf(:,:)
!$OMP THREADPRIVATE(pctsrf)
      REAL, ALLOCATABLE, SAVE :: ftsol(:,:)
!$OMP THREADPRIVATE(ftsol)
      REAL, ALLOCATABLE, SAVE :: beta_aridity(:,:)
!$OMP THREADPRIVATE(beta_aridity)
      REAL,ALLOCATABLE,SAVE :: qsol(:),fevap(:,:),z0m(:,:),z0h(:,:),agesno(:,:)
!$OMP THREADPRIVATE(qsol,fevap,z0m,z0h,agesno)
!FC drag des arbres
      REAL, ALLOCATABLE, SAVE :: treedrg(:,:,:)
!$OMP THREADPRIVATE(treedrg)

!      character(len=6), SAVE :: ocean
!!!!!!$OMP THREADPRIVATE(ocean)
!      logical, SAVE :: ok_veget 
!!!!!!$OMP THREADPRIVATE(ok_veget)
      REAL, ALLOCATABLE, SAVE :: falb1(:,:), falb2(:,:)
!$OMP THREADPRIVATE(falb1, falb2)

!albedo SB >>>
      REAL, ALLOCATABLE, SAVE :: falb_dif(:,:,:), falb_dir(:,:,:)
      REAL, ALLOCATABLE, SAVE :: chl_con(:)
!$OMP THREADPRIVATE(falb_dir,falb_dif,chl_con)
!albedo SB <<<


      REAL, ALLOCATABLE, SAVE :: rain_fall(:), snow_fall(:), bs_fall(:)
!$OMP THREADPRIVATE( rain_fall, snow_fall, bs_fall)
      REAL, ALLOCATABLE, SAVE :: solsw(:), solswfdiff(:), sollw(:)
!$OMP THREADPRIVATE(solsw, solswfdiff, sollw)
      REAL, ALLOCATABLE, SAVE :: radsol(:)
!$OMP THREADPRIVATE(radsol)
      REAL, ALLOCATABLE, SAVE :: swradcorr(:)
!$OMP THREADPRIVATE(swradcorr)
#ifdef ISO
      REAL,ALLOCATABLE,SAVE :: xtsol(:,:),fxtevap(:,:,:)
!$OMP THREADPRIVATE(xtsol,fxtevap)
      REAL, ALLOCATABLE, SAVE :: xtrain_fall(:,:), xtsnow_fall(:,:)
!$OMP THREADPRIVATE(xtrain_fall,xtsnow_fall)
#endif

!clesphy0 param physiq
!
! Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
!
      REAL, ALLOCATABLE, SAVE :: zmea(:), zstd(:), zsig(:), zgam(:)
!$OMP THREADPRIVATE(zmea, zstd, zsig, zgam)
      REAL, ALLOCATABLE, SAVE :: zthe(:), zpic(:), zval(:)
!$OMP THREADPRIVATE(zthe, zpic, zval)
!     REAL tabcntr0(100)
      REAL, ALLOCATABLE, SAVE :: rugoro(:)
!$OMP THREADPRIVATE(rugoro)
      REAL, ALLOCATABLE, SAVE :: t_ancien(:,:), q_ancien(:,:)
!$OMP THREADPRIVATE(t_ancien, q_ancien)
      REAL, ALLOCATABLE, SAVE :: ql_ancien(:,:), qs_ancien(:,:), qbs_ancien(:,:)
!$OMP THREADPRIVATE(ql_ancien, qs_ancien, qbs_ancien)
      REAL, ALLOCATABLE, SAVE :: prw_ancien(:), prlw_ancien(:), prsw_ancien(:), prbsw_ancien(:)
!$OMP THREADPRIVATE(prw_ancien, prlw_ancien, prsw_ancien, prbsw_ancien)
#ifdef ISO
      REAL, ALLOCATABLE, SAVE :: xt_ancien(:,:,:),xtl_ancien(:,:,:),xts_ancien(:,:,:)
!$OMP THREADPRIVATE(xt_ancien,xtl_ancien,xts_ancien)
#endif
      REAL, ALLOCATABLE, SAVE :: u_ancien(:,:), v_ancien(:,:)
!$OMP THREADPRIVATE(u_ancien, v_ancien)
!!! RomP >>>
      REAL, ALLOCATABLE, SAVE :: tr_ancien(:,:,:)
!$OMP THREADPRIVATE(tr_ancien)
!!! RomP <<<
      LOGICAL, SAVE :: ancien_ok
!$OMP THREADPRIVATE(ancien_ok)
      REAL, ALLOCATABLE, SAVE :: clwcon(:,:),rnebcon(:,:)
!$OMP THREADPRIVATE(clwcon,rnebcon)
      REAL, ALLOCATABLE, SAVE :: rneb_ancien(:,:)
!$OMP THREADPRIVATE(rneb_ancien)
      REAL, ALLOCATABLE, SAVE :: qtc_cv(:,:),sigt_cv(:,:),detrain_cv(:,:),fm_cv(:,:)
!$OMP THREADPRIVATE(qtc_cv,sigt_cv,detrain_cv,fm_cv)
      REAL, ALLOCATABLE, SAVE :: ratqs(:,:)
!$OMP THREADPRIVATE(ratqs)
      REAL, ALLOCATABLE, SAVE :: pbl_tke(:,:,:) ! turb kinetic energy
      REAL, ALLOCATABLE, SAVE :: coefh(:,:,:) ! Kz enthalpie
      REAL, ALLOCATABLE, SAVE :: coefm(:,:,:) ! Kz momentum
!$OMP THREADPRIVATE(pbl_tke, coefh,coefm)
      REAL, ALLOCATABLE, SAVE :: zmax0(:), f0(:) ! 
!$OMP THREADPRIVATE(zmax0,f0)
      REAL, ALLOCATABLE, SAVE :: sig1(:,:), w01(:,:)
!$OMP THREADPRIVATE(sig1,w01)
      REAL, ALLOCATABLE, SAVE :: entr_therm(:,:), fm_therm(:,:)
!$OMP THREADPRIVATE(entr_therm,fm_therm)
      REAL, ALLOCATABLE, SAVE :: detr_therm(:,:)
!$OMP THREADPRIVATE(detr_therm)
!IM 150408
!     pour phsystoke avec thermiques
      REAL,ALLOCATABLE,SAVE :: clwcon0th(:,:),rnebcon0th(:,:)
!$OMP THREADPRIVATE(clwcon0th,rnebcon0th)
! radiation outputs
      REAL,ALLOCATABLE,SAVE :: swdnc0(:,:), swdn0(:,:), swdn(:,:)
!$OMP THREADPRIVATE(swdnc0,swdn0,swdn)
      REAL,ALLOCATABLE,SAVE :: swupc0(:,:), swup0(:,:), swup(:,:)
!$OMP THREADPRIVATE(swupc0, swup0,swup)
      REAL,ALLOCATABLE,SAVE :: SWdn200clr(:), SWdn200(:)
!$OMP THREADPRIVATE(SWdn200clr,SWdn200)
      REAL,ALLOCATABLE,SAVE :: SWup200clr(:), SWup200(:)
!$OMP THREADPRIVATE(SWup200clr,SWup200)
      REAL,ALLOCATABLE,SAVE :: lwdnc0(:,:), lwdn0(:,:), lwdn(:,:)
!$OMP THREADPRIVATE(lwdnc0,lwdn0,lwdn)
      REAL,ALLOCATABLE,SAVE :: lwupc0(:,:), lwup0(:,:), lwup(:,:)
!$OMP THREADPRIVATE(lwupc0,lwup0,lwup)
      REAL,ALLOCATABLE,SAVE :: LWdn200clr(:), LWdn200(:)
!$OMP THREADPRIVATE(LWdn200clr,LWdn200)
      REAL,ALLOCATABLE,SAVE :: LWup200clr(:), LWup200(:)
!$OMP THREADPRIVATE(LWup200clr,LWup200)
      REAL,ALLOCATABLE,SAVE :: LWdnTOA(:), LWdnTOAclr(:)
!$OMP THREADPRIVATE(LWdnTOA,LWdnTOAclr)
! pressure level
      REAL,ALLOCATABLE,SAVE :: tsumSTD(:,:,:)
!$OMP THREADPRIVATE(tsumSTD)
      REAL,ALLOCATABLE,SAVE :: usumSTD(:,:,:), vsumSTD(:,:,:)
!$OMP THREADPRIVATE(usumSTD,vsumSTD)
      REAL,ALLOCATABLE,SAVE :: wsumSTD(:,:,:), phisumSTD(:,:,:)
!$OMP THREADPRIVATE(wsumSTD,phisumSTD)
      REAL,ALLOCATABLE,SAVE :: qsumSTD(:,:,:), rhsumSTD(:,:,:)
!$OMP THREADPRIVATE(qsumSTD,rhsumSTD)
      REAL,ALLOCATABLE,SAVE :: tnondef(:,:,:) 
!$OMP THREADPRIVATE(tnondef)
      REAL,ALLOCATABLE,SAVE :: uvsumSTD(:,:,:)
!$OMP THREADPRIVATE(uvsumSTD)
      REAL,ALLOCATABLE,SAVE :: vqsumSTD(:,:,:)
!$OMP THREADPRIVATE(vqsumSTD)
      REAL,ALLOCATABLE,SAVE :: vTsumSTD(:,:,:)
!$OMP THREADPRIVATE(vTsumSTD)
      REAL,ALLOCATABLE,SAVE :: wqsumSTD(:,:,:)
!$OMP THREADPRIVATE(wqsumSTD)
      REAL,ALLOCATABLE,SAVE :: vphisumSTD(:,:,:)
!$OMP THREADPRIVATE(vphisumSTD)
      REAL,ALLOCATABLE,SAVE :: wTsumSTD(:,:,:)
!$OMP THREADPRIVATE(wTsumSTD)
      REAL,ALLOCATABLE,SAVE :: u2sumSTD(:,:,:)
!$OMP THREADPRIVATE(u2sumSTD)
      REAL,ALLOCATABLE,SAVE :: v2sumSTD(:,:,:)
!$OMP THREADPRIVATE(v2sumSTD)
      REAL,ALLOCATABLE,SAVE :: T2sumSTD(:,:,:)
!$OMP THREADPRIVATE(T2sumSTD)
      REAL,ALLOCATABLE,SAVE :: O3sumSTD(:,:,:), O3daysumSTD(:,:,:)
!$OMP THREADPRIVATE(O3sumSTD,O3daysumSTD) 
!IM begin
      REAL,ALLOCATABLE,SAVE :: wlevSTD(:,:), ulevSTD(:,:), vlevSTD(:,:)
!$OMP THREADPRIVATE(wlevSTD,ulevSTD,vlevSTD)
      REAL,ALLOCATABLE,SAVE :: tlevSTD(:,:), qlevSTD(:,:), rhlevSTD(:,:)
!$OMP THREADPRIVATE(tlevSTD,qlevSTD,rhlevSTD)
      REAL,ALLOCATABLE,SAVE :: philevSTD(:,:)
!$OMP THREADPRIVATE(philevSTD)
      REAL,ALLOCATABLE,SAVE :: uvSTD(:,:)
!$OMP THREADPRIVATE(uvSTD)
      REAL,ALLOCATABLE,SAVE :: vqSTD(:,:)
!$OMP THREADPRIVATE(vqSTD)
      REAL,ALLOCATABLE,SAVE :: vTSTD(:,:)
!$OMP THREADPRIVATE(vTSTD)
      REAL,ALLOCATABLE,SAVE :: wqSTD(:,:)
!$OMP THREADPRIVATE(wqSTD)
      REAL,ALLOCATABLE,SAVE :: vphiSTD(:,:)
!$OMP THREADPRIVATE(vphiSTD)
      REAL,ALLOCATABLE,SAVE :: wTSTD(:,:)
!$OMP THREADPRIVATE(wTSTD)
      REAL,ALLOCATABLE,SAVE :: u2STD(:,:)
!$OMP THREADPRIVATE(u2STD)
      REAL,ALLOCATABLE,SAVE :: v2STD(:,:) 
!$OMP THREADPRIVATE(v2STD)
      REAL,ALLOCATABLE,SAVE :: T2STD(:,:)
!$OMP THREADPRIVATE(T2STD)
      REAL,ALLOCATABLE,SAVE :: O3STD(:,:), O3daySTD(:,:)
!$OMP THREADPRIVATE(O3STD,O3daySTD)
!IM end
      INTEGER,ALLOCATABLE,SAVE :: seed_old(:,:)
!$OMP THREADPRIVATE(seed_old)
      REAL,ALLOCATABLE,SAVE :: zuthe(:),zvthe(:)
!$OMP THREADPRIVATE(zuthe,zvthe)
      REAL,ALLOCATABLE,SAVE :: alb_neig(:)
!$OMP THREADPRIVATE(alb_neig)
!cloud base mass flux
      REAL,ALLOCATABLE,SAVE :: ema_cbmf(:)
!$OMP THREADPRIVATE(ema_cbmf)
!cloud base pressure & cloud top pressure
      REAL,ALLOCATABLE,SAVE :: ema_pcb(:), ema_pct(:)
!$OMP THREADPRIVATE(ema_pcb,ema_pct)
      REAL,ALLOCATABLE,SAVE :: Mipsh(:,:)     ! mass flux shed from  adiab. ascents
!$OMP THREADPRIVATE(Mipsh)
      REAL,ALLOCATABLE,SAVE :: Ma(:,:)       ! undilute upward mass flux
!$OMP THREADPRIVATE(Ma)
      REAL,ALLOCATABLE,SAVE :: qcondc(:,:)    ! in-cld water content from convect
!$OMP THREADPRIVATE(qcondc)
      REAL,ALLOCATABLE,SAVE :: wd(:) ! sb
!$OMP THREADPRIVATE(wd)
      REAL,ALLOCATABLE,SAVE :: sigd(:)
!$OMP THREADPRIVATE(sigd)
!
      REAL,ALLOCATABLE,SAVE :: cin(:)
!$OMP THREADPRIVATE(cin)
! ftd : convective heating due to unsaturated downdraughts
      REAL,ALLOCATABLE,SAVE :: ftd(:,:)
!$OMP THREADPRIVATE(ftd)
! fqd : convective moistening due to unsaturated downdraughts
      REAL,ALLOCATABLE,SAVE :: fqd(:,:),fqcomp(:,:)     
!$OMP THREADPRIVATE(fqd,fqcomp)
#ifdef ISO
      REAL, ALLOCATABLE, SAVE :: fxtd(:,:,:)
!$OMP THREADPRIVATE(fxtd)
#endif
!34EK
! -- Variables de controle de ALE et ALP
!ALE : Energie disponible pour soulevement : utilisee par la 
!      convection d'Emanuel pour le declenchement et la regulation
      REAL,ALLOCATABLE,SAVE :: ALE(:)
!$OMP THREADPRIVATE(ALE)
!ALP : Puissance  disponible pour soulevement
      REAL,ALLOCATABLE,SAVE :: ALP(:)
!$OMP THREADPRIVATE(ALP)
!
! nouvelles variables pour le couplage convection-couche limite
      REAL,ALLOCATABLE,SAVE :: Ale_bl(:)
!$OMP THREADPRIVATE(Ale_bl)
      REAL,ALLOCATABLE,SAVE :: Alp_bl(:)
!$OMP THREADPRIVATE(Alp_bl)
      INTEGER,ALLOCATABLE,SAVE :: lalim_conv(:)
!$OMP THREADPRIVATE(lalim_conv)
      REAL,ALLOCATABLE,SAVE :: wght_th(:,:)
!$OMP THREADPRIVATE(wght_th)
      REAL,ALLOCATABLE,SAVE    :: ale_wake(:)
!$OMP THREADPRIVATE(ale_wake)
      REAL,ALLOCATABLE,SAVE    :: ale_bl_stat(:)
!$OMP THREADPRIVATE(ale_bl_stat)
!
! variables de la wake
! wake_deltat : ecart de temperature avec la zone non perturbee
! wake_deltaq : ecart d'humidite avec la zone non perturbee
! wake_s      : fraction surfacique occupee par la poche froide
! awake_dens  : number of active wakes per unit area
! wake_dens   : number of wakes per unit area
! cv_gen      : birth rate of cumulonimbus per unit area.
! wake_occ    : occurence of wakes (= 1 if wakes occur, =0 otherwise)
! wake_Cstar  : vitesse d'etalement de la poche
! wake_pe     : wake potential energy - WAPE
! wake_fip    : Gust Front Impinging power - ALP
      REAL,ALLOCATABLE,SAVE :: wake_deltat(:,:)
!$OMP THREADPRIVATE(wake_deltat)
      REAL,ALLOCATABLE,SAVE :: wake_deltaq(:,:)
!$OMP THREADPRIVATE(wake_deltaq)
#ifdef ISO
      REAL, ALLOCATABLE, SAVE :: wake_deltaxt(:,:,:)
!$OMP THREADPRIVATE(wake_deltaxt)
#endif
      REAL,ALLOCATABLE,SAVE :: wake_s(:)
!$OMP THREADPRIVATE(wake_s)
      REAL,ALLOCATABLE,SAVE :: awake_dens(:), wake_dens(:)
!$OMP THREADPRIVATE(awake_dens, wake_dens)
      REAL,ALLOCATABLE,SAVE :: cv_gen(:)
!$OMP THREADPRIVATE(cv_gen)
      REAL,ALLOCATABLE,SAVE :: wake_Cstar(:)
!$OMP THREADPRIVATE(wake_Cstar)
      REAL,ALLOCATABLE,SAVE :: wake_pe(:)
!$OMP THREADPRIVATE(wake_pe)
      REAL,ALLOCATABLE,SAVE :: wake_fip(:)
!$OMP THREADPRIVATE(wake_fip)
!
!jyg<
! variables related to the spitting of the PBL between wake and 
! off-wake regions.
! wake_delta_pbl_TKE : difference TKE_w - TKE_x
      REAL,ALLOCATABLE,SAVE :: wake_delta_pbl_TKE(:,:,:)
!$OMP THREADPRIVATE(wake_delta_pbl_TKE)
!nrlmd<
      REAL, ALLOCATABLE, SAVE :: delta_tsurf(:,:) ! Surface temperature difference inside-outside cold pool
!$OMP THREADPRIVATE(delta_tsurf)
!>nrlmd
!>jyg
!
! pfrac_impa : Produits des coefs lessivage impaction
! pfrac_nucl : Produits des coefs lessivage nucleation
! pfrac_1nucl: Produits des coefs lessi nucl (alpha = 1) 
      REAL,ALLOCATABLE,SAVE :: pfrac_impa(:,:), pfrac_nucl(:,:)
!$OMP THREADPRIVATE(pfrac_impa,pfrac_nucl)
      REAL,ALLOCATABLE,SAVE :: pfrac_1nucl(:,:)
!$OMP THREADPRIVATE(pfrac_1nucl)
!
      REAL,ALLOCATABLE,SAVE :: total_rain(:), nday_rain(:)  
!$OMP THREADPRIVATE(total_rain,nday_rain)
      REAL,ALLOCATABLE,SAVE :: ndayrain_mth(:)
!$OMP THREADPRIVATE(ndayrain_mth)
      REAL,ALLOCATABLE,SAVE :: paire_ter(:)
!$OMP THREADPRIVATE(paire_ter)
! albsol1: albedo du sol total pour SW visible
! albsol2: albedo du sol total pour SW proche IR
      REAL,ALLOCATABLE,SAVE :: albsol1(:), albsol2(:)
!$OMP THREADPRIVATE(albsol1,albsol2)

!albedo SB >>>
      REAL,ALLOCATABLE,SAVE :: albsol_dif(:,:),albsol_dir(:,:)
!$OMP THREADPRIVATE(albsol_dif,albsol_dir)
!albedo SB <<<


      REAL, ALLOCATABLE, SAVE:: wo(:, :, :)
      ! column-density of ozone in a layer, in kilo-Dobsons
      ! Third dimension has size 1 or 2.
      ! "wo(:, :, 1)" is for the average day-night field, 
      ! "wo(:, :, 2)" is for daylight time.
      !$OMP THREADPRIVATE(wo)

! heat : chauffage solaire
! heat0: chauffage solaire ciel clair
! cool : refroidissement infrarouge
! cool0 : refroidissement infrarouge ciel clair
! sollwdown : downward LW flux at surface
! sollwdownclr : downward CS LW flux at surface
! toplwdown : downward CS LW flux at TOA
! toplwdownclr : downward CS LW flux at TOA
! heat_volc : chauffage solaire du au volcanisme
! cool_volc : refroidissement infrarouge du au volcanisme
      REAL,ALLOCATABLE,SAVE :: clwcon0(:,:),rnebcon0(:,:)
!$OMP THREADPRIVATE(clwcon0,rnebcon0)
      REAL,ALLOCATABLE,SAVE :: heat(:,:)   
!$OMP THREADPRIVATE(heat)
      REAL,ALLOCATABLE,SAVE :: heat0(:,:)
!$OMP THREADPRIVATE(heat0)
      REAL,ALLOCATABLE,SAVE :: cool(:,:)
!$OMP THREADPRIVATE(cool)
      REAL,ALLOCATABLE,SAVE :: cool0(:,:)
!$OMP THREADPRIVATE(cool0)
      REAL,ALLOCATABLE,SAVE :: heat_volc(:,:)   
!$OMP THREADPRIVATE(heat_volc)
      REAL,ALLOCATABLE,SAVE :: cool_volc(:,:)
!$OMP THREADPRIVATE(cool_volc)
      REAL,ALLOCATABLE,SAVE :: topsw(:), toplw(:)
!$OMP THREADPRIVATE(topsw,toplw)
      REAL,ALLOCATABLE,SAVE :: sollwdown(:)
!$OMP THREADPRIVATE(sollwdown)
      REAL,ALLOCATABLE,SAVE :: gustiness(:)
!$OMP THREADPRIVATE(gustiness)
      REAL,ALLOCATABLE,SAVE :: sollwdownclr(:)
!$OMP THREADPRIVATE(sollwdownclr)
      REAL,ALLOCATABLE,SAVE :: toplwdown(:)
!$OMP THREADPRIVATE(toplwdown)
      REAL,ALLOCATABLE,SAVE :: toplwdownclr(:)
!$OMP THREADPRIVATE(toplwdownclr)
      REAL,ALLOCATABLE,SAVE :: topsw0(:),toplw0(:),solsw0(:),sollw0(:)
!$OMP THREADPRIVATE(topsw0,toplw0,solsw0,sollw0)
      REAL,ALLOCATABLE,SAVE :: albpla(:)
!$OMP THREADPRIVATE(albpla)

!IM ajout variables CFMIP2/CMIP5
      REAL,ALLOCATABLE,SAVE :: heatp(:,:), coolp(:,:)
!$OMP THREADPRIVATE(heatp, coolp)
      REAL,ALLOCATABLE,SAVE :: heat0p(:,:), cool0p(:,:)
!$OMP THREADPRIVATE(heat0p, cool0p)
      REAL,ALLOCATABLE,SAVE :: radsolp(:), topswp(:), toplwp(:)
!$OMP THREADPRIVATE(radsolp, topswp, toplwp)
      REAL,ALLOCATABLE,SAVE :: albplap(:)
!$OMP THREADPRIVATE(albplap)
      REAL,ALLOCATABLE,SAVE :: solswp(:), solswfdiffp(:), sollwp(:)
!$OMP THREADPRIVATE(solswp, solswfdiffp, sollwp)
      REAL,ALLOCATABLE,SAVE :: sollwdownp(:)
!$OMP THREADPRIVATE(sollwdownp)
      REAL,ALLOCATABLE,SAVE :: topsw0p(:),toplw0p(:)
      REAL,ALLOCATABLE,SAVE :: solsw0p(:),sollw0p(:)
!$OMP THREADPRIVATE(topsw0p,toplw0p,solsw0p,sollw0p)
      REAL,ALLOCATABLE,SAVE :: lwdnc0p(:,:), lwdn0p(:,:), lwdnp(:,:)
      REAL,ALLOCATABLE,SAVE :: lwupc0p(:,:), lwup0p(:,:), lwupp(:,:)
!$OMP THREADPRIVATE(lwdnc0p, lwdn0p, lwdnp, lwupc0p, lwup0p, lwupp)
      REAL,ALLOCATABLE,SAVE :: swdnc0p(:,:), swdn0p(:,:), swdnp(:,:)
      REAL,ALLOCATABLE,SAVE :: swupc0p(:,:), swup0p(:,:), swupp(:,:)
!$OMP THREADPRIVATE(swdnc0p, swdn0p, swdnp, swupc0p, swup0p, swupp)

!AI ajout variables double appel Ecrad (3Deffect)
      REAL,ALLOCATABLE,SAVE :: heat_s2(:,:), cool_s2(:,:)
!$OMP THREADPRIVATE(heat_s2, cool_s2)
      REAL,ALLOCATABLE,SAVE :: heat0_s2(:,:), cool0_s2(:,:)
!$OMP THREADPRIVATE(heat0_s2, cool0_s2)
      REAL,ALLOCATABLE,SAVE :: radsol_s2(:), topsw_s2(:), toplw_s2(:)
!$OMP THREADPRIVATE(radsol_s2, topsw_s2, toplw_s2)
      REAL,ALLOCATABLE,SAVE :: albpla_s2(:)
!$OMP THREADPRIVATE(albpla_s2)
      REAL,ALLOCATABLE,SAVE :: solsw_s2(:), solswfdiff_s2(:), sollw_s2(:)
!$OMP THREADPRIVATE(solsw_s2, solswfdiff_s2, sollw_s2)
      REAL,ALLOCATABLE,SAVE :: sollwdown_s2(:)
!$OMP THREADPRIVATE(sollwdown_s2)
      REAL,ALLOCATABLE,SAVE :: topsw0_s2(:),toplw0_s2(:)
      REAL,ALLOCATABLE,SAVE :: solsw0_s2(:),sollw0_s2(:)
!$OMP THREADPRIVATE(topsw0_s2,toplw0_s2,solsw0_s2,sollw0_s2)
      REAL,ALLOCATABLE,SAVE :: lwdnc0_s2(:,:), lwdn0_s2(:,:), lwdn_s2(:,:)
      REAL,ALLOCATABLE,SAVE :: lwupc0_s2(:,:), lwup0_s2(:,:), lwup_s2(:,:)
!$OMP THREADPRIVATE(lwdnc0_s2,lwdn0_s2,lwdn_s2,lwupc0_s2,lwup0_s2,lwup_s2)       
      REAL,ALLOCATABLE,SAVE :: swdnc0_s2(:,:), swdn0_s2(:,:), swdn_s2(:,:)
      REAL,ALLOCATABLE,SAVE :: swupc0_s2(:,:), swup0_s2(:,:), swup_s2(:,:)
!$OMP THREADPRIVATE(swdnc0_s2, swdn0_s2, swdn_s2, swupc0_s2, swup0_s2, swup_s2)

! pbase : cloud base pressure
! bbase : cloud base buoyancy
      REAL,ALLOCATABLE,SAVE :: cape(:)
!$OMP THREADPRIVATE(cape)
      REAL,ALLOCATABLE,SAVE :: pbase(:)
!$OMP THREADPRIVATE(pbase)
      REAL,ALLOCATABLE,SAVE :: bbase(:)
!$OMP THREADPRIVATE(bbase)
!
      REAL,SAVE,ALLOCATABLE :: zqasc(:,:)
!$OMP THREADPRIVATE( zqasc)
      INTEGER,ALLOCATABLE,SAVE :: ibas_con(:), itop_con(:)
!$OMP THREADPRIVATE(ibas_con,itop_con)
      REAL,SAVE,ALLOCATABLE :: rain_con(:)
!$OMP THREADPRIVATE(rain_con)
      REAL,SAVE,ALLOCATABLE :: snow_con(:)
!$OMP THREADPRIVATE(snow_con)
!
#ifdef ISO
      REAL,SAVE,ALLOCATABLE :: xtrain_con(:,:)
!$OMP THREADPRIVATE(xtrain_con)
      REAL,SAVE,ALLOCATABLE :: xtsnow_con(:,:)
!$OMP THREADPRIVATE(xtsnow_con)
#endif
      REAL,SAVE,ALLOCATABLE :: rlonPOS(:)
!$OMP THREADPRIVATE(rlonPOS)
      REAL,SAVE,ALLOCATABLE :: newsst(:)
!$OMP THREADPRIVATE(newsst)
      REAL,SAVE,ALLOCATABLE :: ustar(:,:),u10m(:,:), v10m(:,:),wstar(:,:)
!$OMP THREADPRIVATE(ustar,u10m,v10m,wstar)
!
! ok_ade=T -ADE=topswad-topsw
! ok_aie=T ->
!       ok_ade=T -AIE=topswai-topswad
!       ok_ade=F -AIE=topswai-topsw
!
!topswad, solswad : Aerosol direct effect
      REAL,SAVE,ALLOCATABLE :: topswad(:), solswad(:)
!$OMP THREADPRIVATE(topswad,solswad)
!topswai, solswai : Aerosol indirect effect
      REAL,SAVE,ALLOCATABLE :: topswai(:), solswai(:)
!$OMP THREADPRIVATE(topswai,solswai)

      REAL,SAVE,ALLOCATABLE :: tau_aero(:,:,:,:), piz_aero(:,:,:,:), cg_aero(:,:,:,:)
!$OMP THREADPRIVATE(tau_aero, piz_aero, cg_aero)
      REAL,SAVE,ALLOCATABLE :: tau_aero_sw_rrtm(:,:,:,:), piz_aero_sw_rrtm(:,:,:,:), cg_aero_sw_rrtm(:,:,:,:)
!$OMP THREADPRIVATE(tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm)
      REAL,SAVE,ALLOCATABLE :: tau_aero_lw_rrtm(:,:,:,:), piz_aero_lw_rrtm(:,:,:,:), cg_aero_lw_rrtm(:,:,:,:)
!$OMP THREADPRIVATE(tau_aero_lw_rrtm, piz_aero_lw_rrtm, cg_aero_lw_rrtm)
      REAL,SAVE,ALLOCATABLE :: ccm(:,:,:)
!$OMP THREADPRIVATE(ccm)

      REAL,SAVE,ALLOCATABLE :: ale_bl_trig(:)
!$OMP THREADPRIVATE(ale_bl_trig)

      REAL,SAVE,ALLOCATABLE :: ratqs_inter_(:,:)
!$OMP THREADPRIVATE(ratqs_inter_)

#ifdef ISO
#ifdef ISOTRAC
      INTEGER,SAVE,ALLOCATABLE :: bassin_map(:)
!$OMP THREADPRIVATE(bassin_map)
      INTEGER,SAVE,ALLOCATABLE :: boite_map(:,:)
!$OMP THREADPRIVATE(boite_map)
#endif   
#endif
      REAL, ALLOCATABLE, SAVE:: du_gwd_rando(:, :), du_gwd_front(:, :)
      !$OMP THREADPRIVATE(du_gwd_rando, du_gwd_front)
      ! tendencies on wind due to gravity waves

      LOGICAL,SAVE :: is_initialized=.FALSE.
!$OMP THREADPRIVATE(is_initialized)    

      ! Ocean-atmosphere interface:

      REAL, ALLOCATABLE, SAVE:: ds_ns(:) ! (klon)
      ! "delta salinity near surface". Salinity variation in the
      ! near-surface turbulent layer. That is subskin salinity minus
      ! foundation salinity. In ppt.

      REAL, ALLOCATABLE, SAVE:: dt_ns(:) ! (klon)
      ! "delta temperature near surface". Temperature variation in the
      ! near-surface turbulent layer. That is subskin temperature
      ! minus foundation temperature. (Can be negative.) In K.
      
      REAL, ALLOCATABLE, SAVE:: delta_sst(:) ! (klon)
      ! Ocean-air interface temperature minus bulk SST, in
      ! K. Allocated and defined only if activate_ocean_skin >= 1.

      REAL, ALLOCATABLE, SAVE:: delta_sal(:) ! (klon)
      ! Ocean-air interface salinity minus bulk salinity, in ppt
      
      REAL, ALLOCATABLE, SAVE:: dter(:) ! (klon)
      ! Temperature variation in the diffusive microlayer, that is
      ! ocean-air interface temperature minus subskin temperature. In K.

      REAL, SAVE, ALLOCATABLE:: dser(:) ! (klon)
      ! Salinity variation in the diffusive microlayer, that is
      ! ocean-air interface salinity minus subskin salinity. In ppt.

      real, SAVE, ALLOCATABLE:: dt_ds(:) ! (klon)
      ! (tks / tkt) * dTer, in K

      !$OMP THREADPRIVATE(delta_sal, ds_ns, dt_ns, delta_sst, dter, dser, dt_ds)

    CONTAINS

!======================================================================
SUBROUTINE phys_state_var_init(read_climoz)
USE dimphy
USE aero_mod
USE infotrac_phy, ONLY : nbtr
#ifdef ISO
USE infotrac_phy, ONLY : ntraciso=>ntiso,niso
#endif
USE indice_sol_mod
use config_ocean_skin_m, only: activate_ocean_skin
use surface_data, only: type_ocean
IMPLICIT NONE

integer, intent(in)::  read_climoz
! read ozone climatology
! Allowed values are 0, 1 and 2
! 0: do not read an ozone climatology
! 1: read a single ozone climatology that will be used day and night
! 2: read two ozone climatologies, the average day and night
! climatology and the daylight climatology

include "clesphys.h"

      print*, 'is_initialized', is_initialized
      IF (is_initialized) RETURN
      is_initialized=.TRUE.
      ALLOCATE(pctsrf(klon,nbsrf))
      ALLOCATE(ftsol(klon,nbsrf))
      ALLOCATE(beta_aridity(klon,nbsrf))
      ALLOCATE(qsol(klon),fevap(klon,nbsrf))
      ALLOCATE(z0m(klon,nbsrf+1),z0h(klon,nbsrf+1),agesno(klon,nbsrf))
!FC
      ALLOCATE(treedrg(klon,klev,nbsrf))
      ALLOCATE(falb1(klon,nbsrf))
      ALLOCATE(falb2(klon,nbsrf))
!albedo SB >>> 
      print*, 'allocate falb'
      ALLOCATE(falb_dir(klon,nsw,nbsrf),falb_dif(klon,nsw,nbsrf))
!!      print*, 'allocate falb good', falb_dir(1,1,1)
      ALLOCATE(chl_con(klon))
!albedo SB <<<
      ALLOCATE(rain_fall(klon))
      ALLOCATE(snow_fall(klon))
      ALLOCATE(bs_fall(klon))
      ALLOCATE(solsw(klon), solswfdiff(klon), sollw(klon))
      sollw=0.0
      ALLOCATE(radsol(klon))
      ALLOCATE(swradcorr(klon))
      ALLOCATE(zmea(klon), zstd(klon), zsig(klon), zgam(klon))
      ALLOCATE(zthe(klon), zpic(klon), zval(klon))

      ALLOCATE(rugoro(klon))
      ALLOCATE(t_ancien(klon,klev), q_ancien(klon,klev))
      ALLOCATE(ql_ancien(klon,klev), qs_ancien(klon,klev), qbs_ancien(klon,klev))
      ALLOCATE(prw_ancien(klon), prlw_ancien(klon), prsw_ancien(klon), prbsw_ancien(klon))
      ALLOCATE(u_ancien(klon,klev), v_ancien(klon,klev))
!!! Rom P >>>
      ALLOCATE(tr_ancien(klon,klev,nbtr))
!!! Rom P <<<
      ALLOCATE(clwcon(klon,klev),rnebcon(klon,klev))
      ALLOCATE(rneb_ancien(klon,klev))
      ALLOCATE(qtc_cv(klon,klev),sigt_cv(klon,klev),detrain_cv(klon,klev),fm_cv(klon,klev))
      ALLOCATE(ratqs(klon,klev))
      ALLOCATE(pbl_tke(klon,klev+1,nbsrf+1))
!nrlmd<
      ALLOCATE(delta_tsurf(klon,nbsrf))
!>nrlmd
      ALLOCATE(coefh(klon,klev+1,nbsrf+1))
      ALLOCATE(coefm(klon,klev+1,nbsrf+1))
      ! initialize cleanly coefh,coefm
      ! (most of the time in the code these are assumed to be on klev levels)
      coefh(:,:,:)=0
      coefm(:,:,:)=0
      ALLOCATE(zmax0(klon), f0(klon))
      ALLOCATE(sig1(klon,klev), w01(klon,klev))
      ALLOCATE(entr_therm(klon,klev), fm_therm(klon,klev+1))
      ALLOCATE(detr_therm(klon,klev))
!     pour phsystoke avec thermiques
      ALLOCATE(clwcon0th(klon,klev),rnebcon0th(klon,klev))
! radiation outputs
      ALLOCATE(swdnc0(klon,klevp1), swdn0(klon,klevp1), swdn(klon,klevp1))
      ALLOCATE(swupc0(klon,klevp1), swup0(klon,klevp1), swup(klon,klevp1))
      ALLOCATE(lwdnc0(klon,klevp1), lwdn0(klon,klevp1), lwdn(klon,klevp1))
      ALLOCATE(lwupc0(klon,klevp1), lwup0(klon,klevp1), lwup(klon,klevp1))
      ALLOCATE(SWdn200clr(klon), SWdn200(klon))
      ALLOCATE(SWup200clr(klon), SWup200(klon))
      ALLOCATE(LWdn200clr(klon), LWdn200(klon))
      ALLOCATE(LWup200clr(klon), LWup200(klon))
      ALLOCATE(LWdnTOA(klon), LWdnTOAclr(klon))
! pressure level
      ALLOCATE(tsumSTD(klon,nlevSTD,nout))
      ALLOCATE(usumSTD(klon,nlevSTD,nout), vsumSTD(klon,nlevSTD,nout))
      ALLOCATE(wsumSTD(klon,nlevSTD,nout), phisumSTD(klon,nlevSTD,nout))
      ALLOCATE(qsumSTD(klon,nlevSTD,nout), rhsumSTD(klon,nlevSTD,nout))
      ALLOCATE(tnondef(klon,nlevSTD,nout))
      ALLOCATE(uvsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vqsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vTsumSTD(klon,nlevSTD,nout))
      ALLOCATE(wqsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vphisumSTD(klon,nlevSTD,nout))
      ALLOCATE(wTsumSTD(klon,nlevSTD,nout))
      ALLOCATE(u2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(v2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(T2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(O3sumSTD(klon,nlevSTD,nout))
      ALLOCATE(O3daysumSTD(klon,nlevSTD,nout))
!IM beg
      ALLOCATE(wlevSTD(klon,nlevSTD), ulevSTD(klon,nlevSTD), vlevSTD(klon,nlevSTD))
      ALLOCATE(tlevSTD(klon,nlevSTD), qlevSTD(klon,nlevSTD), rhlevSTD(klon,nlevSTD))
      ALLOCATE(philevSTD(klon,nlevSTD))
      ALLOCATE(uvSTD(klon,nlevSTD),vqSTD(klon,nlevSTD))
      ALLOCATE(vTSTD(klon,nlevSTD),wqSTD(klon,nlevSTD))
      ALLOCATE(vphiSTD(klon,nlevSTD),wTSTD(klon,nlevSTD))
      ALLOCATE(u2STD(klon,nlevSTD),v2STD(klon,nlevSTD))
      ALLOCATE(T2STD(klon,nlevSTD))
      ALLOCATE(O3STD(klon,nlevSTD))
      ALLOCATE(O3daySTD(klon,nlevSTD))
!IM end
      ALLOCATE(seed_old(klon,napisccp))
      ALLOCATE(zuthe(klon),zvthe(klon))
      ALLOCATE(alb_neig(klon))
!cloud base mass flux
      ALLOCATE(ema_cbmf(klon))
!cloud base pressure & cloud top pressure
      ALLOCATE(ema_pcb(klon), ema_pct(klon))
!
      ALLOCATE(Mipsh(klon,klev))
      ALLOCATE(Ma(klon,klev))
      ALLOCATE(qcondc(klon,klev))
      ALLOCATE(wd(klon))
      ALLOCATE(sigd(klon))
      ALLOCATE(cin(klon), ALE(klon), ALP(klon))
      ALLOCATE(ftd(klon,klev), fqd(klon,klev),fqcomp(klon,klev))
      ALLOCATE(Ale_bl(klon))
      ALLOCATE(ale_wake(klon))
      ALLOCATE(ale_bl_stat(klon))
      ale_bl_stat(:)=0
      ALLOCATE(Alp_bl(klon))
      ALLOCATE(lalim_conv(klon))
      ALLOCATE(wght_th(klon,klev))
      ALLOCATE(wake_deltat(klon,klev), wake_deltaq(klon,klev))
      ALLOCATE(wake_s(klon), awake_dens(klon), wake_dens(klon))
!!      awake_dens = 0.  ! initialized in phyetat0
      ALLOCATE(cv_gen(klon))
      ALLOCATE(wake_Cstar(klon))
      ALLOCATE(wake_pe(klon), wake_fip(klon))
!jyg<
      ALLOCATE(wake_delta_pbl_TKE(klon,klev+1,nbsrf+1))
!>jyg
      ALLOCATE(pfrac_impa(klon,klev), pfrac_nucl(klon,klev))
      ALLOCATE(pfrac_1nucl(klon,klev))
      ALLOCATE(total_rain(klon), nday_rain(klon))
      ALLOCATE(ndayrain_mth(klon))
      ALLOCATE(paire_ter(klon))
      ALLOCATE(albsol1(klon), albsol2(klon))
!albedo SB >>>
      ALLOCATE(albsol_dir(klon,nsw),albsol_dif(klon,nsw))
!albedo SB <<<

      if (read_climoz <= 1) then
         ALLOCATE(wo(klon,klev, 1))
      else
         ! read_climoz == 2
         ALLOCATE(wo(klon,klev, 2))
      end if
      
      ALLOCATE(clwcon0(klon,klev),rnebcon0(klon,klev))
      ALLOCATE(heat(klon,klev), heat0(klon,klev)) 
      ALLOCATE(cool(klon,klev), cool0(klon,klev))
      ALLOCATE(heat_volc(klon,klev), cool_volc(klon,klev)) 
      ALLOCATE(topsw(klon), toplw(klon))
      ALLOCATE(sollwdown(klon), sollwdownclr(klon))
      sollwdown = 0.
      ALLOCATE(toplwdown(klon), toplwdownclr(klon))
      ALLOCATE(topsw0(klon),toplw0(klon),solsw0(klon),sollw0(klon))
      sollw0 = 0.
      ALLOCATE(albpla(klon))
!IM ajout variables CFMIP2/CMIP5
      ALLOCATE(heatp(klon,klev), coolp(klon,klev))
      ALLOCATE(heat0p(klon,klev), cool0p(klon,klev))
      ALLOCATE(radsolp(klon), topswp(klon), toplwp(klon))
      ALLOCATE(albplap(klon))
      ALLOCATE(solswp(klon), solswfdiffp(klon), sollwp(klon))
      ALLOCATE(gustiness(klon))
      ALLOCATE(sollwdownp(klon))
      ALLOCATE(topsw0p(klon),toplw0p(klon))
      ALLOCATE(solsw0p(klon),sollw0p(klon))
      ALLOCATE(lwdnc0p(klon,klevp1), lwdn0p(klon,klevp1), lwdnp(klon,klevp1))
      ALLOCATE(lwupc0p(klon,klevp1), lwup0p(klon,klevp1), lwupp(klon,klevp1))
      ALLOCATE(swdnc0p(klon,klevp1), swdn0p(klon,klevp1), swdnp(klon,klevp1))
      ALLOCATE(swupc0p(klon,klevp1), swup0p(klon,klevp1), swupp(klon,klevp1))

!AI Ajout pour Ecrad (3Deffect)       
      ALLOCATE(heat_s2(klon,klev), cool_s2(klon,klev))
      ALLOCATE(heat0_s2(klon,klev), cool0_s2(klon,klev))
      ALLOCATE(radsol_s2(klon), topsw_s2(klon), toplw_s2(klon))
      ALLOCATE(albpla_s2(klon))
      ALLOCATE(solsw_s2(klon), solswfdiff_s2(klon), sollw_s2(klon))
      ALLOCATE(sollwdown_s2(klon))
      ALLOCATE(topsw0_s2(klon),toplw0_s2(klon))
      ALLOCATE(solsw0_s2(klon),sollw0_s2(klon))
      ALLOCATE(lwdnc0_s2(klon,klevp1), lwdn0_s2(klon,klevp1), lwdn_s2(klon,klevp1))
      ALLOCATE(lwupc0_s2(klon,klevp1), lwup0_s2(klon,klevp1), lwup_s2(klon,klevp1))
      ALLOCATE(swdnc0_s2(klon,klevp1), swdn0_s2(klon,klevp1), swdn_s2(klon,klevp1))
      ALLOCATE(swupc0_s2(klon,klevp1), swup0_s2(klon,klevp1), swup_s2(klon,klevp1))

      ALLOCATE(cape(klon))
      ALLOCATE(pbase(klon),bbase(klon))
      ALLOCATE(zqasc(klon,klev))
      ALLOCATE(ibas_con(klon), itop_con(klon))
      ALLOCATE(rain_con(klon), snow_con(klon))
      ALLOCATE(rlonPOS(klon))
      ALLOCATE(newsst(klon))
      ALLOCATE(ustar(klon,nbsrf),u10m(klon,nbsrf), v10m(klon,nbsrf),wstar(klon,nbsrf+1))
      ALLOCATE(topswad(klon), solswad(klon))
      ALLOCATE(topswai(klon), solswai(klon))
      ALLOCATE(tau_aero(klon,klev,naero_grp,nbands),piz_aero(klon,klev,naero_grp,nbands),cg_aero(klon,klev,naero_grp,nbands))
      ALLOCATE(tau_aero_sw_rrtm(klon,klev,2,nbands_sw_rrtm),piz_aero_sw_rrtm(klon,klev,2,nbands_sw_rrtm))
      ALLOCATE(cg_aero_sw_rrtm(klon,klev,2,nbands_sw_rrtm))
      ALLOCATE(tau_aero_lw_rrtm(klon,klev,2,nbands_lw_rrtm),piz_aero_lw_rrtm(klon,klev,2,nbands_lw_rrtm))
      ALLOCATE(cg_aero_lw_rrtm(klon,klev,2,nbands_lw_rrtm))
      ALLOCATE(ccm(klon,klev,nbands))

#ifdef ISO
      ALLOCATE(xtsol(niso,klon),fxtevap(ntraciso,klon,nbsrf))
      ALLOCATE(fxtd(ntraciso,klon,klev))
      ALLOCATE(wake_deltaxt(ntraciso,klon,klev))
      ALLOCATE(xt_ancien(ntraciso,klon,klev))
      ALLOCATE(xtl_ancien(ntraciso,klon,klev))
      ALLOCATE(xts_ancien(ntraciso,klon,klev))
      ALLOCATE(xtrain_fall(ntraciso,klon))
      ALLOCATE(xtsnow_fall(ntraciso,klon))
      ALLOCATE(xtrain_con(ntraciso,klon))
      ALLOCATE(xtsnow_con(ntraciso,klon))
#ifdef ISOTRAC
      ALLOCATE(bassin_map(klon))
      ALLOCATE(boite_map(klon,klev))  
#endif      
#endif

      ALLOCATE(ale_bl_trig(klon))
      ALLOCATE(ratqs_inter_(klon,klev))
      IF (ok_gwd_rando) THEN
        ALLOCATE(du_gwd_rando(klon, klev))
        du_gwd_rando(:,:)=0.
      ENDIF
      IF (.not. ok_hines .and. ok_gwd_rando) THEN
        ALLOCATE(du_gwd_front(klon, klev))
        du_gwd_front(:,:) = 0 !ym missing init   
      ENDIF

      if (activate_ocean_skin >= 1) then
         ALLOCATE(delta_sal(klon), ds_ns(klon), dt_ns(klon), delta_sst(klon), &
              dter(klon), dser(klon))
         if (activate_ocean_skin == 2 .and. type_ocean == "couple") &
              allocate(dt_ds(klon))
      end if

    END SUBROUTINE phys_state_var_init

!======================================================================
    SUBROUTINE phys_state_var_end
      ! Useful only for lmdz1d.
!USE dimphy
USE indice_sol_mod
use config_ocean_skin_m, only: activate_ocean_skin
use surface_data, only: type_ocean
IMPLICIT NONE
include "clesphys.h"

      DEALLOCATE(pctsrf, ftsol, falb1, falb2)
      DEALLOCATE(beta_aridity)
      DEALLOCATE(qsol,fevap,z0m,z0h,agesno)
!FC
      DEALLOCATE(treedrg)
      DEALLOCATE(rain_fall, snow_fall, bs_fall,solsw, solswfdiff, sollw, radsol, swradcorr)
      DEALLOCATE(zmea, zstd, zsig, zgam)
      DEALLOCATE(zthe, zpic, zval)
      DEALLOCATE(rugoro, t_ancien, q_ancien, clwcon, rnebcon)
      DEALLOCATE(qs_ancien, ql_ancien, qbs_ancien, rneb_ancien)
      DEALLOCATE(prw_ancien, prlw_ancien, prsw_ancien, prbsw_ancien)
      DEALLOCATE(qtc_cv,sigt_cv,detrain_cv,fm_cv)
      DEALLOCATE(u_ancien, v_ancien)
      DEALLOCATE(tr_ancien)                           !RomP
      DEALLOCATE(ratqs, pbl_tke,coefh,coefm)
      DEALLOCATE(zmax0, f0)
      DEALLOCATE(sig1, w01)
      DEALLOCATE(entr_therm, fm_therm)
      DEALLOCATE(detr_therm)
      DEALLOCATE(clwcon0th, rnebcon0th)
! radiation outputs
      DEALLOCATE(swdnc0, swdn0, swdn)
      DEALLOCATE(swupc0, swup0, swup)
      DEALLOCATE(lwdnc0, lwdn0, lwdn)
      DEALLOCATE(lwupc0, lwup0, lwup)
      DEALLOCATE(SWdn200clr, SWdn200)
      DEALLOCATE(SWup200clr, SWup200)
      DEALLOCATE(LWdn200clr, LWdn200)
      DEALLOCATE(LWup200clr, LWup200)
      DEALLOCATE(LWdnTOA, LWdnTOAclr)
! pressure level
      DEALLOCATE(tsumSTD)
      DEALLOCATE(usumSTD, vsumSTD)
      DEALLOCATE(wsumSTD, phisumSTD)
      DEALLOCATE(tnondef)
      DEALLOCATE(qsumSTD, rhsumSTD)
      DEALLOCATE(uvsumSTD)
      DEALLOCATE(vqsumSTD)
      DEALLOCATE(vTsumSTD)
      DEALLOCATE(wqsumSTD)
      DEALLOCATE(vphisumSTD)
      DEALLOCATE(wTsumSTD)
      DEALLOCATE(u2sumSTD)
      DEALLOCATE(v2sumSTD)
      DEALLOCATE(T2sumSTD)
      DEALLOCATE(O3sumSTD)
      DEALLOCATE(O3daysumSTD)
!IM beg
      DEALLOCATE(wlevSTD,ulevSTD,vlevSTD,tlevSTD,qlevSTD,rhlevSTD,philevSTD)
      DEALLOCATE(uvSTD,vqSTD,vTSTD,wqSTD,vphiSTD,wTSTD,u2STD,v2STD,T2STD,O3STD,O3daySTD)
!IM end
      DEALLOCATE(seed_old)
      DEALLOCATE(zuthe, zvthe)
      DEALLOCATE(alb_neig)
      DEALLOCATE(ema_cbmf)
      DEALLOCATE(ema_pcb, ema_pct)
      DEALLOCATE(Mipsh, Ma, qcondc)
      DEALLOCATE(wd, sigd)
      DEALLOCATE(cin, ALE, ALP)
      DEALLOCATE(ftd, fqd, fqcomp)
      DEALLOCATE(Ale_bl, Alp_bl)
      DEALLOCATE(ale_wake)
      DEALLOCATE(ale_bl_stat)
      DEALLOCATE(lalim_conv, wght_th)
      DEALLOCATE(wake_deltat, wake_deltaq)
      DEALLOCATE(wake_s, awake_dens, wake_dens)
      DEALLOCATE(cv_gen)
      DEALLOCATE(wake_Cstar, wake_pe, wake_fip)
!jyg<
      DEALLOCATE(wake_delta_pbl_TKE)
!nrlmd<
      DEALLOCATE(delta_tsurf)
!>nrlmd
!>jyg
      DEALLOCATE(pfrac_impa, pfrac_nucl)
      DEALLOCATE(pfrac_1nucl)
      DEALLOCATE(total_rain, nday_rain)
      DEALLOCATE(ndayrain_mth)
      DEALLOCATE(paire_ter)
      DEALLOCATE(albsol1, albsol2)
!albedo SB >>>
      DEALLOCATE(albsol_dir,albsol_dif,falb_dir,falb_dif,chl_con)
!albedo SB <<<
      DEALLOCATE(wo)
      DEALLOCATE(clwcon0,rnebcon0)
      DEALLOCATE(heat, heat0) 
      DEALLOCATE(cool, cool0)
      DEALLOCATE(heat_volc, cool_volc) 
      DEALLOCATE(topsw, toplw)
      DEALLOCATE(sollwdown, sollwdownclr)
      DEALLOCATE(gustiness)
      DEALLOCATE(toplwdown, toplwdownclr)
      DEALLOCATE(topsw0,toplw0,solsw0,sollw0)
      DEALLOCATE(albpla)

!AI Ajout pour Ecrad (3Deffect)
      DEALLOCATE(heat_s2, cool_s2)
      DEALLOCATE(heat0_s2, cool0_s2)
      DEALLOCATE(radsol_s2, topsw_s2, toplw_s2)
      DEALLOCATE(albpla_s2)
      DEALLOCATE(solsw_s2, solswfdiff_s2, sollw_s2)
      DEALLOCATE(sollwdown_s2)
      DEALLOCATE(topsw0_s2,toplw0_s2)
      DEALLOCATE(solsw0_s2,sollw0_s2)
      DEALLOCATE(lwdnc0_s2, lwdn0_s2, lwdn_s2)
      DEALLOCATE(lwupc0_s2, lwup0_s2, lwup_s2)
      DEALLOCATE(swdnc0_s2, swdn0_s2, swdn_s2)
      DEALLOCATE(swupc0_s2, swup0_s2, swup_s2)

!IM ajout variables CFMIP2/CMIP5
      DEALLOCATE(heatp, coolp)
      DEALLOCATE(heat0p, cool0p)
      DEALLOCATE(radsolp, topswp, toplwp)
      DEALLOCATE(albplap)
      DEALLOCATE(solswp, solswfdiffp, sollwp)
      DEALLOCATE(sollwdownp)
      DEALLOCATE(topsw0p,toplw0p)
      DEALLOCATE(solsw0p,sollw0p)
      DEALLOCATE(lwdnc0p, lwdn0p, lwdnp)
      DEALLOCATE(lwupc0p, lwup0p, lwupp)
      DEALLOCATE(swdnc0p, swdn0p, swdnp)
      DEALLOCATE(swupc0p, swup0p, swupp)
      DEALLOCATE(cape)
      DEALLOCATE(pbase,bbase)
      DEALLOCATE(zqasc)
      DEALLOCATE(ibas_con, itop_con)
      DEALLOCATE(rain_con, snow_con)
      DEALLOCATE(rlonPOS)
      DEALLOCATE(newsst)
      DEALLOCATE(ustar,u10m, v10m,wstar)
      DEALLOCATE(topswad, solswad)
      DEALLOCATE(topswai, solswai)
      DEALLOCATE(tau_aero,piz_aero,cg_aero)
      DEALLOCATE(tau_aero_sw_rrtm,piz_aero_sw_rrtm,cg_aero_sw_rrtm)
      DEALLOCATE(tau_aero_lw_rrtm,piz_aero_lw_rrtm,cg_aero_lw_rrtm)
      DEALLOCATE(ccm)
      if (ok_gwd_rando) DEALLOCATE(du_gwd_rando)
      if (.not. ok_hines .and. ok_gwd_rando) DEALLOCATE(du_gwd_front)
      DEALLOCATE(ale_bl_trig)
      DEALLOCATE(ratqs_inter_)

      if (activate_ocean_skin >= 1) then
         deALLOCATE(delta_sal, ds_ns, dt_ns, delta_sst, dter, dser)
         if (activate_ocean_skin == 2 .and. type_ocean == "couple") &
              deALLOCATE(dt_ds)
      end if

#ifdef ISO    
      DEALLOCATE(xtsol,fxtevap)  
      DEALLOCATE(xt_ancien,xtl_ancien,xts_ancien, fxtd, wake_deltaxt)
      DEALLOCATE(xtrain_fall, xtsnow_fall, xtrain_con, xtsnow_con)
#ifdef ISOTRAC 
      DEALLOCATE(bassin_map,boite_map) 
#endif        
#endif
      is_initialized=.FALSE.
      
END SUBROUTINE phys_state_var_end

      END MODULE phys_state_var_mod
