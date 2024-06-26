!$Id $
!
MODULE traclmdz_mod

! 
! In this module all tracers specific to LMDZ are treated. This module is used 
! only if running without any other chemestry model as INCA or REPROBUS.  
!
  IMPLICIT NONE

  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: masktr   ! Masque reservoir de sol traceur
!$OMP THREADPRIVATE(masktr)                        ! Masque de l'echange avec la surface (1 = reservoir) ou (possible >= 1 )
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: fshtr    ! Flux surfacique dans le reservoir de sol
!$OMP THREADPRIVATE(fshtr)
  REAL,DIMENSION(:),ALLOCATABLE,SAVE   :: hsoltr   ! Epaisseur equivalente du reservoir de sol
!$OMP THREADPRIVATE(hsoltr)
!
!Radioelements:
!--------------
!
  REAL,DIMENSION(:),ALLOCATABLE,SAVE   :: tautr    ! Constante de decroissance radioactive
!$OMP THREADPRIVATE(tautr)
  REAL,DIMENSION(:),ALLOCATABLE,SAVE   :: vdeptr   ! Vitesse de depot sec dans la couche Brownienne
!$OMP THREADPRIVATE(vdeptr)
  REAL,DIMENSION(:),ALLOCATABLE,SAVE   :: scavtr   ! Coefficient de lessivage
!$OMP THREADPRIVATE(scavtr)
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: srcbe    ! Production du beryllium7 dans l atmosphere (U/s/kgA)
!$OMP THREADPRIVATE(srcbe)

  LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: radio    ! radio(it)   = true  => decroisssance radioactive
!$OMP THREADPRIVATE(radio)  

  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: trs     ! Conc. radon ds le sol
!$OMP THREADPRIVATE(trs)

  INTEGER,SAVE :: id_aga      ! Identification number for tracer : Age of stratospheric air
!$OMP THREADPRIVATE(id_aga)
  INTEGER,SAVE :: lev_1p5km   ! Approximative vertical layer number at 1.5km above surface, used for calculation of the age of air. The result shouldn't be that sensible to the exactness of this value as long as it is in the lower troposphere. 
!$OMP THREADPRIVATE(lev_1p5km)

  INTEGER,SAVE :: id_rn, id_pb ! Identification number for tracer : radon (Rn222), lead (Pb210)
!$OMP THREADPRIVATE(id_rn, id_pb)

  INTEGER,SAVE :: id_be       ! Activation et position du traceur Be7 [ id_be=0 -> desactive ] 
!$OMP THREADPRIVATE(id_be)

  INTEGER,SAVE :: id_pcsat, id_pcocsat, id_pcq ! traceurs pseudo-vapeur CL qsat, qsat_oc, q
!$OMP THREADPRIVATE(id_pcsat, id_pcocsat, id_pcq)
  INTEGER,SAVE :: id_pcs0, id_pcos0, id_pcq0   ! traceurs pseudo-vapeur CL qsat, qsat_oc, q
!                                              ! qui ne sont pas transportes par la convection
!$OMP THREADPRIVATE(id_pcs0, id_pcos0, id_pcq0)

  INTEGER, SAVE:: id_o3
!$OMP THREADPRIVATE(id_o3)
! index of ozone tracer with Cariolle parameterization
! 0 means no ozone tracer

  LOGICAL,SAVE :: rnpb=.FALSE. ! Presence du couple Rn222, Pb210
!$OMP THREADPRIVATE(rnpb)


CONTAINS

  SUBROUTINE traclmdz_from_restart(trs_in)
    ! This subroutine initialize the module saved variable trs with values from restart file (startphy.nc). 
    ! This subroutine is called from phyetat0 after the field trs_in has been read.
    
    USE dimphy
    USE infotrac_phy, ONLY: nbtr
    
    ! Input argument
    REAL,DIMENSION(klon,nbtr), INTENT(IN) :: trs_in 
    
    ! Local variables
    INTEGER :: ierr
    
    ! Allocate restart variables trs
    ALLOCATE( trs(klon,nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_from_restart', 'pb in allocation 1',1)
    
    ! Initialize trs with values read from restart file 
    trs(:,:) = trs_in(:,:)
    
  END SUBROUTINE traclmdz_from_restart


  SUBROUTINE traclmdz_init(pctsrf, xlat, xlon, ftsol, tr_seri, t_seri, pplay, sh, pdtphys, aerosol, lessivage)
    ! This subroutine allocates and initialize module variables and control variables.
    ! Initialization of the tracers should be done here only for those not found in the restart file.
    USE dimphy
    USE infotrac_phy, ONLY: nbtr, nqtot, tracers, pbl_flg, conv_flg
    USE regr_pr_comb_coefoz_m, ONLY: alloc_coefoz
    USE press_coefoz_m, ONLY: press_coefoz
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
    USE indice_sol_mod
    USE print_control_mod, ONLY: lunout

! Input variables
    REAL,DIMENSION(klon,nbsrf),INTENT(IN)     :: pctsrf ! Pourcentage de sol f(nature du sol)
    REAL,DIMENSION(klon),INTENT(IN)           :: xlat   ! latitudes en degres pour chaque point 
    REAL,DIMENSION(klon),INTENT(IN)           :: xlon   ! longitudes en degres pour chaque point 
    REAL,DIMENSION(klon,nbsrf),INTENT(IN)     :: ftsol  ! Temperature du sol (surf)(Kelvin)
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT) :: tr_seri! Concentration Traceur [U/KgA]  
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
    REAL,INTENT(IN)                        :: pdtphys ! Pas d'integration pour la physique (seconde)  

! Output variables
    LOGICAL,DIMENSION(nbtr), INTENT(OUT) :: aerosol
    LOGICAL,INTENT(OUT)                  :: lessivage
        
! Local variables    
    INTEGER :: ierr, it, iq, i, k
    REAL, DIMENSION(klon_glo,klev) :: varglo ! variable temporaire sur la grille global    
    REAL, DIMENSION(klev)          :: mintmp, maxtmp
    LOGICAL                        :: zero
! RomP >>> profil initial Be7
      integer ilesfil
      parameter (ilesfil=1)
      integer  irr,kradio
      real     beryllium(klon,klev)
! profil initial Pb210
      integer ilesfil2
      parameter (ilesfil2=1)
      integer  irr2,kradio2
      real     plomb(klon,klev)
!! RomP <<<
! --------------------------------------------
! Allocation
! --------------------------------------------
    ALLOCATE( scavtr(nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 9',1)
    scavtr(:)=1.
    
    ALLOCATE( radio(nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 11',1)
    radio(:) = .false.    ! Par defaut pas decroissance radioactive
    
    ALLOCATE( masktr(klon,nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 2',1)
    
    ALLOCATE( fshtr(klon,nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 3',1)
    
    ALLOCATE( hsoltr(nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 4',1)
    
    ALLOCATE( tautr(nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 5',1)
    tautr(:)  = 0.
    
    ALLOCATE( vdeptr(nbtr), stat=ierr)
    IF (ierr /= 0) CALL abort_physic('traclmdz_init', 'pb in allocation 6',1)
    vdeptr(:) = 0.


    lessivage  = .TRUE.
!!jyg(20130206) : le choix d activation du lessivage est fait dans phytrac avec iflag_lscav
!!    call getin('lessivage',lessivage)
!!    if(lessivage) then
!!     print*,'lessivage lsc ON'
!!    else
!!     print*,'lessivage lsc OFF'
!!    endif
    aerosol(:) = .FALSE.  ! Tous les traceurs sont des gaz par defaut
    
!
! Recherche des traceurs connus : Be7, O3, CO2,...
! --------------------------------------------
    id_rn=0; id_pb=0; id_aga=0; id_be=0; id_o3=0
    id_pcsat=0; id_pcocsat=0; id_pcq=0; id_pcs0=0; id_pcos0=0; id_pcq0=0
    it = 0
    DO iq = 1, nqtot
       IF(.NOT.(tracers(iq)%isInPhysics)) CYCLE
       it = it+1
       SELECT CASE(tracers(iq)%name)
         CASE("RN");                  id_rn     = it ! radon
         CASE("PB");                  id_pb     = it ! plomb
         CASE("Aga","AGA");           id_aga    = it ! Age of stratospheric air
         CASE("Be","BE","Be7","BE7"); id_be     = it ! Recherche du Beryllium 7
         CASE("o3","O3");             id_o3     = it ! Recherche de l'ozone
         CASE("pcsat",  "Pcsat");     id_pcsat  = it
         CASE("pcocsat","Pcocsat");   id_pcocsat= it
         CASE("pcq",    "Pcq");       id_pcq    = it
         CASE("pcs0",   "Pcs0");      id_pcs0   = it
         CASE("pcos0",  "Pcos0");     id_pcos0  = it
         CASE("pcq0",   "Pcq0");      id_pcq0   = it
         CASE DEFAULT
           WRITE(lunout,*) 'This is an unknown tracer in LMDZ : ', trim(tracers(iq)%name)
       END SELECT

       SELECT CASE(tracers(iq)%name)
         CASE("PB")                        !--- RomP >>> profil initial de PB210 
           OPEN(ilesfil2,file='prof.pb210',status='old',iostat=irr2)
           IF(irr2 == 0) THEN 
             READ(ilesfil2,*) kradio2
             WRITE(lunout,*)'number of levels for pb210 profile ',kradio2
             DO k=kradio2,1,-1; READ (ilesfil2,*) plomb(:,k); END DO
             CLOSE(ilesfil2)
             tr_seri(:,:,id_pb) = plomb(:,:)
           ELSE
             WRITE(lunout,*)'prof.pb210 does not exist: use restart values'
           END IF
         CASE("Aga","AGA")
           radio(id_aga) = .FALSE.
           aerosol(id_aga) = .FALSE.
           pbl_flg(id_aga) = 0 
           ! Find the first model layer above 1.5km from the surface
           IF (klev>=30) THEN
              lev_1p5km=6                  !--- NB: This value is for klev=39
           ELSE IF (klev>=10) THEN
              lev_1p5km=5                  !--- NB: This value is for klev=19
           ELSE
              lev_1p5km=klev/2
           END IF
         CASE("Be","BE","Be7","BE7")
           ALLOCATE( srcbe(klon,klev) )
           radio(id_be) = .TRUE.
           aerosol(id_be) = .TRUE.         !--- Le Be est un aerosol
           CALL init_be(pctsrf,pplay,masktr(:,id_be),tautr(id_be),vdeptr(id_be),scavtr(id_be),srcbe)
           WRITE(lunout,*) 'Initialisation srcBe: OK'
                                           !--- RomP >>> profil initial de Be7
           OPEN(ilesfil,file='prof.be7',status='old',iostat=irr)
           IF(irr == 0) THEN
             READ(ilesfil,*) kradio
             WRITE(lunout,*)'number of levels for Be7 profile ',kradio
             DO k=kradio,1,-1; READ(ilesfil,*) beryllium(:,k); END DO
             CLOSE(ilesfil)
             tr_seri(:,:,id_be)=beryllium(:,:)
           ELSE
             WRITE(lunout,*)'Prof. Be7 does not exist: use restart values'
           END IF
         CASE("o3","O3")                    !--- Parametrisation par la chimie de Cariolle
           CALL alloc_coefoz                !--- Allocate ozone coefficients
           CALL press_coefoz                !--- Read input pressure levels
         CASE("pcs0","Pcs0", "pcos0","Pcos0", "pcq0","Pcq0")
           conv_flg(it)=0                   !--- No transport by convection for this tracer
       END SELECT
    END DO

!
! Valeurs specifiques pour les traceurs Rn222 et Pb210
! ----------------------------------------------
    IF ( id_rn/=0 .AND. id_pb/=0 ) THEN
       rnpb = .TRUE.
       radio(id_rn)= .TRUE.
       radio(id_pb)= .TRUE.
       pbl_flg(id_rn) = 0 ! au lieu de clsol=true ! CL au sol calcule
       pbl_flg(id_pb) = 0 ! au lieu de clsol=true
       aerosol(id_rn) = .FALSE.
       aerosol(id_pb) = .TRUE. ! le Pb est un aerosol
       
       CALL initrrnpb (ftsol,pctsrf,masktr,fshtr,hsoltr,tautr,vdeptr,scavtr)
    END IF

!
! Check if all tracers have restart values
! ----------------------------------------------
    it = 0
    DO iq = 1, nqtot
       IF(.NOT.(tracers(iq)%isAdvected .AND. tracers(iq)%isInPhysics)) CYCLE
       it = it+1
       ! Test if tracer is zero everywhere. 
       ! Done by master process MPI and master thread OpenMP
       CALL gather(tr_seri(:,:,it),varglo)
       IF (is_mpi_root .AND. is_omp_root) THEN
          mintmp=MINVAL(varglo,dim=1)
          maxtmp=MAXVAL(varglo,dim=1)
          IF (MINVAL(mintmp,dim=1)==0. .AND. MAXVAL(maxtmp,dim=1)==0.) THEN
             ! Tracer is zero everywhere
             zero=.TRUE.
          ELSE
             zero=.FALSE.
          END IF
       END IF

       ! Distribute variable at all processes
       CALL bcast(zero)

       ! Initalize tracer that was not found in restart file.
       IF (zero) THEN
          ! The tracer was not found in restart file or it was equal zero everywhere.
          WRITE(lunout,*) "The tracer ",trim(tracers(iq)%name)," will be initialized"
          IF (it==id_pcsat .OR. it==id_pcq .OR. &
               it==id_pcs0 .OR. it==id_pcq0) THEN
             tr_seri(:,:,it) = 100.
          ELSE IF (it==id_pcocsat .OR. it==id_pcos0) THEN
             DO i = 1, klon
                IF ( pctsrf (i, is_oce) == 0. ) THEN
                   tr_seri(i,:,it) = 0.
                ELSE
                   tr_seri(i,:,it) = 100.
                END IF
             END DO
          ELSE
             ! No specific initialization exist for this tracer
             tr_seri(:,:,it) = 0.
          END IF
       END IF
    END DO

  END SUBROUTINE traclmdz_init

  SUBROUTINE traclmdz(nstep, julien, gmtime, pdtphys, t_seri, paprs, pplay, &
       cdragh,  coefh, yu1, yv1, ftsol, pctsrf, xlat, xlon, couchelimite, sh, &
       rh, pphi, ustar, wstar, ale_bl, ale_wake,  zu10m, zv10m, &
       tr_seri, source, d_tr_cl,d_tr_dec, zmasse)               !RomP
    
    USE dimphy
    USE infotrac_phy, ONLY: nbtr, pbl_flg
    USE strings_mod,  ONLY: int2str
    USE regr_pr_comb_coefoz_m, ONLY: regr_pr_comb_coefoz
    USE o3_chem_m, ONLY: o3_chem
    USE indice_sol_mod

    INCLUDE "YOMCST.h"

!==========================================================================
!                   -- DESCRIPTION DES ARGUMENTS --
!==========================================================================

! Input arguments
!
!Configuration grille,temps:
    INTEGER,INTENT(IN) :: nstep      ! nombre d'appels de la physiq
    INTEGER,INTENT(IN) :: julien     ! Jour julien
    REAL,INTENT(IN)    :: gmtime
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)  
    REAL,DIMENSION(klon),INTENT(IN) :: xlat    ! latitudes pour chaque point 
    REAL, INTENT(IN):: xlon(:) ! dim(klon) longitude

!
!Physique: 
!--------
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,intent(in):: zmasse (:, :)   ! dim(klon,klev) density of air, in kg/m2


!Couche limite:
!--------------
!
    REAL,DIMENSION(klon),INTENT(IN)      :: cdragh     ! coeff drag pour T et Q
    REAL,DIMENSION(klon,klev),INTENT(IN) :: coefh      ! coeff melange CL (m**2/s)
    REAL,DIMENSION(klon),INTENT(IN)      :: yu1        ! vents au premier niveau
    REAL,DIMENSION(klon),INTENT(IN)      :: yv1        ! vents au premier niveau
    LOGICAL,INTENT(IN)                   :: couchelimite
    REAL,DIMENSION(klon,klev),INTENT(IN) :: sh         ! humidite specifique
    REAL,DIMENSION(klon,klev),INTENT(IN) :: rh      ! Humidite relative
    REAL,DIMENSION(klon,klev),INTENT(IN) :: pphi    ! geopotentie
    REAL,DIMENSION(klon),INTENT(IN)      :: ustar   ! ustar (m/s)
    REAL,DIMENSION(klon),INTENT(IN)      :: wstar,ale_bl,ale_wake   ! wstar (m/s) and Avail. Lifti. Energ.
    REAL,DIMENSION(klon),INTENT(IN)      :: zu10m   ! vent zonal 10m (m/s)
    REAL,DIMENSION(klon),INTENT(IN)      :: zv10m   ! vent zonal 10m (m/s)

! Arguments necessaires pour les sources et puits de traceur:
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol  ! Temperature du sol (surf)(Kelvin)
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf ! Pourcentage de sol f(nature du sol)

! InOutput argument
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]  

! Output argument
    REAL,DIMENSION(klon,nbtr), INTENT(OUT)        :: source  ! a voir lorsque le flux de surface est prescrit 
    REAL,DIMENSION(klon,klev,nbtr), INTENT(OUT)   :: d_tr_cl ! Td couche limite/traceur

!=======================================================================================
!                        -- VARIABLES LOCALES TRACEURS --
!=======================================================================================

    INTEGER :: i, k, it
    INTEGER :: lmt_pas ! number of time steps of "physics" per day

    REAL,DIMENSION(klon)           :: d_trs    ! Td dans le reservoir
    REAL,DIMENSION(klon,klev)      :: qsat     ! pression de la vapeur a saturation
    REAL,DIMENSION(klon,klev,nbtr) :: d_tr_dec ! Td radioactive
    REAL                           :: zrho     ! Masse Volumique de l'air KgA/m3
    REAL                           :: amn, amx
!
!=================================================================
!  Ajout de la production en  Be7 (Beryllium) srcbe U/s/kgA
!=================================================================
!
    IF ( id_be /= 0 ) THEN
       DO k = 1, klev
          DO i = 1, klon
             tr_seri(i,k,id_be) = tr_seri(i,k,id_be)+srcbe(i,k)*pdtphys
          END DO
       END DO
       WRITE(*,*) 'Ajout srcBe dans tr_seri: OK'
    END IF
  

!=================================================================
! Update pseudo-vapor tracers 
!=================================================================

    CALL q_sat(klon*klev,t_seri,pplay,qsat)

    IF ( id_pcsat /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            tr_seri(i,k,id_pcsat) = qsat(i,k)
         ELSE
            tr_seri(i,k,id_pcsat) = MIN (qsat(i,k), tr_seri(i,k,id_pcsat))            
         END IF
      END DO 
     END DO
    END IF

    IF ( id_pcocsat /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            IF ( pctsrf (i, is_oce) > 0. ) THEN
               tr_seri(i,k,id_pcocsat) = qsat(i,k)
            ELSE
               tr_seri(i,k,id_pcocsat) = 0.
          END IF
       ELSE 
          tr_seri(i,k,id_pcocsat) = MIN (qsat(i,k), tr_seri(i,k,id_pcocsat))
       END IF
      END DO
     END DO
    END IF

    IF ( id_pcq /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            tr_seri(i,k,id_pcq) = sh(i,k)
         ELSE
            tr_seri(i,k,id_pcq) = MIN (qsat(i,k), tr_seri(i,k,id_pcq))    
         END IF
      END DO
     END DO
    END IF


    IF ( id_pcs0 /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            tr_seri(i,k,id_pcs0) = qsat(i,k)
         ELSE
            tr_seri(i,k,id_pcs0) = MIN (qsat(i,k), tr_seri(i,k,id_pcs0))    
         END IF
      END DO
     END DO
    END IF


    IF ( id_pcos0 /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            IF ( pctsrf (i, is_oce) > 0. ) THEN
               tr_seri(i,k,id_pcos0) = qsat(i,k)
            ELSE
               tr_seri(i,k,id_pcos0) = 0.
            END IF
         ELSE
            tr_seri(i,k,id_pcos0) = MIN (qsat(i,k), tr_seri(i,k,id_pcos0))
         END IF
      END DO
     END DO
    END IF


    IF ( id_pcq0 /= 0 ) THEN
     DO k = 1, klev
      DO i = 1, klon
         IF ( pplay(i,k).GE.85000.) THEN
            tr_seri(i,k,id_pcq0) = sh(i,k)
         ELSE 
            tr_seri(i,k,id_pcq0) = MIN (qsat(i,k), tr_seri(i,k,id_pcq0))
         END IF
      END DO
     END DO
    END IF

!=================================================================
! Update tracer : Age of stratospheric air 
!=================================================================
    IF (id_aga/=0) THEN
       
       ! Bottom layers
       DO k = 1, lev_1p5km
          tr_seri(:,k,id_aga) = 0.0
       END DO
       
       ! Layers above 1.5km
       DO k = lev_1p5km+1,klev-1
          tr_seri(:,k,id_aga) = tr_seri(:,k,id_aga) + pdtphys
       END DO
       
       ! Top layer
       tr_seri(:,klev,id_aga) = tr_seri(:,klev-1,id_aga)
       
    END IF

!======================================================================
!     -- Calcul de l'effet de la couche limite --
!======================================================================

    IF (couchelimite) THEN             
       source(:,:) = 0.0

       IF (id_be /=0) THEN
          DO i=1, klon
             zrho = pplay(i,1)/t_seri(i,1)/RD
             source(i,id_be) = - vdeptr(id_be)*tr_seri(i,1,id_be)*zrho
          END DO
       END IF

    END IF
    
    DO it=1, nbtr
       IF (couchelimite .AND. pbl_flg(it) == 0 .AND. (it==id_rn .OR. it==id_pb)) THEN 
          ! couche limite avec quantite dans le sol calculee
          CALL cltracrn(it, pdtphys, yu1, yv1,     &
               cdragh, coefh,t_seri,ftsol,pctsrf,  &
               tr_seri(:,:,it),trs(:,it),          &
               paprs, pplay, zmasse * rg, &
               masktr(:,it),fshtr(:,it),hsoltr(it),&
               tautr(it),vdeptr(it),               &
               xlat,d_tr_cl(:,:,it),d_trs)
          
          DO k = 1, klev
             DO i = 1, klon
                tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_cl(i,k,it)
             END DO
          END DO
        
          ! Traceur dans le reservoir sol
          DO i = 1, klon
             trs(i,it) = trs(i,it) + d_trs(i)
          END DO
       END IF
    END DO
          

!======================================================================
!   Calcul de l'effet du puits radioactif
!======================================================================
    CALL radio_decay (radio,rnpb,pdtphys,tautr,tr_seri,d_tr_dec)

    DO it=1,nbtr
       IF(radio(it)) then     
          DO k = 1, klev
             DO i = 1, klon
                tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_dec(i,k,it)
             END DO
          END DO
          CALL minmaxqfi(tr_seri(:,:,it),0.,1.e33,'puits rn it='//TRIM(int2str(it)))
       END IF
    END DO

!======================================================================
!   Parameterization of ozone chemistry
!======================================================================

    IF (id_o3 /= 0) then
       lmt_pas = NINT(86400./pdtphys)
       IF (MOD(nstep - 1, lmt_pas) == 0) THEN
          ! Once per day, update the coefficients for ozone chemistry:
          CALL regr_pr_comb_coefoz(julien, xlat, paprs, pplay)
       END IF
       CALL o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, xlat, &
            xlon, tr_seri(:, :, id_o3))
    END IF

  END SUBROUTINE traclmdz


  SUBROUTINE traclmdz_to_restart(trs_out)
    ! This subroutine is called from phyredem.F where the module 
    ! variable trs is written to restart file (restartphy.nc)
    USE dimphy
    USE infotrac_phy, ONLY: nbtr
    
    REAL,DIMENSION(klon,nbtr), INTENT(OUT) :: trs_out
    INTEGER :: ierr

    IF ( ALLOCATED(trs) ) THEN
       trs_out(:,:) = trs(:,:)
    ELSE
       ! No previous allocate of trs. This is the case for create_etat0_limit.
       trs_out(:,:) = 0.0
    END IF
    
  END SUBROUTINE traclmdz_to_restart
  

END MODULE traclmdz_mod
