!
! $Id $
!
SUBROUTINE cvltr_scav(pdtime, da, phi,phi2,d1a,dam, mpIN,epIN, &
     sigd,sij,wght_cvfd,clw,elij,epmlmMm,eplaMm,    &
     pmflxrIN,pmflxsIN,ev,te,wdtrainA,wdtrainM,     &
     paprs,it,tr,upd,dnd,inb,icb,                   &
     ccntrAA_3d,ccntrENV_3d,coefcoli_3d,            &
     dtrcv,trsptd,dtrSscav,dtrsat,dtrUscav,qDi,qPr, &
     qPa,qMel,qTrdi,dtrcvMA,Mint,                   &
     zmfd1a,zmfphi2,zmfdam)
  !
  USE IOIPSL 
  USE dimphy
  USE infotrac_phy, ONLY : nbtr
  IMPLICIT NONE 
  !=====================================================================
  ! Objet : convection des traceurs / KE
  ! Auteurs: M-A Filiberti and J-Y Grandpeix
  ! modifiee par R Pilon : lessivage des traceurs / KE
  !=====================================================================

  include "YOMCST.h"
  include "YOECUMF.h"
  include "conema3.h"
  include "chem.h"

  ! Entree
  REAL,INTENT(IN)                           :: pdtime
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: da
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi
  ! RomP
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: d1a,dam ! matrices pour simplifier
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi2    ! l'ecriture des tendances
  !
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: mpIN
  REAL,DIMENSION(klon,klev+1),INTENT(IN)    :: paprs  ! pression aux 1/2 couches (bas en haut)
  INTEGER,INTENT(IN)                        :: it     ! numero du traceur
  REAL,DIMENSION(klon,klev,nbtr),INTENT(IN) :: tr     ! q de traceur (bas en haut)
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: upd   ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: dnd   ! saturated downdraft mass flux
  !
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainA   ! masses precipitantes de l'asc adiab
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainM   ! masses precipitantes des melanges
  !JE  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pmflxrIN   ! vprecip: eau
  REAL,DIMENSION(klon,klev+1),INTENT(IN)      :: pmflxrIN   ! vprecip: eau
  !JE  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pmflxsIN   ! vprecip: neige
  REAL,DIMENSION(klon,klev+1),INTENT(IN)      :: pmflxsIN   ! vprecip: neige
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: ev         ! evaporation cv30_routine
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: epIN
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: te 
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: sij        ! fraction dair de lenv
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wght_cvfd  ! weights of the layers feeding convection
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: elij       ! contenu en eau condensée spécifique/conc deau condensée massique
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: epmlmMm    ! eau condensee precipitee dans mel masse dair sat
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: eplaMm    ! eau condensee precipitee dans aa masse dair sat

  REAL,DIMENSION(klon,klev),INTENT(IN)      :: clw        ! contenu en eau condensée dans lasc adiab
  REAL,DIMENSION(klon),INTENT(IN)           :: sigd
  INTEGER,DIMENSION(klon),INTENT(IN)        :: icb,inb
  !
  REAL,DIMENSION(klon,klev),INTENT(IN) :: ccntrAA_3d
  REAL,DIMENSION(klon,klev),INTENT(IN) :: ccntrENV_3d
  REAL,DIMENSION(klon,klev),INTENT(IN) :: coefcoli_3d
  !
  ! Sortie
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: dtrcv     ! tendance totale (bas en haut)
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: dtrcvMA ! M-A Filiberti
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: trsptd    ! tendance du transport 
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: dtrSscav  ! tendance du lessivage courant sat
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: dtrsat    ! tendance trsp+sat scav
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)     :: dtrUscav  ! tendance du lessivage courant unsat
  !
  ! Variables locales
  INTEGER                           :: i,j,k
  REAL,DIMENSION(klon,klev)         :: dxpres   ! difference de pression entre niveau (j+1) et (j)
  REAL                              :: pdtimeRG ! pas de temps * gravite
  REAL,DIMENSION(klon,nbtr)         :: qfeed     ! tracer concentration feeding convection
  ! variables pour les courants satures
  REAL,DIMENSION(klon,klev,klev)    :: zmd
  REAL,DIMENSION(klon,klev,klev)    :: za
  REAL,DIMENSION(klon,klev,nbtr)    :: zmfd,zmfa
  REAL,DIMENSION(klon,klev,nbtr)    :: zmfp,zmfu

  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: zmfd1a
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: zmfdam
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: zmfphi2

  ! RomP ! les variables sont nettoyees des valeurs aberrantes
  REAL,DIMENSION(klon,klev)         :: Pa, Pm  ! pluie AA et mélanges, var temporaire
  REAL,DIMENSION(klon,klev)         :: pmflxs,pmflxr ! pmflxrIN,pmflxsIN sans valeur aberante
  REAL,DIMENSION(klon,klev)         :: mp            ! flux de masse
  REAL,DIMENSION(klon,klev)         :: ep            ! fraction d'eau convertie en precipitation 
  REAL,DIMENSION(klon,klev)         :: evap          ! evaporation : variable temporaire 
  REAL,DIMENSION(klon,klev)         :: rho    !environmental density

  REAL,DIMENSION(klon,klev)         :: kappa ! denominateur du au calcul de la matrice 
  ! pour obtenir qd et qp
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: qTrdi ! traceurs descente air insature transport MA
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: qDi  ! traceurs descente insaturees
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: qPr  ! traceurs colonne precipitante
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: qPa  ! traceurs dans les precip issues lasc. adiab.
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)    :: qMel ! traceurs dans les precip issues des melanges
  REAL,DIMENSION(klon,klev,nbtr)                :: qMeltmp ! variable temporaire 
  REAL,DIMENSION(klon,klev,nbtr)                :: qpmMint
  REAL,DIMENSION(klon,klev),INTENT(OUT)         :: Mint
  ! tendances
  REAL                              :: tdcvMA           ! terme de transport de traceur (schema Marie Angele)
  REAL                              :: trsptrac         ! terme de transport de traceur par l'air
  REAL                              :: scavtrac         ! terme de lessivage courant sature
  REAL                              :: uscavtrac        ! terme de lessivage courant insature
  ! impaction
!!!       Correction apres discussion Romain P. / Olivier B. 
!!!  REAL,PARAMETER                    :: rdrop=2.5e-3     ! rayon des gouttes d'eau
  REAL,PARAMETER                    :: rdrop=1.e-3     ! rayon des gouttes d'eau 
!!!
  REAL,DIMENSION(klon,klev)         :: imp              ! coefficient d'impaction
  !
  LOGICAL,DIMENSION(klon,klev) :: NO_precip
  ! var tmp tests
  REAL                           :: conserv
  real                           :: conservMA

!jyg<
!!  ! ======================================================
!!  ! calcul de l'impaction
!!  ! ======================================================
!!
!!  ! impaction sur la surface de la colonne de la descente insaturee
!!  ! On prend la moyenne des precip entre le niveau i+1 et i
!!  ! I=3/4* (P(1+1)+P(i))/2 / (sigd*r*rho_l)
!!  ! 1000kg/m3= densite de l'eau
!!  ! 0.75e-3 = 3/4 /1000
!!  ! Par la suite, I est tout le temps multiplie par sig_d pour avoir l'impaction sur la surface de la maille
!!!!  ! on le neglige ici pour simplifier le code
!!
!!  DO j=1,klev-1
!!     DO i=1,klon
!!        imp(i,j) = coefcoli_3d(i,j)*0.75e-3/rdrop *&
!!             0.5*(pmflxr(i,j+1)+pmflxs(i,j+1)+pmflxr(i,j)+pmflxs(i,j))
!!     ENDDO
!!  ENDDO
!>jyg
  !
  ! initialisation pour flux de traceurs, td et autre
  !
  trsptrac = 0.
  scavtrac = 0.
  uscavtrac = 0.
  qfeed(:,it) = 0.              !RL
  DO j=1,klev
     DO i=1,klon
        zmfd(i,j,it)=0.
        zmfa(i,j,it)=0.
        zmfu(i,j,it)=0.
        zmfp(i,j,it)=0.
        zmfphi2(i,j,it)=0.
        zmfd1a(i,j,it)=0.
        zmfdam(i,j,it)=0.
        qDi(i,j,it)=0.
        qPr(i,j,it)=0.
        qPa(i,j,it)=0.
        qMel(i,j,it)=0.
        qMeltmp(i,j,it)=0.
        qTrdi(i,j,it)=0.
        kappa(i,j)=0.
        trsptd(i,j,it)=0.
        dtrsat(i,j,it)=0.
        dtrSscav(i,j,it)=0.
        dtrUscav(i,j,it)=0.
        dtrcv(i,j,it)=0.
        dtrcvMA(i,j,it)=0.
        evap(i,j)=0.
        dxpres(i,j)=0.
        qpmMint(i,j,it)=0.
        Mint(i,j)=0.
     END DO
  END DO

  ! suppression des valeurs très faibles (~1e-320)
  ! multiplication de levaporation pour lavoir par unite de temps
  ! et par unite de surface de la maille
  ! -> cv30_unsat : evap : masse evaporee/s/(m2 de la descente)
  DO j=1,klev
     DO i=1,klon
        IF(ev(i,j).lt.1.e-16) THEN
           evap(i,j)=0.
        ELSE
           evap(i,j)=ev(i,j)*sigd(i)
        ENDIF
     END DO
  END DO

  DO j=1,klev
     DO i=1,klon
        IF(j.LT.klev) THEN
           IF(epIN(i,j).LT.1.e-32) THEN
              ep(i,j)=0.
           ELSE
              ep(i,j)=epIN(i,j)
           ENDIF
        ELSE
           ep(i,j)=epmax
        ENDIF
        IF(mpIN(i,j).LT.1.e-32) THEN
           mp(i,j)=0.
        ELSE
           mp(i,j)=mpIN(i,j)
        ENDIF
        IF(pmflxsIN(i,j).LT.1.e-32) THEN
           pmflxs(i,j)=0.
        ELSE
           pmflxs(i,j)=pmflxsIN(i,j)
        ENDIF
        IF(pmflxrIN(i,j).LT.1.e-32) THEN
           pmflxr(i,j)=0.
        ELSE
           pmflxr(i,j)=pmflxrIN(i,j)
        ENDIF
        IF(wdtrainA(i,j).LT.1.e-32) THEN
           Pa(i,j)=0.
        ELSE
           Pa(i,j)=wdtrainA(i,j)
        ENDIF
        IF(wdtrainM(i,j).LT.1.e-32) THEN
           Pm(i,j)=0.
        ELSE
           Pm(i,j)=wdtrainM(i,j)
        ENDIF
     END DO
  END DO

  !==========================================
  DO j = klev-1,1,-1
     DO i = 1,klon
        NO_precip(i,j) = (pmflxr(i,j+1)+pmflxs(i,j+1)).LT.1.e-10&
             .AND.Pa(i,j).LT.1.e-10.AND.Pm(i,j).LT.1.e-10
     END DO
  END DO

!jyg<
  ! ======================================================
  ! calcul de l'impaction
  ! ======================================================

  ! impaction sur la surface de la colonne de la descente insaturee
  ! On prend la moyenne des precip entre le niveau i+1 et i
  ! I=3/4* (P(1+1)+P(i))/2 / (sigd*r*rho_l)
  ! 1000kg/m3= densite de l'eau
  ! 0.75e-3 = 3/4 /1000
  ! Par la suite, I est tout le temps multiplie par sig_d pour avoir l'impaction sur la surface de la maille
  ! on le neglige ici pour simplifier le code

  DO j=1,klev-1
     DO i=1,klon
        imp(i,j) = coefcoli_3d(i,j)*0.75e-3/rdrop *&
             0.5*(pmflxr(i,j+1)+pmflxs(i,j+1)+pmflxr(i,j)+pmflxs(i,j))
     ENDDO
  ENDDO
!>jyg
  ! =========================================
  ! calcul des tendances liees au downdraft
  ! =========================================
  !cdir collapse
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmd(i,j,k)=0.
           za (i,j,k)=0.
        END DO
     END DO
  END DO
  ! calcul de la matrice d echange
  ! matrice de distribution de la masse entrainee en k
  ! commmentaire RomP : mp > 0
  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))   ! ~ mk(k)
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
           IF(mp(i,j+1).GT.1.e-10) THEN
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1)) !det ~ mk(j)=mk(j+1)*mp(i,j)/mp(i,j+1)
           ENDIF
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!!!!! quantite  de traceur dans la descente d'air insaturee  :   4 juin 2012
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           IF(mp(i,j+1).GT.1.e-10) THEN
              qTrdi(i,j+1,it)=qTrdi(i,j+1,it)+(zmd(i,j+1,k)/mp(i,j+1))*tr(i,k,it)
           ELSE
              qTrdi(i,j,it)=0.!tr(i,j,it)
           ENDIF
        ENDDO
     ENDDO
  ENDDO
!!!!!
  !
  ! rajout du terme lie a l ascendance induite
  !
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
  !
  ! tendance courants insatures  ! sans lessivage ancien schema
  !
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd(i,j,it)=zmfd(i,j,it)+za(i,j,k)*(tr(i,k,it)-tr(i,j,it))
        END DO
     END DO
  END DO
  !
  ! =========================================
  ! calcul des tendances liees aux courants satures   j <-> z ; k <-> z'
  ! =========================================
  !RL
  !  Feeding concentrations
  DO j=1,klev
     DO i=1,klon
        qfeed(i,it)=qfeed(i,it)+wght_cvfd(i,j)*tr(i,j,it)
     END DO
  END DO
  !RL
  !
  DO j=1,klev
     DO i=1,klon
        !RL
        !!        zmfa(i,j,it)=da(i,j)*(tr(i,1,it)-tr(i,j,it))                     ! da
        zmfa(i,j,it)=da(i,j)*(qfeed(i,it)-tr(i,j,it))                     ! da
        !RL
     END DO
  END DO
  !
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp(i,j,it)=zmfp(i,j,it)+phi(i,j,k)*(tr(i,k,it)-tr(i,j,it))  ! phi
        END DO
     END DO
  END DO
  ! RomP ajout des matrices liees au lessivage
  DO j=1,klev
     DO i=1,klon
        zmfd1a(i,j,it)=d1a(i,j)*tr(i,1,it)   ! da1
        zmfdam(i,j,it)=dam(i,j)*tr(i,1,it)   ! dam
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfphi2(i,j,it)=zmfphi2(i,j,it)+phi2(i,j,k)*tr(i,k,it)  ! psi
        END DO
     END DO
  END DO
  DO j=1,klev-1
     DO i=1,klon
        zmfu(i,j,it)=max(0.,upd(i,j+1)+dnd(i,j+1))*(tr(i,j+1,it)-tr(i,j,it))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu(i,j,it)=zmfu(i,j,it)+min(0.,upd(i,j)+dnd(i,j))*(tr(i,j,it)-tr(i,j-1,it))
     END DO
  END DO
  ! ===================================================
  ! calcul des tendances liees aux courants insatures
  ! ===================================================
  !  pression  
  DO k=1, klev
     DO i=1, klon
        dxpres(i,k)=paprs(i,k)-paprs(i,k+1)
     ENDDO
  ENDDO
  pdtimeRG=pdtime*RG

  ! q_pa et q_pm    traceurs issues des courants satures se retrouvant dans les precipitations
  DO j=1,klev
     DO i=1,klon
        IF(j.GE.icb(i).AND.j.LE.inb(i)) THEN
           IF(clw(i,j).GT.1.e-16) THEN
              !JE           qPa(i,j,it)=ccntrAA_coef*tr(i,1,it)/clw(i,j)
              qPa(i,j,it)=ccntrAA_3d(i,j)*tr(i,1,it)/clw(i,j)
           ELSE
              qPa(i,j,it)=0.
           ENDIF
        ENDIF
     END DO
  END DO

  ! calcul de q_pm en 2 parties :
  ! 1) calcul de sa valeur pour un niveau z' donne
  ! 2) integration sur la verticale sur z'
  DO j=1,klev
     DO k=1,j-1
        DO i=1,klon
           IF(k.GE.icb(i).AND.k.LE.inb(i).AND.& 
                j.LE.inb(i)) THEN
              IF(elij(i,k,j).GT.1.e-16) THEN
                 !JE             qMeltmp(i,j,it)=((1-ep(i,k))*ccntrAA_coef*tr(i,1,it)&
                 !JE                         *(1.-sij(i,k,j))  +ccntrENV_coef&
                 !JE                         *tr(i,k,it)*sij(i,k,j)) / elij(i,k,j)
                 qMeltmp(i,j,it)=((1-ep(i,k))*ccntrAA_3d(i,k)*tr(i,1,it)&
                      *(1.-sij(i,k,j))  +ccntrENV_3d(i,k)&
                      *tr(i,k,it)*sij(i,k,j)) / elij(i,k,j)
              ELSE
                 qMeltmp(i,j,it)=0.
              ENDIF
              qpmMint(i,j,it)=qpmMint(i,j,it) + qMeltmp(i,j,it)*epmlmMm(i,j,k)
              Mint(i,j)=Mint(i,j) + epmlmMm(i,j,k)
           ENDIF ! end if dans nuage
        END DO
     END DO
  END DO

  DO j=1,klev
     DO i=1,klon
        IF(Mint(i,j).GT.1.e-16) THEN
           qMel(i,j,it)=qpmMint(i,j,it)/Mint(i,j)
        ELSE
           qMel(i,j,it)=0.
        ENDIF
     END DO
  END DO

  ! calcul de q_d et q_p    traceurs de la descente precipitante 
  DO j=klev-1,1,-1
     DO i=1,klon
        IF(mp(i,j+1).GT.mp(i,j).AND.mp(i,j+1).GT.1.e-10) THEN  ! detrainement
           kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                (-mp(i,j+1)-imp(i,j)/RG*dxpres(i,j))&
                + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))

        ELSEIF(mp(i,j).GT.mp(i,j+1).AND.mp(i,j).GT.1.e-10) THEN! entrainement
           IF(j.eq.1) THEN
              kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                   (-mp(i,2)-imp(i,j)/RG*dxpres(i,j))&
                   + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))
           ELSE
              kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                   (-mp(i,j)-imp(i,j)/RG*dxpres(i,j))&
                   + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))
           ENDIF
        ELSE
           kappa(i,j)=1.
        ENDIF
     ENDDO
  ENDDO

  DO j=klev-1,1,-1
     DO i=1,klon
        IF (abs(kappa(i,j)).LT.1.e-25) THEN    !si denominateur nul (il peut y avoir des mp!=0)
           kappa(i,j)=1.
           IF(j.eq.1) THEN
              qDi(i,j,it)=qDi(i,j+1,it) !orig tr(i,j,it)   ! mp(1)=0 donc tout vient de la couche supérieure
           ELSEIF(mp(i,j+1).GT.mp(i,j).AND.mp(i,j+1).GT.1.e-10) THEN
              qDi(i,j,it)=qDi(i,j+1,it)
           ELSEIF(mp(i,j).GT.mp(i,j+1).AND.mp(i,j).GT.1.e-10) THEN! entrainement
              qDi(i,j,it)=(-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it))/(-mp(i,j))
           ELSE  ! si  mp (i)=0 et mp(j+1)=0
              qDi(i,j,it)=tr(i,j,it) ! orig 0.
           ENDIF

           IF(NO_precip(i,j)) THEN
              qPr(i,j,it)=0.
           ELSE
              qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                   Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                   +imp(i,j)/RG*dxpres(i,j)*qDi(i,j,it))/&
                   (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
           ENDIF
        ELSE   !     denominateur non nul
           kappa(i,j)=1./kappa(i,j)     
           ! calcul de qd et qp
           !!jyg  (20130119) correction pour le sommet du nuage
           !!     if(j.GE.inb(i)) THEN       !au-dessus du nuage, sommet inclu
           if(j.GT.inb(i)) THEN       !au-dessus du nuage
              qDi(i,j,it)=tr(i,j,it)   ! pas de descente => environnement = descente insaturee
              qPr(i,j,it)=0.

              !  vvv premiere couche du modele ou mp(1)=0  ! det tout le temps  vvv
           ELSEIF(j.eq.1) THEN 
              if(mp(i,2).GT.1.e-10) THEN !mp(2) non nul -> detrainement (car mp(1) = 0) !ent pas possible
                 if(NO_precip(i,j)) THEN !pas de precip en (i)
                    qDi(i,j,it)=qDi(i,j+1,it)
                    qPr(i,j,it)=0.
                 ELSE
                    qDi(i,j,it)=kappa(i,j)*(&
                         (-evap(i,j)/RG*dxpres(i,j))*((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)) +&
                         (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                         (-mp(i,j+1)*qDi(i,j+1,it)))

                    qPr(i,j,it)=kappa(i,j)*(&
                         (-mp(i,j+1)-imp(i,j)/RG*dxpres(i,j))*&
                         ((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it))&
                         +(-mp(i,j+1)*qDi(i,j+1,it)) * (imp(i,j)/RG*dxpres(i,j)))
                 ENDIF

              ELSE !mp(2) nul -> plus de descente insaturee -> pluie agit sur environnement
                 qDi(i,j,it)=tr(i,j,it) ! orig 0.
                 if(NO_precip(i,j)) THEN
                    qPr(i,j,it)=0.
                 ELSE
                    qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                         +imp(i,j)/RG*dxpres(i,j)*tr(i,j,it))/&
                         (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
                 ENDIF

              ENDIF  !mp(2) nul ou non

              ! vvv  (j!=1.AND.j.LT.inb(i))  en-dessous du sommet nuage   vvv
           ELSE 
              !------------------------------------------------------------- detrainement 
              if(mp(i,j+1).GT.mp(i,j).AND.mp(i,j+1).GT.1.e-10) THEN !mp(i,j).GT.1.e-10) THEN
                 if(NO_precip(i,j)) THEN 
                    qDi(i,j,it)=qDi(i,j+1,it)
                    qPr(i,j,it)=0.
                 ELSE
                    qDi(i,j,it)=kappa(i,j)*(&
                         (-evap(i,j)/RG*dxpres(i,j))*((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)) +&
                         (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                         (-mp(i,j+1)*qDi(i,j+1,it)))
                    !
                    qPr(i,j,it)=kappa(i,j)*(&
                         (-mp(i,j+1)-imp(i,j)/RG*dxpres(i,j))*&
                         ((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it))&
                         +(-mp(i,j+1)*qDi(i,j+1,it)) * (imp(i,j)/RG*dxpres(i,j)))
                 ENDIF !precip
                 !------------------------------------------------------------- entrainement 
              ELSEIF(mp(i,j).GT.mp(i,j+1).AND.mp(i,j).GT.1.e-10) THEN 
                 if(NO_precip(i,j)) THEN
                    qDi(i,j,it)=(-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it))/(-mp(i,j))
                    qPr(i,j,it)=0.
                 ELSE
                    qDi(i,j,it)=kappa(i,j)*(&
                         (-evap(i,j)/RG*dxpres(i,j))*((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)) +&
                         (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                         (-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it)))
                    !
                    qPr(i,j,it)=kappa(i,j)*(&
                         (-mp(i,j)-imp(i,j)/RG*dxpres(i,j))*&
                         ((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it))&
                         +(-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it))*&
                         (imp(i,j)/RG*dxpres(i,j)))
                 ENDIF !precip
                 !------------------------------------------------------------- ENDIF ! ent/det
              ELSE  !mp nul
                 qDi(i,j,it)=tr(i,j,it) ! orig 0.
                 if(NO_precip(i,j)) THEN
                    qPr(i,j,it)=0.
                 ELSE
                    qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                         Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                         +imp(i,j)/RG*dxpres(i,j)*tr(i,j,it))/&
                         (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
                 ENDIF
              ENDIF ! mp nul ou non
           ENDIF ! condition sur j
        ENDIF ! kappa
     ENDDO
  ENDDO

  !! print test descente insaturee
  !  DO j=klev,1,-1
  !   DO i=1,klon
  !     if(it.eq.3) THEN
  !   write(*,'(I2,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12)') j,&
  !!         'zmfdam',zmfdam(i,j,it),'zmfpsi',zmfphi2(i,j,it),&
  !          'zmfdam+zmfpsi',zmfdam(i,j,it)+zmfphi2(i,j,it),'qpmMint',qpmMint(i,j,it),&
  !          'Pm',Pm(i,j),'Mint',Mint(i,j),&
  !!          'zmfa',zmfa(i,j,it),'zmfp',zmfp(i,j,it),&
  !        'zmfdam',zmfdam(i,j,it),'zmfpsi',zmfphi2(i,j,it),'zmfd1a',zmfd1a(i,j,it)
  !!          'Pa',Pa(i,j),'eplaMm',eplaMm(i,j)
  !!         'zmfd1a=da1*qa',zmfd1a(i,j,it),'Pa*qPa',wdtrainA(i,j)*qPa(i,j,it),'da1',d1a(i,j)
  !     ENDIF
  !   ENDDO
  !  ENDDO


  ! ===================================================
  ! calcul final des tendances
  ! ===================================================

  DO k=klev-1,1,-1
     DO i=1, klon
        ! transport
        tdcvMA=zmfd(i,k,it)+zmfu(i,k,it)+zmfa(i,k,it)+zmfp(i,k,it)   ! double comptage des downdraft insatures
        trsptrac=zmfu(i,k,it)+zmfa(i,k,it)+zmfp(i,k,it)
        ! lessivage courants satures
        !JE     scavtrac=-ccntrAA_coef*zmfd1a(i,k,it)&
        !JE               -zmfphi2(i,k,it)*ccntrENV_coef&
        !JE               -zmfdam(i,k,it)*ccntrAA_coef
        scavtrac=-ccntrAA_3d(i,k)*zmfd1a(i,k,it)&
             -zmfphi2(i,k,it)*ccntrENV_3d(i,k)&
             -zmfdam(i,k,it)*ccntrAA_3d(i,k)
        ! lessivage courants insatures
        if(k.LE.inb(i).AND.k.GT.1) THEN   ! tendances dans le nuage
           !------------------------------------------------------------- detrainement
           if(mp(i,k+1).GT.mp(i,k).AND.mp(i,k+1).GT.1.e-10) THEN
              uscavtrac= (-mp(i,k)+mp(i,k+1))*(qDi(i,k,it)-tr(i,k,it))&
                   + mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))
              !
              !        if(it.eq.3) write(*,'(I2,1X,a,5X,e20.12,82X,a,e20.12)')k,' det incloud',&
              !                    (-mp(i,k)+mp(i,k+1))*(qDi(i,k,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k)+&
              !                    mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k),&
              !                    'mp',mp(i,k)
              !------------------------------------------------------------- entrainement 
           ELSEIF(mp(i,k).GT.mp(i,k+1).AND.mp(i,k).GT.1.e-10) THEN
              uscavtrac= mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))
              !
              !         if(it.eq.3) write(*,'(I2,1X,a,5X,e20.12,82X,a,e20.12)')k,' ent incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
              !=!------------------------------------------------------------- end ent/det 
           ELSE !        mp(i,k+1)=0. et mp(i,k)=0.    pluie directement sur l environnement

              if(NO_precip(i,k)) THEN
                 uscavtrac=0.
                 !         if(it.eq.3) write(*,'(I2,1X,a,e20.12,82X,a,e20.12)')k,' no P ent incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
              ELSE
                 uscavtrac=-imp(i,k)*tr(i,k,it)*dxpres(i,k)/RG+evap(i,k)*qPr(i,k,it)*dxpres(i,k)/RG
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
                 !!JE adds
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac, 'imp',imp(i,k)
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac, 'tr',tr(i,k,it)
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac, 'evap',evap(i,k)
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac, 'qPr',qPr(i,k,it)
                 !         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac, 'dxpres',dxpres(i,k)
                 !!Je end

              ENDIF
           ENDIF  ! mp/det/ent
           !------------------------------------------------------------- premiere couche 
        ELSEIF(k.eq.1) THEN  !                                      mp(1)=0.
           if(mp(i,2).GT.1.e-10) THEN  !detrainement
              uscavtrac= (-0.+mp(i,2))*(qDi(i,k,it)-tr(i,k,it)) !&
              !                   + mp(i,2)*(0.-tr(i,k,it))
              !
              !       if(it.eq.3) write(*,'(I2,1X,a,e20.12,84X,a,e20.12)')k,' 1 det',&
              !                  (-0.+mp(i,2))*(qDi(i,k,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k)+&
              !                  mp(i,2)*(0.-tr(i,k,it))*pdtimeRG/dxpres(i,k),&
              !                   'mp',mp(i,k)
           ELSE   ! mp(2) = 0 = mp(1) pas de descente insaturee, rien ne se passe s'il ne pleut pas, sinon pluie->env
              if(NO_precip(i,1)) THEN
                 uscavtrac=0.
              ELSE
                 uscavtrac=-imp(i,k)*tr(i,k,it)*dxpres(i,k)/RG+evap(i,k)*qPr(i,k,it)*dxpres(i,k)/RG
              ENDIF
              !         if(it.eq.3) write(*,'(I2,1X,a,2X,e20.12,82X,a,e20.12)')k,'1 P env incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
           ENDIF

        ELSE   ! k > INB  au-dessus du nuage
           uscavtrac=0.
        ENDIF

        ! =====    tendances finales  ======
        trsptd(i,k,it)=trsptrac*pdtimeRG/dxpres(i,k)              ! td transport sans eau dans courants satures
        dtrSscav(i,k,it)=scavtrac*pdtimeRG/dxpres(i,k)            ! td du lessivage dans courants satures
        dtrUscav(i,k,it)=uscavtrac*pdtimeRG/dxpres(i,k)           ! td courant insat
        dtrsat(i,k,it)=(trsptrac+scavtrac)*pdtimeRG/dxpres(i,k)   ! td courant sat
        dtrcv(i,k,it)=(trsptrac+scavtrac+uscavtrac)*pdtimeRG/dxpres(i,k)!dtrsat(i,k,it)+dtrUscav(i,k,it)         td conv
!!!!!!
        dtrcvMA(i,k,it)=tdcvMA*pdtimeRG/dxpres(i,k) ! MA tendance convection
     ENDDO
  ENDDO

  ! test de conservation du traceur
  !print*,"_____________________________________________________________"
  !print*,"                                                             "
  !      conserv=0.
  !      conservMA=0.
  !      DO k= klev-1,1,-1
  !        DO i=1, klon
  !         conserv=conserv+dtrcv(i,k,it)*   &
  !        (paprs(i,k)-paprs(i,k+1))/RG
  !         conservMA=conservMA+dtrcvMA(i,k,it)*   &
  !        (paprs(i,k)-paprs(i,k+1))/RG
  !
  !      if(it.eq.3)  write(*,'(I2,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12)') k,&
  !         'MA td ',dtrcvMA(i,k,it)*dxpres(i,k)/RG,&
  !         ' td',dtrcv(i,k,it)*dxpres(i,k)/RG,'         conservMA ',conservMA,'conserv ',conserv
  !!
  !        ENDDO
  !      ENDDO
  !       if(it.eq.3) print *,'it',it,'conserv ',conserv,'conservMA ',conservMA

END SUBROUTINE cvltr_scav
