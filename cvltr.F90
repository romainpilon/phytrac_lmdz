!
! $Id $
!
SUBROUTINE cvltr(pdtime, da, phi,phi2,d1a,dam, mpIN,epIN, &
!!           sigd,sij,clw,elij,epmlmMm,eplaMm,              &   !RL
           sigd,sij,wght_cvfd,clw,elij,epmlmMm,eplaMm,    &     !RL
           pmflxrIN,pmflxsIN,ev,te,wdtrainA,wdtrainM,     &
           paprs,it,tr,upd,dnd,inb,icb,                   &
           dtrcv,trsptd,dtrSscav,dtrsat,dtrUscav,qDi,qPr, &
           qPa,qMel,qTrdi,dtrcvMA,Mint,                   &
           zmfd1a,zmfphi2,zmfdam)
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
!  REAL,DIMENSION(klon,klev),INTENT(IN)    :: pplay ! pression aux 1/2 couches (bas en haut)
  REAL,DIMENSION(klon,klev,nbtr),INTENT(IN) :: tr     ! q de traceur (bas en haut)
  INTEGER,INTENT(IN)                        :: it 
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: upd   ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: dnd   ! saturated downdraft mass flux
!
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainA   ! masses precipitantes de l'asc adiab
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: wdtrainM   ! masses precipitantes des melanges
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pmflxrIN   ! vprecip: eau
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pmflxsIN   ! vprecip: neige
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
! variables pour les courants satures
  REAL,DIMENSION(klon,klev,klev)    :: zmd
  REAL,DIMENSION(klon,klev,klev)    :: za
  REAL,DIMENSION(klon,klev,nbtr)    :: zmfd,zmfa
  REAL,DIMENSION(klon,klev,nbtr)    :: zmfp,zmfu
  REAL,DIMENSION(klon,nbtr)         :: qfeed     ! tracer concentration feeding convection

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
! parametres lessivage
  REAL                    :: ccntrAA_coef        ! \alpha_a : fract aerosols de l'AA convertis en CCN
  REAL                    :: ccntrENV_coef       ! \beta_m  : fract aerosols de l'env convertis en CCN
  REAL                    :: coefcoli            ! coefficient de collision des gouttes sur les aerosols
!
  LOGICAL,DIMENSION(klon,klev) :: NO_precip
!  LOGICAL                      :: scavON
! var tmp tests
  REAL                           :: conserv
  real                           :: conservMA

! coefficient lessivage
   ccntrAA_coef     = 0.
   ccntrENV_coef    = 0.
   coefcoli         = 0.

!$OMP MASTER
  call getin('ccntrAA_coef',ccntrAA_coef)
  call getin('ccntrENV_coef',ccntrENV_coef)
  call getin('coefcoli',coefcoli)
!$OMP END MASTER
!$OMP BARRIER
  print*,'cvltr coef lessivage convectif', ccntrAA_coef,ccntrENV_coef,coefcoli

!  scavON=.TRUE.
!  if(scavON) then
!   ccntrAA_coef     = 1.
!   ccntrENV_coef    = 1.
!   coefcoli         = 1.
!  else
!   ccntrAA_coef     = 0.
!   ccntrENV_coef    = 0.
!   coefcoli         = 0.
!  endif

! ======================================================
! calcul de l'impaction
! ======================================================
!initialisation
  do j=1,klev
   do i=1,klon
     imp(i,j)=0.
   enddo
  enddo
! impaction sur la surface de la colonne de la descente insaturee
! On prend la moyenne des precip entre le niveau i+1 et i
! I=3/4* (P(1+1)+P(i))/2 / (sigd*r*rho_l)
!  1000kg/m3= densité de l'eau
! 0.75e-3 = 3/4 /1000
! Par la suite, I est tout le temps multiplié par sig_d pour avoir l'impaction sur la surface de la maille
! on le néglige ici pour simplifier le code
  do j=1,klev-1
   do i=1,klon
    imp(i,j) = coefcoli*0.75e-3/rdrop *&
             0.5*(pmflxr(i,j+1)+pmflxs(i,j+1)+pmflxr(i,j)+pmflxs(i,j))
!    rho(i,j)=pplay(i,j)/(rd*te(i,j))
   enddo
  enddo
!
! initialisation pour flux de traceurs, td et autre
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
    if(ev(i,j).lt.1.e-16) then
     evap(i,j)=0.
    else
     evap(i,j)=ev(i,j)*sigd(i)
    endif
   END DO
  END DO

  DO j=1,klev
   DO i=1,klon
   if(j.lt.klev) then
    if(epIN(i,j).lt.1.e-32) then
     ep(i,j)=0.
    else
     ep(i,j)=epIN(i,j)
    endif
   else
    ep(i,j)=epmax
   endif
    if(mpIN(i,j).lt.1.e-32) then
     mp(i,j)=0.
    else
     mp(i,j)=mpIN(i,j)
    endif
    if(pmflxsIN(i,j).lt.1.e-32) then
     pmflxs(i,j)=0.
    else
     pmflxs(i,j)=pmflxsIN(i,j)
    endif
    if(pmflxrIN(i,j).lt.1.e-32) then
     pmflxr(i,j)=0.
    else
     pmflxr(i,j)=pmflxrIN(i,j)
    endif
    if(wdtrainA(i,j).lt.1.e-32) then
     Pa(i,j)=0.
    else
     Pa(i,j)=wdtrainA(i,j)
    endif
    if(wdtrainM(i,j).lt.1.e-32) then
     Pm(i,j)=0.
    else
     Pm(i,j)=wdtrainM(i,j)
    endif
   END DO
  END DO

!==========================================
  DO j = klev-1,1,-1
   DO i = 1,klon
     NO_precip(i,j) = (pmflxr(i,j+1)+pmflxs(i,j+1)).lt.1.e-10&
                    .and.Pa(i,j).lt.1.e-10.and.Pm(i,j).lt.1.e-10
   END DO
  END DO

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
           if(mp(i,j+1).gt.1.e-10) then
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1)) !det ~ mk(j)=mk(j+1)*mp(i,j)/mp(i,j+1)
           ENDif
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
        if(mp(i,j+1).gt.1.e-10) then
          qTrdi(i,j+1,it)=qTrdi(i,j+1,it)+(zmd(i,j+1,k)/mp(i,j+1))*tr(i,k,it)
        else
          qTrdi(i,j,it)=0.!tr(i,j,it)
        endif
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
!
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
        if(j.ge.icb(i).and.j.le.inb(i)) then
          if(clw(i,j).gt.1.e-16) then
           qPa(i,j,it)=ccntrAA_coef*tr(i,1,it)/clw(i,j)
          else
           qPa(i,j,it)=0.
          endif
        endif
     END DO
  END DO

! calcul de q_pm en 2 parties :
! 1) calcul de sa valeur pour un niveau z' donne
! 2) integration sur la verticale sur z'
     DO j=1,klev
      DO k=1,j-1
        DO i=1,klon
        if(k.ge.icb(i).and.k.le.inb(i).and.& 
           j.le.inb(i)) then
            if(elij(i,k,j).gt.1.e-16) then
             qMeltmp(i,j,it)=((1-ep(i,k))*ccntrAA_coef*tr(i,1,it)&
                         *(1.-sij(i,k,j))  +ccntrENV_coef&
                         *tr(i,k,it)*sij(i,k,j)) / elij(i,k,j)
            else
             qMeltmp(i,j,it)=0.
            endif
          qpmMint(i,j,it)=qpmMint(i,j,it) + qMeltmp(i,j,it)*epmlmMm(i,j,k)
          Mint(i,j)=Mint(i,j) + epmlmMm(i,j,k)
        endif ! end if dans nuage
        END DO
     END DO
  END DO

     DO j=1,klev
        DO i=1,klon
          if(Mint(i,j).gt.1.e-16) then
           qMel(i,j,it)=qpmMint(i,j,it)/Mint(i,j)
          else
           qMel(i,j,it)=0.
          endif
     END DO
  END DO

! calcul de q_d et q_p    traceurs de la descente precipitante 
   DO j=klev-1,1,-1
    DO i=1,klon
     if(mp(i,j+1).gt.mp(i,j).and.mp(i,j+1).gt.1.e-10) then  ! detrainement
       kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                (-mp(i,j+1)-imp(i,j)/RG*dxpres(i,j))&
                + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))
             
     elseif(mp(i,j).gt.mp(i,j+1).and.mp(i,j).gt.1.e-10) then! entrainement
       if(j.eq.1) then
        kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                 (-mp(i,2)-imp(i,j)/RG*dxpres(i,j))&
                 + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))
       else
        kappa(i,j)=((pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))*&
                 (-mp(i,j)-imp(i,j)/RG*dxpres(i,j))&
                 + (imp(i,j)/RG*dxpres(i,j))*(evap(i,j)/RG*dxpres(i,j)))
       endif
      else
        kappa(i,j)=1.
      endif
    ENDDO
   ENDDO

  DO j=klev-1,1,-1
   DO i=1,klon
    if (abs(kappa(i,j)).lt.1.e-25) then    !si denominateur nul (il peut y avoir des mp!=0)
     kappa(i,j)=1.
     if(j.eq.1) then
       qDi(i,j,it)=qDi(i,j+1,it) !orig tr(i,j,it)   ! mp(1)=0 donc tout vient de la couche supérieure
     elseif(mp(i,j+1).gt.mp(i,j).and.mp(i,j+1).gt.1.e-10) then
       qDi(i,j,it)=qDi(i,j+1,it)
     elseif(mp(i,j).gt.mp(i,j+1).and.mp(i,j).gt.1.e-10) then! entrainement
       qDi(i,j,it)=(-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it))/(-mp(i,j))
     else  ! si  mp (i)=0 et mp(j+1)=0
       qDi(i,j,it)=tr(i,j,it) ! orig 0.
     endif 

      if(NO_precip(i,j)) then
       qPr(i,j,it)=0.
      else
       qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                   Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                   +imp(i,j)/RG*dxpres(i,j)*qDi(i,j,it))/&
                  (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
      endif
    else   !     denominateur non nul
     kappa(i,j)=1./kappa(i,j)     
! calcul de qd et qp
!!jyg  (20130119) correction pour le sommet du nuage
!!     if(j.ge.inb(i)) then       !au-dessus du nuage, sommet inclu
     if(j.gt.inb(i)) then       !au-dessus du nuage
       qDi(i,j,it)=tr(i,j,it)   ! pas de descente => environnement = descente insaturee
       qPr(i,j,it)=0.

!  vvv premiere couche du modele ou mp(1)=0  ! det tout le temps  vvv
     elseif(j.eq.1) then 
      if(mp(i,2).gt.1.e-10) then !mp(2) non nul -> detrainement (car mp(1) = 0) !ent pas possible
       if(NO_precip(i,j)) then !pas de precip en (i)
        qDi(i,j,it)=qDi(i,j+1,it)
        qPr(i,j,it)=0.
       else
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
       endif

      else !mp(2) nul -> plus de descente insaturee -> pluie agit sur environnement
        qDi(i,j,it)=tr(i,j,it) ! orig 0.
        if(NO_precip(i,j)) then
         qPr(i,j,it)=0.
        else
         qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                   Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                   +imp(i,j)/RG*dxpres(i,j)*tr(i,j,it))/&
                  (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
        endif

      endif  !mp(2) nul ou non

! vvv  (j!=1.and.j.lt.inb(i))  en-dessous du sommet nuage   vvv
     else 
!------------------------------------------------------------- detrainement 
      if(mp(i,j+1).gt.mp(i,j).and.mp(i,j+1).gt.1.e-10) then !mp(i,j).gt.1.e-10) then
         if(NO_precip(i,j)) then 
          qDi(i,j,it)=qDi(i,j+1,it)
          qPr(i,j,it)=0.
         else
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
         endif !precip
!------------------------------------------------------------- entrainement 
      elseif(mp(i,j).gt.mp(i,j+1).and.mp(i,j).gt.1.e-10) then 
         if(NO_precip(i,j)) then
          qDi(i,j,it)=(-mp(i,j+1)*(qDi(i,j+1,it)-tr(i,j,it))-mp(i,j)*tr(i,j,it))/(-mp(i,j))
          qPr(i,j,it)=0.
         else
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
         endif !precip
!------------------------------------------------------------- endif ! ent/det
      else  !mp nul
        qDi(i,j,it)=tr(i,j,it) ! orig 0.
        if(NO_precip(i,j)) then
         qPr(i,j,it)=0.
        else
         qPr(i,j,it)=((pmflxr(i,j+1)+pmflxs(i,j+1))*qPr(i,j+1,it)+&
                   Pa(i,j)*qPa(i,j,it)+Pm(i,j)*qMel(i,j,it)&
                   +imp(i,j)/RG*dxpres(i,j)*tr(i,j,it))/&
                  (pmflxr(i,j+1)+pmflxs(i,j+1)+Pa(i,j)+Pm(i,j))
        endif
      endif ! mp nul ou non
     endif ! condition sur j
    endif ! kappa
   ENDDO
  ENDDO

!! print test descente insaturee
!  DO j=klev,1,-1
!   DO i=1,klon
!     if(it.eq.3) then
!   write(*,'(I2,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12)') j,&
!!         'zmfdam',zmfdam(i,j,it),'zmfpsi',zmfphi2(i,j,it),&
!          'zmfdam+zmfpsi',zmfdam(i,j,it)+zmfphi2(i,j,it),'qpmMint',qpmMint(i,j,it),&
!          'Pm',Pm(i,j),'Mint',Mint(i,j),&
!!          'zmfa',zmfa(i,j,it),'zmfp',zmfp(i,j,it),&
!        'zmfdam',zmfdam(i,j,it),'zmfpsi',zmfphi2(i,j,it),'zmfd1a',zmfd1a(i,j,it)
!!          'Pa',Pa(i,j),'eplaMm',eplaMm(i,j)
!!         'zmfd1a=da1*qa',zmfd1a(i,j,it),'Pa*qPa',wdtrainA(i,j)*qPa(i,j,it),'da1',d1a(i,j)
!     endif
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
     scavtrac=-ccntrAA_coef*zmfd1a(i,k,it)&
               -zmfphi2(i,k,it)*ccntrENV_coef&
               -zmfdam(i,k,it)*ccntrAA_coef
! lessivage courants insatures
   if(k.le.inb(i).and.k.gt.1) then   ! tendances dans le nuage
!------------------------------------------------------------- detrainement
      if(mp(i,k+1).gt.mp(i,k).and.mp(i,k+1).gt.1.e-10) then
        uscavtrac= (-mp(i,k)+mp(i,k+1))*(qDi(i,k,it)-tr(i,k,it))&
                   + mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))
!
!        if(it.eq.3) write(*,'(I2,1X,a,5X,e20.12,82X,a,e20.12)')k,' det incloud',&
!                    (-mp(i,k)+mp(i,k+1))*(qDi(i,k,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k)+&
!                    mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k),&
!                    'mp',mp(i,k)
!------------------------------------------------------------- entrainement 
      elseif(mp(i,k).gt.mp(i,k+1).and.mp(i,k).gt.1.e-10) then
         uscavtrac= mp(i,k)*(tr(i,k-1,it)-tr(i,k,it))
!
!         if(it.eq.3) write(*,'(I2,1X,a,5X,e20.12,82X,a,e20.12)')k,' ent incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
!=!------------------------------------------------------------- end ent/det 
      else !        mp(i,k+1)=0. et mp(i,k)=0.    pluie directement sur l environnement

        if(NO_precip(i,k)) then
          uscavtrac=0.
!         if(it.eq.3) write(*,'(I2,1X,a,e20.12,82X,a,e20.12)')k,' no P ent incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
        else
          uscavtrac=-imp(i,k)*tr(i,k,it)*dxpres(i,k)/RG+evap(i,k)*qPr(i,k,it)*dxpres(i,k)/RG
!         if(it.eq.3) write(*,'(I2,1X,a,3X,e20.12,82X,a,e20.12)')k,' P env incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
        endif
      endif  ! mp/det/ent
!------------------------------------------------------------- premiere couche 
   elseif(k.eq.1) then  !                                      mp(1)=0.
      if(mp(i,2).gt.1.e-10) then  !detrainement
         uscavtrac= (-0.+mp(i,2))*(qDi(i,k,it)-tr(i,k,it)) !&
!                   + mp(i,2)*(0.-tr(i,k,it))
!
!       if(it.eq.3) write(*,'(I2,1X,a,e20.12,84X,a,e20.12)')k,' 1 det',&
!                  (-0.+mp(i,2))*(qDi(i,k,it)-tr(i,k,it))*pdtimeRG/dxpres(i,k)+&
!                  mp(i,2)*(0.-tr(i,k,it))*pdtimeRG/dxpres(i,k),&
!                   'mp',mp(i,k)
      else   ! mp(2) = 0 = mp(1) pas de descente insaturee, rien ne se passe s'il ne pleut pas, sinon pluie->env
        if(NO_precip(i,1)) then
          uscavtrac=0.
        else
          uscavtrac=-imp(i,k)*tr(i,k,it)*dxpres(i,k)/RG+evap(i,k)*qPr(i,k,it)*dxpres(i,k)/RG
        endif
!         if(it.eq.3) write(*,'(I2,1X,a,2X,e20.12,82X,a,e20.12)')k,'1 P env incloud',uscavtrac*pdtimeRG/dxpres(i,k), 'mp',mp(i,k)
      endif

   else   ! k > INB  au-dessus du nuage
    uscavtrac=0.
   endif

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

END SUBROUTINE cvltr
