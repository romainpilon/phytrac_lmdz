!$Id $

SUBROUTINE lsc_scav(pdtime,it,iflag_lscav, aerosol,  &
                    oliq,flxr,flxs,rneb,beta_fisrt,  &
                    beta_v1,pplay,paprs,t,tr_seri,   & 
                    d_tr_insc,d_tr_bcscav,d_tr_evap,qPrls) 
  USE ioipsl
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE traclmdz_mod
  USE infotrac_phy,ONLY : nbtr
  USE iophy
  IMPLICIT NONE 
!=====================================================================
! Objet : depot humide (lessivage et evaporation) de traceurs
! Inspired by routines of Olivier Boucher (mars 1998)
! author R. Pilon 10 octobre 2012
! last modification 16/01/2013 (reformulation partie evaporation)
!=====================================================================

  include "chem.h"
  include "YOMCST.h"
  include "YOECUMF.h" 

! inputs
  REAL,INTENT(IN)                        :: pdtime ! time step (s)
  INTEGER,INTENT(IN)                     :: it     ! tracer number
  INTEGER,INTENT(IN)                     :: iflag_lscav ! LS scavenging param: 3=Reddy_Boucher2004, 4=3+RPilon.
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: flxr     ! flux precipitant de pluie
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: flxs     ! flux precipitant de neige
  REAL,INTENT(IN)                        :: oliq ! contenu en eau liquide dans le nuage (kg/kg)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: rneb
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay    ! pression
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs    ! pression
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: t        ! temperature
! tracers
  LOGICAL,DIMENSION(nbtr), INTENT(IN)         :: aerosol
  REAL,DIMENSION(klon,klev,nbtr),INTENT(IN)   :: tr_seri        ! q de traceur  
  REAL,DIMENSION(klon,klev),INTENT(IN)        :: beta_fisrt     ! taux de conversion de l'eau cond
  REAL,DIMENSION(klon,klev),INTENT(OUT)       :: beta_v1        ! -- (originale version)
  REAL,DIMENSION(klon)                        :: his_dh         ! tendance de traceur integre verticalement
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)  :: d_tr_insc      ! tendance du traceur 
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)  :: d_tr_bcscav  ! tendance de traceur
  REAL,DIMENSION(klon,klev,nbtr),INTENT(OUT)  :: d_tr_evap
  REAL,DIMENSION(klon,nbtr),INTENT(OUT)       :: qPrls      !jyg: concentration tra dans pluie LS a la surf.
  REAL :: dxin,dxev                              ! tendance temporaire de traceur incloud
  REAL,DIMENSION(klon,klev) :: dxbc       ! tendance temporaire de traceur bc

!  variables locales      
 LOGICAL,SAVE :: debut=.TRUE.
!$OMP THREADPRIVATE(debut)
!
  REAL,PARAMETER :: henry=1.4  ! constante de Henry en mol/l/atm ~1.4 for gases
  REAL           :: henry_t    !  constante de Henry a T t  (mol/l/atm)
  REAL,PARAMETER :: kk=2900.   ! coefficient de dependence en T (K)
  REAL :: f_a     !  rapport de la phase aqueuse a la phase gazeuse
  REAL :: beta    !  taux de conversion de l'eau en pluie

  INTEGER :: i, k
  REAL,DIMENSION(klon,klev)    :: scav  !  water liquid content / fraction aqueuse du constituant
  REAL,DIMENSION(klon,klev)    :: zrho
  REAL,DIMENSION(klon,klev)    :: zdz
  REAL,DIMENSION(klon,klev)    :: zmass ! layer mass

  REAL           :: frac_ev       ! cste pour la reevaporation : dropplet shrinking
!  frac_ev = frac_gas ou frac_aer
  REAL,PARAMETER :: frac_gas=1.0  ! cste pour la reevaporation pour les gaz
  REAL           :: frac_aer      ! cste pour la reevaporation pour les particules
  REAL,DIMENSION(klon,klev) :: deltaP     ! P(i+1)-P(i)
  REAL,DIMENSION(klon,klev) :: beta_ev    !  dP/P(i+1)

!  101.325  m3/l x Pa/atm
!  R        Pa.m3/mol/K
!   cste de dissolution pour le depot humide
  REAL,SAVE :: frac_fine_scav
  REAL,SAVE :: frac_coar_scav
!$OMP THREADPRIVATE(frac_fine_scav, frac_coar_scav)

! below-cloud scav variables
! aerosol : alpha_r=0.001, gas 0.001  (Pruppacher & Klett 1967)
  REAL,SAVE :: alpha_r  !  coefficient d'impaction pour la pluie
  REAL,SAVE :: alpha_s  !  coefficient d'impaction pour la neige  
  REAL,SAVE :: R_r      !  mean raindrop radius (m) 
  REAL,SAVE :: R_s      !  mean snow crystal radius (m)
!$OMP THREADPRIVATE(alpha_r, alpha_s, R_r, R_s)
  REAL           :: pr, ps, ice, water
!  REAL :: conserv
!
!
  IF (debut) THEN
!
      alpha_r=0.001        !  coefficient d'impaction pour la pluie
      alpha_s=0.01         !  coefficient d'impaction pour la neige  
      R_r=0.001            !  mean raindrop radius (m) 
      R_s=0.001            !  mean snow crystal radius (m)
      frac_fine_scav=0.7
      frac_coar_scav=0.7
!     Droplet size shrinks by evap
      frac_aer=0.5
      debut=.FALSE.
!
      OPEN(99,file='lsc_scav_param.data',status='old', &
                form='formatted',err=9999)
      READ(99,*,end=9998)  alpha_r
      READ(99,*,end=9998)  alpha_s
      READ(99,*,end=9998)  R_r
      READ(99,*,end=9998)  R_s
      READ(99,*,end=9998)  frac_fine_scav
      READ(99,*,end=9998)  frac_coar_scav
      READ(99,*,end=9998)  frac_aer
9998  CONTINUE
      CLOSE(99)
9999  CONTINUE

   print*,'alpha_r',alpha_r
   print*,'alpha_s',alpha_s
   print*,'R_r',R_r
   print*,'R_s',R_s
   print*,'frac_fine_scav',frac_fine_scav
   print*,'frac_coar_scav',frac_coar_scav
   print*,'frac_aer ev',frac_aer
!
  ENDIF !(debut)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialization
  dxin=0.
  dxev=0.
  beta_ev=0.

  DO i=1,klon
   his_dh(i)=0.
  ENDDO

  DO k=1,klev
   DO i=1, klon 
    dxbc(i,k)=0.
    beta_v1(i,k)=0.
    deltaP(i,k)=0.
   ENDDO
  ENDDO

!  pressure and size of the layer
  DO k=klev, 1, -1
   DO i=1, klon 
     zrho(i,k)=pplay(i,k)/t(i,k)/RD   
     zdz(i,k)=(paprs(i,k)-paprs(i,k+1))/zrho(i,k)/RG
     zmass(i,k)=(paprs(i,k)-paprs(i,k+1))/RG
   ENDDO
  ENDDO

!jyg<
!! Temporary correction: all non-aerosol tracers are dealt with in the same way.
!! Should be updated once it has been decided how gases should be dealt with.
  IF (aerosol(it)) THEN
      frac_ev=frac_aer
  ELSE                                                !  gas
      frac_ev=frac_gas
  ENDIF
!
!jyg<
  IF (aerosol(it)) THEN ! aerosol
    DO k=1, klev
      DO i=1, klon 
       scav(i,k)=frac_fine_scav
      ENDDO 
    ENDDO
  ELSE                  ! gas
    DO k=1, klev
      DO i=1, klon 
       henry_t=henry*exp(-kk*(1./298.-1./t(i,k)))    !  mol/l/atm
       f_a=henry_t/101.325*R*t(i,k)*oliq*zrho(i,k)/rho_water
       scav(i,k)=f_a/(1.+f_a)
      ENDDO 
    ENDDO
  ENDIF

  DO k=klev-1, 1, -1
    DO i=1, klon 
!  incloud scavenging
      IF (iflag_lscav .EQ. 4) THEN
        beta=beta_fisrt(i,k)*rneb(i,k)
      ELSE
        beta=flxr(i,k)-flxr(i,k+1)+flxs(i,k)-flxs(i,k+1)
        beta=beta/zmass(i,k)/oliq
        beta=MAX(0.,beta)
      ENDIF ! (iflag_lscav .eq. 4)
      beta_v1(i,k)=beta    !! for output
!
      dxin=tr_seri(i,k,it)*(exp(-scav(i,k)*beta*pdtime)-1.)
      his_dh(i)=his_dh(i)-dxin*zmass(i,k)/pdtime !  kg/m2/s
      d_tr_insc(i,k,it)=dxin                     !  kg/kg/timestep

!  below-cloud impaction
!jyg<
      IF (.NOT.aerosol(it)) THEN
        d_tr_bcscav(i,k,it)=0.
      ELSE
        pr=0.5*(flxr(i,k)+flxr(i,k+1))
        ps=0.5*(flxs(i,k)+flxs(i,k+1))
        water=pr*alpha_r/R_r/rho_water
        ice=ps*alpha_s/R_s/rho_ice
        dxbc(i,k)=-3./4.*tr_seri(i,k,it)*pdtime*(water+ice)
!      add tracers from below cloud scav in his_dh
        his_dh(i)=his_dh(i)-dxbc(i,k)*zmass(i,k)/pdtime !  kg/m2/s
        d_tr_bcscav(i,k,it)=dxbc(i,k)                   !  kg/kg/timestep
      ENDIF

!  reevaporation
      deltaP(i,k)=flxr(i,k+1)+flxs(i,k+1)-flxr(i,k)-flxs(i,k)
      deltaP(i,k)=max(deltaP(i,k),0.)

      IF (flxr(i,k+1)+flxs(i,k+1).GT.1.e-16) THEn
        beta_ev(i,k)=deltaP(i,k)/(flxr(i,k+1)+flxs(i,k+1))
      ELSE
        beta_ev(i,k)=0.
      ENDIF

      beta_ev(i,k)=max(min(1.,beta_ev(i,k)),0.)

!jyg
      IF (ABS(1.-(1.-frac_ev)*beta_ev(i,k)).GT.1.e-16) THEN
! remove tracers from precipitation owing to release by evaporation in his_dh
        dxev=frac_ev*beta_ev(i,k)*his_dh(i)*pdtime/zmass(i,k)/(1.-(1.-frac_ev)*beta_ev(i,k))
        his_dh(i)=his_dh(i)*(1.-frac_ev*beta_ev(i,k)/(1.-(1.-frac_ev)*beta_ev(i,k)))
      ELSE
        dxev=his_dh(i)*pdtime/zmass(i,k)
        his_dh(i)=0.
      ENDIF
!
!      print*,  k, 'beta_ev',beta_ev
! remove tracers from precipitation owing to release by evaporation in his_dh
!      dxev=frac_ev*deltaP(i,k)*pdtime * his_dh(i) /(zrho(i,k)*zdz(i,k))
!rplmd
!      dxev=frac_ev*deltaP(i,k)*his_dh(i) *pdtime/(zrho(i,k)*zdz(i,k))/max(flxr(i,k)+flxs(i,k),1.e-16)

      d_tr_evap(i,k,it)=dxev       !  kg/kg/timestep
!
    ENDDO
  ENDDO
!
  DO i = 1,klon
     qPrls(i,it) = his_dh(i)/max(flxr(i,1)+flxs(i,1),1.e-16)
  ENDDO
!
! test de conservation
!      conserv=0.
!      DO k= klev,1,-1
!        DO i=1, klon
!         conserv=conserv+d_tr_insc(i,k,it)*(paprs(i,k)-paprs(i,k+1))/RG &
!                +d_tr_bcscav(i,k,it)*(paprs(i,k)-paprs(i,k+1))/RG  &
!                +d_tr_evap(i,k,it)*(paprs(i,k)-paprs(i,k+1))/RG
!      if(it.eq.3) write(*,'(I2,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12,2X,a,e20.12)'),&
!      k,'lsc conserv ',conserv,'insc',d_tr_insc(i,k,it),'bc',d_tr_bcscav(i,k,it),'ev',d_tr_evap(i,k,it)
!       ENDDO
!     ENDDO

END SUBROUTINE lsc_scav
