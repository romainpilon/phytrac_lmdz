SUBROUTINE concvl(iflag_clos, &
                  dtime, paprs, pplay, k_upper_cv, &
                  t, q, t_wake, q_wake, s_wake, u, v, tra, ntra, &
                  Ale, Alp, sig1, w01, &
                  d_t, d_q, d_qcomp, d_u, d_v, d_tra, &
                  rain, snow, kbas, ktop, sigd, &
                  cbmf, plcl, plfc, wbeff, convoccur, &
                  upwd, dnwd, dnwdbis, &
                  Ma, mip, Vprecip, &
                  cape, cin, tvp, Tconv, iflag, &
                  pbase, bbase, dtvpdt1, dtvpdq1, dplcldt, dplcldr, &
                  qcondc, wd, pmflxr, pmflxs, &
!RomP >>>
!!     .             da,phi,mp,dd_t,dd_q,lalim_conv,wght_th)
                  da, phi, mp, phii, d1a, dam, sij, qta, clw, elij, &! RomP
                  dd_t, dd_q, lalim_conv, wght_th, &                 ! RomP
                  evap, ep, epmlmMm, eplaMm, &                       ! RomP
                  wdtrainA, wdtrainS, wdtrainM, wght, qtc, sigt, detrain, &
                  tau_cld_cv, coefw_cld_cv, &                           ! RomP+RL, AJ
!RomP <<<
                  epmax_diag) ! epmax_cape
! **************************************************************
! *
! CONCVL                                                      *
! *
! *
! written by   : Sandrine Bony-Lena, 17/05/2003, 11.16.04    *
! modified by :                                               *
! **************************************************************


  USE dimphy
  USE infotrac_phy, ONLY: nbtr
  USE phys_local_var_mod, ONLY: omega
  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE
! ======================================================================
! Auteur(s): S. Bony-Lena (LMD/CNRS) date: ???
! Objet: schema de convection de Emanuel (1991) interface
! ======================================================================
! Arguments:
! dtime--input-R-pas d'integration (s)
! s-------input-R-la vAleur "s" pour chaque couche
! sigs----input-R-la vAleur "sigma" de chaque couche
! sig-----input-R-la vAleur de "sigma" pour chaque niveau
! psolpa--input-R-la pression au sol (en Pa)
! pskapa--input-R-exponentiel kappa de psolpa
! h-------input-R-enthAlpie potentielle (Cp*T/P**kappa)
! q-------input-R-vapeur d'eau (en kg/kg)

! work*: input et output: deux variables de travail,
! on peut les mettre a 0 au debut
! ALE--------input-R-energie disponible pour soulevement
! ALP--------input-R-puissance disponible pour soulevement

! d_h--------output-R-increment de l'enthAlpie potentielle (h)
! d_q--------output-R-increment de la vapeur d'eau
! rain-------output-R-la pluie (mm/s)
! snow-------output-R-la neige (mm/s)
! upwd-------output-R-saturated updraft mass flux (kg/m**2/s)
! dnwd-------output-R-saturated downdraft mass flux (kg/m**2/s)
! dnwd0------output-R-unsaturated downdraft mass flux (kg/m**2/s)
! Ma---------output-R-adiabatic ascent mass flux (kg/m2/s)
! mip--------output-R-mass flux shed by adiabatic ascent (kg/m2/s)
! Vprecip----output-R-vertical profile of total precipitation (kg/m2/s)
! Tconv------output-R-environment temperature seen by convective scheme (K)
! Cape-------output-R-CAPE (J/kg)
! Cin -------output-R-CIN  (J/kg)
! Tvp--------output-R-Temperature virtuelle d'une parcelle soulevee
! adiabatiquement a partir du niveau 1 (K)
! deltapb----output-R-distance entre LCL et base de la colonne (<0 ; Pa)
! Ice_flag---input-L-TRUE->prise en compte de la thermodynamique de la glace
! dd_t-------output-R-increment de la temperature du aux descentes precipitantes
! dd_q-------output-R-increment de la vapeur d'eau du aux desc precip
! lalim_conv-
! wght_th----
! evap-------output-R
! ep---------output-R
! epmlmMm----output-R
! eplaMm-----output-R
! wdtrainA---output-R
! wdtrainS---output-R
! wdtrainM---output-R
! wght-------output-R
! ======================================================================


  include "clesphys.h"

  INTEGER, INTENT(IN)                           :: iflag_clos
  REAL, INTENT(IN)                              :: dtime
  REAL, DIMENSION(klon,klev),   INTENT(IN)      :: pplay
  REAL, DIMENSION(klon,klev+1), INTENT(IN)      :: paprs
  INTEGER,                      INTENT(IN)      :: k_upper_cv
  REAL, DIMENSION(klon,klev),   INTENT(IN)      :: t, q, u, v
  REAL, DIMENSION(klon,klev),   INTENT(IN)      :: t_wake, q_wake
  REAL, DIMENSION(klon),        INTENT(IN)      :: s_wake
  REAL, DIMENSION(klon,klev, nbtr),INTENT(IN)   :: tra
  INTEGER,                      INTENT(IN)      :: ntra
  REAL, DIMENSION(klon),        INTENT(IN)      :: Ale, Alp
!CR:test: on passe lentr et alim_star des thermiques
  INTEGER, DIMENSION(klon),     INTENT(IN)      :: lalim_conv
  REAL, DIMENSION(klon,klev),   INTENT(IN)      :: wght_th

  REAL, DIMENSION(klon,klev),   INTENT(INOUT)   :: sig1, w01

  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: d_t, d_q, d_qcomp, d_u, d_v
  REAL, DIMENSION(klon,klev, nbtr),INTENT(OUT)  :: d_tra
  REAL, DIMENSION(klon),        INTENT(OUT)     :: rain, snow

  INTEGER, DIMENSION(klon),     INTENT(OUT)     :: kbas, ktop
  REAL, DIMENSION(klon),        INTENT(OUT)     :: sigd
  REAL, DIMENSION(klon),        INTENT(OUT)     :: cbmf, plcl, plfc, wbeff
  REAL, DIMENSION(klon),        INTENT(OUT)     :: convoccur
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: upwd, dnwd, dnwdbis

!!       REAL Ma(klon,klev), mip(klon,klev),Vprecip(klon,klev)                    !jyg
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: Ma, mip 
  REAL, DIMENSION(klon,klev+1), INTENT(OUT)     :: Vprecip                        !jyg
  REAL, DIMENSION(klon),        INTENT(OUT)     :: cape, cin
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: tvp
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: Tconv
  INTEGER, DIMENSION(klon),     INTENT(OUT)     :: iflag
  REAL, DIMENSION(klon),        INTENT(OUT)     :: pbase, bbase
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: dtvpdt1, dtvpdq1
  REAL, DIMENSION(klon),        INTENT(OUT)     :: dplcldt, dplcldr
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: qcondc
  REAL, DIMENSION(klon),        INTENT(OUT)     :: wd
  REAL, DIMENSION(klon,klev+1), INTENT(OUT)     :: pmflxr, pmflxs

  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: da, mp
  REAL, DIMENSION(klon,klev,klev),INTENT(OUT)   :: phi
! RomP >>>
  REAL, DIMENSION(klon,klev,klev),INTENT(OUT)   :: phii
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: d1a, dam
  REAL, DIMENSION(klon,klev,klev),INTENT(OUT)   :: sij, elij
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: qta
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: clw
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: dd_t, dd_q
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: evap, ep
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: eplaMm
  REAL, DIMENSION(klon,klev,klev), INTENT(OUT)  :: epmlmMm
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: wdtrainA, wdtrainS, wdtrainM
! RomP <<<
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: wght                       !RL
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: qtc
  REAL, DIMENSION(klon,klev),   INTENT(OUT)     :: sigt, detrain
  REAL,                         INTENT(OUT)     :: tau_cld_cv, coefw_cld_cv
  REAL, DIMENSION(klon),        INTENT(OUT)     :: epmax_diag                ! epmax_cape

!
!  Local
!  ----
  REAL, DIMENSION(klon,klev)                    :: em_p
  REAL, DIMENSION(klon,klev+1)                  :: em_ph
  REAL                                          :: em_sig1feed ! sigma at lower bound of feeding layer
  REAL                                          :: em_sig2feed ! sigma at upper bound of feeding layer
  REAL, DIMENSION(klev)                         :: em_wght ! weight density determining the feeding mixture
  REAL, DIMENSION(klon,klev+1)                  :: Vprecipi                       !jyg
!on enleve le save
! SAVE em_sig1feed,em_sig2feed,em_wght

  REAL, DIMENSION(klon)                         :: rflag
  REAL, DIMENSION(klon)                         :: plim1, plim2
  REAL, DIMENSION(klon)                         :: ptop2
  REAL, DIMENSION(klon,klev)                    :: asupmax
  REAL, DIMENSION(klon)                         :: supmax0, asupmaxmin
  REAL                                          :: zx_t, zdelta, zx_qs, zcor
!
!   INTEGER iflag_mix
!   SAVE iflag_mix
  INTEGER                                       :: noff, minorig
  INTEGER                                       :: i,j, k, itra
  REAL, DIMENSION(klon,klev)                    :: qs, qs_wake
!LF          SAVE cbmf
!IM/JYG      REAL, SAVE, ALLOCATABLE :: cbmf(:)
!!!$OMP THREADPRIVATE(cbmf)!
  REAL, DIMENSION(klon)                         :: cbmflast


! Variables supplementaires liees au bilan d'energie
! Real paire(klon)
!LF      Real ql(klon,klev)
! Save paire
!LF      Save ql
!LF      Real t1(klon,klev),q1(klon,klev)
!LF      Save t1,q1
! Data paire /1./
  REAL, SAVE, ALLOCATABLE :: ql(:, :), q1(:, :), t1(:, :)
!$OMP THREADPRIVATE(ql, q1, t1)

! Variables liees au bilan d'energie et d'enthAlpi
  REAL ztsol(klon)
  REAL        h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot, &
              h_qs_tot, qw_tot, ql_tot, qs_tot, ec_tot
  SAVE        h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot, &
              h_qs_tot, qw_tot, ql_tot, qs_tot, ec_tot
!$OMP THREADPRIVATE(h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot)
!$OMP THREADPRIVATE(h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot)
  REAL        d_h_vcol, d_h_dair, d_qt, d_qw, d_ql, d_qs, d_ec
  REAL        d_h_vcol_phy
  REAL        fs_bound, fq_bound
  SAVE        d_h_vcol_phy
!$OMP THREADPRIVATE(d_h_vcol_phy)
  REAL        zero_v(klon)
  CHARACTER *15 ztit
  INTEGER     ip_ebil ! PRINT level for energy conserv. diag.
  SAVE        ip_ebil
  DATA        ip_ebil/2/
!$OMP THREADPRIVATE(ip_ebil)
  INTEGER     if_ebil ! level for energy conserv. dignostics
  SAVE        if_ebil
  DATA        if_ebil/2/
!$OMP THREADPRIVATE(if_ebil)
!+jld ec_conser
  REAL d_t_ec(klon, klev) ! tendance du a la conersion Ec -> E thermique
  REAL zrcpd
!-jld ec_conser
!LF
  INTEGER nloc
  LOGICAL, SAVE            :: first = .TRUE.
!$OMP THREADPRIVATE(first)
  INTEGER, SAVE            :: itap, igout
!$OMP THREADPRIVATE(itap, igout)


  include "YOMCST.h"
  include "YOMCST2.h"
  include "YOETHF.h"
  include "FCTTRE.h"
!jyg<
  include "conema3.h"
!>jyg

  IF (first) THEN
! Allocate some variables LF 04/2008

!IM/JYG allocate(cbmf(klon))
    ALLOCATE (ql(klon,klev))
    ALLOCATE (t1(klon,klev))
    ALLOCATE (q1(klon,klev))
!
    convoccur(:) = 0.
!
    itap = 0
    igout = klon/2 + 1/klon
  END IF
! Incrementer le compteur de la physique
  itap = itap + 1

! Copy T into Tconv
  DO k = 1, klev
    DO i = 1, klon
      Tconv(i, k) = t(i, k)
    END DO
  END DO

  IF (if_ebil>=1) THEN
    DO i = 1, klon
      ztsol(i) = t(i, 1)
      zero_v(i) = 0.
      DO k = 1, klev
        ql(i, k) = 0.
      END DO
    END DO
  END IF

! ym
  snow(:) = 0

  IF (first) THEN
    first = .FALSE.

! ===========================================================================
! READ IN PARAMETERS FOR THE CLOSURE AND THE MIXING DISTRIBUTION
! ===========================================================================

    IF (iflag_con==3) THEN
!      CALL cv3_inicp()
      CALL cv3_inip()
    END IF

! ===========================================================================
! READ IN PARAMETERS FOR CONVECTIVE INHIBITION BY TROPOS. DRYNESS
! ===========================================================================

! c$$$         open (56,file='supcrit.data')
! c$$$         read (56,*) Supcrit1, Supcrit2
! c$$$         close (56)

    IF (prt_level>=10) WRITE (lunout, *) 'supcrit1, supcrit2', supcrit1, supcrit2

! ===========================================================================
! Initialisation pour les bilans d'eau et d'energie
! ===========================================================================
    IF (if_ebil>=1) d_h_vcol_phy = 0.

    DO i = 1, klon
      cbmf(i) = 0.
!!          plcl(i) = 0.
      sigd(i) = 0.
    END DO
  END IF !(first)

! Initialisation a chaque pas de temps
  plfc(:) = 0.
  wbeff(:) = 100.
  plcl(:) = 0.

  DO k = 1, klev + 1
    DO i = 1, klon
      em_ph(i, k) = paprs(i, k)/100.0
      pmflxr(i, k) = 0.
      pmflxs(i, k) = 0.
    END DO
  END DO

  DO k = 1, klev
    DO i = 1, klon
      em_p(i, k) = pplay(i, k)/100.0
    END DO
  END DO


! Feeding layer

  em_sig1feed = 1.
!jyg<
!  em_sig2feed = 0.97
  em_sig2feed = cvl_sig2feed
!>jyg
! em_sig2feed = 0.8
! Relative Weight densities
  DO k = 1, klev
    em_wght(k) = 1.
  END DO
!CRtest: couche alim des tehrmiques ponderee par a*
! DO i = 1, klon
! do k=1,lalim_conv(i)
! em_wght(k)=wght_th(i,k)
! print*,'em_wght=',em_wght(k),wght_th(i,k)
! end do
! END DO

  IF (iflag_con==4) THEN
    DO k = 1, klev
      DO i = 1, klon
        zx_t = t(i, k)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zx_qs = min(0.5, r2es*foeew(zx_t,zdelta)/em_p(i,k)/100.0)
        zcor = 1./(1.-retv*zx_qs)
        qs(i, k) = zx_qs*zcor
      END DO
      DO i = 1, klon
        zx_t = t_wake(i, k)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zx_qs = min(0.5, r2es*foeew(zx_t,zdelta)/em_p(i,k)/100.0)
        zcor = 1./(1.-retv*zx_qs)
        qs_wake(i, k) = zx_qs*zcor
      END DO
    END DO
  ELSE ! iflag_con=3 (modif de puristes qui fait la diffce pour la convergence numerique)
    DO k = 1, klev
      DO i = 1, klon
        zx_t = t(i, k)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zx_qs = r2es*foeew(zx_t, zdelta)/em_p(i, k)/100.0
        zx_qs = min(0.5, zx_qs)
        zcor = 1./(1.-retv*zx_qs)
        zx_qs = zx_qs*zcor
        qs(i, k) = zx_qs
      END DO
      DO i = 1, klon
        zx_t = t_wake(i, k)
        zdelta = max(0., sign(1.,rtt-zx_t))
        zx_qs = r2es*foeew(zx_t, zdelta)/em_p(i, k)/100.0
        zx_qs = min(0.5, zx_qs)
        zcor = 1./(1.-retv*zx_qs)
        zx_qs = zx_qs*zcor
        qs_wake(i, k) = zx_qs
      END DO
    END DO
  END IF ! iflag_con

! ------------------------------------------------------------------

! Main driver for convection:
!                   iflag_con=3 -> nvlle version de KE (JYG)
!                   iflag_con = 30  -> equivAlent to convect3
!                   iflag_con = 4  -> equivAlent to convect1/2


  IF (iflag_con==30) THEN

! print *, '-> cv_driver'      !jyg
    CALL cv_driver(klon, klev, klevp1, ntra, iflag_con, &
                   t, q, qs, u, v, tra, &
                   em_p, em_ph, iflag, &
                   d_t, d_q, d_u, d_v, d_tra, rain, &
                   Vprecip, cbmf, sig1, w01, & !jyg
                   kbas, ktop, &
                   dtime, Ma, upwd, dnwd, dnwdbis, qcondc, wd, cape, &
                   da, phi, mp, phii, d1a, dam, sij, clw, elij, &       !RomP
                   evap, ep, epmlmMm, eplaMm, &                         !RomP
                   wdtrainA, wdtrainM, &                                !RomP
                   epmax_diag) ! epmax_cape
!           print *, 'cv_driver ->'      !jyg

    DO i = 1, klon
      cbmf(i) = Ma(i, kbas(i))
    END DO

!RL
    wght(:, :) = 0.
    DO i = 1, klon
      wght(i, 1) = 1.
    END DO
!RL

  ELSE

!LF   necessary for gathered fields
    nloc = klon
    CALL cva_driver(klon, klev, klev+1, ntra, nloc, k_upper_cv, &
                    iflag_con, iflag_mix, iflag_ice_thermo, &
                    iflag_clos, ok_conserv_q, dtime, cvl_comp_threshold, &
                    t, q, qs, t_wake, q_wake, qs_wake, s_wake, u, v, tra, &
                    em_p, em_ph, &
                    Ale, Alp, omega, &
                    em_sig1feed, em_sig2feed, em_wght, &
                    iflag, d_t, d_q, d_qcomp, d_u, d_v, d_tra, rain, kbas, ktop, &
                    cbmf, plcl, plfc, wbeff, sig1, w01, ptop2, sigd, &
                    Ma, mip, Vprecip, Vprecipi, upwd, dnwd, dnwdbis, qcondc, wd, &
                    cape, cin, tvp, &
                    dd_t, dd_q, plim1, plim2, asupmax, supmax0, &
                    asupmaxmin, lalim_conv, &
!AC!+!RomP+jyg
!!                   da,phi,mp,phii,d1a,dam,sij,clw,elij, &               ! RomP
!!                   evap,ep,epmlmMm,eplaMm,                              ! RomP
                    da, phi, mp, phii, d1a, dam, sij, wght, &           ! RomP+RL
                    qta, clw, elij, evap, ep, epmlmMm, eplaMm, &        ! RomP+RL
                    wdtrainA, wdtrainS, wdtrainM, qtc, sigt, detrain, &
                    tau_cld_cv, coefw_cld_cv, &                         ! RomP,AJ
!AC!+!RomP+jyg
                    epmax_diag) ! epmax_cape
  END IF
! ------------------------------------------------------------------
  IF (prt_level>=10) WRITE (lunout, *) ' cva_driver -> cbmf,plcl,plfc,wbeff, d_t, d_q ', &
                                         cbmf(1), plcl(1), plfc(1), wbeff(1), d_t(1,1), d_q(1,1)

  DO i = 1, klon
    rain(i) = rain(i)/86400.
    rflag(i) = iflag(i)
  END DO

  DO k = 1, klev
    DO i = 1, klon
      d_t(i, k) = dtime*d_t(i, k)
      d_q(i, k) = dtime*d_q(i, k)
      d_u(i, k) = dtime*d_u(i, k)
      d_v(i, k) = dtime*d_v(i, k)
    END DO
  END DO

  IF (iflag_con==3) THEN
    DO i = 1,klon
      IF (wbeff(i) > 100. .OR. wbeff(i) == 0 .OR. iflag(i) > 3) THEN
        wbeff(i) = 0.
        convoccur(i) = 0.  
      ELSE
        convoccur(i) = 1.
      ENDIF
    ENDDO
  ENDIF

  IF (iflag_con==30) THEN
    DO itra = 1, ntra
      DO k = 1, klev
        DO i = 1, klon
!RL!            d_tra(i,k,itra) =dtime*d_tra(i,k,itra)
          d_tra(i, k, itra) = 0.
        END DO
      END DO
    END DO
  END IF

!!AC!
  IF (iflag_con==3) THEN
    DO itra = 1, ntra
      DO k = 1, klev
        DO i = 1, klon
!RL!            d_tra(i,k,itra) =dtime*d_tra(i,k,itra)
          d_tra(i, k, itra) = 0.
        END DO
      END DO
    END DO
  END IF
!!AC!

  DO k = 1, klev
    DO i = 1, klon
      t1(i, k) = t(i, k) + d_t(i, k)
      q1(i, k) = q(i, k) + d_q(i, k)
    END DO
  END DO
!                                                     !jyg
  IF (iflag_con == 30 .OR. iflag_ice_thermo ==0) THEN
! --Separation neige/pluie (pour diagnostics)         !jyg
    DO k = 1, klev                                    !jyg
      DO i = 1, klon                                  !jyg
        IF (t1(i,k)<rtt) THEN                         !jyg
          pmflxs(i, k) = Vprecip(i, k)                !jyg
        ELSE                                          !jyg
          pmflxr(i, k) = Vprecip(i, k)                !jyg
        END IF                                        !jyg
      END DO                                          !jyg
    END DO                                            !jyg
  ELSE
    DO k = 1, klev                                    !jyg
      DO i = 1, klon                                  !jyg
        pmflxs(i, k) = Vprecipi(i, k)                 !jyg
        pmflxr(i, k) = Vprecip(i, k)-Vprecipi(i, k)   !jyg
      END DO                                          !jyg
    END DO                                            !jyg
  ENDIF

! c      IF (if_ebil.ge.2) THEN
! c        ztit='after convect'
! c        CALL diagetpq(paire,ztit,ip_ebil,2,2,dtime
! c     e      , t1,q1,ql,qs,u,v,paprs,pplay
! c     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
! c         call diagphy(paire,ztit,ip_ebil
! c     e      , zero_v, zero_v, zero_v, zero_v, zero_v
! c     e      , zero_v, rain, zero_v, ztsol
! c     e      , d_h_vcol, d_qt, d_ec
! c     s      , fs_bound, fq_bound )
! c      END IF


! les traceurs ne sont pas mis dans cette version de convect4:
  IF (iflag_con==4) THEN
    DO itra = 1, ntra
      DO k = 1, klev
        DO i = 1, klon
          d_tra(i, k, itra) = 0.
        END DO
      END DO
    END DO
  END IF
! print*, 'concvl->: dd_t,dd_q ',dd_t(1,1),dd_q(1,1)

  DO k = 1, klev
    DO i = 1, klon
      dtvpdt1(i, k) = 0.
      dtvpdq1(i, k) = 0.
    END DO
  END DO
  DO i = 1, klon
    dplcldt(i) = 0.
    dplcldr(i) = 0.
  END DO

  IF (prt_level>=20) THEN
    DO k = 1, klev
! print*,'physiq apres_add_con i k it d_u d_v d_t d_q qdl0',igout, &
!         k,itap,d_u_con(igout,k) ,d_v_con(igout,k), d_t_con(igout,k), &
!         d_q_con(igout,k),dql0(igout,k)
! print*,'phys apres_add_con itap Ma cin ALE ALP wak t q undi t q', &
!         itap,Ma(igout,k),cin(igout),ALE(igout), ALP(igout), &
!         t_wake(igout,k),q_wake(igout,k),t_undi(igout,k),q_undi(igout,k)
! print*,'phy apres_add_con itap CON rain snow EMA wk1 wk2 Vpp mip', &
!         itap,rain_con(igout),snow_con(igout),ema_work1(igout,k), &
!         ema_work2(igout,k),Vprecip(igout,k), mip(igout,k)
! print*,'phy apres_add_con itap upwd dnwd dnwd0 cape tvp Tconv ', &
!         itap,upwd(igout,k),dnwd(igout,k),dnwd0(igout,k),cape(igout), &
!         tvp(igout,k),Tconv(igout,k)
! print*,'phy apres_add_con itap dtvpdt dtvdq dplcl dplcldr qcondc', &
!         itap,dtvpdt1(igout,k),dtvpdq1(igout,k),dplcldt(igout), &
!         dplcldr(igout),qcondc(igout,k)
! print*,'phy apres_add_con itap wd pmflxr Kpmflxr Kp1 Kpmflxs Kp1', &
!         itap,wd(igout),pmflxr(igout,k),pmflxr(igout,k+1),pmflxs(igout,k), &
!         pmflxs(igout,k+1)
! print*,'phy apres_add_con itap da phi mp ftd fqd lalim wgth', &
!         itap,da(igout,k),phi(igout,k,k),mp(igout,k),ftd(igout,k), &
!         fqd(igout,k),lalim_conv(igout),wght_th(igout,k)
    END DO
  END IF !(prt_level.EQ.20) THEN

  RETURN
END SUBROUTINE concvl

