SUBROUTINE cv3a_uncompress(nloc, len, ncum, nd, ntra, idcum, compress, &
                           iflag, kbas, ktop, &
                           precip, cbmf, plcl, plfc, wbeff, sig, w0, ptop2, &
                           ft, fq, fqcomp, fu, fv, ftra,  &
                           sigd, ma, mip, vprecip, vprecipi, upwd, dnwd, dnwd0, &
                           qcondc, wd, cape, cin, &
                           tvp, &
                           ftd, fqd, &
                           plim1, plim2, asupmax, supmax0, &
                           asupmaxmin, &
                           da, phi, mp, phi2, d1a, dam, sigij, &                ! RomP+AC+jyg
                           qta, clw, elij, evap, ep, epmlmMm, eplaMm, &         ! RomP+jyg
                           wdtrainA, wdtrainS, wdtrainM, &                      ! RomP
                           qtc, sigt, detrain,         &
                           epmax_diag, & ! epmax_cape
                           iflag1, kbas1, ktop1, &
                           precip1, cbmf1, plcl1, plfc1, wbeff1, sig1, w01, ptop21, &
                           ft1, fq1, fqcomp1, fu1, fv1, ftra1, &
                           sigd1, ma1, mip1, vprecip1, vprecipi1, upwd1, dnwd1, dnwd01, &
                           qcondc1, wd1, cape1, cin1, &
                           tvp1, &
                           ftd1, fqd1, &
                           plim11, plim21, asupmax1, supmax01, &
                           asupmaxmin1, &
                           da1, phi1, mp1, phi21, d1a1, dam1, sigij1, &         ! RomP+AC+jyg
                           qta1, clw1, elij1, evap1, ep1, epmlmMm1, eplaMm1, &  ! RomP+jyg
                           wdtrainA1, wdtrainS1, wdtrainM1, &                   ! RomP
                           qtc1, sigt1, detrain1, &
                           epmax_diag1) ! epmax_cape

  ! **************************************************************
  ! *
  ! CV3A_UNCOMPRESS                                             *
  ! *
  ! *
  ! written by   : Sandrine Bony-Lena , 17/05/2003, 11.22.15    *
  ! modified by  : Jean-Yves Grandpeix, 23/06/2003, 10.36.17    *
  ! **************************************************************

  IMPLICIT NONE

  include "cv3param.h"

  ! inputs:
  INTEGER, INTENT (IN)                               :: nloc, len, ncum, nd, ntra
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: idcum(nloc)
!jyg<
  LOGICAL, INTENT (IN)                               :: compress
!>jyg
  INTEGER, DIMENSION (nloc), INTENT (IN)             ::iflag, kbas, ktop
  REAL, DIMENSION (nloc), INTENT (IN)                :: precip, cbmf, plcl, plfc
  REAL, DIMENSION (nloc), INTENT (IN)                :: wbeff
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: sig, w0
  REAL, DIMENSION (nloc), INTENT (IN)                :: ptop2
  REAL, DIMENSION (nloc), INTENT (IN)                :: epmax_diag
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: ft, fq, fqcomp, fu, fv
  REAL, DIMENSION (nloc, nd, ntra), INTENT (IN)      :: ftra
  REAL, DIMENSION (nloc), INTENT (IN)                :: sigd
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: ma, mip
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: vprecip
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: vprecipi
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: upwd, dnwd, dnwd0
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: qcondc
  REAL, DIMENSION (nloc), INTENT (IN)                :: wd, cape, cin
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: tvp
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: ftd, fqd
  REAL, DIMENSION (nloc), INTENT (IN)                :: plim1, plim2
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: asupmax
  REAL, DIMENSION (nloc), INTENT (IN)                :: supmax0, asupmaxmin

  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: da
  REAL, DIMENSION (nloc, nd, nd), INTENT (IN)        :: phi                    !AC!
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: mp                     !RomP
  REAL, DIMENSION (nloc, nd, nd), INTENT (IN)        :: phi2                   !RomP
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: d1a, dam               !RomP
  REAL, DIMENSION (nloc, nd, nd), INTENT (IN)        :: sigij                  !RomP
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: qta                    !jyg
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: clw                    !RomP
  REAL, DIMENSION (nloc, nd, nd), INTENT (IN)        :: elij                   !RomP
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: evap, ep               !RomP
  REAL, DIMENSION (nloc, nd, nd), INTENT (IN)        :: epmlmMm                !RomP+jyg
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: eplamM                 !RomP+jyg
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: qtc, sigt, detrain              !RomP
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: wdtrainA, wdtrainS, wdtrainM     !RomP

  ! outputs:
  INTEGER, DIMENSION (len), INTENT (OUT)             :: iflag1, kbas1, ktop1
  REAL, DIMENSION (len), INTENT (OUT)                :: precip1, cbmf1, plcl1, plfc1
  REAL, DIMENSION (len), INTENT (OUT)                :: wbeff1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: sig1, w01
  REAL, DIMENSION (len), INTENT (OUT)                :: epmax_diag1 ! epmax_cape
  REAL, DIMENSION (len), INTENT (OUT)                :: ptop21
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ft1, fq1, fqcomp1, fu1, fv1
  REAL, DIMENSION (len, nd, ntra), INTENT (OUT)      :: ftra1
  REAL, DIMENSION (len), INTENT (OUT)                :: sigd1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ma1, mip1
  REAL, DIMENSION (len, nd+1), INTENT (OUT)          :: vprecip1
  REAL, DIMENSION (len, nd+1), INTENT (OUT)          :: vprecipi1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: upwd1, dnwd1, dnwd01
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qcondc1
  REAL, DIMENSION (len), INTENT (OUT)                :: wd1, cape1, cin1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: tvp1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ftd1, fqd1
  REAL, DIMENSION (len), INTENT (OUT)                :: plim11, plim21
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: asupmax1
  REAL, DIMENSION (len), INTENT (OUT)                :: supmax01, asupmaxmin1
                                                    
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: da1
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: phi1                   !AC!
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: mp1                    !RomP
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: phi21                  !RomP
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: d1a1, dam1 !RomP       !RomP
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: sigij1                 !RomP
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qta1                   !jyg
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: clw1                   !RomP
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: elij1                  !RomP
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: evap1, ep1             !RomP
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: epmlmMm1               !RomP+jyg
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: eplamM1                !RomP+jyg
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qtc1, sigt1, detrain1            !RomP
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: wdtrainA1, wdtrainS1, wdtrainM1   !RomP


  ! local variables:
  INTEGER i, k, j
  INTEGER jdcum
  ! c    integer k1,k2

!jyg<
  IF (compress) THEN
!>jyg
    DO i = 1, ncum
      sig1(idcum(i), nd) = sig(i, nd)
      ptop21(idcum(i)) = ptop2(i)
      sigd1(idcum(i)) = sigd(i)
      precip1(idcum(i)) = precip(i)
      cbmf1(idcum(i)) = cbmf(i)
      plcl1(idcum(i)) = plcl(i)
      plfc1(idcum(i)) = plfc(i)
      wbeff1(idcum(i)) = wbeff(i)
      iflag1(idcum(i)) = iflag(i)
      kbas1(idcum(i)) = kbas(i)
      ktop1(idcum(i)) = ktop(i)
      wd1(idcum(i)) = wd(i)
      cape1(idcum(i)) = cape(i)
      cin1(idcum(i)) = cin(i)
      plim11(idcum(i)) = plim1(i)
      plim21(idcum(i)) = plim2(i)
      supmax01(idcum(i)) = supmax0(i)
      asupmaxmin1(idcum(i)) = asupmaxmin(i)
      epmax_diag1(idcum(i)) = epmax_diag(i)
    END DO
   
    DO k = 1, nl
      DO i = 1, ncum
        sig1(idcum(i), k) = sig(i, k)
        w01(idcum(i), k) = w0(i, k)
        ft1(idcum(i), k) = ft(i, k)
        fq1(idcum(i), k) = fq(i, k)
        fqcomp1(idcum(i), k) = fqcomp(i, k)
        fu1(idcum(i), k) = fu(i, k)
        fv1(idcum(i), k) = fv(i, k)
        ma1(idcum(i), k) = ma(i, k)
        mip1(idcum(i), k) = mip(i, k)
        vprecip1(idcum(i), k) = vprecip(i, k)
        vprecipi1(idcum(i), k) = vprecipi(i, k)
        upwd1(idcum(i), k) = upwd(i, k)
        dnwd1(idcum(i), k) = dnwd(i, k)
        dnwd01(idcum(i), k) = dnwd0(i, k)
        qcondc1(idcum(i), k) = qcondc(i, k)
        tvp1(idcum(i), k) = tvp(i, k)
        ftd1(idcum(i), k) = ftd(i, k)
        fqd1(idcum(i), k) = fqd(i, k)
        asupmax1(idcum(i), k) = asupmax(i, k)
   
        da1(idcum(i), k) = da(i, k) !AC!
        mp1(idcum(i), k) = mp(i, k) !RomP
        d1a1(idcum(i), k) = d1a(i, k) !RomP
        dam1(idcum(i), k) = dam(i, k) !RomP
        qta1(idcum(i), k) = qta(i, k) !jyg
        clw1(idcum(i), k) = clw(i, k) !RomP
        evap1(idcum(i), k) = evap(i, k) !RomP
        ep1(idcum(i), k) = ep(i, k) !RomP
        eplamM1(idcum(i), k) = eplamM(i, k) !RomP+jyg
        wdtrainA1(idcum(i), k) = wdtrainA(i, k) !RomP
        wdtrainS1(idcum(i), k) = wdtrainS(i, k) !RomP
        wdtrainM1(idcum(i), k) = wdtrainM(i, k) !RomP
        qtc1(idcum(i), k) = qtc(i, k)
        sigt1(idcum(i), k) = sigt(i, k)
        detrain1(idcum(i), k) = detrain(i, k)
   
      END DO
    END DO

! Fluxes are defined on a staggered grid and extend up to nl+1
    DO i = 1, ncum
      ma1(idcum(i), nlp) = 0.
      vprecip1(idcum(i), nlp) = 0.
      vprecipi1(idcum(i), nlp) = 0.
      upwd1(idcum(i), nlp) = 0.
      dnwd1(idcum(i), nlp) = 0.
      dnwd01(idcum(i), nlp) = 0.
    END DO
   
    ! AC!        do 2100 j=1,ntra
    ! AC!c oct3         do 2110 k=1,nl
    ! AC!         do 2110 k=1,nd ! oct3
    ! AC!          do 2120 i=1,ncum
    ! AC!            ftra1(idcum(i),k,j)=ftra(i,k,j)
    ! AC! 2120     continue
    ! AC! 2110    continue
    ! AC! 2100   continue
   
    ! AC!
!jyg<
!  Essais pour gagner du temps en diminuant l'adressage indirect 
!!    DO j = 1, nd
!!      DO k = 1, nd
!!        DO i = 1, ncum
!!          phi1(idcum(i), k, j) = phi(i, k, j) !AC!
!!          phi21(idcum(i), k, j) = phi2(i, k, j) !RomP
!!          sigij1(idcum(i), k, j) = sigij(i, k, j) !RomP
!!          elij1(idcum(i), k, j) = elij(i, k, j) !RomP
!!          epmlmMm(idcum(i), k, j) = epmlmMm(i, k, j) !RomP+jyg
!!        END DO
!!      END DO
!!    END DO

!!      DO i = 1, ncum
!!        jdcum=idcum(i)
!!        phi1    (jdcum, 1:nl+1, 1:nl+1) = phi    (i, 1:nl+1, 1:nl+1)          !AC!
!!        phi21   (jdcum, 1:nl+1, 1:nl+1) = phi2   (i, 1:nl+1, 1:nl+1)          !RomP
!!        sigij1  (jdcum, 1:nl+1, 1:nl+1) = sigij  (i, 1:nl+1, 1:nl+1)          !RomP
!!        elij1   (jdcum, 1:nl+1, 1:nl+1) = elij   (i, 1:nl+1, 1:nl+1)          !RomP
!!        epmlmMm1(jdcum, 1:nl+1, 1:nl+1) = epmlmMm(i, 1:nl+1, 1:nl+1)          !RomP+jyg
!!      END DO
!  These tracer associated arrays are defined up to nl, not nl+1
  DO i = 1, ncum
    jdcum=idcum(i)
    DO k = 1,nl
      DO j = 1,nl
        phi1    (jdcum, j, k) = phi    (i, j, k)          !AC!
        phi21   (jdcum, j, k) = phi2   (i, j, k)          !RomP
        sigij1  (jdcum, j, k) = sigij  (i, j, k)          !RomP
        elij1   (jdcum, j, k) = elij   (i, j, k)          !RomP
        epmlmMm1(jdcum, j, k) = epmlmMm(i, j, k)          !RomP+jyg
      END DO
    ENDDO
  ENDDO
!>jyg
    ! AC!
   
   
    ! do 2220 k2=1,nd
    ! do 2210 k1=1,nd
    ! do 2200 i=1,ncum
    ! ment1(idcum(i),k1,k2) = ment(i,k1,k2)
    ! sigij1(idcum(i),k1,k2) = sigij(i,k1,k2)
    ! 2200      enddo
    ! 2210     enddo
    ! 2220    enddo
!
!jyg<
  ELSE  !(compress)
!
      sig1(:,nd) = sig(:,nd)
      ptop21(:) = ptop2(:)
      sigd1(:) = sigd(:)
      precip1(:) = precip(:)
      cbmf1(:) = cbmf(:)
      plcl1(:) = plcl(:)
      plfc1(:) = plfc(:)
      wbeff1(:) = wbeff(:)
      iflag1(:) = iflag(:)
      kbas1(:) = kbas(:)
      ktop1(:) = ktop(:)
      wd1(:) = wd(:)
      cape1(:) = cape(:)
      cin1(:) = cin(:)
      plim11(:) = plim1(:)
      plim21(:) = plim2(:)
      supmax01(:) = supmax0(:)
      asupmaxmin1(:) = asupmaxmin(:)
!
      sig1(:, 1:nl) = sig(:, 1:nl)
      w01(:, 1:nl) = w0(:, 1:nl)
      ft1(:, 1:nl) = ft(:, 1:nl)
      fq1(:, 1:nl) = fq(:, 1:nl)
      fqcomp1(:, 1:nl) = fqcomp(:, 1:nl)
      fu1(:, 1:nl) = fu(:, 1:nl)
      fv1(:, 1:nl) = fv(:, 1:nl)
      ma1(:, 1:nl) = ma(:, 1:nl)
      mip1(:, 1:nl) = mip(:, 1:nl)
      vprecip1(:, 1:nl) = vprecip(:, 1:nl)
      vprecipi1(:, 1:nl) = vprecipi(:, 1:nl)
      upwd1(:, 1:nl) = upwd(:, 1:nl)
      dnwd1(:, 1:nl) = dnwd(:, 1:nl)
      dnwd01(:, 1:nl) = dnwd0(:, 1:nl)
      qcondc1(:, 1:nl) = qcondc(:, 1:nl)
      tvp1(:, 1:nl) = tvp(:, 1:nl)
      ftd1(:, 1:nl) = ftd(:, 1:nl)
      fqd1(:, 1:nl) = fqd(:, 1:nl)
      asupmax1(:, 1:nl) = asupmax(:, 1:nl)

      da1(:, 1:nl) = da(:, 1:nl)              !AC!
      mp1(:, 1:nl) = mp(:, 1:nl)              !RomP
      d1a1(:, 1:nl) = d1a(:, 1:nl)            !RomP
      dam1(:, 1:nl) = dam(:, 1:nl)            !RomP
      qta1(:, 1:nl) = qta(:, 1:nl)            !jyg
      clw1(:, 1:nl) = clw(:, 1:nl)            !RomP
      evap1(:, 1:nl) = evap(:, 1:nl)          !RomP
      ep1(:, 1:nl) = ep(:, 1:nl)              !RomP
      eplamM1(:, 1:nl) = eplamM(:, 1:nl)       !RomP+jyg
      wdtrainA1(:, 1:nl) = wdtrainA(:, 1:nl)  !RomP
      wdtrainS1(:, 1:nl) = wdtrainS(:, 1:nl)  !RomP
      wdtrainM1(:, 1:nl) = wdtrainM(:, 1:nl)  !RomP
      qtc1(:, 1:nl) = qtc(:, 1:nl)
      sigt1(:, 1:nl) = sigt(:, 1:nl)
      detrain1(:, 1:nl) = detrain(:, 1:nl)
!
      ma1(:, nlp) = 0.
      vprecip1(:, nlp) = 0.
      vprecipi1(:, nlp) = 0.
      upwd1(:, nlp) = 0.
      dnwd1(:, nlp) = 0.
      dnwd01(:, nlp) = 0.

!
      phi1    (:, 1:nl, 1:nl) = phi    (:, 1:nl, 1:nl)  !AC!
      phi21   (:, 1:nl, 1:nl) = phi2   (:, 1:nl, 1:nl)  !RomP
      sigij1  (:, 1:nl, 1:nl) = sigij  (:, 1:nl, 1:nl)  !RomP
      elij1   (:, 1:nl, 1:nl) = elij   (:, 1:nl, 1:nl)  !RomP
      epmlmMm1(:, 1:nl, 1:nl) = epmlmMm(:, 1:nl, 1:nl)  !RomP+jyg
  ENDIF !(compress)
!>jyg

  RETURN
END SUBROUTINE cv3a_uncompress

