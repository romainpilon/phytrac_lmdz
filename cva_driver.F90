
! $Id: cva_driver.F90 4613 2023-07-07 02:01:37Z fhourdin $

SUBROUTINE cva_driver(len, nd, ndp1, ntra, nloc, k_upper, &
                      iflag_con, iflag_mix, iflag_ice_thermo, iflag_clos, ok_conserv_q, &
!!                      delt, t1, q1, qs1, t1_wake, q1_wake, qs1_wake, s1_wake, &  ! jyg
                      delt, comp_threshold, &                                      ! jyg
                      t1, q1, qs1, t1_wake, q1_wake, qs1_wake, s1_wake, &          ! jyg
                      u1, v1, tra1, &
                      p1, ph1, &
                      Ale1, Alp1, omega1, &
                      sig1feed1, sig2feed1, wght1, &
                      iflag1, ft1, fq1, fqcomp1, fu1, fv1, ftra1, &
                      precip1, kbas1, ktop1, &
                      cbmf1, plcl1, plfc1, wbeff1, &
                      sig1, w01, & !input/output
                      ptop21, sigd1, &
                      ma1, mip1, Vprecip1, Vprecipi1, upwd1, dnwd1, dnwd01, &      ! jyg
                      qcondc1, wd1, &
                      cape1, cin1, tvp1, &
                      ftd1, fqd1, &
                      Plim11, Plim21, asupmax1, supmax01, asupmaxmin1, &
                      lalim_conv1, & 
!!                      da1,phi1,mp1,phi21,d1a1,dam1,sigij1,clw1, &        ! RomP
!!                      elij1,evap1,ep1,epmlmMm1,eplaMm1, &                ! RomP
                      da1, phi1, mp1, phi21, d1a1, dam1, sigij1, wghti1, & ! RomP, RL
                      qta1, clw1, elij1, evap1, ep1, epmlmMm1, eplaMm1, &  ! RomP, RL
                      wdtrainA1, wdtrainS1, wdtrainM1, qtc1, sigt1, detrain1, tau_cld_cv, &     !!jygprl
                      coefw_cld_cv, &                                      ! RomP, AJ
                      epmax_diag1)  ! epmax_cape
! **************************************************************
! *
! CV_DRIVER                                                   *
! *
! *
! written by   : Sandrine Bony-Lena , 17/05/2003, 11.19.41    *
! modified by :                                               *
! **************************************************************
! **************************************************************

  USE print_control_mod, ONLY: prt_level, lunout
  USE add_phys_tend_mod, ONLY: fl_cor_ebil
  IMPLICIT NONE

! .............................START PROLOGUE............................


! All argument names (except len,nd,ntra,nloc,delt and the flags) have a "1" appended.
! The "1" is removed for the corresponding compressed variables.
! PARAMETERS:
! Name            Type         Usage            Description
! ----------      ----------     -------  ----------------------------

! len           Integer        Input        first (i) dimension
! nd            Integer        Input        vertical (k) dimension
! ndp1          Integer        Input        nd + 1
! ntra          Integer        Input        number of tracors
! nloc          Integer        Input        dimension of arrays for compressed fields
! k_upper       Integer        Input        upmost level for vertical loops
! iflag_con     Integer        Input        version of convect (3/4)
! iflag_mix     Integer        Input        version of mixing  (0/1/2)
! iflag_ice_thermo Integer        Input        accounting for ice thermodynamics (0/1)
! iflag_clos    Integer        Input        version of closure (0/1)
! tau_cld_cv    Real           Input        characteristic time of dissipation of mixing fluxes
! coefw_cld_cv  Real           Input        coefficient for updraft velocity in convection
! ok_conserv_q  Logical        Input        when true corrections for water conservation are swtiched on
! delt          Real           Input        time step
! comp_threshold Real           Input       threshold on the fraction of convective points below which
!                                            fields  are compressed
! t1            Real           Input        temperature (sat draught envt)
! q1            Real           Input        specific hum (sat draught envt)
! qs1           Real           Input        sat specific hum (sat draught envt)
! t1_wake       Real           Input        temperature (unsat draught envt)
! q1_wake       Real           Input        specific hum(unsat draught envt)
! qs1_wake      Real           Input        sat specific hum(unsat draughts envt)
! s1_wake       Real           Input        fractionnal area covered by wakes
! u1            Real           Input        u-wind
! v1            Real           Input        v-wind
! tra1          Real           Input        tracors
! p1            Real           Input        full level pressure
! ph1           Real           Input        half level pressure
! ALE1          Real           Input        Available lifting Energy
! ALP1          Real           Input        Available lifting Power
! sig1feed1     Real           Input        sigma coord at lower bound of feeding layer
! sig2feed1     Real           Input        sigma coord at upper bound of feeding layer
! wght1         Real           Input        weight density determining the feeding mixture
! iflag1        Integer        Output       flag for Emanuel conditions
! ft1           Real           Output       temp tend
! fq1           Real           Output       spec hum tend
! fqcomp1       Real           Output       spec hum tend (only mixed draughts)
! fu1           Real           Output       u-wind tend
! fv1           Real           Output       v-wind tend
! ftra1         Real           Output       tracor tend
! precip1       Real           Output       precipitation
! kbas1         Integer        Output       cloud base level
! ktop1         Integer        Output       cloud top level
! cbmf1         Real           Output       cloud base mass flux
! sig1          Real           In/Out       section adiabatic updraft
! w01           Real           In/Out       vertical velocity within adiab updraft
! ptop21        Real           In/Out       top of entraining zone
! Ma1           Real           Output       mass flux adiabatic updraft
! mip1          Real           Output       mass flux shed by the adiabatic updraft
! Vprecip1      Real           Output       vertical profile of total precipitation
! Vprecipi1     Real           Output       vertical profile of ice precipitation
! upwd1         Real           Output       total upward mass flux (adiab+mixed)
! dnwd1         Real           Output       saturated downward mass flux (mixed)
! dnwd01        Real           Output       unsaturated downward mass flux
! qcondc1       Real           Output       in-cld mixing ratio of condensed water
! wd1           Real           Output       downdraft velocity scale for sfc fluxes
! cape1         Real           Output       CAPE
! cin1          Real           Output       CIN
! tvp1          Real           Output       adiab lifted parcell virt temp
! ftd1          Real           Output       precip temp tend
! fqt1          Real           Output       precip spec hum tend
! Plim11        Real           Output
! Plim21        Real           Output
! asupmax1      Real           Output
! supmax01      Real           Output
! asupmaxmin1   Real           Output

! ftd1          Real           Output  Array of temperature tendency due to precipitations (K/s) of dimension ND,
!                                      defined at same grid levels as T, Q, QS and P.

! fqd1          Real           Output  Array of specific humidity tendencies due to precipitations ((gm/gm)/s)
!                                      of dimension ND, defined at same grid levels as T, Q, QS and P.

! wdtrainA1     Real           Output   precipitation ejected from adiabatic draught;
!                                         should be used in tracer transport (cvltr)
! wdtrainS1     Real           Output   precipitation detrained from shedding of adiabatic draught;
!                                         used in tracer transport (cvltr)
! wdtrainM1     Real           Output   precipitation detrained from mixed draughts;
!                                         used in tracer transport (cvltr)
! da1           Real           Output     used in tracer transport (cvltr)
! phi1          Real           Output     used in tracer transport (cvltr)
! mp1           Real           Output     used in tracer transport (cvltr)
! qtc1          Real           Output     specific humidity in convection
! sigt1         Real           Output     surface fraction in adiabatic updrafts                                         
! detrain1      Real           Output     detrainment terme klein
! phi21         Real           Output     used in tracer transport (cvltr)
                                         
! d1a1          Real           Output     used in tracer transport (cvltr)
! dam1          Real           Output     used in tracer transport (cvltr)
                                         
! epmlmMm1      Real           Output     used in tracer transport (cvltr)
! eplaMm1       Real           Output     used in tracer transport (cvltr)
                                         
! evap1         Real           Output    
! ep1           Real           Output    
! sigij1        Real           Output     used in tracer transport (cvltr)
! clw1          Real           Output   condensed water content of the adiabatic updraught
! elij1         Real           Output
! wghti1        Real           Output   final weight of the feeding layers,
!                                         used in tracer transport (cvltr)


! S. Bony, Mar 2002:
! * Several modules corresponding to different physical processes
! * Several versions of convect may be used:
!         - iflag_con=3: version lmd  (previously named convect3)
!         - iflag_con=4: version 4.3b (vect. version, previously convect1/2)
! + tard: - iflag_con=5: version lmd with ice (previously named convectg)
! S. Bony, Oct 2002:
! * Vectorization of convect3 (ie version lmd)

! ..............................END PROLOGUE.............................



! Input
  INTEGER, INTENT (IN)                               :: len
  INTEGER, INTENT (IN)                               :: nd
  INTEGER, INTENT (IN)                               :: ndp1
  INTEGER, INTENT (IN)                               :: ntra
  INTEGER, INTENT(IN)                                :: nloc ! (nloc=len)  pour l'instant
  INTEGER, INTENT (IN)                               :: k_upper
  INTEGER, INTENT (IN)                               :: iflag_con
  INTEGER, INTENT (IN)                               :: iflag_mix
  INTEGER, INTENT (IN)                               :: iflag_ice_thermo
  INTEGER, INTENT (IN)                               :: iflag_clos
  LOGICAL, INTENT (IN)                               :: ok_conserv_q
  REAL, INTENT (IN)                                  :: tau_cld_cv
  REAL, INTENT (IN)                                  :: coefw_cld_cv
  REAL, INTENT (IN)                                  :: delt
  REAL, INTENT (IN)                                  :: comp_threshold
  REAL, DIMENSION (len, nd), INTENT (IN)             :: t1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: q1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: qs1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: t1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: q1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: qs1_wake
  REAL, DIMENSION (len), INTENT (IN)                 :: s1_wake
  REAL, DIMENSION (len, nd), INTENT (IN)             :: u1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: v1
  REAL, DIMENSION (len, nd, ntra), INTENT (IN)       :: tra1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: p1
  REAL, DIMENSION (len, ndp1), INTENT (IN)           :: ph1
  REAL, DIMENSION (len), INTENT (IN)                 :: Ale1
  REAL, DIMENSION (len), INTENT (IN)                 :: Alp1
  REAL, DIMENSION (len, nd), INTENT (IN)             :: omega1
  REAL, INTENT (IN)                                  :: sig1feed1 ! pressure at lower bound of feeding layer
  REAL, INTENT (IN)                                  :: sig2feed1 ! pressure at upper bound of feeding layer
  REAL, DIMENSION (nd), INTENT (IN)                  :: wght1     ! weight density determining the feeding mixture
  INTEGER, DIMENSION (len), INTENT (IN)              :: lalim_conv1

! Input/Output
  REAL, DIMENSION (len, nd), INTENT (INOUT)          :: sig1
  REAL, DIMENSION (len, nd), INTENT (INOUT)          :: w01

! Output
  INTEGER, DIMENSION (len), INTENT (OUT)             :: iflag1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ft1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: fq1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: fqcomp1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: fu1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: fv1
  REAL, DIMENSION (len, nd, ntra), INTENT (OUT)      :: ftra1
  REAL, DIMENSION (len), INTENT (OUT)                :: precip1
  INTEGER, DIMENSION (len), INTENT (OUT)             :: kbas1
  INTEGER, DIMENSION (len), INTENT (OUT)             :: ktop1
  REAL, DIMENSION (len), INTENT (OUT)                :: cbmf1
  REAL, DIMENSION (len), INTENT (OUT)                :: plcl1
  REAL, DIMENSION (len), INTENT (OUT)                :: plfc1
  REAL, DIMENSION (len), INTENT (OUT)                :: wbeff1
  REAL, DIMENSION (len), INTENT (OUT)                :: ptop21
  REAL, DIMENSION (len), INTENT (OUT)                :: sigd1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ma1        ! adiab. asc. mass flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: mip1       ! mass flux shed from adiab. ascent (extensive)
! real Vprecip1(len,nd)
  REAL, DIMENSION (len, ndp1), INTENT (OUT)          :: vprecip1   ! tot precipitation flux (staggered grid)
  REAL, DIMENSION (len, ndp1), INTENT (OUT)          :: vprecipi1  ! ice precipitation flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: upwd1      ! upwd sat. mass flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: dnwd1      ! dnwd sat. mass flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: dnwd01     ! unsat. mass flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qcondc1    ! max cloud condensate (intensive)  ! cld
  REAL, DIMENSION (len), INTENT (OUT)                :: wd1             ! gust
  REAL, DIMENSION (len), INTENT (OUT)                :: cape1
  REAL, DIMENSION (len), INTENT (OUT)                :: cin1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: tvp1       ! Virt. temp. in the adiab. ascent

!AC!
!!      real da1(len,nd),phi1(len,nd,nd)
!!      real da(len,nd),phi(len,nd,nd)
!AC!
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ftd1       ! Temp. tendency due to the sole unsat. drafts
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: fqd1       ! Moist. tendency due to the sole unsat. drafts
  REAL, DIMENSION (len), INTENT (OUT)                :: Plim11
  REAL, DIMENSION (len), INTENT (OUT)                :: Plim21
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: asupmax1   ! Highest mixing fraction of mixed updraughts
  REAL, DIMENSION (len), INTENT (OUT)                :: supmax01
  REAL, DIMENSION (len), INTENT (OUT)                :: asupmaxmin1
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qtc1    ! in cloud water content (intensive)   ! cld
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: sigt1   ! fract. cloud area (intensive)        ! cld
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: detrain1   ! detrainement term of mixed draughts in environment

! RomP >>>
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: wdtrainA1, wdtrainS1, wdtrainM1 ! precipitation sources (extensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: mp1  ! unsat. mass flux (staggered grid)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: da1  ! detrained mass flux of adiab. asc. air (extensive)
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: phi1 ! mass flux of envt. air in mixed draughts (extensive)
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: epmlmMm1  ! (extensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: eplaMm1   ! (extensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: evap1 ! evaporation rate in precip. downdraft. (intensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: ep1
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: sigij1 ! mass fraction of env. air in mixed draughts (intensive) 
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: elij1! cond. water per unit mass of mixed draughts (intensive) 
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: qta1 ! total water per unit mass of the adiab. asc. (intensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: clw1 ! cond. water per unit mass of the adiab. asc. (intensive) 
!JYG,RL
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: wghti1   ! final weight of the feeding layers (extensive)
!JYG,RL
  REAL, DIMENSION (len, nd, nd), INTENT (OUT)        :: phi21    ! (extensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: d1a1     ! (extensive)
  REAL, DIMENSION (len, nd), INTENT (OUT)            :: dam1     ! (extensive)
! RomP <<<
  REAL, DIMENSION (len ), INTENT (OUT)               :: epmax_diag1      

! -------------------------------------------------------------------
! Prolog by Kerry Emanuel.
! -------------------------------------------------------------------
! --- ARGUMENTS
! -------------------------------------------------------------------
! --- On input:

! t:   Array of absolute temperature (K) of dimension ND, with first
! index corresponding to lowest model level. Note that this array
! will be altered by the subroutine if dry convective adjustment
! occurs and if IPBL is not equal to 0.

! q:   Array of specific humidity (gm/gm) of dimension ND, with first
! index corresponding to lowest model level. Must be defined
! at same grid levels as T. Note that this array will be altered
! if dry convective adjustment occurs and if IPBL is not equal to 0.

! qs:  Array of saturation specific humidity of dimension ND, with first
! index corresponding to lowest model level. Must be defined
! at same grid levels as T. Note that this array will be altered
! if dry convective adjustment occurs and if IPBL is not equal to 0.

! t_wake: Array of absolute temperature (K), seen by unsaturated draughts,
! of dimension ND, with first index corresponding to lowest model level.

! q_wake: Array of specific humidity (gm/gm), seen by unsaturated draughts,
! of dimension ND, with first index corresponding to lowest model level.
! Must be defined at same grid levels as T.

! qs_wake: Array of saturation specific humidity, seen by unsaturated draughts,
! of dimension ND, with first index corresponding to lowest model level.
! Must be defined at same grid levels as T.

! s_wake: Array of fractionnal area occupied by the wakes.

! u:   Array of zonal wind velocity (m/s) of dimension ND, witth first
! index corresponding with the lowest model level. Defined at
! same levels as T. Note that this array will be altered if
! dry convective adjustment occurs and if IPBL is not equal to 0.

! v:   Same as u but for meridional velocity.

! tra: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
! where NTRA is the number of different tracers. If no
! convective tracer transport is needed, define a dummy
! input array of dimension (ND,1). Tracers are defined at
! same vertical levels as T. Note that this array will be altered
! if dry convective adjustment occurs and if IPBL is not equal to 0.

! p:   Array of pressure (mb) of dimension ND, with first
! index corresponding to lowest model level. Must be defined
! at same grid levels as T.

! ph:  Array of pressure (mb) of dimension ND+1, with first index
! corresponding to lowest level. These pressures are defined at
! levels intermediate between those of P, T, Q and QS. The first
! value of PH should be greater than (i.e. at a lower level than)
! the first value of the array P.

! ALE:  Available lifting Energy

! ALP:  Available lifting Power

! nl:  The maximum number of levels to which convection can penetrate, plus 1.
!       NL MUST be less than or equal to ND-1.

! delt: The model time step (sec) between calls to CONVECT

! ----------------------------------------------------------------------------
! ---   On Output:

! iflag: An output integer whose value denotes the following:
!       VALUE   INTERPRETATION
!       -----   --------------
!         0     Moist convection occurs.
!         1     Moist convection occurs, but a CFL condition
!               on the subsidence warming is violated. This
!               does not cause the scheme to terminate.
!         2     Moist convection, but no precip because ep(inb) lt 0.0001
!         3     No moist convection because new cbmf is 0 and old cbmf is 0.
!         4     No moist convection; atmosphere is not
!               unstable
!         6     No moist convection because ihmin le minorig.
!         7     No moist convection because unreasonable
!               parcel level temperature or specific humidity.
!         8     No moist convection: lifted condensation
!               level is above the 200 mb level.
!         9     No moist convection: cloud base is higher
!               then the level NL-1.
!        10     No moist convection: cloud top is too warm.
!        14     No moist convection; atmosphere is very
!               stable (=> no computation)
!

! ft:   Array of temperature tendency (K/s) of dimension ND, defined at same
!       grid levels as T, Q, QS and P.

! fq:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!       defined at same grid levels as T, Q, QS and P.

! fu:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!      defined at same grid levels as T.

! fv:   Same as FU, but for forcing of meridional velocity.

! ftra: Array of forcing of tracer content, in tracer mixing ratio per
!       second, defined at same levels as T. Dimensioned (ND,NTRA).

! precip: Scalar convective precipitation rate (mm/day).

! wd:   A convective downdraft velocity scale. For use in surface
!       flux parameterizations. See convect.ps file for details.

! tprime: A convective downdraft temperature perturbation scale (K).
!         For use in surface flux parameterizations. See convect.ps
!         file for details.

! qprime: A convective downdraft specific humidity
!         perturbation scale (gm/gm).
!         For use in surface flux parameterizations. See convect.ps
!         file for details.

! cbmf: The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!       BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!       ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!       by the calling program between calls to CONVECT.

! det:   Array of detrainment mass flux of dimension ND.
! -------------------------------------------------------------------

! Local (non compressed) arrays


  INTEGER i, k, il
  INTEGER nword1, nword2, nword3, nword4
  INTEGER icbmax
  INTEGER nk1(len)
  INTEGER icb1(len)
  INTEGER icbs1(len)

  LOGICAL ok_inhib ! True => possible inhibition of convection by dryness
  LOGICAL, SAVE :: debut = .TRUE.
!$OMP THREADPRIVATE(debut)

  REAL coef_convective(len)   ! = 1 for convective points, = 0 otherwise
  REAL tnk1(len)
  REAL thnk1(len)
  REAL qnk1(len)
  REAL gznk1(len)
  REAL qsnk1(len)
  REAL unk1(len)
  REAL vnk1(len)
  REAL cpnk1(len)
  REAL hnk1(len)
  REAL pbase1(len)
  REAL buoybase1(len)

  REAL lf1(len, nd), lf1_wake(len, nd)
  REAL lv1(len, nd), lv1_wake(len, nd)
  REAL cpn1(len, nd), cpn1_wake(len, nd)
  REAL tv1(len, nd), tv1_wake(len, nd)
  REAL gz1(len, nd), gz1_wake(len, nd)
  REAL hm1(len, nd)
  REAL h1(len, nd), h1_wake(len, nd)
  REAL tp1(len, nd)
  REAL th1(len, nd), th1_wake(len, nd)

  REAL bid(len, nd) ! dummy array

  INTEGER ncum

  REAL p1feed1(len) ! pressure at lower bound of feeding layer
  REAL p2feed1(len) ! pressure at upper bound of feeding layer
!JYG,RL
!!      real wghti1(len,nd) ! weights of the feeding layers
!JYG,RL

! (local) compressed fields:


  INTEGER idcum(nloc)
!jyg<
  LOGICAL compress    ! True if compression occurs
!>jyg
  INTEGER iflag(nloc), nk(nloc), icb(nloc)
  INTEGER nent(nloc, nd)
  INTEGER icbs(nloc)
  INTEGER inb(nloc), inbis(nloc)

  REAL cbmf(nloc), plcl(nloc), plfc(nloc), wbeff(nloc)
  REAL t(nloc, nd), q(nloc, nd), qs(nloc, nd)
  REAL t_wake(nloc, nd), q_wake(nloc, nd), qs_wake(nloc, nd)
  REAL s_wake(nloc)
  REAL u(nloc, nd), v(nloc, nd)
  REAL gz(nloc, nd), h(nloc, nd)
  REAL h_wake(nloc, nd)
  REAL lv(nloc, nd), lf(nloc, nd), cpn(nloc, nd)
  REAL lv_wake(nloc, nd), lf_wake(nloc, nd), cpn_wake(nloc, nd)
  REAL p(nloc, nd), ph(nloc, nd+1), tv(nloc, nd), tp(nloc, nd)
  REAL tv_wake(nloc, nd)
  REAL clw(nloc, nd)
  REAL, DIMENSION(nloc, nd)    :: qta, qpreca                       !!jygprl
  REAL dph(nloc, nd)
  REAL pbase(nloc), buoybase(nloc), th(nloc, nd)
  REAL th_wake(nloc, nd)
  REAL tvp(nloc, nd)
  REAL sig(nloc, nd), w0(nloc, nd), ptop2(nloc)
  REAL hp(nloc, nd), ep(nloc, nd), sigp(nloc, nd)
  REAL buoy(nloc, nd)
  REAL cape(nloc)
  REAL cin(nloc)
  REAL m(nloc, nd)
  REAL mm(nloc, nd)
  REAL ment(nloc, nd, nd), sigij(nloc, nd, nd)
  REAL qent(nloc, nd, nd)
  REAL hent(nloc, nd, nd)
  REAL uent(nloc, nd, nd), vent(nloc, nd, nd)
  REAL ments(nloc, nd, nd), qents(nloc, nd, nd)
  REAL elij(nloc, nd, nd)
  REAL supmax(nloc, nd)
  REAL Ale(nloc), Alp(nloc), coef_clos(nloc)
  REAL omega(nloc,nd)
  REAL sigd(nloc)
! real mp(nloc,nd), qp(nloc,nd), up(nloc,nd), vp(nloc,nd)
! real wt(nloc,nd), water(nloc,nd), evap(nloc,nd), ice(nloc,nd)
! real b(nloc,nd), sigd(nloc)
! save mp,qp,up,vp,wt,water,evap,b
  REAL, DIMENSION(len,nd)     :: mp, qp, up, vp
  REAL, DIMENSION(len,nd)     :: wt, water, evap
  REAL, DIMENSION(len,nd)     :: ice, fondue, b
  REAL, DIMENSION(len,nd)     :: frac_a, frac_s, faci               !!jygprl
  REAL ft(nloc, nd), fq(nloc, nd), fqcomp(nloc, nd)
  REAL ftd(nloc, nd), fqd(nloc, nd)
  REAL fu(nloc, nd), fv(nloc, nd)
  REAL upwd(nloc, nd), dnwd(nloc, nd), dnwd0(nloc, nd)
  REAL ma(nloc, nd), mip(nloc, nd)
!!  REAL tls(nloc, nd), tps(nloc, nd)                 ! unused . jyg
  REAL qprime(nloc), tprime(nloc)
  REAL precip(nloc)
! real Vprecip(nloc,nd)
  REAL vprecip(nloc, nd+1)
  REAL vprecipi(nloc, nd+1)
  REAL tra(nloc, nd, ntra), trap(nloc, nd, ntra)
  REAL ftra(nloc, nd, ntra), traent(nloc, nd, nd, ntra)
  REAL qcondc(nloc, nd)      ! cld
  REAL wd(nloc)                ! gust
  REAL Plim1(nloc), plim2(nloc)
  REAL asupmax(nloc, nd)
  REAL supmax0(nloc)
  REAL asupmaxmin(nloc)

  REAL tnk(nloc), qnk(nloc), gznk(nloc)
  REAL wghti(nloc, nd)
  REAL hnk(nloc), unk(nloc), vnk(nloc)

  REAL qtc(nloc, nd)         ! cld
  REAL sigt(nloc, nd)        ! cld
  REAL detrain(nloc, nd)     ! cld
 
! RomP >>>
  REAL wdtrainA(nloc, nd), wdtrainS(nloc, nd), wdtrainM(nloc, nd)   !!jygprl
  REAL da(len, nd), phi(len, nd, nd)
  REAL epmlmMm(nloc, nd, nd), eplaMm(nloc, nd)
  REAL phi2(len, nd, nd)
  REAL d1a(len, nd), dam(len, nd)
! RomP <<<
  REAL epmax_diag(nloc) ! epmax_cape

  CHARACTER (LEN=20) :: modname = 'cva_driver'
  CHARACTER (LEN=80) :: abort_message

  REAL, PARAMETER    :: Cin_noconv = -100000.
  REAL, PARAMETER    :: Cape_noconv = -1.

  INTEGER,SAVE                                       :: igout=1
!$OMP THREADPRIVATE(igout)


! print *, 't1, t1_wake ',(k,t1(1,k),t1_wake(1,k),k=1,nd)
! print *, 'q1, q1_wake ',(k,q1(1,k),q1_wake(1,k),k=1,nd)

! -------------------------------------------------------------------
! --- SET CONSTANTS AND PARAMETERS
! -------------------------------------------------------------------

! -- set simulation flags:
! (common cvflag)

  CALL cv_flag(iflag_ice_thermo)

! -- set thermodynamical constants:
! (common cvthermo)

  CALL cv_thermo(iflag_con)

! -- set convect parameters

! includes microphysical parameters and parameters that
! control the rate of approach to quasi-equilibrium)
! (common cvparam)

  IF (iflag_con==3) THEN
    CALL cv3_param(nd, k_upper, delt)

  END IF

  IF (iflag_con==4) THEN
    CALL cv_param(nd)
  END IF

! ---------------------------------------------------------------------
! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
! ---------------------------------------------------------------------
  nword1 = len
  nword2 = len*nd
  nword3 = len*nd*ntra
  nword4 = len*nd*nd

  iflag1(:) = 0
  ktop1(:) = 0
  kbas1(:) = 0
  ft1(:, :) = 0.0
  fq1(:, :) = 0.0
  fqcomp1(:, :) = 0.0
  fu1(:, :) = 0.0
  fv1(:, :) = 0.0
  ftra1(:, :, :) = 0.
  precip1(:) = 0.
  cbmf1(:) = 0.
  plcl1(:) = 0.
  plfc1(:) = 0.
  wbeff1(:) = 0.
  ptop21(:) = 0.
  sigd1(:) = 0.
  ma1(:, :) = 0.
  mip1(:, :) = 0.
  vprecip1(:, :) = 0.
  vprecipi1(:, :) = 0.
  upwd1(:, :) = 0.
  dnwd1(:, :) = 0.
  dnwd01(:, :) = 0.
  qcondc1(:, :) = 0.
  wd1(:) = 0.
  cape1(:) = 0.
  cin1(:) = 0.
  tvp1(:, :) = 0.
  ftd1(:, :) = 0.
  fqd1(:, :) = 0.
  Plim11(:) = 0.
  Plim21(:) = 0.
  asupmax1(:, :) = 0.
  supmax01(:) = 0.
  asupmaxmin1(:) = 0.

  tvp(:, :) = 0. !ym missing init, need to have a look by developpers
  tv(:, :) = 0. !ym missing init, need to have a look by developpers

  DO il = 1, len
!!    cin1(il) = -100000.
!!    cape1(il) = -1.
    cin1(il) = Cin_noconv
    cape1(il) = Cape_noconv
  END DO

!!  IF (iflag_con==3) THEN
!!    DO il = 1, len
!!      sig1(il, nd) = sig1(il, nd) + 1.
!!      sig1(il, nd) = amin1(sig1(il,nd), 12.1)
!!    END DO
!!  END IF

  IF (iflag_con==3) THEN
      CALL cv3_incrcount(len,nd,delt,sig1)
  END IF  ! (iflag_con==3)

! RomP >>>
  sigt1(:, :) = 0.
  detrain1(:, :) = 0.
  qtc1(:, :) = 0.
  wdtrainA1(:, :) = 0.
  wdtrainS1(:, :) = 0.
  wdtrainM1(:, :) = 0.
  da1(:, :) = 0.
  phi1(:, :, :) = 0.
  epmlmMm1(:, :, :) = 0.
  eplaMm1(:, :) = 0.
  mp1(:, :) = 0.
  evap1(:, :) = 0.
  ep1(:, :) = 0.
  sigij1(:, :, :) = 0.
  elij1(:, :, :) = 0.
  qta1(:,:) = 0.
  clw1(:,:) = 0.
  wghti1(:,:) = 0.
  phi21(:, :, :) = 0.
  d1a1(:, :) = 0.
  dam1(:, :) = 0.
! RomP <<<
! ---------------------------------------------------------------------
! --- INITIALIZE LOCAL ARRAYS AND PARAMETERS
! ---------------------------------------------------------------------

  DO il = 1, nloc
    coef_clos(il) = 1.
  END DO

! --------------------------------------------------------------------
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
! --------------------------------------------------------------------

  IF (iflag_con==3) THEN

    IF (debut) THEN
      PRINT *, 'Emanuel version 3 nouvelle'
    END IF
! print*,'t1, q1 ',t1,q1
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_prelim'
    CALL cv3_prelim(len, nd, ndp1, t1, q1, p1, ph1, &           ! nd->na
                    lv1, lf1, cpn1, tv1, gz1, h1, hm1, th1)


        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_prelim'
    CALL cv3_prelim(len, nd, ndp1, t1_wake, q1_wake, p1, ph1, & ! nd->na
                    lv1_wake, lf1_wake, cpn1_wake, tv1_wake, gz1_wake, &
                    h1_wake, bid, th1_wake)

  END IF

  IF (iflag_con==4) THEN
    PRINT *, 'Emanuel version 4 '
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_prelim'
    CALL cv_prelim(len, nd, ndp1, t1, q1, p1, ph1, &
                   lv1, cpn1, tv1, gz1, h1, hm1)
  END IF

! --------------------------------------------------------------------
! --- CONVECTIVE FEED
! --------------------------------------------------------------------

! compute feeding layer potential temperature and mixing ratio :

! get bounds of feeding layer

! test niveaux couche alimentation KE
  IF (sig1feed1==sig2feed1) THEN
    WRITE (lunout, *) 'impossible de choisir sig1feed=sig2feed'
    WRITE (lunout, *) 'changer la valeur de sig2feed dans physiq.def'
    abort_message = ''
    CALL abort_physic(modname, abort_message, 1)
  END IF

  DO i = 1, len
    p1feed1(i) = sig1feed1*ph1(i, 1)
    p2feed1(i) = sig2feed1*ph1(i, 1)
!test maf
!   p1feed1(i)=ph1(i,1)
!   p2feed1(i)=ph1(i,2)
!   p2feed1(i)=ph1(i,3)
!testCR: on prend la couche alim des thermiques
!   p2feed1(i)=ph1(i,lalim_conv1(i)+1)
!   print*,'lentr=',lentr(i),ph1(i,lentr(i)+1),ph1(i,2)
  END DO

  IF (iflag_con==3) THEN
  END IF
  DO i = 1, len
! print*,'avant cv3_feed Plim',p1feed1(i),p2feed1(i)
  END DO
  IF (iflag_con==3) THEN

! print*, 'IFLAG1 avant cv3_feed'
! print*,'len,nd',len,nd
! write(*,'(64i1)') iflag1(2:len-1)

        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_feed'
    CALL cv3_feed(len, nd, ok_conserv_q, &                 ! nd->na
                  t1, q1, u1, v1, p1, ph1, h1, gz1, & 
                  p1feed1, p2feed1, wght1, &
                  wghti1, tnk1, thnk1, qnk1, qsnk1, unk1, vnk1, &
                  cpnk1, hnk1, nk1, icb1, icbmax, iflag1, gznk1, plcl1)
  END IF

! print*, 'IFLAG1 apres cv3_feed'
! print*,'len,nd',len,nd
! write(*,'(64i1)') iflag1(2:len-1)

  IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_feed'
    CALL cv_feed(len, nd, t1, q1, qs1, p1, hm1, gz1, &
                 nk1, icb1, icbmax, iflag1, tnk1, qnk1, gznk1, plcl1)
  END IF

! print *, 'cv3_feed-> iflag1, plcl1 ',iflag1(1),plcl1(1)

! --------------------------------------------------------------------
! --- UNDILUTE (ADIABATIC) UPDRAFT / 1st part
! (up through ICB for convect4, up through ICB+1 for convect3)
! Calculates the lifted parcel virtual temperature at nk, the
! actual temperature, and the adiabatic liquid water content.
! --------------------------------------------------------------------

  IF (iflag_con==3) THEN

        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_undilute1'
    CALL cv3_undilute1(len, nd, t1, qs1, gz1, plcl1, p1, icb1, tnk1, qnk1, & ! nd->na
                       gznk1, tp1, tvp1, clw1, icbs1)
  END IF


  IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_undilute1'
    CALL cv_undilute1(len, nd, t1, q1, qs1, gz1, p1, nk1, icb1, icbmax, &
                      tp1, tvp1, clw1)
  END IF

! -------------------------------------------------------------------
! --- TRIGGERING
! -------------------------------------------------------------------

! print *,' avant triggering, iflag_con ',iflag_con

  IF (iflag_con==3) THEN

        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_trigger'
    CALL cv3_trigger(len, nd, icb1, plcl1, p1, th1, tv1, tvp1, thnk1, & ! nd->na
                      pbase1, buoybase1, iflag1, sig1, w01)


! print*, 'IFLAG1 apres cv3_triger'
! print*,'len,nd',len,nd
! write(*,'(64i1)') iflag1(2:len-1)

! call dump2d(iim,jjm-1,sig1(2)
  END IF

  IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_trigger'
    CALL cv_trigger(len, nd, icb1, cbmf1, tv1, tvp1, iflag1)
  END IF


! =====================================================================
! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
! =====================================================================

!  Determine the number "ncum" of convective gridpoints, the list "idcum" of convective
!  gridpoints and the weights "coef_convective" (= 1. for convective gridpoints and = 0.
!  elsewhere).
  ncum = 0
  coef_convective(:) = 0.
  DO i = 1, len
    IF (iflag1(i)==0) THEN
      coef_convective(i) = 1.
      ncum = ncum + 1
      idcum(ncum) = i
    END IF
  END DO

! print*,'len, ncum = ',len,ncum

  IF (ncum>0) THEN

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- COMPRESS THE FIELDS
!       (-> vectorization over convective gridpoints)
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IF (iflag_con==3) THEN
! print*,'ncum tv1 ',ncum,tv1
! print*,'tvp1 ',tvp1
!jyg<
!   If the fraction of convective points is larger than comp_threshold, then compression
!   is assumed useless.
!
  compress = ncum .lt. len*comp_threshold
!
  IF (.not. compress) THEN
    DO i = 1,len
      idcum(i) = i
    ENDDO
  ENDIF
!
!>jyg
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3a_compress'
      CALL cv3a_compress(len, nloc, ncum, nd, ntra, compress, &
                         iflag1, nk1, icb1, icbs1, &
                         plcl1, tnk1, qnk1, gznk1, hnk1, unk1, vnk1, &
                         wghti1, pbase1, buoybase1, &
                         t1, q1, qs1, t1_wake, q1_wake, qs1_wake, s1_wake, &
                         u1, v1, gz1, th1, th1_wake, &
                         tra1, &
                         h1, lv1, lf1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
                         h1_wake, lv1_wake, lf1_wake, cpn1_wake, tv1_wake, &
                         sig1, w01, ptop21, &
                         Ale1, Alp1, omega1, &
                         iflag, nk, icb, icbs, &
                         plcl, tnk, qnk, gznk, hnk, unk, vnk, &
                         wghti, pbase, buoybase, &
                         t, q, qs, t_wake, q_wake, qs_wake, s_wake, &
                         u, v, gz, th, th_wake, &
                         tra, &
                         h, lv, lf, cpn, p, ph, tv, tp, tvp, clw, &
                         h_wake, lv_wake, lf_wake, cpn_wake, tv_wake, &
                         sig, w0, ptop2, &
                         Ale, Alp, omega)

! print*,'tv ',tv
! print*,'tvp ',tvp

    END IF

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_compress'
      CALL cv_compress(len, nloc, ncum, nd, &
                       iflag1, nk1, icb1, &
                       cbmf1, plcl1, tnk1, qnk1, gznk1, &
                       t1, q1, qs1, u1, v1, gz1, &
                       h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
                       iflag, nk, icb, &
                       cbmf, plcl, tnk, qnk, gznk, &
                       t, q, qs, u, v, gz, h, lv, cpn, p, ph, tv, tp, tvp, clw, &
                       dph)
    END IF

! -------------------------------------------------------------------
! --- UNDILUTE (ADIABATIC) UPDRAFT / second part :
! ---   FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
! ---   &
! ---   COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
! ---   FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
! ---   &
! ---   FIND THE LEVEL OF NEUTRAL BUOYANCY
! -------------------------------------------------------------------

    IF (iflag_con==3) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_undilute2'
      CALL cv3_undilute2(nloc, ncum, nd, iflag, icb, icbs, nk, &        !na->nd
                         tnk, qnk, gznk, hnk, t, q, qs, gz, &
                         p, ph, h, tv, lv, lf, pbase, buoybase, plcl, &
                         inb, tp, tvp, clw, hp, ep, sigp, buoy, &
                         frac_a, frac_s, qpreca, qta)                        !!jygprl
    END IF

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_undilute2'
      CALL cv_undilute2(nloc, ncum, nd, icb, nk, &
                        tnk, qnk, gznk, t, q, qs, gz, &
                        p, dph, h, tv, lv, &
                        inb, inbis, tp, tvp, clw, hp, ep, sigp, frac_s)
    END IF

    ! epmax_cape
    ! on recalcule ep et hp    
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_epmax_cape'
    call cv3_epmax_fn_cape(nloc,ncum,nd &
                , ep,hp,icb,inb,clw,nk,t,h,hnk,lv,lf,frac_s &
                , pbase, p, ph, tv, buoy, sig, w0,iflag &
                , epmax_diag)

! -------------------------------------------------------------------
! --- MIXING(1)   (if iflag_mix .ge. 1)
! -------------------------------------------------------------------
    IF (iflag_con==3) THEN
!      IF ((iflag_ice_thermo==1) .AND. (iflag_mix/=0)) THEN
!        WRITE (*, *) ' iflag_ice_thermo==1 requires iflag_mix==0', ' but iflag_mix=', iflag_mix, &
!          '. Might as well stop here.'
!        STOP
!      END IF
      IF (iflag_mix>=1) THEN
        CALL zilch(supmax, nloc*nd)
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3p_mixing'
        CALL cv3p_mixing(nloc, ncum, nd, nd, ntra, icb, nk, inb, &           ! na->nd
!!                         ph, t, q, qs, u, v, tra, h, lv, lf, frac, qnk, &
                         ph, t, q, qs, u, v, tra, h, lv, lf, frac_s, qta, &      !!jygprl
                         unk, vnk, hp, tv, tvp, ep, clw, sig, &
                         ment, qent, hent, uent, vent, nent, &
                         sigij, elij, supmax, ments, qents, traent)
! print*, 'cv3p_mixing-> supmax ', (supmax(1,k), k=1,nd)

      ELSE
        CALL zilch(supmax, nloc*nd)
      END IF
    END IF
! -------------------------------------------------------------------
! --- CLOSURE
! -------------------------------------------------------------------


    IF (iflag_con==3) THEN
      IF (iflag_clos==0) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_closure'
        CALL cv3_closure(nloc, ncum, nd, icb, inb, &           ! na->nd
                         pbase, p, ph, tv, buoy, &
                         sig, w0, cape, m, iflag)
      END IF   ! iflag_clos==0

      ok_inhib = iflag_mix == 2

      IF (iflag_clos==1) THEN
        PRINT *, ' pas d appel cv3p_closure'
! c        CALL cv3p_closure(nloc,ncum,nd,icb,inb              ! na->nd
! c    :                       ,pbase,plcl,p,ph,tv,tvp,buoy
! c    :                       ,supmax
! c    o                       ,sig,w0,ptop2,cape,cin,m)
      END IF   ! iflag_clos==1

      IF (iflag_clos==2) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3p1_closure'
        CALL cv3p1_closure(nloc, ncum, nd, icb, inb, &         ! na->nd
                           pbase, plcl, p, ph, tv, tvp, buoy, &
                           supmax, ok_inhib, Ale, Alp, omega, &
                           sig, w0, ptop2, cape, cin, m, iflag, coef_clos, &
                           Plim1, plim2, asupmax, supmax0, &
                           asupmaxmin, cbmf, plfc, wbeff)
        if (prt_level >= 10) &
             PRINT *, 'cv3p1_closure-> plfc,wbeff ', plfc(1), wbeff(1)
      END IF   ! iflag_clos==2

      IF (iflag_clos==3) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3p2_closure'
        CALL cv3p2_closure(nloc, ncum, nd, icb, inb, &         ! na->nd
                           pbase, plcl, p, ph, tv, tvp, buoy, &
                           supmax, ok_inhib, Ale, Alp, omega, &
                           sig, w0, ptop2, cape, cin, m, iflag, coef_clos, &
                           Plim1, plim2, asupmax, supmax0, &
                           asupmaxmin, cbmf, plfc, wbeff)
        if (prt_level >= 10) &
             PRINT *, 'cv3p2_closure-> plfc,wbeff ', plfc(1), wbeff(1)
      END IF   ! iflag_clos==3
    END IF ! iflag_con==3

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_closure'
      CALL cv_closure(nloc, ncum, nd, nk, icb, &
                         tv, tvp, p, ph, dph, plcl, cpn, &
                         iflag, cbmf)
    END IF

! print *,'cv_closure-> cape ',cape(1)

! -------------------------------------------------------------------
! --- MIXING(2)
! -------------------------------------------------------------------

    IF (iflag_con==3) THEN
      IF (iflag_mix==0) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_mixing'
        CALL cv3_mixing(nloc, ncum, nd, nd, ntra, icb, nk, inb, &             ! na->nd
                        ph, t, q, qs, u, v, tra, h, lv, lf, frac_s, qnk, &
                        unk, vnk, hp, tv, tvp, ep, clw, m, sig, &
                        ment, qent, uent, vent, nent, sigij, elij, ments, qents, traent)
        CALL zilch(hent, nloc*nd*nd)
      ELSE
!!jyg:  Essais absurde pour voir
!!        mm(:,1) = 0.
!!        DO  i = 2,nd
!!          mm(:,i) = m(:,i)*(1.-qta(:,i-1))
!!        ENDDO
        mm(:,:) = m(:,:)
        CALL cv3_mixscale(nloc, ncum, nd, ment, mm)
        IF (debut) THEN
          PRINT *, ' cv3_mixscale-> '
        END IF !(debut) THEN
      END IF
    END IF

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_mixing'
      CALL cv_mixing(nloc, ncum, nd, icb, nk, inb, inbis, &
                     ph, t, q, qs, u, v, h, lv, qnk, &
                     hp, tv, tvp, ep, clw, cbmf, &
                     m, ment, qent, uent, vent, nent, sigij, elij)
    END IF                                                                                         

    IF (debut) THEN
      PRINT *, ' cv_mixing ->'
    END IF !(debut) THEN
! do i = 1,nd
! print*,'cv_mixing-> i,ment ',i,(ment(1,i,j),j=1,nd)
! enddo

! -------------------------------------------------------------------
! --- UNSATURATED (PRECIPITATING) DOWNDRAFTS
! -------------------------------------------------------------------
    IF (iflag_con==3) THEN
      IF (debut) THEN
        PRINT *, ' cva_driver -> cv3_unsat '
      END IF !(debut) THEN

        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_unsat'
      CALL cv3_unsat(nloc, ncum, nd, nd, ntra, icb, inb, iflag, &              ! na->nd
                     t_wake, q_wake, qs_wake, gz, u, v, tra, p, ph, &
                     th_wake, tv_wake, lv_wake, lf_wake, cpn_wake, &
                     ep, sigp, clw, frac_s, qpreca, frac_a, qta, &                    !!jygprl
                     m, ment, elij, delt, plcl, coef_clos, &
                     mp, qp, up, vp, trap, wt, water, evap, fondue, ice, &
                     faci, b, sigd, &
!!                     wdtrainA, wdtrainM)                                       ! RomP
                     wdtrainA, wdtrainS, wdtrainM)                               !!jygprl
!
      IF (prt_level >= 10) THEN
        Print *, 'cva_driver after cv3_unsat:mp , water, ice, evap, fondue '
        DO k = 1,nd
        write (6, '(i4,5(1x,e13.6))'), &
          k, mp(igout,k), water(igout,k), ice(igout,k), &
           evap(igout,k), fondue(igout,k)
        ENDDO
        Print *, 'cva_driver after cv3_unsat: wdtrainA, wdtrainS, wdtrainM '     !!jygprl
        DO k = 1,nd
        write (6, '(i4,3(1x,e13.6))'), &
           k, wdtrainA(igout,k), wdtrainS(igout,k), wdtrainM(igout,k)            !!jygprl
        ENDDO
      ENDIF
!
    END IF  !(iflag_con==3)

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_unsat'
      CALL cv_unsat(nloc, ncum, nd, inb, t, q, qs, gz, u, v, p, ph, &
                     h, lv, ep, sigp, clw, m, ment, elij, &
                     iflag, mp, qp, up, vp, wt, water, evap)
    END IF

    IF (debut) THEN
      PRINT *, 'cv_unsat-> '
    END IF !(debut) THEN

! print *,'cv_unsat-> mp ',mp
! print *,'cv_unsat-> water ',water
! -------------------------------------------------------------------
! --- YIELD
! (tendencies, precipitation, variables of interface with other
! processes, etc)
! -------------------------------------------------------------------

    IF (iflag_con==3) THEN

        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_yield'
      CALL cv3_yield(nloc, ncum, nd, nd, ntra, ok_conserv_q, &                      ! na->nd
                     icb, inb, delt, &
                     t, q, t_wake, q_wake, s_wake, u, v, tra, &
                     gz, p, ph, h, hp, lv, lf, cpn, th, th_wake, &
                     ep, clw, qpreca, m, tp, mp, qp, up, vp, trap, &
                     wt, water, ice, evap, fondue, faci, b, sigd, &
                     ment, qent, hent, iflag_mix, uent, vent, &
                     nent, elij, traent, sig, &
                     tv, tvp, wghti, &
                     iflag, precip, Vprecip, Vprecipi, ft, fq, fqcomp, fu, fv, ftra, &      ! jyg
                     cbmf, upwd, dnwd, dnwd0, ma, mip, &
!!                     tls, tps, &                            ! useless . jyg
                     qcondc, wd, &
!!                     ftd, fqd, qnk, qtc, sigt, tau_cld_cv, coefw_cld_cv)
                     ftd, fqd, qta, qtc, sigt, detrain, tau_cld_cv, coefw_cld_cv)         !!jygprl
!
!         Test conseravtion de l'eau
!
      IF (debut) THEN
        PRINT *, ' cv3_yield -> fqd(1) = ', fqd(igout, 1)
      END IF !(debut) THEN
!   
      IF (prt_level >= 10) THEN
        Print *, 'cva_driver after cv3_yield:ft(1) , ftd(1) ', &
                    ft(igout,1), ftd(igout,1)
        Print *, 'cva_driver after cv3_yield:fq(1) , fqd(1) ', &
                    fq(igout,1), fqd(igout,1)
      ENDIF
!   
    END IF

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_yield'
      CALL cv_yield(nloc, ncum, nd, nk, icb, inb, delt, &
                     t, q, u, v, &
                     gz, p, ph, h, hp, lv, cpn, &
                     ep, clw, frac_s, m, mp, qp, up, vp, &
                     wt, water, evap, &
                     ment, qent, uent, vent, nent, elij, &
                     tv, tvp, &
                     iflag, wd, qprime, tprime, &
                     precip, cbmf, ft, fq, fu, fv, ma, qcondc)
    END IF

!AC!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!--- passive tracers
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IF (iflag_con==3) THEN
!RomP >>>
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3_tracer'
      CALL cv3_tracer(nloc, len, ncum, nd, nd, &
                     ment, sigij, da, phi, phi2, d1a, dam, &
                     ep, vprecip, elij, clw, epmlmMm, eplaMm, &
                     icb, inb)
!RomP <<<
    END IF

!AC!

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! --- UNCOMPRESS THE FIELDS
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    IF (iflag_con==3) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv3a_uncompress'
      CALL cv3a_uncompress(nloc, len, ncum, nd, ntra, idcum, compress, &
                           iflag, icb, inb, &
                           precip, cbmf, plcl, plfc, wbeff, sig, w0, ptop2, &
                           ft, fq, fqcomp, fu, fv, ftra, &
                           sigd, ma, mip, vprecip, vprecipi, upwd, dnwd, dnwd0, &
                           qcondc, wd, cape, cin, &
                           tvp, &
                           ftd, fqd, &
                           Plim1, plim2, asupmax, supmax0, &
                           asupmaxmin, &
                           da, phi, mp, phi2, d1a, dam, sigij, &         ! RomP
                           qta, clw, elij, evap, ep, epmlmMm, eplaMm, &  ! RomP
                           wdtrainA, wdtrainS, wdtrainM, &                         ! RomP
                           qtc, sigt, detrain, epmax_diag, & ! epmax_cape
                           iflag1, kbas1, ktop1, &
                           precip1, cbmf1, plcl1, plfc1, wbeff1, sig1, w01, ptop21, &
                           ft1, fq1, fqcomp1, fu1, fv1, ftra1, &
                           sigd1, ma1, mip1, vprecip1, vprecipi1, upwd1, dnwd1, dnwd01, &
                           qcondc1, wd1, cape1, cin1, &
                           tvp1, &
                           ftd1, fqd1, &
                           Plim11, plim21, asupmax1, supmax01, &
                           asupmaxmin1, &
                           da1, phi1, mp1, phi21, d1a1, dam1, sigij1,  &       ! RomP
                           qta1, clw1, elij1, evap1, ep1, epmlmMm1, eplaMm1, & ! RomP
                           wdtrainA1, wdtrainS1, wdtrainM1,                  & ! RomP
                           qtc1, sigt1, detrain1, epmax_diag1) ! epmax_cape
!   
      IF (prt_level >= 10) THEN
        Print *, 'cva_driver after cv3_uncompress:ft1(1) , ftd1(1) ', &
                    ft1(igout,1), ftd1(igout,1)
        Print *, 'cva_driver after cv3_uncompress:fq1(1) , fqd1(1) ', &
                    fq1(igout,1), fqd1(igout,1)
      ENDIF
!   
    END IF

    IF (iflag_con==4) THEN
        if (prt_level >= 9) &
             PRINT *, 'cva_driver -> cv_uncompress'
      CALL cv_uncompress(nloc, len, ncum, nd, idcum, &
                           iflag, &
                           precip, cbmf, &
                           ft, fq, fu, fv, &
                           ma, qcondc, &
                           iflag1, &
                           precip1,cbmf1, &
                           ft1, fq1, fu1, fv1, &
                           ma1, qcondc1)
    END IF

  END IF ! ncum>0
!
!
  DO i = 1,len
    IF (iflag1(i) == 14) THEN
      Cin1(i) = Cin_noconv
      Cape1(i) = Cape_noconv
    ENDIF
  ENDDO

!
! In order take into account the possibility of changing the compression,
! reset m, sig and w0 to zero for non-convective points.
  DO k = 1,nd-1
        sig1(:, k) = sig1(:, k)*coef_convective(:)
        w01(:, k)  = w01(:, k)*coef_convective(:)
  ENDDO

  IF (debut) THEN
    PRINT *, ' cv_uncompress -> '
    debut = .FALSE.
  END IF  !(debut) THEN


  RETURN
END SUBROUTINE cva_driver
