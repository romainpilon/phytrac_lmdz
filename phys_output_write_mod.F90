!
! $Id: phys_output_write_mod.F90 4703 2023-09-20 13:09:14Z fairhead $
!
MODULE phys_output_write_mod

  USE phytrac_mod, ONLY : d_tr_cl, d_tr_th, d_tr_cv, d_tr_lessi_impa, &
       d_tr_lessi_nucl, d_tr_insc, d_tr_bcscav, d_tr_evapls, d_tr_ls,  &
       d_tr_trsp, d_tr_sscav, d_tr_sat, d_tr_uscav, flux_tr_dry

  ! Author: Abderrahmane IDELKADI (original include file)
  ! Author: Laurent FAIRHEAD (transformation to module/subroutine)
  ! Author: Ulysse GERARD (effective implementation)

CONTAINS 

  ! ug Routine pour définir (lors du premier passage) ET sortir les variables
  SUBROUTINE phys_output_write(itap, pdtphys, paprs, pphis, &
       pplay, lmax_th, aerosol_couple,         &
       ok_ade, ok_aie, ok_volcan, ivap, iliq, isol, ibs, ok_sync, &
       ptconv, read_climoz, clevSTD, ptconvth, &
       d_u, d_t, qx, d_qx, zmasse, flag_aerosol, flag_aerosol_strat, ok_cdnc, t, u1, v1)

    ! This subroutine does the actual writing of diagnostics that were
    ! defined and initialised in phys_output_mod.F90

    USE dimphy, ONLY: klon, klev, klevp1
    USE infotrac_phy, ONLY: nbtr, nqtot, nqo, type_trac, tracers, niso, ntiso
    USE strings_mod,  ONLY: maxlen
    USE mod_phys_lmdz_para, ONLY: is_north_pole_phy,is_south_pole_phy
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
    USE time_phylmdz_mod, ONLY: day_step_phy, start_time, itau_phy
    USE vertical_layers_mod, ONLY : ap, bp, aps, bps
    USE phystokenc_mod, ONLY: offline
    USE phys_output_ctrlout_mod, ONLY: o_phis, o_aire, is_ter, is_lic, is_oce, &
         o_longitude, o_latitude, &
         o_Ahyb, o_Bhyb,o_Ahyb_bounds, o_Bhyb_bounds, & 
         o_Ahyb_mid, o_Bhyb_mid,o_Ahyb_mid_bounds, o_Bhyb_mid_bounds, & 
         is_ave, is_sic, o_contfracATM, o_contfracOR, &
         o_aireTER, o_flat, o_slp, o_ptstar, o_pt0, o_tsol, &
         o_t2m, o_t2m_min, o_t2m_max, &
         o_t2m_min_mon, o_t2m_max_mon, &
         o_nt2mout, o_nt2moutfg, &
         o_nq2mout, o_nq2moutfg, &
         o_nu2mout, o_nu2moutfg, &
         o_q2m, o_ustar, o_u10m, o_v10m, &
         o_wind10m, o_wind10max, o_wind100m, o_gusts, o_sicf, &
         o_loadfactor_wind_onshore, o_loadfactor_wind_offshore, &
         o_psol, o_mass, o_qsurf, o_qsol, &
         o_precip, o_rain_fall, o_rain_con, o_ndayrain, o_plul, o_pluc, o_plun, &
         o_snow, o_msnow, o_fsnow, o_evap, o_snowerosion, o_ustart_lic, o_rhosnow_lic, o_bsfall, & 
         o_ep,o_epmax_diag, & ! epmax_cape
         o_tops, o_tops0, o_topl, o_topl0, &
         o_SWupTOA, o_SWupTOAclr, o_SWupTOAcleanclr, o_SWdnTOA, o_fdiffSWdnSFC, &
         o_SWdnTOAclr, o_nettop, o_SWup200, &
         o_SWup200clr, o_SWdn200, o_SWdn200clr, &
         o_LWup200, o_LWup200clr, o_LWdn200, &
         o_LWdn200clr, o_sols, o_sols0, &
         o_soll, o_radsol, o_soll0, o_SWupSFC, &
         o_SWupSFCclr, o_SWupSFCcleanclr, o_SWdnSFC, o_SWdnSFCclr, o_SWdnSFCcleanclr, &
         o_LWupSFC, o_LWdnSFC, o_LWupSFCclr, &
         o_LWdnSFCclr, o_LWupTOAcleanclr, o_LWdnSFCcleanclr, o_bils, o_bils_diss, &
         o_bils_ec,o_bils_ech, o_bils_tke, o_bils_kinetic, &
         o_bils_latent, o_bils_enthalp, o_sens, &
         o_fder, o_ffonte, o_fqcalving, o_fqfonte, o_mrroli, o_runofflic, &
         o_taux, o_tauy, o_snowsrf, o_qsnow, &
         o_snowhgt, o_toice, o_sissnow, o_runoff, &
         o_albslw3, o_pourc_srf, o_fract_srf, &
         o_taux_srf, o_tauy_srf, o_tsol_srf, &
         o_evappot_srf, o_ustar_srf, o_u10m_srf, &
         o_v10m_srf, o_t2m_srf, o_evap_srf, &
         o_sens_srf, o_lat_srf, o_flw_srf, &
         o_fsw_srf, o_wbils_srf, o_wbilo_srf, &
         o_wevap_srf, o_wrain_srf, o_wsnow_srf, &
         o_tke_srf, o_tke_max_srf,o_dltpbltke_srf, o_wstar, &
         o_l_mixmin,o_l_mix, &
         o_cdrm, o_cdrh, o_cldl, o_cldm, o_cldh, &
         o_cldt, o_JrNt, o_cldljn, o_cldmjn, &
         o_cldhjn, o_cldtjn, o_cldq, o_lwp, o_iwp, &
         o_ue, o_ve, o_uq, o_vq, o_cape, o_pbase, &
         o_uwat, o_vwat, &
         o_ptop, o_fbase, o_plcl, o_plfc, &
         o_wbeff, o_convoccur, o_cape_max, o_upwd, o_ep,o_epmax_diag, &
         o_Mipsh, o_Ma, &
         o_dnwd, o_dnwd0, o_ftime_deepcv, o_ftime_con, o_mc, &
         o_prw, o_prlw, o_prsw, o_prbsw, o_s_pblh, o_s_pblt, o_s_lcl, &
         o_s_therm, o_uSTDlevs, o_vSTDlevs, &
         o_wSTDlevs, o_zSTDlevs, o_qSTDlevs, &
         o_tSTDlevs, epsfra, o_t_oce_sic, &
         o_ale_bl, o_alp_bl, o_ale_wk, o_alp_wk, &
         o_dtvdf_x    , o_dtvdf_w    , o_dqvdf_x    , o_dqvdf_w    , &
         o_sens_x     , o_sens_w     , o_flat_x     , o_flat_w     , &
         o_delta_tsurf, o_delta_tsurf_srf, &
         o_cdragh_x   , o_cdragh_w   , o_cdragm_x   , o_cdragm_w   , &
         o_kh         , o_kh_x       , o_kh_w       , &
         o_ale, o_alp, o_cin, o_WAPE, o_wake_h, o_cv_gen, o_wake_dens, &
         o_wake_s, o_wake_deltat, o_wake_deltaq, &
         o_wake_omg, o_dtwak, o_dqwak, o_dqwak2d, o_Vprecip, &
         o_qtaa, o_Clwaa, &
         o_ftd, o_fqd, o_wdtrainA, o_wdtrainS, o_wdtrainM, &
         o_n2, o_s2, o_proba_notrig, &
         o_random_notrig, o_ale_bl_stat, &
         o_ale_bl_trig, o_alp_bl_det, &
         o_alp_bl_fluct_m, o_alp_bl_fluct_tke, &
         o_alp_bl_conv, o_alp_bl_stat, &
         o_slab_qflux, o_tslab, o_slab_bils, &
         o_slab_bilg, o_slab_sic, o_slab_tice, &
         o_slab_hdiff, o_slab_ekman, o_slab_gm,  &
         o_weakinv, o_dthmin, o_cldtau, &
         o_cldemi, o_pr_con_l, o_pr_con_i, &
         o_pr_lsc_l, o_pr_lsc_i, o_pr_bs, o_re, o_fl, &
         o_rh2m, &
         o_qsat2m, o_tpot, o_tpote, o_SWnetOR, &
         o_LWdownOR, o_snowl, &
         o_solldown, o_dtsvdfo, o_dtsvdft, &
         o_dtsvdfg, o_dtsvdfi, o_z0m, o_z0h,  o_od443aer, o_od550aer, &
         o_dryod550aer, o_od865aer, o_abs550aer, o_od550lt1aer, &
         o_sconcso4, o_sconcno3, o_sconcoa, o_sconcbc, &
         o_sconcss, o_sconcdust, o_concso4, o_concno3, &
         o_concoa, o_concbc, o_concss, o_concdust, &
         o_loadso4, o_loadoa, o_loadbc, o_loadss, &
         o_loaddust, o_loadno3, o_tausumaero, & 
         o_drytausumaero, o_tausumaero_lw, &
         o_topswad, o_topswad0, o_solswad, o_solswad0, &
         o_toplwad, o_toplwad0, o_sollwad, o_sollwad0, &
         o_swtoaas_nat, o_swsrfas_nat, &
         o_swtoacs_nat, o_swtoaas_ant, &
         o_swsrfas_ant, o_swtoacs_ant, &
         o_swsrfcs_ant, o_swtoacf_nat, &
         o_swsrfcf_nat, o_swtoacf_ant, &
         o_swsrfcs_nat, o_swsrfcf_ant, &
         o_swtoacf_zero, o_swsrfcf_zero, &
         o_topswai, o_solswai, o_toplwai, o_sollwai, o_scdnc, &
         o_cldncl, o_reffclws, o_reffclwc, o_solbnd, o_stratomask,&
         o_cldnvi, o_lcc, o_lcc3d, o_lcc3dcon, &
         o_lcc3dstra, o_icc3dcon, o_icc3dstra, &
         o_cldicemxrat, o_cldwatmxrat, o_reffclwtop, o_ec550aer, &
         o_lwcon, o_iwcon, o_temp, o_theta, &
         o_ovapinit, o_ovap, o_oliq, o_ocond, o_geop,o_qbs, &
         o_vitu, o_vitv, o_vitw, o_pres, o_paprs, &
         o_zfull, o_zhalf, o_rneb, o_rnebjn, o_rnebcon, &
         o_rnebls, o_rneblsvol, o_rhum, o_rhl, o_rhi, o_ozone, o_ozone_light, &
         o_pfraclr, o_pfracld, &
         o_duphy, o_dtphy, o_dqphy, o_dqphy2d, o_dqlphy, o_dqlphy2d, &
         o_dqsphy, o_dqsphy2d, o_dqbsphy, o_dqbsphy2d, o_albe_srf, o_z0m_srf, o_z0h_srf, &
         o_ages_srf, o_snow_srf, o_alb1, o_alb2, o_tke, o_tke_dissip, &
         o_tke_max, o_kz, o_kz_max, o_clwcon, &
         o_dtdyn, o_dqdyn, o_dqdyn2d, o_dqldyn, o_dqldyn2d, & 
         o_dqsdyn, o_dqsdyn2d, o_dqbsdyn, o_dqbsdyn2d, o_dudyn, o_dvdyn, &
         o_dtcon, o_tntc, o_ducon, o_dvcon, &
         o_dqcon, o_dqcon2d, o_tnhusc, o_tnhusc, o_dtlsc, &
         o_dtlschr, o_dqlsc, o_dqlsc2d, o_beta_prec, &
         o_dtlscth, o_dtlscst, o_dqlscth, o_dqlscth2d, &
         o_dqlscst, o_dqlscst2d, o_plulth, o_plulst, &
         o_ptconvth, o_lmaxth, o_dtvdf, &
         o_dtdis, o_dqvdf, o_dqvdf2d, o_dteva, o_dqeva, o_dqeva2d, &
         o_dqbsvdf, o_dtbs, o_dqbs, o_dqbsbs, &
         o_ptconv, o_ratqs, o_dtthe, & 
         o_duthe, o_dvthe, o_ftime_th, &
         o_f_th, o_e_th, o_w_th, o_q_th, &
         o_a_th, o_cloudth_sth, o_cloudth_senv, &
         o_cloudth_sigmath, o_cloudth_sigmaenv, &
         o_d_th, o_f0_th, o_zmax_th, &
         o_dqthe, o_dqthe2d, o_dtajs, o_dqajs, o_dqajs2d, o_dtswr, &
         o_dtsw0, o_dtlwr, o_dtlw0, o_dtec, &
         o_duvdf, o_dvvdf, o_duoro, o_dvoro, &
         o_dtoro, o_dulif, o_dvlif, o_dtlif, &
         o_du_gwd_hines, o_dv_gwd_hines, o_dthin, o_dqch4, o_rsu, &
         o_du_gwd_front, o_dv_gwd_front, &
         o_east_gwstress, o_west_gwstress, &
         o_rsd, o_rlu, o_rld, o_rsucs, o_rsdcs, o_rsucsaf, o_rsdcsaf, &
         o_rlucs, o_rldcs, o_tnt, o_tntr, &
         o_tntscpbl, o_tnhus, o_tnhusscpbl, &
         o_evu, o_h2o, o_mcd, o_dmc, o_ref_liq, &
         o_ref_ice, o_rsut4co2, o_rlut4co2, &
         o_rsutcs4co2, o_rlutcs4co2, o_rsu4co2, &
         o_rlu4co2, o_rsucs4co2, o_rlucs4co2, &
         o_rsd4co2, o_rld4co2, o_rsdcs4co2, &
         o_rldcs4co2, o_tnondef, o_ta, o_zg, &
         o_hus, o_hur, o_ua, o_va, o_wap, &
         o_psbg, o_tro3, o_tro3_daylight, &
         o_uxv, o_vxq, o_vxT, o_wxq, o_vxphi, &
         o_wxT, o_uxu, o_vxv, o_TxT, o_trac, &
         o_dtr_vdf, o_dtr_the, o_dtr_con, &
         o_dtr_lessi_impa, o_dtr_lessi_nucl, &
         o_dtr_insc, o_dtr_bcscav, o_dtr_evapls, &
         o_dtr_ls, o_dtr_trsp, o_dtr_sscav, o_dtr_dry, &
         o_dtr_sat, o_dtr_uscav, o_trac_cum, o_du_gwd_rando, o_dv_gwd_rando, &
         o_ustr_gwd_hines,o_vstr_gwd_hines,o_ustr_gwd_rando,o_vstr_gwd_rando, &
         o_ustr_gwd_front,o_vstr_gwd_front, &
         o_sens_prec_liq_oce, o_sens_prec_liq_sic, &
         o_sens_prec_sol_oce, o_sens_prec_sol_sic, &
         o_lat_prec_liq_oce, o_lat_prec_liq_sic, &
         o_lat_prec_sol_oce, o_lat_prec_sol_sic, &
         o_sza, &
! Marine
         o_map_prop_hc, o_map_prop_hist, o_map_emis_hc, o_map_iwp_hc, &
         o_map_deltaz_hc, o_map_pcld_hc, o_map_tcld_hc, &
         o_map_emis_hist, o_map_iwp_hist, o_map_deltaz_hist, &
         o_map_rad_hist, &
         o_map_emis_Cb, o_map_pcld_Cb, o_map_tcld_Cb, &
         o_map_emis_ThCi, o_map_pcld_ThCi, o_map_tcld_ThCi, &
         o_map_emis_Anv, o_map_pcld_Anv, o_map_tcld_Anv, &
         o_map_ntot, o_map_hc,o_map_hist,o_map_Cb,o_map_ThCi,o_map_Anv, &
#ifdef ISO
! Isotopes
         o_xtprecip,o_xtplul,o_xtpluc,o_xtovap,o_xtoliq,o_xtcond, &
         o_xtevap,o_dxtdyn,o_dxtldyn,o_dxtcon,o_dxtlsc,o_dxteva, &
         o_dxtajs,o_dxtvdf,o_dxtthe, o_dxtch4, &
         o_dxtprod_nucl,o_dxtcosmo,o_dxtdecroiss, &
         o_xtevap_srf, &
#endif
! Tropopause
         o_alt_tropo, & 
         o_p_tropopause, o_z_tropopause, o_t_tropopause,  &
         o_col_O3_strato, o_col_O3_tropo,                 & 
!--aviation & supersaturation
         o_oclr, o_ocld, o_oss, o_ovc, o_rnebss, o_rnebclr, o_rnebseri, o_gammass, &
         o_N1_ss, o_N2_ss, o_qsatl, o_qsats, &
         o_drnebsub, o_drnebcon, o_drnebtur, o_drnebavi, o_flight_m, o_flight_h2o, &
         o_Tcontr, o_qcontr, o_qcontr2, o_fcontrN, o_fcontrP, &
!--interactive CO2
         o_flx_co2_ocean, o_flx_co2_ocean_cor, &
         o_flx_co2_land, o_flx_co2_land_cor, & 
         o_flx_co2_ff, o_flx_co2_bb, & 
         o_delta_sst, o_delta_sal, o_ds_ns, o_dt_ns, o_dter, o_dser, o_tkt, &
         o_tks, o_taur, o_sss, &
!FC
         o_zxfluxt,o_zxfluxq

#ifdef CPP_ECRAD
    USE phys_output_ctrlout_mod, ONLY:  &
         o_soll0_s2,o_soll_s2,o_sols0_s2,o_sols_s2, &
         o_topl0_s2,o_topl_s2,o_tops0_s2,o_tops_s2   
#endif

#ifdef CPP_StratAer
    USE phys_output_ctrlout_mod, ONLY:  & 
         o_budg_3D_nucl, o_budg_3D_cond_evap, o_budg_3D_ocs_to_so2, o_budg_3D_so2_to_h2so4, &
         o_budg_sed_part, o_R2SO4, o_OCS_lifetime, o_SO2_lifetime, &
         o_budg_3D_backgr_ocs, o_budg_3D_backgr_so2, & 
         o_budg_dep_dry_ocs, o_budg_dep_wet_ocs, &
         o_budg_dep_dry_so2, o_budg_dep_wet_so2, &
         o_budg_dep_dry_h2so4, o_budg_dep_wet_h2so4, &
         o_budg_dep_dry_part, o_budg_dep_wet_part, & 
         o_budg_emi_ocs, o_budg_emi_so2, o_budg_emi_h2so4, o_budg_emi_part, & 
         o_budg_ocs_to_so2, o_budg_so2_to_h2so4, o_budg_h2so4_to_part, &
         o_surf_PM25_sulf, o_ext_strat_550, o_tau_strat_550, &
         o_vsed_aer, o_tau_strat_1020, o_ext_strat_1020, o_f_r_wet
#endif

    USE ice_sursat_mod, ONLY: flight_m, flight_h2o
    
    USE phys_output_ctrlout_mod, ONLY: o_heat_volc, o_cool_volc !NL
    USE phys_state_var_mod, ONLY: heat_volc, cool_volc !NL

    USE phys_state_var_mod, ONLY: pctsrf, rain_fall, snow_fall, bs_fall,&
         qsol, z0m, z0h, fevap, agesno, &
         nday_rain, ndayrain_mth, rain_con, snow_con, &
         topsw, toplw, toplw0, swup, swdn, solswfdiff, &
         topsw0, swupc0, swdnc0, swup0, swdn0, SWup200, SWup200clr, &
         SWdn200, SWdn200clr, LWup200, LWup200clr, &
         LWdn200, LWdn200clr, solsw, solsw0, sollw, &
         radsol, swradcorr, sollw0, sollwdown, sollw, gustiness, &
         sollwdownclr, lwdnc0, lwdn0, ftsol, ustar, u10m, &
         v10m, pbl_tke, wake_delta_pbl_TKE, &
         delta_tsurf, &
         wstar, cape, ema_pcb, ema_pct, &
         ema_cbmf, Mipsh, Ma, fm_therm, ale_bl, alp_bl, ale, &
         alp, cin, wake_pe, wake_dens, cv_gen, wake_s, wake_deltat, &
         wake_deltaq, ftd, fqd, ale_bl_trig, albsol1, &
         ale_wake, ale_bl_stat, &
         rnebcon, wo, falb1, albsol2, coefh, clwcon0, &
         ratqs, entr_therm, zqasc, detr_therm, f0, &
         lwup, lwdn, lwupc0, lwup0, coefm, &
         swupp, lwupp, swupc0p, swup0p, lwupc0p, lwup0p, swdnp, lwdnp, &
         swdnc0p, swdn0p, lwdnc0p, lwdn0p, tnondef, O3sumSTD, uvsumSTD, &
         vqsumSTD, vTsumSTD, O3daysumSTD, wqsumSTD, &
         vphisumSTD, wTsumSTD, u2sumSTD, v2sumSTD, &
         T2sumSTD, nlevSTD, du_gwd_rando, du_gwd_front, &
         ulevSTD, vlevSTD, wlevSTD, philevSTD, qlevSTD, tlevSTD, &
         rhlevSTD, O3STD, O3daySTD, uvSTD, vqSTD, vTSTD, wqSTD, vphiSTD, &
         wTSTD, u2STD, v2STD, T2STD, missing_val_nf90, delta_sal, ds_ns, &
#ifdef ISO
        xtrain_con, xtsnow_con, xtrain_fall, xtsnow_fall, fxtevap, &
#endif
         dt_ns, delta_sst, dter, dser
         
! AI 08 2023 pour ECRAD 3Deffect
#ifdef CPP_ECRAD
    USE phys_state_var_mod, ONLY: &
        sollw0_s2,sollw_s2,solsw0_s2,solsw_s2, &
        toplw0_s2,toplw_s2,topsw0_s2,topsw_s2
#endif


    USE phys_local_var_mod, ONLY: zxfluxlat, slp, ptstar, pt0, zxtsol, zt2m, &
         zn2mout, t2m_min_mon, t2m_max_mon, evap, &
         snowerosion, zxustartlic, zxrhoslic, &
         l_mixmin,l_mix, tke_dissip, &
         zu10m, zv10m, zq2m, zustar, zxqsurf, &
         rain_lsc, rain_num, snow_lsc, bils, sens, fder, &
         zxffonte, zxfqcalving, zxfqfonte, zxrunofflic, fluxu, &
         fluxv, zxsnow, qsnow, snowhgt, to_ice, &
         sissnow, runoff, albsol3_lic, evap_pot, &
         t2m, fluxt, fluxlat, fsollw, fsolsw, &
         wfbils, wfbilo, wfevap, wfrain, wfsnow, &
         cdragm, cdragh, cldl, cldm, &
         cldh, cldt, JrNt,   & ! only output names: cldljn,cldmjn,cldhjn,cldtjn
         cldq, flwp, fiwp, ue, ve, uq, vq, &
         uwat, vwat, &
         rneb_seri, d_rneb_dyn, &
         plcl, plfc, wbeff, convoccur, upwd, dnwd, dnwd0, prw, prlw, prsw,prbsw, &
         s_pblh, s_pblt, s_lcl, s_therm, uwriteSTD, &
         vwriteSTD, wwriteSTD, phiwriteSTD, qwriteSTD, &
         twriteSTD, alp_wake, &
!!         dtvdf_x    ,dtvdf_w    ,dqvdf_x    ,dqvdf_w    , &
         d_t_vdf_x    ,d_t_vdf_w    ,d_q_vdf_x    ,d_q_vdf_w    , &
         sens_x     ,sens_w     ,zxfluxlat_x,zxfluxlat_w, &
         cdragh_x   ,cdragh_w   ,cdragm_x   ,cdragm_w   , &
         kh         ,kh_x       ,kh_w       , &
         wake_h, &
         wake_omg, d_t_wake, d_q_wake, Vprecip, qtaa, Clw, &
         wdtrainA, wdtrainS, wdtrainM, n2, s2, proba_notrig, &
         random_notrig, &
         qclr, qcld, qss, qvc, rnebclr, rnebss, gamma_ss, &
         N1_ss, N2_ss, zqsatl, zqsats, &
         Tcontr, qcontr, qcontr2, fcontrN, fcontrP, &
         drneb_sub, drneb_con, drneb_tur, drneb_avi, &
         alp_bl_det, alp_bl_fluct_m, alp_bl_conv, &
         alp_bl_stat, alp_bl_fluct_tke, slab_wfbils, &
         weak_inversion, dthmin, cldtau, cldemi, &
         pmflxr, pmflxs, prfl, psfl,bsfl, re, fl, rh2m, &
         qsat2m, tpote, tpot, d_ts, od443aer, od550aer, dryod550aer, &
         od865aer, abs550aer, od550lt1aer, sconcso4, sconcno3, &
         sconcoa, sconcbc, sconcss, sconcdust, concso4, concno3, &
         concoa, concbc, concss, concdust, loadso4, &
         loadoa, loadbc, loadss, loaddust, loadno3, tausum_aero, drytausum_aero, &
         topswad_aero, topswad0_aero, solswad_aero, &
         solswad0_aero, topsw_aero, solsw_aero, &
         topsw0_aero, solsw0_aero, topswcf_aero, &
         solswcf_aero, topswai_aero, solswai_aero, &
         toplwad_aero, toplwad0_aero, sollwad_aero, &
         sollwad0_aero, toplwai_aero, sollwai_aero, &
         stratomask,&
         zfice, &
         ec550aer, flwc, fiwc, t_seri, theta, q_seri, &
         ql_seri, qs_seri, qbs_seri, tr_seri, qbs_seri,&
         zphi, u_seri, v_seri, omega, cldfra, &
         rneb, rnebjn, rneblsvol, zx_rh, zx_rhl, zx_rhi, &
         pfraclr, pfracld, d_t_dyn,  & 
         d_q_dyn,  d_ql_dyn, d_qs_dyn, d_qbs_dyn,  &
         d_q_dyn2d,  d_ql_dyn2d, d_qs_dyn2d, d_qbs_dyn2d, &
         d_u_dyn, d_v_dyn, d_t_con, d_t_ajsb, d_t_ajs, &
         d_u_ajs, d_v_ajs, &
         d_u_con, d_v_con, d_q_con, d_q_ajs, d_t_lsc, &
         d_t_lwr,d_t_lw0,d_t_swr,d_t_sw0, &
         d_t_eva, d_q_lsc, beta_prec, d_t_lscth, &
         d_t_lscst, d_q_lscth, d_q_lscst, plul_th, &
         plul_st, d_t_vdf, d_t_diss, d_q_vdf, d_q_eva, &
         d_t_bs, d_q_bs, d_qbs_bs, d_qbs_vdf, &
         zw2, fraca, zmax_th, d_q_ajsb, d_t_ec, d_u_vdf, &
         d_v_vdf, d_u_oro, d_v_oro, d_t_oro, d_u_lif, &
         d_v_lif, d_t_lif, du_gwd_hines, dv_gwd_hines, d_t_hin, &
         dv_gwd_rando, dv_gwd_front, &
         east_gwstress, west_gwstress, &
         d_q_ch4, pmfd, pmfu, ref_liq, ref_ice, rhwriteSTD, &
#ifdef ISO
        xtrain_lsc, xtsnow_lsc, xt_seri, xtl_seri,xts_seri,xtevap, &
        d_xt_dyn,d_xtl_dyn,d_xt_con,d_xt_vdf,d_xt_ajsb, &
        d_xt_lsc,d_xt_eva,d_xt_ch4, &
        d_xt_ajs, d_xt_ajsb, &
        d_xt_prod_nucl,d_xt_cosmo,d_xt_decroiss, &
#endif
         ep, epmax_diag, &  ! epmax_cape
         p_tropopause, t_tropopause, z_tropopause, &
         zxfluxt,zxfluxq, &
! offline 
         da, mp, phi, wght_cvfd
    USE phys_output_var_mod, ONLY: scdnc, cldncl, reffclwtop, lcc, reffclws, &
         reffclwc, cldnvi, lcc3d, lcc3dcon, lcc3dstra, icc3dcon, icc3dstra
    

#ifdef CPP_StratAer 
    USE phys_local_var_mod, ONLY:  &
         budg_3D_nucl, budg_3D_cond_evap, budg_3D_ocs_to_so2, budg_3D_so2_to_h2so4, &
         budg_sed_part, R2SO4, OCS_lifetime, SO2_lifetime, &
         budg_3D_backgr_ocs, budg_3D_backgr_so2, & 
         budg_dep_dry_ocs, budg_dep_wet_ocs, &
         budg_dep_dry_so2, budg_dep_wet_so2, &
         budg_dep_dry_h2so4, budg_dep_wet_h2so4, &
         budg_dep_dry_part, budg_dep_wet_part, &
         budg_emi_ocs, budg_emi_so2, budg_emi_h2so4, budg_emi_part, & 
         budg_ocs_to_so2, budg_so2_to_h2so4, budg_h2so4_to_part, & 
         surf_PM25_sulf, tau_strat_550, tausum_strat, &
         vsed_aer, tau_strat_1020, f_r_wet
#endif

    USE carbon_cycle_mod, ONLY: fco2_ff, fco2_bb, fco2_land, fco2_ocean
    USE carbon_cycle_mod, ONLY: fco2_ocean_cor, fco2_land_cor

    USE phys_output_var_mod, ONLY: vars_defined, snow_o, zfra_o, bils_diss, &
         bils_ec,bils_ech, bils_tke, bils_kinetic, bils_latent, bils_enthalp, &
         itau_con, nfiles, clef_files, nid_files, dryaod_diag, &
         zustr_gwd_hines, zvstr_gwd_hines,zustr_gwd_rando, zvstr_gwd_rando, &
         zustr_gwd_front, zvstr_gwd_front, sza_o,    &
         sens_prec_liq_o, sens_prec_sol_o, lat_prec_liq_o, lat_prec_sol_o, &
         cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv, &
! Marine
         map_prop_hc, map_prop_hist, &
         map_emis_hc,map_iwp_hc,map_deltaz_hc,&
         map_pcld_hc,map_tcld_hc,&
         map_emis_hist,map_iwp_hist,map_deltaz_hist,&
         map_rad_hist,&
         map_ntot,map_hc,map_hist,&
         map_Cb,map_ThCi,map_Anv,&
         map_emis_Cb,map_pcld_Cb,map_tcld_Cb,&
         map_emis_ThCi,map_pcld_ThCi,map_tcld_ThCi,&
         map_emis_Anv,map_pcld_Anv,map_tcld_Anv, &
         alt_tropo, &
!Ionela
         ok_4xCO2atm, tkt, tks, taur, sss

    USE ocean_slab_mod, ONLY: nslay, tslab, slab_bilg, tice, seaice, &
        slab_ekman,slab_hdiff,slab_gm,dt_ekman, dt_hdiff, dt_gm, dt_qflux
    USE pbl_surface_mod, ONLY: snow
    USE indice_sol_mod, ONLY: nbsrf
#ifdef ISO
    USE isotopes_mod, ONLY: iso_HTO
#endif
    USE geometry_mod, ONLY: cell_area, latitude_deg, longitude_deg
    USE surface_data, ONLY: type_ocean, version_ocean, ok_veget, landice_opt
    USE aero_mod, ONLY: naero_tot, id_STRAT_phy
    USE ioipsl, ONLY: histend, histsync
    USE iophy, ONLY: set_itau_iophy, histwrite_phy
    USE netcdf, ONLY: nf90_fill_real
    USE print_control_mod, ONLY: prt_level,lunout
    ! ug Pour les sorties XIOS
    USE lmdz_xios
    USE wxios, ONLY: wxios_closedef, missing_val_xios=>missing_val, wxios_set_context
    USE phys_cal_mod, ONLY : mth_len

#ifdef CPP_RRTM
    USE YOESW, ONLY : RSUN
#endif
    USE tracinca_mod, ONLY: config_inca
    use config_ocean_skin_m, only: activate_ocean_skin

    USE vertical_layers_mod, ONLY: presnivs

    IMPLICIT NONE

    INCLUDE "clesphys.h"
    INCLUDE "alpale.h"
    INCLUDE "compbl.h"
    INCLUDE "YOMCST.h"

    ! Input
    INTEGER :: itap, ivap, iliq, isol, ibs, read_climoz
    INTEGER, DIMENSION(klon) :: lmax_th
    LOGICAL :: aerosol_couple, ok_sync
    LOGICAL :: ok_ade, ok_aie, ok_volcan
    LOGICAL, DIMENSION(klon, klev) :: ptconv, ptconvth
    REAL :: pdtphys
    CHARACTER (LEN=4), DIMENSION(nlevSTD) :: clevSTD
    REAL, DIMENSION(klon,nlevSTD) :: zx_tmp_fi3d_STD
    REAL, DIMENSION(klon) :: pphis
    REAL, DIMENSION(klon, klev) :: pplay, d_u, d_t
    REAL, DIMENSION(klon, klev+1) :: paprs
    REAL, DIMENSION(klon,klev,nqtot) :: qx, d_qx
    REAL, DIMENSION(klon, klev) :: zmasse
    INTEGER :: flag_aerosol_strat
    INTEGER :: flag_aerosol 
    LOGICAL :: ok_cdnc
    REAL, DIMENSION(klon,klev) :: t   ! output for phystoken - offline flux 
    REAL, DIMENSION(klon) :: u1, v1   ! output for phystoken - offline flux
    
    REAL, DIMENSION(3) :: freq_moyNMC

    ! Local
    INTEGER :: itau_w
    INTEGER :: i, iinit, iinitend=1, iff, iq, nsrf, k, ll, naero
    REAL, DIMENSION (klon) :: zx_tmp_fi2d, zpt_conv2d, wind100m
    REAL, DIMENSION (klon,klev) :: zx_tmp_fi3d, zpt_conv
    REAL, DIMENSION (klon,klev+1) :: zx_tmp_fi3d1
    REAL, DIMENSION (klon,NSW) :: zx_tmp_fi3dsp
    CHARACTER (LEN=4)              :: bb2
    INTEGER, DIMENSION(nbp_lon*nbp_lat)  :: ndex2d
    INTEGER, DIMENSION(nbp_lon*nbp_lat*klev) :: ndex3d
    REAL, PARAMETER :: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2
!   REAL, PARAMETER :: missing_val=nf90_fill_real
    REAL, DIMENSION(klev+1,2) :: Ahyb_bounds, Bhyb_bounds
    REAL, DIMENSION(klev,2) :: Ahyb_mid_bounds, Bhyb_mid_bounds
    INTEGER :: ilev
    INTEGER, SAVE :: kmax_100m
!$OMP THREADPRIVATE(kmax_100m)
    REAL :: x
    REAL :: missing_val
    REAL, PARAMETER :: un_jour=86400.
    CHARACTER(len=12) :: nvar    
    INTEGER :: ISW, itr, ixt, it
    CHARACTER*1 ch1
    CHARACTER(LEN=maxlen) :: varname, dn
    REAL, DIMENSION(klon,klev) :: coefh_stok 
    
    
#ifdef CPP_StratAer 
    LOGICAL, PARAMETER :: debug_strataer=.FALSE.
    CHARACTER(LEN=maxlen) :: unt
#endif
    REAL,DIMENSION(klon,klev) :: z, dz
    REAL,DIMENSION(klon)      :: zrho, zt

    ! On calcul le nouveau tau:
    itau_w = itau_phy + itap 
    ! On le donne à iophy pour que les histwrite y aient accès:
    CALL set_itau_iophy(itau_w)

 !   IF (.NOT.vars_defined) THEN
       iinitend = 1
 !   ELSE
 !      iinitend = 1
 !   ENDIF

    IF (using_xios) CALL wxios_set_context

    IF (using_xios) THEN
      missing_val=missing_val_xios
    ELSE
      missing_val=missing_val_nf90
    ENDIF

    IF (.NOT.vars_defined) THEN
      kmax_100m=1
      DO k=1, klev-1
        IF (presnivs(k).GT.0.97*101325.) kmax_100m = k !--finding out max level for 100 m with a good margin
      ENDDO
    ENDIF

    Ahyb_bounds(1,1) = 0.
    Ahyb_bounds(1,2) = aps(1)
    Bhyb_bounds(1,1) = 1.
    Bhyb_bounds(1,2) = bps(1)    
    DO ilev=2,klev
      Ahyb_bounds(ilev,1) = aps(ilev-1)
      Ahyb_bounds(ilev,2) = aps(ilev)
      Bhyb_bounds(ilev,1) = bps(ilev-1)
      Bhyb_bounds(ilev,2) = bps(ilev)
    ENDDO
     Ahyb_bounds(klev+1,1) = aps(klev)
     Ahyb_bounds(klev+1,2) = 0.
     Bhyb_bounds(klev+1,1) = bps(klev)
     Bhyb_bounds(klev+1,2) = 0.

    DO ilev=1, klev
      Ahyb_mid_bounds(ilev,1) = ap(ilev)
      Ahyb_mid_bounds(ilev,2) = ap(ilev+1)
      Bhyb_mid_bounds(ilev,1) = bp(ilev)
      Bhyb_mid_bounds(ilev,2) = bp(ilev+1)
    END DO


    ! ug la boucle qui suit ne sert qu'une fois, pour l'initialisation, sinon il n'y a toujours qu'un seul passage:
    DO iinit=1, iinitend
!      print *,'IFF iinit=', iinit, iinitend 
       IF (using_xios) THEN
         !$OMP MASTER
         IF (vars_defined) THEN
            IF (prt_level >= 10) then
               write(lunout,*)"phys_output_write: call xios_update_calendar, itau_w=",itau_w
            ENDIF
!            CALL xios_update_calendar(itau_w)
            CALL xios_update_calendar(itap)
         ENDIF
         !$OMP END MASTER
         !$OMP BARRIER
       ENDIF

       ! On procède à l'écriture ou à la définition des nombreuses variables:
!!! Champs 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL histwrite_phy(o_phis, pphis)

       zx_tmp_fi2d = cell_area
       IF (is_north_pole_phy) then
         zx_tmp_fi2d(1) = cell_area(1)/nbp_lon
       ENDIF
       IF (is_south_pole_phy) then
         zx_tmp_fi2d(klon) = cell_area(klon)/nbp_lon
       ENDIf
       CALL histwrite_phy(o_aire, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=pctsrf(i,is_ter)+pctsrf(i,is_lic)
          ENDDO
       ENDIF

       CALL histwrite_phy(o_contfracATM, zx_tmp_fi2d)
       CALL histwrite_phy(o_contfracOR, pctsrf(:,is_ter))
!
       IF (using_xios) THEN

         CALL histwrite_phy("R_ecc",R_ecc)
         CALL histwrite_phy("R_peri",R_peri)
         CALL histwrite_phy("R_incl",R_incl)
         CALL histwrite_phy("solaire",solaire)
         CALL histwrite_phy(o_Ahyb, ap)
         CALL histwrite_phy(o_Bhyb, bp)
         CALL histwrite_phy(o_Ahyb_bounds, Ahyb_bounds)
         CALL histwrite_phy(o_Bhyb_bounds, Bhyb_bounds)
         CALL histwrite_phy(o_Ahyb_mid, aps)
         CALL histwrite_phy(o_Bhyb_mid, bps)
         CALL histwrite_phy(o_Ahyb_mid_bounds, Ahyb_mid_bounds)
         CALL histwrite_phy(o_Bhyb_mid_bounds, Bhyb_mid_bounds)
         CALL histwrite_phy(o_longitude, longitude_deg)
         CALL histwrite_phy(o_latitude, latitude_deg)
!
#ifdef CPP_RRTM
        IF (iflag_rrtm.EQ.1) THEN
          DO ISW=1, NSW
            WRITE(ch1,'(i1)') ISW
  !         zx_tmp_0d=RSUN(ISW)
  !         CALL histwrite_phy("rsun"//ch1,zx_tmp_0d)
            CALL histwrite_phy("rsun"//ch1,RSUN(ISW))
          ENDDO
        ENDIF
#endif
!
        CALL histwrite_phy("co2_ppm",co2_ppm)
        CALL histwrite_phy("CH4_ppb",CH4_ppb)
        CALL histwrite_phy("N2O_ppb",N2O_ppb)
        CALL histwrite_phy("CFC11_ppt",CFC11_ppt)
        CALL histwrite_phy("CFC12_ppt",CFC12_ppt)
!
      ENDIF !using_xios

!!! Champs 2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulateur AIRS
     IF (ok_airs) then
      CALL histwrite_phy(o_alt_tropo,alt_tropo)
 
      CALL histwrite_phy(o_map_prop_hc,map_prop_hc)
      CALL histwrite_phy(o_map_prop_hist,map_prop_hist)

      CALL histwrite_phy(o_map_emis_hc,map_emis_hc)
      CALL histwrite_phy(o_map_iwp_hc,map_iwp_hc)
      CALL histwrite_phy(o_map_deltaz_hc,map_deltaz_hc)
      CALL histwrite_phy(o_map_pcld_hc,map_pcld_hc)
      CALL histwrite_phy(o_map_tcld_hc,map_tcld_hc)

      CALL histwrite_phy(o_map_emis_hist,map_emis_hist)
      CALL histwrite_phy(o_map_iwp_hist,map_iwp_hist)
      CALL histwrite_phy(o_map_deltaz_hist,map_deltaz_hist)

      CALL histwrite_phy(o_map_ntot,map_ntot)
      CALL histwrite_phy(o_map_hc,map_hc)
      CALL histwrite_phy(o_map_hist,map_hist)

      CALL histwrite_phy(o_map_Cb,map_Cb)
      CALL histwrite_phy(o_map_ThCi,map_ThCi)
      CALL histwrite_phy(o_map_Anv,map_Anv)

      CALL histwrite_phy(o_map_emis_Cb,map_emis_Cb)
      CALL histwrite_phy(o_map_pcld_Cb,map_pcld_Cb)
      CALL histwrite_phy(o_map_tcld_Cb,map_tcld_Cb)

      CALL histwrite_phy(o_map_emis_ThCi,map_emis_ThCi)
      CALL histwrite_phy(o_map_pcld_ThCi,map_pcld_ThCi)
      CALL histwrite_phy(o_map_tcld_ThCi,map_tcld_ThCi)

      CALL histwrite_phy(o_map_emis_Anv,map_emis_Anv)
      CALL histwrite_phy(o_map_pcld_Anv,map_pcld_Anv)
      CALL histwrite_phy(o_map_tcld_Anv,map_tcld_Anv)
     ENDIF

       CALL histwrite_phy(o_sza, sza_o)
       CALL histwrite_phy(o_flat, zxfluxlat)
       CALL histwrite_phy(o_ptstar, ptstar)
       CALL histwrite_phy(o_pt0, pt0)
       CALL histwrite_phy(o_slp, slp)
       CALL histwrite_phy(o_tsol, zxtsol)
       CALL histwrite_phy(o_t2m, zt2m)
       CALL histwrite_phy(o_t2m_min, zt2m)
       CALL histwrite_phy(o_t2m_max, zt2m)
       CALL histwrite_phy(o_t2m_max_mon, t2m_max_mon)
       CALL histwrite_phy(o_t2m_min_mon, t2m_min_mon)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,1))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nt2mout, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,2))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nt2moutfg, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,3))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nq2mout, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,4))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nq2moutfg, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,5))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nu2mout, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=real(zn2mout(i,6))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_nu2moutfg, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_wind10m, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i=1, klon
             zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
          ENDDO
       ENDIF
       CALL histwrite_phy(o_wind10max, zx_tmp_fi2d)

       CALL histwrite_phy(o_gusts, gustiness)

       IF (vars_defined) THEN
          DO k = 1, kmax_100m                                      !--we could stop much lower
            zrho(:) = pplay(:,k)/t_seri(:,k)/RD                    ! air density in kg/m3
            dz(:,k) = (paprs(:,k)-paprs(:,k+1))/zrho(:)/RG         ! layer thickness in m
            IF (k==1) THEN
              z(:,1) = (paprs(:,1)-pplay(:,1))/zrho(:)/RG          ! altitude middle of first layer in m
              zt(:)  = dz(:,1)                                     ! altitude top of first layer in m
            ELSE
              z(:,k) = zt(:) + (paprs(:,k)-pplay(:,k))/zrho(:)/RG  ! altitude middle of layer k in m
              zt(:)  = zt(:) + dz(:,k)                             ! altitude top of layer k in m
            ENDIF
          ENDDO
          wind100m(:)=missing_val
          DO k=1, kmax_100m-1                                      !--we could stop much lower
            DO i=1,klon
              IF (z(i,k).LT.100..AND.z(i,k+1).GE.100.) THEN
                wind100m(i)=SQRT( (u_seri(i,k)+(100.-z(i,k))/(z(i,k+1)-z(i,k))*(u_seri(i,k+1)-u_seri(i,k)))**2.0 + &
                                  (v_seri(i,k)+(100.-z(i,k))/(z(i,k+1)-z(i,k))*(v_seri(i,k+1)-v_seri(i,k)))**2.0 )
              ENDIF
            ENDDO
          ENDDO
       ENDIF
       CALL histwrite_phy(o_wind100m, wind100m)

       IF (vars_defined) THEN
         !--polynomial fit for 14,Vestas,1074,V136/3450 kW windmill - Olivier
         DO i=1,klon
           IF (pctsrf(i,is_ter).GT.0.05 .AND. wind100m(i).NE.missing_val) THEN
             x=wind100m(i)
             IF (x.LE.3.0 .OR. x.GE.22.5) THEN
               zx_tmp_fi2d(i)=0.0
             ELSE IF (x.GE.10.0) THEN
               zx_tmp_fi2d(i)=1.0
             ELSE
               zx_tmp_fi2d(i)= 10.73 + x*(-14.69 + x*(8.339 + x*(-2.59 + x*(0.4893 + x*(-0.05898 + x*(0.004627 + & 
                               x*(-0.0002352 + x*(7.478e-06 + x*(-1.351e-07 + x*(1.059e-09))))))))))
               zx_tmp_fi2d(i)=MIN(MAX(zx_tmp_fi2d(i),0.0),1.0)
             ENDIF
           ELSE
             zx_tmp_fi2d(i)=missing_val
           ENDIF
         ENDDO
       ENDIF
       CALL histwrite_phy(o_loadfactor_wind_onshore, zx_tmp_fi2d)

       IF (vars_defined) THEN
         !--polynomial fit for 14,Vestas,867,V164/8000 kW - Olivier
         DO i=1,klon
           IF (pctsrf(i,is_oce).GT.0.05 .AND. wind100m(i).NE.missing_val) THEN
             x=wind100m(i)
             IF (x.LE.3.0 .OR. x.GE.25.5) THEN
               zx_tmp_fi2d(i)=0.0
             ELSE IF (x.GE.12.5) THEN
               zx_tmp_fi2d(i)=1.0
             ELSE
               zx_tmp_fi2d(i)= 20.59 + x*(-22.39 + x*(10.25 + x*(-2.601 + x*(0.4065 + x*(-0.04099 + x*(0.002716 + & 
                               x*(-0.0001175 + x*(3.195e-06 + x*(-4.959e-08 + x*(3.352e-10))))))))))
               zx_tmp_fi2d(i)=MIN(MAX(zx_tmp_fi2d(i),0.0),1.0)
             ENDIF
           ELSE
             zx_tmp_fi2d(i)=missing_val
           ENDIF
         ENDDO
       ENDIF
       CALL histwrite_phy(o_loadfactor_wind_offshore, zx_tmp_fi2d)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = pctsrf(i,is_sic)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_sicf, zx_tmp_fi2d)
       CALL histwrite_phy(o_q2m, zq2m)
       IF (vars_defined) zx_tmp_fi2d = zustar
       CALL histwrite_phy(o_ustar, zx_tmp_fi2d)
       CALL histwrite_phy(o_u10m, zu10m)
       CALL histwrite_phy(o_v10m, zv10m)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = paprs(i,1)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_psol, zx_tmp_fi2d)
       CALL histwrite_phy(o_mass, zmasse)
       CALL histwrite_phy(o_qsurf, zxqsurf)

       IF (.NOT. ok_veget) THEN
          CALL histwrite_phy(o_qsol, qsol)
       ENDIF

       IF (vars_defined) THEN
          IF (ok_bs) THEN
             DO i = 1, klon
             zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i) + bs_fall(i)
             ENDDO
          ELSE
             DO i = 1, klon
             zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i)
             ENDDO
          ENDIF
       ENDIF

       CALL histwrite_phy(o_precip, zx_tmp_fi2d)
       CALL histwrite_phy(o_rain_fall, rain_fall)
       CALL histwrite_phy(o_ndayrain, ndayrain_mth)

       ! epmax_cape:
!       CALL histwrite_phy(o_epmax_diag, epmax_diag)
       CALL histwrite_phy(o_ep, ep)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = rain_lsc(i) + snow_lsc(i)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_plul, zx_tmp_fi2d)
       CALL histwrite_phy(o_plun, rain_num)

       IF (vars_defined) THEN
          DO i = 1, klon
             zx_tmp_fi2d(i) = rain_con(i) + snow_con(i)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_pluc, zx_tmp_fi2d)
       CALL histwrite_phy(o_rain_con, rain_con)
       CALL histwrite_phy(o_snow, snow_fall)
       CALL histwrite_phy(o_msnow, zxsnow)
       CALL histwrite_phy(o_fsnow, zfra_o)
       CALL histwrite_phy(o_evap, evap)

       IF (ok_bs) THEN
           CALL histwrite_phy(o_bsfall, bs_fall)     
           CALL histwrite_phy(o_snowerosion, snowerosion)
           CALL histwrite_phy(o_ustart_lic, zxustartlic)
           CALL histwrite_phy(o_rhosnow_lic, zxrhoslic)
       ENDIF

       IF (vars_defined) THEN
         zx_tmp_fi2d = topsw*swradcorr
       ENDIF
       CALL histwrite_phy(o_tops, zx_tmp_fi2d)

       IF (vars_defined) THEN
         zx_tmp_fi2d = topsw0*swradcorr
       ENDIF
       CALL histwrite_phy(o_tops0, zx_tmp_fi2d)

       CALL histwrite_phy(o_topl, toplw)
       CALL histwrite_phy(o_topl0, toplw0)

! offline 
       IF (using_xios) THEN
         IF (offline) THEN

            coefh_stok(:,1)      = cdragh(:)
            coefh_stok(:,2:klev) = coefh(:,2:klev, is_ave)
          
            CALL histwrite_phy('upwd_stok', upwd)
            CALL histwrite_phy('t_stok', t)
            CALL histwrite_phy('fm_th_stok', fm_therm(:,1:klev))
            CALL histwrite_phy('en_th_stok', entr_therm)
            CALL histwrite_phy('da_stok',da )
            CALL histwrite_phy('mp_stok',mp )
            CALL histwrite_phy('dnwd_stok', dnwd)
            CALL histwrite_phy('wght_stok', wght_cvfd)
            CALL histwrite_phy('coefh_stok', coefh_stok)
            CALL histwrite_phy('yu1_stok', u1) 
            CALL histwrite_phy('yv1_stok', v1)

            DO k=1,klev
               IF (k<10) THEN
                  WRITE(nvar,'(i1)') k
               ELSE IF (k<100) THEN
                  WRITE(nvar,'(i2)') k
               ELSE
                  WRITE(nvar,'(i3)') k
               END IF
               nvar='phi_lev'//trim(nvar)
               CALL histwrite_phy(nvar,phi(:,:,k))
            END DO
          
         ENDIF
       ENDIF



       
       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swup(:,klevp1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWupTOA, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swup0(:,klevp1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWupTOAclr, zx_tmp_fi2d)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             zx_tmp_fi2d(:) = swupc0(:,klevp1)*swradcorr(:)
          ENDIF
          CALL histwrite_phy(o_SWupTOAcleanclr, zx_tmp_fi2d)
       ENDIF

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swdn(:,klevp1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWdnTOA, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swdn0(:,klevp1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWdnTOAclr, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = topsw(:)*swradcorr(:)-toplw(:)
       ENDIF
       CALL histwrite_phy(o_nettop, zx_tmp_fi2d)
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = SWup200*swradcorr
       ENDIF
       CALL histwrite_phy(o_SWup200, zx_tmp_fi2d)
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = SWup200clr*swradcorr
       ENDIF
       CALL histwrite_phy(o_SWup200clr, zx_tmp_fi2d)
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = SWdn200*swradcorr
       ENDIF
       CALL histwrite_phy(o_SWdn200, zx_tmp_fi2d)
       
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = SWdn200clr*swradcorr
       ENDIF
       CALL histwrite_phy(o_SWdn200clr, zx_tmp_fi2d)
       
       CALL histwrite_phy(o_LWup200, LWup200)
       CALL histwrite_phy(o_LWup200clr, LWup200clr)
       CALL histwrite_phy(o_LWdn200, LWdn200)
       CALL histwrite_phy(o_LWdn200clr, LWdn200clr)
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = solsw*swradcorr
       ENDIF
       CALL histwrite_phy(o_sols, zx_tmp_fi2d)
       
       
       IF (vars_defined) THEN
          zx_tmp_fi2d = solsw0*swradcorr
       ENDIF
       CALL histwrite_phy(o_sols0, zx_tmp_fi2d)
       CALL histwrite_phy(o_soll, sollw)
       CALL histwrite_phy(o_soll0, sollw0)
       CALL histwrite_phy(o_radsol, radsol)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swup(:,1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWupSFC, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swup0(:,1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWupSFCclr, zx_tmp_fi2d)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             zx_tmp_fi2d(:) = swupc0(:,1)*swradcorr(:)
          ENDIF
          CALL histwrite_phy(o_SWupSFCcleanclr, zx_tmp_fi2d)
       ENDIF

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swdn(:,1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWdnSFC, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = swdn0(:,1)*swradcorr(:)
       ENDIF
       CALL histwrite_phy(o_SWdnSFCclr, zx_tmp_fi2d)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             zx_tmp_fi2d(:) = swdnc0(:,1)*swradcorr(:)
          ENDIF
          CALL histwrite_phy(o_SWdnSFCcleanclr, zx_tmp_fi2d)
       ENDIF

       CALL histwrite_phy(o_fdiffSWdnSFC, solswfdiff)

       IF (vars_defined) THEN
          zx_tmp_fi2d(:)=sollwdown(:)-sollw(:)
       ENDIF
       CALL histwrite_phy(o_LWupSFC, zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdnSFC, sollwdown)

       IF (vars_defined) THEN
          sollwdownclr(1:klon) = -1.*lwdn0(1:klon,1)
          zx_tmp_fi2d(1:klon)=sollwdownclr(1:klon)-sollw0(1:klon)
       ENDIF
       CALL histwrite_phy(o_LWupSFCclr, zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdnSFCclr, sollwdownclr)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             zx_tmp_fi2d(:) = lwupc0(:,klevp1)
          ENDIF
          CALL histwrite_phy(o_LWupTOAcleanclr, zx_tmp_fi2d)
       ENDIF

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             zx_tmp_fi2d(:) = -1.*lwdnc0(:,1)
          ENDIF
          CALL histwrite_phy(o_LWdnSFCcleanclr, zx_tmp_fi2d)
       ENDIF

!AI 08 2023 Ecrad 3Deffect
#ifdef CPP_ECRAD
     if (ok_3Deffect) then
        IF (vars_defined) THEN
          zx_tmp_fi2d = solsw_s2*swradcorr
       ENDIF
       CALL histwrite_phy(o_sols_s2, zx_tmp_fi2d)
       IF (vars_defined) THEN
          zx_tmp_fi2d = solsw0_s2*swradcorr
       ENDIF
       CALL histwrite_phy(o_sols0_s2, zx_tmp_fi2d)
       CALL histwrite_phy(o_soll_s2, sollw_s2)
       CALL histwrite_phy(o_soll0_s2, sollw0_s2)
       IF (vars_defined) THEN
         zx_tmp_fi2d = topsw_s2*swradcorr
       ENDIF
       CALL histwrite_phy(o_tops_s2, zx_tmp_fi2d)

       IF (vars_defined) THEN
         zx_tmp_fi2d = topsw0_s2*swradcorr
       ENDIF
       CALL histwrite_phy(o_tops0_s2, zx_tmp_fi2d)

       CALL histwrite_phy(o_topl_s2, toplw_s2)
       CALL histwrite_phy(o_topl0_s2, toplw0_s2)
     endif
#endif       

       CALL histwrite_phy(o_bils, bils)
       CALL histwrite_phy(o_bils_diss, bils_diss)
       CALL histwrite_phy(o_bils_ec, bils_ec)
       CALL histwrite_phy(o_bils_ech, bils_ech)
       CALL histwrite_phy(o_bils_tke, bils_tke)
       CALL histwrite_phy(o_bils_kinetic, bils_kinetic)
       CALL histwrite_phy(o_bils_latent, bils_latent)
       CALL histwrite_phy(o_bils_enthalp, bils_enthalp)

       IF (vars_defined) THEN
          zx_tmp_fi2d(1:klon)=-1*sens(1:klon)
       ENDIF
       CALL histwrite_phy(o_sens, zx_tmp_fi2d)
       CALL histwrite_phy(o_fder, fder)
       CALL histwrite_phy(o_ffonte, zxffonte)
       CALL histwrite_phy(o_fqcalving, zxfqcalving)
       CALL histwrite_phy(o_fqfonte, zxfqfonte)
       IF (vars_defined) THEN
          zx_tmp_fi2d(1:klon)=(zxfqfonte(1:klon)+rain_fall(1:klon))*pctsrf(1:klon,is_lic)
       ENDIF
       CALL histwrite_phy(o_mrroli, zx_tmp_fi2d)
       CALL histwrite_phy(o_runofflic, zxrunofflic)
       IF (vars_defined) THEN
          zx_tmp_fi2d=0.
          DO nsrf=1,nbsrf
             zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxu(:,1,nsrf)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_taux, zx_tmp_fi2d)

       IF (vars_defined) THEN
          zx_tmp_fi2d=0.
          DO nsrf=1,nbsrf
             zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxv(:,1,nsrf)
          ENDDO
       ENDIF
       CALL histwrite_phy(o_tauy, zx_tmp_fi2d)

       DO nsrf = 1, nbsrf

          IF (vars_defined)             zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)*100.
          CALL histwrite_phy(o_pourc_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)           zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)
          CALL histwrite_phy(o_fract_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxu( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_taux_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxv( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_tauy_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = ftsol( 1 : klon, nsrf)
          CALL histwrite_phy(o_tsol_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = evap_pot( 1 : klon, nsrf)
          CALL histwrite_phy(o_evappot_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = ustar(1 : klon, nsrf)
          CALL histwrite_phy(o_ustar_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = u10m(1 : klon, nsrf)
          CALL histwrite_phy(o_u10m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = v10m(1 : klon, nsrf)
          CALL histwrite_phy(o_v10m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = t2m(1 : klon, nsrf)
          CALL histwrite_phy(o_t2m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)       zx_tmp_fi2d(1 : klon) = fevap(1 : klon, nsrf)
          CALL histwrite_phy(o_evap_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)        zx_tmp_fi2d(1 : klon) = fluxt( 1 : klon, 1, nsrf)
          CALL histwrite_phy(o_sens_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fluxlat( 1 : klon, nsrf)
          CALL histwrite_phy(o_lat_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fsollw( 1 : klon, nsrf)
          CALL histwrite_phy(o_flw_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, nsrf)
          CALL histwrite_phy(o_fsw_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfbils( 1 : klon, nsrf)
          CALL histwrite_phy(o_wbils_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfbilo( 1 : klon, nsrf)
          CALL histwrite_phy(o_wbilo_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfevap( 1 : klon, nsrf)
          CALL histwrite_phy(o_wevap_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfrain( 1 : klon, nsrf)
          CALL histwrite_phy(o_wrain_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined)         zx_tmp_fi2d(1 : klon) = wfsnow( 1 : klon, nsrf)
          CALL histwrite_phy(o_wsnow_srf(nsrf), zx_tmp_fi2d)

          IF (iflag_pbl > 1) THEN
             CALL histwrite_phy(o_tke_srf(nsrf),  pbl_tke(:,1:klev,nsrf))
             !CALL histwrite_phy(o_l_mix(nsrf),  l_mix(:,1:klev,nsrf))
             CALL histwrite_phy(o_l_mixmin(nsrf),  l_mixmin(:,1:klev,nsrf))
             CALL histwrite_phy(o_tke_max_srf(nsrf),  pbl_tke(:,1:klev,nsrf))

                         
          ENDIF
!jyg<
          IF (iflag_pbl > 1 .AND. iflag_wake>=1  .AND. iflag_pbl_split >=1) THEN
             CALL histwrite_phy(o_dltpbltke_srf(nsrf), wake_delta_pbl_TKE(:,1:klev,nsrf))
          ENDIF
!>jyg
!          IF (iflag_pbl > 1 .AND. ifl_pbltree  >=1 ) THEN
!       CALL histwrite_phy(o_treedrg_srf(nsrf), treedrg(:,1:klev,nsrf))
!            ENDIF

       ENDDO
       
                
        IF (iflag_pbl > 1) THEN
          zx_tmp_fi3d=0.
          IF (vars_defined) THEN
             DO nsrf=1,nbsrf
                DO k=1,klev
                   zx_tmp_fi3d(:,k)=zx_tmp_fi3d(:,k) &
                        +pctsrf(:,nsrf)*tke_dissip(:,k,nsrf)
                ENDDO
             ENDDO
          ENDIF
          
          CALL histwrite_phy(o_tke_dissip, zx_tmp_fi3d)    
       ENDIF

       IF (vars_defined) zx_tmp_fi2d(1 : klon) = sens_prec_liq_o(1 : klon, 1)
       CALL histwrite_phy(o_sens_prec_liq_oce, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = sens_prec_liq_o(1 : klon, 2)
       CALL histwrite_phy(o_sens_prec_liq_sic, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = sens_prec_sol_o(1 : klon, 1)
       CALL histwrite_phy(o_sens_prec_sol_oce, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = sens_prec_sol_o(1 : klon, 2)
       CALL histwrite_phy(o_sens_prec_sol_sic, zx_tmp_fi2d)       

       IF (vars_defined) zx_tmp_fi2d(1 : klon) = lat_prec_liq_o(1 : klon, 1)
       CALL histwrite_phy(o_lat_prec_liq_oce, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = lat_prec_liq_o(1 : klon, 2)
       CALL histwrite_phy(o_lat_prec_liq_sic, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = lat_prec_sol_o(1 : klon, 1)
       CALL histwrite_phy(o_lat_prec_sol_oce, zx_tmp_fi2d)       
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = lat_prec_sol_o(1 : klon, 2)
       CALL histwrite_phy(o_lat_prec_sol_sic, zx_tmp_fi2d)       

       DO nsrf=1,nbsrf+1
          CALL histwrite_phy(o_wstar(nsrf), wstar(1 : klon, nsrf))
       ENDDO

       CALL histwrite_phy(o_cdrm, cdragm)
       CALL histwrite_phy(o_cdrh, cdragh)
       CALL histwrite_phy(o_cldl, cldl)
       CALL histwrite_phy(o_cldm, cldm)
       CALL histwrite_phy(o_cldh, cldh)
       CALL histwrite_phy(o_cldt, cldt)
       CALL histwrite_phy(o_JrNt, JrNt)
       
       IF (vars_defined)  zx_tmp_fi2d=cldl*JrNt      
       CALL histwrite_phy(o_cldljn, zx_tmp_fi2d)
       
       IF (vars_defined)  zx_tmp_fi2d=cldm*JrNt      
       CALL histwrite_phy(o_cldmjn, zx_tmp_fi2d)
       
       IF (vars_defined)  zx_tmp_fi2d=cldh*JrNt
       CALL histwrite_phy(o_cldhjn, zx_tmp_fi2d)
       
       IF (vars_defined)  zx_tmp_fi2d=cldt*JrNt
       CALL histwrite_phy(o_cldtjn, zx_tmp_fi2d)
       
       CALL histwrite_phy(o_cldq, cldq)
       IF (vars_defined)       zx_tmp_fi2d(1:klon) = flwp(1:klon)
       CALL histwrite_phy(o_lwp, zx_tmp_fi2d)
       IF (vars_defined)       zx_tmp_fi2d(1:klon) = fiwp(1:klon)
       CALL histwrite_phy(o_iwp, zx_tmp_fi2d)
       CALL histwrite_phy(o_ue, ue)
       CALL histwrite_phy(o_ve, ve)
       CALL histwrite_phy(o_uq, uq)
       CALL histwrite_phy(o_vq, vq)
       CALL histwrite_phy(o_uwat, uwat)
       CALL histwrite_phy(o_vwat, vwat)
       IF (iflag_con.GE.3) THEN ! sb
          CALL histwrite_phy(o_cape, cape)
          CALL histwrite_phy(o_pbase, ema_pcb)
          CALL histwrite_phy(o_ptop, ema_pct)
          CALL histwrite_phy(o_fbase, ema_cbmf)
          IF (iflag_con /= 30) THEN
             CALL histwrite_phy(o_plcl, plcl)
             CALL histwrite_phy(o_plfc, plfc)
             CALL histwrite_phy(o_wbeff, wbeff)
             CALL histwrite_phy(o_convoccur, convoccur)
          ENDIF

          CALL histwrite_phy(o_cape_max, cape)

          CALL histwrite_phy(o_upwd, upwd)
          CALL histwrite_phy(o_Ma, Ma)
          CALL histwrite_phy(o_dnwd, dnwd)
          CALL histwrite_phy(o_dnwd0, dnwd0)
          !! The part relative to the frequency of occurence of convection
          !! is now grouped with the part relative to thermals and shallow
          !! convection (output of the 3 fields: ftime_deepcv, ftime_th and
          !!  ftime_con).
          IF (vars_defined) THEN
             IF (iflag_thermals>=1)THEN
                zx_tmp_fi3d=dnwd+dnwd0+upwd+fm_therm(:,1:klev)
             ELSE
                zx_tmp_fi3d=dnwd+dnwd0+upwd
             ENDIF
          ENDIF
          CALL histwrite_phy(o_mc, zx_tmp_fi3d)
       ENDIF !iflag_con .GE. 3
       CALL histwrite_phy(o_prw, prw)
       CALL histwrite_phy(o_prlw, prlw)
       CALL histwrite_phy(o_prsw, prsw)
       IF (ok_bs) THEN
       CALL histwrite_phy(o_prbsw, prbsw)
       ENDIF
       CALL histwrite_phy(o_s_pblh, s_pblh)
       CALL histwrite_phy(o_s_pblt, s_pblt)
       CALL histwrite_phy(o_s_lcl, s_lcl)
       CALL histwrite_phy(o_s_therm, s_therm)
       !IM : Les champs suivants (s_capCL, s_oliqCL, s_cteiCL, s_trmb1, s_trmb2, s_trmb3) ne sont pas definis dans HBTM.F
       !       IF (o_s_capCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_capCL%name,itau_w,s_capCL)
       !       ENDIF
       !       IF (o_s_oliqCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_oliqCL%name,itau_w,s_oliqCL)
       !       ENDIF
       !       IF (o_s_cteiCL%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_cteiCL%name,itau_w,s_cteiCL)
       !       ENDIF
       !       IF (o_s_trmb1%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb1%name,itau_w,s_trmb1)
       !       ENDIF
       !       IF (o_s_trmb2%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb2%name,itau_w,s_trmb2)
       !       ENDIF
       !       IF (o_s_trmb3%flag(iff)<=lev_files(iff)) THEN
       !     CALL histwrite_phy(nid_files(iff),clef_stations(iff),
       !    $o_s_trmb3%name,itau_w,s_trmb3)
       !       ENDIF

#ifdef CPP_IOIPSL
       IF (.NOT. using_xios) THEN
         IF (.NOT.ok_all_xml) THEN
           ! ATTENTION, LES ANCIENS HISTWRITE ONT ETES CONSERVES EN ATTENDANT MIEUX:
           ! Champs interpolles sur des niveaux de pression
            DO iff=1, nfiles
              ll=0
              DO k=1, nlevSTD
                bb2=clevSTD(k) 
                  IF (bb2.EQ."850".OR.bb2.EQ."700".OR. &
                       bb2.EQ."500".OR.bb2.EQ."200".OR. &
                       bb2.EQ."100".OR. &
                       bb2.EQ."50".OR.bb2.EQ."10") THEN

                      ! a refaire correctement !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ll=ll+1
                      CALL histwrite_phy(o_uSTDlevs(ll),uwriteSTD(:,k,iff), iff)
                      CALL histwrite_phy(o_vSTDlevs(ll),vwriteSTD(:,k,iff), iff)
                      CALL histwrite_phy(o_wSTDlevs(ll),wwriteSTD(:,k,iff), iff)
                      CALL histwrite_phy(o_zSTDlevs(ll),phiwriteSTD(:,k,iff), iff)
                      CALL histwrite_phy(o_qSTDlevs(ll),qwriteSTD(:,k,iff), iff)
                      CALL histwrite_phy(o_tSTDlevs(ll),twriteSTD(:,k,iff), iff)

                  ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR.
              ENDDO
            ENDDO
         ENDIF
       ENDIF
#endif

       IF (using_xios) THEN
         IF (ok_all_xml) THEN
           !XIOS  CALL xios_get_field_attr("u850",default_value=missing_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ll=0
            DO k=1, nlevSTD
              bb2=clevSTD(k) 
              IF (bb2.EQ."850".OR.bb2.EQ."700".OR. &
                  bb2.EQ."500".OR.bb2.EQ."200".OR. &
                  bb2.EQ."100".OR. &
                  bb2.EQ."50".OR.bb2.EQ."10") THEN
                  ll=ll+1
                  CALL histwrite_phy(o_uSTDlevs(ll),ulevSTD(:,k))
                  CALL histwrite_phy(o_vSTDlevs(ll),vlevSTD(:,k))
                  CALL histwrite_phy(o_wSTDlevs(ll),wlevSTD(:,k))
                  CALL histwrite_phy(o_zSTDlevs(ll),philevSTD(:,k))
                  CALL histwrite_phy(o_qSTDlevs(ll),qlevSTD(:,k))
                  CALL histwrite_phy(o_tSTDlevs(ll),tlevSTD(:,k))
              ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR.
            ENDDO
         ENDIF
       ENDIF

       IF (vars_defined) THEN
          DO i=1, klon
             IF (pctsrf(i,is_oce).GT.epsfra.OR. &
                  pctsrf(i,is_sic).GT.epsfra) THEN
                zx_tmp_fi2d(i) = (ftsol(i, is_oce) * pctsrf(i,is_oce)+ &
                     ftsol(i, is_sic) * pctsrf(i,is_sic))/ &
                     (pctsrf(i,is_oce)+pctsrf(i,is_sic))
             ELSE
                zx_tmp_fi2d(i) = 273.15
             ENDIF
          ENDDO
       ENDIF
       CALL histwrite_phy(o_t_oce_sic, zx_tmp_fi2d)

       ! Couplage convection-couche limite
       IF (iflag_con.GE.3) THEN
          IF (iflag_coupl>=1) THEN
             CALL histwrite_phy(o_ale_bl, ale_bl)
             CALL histwrite_phy(o_alp_bl, alp_bl)
          ENDIF !iflag_coupl>=1
       ENDIF !(iflag_con.GE.3)
       ! Wakes
       IF (iflag_con.EQ.3) THEN
          CALL histwrite_phy(o_Mipsh, Mipsh)
          IF (iflag_wake>=1) THEN
             CALL histwrite_phy(o_ale_wk, ale_wake)
             CALL histwrite_phy(o_alp_wk, alp_wake)
             IF (iflag_pbl_split>=1) THEN
!!               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=dtvdf_x(1:klon,1:klev)/pdtphys
!!               CALL histwrite_phy(o_dtvdf_x    ,zx_tmp_fi3d) 
!!               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=dtvdf_w(1:klon,1:klev)/pdtphys
!!               CALL histwrite_phy(o_dtvdf_w    ,zx_tmp_fi3d) 
!!               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=dqvdf_x(1:klon,1:klev)/pdtphys
!!               CALL histwrite_phy(o_dqvdf_x    ,zx_tmp_fi3d) 
!!               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=dqvdf_w(1:klon,1:klev)/pdtphys
!
               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_vdf_x(1:klon,1:klev)/pdtphys
               CALL histwrite_phy(o_dtvdf_x    ,zx_tmp_fi3d) 
               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_vdf_w(1:klon,1:klev)/pdtphys
               CALL histwrite_phy(o_dtvdf_w    ,zx_tmp_fi3d) 
               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_vdf_x(1:klon,1:klev)/pdtphys
               CALL histwrite_phy(o_dqvdf_x    ,zx_tmp_fi3d) 
               IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_vdf_w(1:klon,1:klev)/pdtphys
!
               CALL histwrite_phy(o_dqvdf_w    ,zx_tmp_fi3d) 
       IF (vars_defined)  zx_tmp_fi2d(1:klon)=-1*sens_x(1:klon)
               CALL histwrite_phy(o_sens_x     ,zx_tmp_fi2d) 
       IF (vars_defined)  zx_tmp_fi2d(1:klon)=-1*sens_w(1:klon)
               CALL histwrite_phy(o_sens_w     ,zx_tmp_fi2d) 
               CALL histwrite_phy(o_flat_x     ,zxfluxlat_x) 
               CALL histwrite_phy(o_flat_w     ,zxfluxlat_w) 
          zx_tmp_fi2d=0.
          IF (vars_defined) THEN
             DO nsrf=1,nbsrf
                   zx_tmp_fi2d(:)=zx_tmp_fi2d(:) &
                        +pctsrf(:,nsrf)*delta_tsurf(:,nsrf)
             ENDDO
          ENDIF
               CALL histwrite_phy(o_delta_tsurf,zx_tmp_fi2d) 
               CALL histwrite_phy(o_cdragh_x   ,cdragh_x   ) 
               CALL histwrite_phy(o_cdragh_w   ,cdragh_w   ) 
               CALL histwrite_phy(o_cdragm_x   ,cdragm_x   ) 
               CALL histwrite_phy(o_cdragm_w   ,cdragm_w   ) 
               CALL histwrite_phy(o_kh         ,kh         ) 
               CALL histwrite_phy(o_kh_x       ,kh_x       ) 
               CALL histwrite_phy(o_kh_w       ,kh_w       ) 
             ENDIF   ! (iflag_pbl_split>=1)
             CALL histwrite_phy(o_ale, ale)
             CALL histwrite_phy(o_alp, alp)
             CALL histwrite_phy(o_cin, cin)
             CALL histwrite_phy(o_WAPE, wake_pe)
             CALL histwrite_phy(o_cv_gen, cv_gen)
             CALL histwrite_phy(o_wake_h, wake_h)
             CALL histwrite_phy(o_wake_dens, wake_dens)
             CALL histwrite_phy(o_wake_s, wake_s)
             CALL histwrite_phy(o_wake_deltat, wake_deltat)
             CALL histwrite_phy(o_wake_deltaq, wake_deltaq)
             CALL histwrite_phy(o_wake_omg, wake_omg)
             IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_wake(1:klon,1:klev) &
                  /pdtphys
             CALL histwrite_phy(o_dtwak, zx_tmp_fi3d)
             IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_wake(1:klon,1:klev)/pdtphys
             CALL histwrite_phy(o_dqwak, zx_tmp_fi3d)
             IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
             CALL histwrite_phy(o_dqwak2d, zx_tmp_fi2d)
          ENDIF ! iflag_wake>=1
          CALL histwrite_phy(o_ftd, ftd)
          CALL histwrite_phy(o_fqd, fqd)
       ENDIF !(iflag_con.EQ.3)
       IF (iflag_con.EQ.3.OR.iflag_con.EQ.30) THEN
          ! sortie RomP convection descente insaturee iflag_con=30
          ! etendue a iflag_con=3 (jyg)
          CALL histwrite_phy(o_Vprecip, Vprecip)
          CALL histwrite_phy(o_qtaa, qtaa)
          CALL histwrite_phy(o_clwaa, clw)
          CALL histwrite_phy(o_wdtrainA, wdtrainA)
          CALL histwrite_phy(o_wdtrainS, wdtrainS)
          CALL histwrite_phy(o_wdtrainM, wdtrainM)
       ENDIF !(iflag_con.EQ.3.or.iflag_con.EQ.30)
!!! nrlmd le 10/04/2012
       IF (iflag_trig_bl>=1) THEN
          CALL histwrite_phy(o_n2, n2)
          CALL histwrite_phy(o_s2, s2)
          CALL histwrite_phy(o_proba_notrig, proba_notrig)
          CALL histwrite_phy(o_random_notrig, random_notrig)
          CALL histwrite_phy(o_ale_bl_stat, ale_bl_stat)
          CALL histwrite_phy(o_ale_bl_trig, ale_bl_trig)
       ENDIF  !(iflag_trig_bl>=1)
       IF (iflag_clos_bl>=1) THEN
          CALL histwrite_phy(o_alp_bl_det, alp_bl_det)
          CALL histwrite_phy(o_alp_bl_fluct_m, alp_bl_fluct_m)
          CALL histwrite_phy(o_alp_bl_fluct_tke,  &
               alp_bl_fluct_tke)
          CALL histwrite_phy(o_alp_bl_conv, alp_bl_conv)
          CALL histwrite_phy(o_alp_bl_stat, alp_bl_stat)
       ENDIF  !(iflag_clos_bl>=1)
!!! fin nrlmd le 10/04/2012
       ! Output of slab ocean variables
       IF (type_ocean=='slab ') THEN
          CALL histwrite_phy(o_slab_bils, slab_wfbils)
          IF (nslay.EQ.1) THEN
              IF (vars_defined) zx_tmp_fi2d(:)=tslab(:,1)
              CALL histwrite_phy(o_tslab, zx_tmp_fi2d)
              IF (vars_defined) zx_tmp_fi2d(:)=dt_qflux(:,1)
              CALL histwrite_phy(o_slab_qflux, zx_tmp_fi2d)
          ELSE
              CALL histwrite_phy(o_tslab, tslab(:,1:nslay))
              CALL histwrite_phy(o_slab_qflux, dt_qflux(:,1:nslay))
          ENDIF
          IF (version_ocean=='sicINT') THEN
              CALL histwrite_phy(o_slab_bilg, slab_bilg)
              CALL histwrite_phy(o_slab_tice, tice)
              CALL histwrite_phy(o_slab_sic, seaice)
          ENDIF
          IF (slab_gm) THEN
             CALL histwrite_phy(o_slab_gm, dt_gm(:,1:nslay))
          ENDIF
          IF (slab_hdiff) THEN
            IF (nslay.EQ.1) THEN
                IF (vars_defined) zx_tmp_fi2d(:)=dt_hdiff(:,1)
                CALL histwrite_phy(o_slab_hdiff, zx_tmp_fi2d)
            ELSE
                CALL histwrite_phy(o_slab_hdiff, dt_hdiff(:,1:nslay))
            ENDIF
          ENDIF
          IF (slab_ekman.GT.0) THEN
            IF (nslay.EQ.1) THEN
                IF (vars_defined) zx_tmp_fi2d(:)=dt_ekman(:,1)
                CALL histwrite_phy(o_slab_ekman, zx_tmp_fi2d)
            ELSE
                CALL histwrite_phy(o_slab_ekman, dt_ekman(:,1:nslay))
            ENDIF
          ENDIF
       ENDIF !type_ocean == force/slab
       CALL histwrite_phy(o_weakinv, weak_inversion)
       CALL histwrite_phy(o_dthmin, dthmin)
       CALL histwrite_phy(o_cldtau, cldtau)
       CALL histwrite_phy(o_cldemi, cldemi)
       CALL histwrite_phy(o_pr_con_l, pmflxr(:,1:klev))
       CALL histwrite_phy(o_pr_con_i, pmflxs(:,1:klev))
       CALL histwrite_phy(o_pr_lsc_l, prfl(:,1:klev))
       CALL histwrite_phy(o_pr_lsc_i, psfl(:,1:klev))
       CALL histwrite_phy(o_re, re)
       CALL histwrite_phy(o_fl, fl)

       IF (ok_bs) THEN
         CALL histwrite_phy(o_pr_bs, bsfl(:,1:klev))
       ENDIF

       IF (vars_defined) THEN
          DO i=1, klon
             IF (zt2m(i).LE.273.15) then
                zx_tmp_fi2d(i)=MAX(0.,rh2m(i)*100.)
             ELSE
                zx_tmp_fi2d(i)=MAX(0.,MIN(100.,rh2m(i)*100.))
             ENDIF
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rh2m, zx_tmp_fi2d)

!       IF (vars_defined) THEN
!          DO i=1, klon
!             zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
!          ENDDO
!       ENDIF
!       CALL histwrite_phy(o_rh2m_min, zx_tmp_fi2d)

!       IF (vars_defined) THEN
!          DO i=1, klon
!             zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
!          ENDDO
!       ENDIF
!       CALL histwrite_phy(o_rh2m_max, zx_tmp_fi2d)

       CALL histwrite_phy(o_qsat2m, qsat2m)
       CALL histwrite_phy(o_tpot, tpot)
       CALL histwrite_phy(o_tpote, tpote)
       IF (vars_defined) zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, is_ter)
       CALL histwrite_phy(o_SWnetOR,  zx_tmp_fi2d)
       CALL histwrite_phy(o_LWdownOR, sollwdown)
       CALL histwrite_phy(o_snowl, snow_lsc)
       CALL histwrite_phy(o_solldown, sollwdown)
       CALL histwrite_phy(o_dtsvdfo, d_ts(:,is_oce))
       CALL histwrite_phy(o_dtsvdft, d_ts(:,is_ter))
       CALL histwrite_phy(o_dtsvdfg,  d_ts(:,is_lic))
       CALL histwrite_phy(o_dtsvdfi, d_ts(:,is_sic))
       CALL histwrite_phy(o_z0m, z0m(:,nbsrf+1))
       CALL histwrite_phy(o_z0h, z0h(:,nbsrf+1))

       ! od550 per species
!--OLIVIER
!This is warranted by treating INCA aerosols as offline aerosols
#ifndef CPP_ECRAD
       IF (flag_aerosol.GT.0) THEN
          IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 

             CALL histwrite_phy(o_od443aer, od443aer)
             CALL histwrite_phy(o_od550aer, od550aer)
             CALL histwrite_phy(o_od865aer, od865aer)
             CALL histwrite_phy(o_abs550aer, abs550aer)
             CALL histwrite_phy(o_od550lt1aer, od550lt1aer)
             CALL histwrite_phy(o_sconcso4, sconcso4)
             CALL histwrite_phy(o_sconcno3, sconcno3)
             CALL histwrite_phy(o_sconcoa, sconcoa)
             CALL histwrite_phy(o_sconcbc, sconcbc)
             CALL histwrite_phy(o_sconcss, sconcss)
             CALL histwrite_phy(o_sconcdust, sconcdust)
             CALL histwrite_phy(o_concso4, concso4)
             CALL histwrite_phy(o_concno3, concno3)
             CALL histwrite_phy(o_concoa, concoa)
             CALL histwrite_phy(o_concbc, concbc)
             CALL histwrite_phy(o_concss, concss)
             CALL histwrite_phy(o_concdust, concdust)
             CALL histwrite_phy(o_loadso4, loadso4)
             CALL histwrite_phy(o_loadoa, loadoa)
             CALL histwrite_phy(o_loadbc, loadbc)
             CALL histwrite_phy(o_loadss, loadss)
             CALL histwrite_phy(o_loaddust, loaddust)
             CALL histwrite_phy(o_loadno3, loadno3)
             CALL histwrite_phy(o_dryod550aer, dryod550aer)
             DO naero = 1, naero_tot-1
                CALL histwrite_phy(o_drytausumaero(naero),drytausum_aero(:,naero))
             END DO
          ENDIF
       ENDIF
       !--STRAT AER
       IF (flag_aerosol.GT.0.OR.flag_aerosol_strat.GT.0) THEN
          DO naero = 1, naero_tot
             CALL histwrite_phy(o_tausumaero(naero),tausum_aero(:,2,naero))
          END DO
       ENDIF
       IF (flag_aerosol_strat.GT.0) THEN
          CALL histwrite_phy(o_tausumaero_lw,tausum_aero(:,6,id_STRAT_phy))
       ENDIF

       CALL histwrite_phy(o_p_tropopause, p_tropopause)
       CALL histwrite_phy(o_t_tropopause, t_tropopause)
       CALL histwrite_phy(o_z_tropopause, z_tropopause)

! ThL -- In the following, we assume read_climoz == 1
       IF (vars_defined) THEN
         zx_tmp_fi2d = 0.0    ! Computation for strato, added ThL
         DO k=1, klev
            zx_tmp_fi2d(:) = zx_tmp_fi2d(:) + wo(:,k,1) * stratomask(:,k) * 1.e3
         END DO
       ENDIF
       CALL histwrite_phy(o_col_O3_strato, zx_tmp_fi2d) ! Added ThL

       IF (vars_defined) THEN
         zx_tmp_fi2d = 0.0    ! Computation for tropo, added ThL
         DO k=1, klev
            zx_tmp_fi2d(:) = zx_tmp_fi2d(:) + wo(:,k,1) * (1.0-stratomask(:,k)) * 1.e3
         END DO
       ENDIF
       CALL histwrite_phy(o_col_O3_tropo, zx_tmp_fi2d)   ! Added ThL
! end add ThL

#ifdef CPP_StratAer
       IF (type_trac=='coag') THEN
          CALL histwrite_phy(o_R2SO4, R2SO4)
          CALL histwrite_phy(o_OCS_lifetime, OCS_lifetime)
          CALL histwrite_phy(o_SO2_lifetime, SO2_lifetime)
          CALL histwrite_phy(o_budg_3D_backgr_ocs,   budg_3D_backgr_ocs)
          CALL histwrite_phy(o_budg_3D_backgr_so2,   budg_3D_backgr_so2)
          CALL histwrite_phy(o_budg_3D_ocs_to_so2,   budg_3D_ocs_to_so2)
          CALL histwrite_phy(o_budg_3D_so2_to_h2so4, budg_3D_so2_to_h2so4)
          CALL histwrite_phy(o_budg_3D_nucl,         budg_3D_nucl)
          CALL histwrite_phy(o_budg_3D_cond_evap,    budg_3D_cond_evap)
          CALL histwrite_phy(o_budg_dep_dry_ocs,     budg_dep_dry_ocs)
          CALL histwrite_phy(o_budg_dep_wet_ocs,     budg_dep_wet_ocs)
          CALL histwrite_phy(o_budg_dep_dry_so2,     budg_dep_dry_so2)
          CALL histwrite_phy(o_budg_dep_wet_so2,     budg_dep_wet_so2)
          CALL histwrite_phy(o_budg_dep_dry_h2so4,   budg_dep_dry_h2so4)
          CALL histwrite_phy(o_budg_dep_wet_h2so4,   budg_dep_wet_h2so4)
          CALL histwrite_phy(o_budg_dep_dry_part,    budg_dep_dry_part)
          CALL histwrite_phy(o_budg_dep_wet_part,    budg_dep_wet_part)
          CALL histwrite_phy(o_budg_emi_ocs,         budg_emi_ocs)
          CALL histwrite_phy(o_budg_emi_so2,         budg_emi_so2)
          CALL histwrite_phy(o_budg_emi_h2so4,       budg_emi_h2so4)
          CALL histwrite_phy(o_budg_emi_part,        budg_emi_part)
          CALL histwrite_phy(o_budg_ocs_to_so2,      budg_ocs_to_so2)
          CALL histwrite_phy(o_budg_so2_to_h2so4,    budg_so2_to_h2so4)
          CALL histwrite_phy(o_budg_h2so4_to_part,   budg_h2so4_to_part)
          CALL histwrite_phy(o_budg_sed_part,        budg_sed_part)
          CALL histwrite_phy(o_surf_PM25_sulf, surf_PM25_sulf)
          CALL histwrite_phy(o_vsed_aer, vsed_aer)
          CALL histwrite_phy(o_f_r_wet, f_r_wet)
          CALL histwrite_phy(o_ext_strat_550, tau_strat_550)
          CALL histwrite_phy(o_ext_strat_1020, tau_strat_1020)
          CALL histwrite_phy(o_tau_strat_550, tausum_strat(:,1))
          CALL histwrite_phy(o_tau_strat_1020, tausum_strat(:,2))
       ENDIF
#endif
       !NL
       IF (ok_volcan .AND. ok_ade) THEN
          DO k=1, klev
             IF (vars_defined) zx_tmp_fi3d(:,k)=heat_volc(:,k)*swradcorr(:)
          ENDDO
          CALL histwrite_phy(o_heat_volc, zx_tmp_fi3d)
          DO k=1, klev
             IF (vars_defined) zx_tmp_fi3d(:,k)=cool_volc(:,k)
          ENDDO
          CALL histwrite_phy(o_cool_volc, zx_tmp_fi3d)
       ENDIF
       IF (ok_ade) THEN
          IF (vars_defined) zx_tmp_fi2d(:)=topswad_aero*swradcorr
          CALL histwrite_phy(o_topswad, zx_tmp_fi2d)
          
          IF (vars_defined) zx_tmp_fi2d(:)=topswad0_aero*swradcorr
          CALL histwrite_phy(o_topswad0, zx_tmp_fi2d)
                    
          IF (vars_defined) zx_tmp_fi2d(:)=solswad_aero*swradcorr
          CALL histwrite_phy(o_solswad, zx_tmp_fi2d)
                    
          IF (vars_defined) zx_tmp_fi2d(:)=solswad0_aero*swradcorr
          CALL histwrite_phy(o_solswad0, zx_tmp_fi2d)
          
          IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 

             CALL histwrite_phy(o_toplwad, toplwad_aero)
             CALL histwrite_phy(o_toplwad0, toplwad0_aero)
             CALL histwrite_phy(o_sollwad, sollwad_aero)
             CALL histwrite_phy(o_sollwad0, sollwad0_aero)
          ENDIF
          !====MS forcing diagnostics
          !ym warning : topsw_aero, solsw_aero, topsw0_aero, solsw0_aero are not defined by model
          !ym => init to 0 in radlwsw_m.F90 ztopsw_aero, zsolsw_aero, ztopsw0_aero, zsolsw0_aero

          IF (vars_defined) zx_tmp_fi2d(:)=topsw_aero(:,1)*swradcorr(:)
          CALL histwrite_phy(o_swtoaas_nat,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=solsw_aero(:,1)*swradcorr(:)
          CALL histwrite_phy(o_swsrfas_nat,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=topsw0_aero(:,1)*swradcorr(:)
          CALL histwrite_phy(o_swtoacs_nat,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=solsw0_aero(:,1)*swradcorr(:)
          CALL histwrite_phy(o_swsrfcs_nat,zx_tmp_fi2d)
          !ant
          IF (vars_defined) zx_tmp_fi2d(:)=topsw_aero(:,2)*swradcorr(:)
          CALL histwrite_phy(o_swtoaas_ant,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=solsw_aero(:,2)*swradcorr(:)
          CALL histwrite_phy(o_swsrfas_ant,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=topsw0_aero(:,2)*swradcorr(:)
          CALL histwrite_phy(o_swtoacs_ant,zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:)=solsw0_aero(:,2)*swradcorr(:)
          CALL histwrite_phy(o_swsrfcs_ant,zx_tmp_fi2d)
          !cf
          IF (.not. aerosol_couple) THEN
             IF (vars_defined) zx_tmp_fi2d(:)=topswcf_aero(:,1)*swradcorr(:)
             CALL histwrite_phy(o_swtoacf_nat,zx_tmp_fi2d)
             IF (vars_defined) zx_tmp_fi2d(:)=solswcf_aero(:,1)*swradcorr(:)
             CALL histwrite_phy(o_swsrfcf_nat,zx_tmp_fi2d)
             IF (vars_defined) zx_tmp_fi2d(:)=topswcf_aero(:,2)*swradcorr(:)
             CALL histwrite_phy(o_swtoacf_ant,zx_tmp_fi2d)
             IF (vars_defined) zx_tmp_fi2d(:)=solswcf_aero(:,2)*swradcorr(:)
             CALL histwrite_phy(o_swsrfcf_ant,zx_tmp_fi2d)
             IF (vars_defined) zx_tmp_fi2d(:)=topswcf_aero(:,3)*swradcorr(:)
             CALL histwrite_phy(o_swtoacf_zero,zx_tmp_fi2d)
             IF (vars_defined) zx_tmp_fi2d(:)=solswcf_aero(:,3)*swradcorr(:)
             CALL histwrite_phy(o_swsrfcf_zero,zx_tmp_fi2d)
          ENDIF
          !====MS forcing diagnostics
       ENDIF
       IF (ok_aie) THEN
          IF (vars_defined) zx_tmp_fi2d(:)= topswai_aero*swradcorr
          CALL histwrite_phy(o_topswai, zx_tmp_fi2d)
          
          IF (vars_defined) zx_tmp_fi2d(:)=toplwai_aero*swradcorr
          CALL histwrite_phy(o_toplwai, zx_tmp_fi2d)
          
          IF (vars_defined) zx_tmp_fi2d(:)=solswai_aero*swradcorr
          CALL histwrite_phy(o_solswai, zx_tmp_fi2d)
          
          IF (vars_defined) zx_tmp_fi2d(:)=sollwai_aero*swradcorr
          CALL histwrite_phy(o_sollwai, zx_tmp_fi2d)
       ENDIF
       IF (flag_aerosol.GT.0.AND.ok_cdnc) THEN
          CALL histwrite_phy(o_scdnc, scdnc)
          CALL histwrite_phy(o_cldncl, cldncl)
          CALL histwrite_phy(o_reffclws, reffclws)
          CALL histwrite_phy(o_reffclwc, reffclwc)
          CALL histwrite_phy(o_cldnvi, cldnvi)
          CALL histwrite_phy(o_lcc, lcc)
          CALL histwrite_phy(o_lcc3d, lcc3d)
          CALL histwrite_phy(o_lcc3dcon, lcc3dcon)
          CALL histwrite_phy(o_lcc3dstra, lcc3dstra)
          CALL histwrite_phy(o_icc3dcon, icc3dcon)
          CALL histwrite_phy(o_icc3dstra, icc3dstra)
          CALL histwrite_phy(o_cldicemxrat, zfice)
          IF (vars_defined) zx_tmp_fi3d(:,:)=1-zfice(:,:)
          CALL histwrite_phy(o_cldwatmxrat, zx_tmp_fi3d)
          CALL histwrite_phy(o_reffclwtop, reffclwtop)
       ENDIF
       ! Champs 3D:
       IF (ok_ade .OR. ok_aie) then
          IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
             CALL histwrite_phy(o_ec550aer, ec550aer)
          ENDIF
       ENDIF

       CALL histwrite_phy(o_lwcon, flwc)
       CALL histwrite_phy(o_iwcon, fiwc)
       CALL histwrite_phy(o_temp, t_seri)
       CALL histwrite_phy(o_theta, theta)
       CALL histwrite_phy(o_ovapinit, qx(:,:,ivap))
       CALL histwrite_phy(o_ovap, q_seri)
       CALL histwrite_phy(o_oliq, ql_seri)
       !FC
       CALL histwrite_phy(o_zxfluxt, zxfluxt)
       CALL histwrite_phy(o_zxfluxq, zx_tmp_fi3d)
       !FC

       IF (vars_defined) zx_tmp_fi3d = ql_seri+qs_seri
       CALL histwrite_phy(o_ocond, zx_tmp_fi3d)

       CALL histwrite_phy(o_geop, zphi)
       CALL histwrite_phy(o_vitu, u_seri)
       CALL histwrite_phy(o_vitv, v_seri)
       CALL histwrite_phy(o_vitw, omega)
       CALL histwrite_phy(o_pres, pplay)
       CALL histwrite_phy(o_paprs, paprs(:,1:klev))
       
       IF (vars_defined) zx_tmp_fi3d = zphi/RG
       CALL histwrite_phy(o_zfull,zx_tmp_fi3d)

       IF (ok_bs) THEN
          CALL histwrite_phy(o_qbs, qbs_seri)
       ENDIF

       IF (using_xios) THEN
!solbnd begin
#ifdef CPP_RRTM
         IF (iflag_rrtm.EQ.1) THEN
           IF (vars_defined) THEN
             DO ISW=1, NSW
               zx_tmp_fi3dsp(:,ISW) = swdn(:,klevp1)*swradcorr(:)*RSUN(ISW)
             ENDDO
             CALL histwrite_phy(o_solbnd, zx_tmp_fi3dsp)
           ENDIF
         ENDIF
#endif
!solbnd end
       ENDIF
#endif

       IF (flag_aerosol_strat.EQ.2) THEN
         CALL histwrite_phy(o_stratomask, stratomask)
       ENDIF
      
       IF (vars_defined)  THEN
        zx_tmp_fi3d(:,1)= pphis(:)/RG
        DO k = 2, klev
         DO i = 1, klon
            zx_tmp_fi3d(i,k) = zphi(i,k-1)/RG + &
                          (zphi(i,k)-zphi(i,k-1))/RG * &
                          (paprs(i,k)-pplay(i,k-1))/(pplay(i,k)-pplay(i,k-1))
         ENDDO
        ENDDO
       ENDIF
       CALL histwrite_phy(o_zhalf, zx_tmp_fi3d)
       CALL histwrite_phy(o_rneb, cldfra)
       CALL histwrite_phy(o_rnebcon, rnebcon)
       CALL histwrite_phy(o_rnebls, rneb)
       CALL histwrite_phy(o_rneblsvol, rneblsvol)
       IF (vars_defined)  THEN
          DO k=1, klev
             DO i=1, klon
                zx_tmp_fi3d(i,k)=cldfra(i,k)*JrNt(i)
             ENDDO
          ENDDO
       ENDIF
       CALL histwrite_phy(o_rnebjn, zx_tmp_fi3d)
       CALL histwrite_phy(o_rhum, zx_rh)
       IF (iflag_ice_thermo .GT. 0) THEN
          IF (vars_defined) zx_tmp_fi3d = zx_rhl * 100.
          CALL histwrite_phy(o_rhl, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d = zx_rhi * 100.
          CALL histwrite_phy(o_rhi, zx_tmp_fi3d)
       ENDIF
     
       IF (ok_new_lscp) THEN
           CALL histwrite_phy(o_pfraclr, pfraclr)
           CALL histwrite_phy(o_pfracld, pfracld)
       ENDIF

!--aviation & supersaturation
       IF (ok_ice_sursat) THEN
         CALL histwrite_phy(o_oclr, qclr)
         CALL histwrite_phy(o_ocld, qcld)
         CALL histwrite_phy(o_oss, qss)
         CALL histwrite_phy(o_ovc, qvc)
         CALL histwrite_phy(o_rnebclr, rnebclr)
         CALL histwrite_phy(o_rnebss, rnebss)
         CALL histwrite_phy(o_rnebseri, rneb_seri)
         CALL histwrite_phy(o_gammass, gamma_ss)
         CALL histwrite_phy(o_N1_ss, N1_ss)
         CALL histwrite_phy(o_N2_ss, N2_ss)
         CALL histwrite_phy(o_drnebsub, drneb_sub)
         CALL histwrite_phy(o_drnebcon, drneb_con)
         CALL histwrite_phy(o_drnebtur, drneb_tur)
         CALL histwrite_phy(o_drnebavi, drneb_avi)
         CALL histwrite_phy(o_qsatl, zqsatl)
         CALL histwrite_phy(o_qsats, zqsats)
         CALL histwrite_phy(o_Tcontr, Tcontr)
         CALL histwrite_phy(o_qcontr, qcontr)
         CALL histwrite_phy(o_qcontr2, qcontr2)
         CALL histwrite_phy(o_fcontrN, fcontrN)
         CALL histwrite_phy(o_fcontrP, fcontrP)
       ENDIF
       IF (ok_plane_contrail) THEN
         CALL histwrite_phy(o_flight_m, flight_m)
       ENDIF
       IF (ok_plane_h2o) THEN
         CALL histwrite_phy(o_flight_h2o, flight_h2o)
       ENDIF
       
       IF (vars_defined) zx_tmp_fi3d = wo(:, :, 1) * dobson_u * 1e3 / zmasse / rmo3 * rmd
       CALL histwrite_phy(o_ozone, zx_tmp_fi3d)

       IF (read_climoz == 2) THEN
         IF (vars_defined) zx_tmp_fi3d = wo(:, :, 2) * dobson_u * 1e3 / zmasse / rmo3 * rmd 
         CALL histwrite_phy(o_ozone_light, zx_tmp_fi3d)
       ENDIF

       CALL histwrite_phy(o_duphy, d_u)

       CALL histwrite_phy(o_dtphy, d_t)

       CALL histwrite_phy(o_dqphy,  d_qx(:,:,ivap))
       IF (vars_defined) CALL water_int(klon,klev,d_qx(:,:,ivap),zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqphy2d,  zx_tmp_fi2d)

       CALL histwrite_phy(o_dqlphy,  d_qx(:,:,iliq))
       IF (vars_defined) CALL water_int(klon,klev,d_qx(:,:,iliq),zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqlphy2d,  zx_tmp_fi2d)

       IF (nqo.EQ.3) THEN
       CALL histwrite_phy(o_dqsphy,  d_qx(:,:,isol))
       IF (vars_defined) CALL water_int(klon,klev,d_qx(:,:,isol),zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqsphy2d,  zx_tmp_fi2d)
       ELSE 
       zx_tmp_fi3d=0.0
       CALL histwrite_phy(o_dqsphy,  zx_tmp_fi3d)
       zx_tmp_fi2d=0.0
       CALL histwrite_phy(o_dqsphy2d,  zx_tmp_fi2d)
       ENDIF


       IF (ok_bs) THEN
       CALL histwrite_phy(o_dqbsphy,  d_qx(:,:,ibs))
       IF (vars_defined) CALL water_int(klon,klev,d_qx(:,:,ibs),zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqbsphy2d,  zx_tmp_fi2d)
       ELSE
       zx_tmp_fi3d=0.0
       CALL histwrite_phy(o_dqbsphy,  zx_tmp_fi3d)
       zx_tmp_fi2d=0.0
       CALL histwrite_phy(o_dqbsphy2d,  zx_tmp_fi2d)
       ENDIF

       DO nsrf=1, nbsrf
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = falb1( 1 : klon, nsrf)
          CALL histwrite_phy(o_albe_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = z0m( 1 : klon, nsrf)
          CALL histwrite_phy(o_z0m_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = z0h( 1 : klon, nsrf)
          CALL histwrite_phy(o_z0h_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = agesno( 1 : klon, nsrf)
          CALL histwrite_phy(o_ages_srf(nsrf), zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(1 : klon) = snow( 1 : klon, nsrf)
          CALL histwrite_phy(o_snow_srf(nsrf), zx_tmp_fi2d)
       ENDDO !nsrf=1, nbsrf
       CALL histwrite_phy(o_alb1, albsol1)
       CALL histwrite_phy(o_alb2, albsol2)
       !FH Sorties pour la couche limite
       IF (iflag_pbl>1) THEN
          zx_tmp_fi3d=0.
          IF (vars_defined) THEN
             DO nsrf=1,nbsrf
                DO k=1,klev
                   zx_tmp_fi3d(:,k)=zx_tmp_fi3d(:,k) &
                        +pctsrf(:,nsrf)*pbl_tke(:,k,nsrf)
                ENDDO
             ENDDO
          ENDIF
          CALL histwrite_phy(o_tke, zx_tmp_fi3d)
          CALL histwrite_phy(o_tke_max, zx_tmp_fi3d)  

       ENDIF

       CALL histwrite_phy(o_kz, coefh(:,:,is_ave))

       CALL histwrite_phy(o_kz_max, coefh(:,:,is_ave))

       CALL histwrite_phy(o_clwcon, clwcon0)
       CALL histwrite_phy(o_dtdyn, d_t_dyn)

       CALL histwrite_phy(o_dqdyn, d_q_dyn)

       CALL histwrite_phy(o_dqdyn2d,d_q_dyn2d)

       CALL histwrite_phy(o_dqldyn, d_ql_dyn)

       CALL histwrite_phy(o_dqldyn2d, d_ql_dyn2d)

       CALL histwrite_phy(o_dqsdyn, d_qs_dyn)

       CALL histwrite_phy(o_dqsdyn2d, d_qs_dyn2d)

       IF (ok_bs) THEN
       CALL histwrite_phy(o_dqbsdyn, d_qbs_dyn)
       CALL histwrite_phy(o_dqbsdyn2d, d_qbs_dyn2d)
       ENDIF

       CALL histwrite_phy(o_dudyn, d_u_dyn)
       CALL histwrite_phy(o_dvdyn, d_v_dyn)

       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dtcon, zx_tmp_fi3d)
       IF (iflag_thermals.EQ.0) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys + &
                  d_t_ajsb(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tntc, zx_tmp_fi3d)
       ELSE IF(iflag_thermals.GE.1.AND.iflag_wake.EQ.1) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys + &
                  d_t_ajs(1:klon,1:klev)/pdtphys + &
                  d_t_wake(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tntc, zx_tmp_fi3d)
       ENDIF
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_ducon, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dvcon, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqcon, zx_tmp_fi3d)
       IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqcon2d, zx_tmp_fi2d)

       IF (iflag_thermals.EQ.0) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_tnhusc, zx_tmp_fi3d)
       ELSE IF (iflag_thermals.GE.1.AND.iflag_wake.EQ.1) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys + &
                  d_q_ajs(1:klon,1:klev)/pdtphys + &
                  d_q_wake(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_tnhusc, zx_tmp_fi3d)
       ENDIF

       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lsc(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlsc, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon, 1:klev)=(d_t_lsc(1:klon,1:klev)+ &
            d_t_eva(1:klon,1:klev))/pdtphys
       CALL histwrite_phy(o_dtlschr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lsc(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqlsc, zx_tmp_fi3d)
       IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqlsc2d, zx_tmp_fi2d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=beta_prec(1:klon,1:klev)
       CALL histwrite_phy(o_beta_prec, zx_tmp_fi3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Sorties specifiques a la separation thermiques/non thermiques
       IF (iflag_thermals>=1) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lscth(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlscth, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lscst(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlscst, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lscth(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqlscth, zx_tmp_fi3d)
          IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
          CALL histwrite_phy(o_dqlscth2d, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_lscst(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqlscst, zx_tmp_fi3d)
          IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
          CALL histwrite_phy(o_dqlscst2d, zx_tmp_fi2d)
          CALL histwrite_phy(o_plulth, plul_th)
          CALL histwrite_phy(o_plulst, plul_st)
          IF (vars_defined) THEN
             DO i=1,klon
                zx_tmp_fi2d(1:klon)=lmax_th(:)
             ENDDO
          ENDIF
          CALL histwrite_phy(o_lmaxth, zx_tmp_fi2d)
          IF (vars_defined) THEN
             DO k=1,klev
                DO i=1,klon
                   IF (ptconvth(i,k)) THEN
                      zx_tmp_fi3d(i,k)=1.
                   ELSE
                      zx_tmp_fi3d(i,k)=0.
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          CALL histwrite_phy(o_ptconvth, zx_tmp_fi3d)
       ENDIF ! iflag_thermals>=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       zpt_conv = 0.
       WHERE (ptconv) zpt_conv = 1.
       CALL histwrite_phy(o_ptconv, zpt_conv)
!!       IF (vars_defined)         zx_tmp_fi2d=float(itau_con)/float(itap)
!!       CALL histwrite_phy(o_ftime_con, zx_tmp_fi2d)
       IF (vars_defined) THEN
          zpt_conv2d(:) = 0.
          DO k=1,klev
            WHERE (ptconv(:,k)) zpt_conv2d(:) = 1.
          ENDDO
       ENDIF
       CALL histwrite_phy(o_ftime_deepcv, zpt_conv2d)
       IF (vars_defined) THEN
          zx_tmp_fi2d(:) = 0.
          DO k=1,klev
            WHERE (ptconvth(:,k)) zx_tmp_fi2d(:) = 1.
          ENDDO
       ENDIF
       CALL histwrite_phy(o_ftime_th, zx_tmp_fi2d)
       IF (vars_defined) THEN
           zx_tmp_fi2d(:) = max(zx_tmp_fi2d(:),zpt_conv2d(:))
       ENDIF
       CALL histwrite_phy(o_ftime_con, zx_tmp_fi2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtvdf, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_diss(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtdis, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqvdf, zx_tmp_fi3d)
       IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqvdf2d, zx_tmp_fi2d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_eva(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dteva, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_eva(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqeva, zx_tmp_fi3d)
       IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqeva2d, zx_tmp_fi2d)
       CALL histwrite_phy(o_ratqs, ratqs)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_ajs(1:klon,1:klev)/pdtphys - &
               d_t_ajsb(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dtthe, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_u_ajs(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_duthe, zx_tmp_fi3d) 
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_v_ajs(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_dvthe, zx_tmp_fi3d)

       IF (ok_bs) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_qbs_vdf(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqbsvdf, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_qbs_bs(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqbsbs, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_bs(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dqbs, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_bs(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtbs, zx_tmp_fi3d)
       ENDIF

       IF (iflag_thermals>=1) THEN
          ! Pour l instant 0 a y reflichir pour les thermiques
          ! regroupe avec ftime_deepcv et ftime_con
          !!zx_tmp_fi2d=0. 
          !!CALL histwrite_phy(o_ftime_th, zx_tmp_fi2d)
          CALL histwrite_phy(o_f_th, fm_therm)
          CALL histwrite_phy(o_e_th, entr_therm)
          CALL histwrite_phy(o_w_th, zw2)
          CALL histwrite_phy(o_q_th, zqasc)
          CALL histwrite_phy(o_a_th, fraca)
          CALL histwrite_phy(o_cloudth_sth, cloudth_sth)
          CALL histwrite_phy(o_cloudth_senv, cloudth_senv)
          CALL histwrite_phy(o_cloudth_sigmath, cloudth_sigmath)
          CALL histwrite_phy(o_cloudth_sigmaenv, cloudth_sigmaenv)
          CALL histwrite_phy(o_d_th, detr_therm)
          CALL histwrite_phy(o_f0_th, f0)
          CALL histwrite_phy(o_zmax_th, zmax_th)
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=d_q_ajs(1:klon,1:klev)/pdtphys - &
                  d_q_ajsb(1:klon,1:klev)/pdtphys
          ENDIF
          CALL histwrite_phy(o_dqthe, zx_tmp_fi3d)
          IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
          CALL histwrite_phy(o_dqthe2d, zx_tmp_fi2d)
       ENDIF !iflag_thermals
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_ajsb(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtajs, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_q_ajsb(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dqajs, zx_tmp_fi3d)
       IF (vars_defined) CALL water_int(klon,klev,zx_tmp_fi3d,zmasse,zx_tmp_fi2d)
       CALL histwrite_phy(o_dqajs2d, zx_tmp_fi2d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_swr(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtswr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_sw0(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtsw0, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lwr(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlwr, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lw0(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtlw0, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_ec(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dtec, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_duvdf, zx_tmp_fi3d)
       IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_vdf(1:klon,1:klev)/pdtphys
       CALL histwrite_phy(o_dvvdf, zx_tmp_fi3d)
       IF (ok_orodr) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_duoro, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dvoro, zx_tmp_fi3d)
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_oro(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtoro, zx_tmp_fi3d)
       ENDIF
       IF (ok_orolf) THEN
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_u_lif(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dulif, zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_v_lif(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dvlif, zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_lif(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dtlif, zx_tmp_fi3d)
       ENDIF

       IF (ok_hines) THEN
          IF (vars_defined) zx_tmp_fi3d=du_gwd_hines/pdtphys
          CALL histwrite_phy(o_du_gwd_hines, zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d= dv_gwd_hines/pdtphys          
          CALL histwrite_phy(o_dv_gwd_hines, zx_tmp_fi3d)
          
          IF (vars_defined) zx_tmp_fi3d(1:klon,1:klev)=d_t_hin(1:klon,1:klev)/pdtphys
          CALL histwrite_phy(o_dthin, zx_tmp_fi3d)
          CALL histwrite_phy(o_ustr_gwd_hines, zustr_gwd_hines)
          CALL histwrite_phy(o_vstr_gwd_hines, zvstr_gwd_hines)
       ENDIF

       IF (.not. ok_hines .and. ok_gwd_rando) THEN
          IF (vars_defined)  zx_tmp_fi3d=du_gwd_front / pdtphys
          CALL histwrite_phy(o_du_gwd_front, zx_tmp_fi3d)
          
          IF (vars_defined)  zx_tmp_fi3d=dv_gwd_front / pdtphys
          CALL histwrite_phy(o_dv_gwd_front, zx_tmp_fi3d)
          
          CALL histwrite_phy(o_ustr_gwd_front, zustr_gwd_front)
          CALL histwrite_phy(o_vstr_gwd_front, zvstr_gwd_front)
       ENDIF

       IF (ok_gwd_rando) THEN
          IF (vars_defined)  zx_tmp_fi3d=du_gwd_rando / pdtphys
          CALL histwrite_phy(o_du_gwd_rando, zx_tmp_fi3d)
          
          IF (vars_defined)  zx_tmp_fi3d=dv_gwd_rando / pdtphys
          CALL histwrite_phy(o_dv_gwd_rando, zx_tmp_fi3d)
          CALL histwrite_phy(o_ustr_gwd_rando, zustr_gwd_rando)
          CALL histwrite_phy(o_vstr_gwd_rando, zvstr_gwd_rando)
          CALL histwrite_phy(o_east_gwstress, east_gwstress )
          CALL histwrite_phy(o_west_gwstress, west_gwstress )
       ENDIF

       IF (ok_qch4) THEN
          IF (vars_defined) zx_tmp_fi3d=d_q_ch4 / pdtphys
          CALL histwrite_phy(o_dqch4, zx_tmp_fi3d)
       ENDIF
       
       IF (vars_defined) THEN
         DO k=1, klevp1
           zx_tmp_fi3d1(:,k)=swup(:,k)*swradcorr(:)
         ENDDO
       ENDIF
       
       CALL histwrite_phy(o_rsu, zx_tmp_fi3d1)


       IF (vars_defined) THEN
         DO k=1, klevp1
           zx_tmp_fi3d1(:,k)=swdn(:,k)*swradcorr(:)
         ENDDO
       ENDIF
       
       CALL histwrite_phy(o_rsd, zx_tmp_fi3d1)

       IF (vars_defined) THEN
         DO k=1, klevp1
           zx_tmp_fi3d1(:,k)=swup0(:,k)*swradcorr(:)
         ENDDO
       ENDIF
       
       CALL histwrite_phy(o_rsucs, zx_tmp_fi3d1)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             DO k=1, klevp1
                zx_tmp_fi3d1(:,k)=swupc0(:,k)*swradcorr(:)
             ENDDO
          ENDIF
          CALL histwrite_phy(o_rsucsaf, zx_tmp_fi3d1)
       ENDIF

       IF (vars_defined) THEN
         DO k=1, klevp1
           zx_tmp_fi3d1(:,k)=swdn0(:,k)*swradcorr(:)
         ENDDO
       ENDIF
       CALL histwrite_phy(o_rsdcs, zx_tmp_fi3d1)

       IF (type_trac/='inca' .OR. config_inca=='aeNP') THEN 
          IF (vars_defined) THEN
             DO k=1, klevp1
                zx_tmp_fi3d1(:,k)=swdnc0(:,k)*swradcorr(:)
             ENDDO
          ENDIF
          CALL histwrite_phy(o_rsdcsaf, zx_tmp_fi3d1)
       ENDIF

       CALL histwrite_phy(o_rlu, lwup)
       CALL histwrite_phy(o_rld, lwdn)
       CALL histwrite_phy(o_rlucs, lwup0)
       CALL histwrite_phy(o_rldcs, lwdn0)

       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t(1:klon,1:klev)+ &
               d_t_dyn(1:klon,1:klev)
       ENDIF
       CALL histwrite_phy(o_tnt, zx_tmp_fi3d)

       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_t_swr(1:klon,1:klev)/pdtphys + &
               d_t_lwr(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_tntr, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)= (d_t_lsc(1:klon,1:klev)+ &
               d_t_eva(1:klon,1:klev)+ &
               d_t_vdf(1:klon,1:klev))/pdtphys
       ENDIF
       CALL histwrite_phy(o_tntscpbl, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_qx(1:klon,1:klev,ivap)+ &
               d_q_dyn(1:klon,1:klev)
       ENDIF
       CALL histwrite_phy(o_tnhus, zx_tmp_fi3d)
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=d_q_lsc(1:klon,1:klev)/pdtphys+ &
               d_q_eva(1:klon,1:klev)/pdtphys
       ENDIF
       CALL histwrite_phy(o_tnhusscpbl, zx_tmp_fi3d)
       CALL histwrite_phy(o_evu, coefm(:,:,is_ave))
       IF (vars_defined) THEN
          zx_tmp_fi3d(1:klon,1:klev)=q_seri(1:klon,1:klev)+ &
               ql_seri(1:klon,1:klev) 
       ENDIF
       CALL histwrite_phy(o_h2o, zx_tmp_fi3d)
       IF (iflag_con >= 3) THEN
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=-1 * (dnwd(1:klon,1:klev)+ &
                  dnwd0(1:klon,1:klev)) 
          ENDIF
          CALL histwrite_phy(o_mcd, zx_tmp_fi3d)
          IF (vars_defined) THEN
             zx_tmp_fi3d(1:klon,1:klev)=upwd(1:klon,1:klev) + &
                  dnwd(1:klon,1:klev)+ dnwd0(1:klon,1:klev) 
          ENDIF
          CALL histwrite_phy(o_dmc, zx_tmp_fi3d)
       ELSE IF (iflag_con == 2) THEN
          CALL histwrite_phy(o_mcd,  pmfd)
          IF (vars_defined) zx_tmp_fi3d = pmfu + pmfd
          CALL histwrite_phy(o_dmc,  zx_tmp_fi3d)
       ENDIF
       CALL histwrite_phy(o_ref_liq, ref_liq)
       CALL histwrite_phy(o_ref_ice, ref_ice)
!
       IF (ok_4xCO2atm) THEN
          IF (vars_defined) zx_tmp_fi2d(:) = swupp(:,klevp1)*swradcorr(:)
          CALL histwrite_phy(o_rsut4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:) = lwupp(:,klevp1)
          CALL histwrite_phy(o_rlut4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:) = swup0p(:,klevp1)*swradcorr(:)
          CALL histwrite_phy(o_rsutcs4co2, zx_tmp_fi2d)
          IF (vars_defined) zx_tmp_fi2d(:) = lwup0p(:,klevp1)
          CALL histwrite_phy(o_rlutcs4co2, zx_tmp_fi2d)
          IF (vars_defined) THEN
            DO k=1, klevp1
              zx_tmp_fi3d1(:,k)=swupp(:,k)*swradcorr(:)
            ENDDO
          ENDIF
          CALL histwrite_phy(o_rsu4co2, zx_tmp_fi3d1)
          IF (vars_defined) THEN
            DO k=1, klevp1
              zx_tmp_fi3d1(:,k)=swup0p(:,k)*swradcorr(:)
            ENDDO
          ENDIF
          CALL histwrite_phy(o_rsucs4co2, zx_tmp_fi3d1)
          IF (vars_defined) THEN
            DO k=1, klevp1
              zx_tmp_fi3d1(:,k)=swdnp(:,k)*swradcorr(:)
            ENDDO
          ENDIF
          CALL histwrite_phy(o_rsd4co2, zx_tmp_fi3d1)
          IF (vars_defined) THEN
            DO k=1, klevp1
              zx_tmp_fi3d1(:,k)=swdn0p(:,k)*swradcorr(:)
            ENDDO
          ENDIF
          CALL histwrite_phy(o_rsdcs4co2, zx_tmp_fi3d1)
          CALL histwrite_phy(o_rlu4co2, lwupp)
          CALL histwrite_phy(o_rlucs4co2, lwup0p)
          CALL histwrite_phy(o_rld4co2, lwdnp)
          CALL histwrite_phy(o_rldcs4co2, lwdn0p)
       ENDIF !ok_4xCO2atm
!!!!!!!!!!!! Sorties niveaux de pression NMC !!!!!!!!!!!!!!!!!!!!
#ifdef CPP_IOIPSL
       IF (.NOT. using_xios) THEN
         IF (.NOT.ok_all_xml) THEN 
           ! ATTENTION, LES ANCIENS HISTWRITE ONT ETES CONSERVES EN ATTENDANT MIEUX:
           ! Champs interpolles sur des niveaux de pression
           DO iff=7, nfiles-1 !--OB: here we deal with files 7,8,9

             CALL histwrite_phy(o_tnondef,tnondef(:,:,iff-6),iff)
             CALL histwrite_phy(o_ta,twriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_zg,phiwriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_hus,qwriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_hur,rhwriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_ua,uwriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_va,vwriteSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_wap,wwriteSTD(:,:,iff-6),iff)
             IF (vars_defined) THEN
               DO k=1, nlevSTD
                  DO i=1, klon
                     IF (tnondef(i,k,iff-6).NE.missing_val) THEN
                       IF (freq_outNMC(iff-6).LT.0) THEN
                          freq_moyNMC(iff-6)=(mth_len*un_jour)/freq_calNMC(iff-6)
                       ELSE
                          freq_moyNMC(iff-6)=freq_outNMC(iff-6)/freq_calNMC(iff-6)
                       ENDIF
                       zx_tmp_fi3d_STD(i,k) = (100.*tnondef(i,k,iff-6))/freq_moyNMC(iff-6)
                     ELSE
                       zx_tmp_fi3d_STD(i,k) = missing_val
                     ENDIF
                  ENDDO
               ENDDO
             ENDIF
             CALL histwrite_phy(o_psbg,zx_tmp_fi3d_STD,iff)
             IF (vars_defined) THEN
               DO k=1, nlevSTD
                  DO i=1, klon
                    IF (O3sumSTD(i,k,iff-6).NE.missing_val) THEN
                       zx_tmp_fi3d_STD(i,k) = O3sumSTD(i,k,iff-6) * 1.e+9
                    ELSE
                       zx_tmp_fi3d_STD(i,k) = missing_val
                    ENDIF
                  ENDDO
               ENDDO !k=1, nlevSTD
             ENDIF
             CALL histwrite_phy(o_tro3,zx_tmp_fi3d_STD,iff)
             IF (read_climoz == 2) THEN
               IF (vars_defined) THEN
                 DO k=1, nlevSTD
                   DO i=1, klon
                      IF (O3daysumSTD(i,k,iff-6).NE.missing_val) THEN
                         zx_tmp_fi3d_STD(i,k) = O3daysumSTD(i,k,iff-6) * 1.e+9
                      ELSE
                         zx_tmp_fi3d_STD(i,k) = missing_val
                      ENDIF
                   ENDDO
                 ENDDO !k=1, nlevSTD
               ENDIF
               CALL histwrite_phy(o_tro3_daylight,zx_tmp_fi3d_STD,iff)
             ENDIF
             CALL histwrite_phy(o_uxv,uvsumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_vxq,vqsumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_vxT,vTsumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_wxq,wqsumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_vxphi,vphisumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_wxT,wTsumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_uxu,u2sumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_vxv,v2sumSTD(:,:,iff-6),iff)
             CALL histwrite_phy(o_TxT,T2sumSTD(:,:,iff-6),iff)
           ENDDO !nfiles
         ENDIF
       ENDIF !.NOT. using_xios
#endif

       IF (using_xios) THEN
         IF (ok_all_xml) THEN 
    !      DO iff=7, nfiles

!         CALL histwrite_phy(o_tnondef,tnondef(:,:,3))
          CALL histwrite_phy(o_ta,tlevSTD(:,:))
          CALL histwrite_phy(o_zg,philevSTD(:,:))
          CALL histwrite_phy(o_hus,qlevSTD(:,:))
          CALL histwrite_phy(o_hur,rhlevSTD(:,:))
          CALL histwrite_phy(o_ua,ulevSTD(:,:))
          CALL histwrite_phy(o_va,vlevSTD(:,:))
          CALL histwrite_phy(o_wap,wlevSTD(:,:))
!         IF (vars_defined) THEN
!            DO k=1, nlevSTD
!               DO i=1, klon
!                  IF (tnondef(i,k,3).NE.missing_val) THEN
!                     IF (freq_outNMC(iff-6).LT.0) THEN
!                        freq_moyNMC(iff-6)=(mth_len*un_jour)/freq_calNMC(iff-6)
!                     ELSE
!                        freq_moyNMC(iff-6)=freq_outNMC(iff-6)/freq_calNMC(iff-6)
!                     ENDIF
!                     zx_tmp_fi3d_STD(i,k) = (100.*tnondef(i,k,3))/freq_moyNMC(iff-6)
!                  ELSE
!                     zx_tmp_fi3d_STD(i,k) = missing_val
!                  ENDIF
!               ENDDO
!            ENDDO
!         ENDIF
!         CALL histwrite_phy(o_psbg,zx_tmp_fi3d_STD)
          IF (vars_defined) THEN
             DO k=1, nlevSTD
                DO i=1, klon
                   IF (O3STD(i,k).NE.missing_val) THEN
                      zx_tmp_fi3d_STD(i,k) = O3STD(i,k) * 1.e+9
                   ELSE
                      zx_tmp_fi3d_STD(i,k) = missing_val
                   ENDIF
                ENDDO
             ENDDO !k=1, nlevSTD
          ENDIF
          CALL histwrite_phy(o_tro3,zx_tmp_fi3d_STD)
          IF (read_climoz == 2) THEN
             IF (vars_defined) THEN
                DO k=1, nlevSTD
                   DO i=1, klon
                      IF (O3daySTD(i,k).NE.missing_val) THEN
                         zx_tmp_fi3d_STD(i,k) = O3daySTD(i,k) * 1.e+9
                      ELSE
                         zx_tmp_fi3d_STD(i,k) = missing_val
                      ENDIF
                   ENDDO
                ENDDO !k=1, nlevSTD
             ENDIF
             CALL histwrite_phy(o_tro3_daylight,zx_tmp_fi3d_STD)
          ENDIF
          CALL histwrite_phy(o_uxv,uvSTD(:,:))
          CALL histwrite_phy(o_vxq,vqSTD(:,:))
          CALL histwrite_phy(o_vxT,vTSTD(:,:))
          CALL histwrite_phy(o_wxq,wqSTD(:,:))
          CALL histwrite_phy(o_vxphi,vphiSTD(:,:))
          CALL histwrite_phy(o_wxT,wTSTD(:,:))
          CALL histwrite_phy(o_uxu,u2STD(:,:))
          CALL histwrite_phy(o_vxv,v2STD(:,:))
          CALL histwrite_phy(o_TxT,T2STD(:,:))
!      ENDDO !nfiles
    ENDIF
  ENDIF !using_xios
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (iflag_phytrac == 1 ) then
!
         IF (type_trac == 'co2i') THEN
           itr = 0
           DO iq = 1, nqtot
             IF(.NOT.tracers(iq)%isInPhysics) CYCLE
             itr = itr + 1
!            write(*,*) 'phys_output_write_mod 2370: itr=',itr
             !--3D fields
             CALL histwrite_phy(o_trac(itr), tr_seri(:,:,itr))
             CALL histwrite_phy(o_dtr_vdf(itr),d_tr_cl(:,:,itr))
             CALL histwrite_phy(o_dtr_the(itr),d_tr_th(:,:,itr))
             CALL histwrite_phy(o_dtr_con(itr),d_tr_cv(:,:,itr))
             !--2D fields
             !--CO2 burden
             zx_tmp_fi2d=0.
             IF (vars_defined) THEN
                DO k=1,klev
                   zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+zmasse(:,k)*tr_seri(:,k,itr)
                ENDDO
             ENDIF
             CALL histwrite_phy(o_trac_cum(itr), zx_tmp_fi2d)
           ENDDO !--iq
           !--CO2 net fluxes
           CALL histwrite_phy(o_flx_co2_land,  fco2_land)
           CALL histwrite_phy(o_flx_co2_ocean, fco2_ocean)
           CALL histwrite_phy(o_flx_co2_ocean_cor, fco2_ocean_cor)
           CALL histwrite_phy(o_flx_co2_land_cor, fco2_land_cor)
           CALL histwrite_phy(o_flx_co2_ff,    fco2_ff)
           CALL histwrite_phy(o_flx_co2_bb,    fco2_bb)

         ELSE IF (type_trac == 'inco') THEN
           itr = 0
           DO iq = 1, nqtot
             IF(.NOT.tracers(iq)%isInPhysics) CYCLE
             itr = itr+1
             IF(tracers(iq)%component /= 'co2i') CYCLE
             !--3D fields
             CALL histwrite_phy(o_trac   (itr),tr_seri(:,:,itr))
             CALL histwrite_phy(o_dtr_vdf(itr),d_tr_cl(:,:,itr))
             CALL histwrite_phy(o_dtr_the(itr),d_tr_th(:,:,itr))
             CALL histwrite_phy(o_dtr_con(itr),d_tr_cv(:,:,itr))
             !--2D fields
             !--CO2 burden
             zx_tmp_fi2d=0.
             IF (vars_defined) THEN
                DO k=1,klev
                   zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+zmasse(:,k)*tr_seri(:,k,itr)
                ENDDO
             ENDIF
             CALL histwrite_phy(o_trac_cum(itr), zx_tmp_fi2d)
           ENDDO !--iq
           !--CO2 net fluxes
           CALL histwrite_phy(o_flx_co2_land,  fco2_land)
           CALL histwrite_phy(o_flx_co2_ocean, fco2_ocean)
           CALL histwrite_phy(o_flx_co2_ocean_cor, fco2_ocean_cor)
           CALL histwrite_phy(o_flx_co2_land_cor, fco2_land_cor)
           CALL histwrite_phy(o_flx_co2_ff,    fco2_ff)
           CALL histwrite_phy(o_flx_co2_bb,    fco2_bb)

         ELSE IF (ANY(type_trac==['lmdz','coag'])) THEN
           itr = 0
           DO iq = 1, nqtot
             IF(.NOT.tracers(iq)%isInPhysics) CYCLE
             itr = itr + 1
!            write(*,*) 'phys_output_write_mod 2337: itr=',itr
             !--3D fields
             CALL histwrite_phy(o_trac(itr), tr_seri(:,:,itr))
             CALL histwrite_phy(o_dtr_vdf(itr),d_tr_cl(:,:,itr))
             CALL histwrite_phy(o_dtr_the(itr),d_tr_th(:,:,itr))
             CALL histwrite_phy(o_dtr_con(itr),d_tr_cv(:,:,itr))
             CALL histwrite_phy(o_dtr_lessi_impa(itr),d_tr_lessi_impa(:,:,itr))
             CALL histwrite_phy(o_dtr_lessi_nucl(itr),d_tr_lessi_nucl(:,:,itr))
             CALL histwrite_phy(o_dtr_insc(itr),d_tr_insc(:,:,itr))
             CALL histwrite_phy(o_dtr_bcscav(itr),d_tr_bcscav(:,:,itr))
             CALL histwrite_phy(o_dtr_evapls(itr),d_tr_evapls(:,:,itr))
             CALL histwrite_phy(o_dtr_ls(itr),d_tr_ls(:,:,itr))
             CALL histwrite_phy(o_dtr_trsp(itr),d_tr_trsp(:,:,itr))
             CALL histwrite_phy(o_dtr_sscav(itr),d_tr_sscav(:,:,itr))
             CALL histwrite_phy(o_dtr_sat(itr),d_tr_sat(:,:,itr))
             CALL histwrite_phy(o_dtr_uscav(itr),d_tr_uscav(:,:,itr))
            !--2D fields
             CALL histwrite_phy(o_dtr_dry(itr), flux_tr_dry(:,itr))
             zx_tmp_fi2d=0.
             IF (vars_defined) THEN
                DO k=1,klev
                   zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+zmasse(:,k)*tr_seri(:,k,itr)
                ENDDO
             ENDIF
             CALL histwrite_phy(o_trac_cum(itr), zx_tmp_fi2d)
           ENDDO !--iq
         ENDIF   !--type_trac
       ENDIF   !(iflag_phytrac==1)

       if (activate_ocean_skin >= 1) then
          CALL histwrite_phy(o_delta_sst, delta_sst)
          CALL histwrite_phy(o_delta_sal, delta_sal)
          CALL histwrite_phy(o_ds_ns, ds_ns)
          CALL histwrite_phy(o_dt_ns, dt_ns)
          CALL histwrite_phy(o_dter, dter)
          CALL histwrite_phy(o_dser, dser)
          CALL histwrite_phy(o_tkt, tkt)
          CALL histwrite_phy(o_tks, tks)
          CALL histwrite_phy(o_taur, taur)
          CALL histwrite_phy(o_sss, sss)
       end if

#ifdef ISO
    do ixt=1,ntiso
!        write(*,*) 'ixt'
        IF (vars_defined) zx_tmp_fi2d(:) = xtrain_fall(ixt,:) + xtsnow_fall(ixt,:)
        CALL histwrite_phy(o_xtprecip(ixt), zx_tmp_fi2d)

        IF (vars_defined) zx_tmp_fi2d(:) = xtrain_lsc(ixt,:) + xtsnow_lsc(ixt,:)
        CALL histwrite_phy(o_xtplul(ixt), zx_tmp_fi2d)

        IF (vars_defined) zx_tmp_fi2d(:) = xtrain_con(ixt,:) + xtsnow_con(ixt,:)
        CALL histwrite_phy(o_xtpluc(ixt), zx_tmp_fi2d)
        CALL histwrite_phy(o_xtevap(ixt),   xtevap(ixt,:))
        CALL histwrite_phy(o_xtovap(ixt),  xt_seri(ixt,:,:))
        CALL histwrite_phy(o_xtoliq(ixt), xtl_seri(ixt,:,:))

        DO nsrf = 1, nbsrf ! ajout Camille 8 mai 2023
        IF (vars_defined)       zx_tmp_fi2d(1 : klon) = fxtevap(ixt,:, nsrf)
        CALL histwrite_phy(o_xtevap_srf(ixt,nsrf), zx_tmp_fi2d)
        ENDDO

        IF (vars_defined) zx_tmp_fi3d(:,:)=xtl_seri(ixt,:,:)+xts_seri(ixt,:,:)
        CALL histwrite_phy(o_xtcond(ixt), zx_tmp_fi3d)
        CALL histwrite_phy(o_dxtdyn(ixt),   d_xt_dyn(ixt,:,:))
        CALL histwrite_phy(o_dxtldyn(ixt), d_xtl_dyn(ixt,:,:))

        IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_con(ixt,:,:)/pdtphys
        CALL histwrite_phy(o_dxtcon(ixt), zx_tmp_fi3d)

        IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_lsc(ixt,:,:)/pdtphys
        CALL histwrite_phy(o_dxtlsc(ixt), zx_tmp_fi3d)

        IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_eva(ixt,:,:)/pdtphys
        CALL histwrite_phy(o_dxteva(ixt), zx_tmp_fi3d)

        IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_vdf(ixt,:,:)/pdtphys
        CALL histwrite_phy(o_dxtvdf(ixt), zx_tmp_fi3d)

        IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_ajsb(ixt,:,:)/pdtphys
        CALL histwrite_phy(o_dxtajs(ixt), zx_tmp_fi3d)

        IF (vars_defined) zx_tmp_fi3d(:,:)=(d_xt_ajs(ixt,:,:)-d_xt_ajsb(ixt,:,:))/pdtphys
        CALL histwrite_phy(o_dxtthe(ixt), zx_tmp_fi3d)

        IF (ok_qch4) THEN
          IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_ch4(ixt,:,:)/pdtphys
          CALL histwrite_phy(o_dxtch4(ixt), zx_tmp_fi3d)
        END IF

        IF (ixt == iso_HTO) THEN
          IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_prod_nucl(ixt,:,:)/pdtphys
          CALL histwrite_phy(o_dxtprod_nucl(ixt), zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_cosmo(ixt,:,:)/pdtphys
          CALL histwrite_phy(o_dxtcosmo(ixt), zx_tmp_fi3d)

          IF (vars_defined) zx_tmp_fi3d(:,:)=d_xt_decroiss(ixt,:,:)/pdtphys
          CALL histwrite_phy(o_dxtdecroiss(ixt), zx_tmp_fi3d)
        END IF

    !write(*,*) 'phys_output_write_mod 2531'
    enddo
#endif

       IF (.NOT.vars_defined) THEN
          !$OMP MASTER
#ifndef CPP_IOIPSL_NO_OUTPUT 
          DO iff=1,nfiles
             IF (clef_files(iff)) THEN
                CALL histend(nid_files(iff))
                ndex2d = 0
                ndex3d = 0
             ENDIF ! clef_files
          ENDDO !  iff
#endif
          !On finalise l'initialisation:
          IF (using_xios) CALL wxios_closedef()

          !$OMP END MASTER
          !$OMP BARRIER
          vars_defined = .TRUE.

       ENDIF !--.NOT.vars_defined

    ENDDO

    IF (vars_defined) THEN
       ! On synchronise les fichiers pour IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
       !$OMP MASTER
       DO iff=1,nfiles
          IF (ok_sync .AND. clef_files(iff)) THEN
             CALL histsync(nid_files(iff))
          ENDIF
       END DO
       !$OMP END MASTER
#endif
    ENDIF

  END SUBROUTINE phys_output_write

END MODULE phys_output_write_mod
