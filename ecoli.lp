\ Model Ecoli_core_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Biomass_Ecoli_core_w_GAM
Subject To
 13dpg[c]: GAPD + PGK = 0
 2pg[c]: - ENO - PGM = 0
 3pg[c]: - 1.496 Biomass_Ecoli_core_w_GAM - PGK + PGM = 0
 6pgc[c]: - GND + PGL = 0
 6pgl[c]: G6PDH2r - PGL = 0
 ac[c]: - ACKr + ACt2r = 0
 ac[e]: - ACt2r - EX_ac(e) = 0
 acald[c]: - ACALD + ACALDt + ALCD2x = 0
 acald[e]: - ACALDt - EX_acald(e) = 0
 accoa[c]: ACALD - 3.7478 Biomass_Ecoli_core_w_GAM - CS - MALS + PDH + PFL
   - PTAr = 0
 acon_C[c]: ACONTa - ACONTb = 0
 actp[c]: ACKr + PTAr = 0
 adp[c]: ACKr + 2 ADK1 + ATPM - ATPS4r + 59.81 Biomass_Ecoli_core_w_GAM
   + GLNS + GLNabc + PFK + PGK + PPCK - PYK + SUCOAS = 0
 akg[c]: - AKGDH + AKGt2r + 4.1182 Biomass_Ecoli_core_w_GAM + GLUDy - GLUSy
   + ICDHyr = 0
 akg[e]: - AKGt2r - EX_akg(e) = 0
 amp[c]: - ADK1 + PPS = 0
 atp[c]: - ACKr - ADK1 - ATPM + ATPS4r - 59.81 Biomass_Ecoli_core_w_GAM
   - GLNS - GLNabc - PFK - PGK - PPCK - PPS + PYK - SUCOAS = 0
 cit[c]: - ACONTa + CS = 0
 co2[c]: AKGDH + CO2t + GND + ICDHyr + ME1 + ME2 + PDH - PPC + PPCK = 0
 co2[e]: - CO2t - EX_co2(e) = 0
 coa[c]: - ACALD - AKGDH + 3.7478 Biomass_Ecoli_core_w_GAM + CS + MALS
   - PDH - PFL + PTAr - SUCOAS = 0
 dhap[c]: FBA - TPI = 0
 e4p[c]: - 0.361 Biomass_Ecoli_core_w_GAM + TALA - TKT2 = 0
 etoh[c]: - ALCD2x + ETOHt2r = 0
 etoh[e]: - ETOHt2r - EX_etoh(e) = 0
 f6p[c]: - 0.0709 Biomass_Ecoli_core_w_GAM + FBP + FRUpts2 - PFK + PGI
   + TALA + TKT2 = 0
 fdp[c]: - FBA - FBP + PFK = 0
 for[c]: FORt2 - FORti + PFL = 0
 for[e]: - EX_for(e) - FORt2 + FORti = 0
 fru[e]: - EX_fru(e) - FRUpts2 = 0
 fum[c]: - FRD7 - FUM + FUMt2_2 + SUCDi = 0
 fum[e]: - EX_fum(e) - FUMt2_2 = 0
 g3p[c]: - 0.129 Biomass_Ecoli_core_w_GAM + FBA - GAPD - TALA + TKT1 + TKT2
   + TPI = 0
 g6p[c]: - 0.205 Biomass_Ecoli_core_w_GAM - G6PDH2r + GLCpts - PGI = 0
 glc_D[e]: - EX_glc(e) - GLCpts = 0
 gln_L[c]: - 0.2557 Biomass_Ecoli_core_w_GAM + GLNS + GLNabc - GLUN - GLUSy
   = 0
 gln_L[e]: - EX_gln_L(e) - GLNabc = 0
 glu_L[c]: - 4.9414 Biomass_Ecoli_core_w_GAM - GLNS - GLUDy + GLUN
   + 2 GLUSy + GLUt2r = 0
 glu_L[e]: - EX_glu_L(e) - GLUt2r = 0
 glx[c]: ICL - MALS = 0
 h2o[c]: ACONTa - ACONTb - ATPM + ATPS4r - 59.81 Biomass_Ecoli_core_w_GAM
   - CS + CYTBD + ENO - FBP - FUM - GLNabc - GLUDy - GLUN + H2Ot - MALS
   - PGL - PPC - PPS = 0
 h2o[e]: - EX_h2o(e) - H2Ot = 0
 h[c]: ACALD + ACt2r + AKGt2r + ALCD2x + ATPM + 3 ATPS4r
   + 59.81 Biomass_Ecoli_core_w_GAM + CS - 2 CYTBD + D_LACt2 + ETOHt2r
   + FORt2 + 2 FUMt2_2 + G6PDH2r + GAPD + GLNS + GLNabc + GLUDy - GLUSy
   + GLUt2r + LDH_D + MALS + 2 MALt2_2 + MDH - 4 NADH16 + PFK + PGL + PIt2r
   + PPC + 2 PPS - PYK + PYRt2r + 2 SUCCt2_2 + SUCCt3 + 2 THD2 = 0
 h[e]: - ACt2r - AKGt2r - 4 ATPS4r + 2 CYTBD - D_LACt2 - ETOHt2r - EX_h(e)
   - FORt2 - 2 FUMt2_2 - GLUt2r - 2 MALt2_2 + 3 NADH16 - PIt2r - PYRt2r
   - 2 SUCCt2_2 - SUCCt3 - 2 THD2 = 0
 icit[c]: ACONTb - ICDHyr - ICL = 0
 lac_D[c]: D_LACt2 - LDH_D = 0
 lac_D[e]: - D_LACt2 - EX_lac_D(e) = 0
 mal_L[c]: FUM + MALS + MALt2_2 - MDH - ME1 - ME2 = 0
 mal_L[e]: - EX_mal_L(e) - MALt2_2 = 0
 nad[c]: - ACALD - AKGDH - ALCD2x - 3.547 Biomass_Ecoli_core_w_GAM - GAPD
   - LDH_D - MDH - ME1 + NADH16 - NADTRHD - PDH + THD2 = 0
 nadh[c]: ACALD + AKGDH + ALCD2x + 3.547 Biomass_Ecoli_core_w_GAM + GAPD
   + LDH_D + MDH + ME1 - NADH16 + NADTRHD + PDH - THD2 = 0
 nadp[c]: 13.0279 Biomass_Ecoli_core_w_GAM - G6PDH2r - GLUDy + GLUSy - GND
   - ICDHyr - ME2 + NADTRHD - THD2 = 0
 nadph[c]: - 13.0279 Biomass_Ecoli_core_w_GAM + G6PDH2r + GLUDy - GLUSy
   + GND + ICDHyr + ME2 - NADTRHD + THD2 = 0
 nh4[c]: - GLNS + GLUDy + GLUN + NH4t = 0
 nh4[e]: - EX_nh4(e) - NH4t = 0
 o2[c]: - 0.5 CYTBD + O2t = 0
 o2[e]: - EX_o2(e) - O2t = 0
 oaa[c]: - 1.7867 Biomass_Ecoli_core_w_GAM - CS + MDH + PPC - PPCK = 0
 pep[c]: - 0.5191 Biomass_Ecoli_core_w_GAM + ENO - FRUpts2 - GLCpts - PPC
   + PPCK + PPS - PYK = 0
 pi[c]: ATPM - ATPS4r + 59.81 Biomass_Ecoli_core_w_GAM + FBP - GAPD + GLNS
   + GLNabc + PIt2r + PPC + PPS - PTAr + SUCOAS = 0
 pi[e]: - EX_pi(e) - PIt2r = 0
 pyr[c]: - 2.8328 Biomass_Ecoli_core_w_GAM + FRUpts2 + GLCpts + LDH_D + ME1
   + ME2 - PDH - PFL - PPS + PYK + PYRt2r = 0
 pyr[e]: - EX_pyr(e) - PYRt2r = 0
 q8[c]: CYTBD + FRD7 - NADH16 - SUCDi = 0
 q8h2[c]: - CYTBD - FRD7 + NADH16 + SUCDi = 0
 r5p[c]: - 0.8977 Biomass_Ecoli_core_w_GAM - RPI - TKT1 = 0
 ru5p_D[c]: GND - RPE + RPI = 0
 s7p[c]: - TALA + TKT1 = 0
 succ[c]: FRD7 + ICL + SUCCt2_2 - SUCCt3 - SUCDi - SUCOAS = 0
 succ[e]: - EX_succ(e) - SUCCt2_2 + SUCCt3 = 0
 succoa[c]: AKGDH + SUCOAS = 0
 xu5p_D[c]: RPE - TKT1 - TKT2 = 0
Bounds
 -1000 <= ACALD <= 1000
 -1000 <= ACALDt <= 1000
 -1000 <= ACKr <= 1000
 -1000 <= ACONTa <= 1000
 -1000 <= ACONTb <= 1000
 -1000 <= ACt2r <= 1000
 -1000 <= ADK1 <= 1000
 AKGDH <= 1000
 -1000 <= AKGt2r <= 1000
 -1000 <= ALCD2x <= 1000
 8.39 <= ATPM <= 1000
 -1000 <= ATPS4r <= 1000
 Biomass_Ecoli_core_w_GAM <= 1000
 -1000 <= CO2t <= 1000
 CS <= 1000
 CYTBD <= 1000
 -1000 <= D_LACt2 <= 1000
 -1000 <= ENO <= 1000
 -1000 <= ETOHt2r <= 1000
 EX_ac(e) <= 1000
 EX_acald(e) <= 1000
 EX_akg(e) <= 1000
 -1000 <= EX_co2(e) <= 1000
 EX_etoh(e) <= 1000
 EX_for(e) <= 1000
 EX_fru(e) <= 1000
 EX_fum(e) <= 1000
 -10 <= EX_glc(e) <= 1000
 EX_gln_L(e) <= 1000
 EX_glu_L(e) <= 1000
 -1000 <= EX_h(e) <= 1000
 -1000 <= EX_h2o(e) <= 1000
 EX_lac_D(e) <= 1000
 EX_mal_L(e) <= 1000
 -1000 <= EX_nh4(e) <= 1000
 -1000 <= EX_o2(e) <= 1000
 -1000 <= EX_pi(e) <= 1000
 EX_pyr(e) <= 1000
 EX_succ(e) <= 1000
 -1000 <= FBA <= 1000
 FBP <= 1000
 FORt2 <= 1000
 FORti <= 1000
 FRD7 <= 1000
 FRUpts2 <= 1000
 -1000 <= FUM <= 1000
 FUMt2_2 <= 1000
 -1000 <= G6PDH2r <= 1000
 -1000 <= GAPD <= 1000
 GLCpts <= 1000
 GLNS <= 1000
 GLNabc <= 1000
 -1000 <= GLUDy <= 1000
 GLUN <= 1000
 GLUSy <= 1000
 -1000 <= GLUt2r <= 1000
 GND <= 1000
 -1000 <= H2Ot <= 1000
 -1000 <= ICDHyr <= 1000
 ICL <= 1000
 -1000 <= LDH_D <= 1000
 MALS <= 1000
 MALt2_2 <= 1000
 -1000 <= MDH <= 1000
 ME1 <= 1000
 ME2 <= 1000
 NADH16 <= 1000
 NADTRHD <= 1000
 -1000 <= NH4t <= 1000
 -1000 <= O2t <= 1000
 PDH <= 1000
 PFK <= 1000
 PFL <= 1000
 -1000 <= PGI <= 1000
 -1000 <= PGK <= 1000
 PGL <= 1000
 -1000 <= PGM <= 1000
 -1000 <= PIt2r <= 1000
 PPC <= 1000
 PPCK <= 1000
 PPS <= 1000
 -1000 <= PTAr <= 1000
 PYK <= 1000
 -1000 <= PYRt2r <= 1000
 -1000 <= RPE <= 1000
 -1000 <= RPI <= 1000
 SUCCt2_2 <= 1000
 SUCCt3 <= 1000
 SUCDi <= 1000
 -1000 <= SUCOAS <= 1000
 -1000 <= TALA <= 1000
 THD2 <= 1000
 -1000 <= TKT1 <= 1000
 -1000 <= TKT2 <= 1000
 -1000 <= TPI <= 1000
End
