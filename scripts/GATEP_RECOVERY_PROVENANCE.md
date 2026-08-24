# CASE-Q013 — Gate-P Globality Recovery Provenance

Frozen candidate: HYP-004C, N_ax=2, n=3, theta1=theta2=2.8, theta'=0, f1/f2 in [0,0.3], log10 zc1 in [3,3.75], log10 zc2 in [3.75,4.5].

Frozen theory: marcobeIIa/mAxiCLASS @ b5a3af9818d04d2fe696e54972545e8c6805f4ea (8 Nov 2025).
Frozen ACT: ACTCollaboration/DR6-ACT-lite @ 0e0cd2c703c62a0e980470b572602233b27750e1.
Frozen SPT: SouthPoleTelescope/spt_candl_data @ 2cec6e762a8c540484dd5acafc529f0035856350.
SPT likelihood: SPT3G_D1_TnE_lite, candl-like==2.0.3, clear_internal_priors=true.
Gate P: ACT DR6 TT/TE/EE ell 600..8500 + SPT-3G D1 TT/TE/EE + tau information. No lensing.
Executed nuisance contract retained: tau~N(0.051,0.006^2), A_act hard [0.5,1.5], P_act hard [0.9,1.1], A_act~N(1,0.003^2), Tcal hard [0.8,1.2], Ecal hard [0.8,1.2], Tcal~N(1,0.0036^2).

K1 Bella, Poulin, Vagnozzi & Knox, "Double the axions, half the tension: multi-field early dark energy eases the Hubble tension", arXiv:2604.13535 (2026).
K2 frozen mAxiCLASS commit above.
K3 frozen ACT DR6 commit above.
K4 act_dr6_cmbonly/tests/test_act.py, official software validation, reference lnL approximately -395.48.
K5 frozen SPT data commit above.
K6 frozen SPT repository publication state documenting candl 2.0.3.
K7 SPT3G_D1_TnE_lite.yaml: Gaussian 52 TT / 72 TE / 72 EE bandpowers and calibrations.
K8 cobaya/SPT3G_D1_TnE/cobaya_ttteee_lite.yaml: clear_internal_priors=true and external prior/bounds.
K9 E. Camphuis et al./SPT-3G, PRD 113, 083504 (2026), arXiv:2506.20707, DOI 10.1103/7wt3-9v2y.
K10 A. R. Khalife et al./SPT-3G, PRD 113, 103546 (2026), arXiv:2507.23355v2, DOI 10.1103/8jjr-7hpb.

R1 gatep-fit-lcdm-32764362888.zip: previous effective chi2 317.4238643933.
R2 gatep-fit-1ax-32764362888.zip: previous local effective chi2 319.1523776616.
R3 gatep-fit-2ax-32764362888.zip: previous local effective chi2 322.2027233657, f2=0.02544272654.
R4 gatep-final-32764362888.zip: both nesting checks false; NOT FINAL.
R5 gatep-act-dr6-data-32764362888.zip: ACT data artifact from successful direct run.

Recovery acceptance: chi2_2ax <= chi2_1ax <= chi2_LCDM within 0.1 numerical tolerance. Delta chi2_P=chi2_2ax-chi2_1ax. Delta AIC_P=Delta chi2_P+4. No Wilks/sigma conversion. Final Gate-P verdict still requires separate theory-precision gate |delta chi2_theory|<=0.1.
