import re
from math import sqrt
from datetime import datetime as dt

from .INCA_i import inca_i
from .INCA_o import inca_o

class inca_api():
    '''class of methods to interact with the INCA matlab interface
    1. output to .m scripts
    2. input from .mat files
    3. internal methods for parsing and writing reactions
    '''
    def __init__(self):
        #modified biomass (INCA does not like exponential terms)
        self.biomass_INCA_iJS2012 = '0.488*ala_DASH_L_c + 0.0004*nadph_c + 0.176*phe_DASH_L_c + 0.582*gly_c + 0.00013*nadp_c + 0.276*ile_DASH_L_c + 0.21*pro_DASH_L_c + 0.25*glu_DASH_L_c + 0.154*glycogen_c + 45.73*atp_c + 0.139*utp_c + 0.00645*clpnEC_e + 0.131*tyr_DASH_L_c + 0.203*gtp_c + 0.00005*nadh_c + 0.0000006*coa_c + 0.00215*nad_c + 0.326*lys_DASH_L_c + 0.25*gln_DASH_L_c + 0.000003*succoa_c + 45.56*h2o_c + 0.205*ser_DASH_L_c + 0.126*ctp_c + 0.001*amp_c + 0.09675*peEC_e + 0.054*trp_DASH_L_c + 0.09*his_DASH_L_c + 0.0276*peptidoEC_e + 0.087*cys_DASH_L_c + 0.0084*lpsEC_e + 0.0247*datp_c + 0.0247*dttp_c + 0.241*thr_DASH_L_c + 0.281*arg_DASH_L_c + 0.00005*accoa_c + 0.402*val_DASH_L_c + 0.007*spmd_c + 0.0254*dgtp_c + 0.0232*pgEC_e + 0.146*met_DASH_L_c + 0.035*ptrc_c + 0.0254*dctp_c + 0.428*leu_DASH_L_c + 0.05*5mthf_c + 0.229*asp_DASH_L_c + 0.229*asn_DASH_L_c + 0.003*g1p_c + 0.0026*psEC_e + 0.00001*fad_c -> 45.56*pi_c + 45.56*h_c + 45.56*adp_c + 0.7332*ppi_c';
        self.biomass_INCA_iJS2012_v1 = '0.488*ala_DASH_L_c (C1:d C2:e C3:f) + 0.0004*nadph_c + 0.176*phe_DASH_L_c (C1:p C2:q C3:r C4:s C5:t C6:u C7:v C8:w C9:x) + 0.582*gly_c (C1:N C2:O) + 0.00013*nadp_c + 0.276*ile_DASH_L_c (C1:2 C2:3 C3:4 C4:5 C5:6 C6:7) + 0.21*pro_DASH_L_c (C1:y C2:z C3:A C4:B C5:C) + 0.25*glu_DASH_L_c (C1:I C2:J C3:K C4:L C5:M) + 0.154*glycogen_c (C1:P C2:Q C3:R C4:S C5:T C6:U) + 45.7318*atp_c + 0.139*utp_c + 0.00645*clpnEC_e + 0.131*tyr_DASH_L_c (C1:b C2:c C3:d C4:e C5:f C6:g C7:h C8:i C9:j) + 0.203*gtp_c + 0.00005*nadh_c + 0.000006*coa_c + 0.00215*nad_c + 0.326*lys_DASH_L_c (C1:e C2:f C3:g C4:h C5:i C6:j) + 0.25*gln_DASH_L_c (C1:D C2:E C3:F C4:G C5:H) + 0.000003*succoa_c (C1:R C2:S C3:T C4:U) + 45.5608*h2o_c + 0.205*ser_DASH_L_c (C1:H C2:I C3:J) + 0.126*ctp_c + 0.001*amp_c + 0.09675*peEC_e + 0.054*trp_DASH_L_c (C1:Z C2:1 C3:2 C4:3 C5:4 C6:5 C7:6 C8:7 C9:8 C10:9 C11:a) + 0.09*his_DASH_L_c (C1:V C2:W C3:X C4:Y C5:Z C6:1) + 0.0276*peptidoEC_e + 0.087*cys_DASH_L_c (C1:u C2:v C3:w) + 0.0084*lpsEC_e + 0.0247*datp_c + 0.0247*dttp_c + 0.241*thr_DASH_L_c (C1:V C2:W C3:X C4:Y) + 0.281*arg_DASH_L_c (C1:g C2:h C3:i C4:j C5:k C6:l) + 0.00005*accoa_c (C1:b C2:c) + 0.402*val_DASH_L_c (C1:k C2:l C3:m C4:n C5:o) + 0.007*spmd_c (C1:K C2:L C3:M C4:N C5:O C6:P C7:Q) + 0.0254*dgtp_c + 0.0232*pgEC_e + 0.146*met_DASH_L_c (C1:k C2:l C3:m C4:n C5:o) + 0.035*ptrc_c (C1:D C2:E C3:F C4:G) + 0.0254*dctp_c + 0.428*leu_DASH_L_c (C1:8 C2:9 C3:a C4:b C5:c C6:d) + 0.05*5mthf_c (C1:a) + 0.229*asp_DASH_L_c (C1:q C2:r C3:s C4:t) + 0.229*asn_DASH_L_c (C1:m C2:n C3:o C4:p) + 0.003*g1p_c (C1:x C2:y C3:z C4:A C5:B C6:C) + 0.0026*psEC_e + 0.00001*fad_c -> 45.5628*pi_c + 45.55735*h_c + 45.5608*adp_c + 0.7332*ppi_c ';
        self.biomass_INCA_iJS2012_v2 = '0.488*ala_DASH_L_c (C1:ala_DASH_L_cd C2:ala_DASH_L_ce C3:ala_DASH_L_cf) + 0.0004*nadph_c + 0.176*phe_DASH_L_c (C1:phe_DASH_L_cp C2:phe_DASH_L_cq C3:phe_DASH_L_cr C4:phe_DASH_L_cs C5:phe_DASH_L_ct C6:phe_DASH_L_cu C7:phe_DASH_L_cv C8:phe_DASH_L_cw C9:phe_DASH_L_cx) + 0.582*gly_c (C1:gly_cN C2:gly_cO) + 0.00013*nadp_c + 0.276*ile_DASH_L_c (C1:ile_DASH_L_c2 C2:ile_DASH_L_c3 C3:ile_DASH_L_c4 C4:ile_DASH_L_c5 C5:ile_DASH_L_c6 C6:ile_DASH_L_c7) + 0.21*pro_DASH_L_c (C1:pro_DASH_L_cy C2:pro_DASH_L_cz C3:pro_DASH_L_cA C4:pro_DASH_L_cB C5:pro_DASH_L_cC) + 0.25*glu_DASH_L_c (C1:glu_DASH_L_cI C2:glu_DASH_L_cJ C3:glu_DASH_L_cK C4:glu_DASH_L_cL C5:glu_DASH_L_cM) + 0.154*glycogen_c (C1:glycogen_cP C2:glycogen_cQ C3:glycogen_cR C4:glycogen_cS C5:glycogen_cT C6:glycogen_cU) + 45.73*atp_c + 0.139*utp_c + 0.00645*clpnEC_e + 0.131*tyr_DASH_L_c (C1:tyr_DASH_L_cb C2:tyr_DASH_L_cc C3:tyr_DASH_L_cd C4:tyr_DASH_L_ce C5:tyr_DASH_L_cf C6:tyr_DASH_L_cg C7:tyr_DASH_L_ch C8:tyr_DASH_L_ci C9:tyr_DASH_L_cj) + 0.203*gtp_c + 0.00005*nadh_c + 0.0000006*coa_c + 0.00215*nad_c + 0.326*lys_DASH_L_c (C1:lys_DASH_L_ce C2:lys_DASH_L_cf C3:lys_DASH_L_cg C4:lys_DASH_L_ch C5:lys_DASH_L_ci C6:lys_DASH_L_cj) + 0.25*gln_DASH_L_c (C1:gln_DASH_L_cD C2:gln_DASH_L_cE C3:gln_DASH_L_cF C4:gln_DASH_L_cG C5:gln_DASH_L_cH) + 0.000003*succoa_c (C1:succoa_cR C2:succoa_cS C3:succoa_cT C4:succoa_cU) + 45.56*h2o_c + 0.205*ser_DASH_L_c (C1:ser_DASH_L_cH C2:ser_DASH_L_cI C3:ser_DASH_L_cJ) + 0.126*ctp_c + 0.001*amp_c + 0.09675*peEC_e + 0.054*trp_DASH_L_c (C1:trp_DASH_L_cZ C2:trp_DASH_L_c1 C3:trp_DASH_L_c2 C4:trp_DASH_L_c3 C5:trp_DASH_L_c4 C6:trp_DASH_L_c5 C7:trp_DASH_L_c6 C8:trp_DASH_L_c7 C9:trp_DASH_L_c8 C10:trp_DASH_L_c9 C11:trp_DASH_L_ca) + 0.09*his_DASH_L_c (C1:his_DASH_L_cV C2:his_DASH_L_cW C3:his_DASH_L_cX C4:his_DASH_L_cY C5:his_DASH_L_cZ C6:his_DASH_L_c1) + 0.0276*peptidoEC_e + 0.087*cys_DASH_L_c (C1:cys_DASH_L_cu C2:cys_DASH_L_cv C3:cys_DASH_L_cw) + 0.0084*lpsEC_e + 0.0247*datp_c + 0.0247*dttp_c + 0.241*thr_DASH_L_c (C1:thr_DASH_L_cV C2:thr_DASH_L_cW C3:thr_DASH_L_cX C4:thr_DASH_L_cY) + 0.281*arg_DASH_L_c (C1:arg_DASH_L_cg C2:arg_DASH_L_ch C3:arg_DASH_L_ci C4:arg_DASH_L_cj C5:arg_DASH_L_ck C6:arg_DASH_L_cl) + 0.00005*accoa_c (C1:accoa_cb C2:accoa_cc) + 0.402*val_DASH_L_c (C1:val_DASH_L_ck C2:val_DASH_L_cl C3:val_DASH_L_cm C4:val_DASH_L_cn C5:val_DASH_L_co) + 0.007*spmd_c (C1:spmd_cK C2:spmd_cL C3:spmd_cM C4:spmd_cN C5:spmd_cO C6:spmd_cP C7:spmd_cQ) + 0.0254*dgtp_c + 0.0232*pgEC_e + 0.146*met_DASH_L_c (C1:met_DASH_L_ck C2:met_DASH_L_cl C3:met_DASH_L_cm C4:met_DASH_L_cn C5:met_DASH_L_co) + 0.035*ptrc_c (C1:ptrc_cD C2:ptrc_cE C3:ptrc_cF C4:ptrc_cG) + 0.0254*dctp_c + 0.428*leu_DASH_L_c (C1:leu_DASH_L_c8 C2:leu_DASH_L_c9 C3:leu_DASH_L_ca C4:leu_DASH_L_cb C5:leu_DASH_L_cc C6:leu_DASH_L_cd) + 0.05*5mthf_c (C1:5mthf_ca) + 0.229*asp_DASH_L_c (C1:asp_DASH_L_cq C2:asp_DASH_L_cr C3:asp_DASH_L_cs C4:asp_DASH_L_ct) + 0.229*asn_DASH_L_c (C1:asn_DASH_L_cm C2:asn_DASH_L_cn C3:asn_DASH_L_co C4:asn_DASH_L_cp) + 0.003*g1p_c (C1:g1p_cx C2:g1p_cy C3:g1p_cz C4:g1p_cA C5:g1p_cB C6:g1p_cC) + 0.0026*psEC_e + 0.00001*fad_c -> 45.56*pi_c + 45.56*h_c + 45.56*adp_c + 0.7332*ppi_c';
        self.biomass_INCA = '0.005707*pg160_c (C1:pg160_c0_C0 C2:pg160_c0_C1 C3:pg160_c0_C2 C4:pg160_c0_C3 C5:pg160_c0_C4 C6:pg160_c0_C5) + 0.000168*coa_c (C1:coa_c0_C0 C2:coa_c0_C1 C3:coa_c0_C2 C4:coa_c0_C3 C5:coa_c0_C4 C6:coa_c0_C5 C7:coa_c0_C6 C8:coa_c0_C7 C9:coa_c0_C8 C10:coa_c0_C9 C11:coa_c0_C10 C12:coa_c0_C11 C13:coa_c0_C12 C14:coa_c0_C13 C15:coa_c0_C14 C16:coa_c0_C15 C17:coa_c0_C16 C18:coa_c0_C17 C19:coa_c0_C18 C20:coa_c0_C19 C21:coa_c0_C20) + 0.000003*lipopb_c + 0.000307*ni2_c + 0.000055*udcpdp_c (C1:udcpdp_c0_C0 C2:udcpdp_c0_C1 C3:udcpdp_c0_C2 C4:udcpdp_c0_C3 C5:udcpdp_c0_C4 C6:udcpdp_c0_C5 C7:udcpdp_c0_C6 C8:udcpdp_c0_C7 C9:udcpdp_c0_C8 C10:udcpdp_c0_C9 C11:udcpdp_c0_C10 C12:udcpdp_c0_C11 C13:udcpdp_c0_C12 C14:udcpdp_c0_C13 C15:udcpdp_c0_C14 C16:udcpdp_c0_C15 C17:udcpdp_c0_C16 C18:udcpdp_c0_C17 C19:udcpdp_c0_C18 C20:udcpdp_c0_C19 C21:udcpdp_c0_C20 C22:udcpdp_c0_C21 C23:udcpdp_c0_C22 C24:udcpdp_c0_C23 C25:udcpdp_c0_C24 C26:udcpdp_c0_C25 C27:udcpdp_c0_C26 C28:udcpdp_c0_C27 C29:udcpdp_c0_C28 C30:udcpdp_c0_C29 C31:udcpdp_c0_C30 C32:udcpdp_c0_C31 C33:udcpdp_c0_C32 C34:udcpdp_c0_C33 C35:udcpdp_c0_C34 C36:udcpdp_c0_C35 C37:udcpdp_c0_C36 C38:udcpdp_c0_C37 C39:udcpdp_c0_C38 C40:udcpdp_c0_C39 C41:udcpdp_c0_C40 C42:udcpdp_c0_C41 C43:udcpdp_c0_C42 C44:udcpdp_c0_C43 C45:udcpdp_c0_C44 C46:udcpdp_c0_C45 C47:udcpdp_c0_C46 C48:udcpdp_c0_C47 C49:udcpdp_c0_C48 C50:udcpdp_c0_C49 C51:udcpdp_c0_C50 C52:udcpdp_c0_C51 C53:udcpdp_c0_C52 C54:udcpdp_c0_C53 C55:udcpdp_c0_C54) + 0.004957*pe181_c (C1:pe181_c0_C0 C2:pe181_c0_C1 C3:pe181_c0_C2 C4:pe181_c0_C3 C5:pe181_c0_C4) + 0.000112*nadp_c (C1:nadp_c0_C0 C2:nadp_c0_C1 C3:nadp_c0_C2 C4:nadp_c0_C3 C5:nadp_c0_C4 C6:nadp_c0_C5 C7:nadp_c0_C6 C8:nadp_c0_C7 C9:nadp_c0_C8 C10:nadp_c0_C9 C11:nadp_c0_C10 C12:nadp_c0_C11 C13:nadp_c0_C12 C14:nadp_c0_C13 C15:nadp_c0_C14 C16:nadp_c0_C15 C17:nadp_c0_C16 C18:nadp_c0_C17 C19:nadp_c0_C18 C20:nadp_c0_C19 C21:nadp_c0_C20) + 0.140101*utp_c (C1:utp_c0_C0 C2:utp_c0_C1 C3:utp_c0_C2 C4:utp_c0_C3 C5:utp_c0_C4 C6:utp_c0_C5 C7:utp_c0_C6 C8:utp_c0_C7 C9:utp_c0_C8) + 0.008253*mg2_c + 0.000024*cobalt2_c + 0.234232*asp_DASH_L_c (C1:asp_DASH_L_c0_C0 C2:asp_DASH_L_c0_C1 C3:asp_DASH_L_c0_C2 C4:asp_DASH_L_c0_C3) + 0.002288*pg181_c (C1:pg181_c0_C0 C2:pg181_c0_C1 C3:pg181_c0_C2 C4:pg181_c0_C3 C5:pg181_c0_C4 C6:pg181_c0_C5) + 0.154187*glycogen_c (C1:glycogen_c0_C0 C2:glycogen_c0_C1 C3:glycogen_c0_C2 C4:glycogen_c0_C3 C5:glycogen_c0_C4 C6:glycogen_c0_C5) + 0.000098*succoa_c (C1:succoa_c0_C0 C2:succoa_c0_C1 C3:succoa_c0_C2 C4:succoa_c0_C3 C5:succoa_c0_C4 C6:succoa_c0_C5 C7:succoa_c0_C6 C8:succoa_c0_C7 C9:succoa_c0_C8 C10:succoa_c0_C9 C11:succoa_c0_C10 C12:succoa_c0_C11 C13:succoa_c0_C12 C14:succoa_c0_C13 C15:succoa_c0_C14 C16:succoa_c0_C15 C17:succoa_c0_C16 C18:succoa_c0_C17 C19:succoa_c0_C18 C20:succoa_c0_C19 C21:succoa_c0_C20 C22:succoa_c0_C21 C23:succoa_c0_C22 C24:succoa_c0_C23 C25:succoa_c0_C24) + 48.752916*h2o_c + 0.031798*pe160_p (C1:pe160_p0_C0 C2:pe160_p0_C1 C3:pe160_p0_C2 C4:pe160_p0_C3 C5:pe160_p0_C4) + 0.000223*gthrd_c (C1:gthrd_c0_C0 C2:gthrd_c0_C1 C3:gthrd_c0_C2 C4:gthrd_c0_C3 C5:gthrd_c0_C4 C6:gthrd_c0_C5 C7:gthrd_c0_C6 C8:gthrd_c0_C7 C9:gthrd_c0_C8 C10:gthrd_c0_C9) + 0.000031*malcoa_c (C1:malcoa_c0_C0 C2:malcoa_c0_C1 C3:malcoa_c0_C2 C4:malcoa_c0_C3 C5:malcoa_c0_C4 C6:malcoa_c0_C5 C7:malcoa_c0_C6 C8:malcoa_c0_C7 C9:malcoa_c0_C8 C10:malcoa_c0_C9 C11:malcoa_c0_C10 C12:malcoa_c0_C11 C13:malcoa_c0_C12 C14:malcoa_c0_C13 C15:malcoa_c0_C14 C16:malcoa_c0_C15 C17:malcoa_c0_C16 C18:malcoa_c0_C17 C19:malcoa_c0_C18 C20:malcoa_c0_C19 C21:malcoa_c0_C20 C22:malcoa_c0_C21 C23:malcoa_c0_C22 C24:malcoa_c0_C23) + 0.209684*ser_DASH_L_c (C1:ser_DASH_L_c0_C0 C2:ser_DASH_L_c0_C1 C3:ser_DASH_L_c0_C2) + 0.234232*asn_DASH_L_c (C1:asn_DASH_L_c0_C0 C2:asn_DASH_L_c0_C1 C3:asn_DASH_L_c0_C2 C4:asn_DASH_L_c0_C3) + 0.000223*amet_c (C1:amet_c0_C0 C2:amet_c0_C1 C3:amet_c0_C2 C4:amet_c0_C3 C5:amet_c0_C4 C6:amet_c0_C5 C7:amet_c0_C6 C8:amet_c0_C7 C9:amet_c0_C8 C10:amet_c0_C9 C11:amet_c0_C10 C12:amet_c0_C11 C13:amet_c0_C12 C14:amet_c0_C13 C15:amet_c0_C14) + 0.595297*gly_c (C1:gly_c0_C0 C2:gly_c0_C1) + 0.000605*murein3px4p_p (C1:murein3px4p_p0_C0 C2:murein3px4p_p0_C1 C3:murein3px4p_p0_C2 C4:murein3px4p_p0_C3 C5:murein3px4p_p0_C4 C6:murein3px4p_p0_C5 C7:murein3px4p_p0_C6 C8:murein3px4p_p0_C7 C9:murein3px4p_p0_C8 C10:murein3px4p_p0_C9 C11:murein3px4p_p0_C10 C12:murein3px4p_p0_C11 C13:murein3px4p_p0_C12 C14:murein3px4p_p0_C13 C15:murein3px4p_p0_C14 C16:murein3px4p_p0_C15 C17:murein3px4p_p0_C16 C18:murein3px4p_p0_C17 C19:murein3px4p_p0_C18 C20:murein3px4p_p0_C19 C21:murein3px4p_p0_C20 C22:murein3px4p_p0_C21 C23:murein3px4p_p0_C22 C24:murein3px4p_p0_C23 C25:murein3px4p_p0_C24 C26:murein3px4p_p0_C25 C27:murein3px4p_p0_C26 C28:murein3px4p_p0_C27 C29:murein3px4p_p0_C28 C30:murein3px4p_p0_C29 C31:murein3px4p_p0_C30 C32:murein3px4p_p0_C31 C33:murein3px4p_p0_C32 C34:murein3px4p_p0_C33 C35:murein3px4p_p0_C34 C36:murein3px4p_p0_C35 C37:murein3px4p_p0_C36 C38:murein3px4p_p0_C37 C39:murein3px4p_p0_C38 C40:murein3px4p_p0_C39 C41:murein3px4p_p0_C40 C42:murein3px4p_p0_C41 C43:murein3px4p_p0_C42 C44:murein3px4p_p0_C43 C45:murein3px4p_p0_C44 C46:murein3px4p_p0_C45 C47:murein3px4p_p0_C46 C48:murein3px4p_p0_C47 C49:murein3px4p_p0_C48 C50:murein3px4p_p0_C49 C51:murein3px4p_p0_C50 C52:murein3px4p_p0_C51 C53:murein3px4p_p0_C52 C54:murein3px4p_p0_C53 C55:murein3px4p_p0_C54 C56:murein3px4p_p0_C55 C57:murein3px4p_p0_C56 C58:murein3px4p_p0_C57 C59:murein3px4p_p0_C58 C60:murein3px4p_p0_C59 C61:murein3px4p_p0_C60 C62:murein3px4p_p0_C61 C63:murein3px4p_p0_C62 C64:murein3px4p_p0_C63 C65:murein3px4p_p0_C64 C66:murein3px4p_p0_C65 C67:murein3px4p_p0_C66 C68:murein3px4p_p0_C67 C69:murein3px4p_p0_C68 C70:murein3px4p_p0_C69 C71:murein3px4p_p0_C70) + 0.055234*trp_DASH_L_c (C1:trp_DASH_L_c0_C0 C2:trp_DASH_L_c0_C1 C3:trp_DASH_L_c0_C2 C4:trp_DASH_L_c0_C3 C5:trp_DASH_L_c0_C4 C6:trp_DASH_L_c0_C5 C7:trp_DASH_L_c0_C6 C8:trp_DASH_L_c0_C7 C9:trp_DASH_L_c0_C8 C10:trp_DASH_L_c0_C9 C11:trp_DASH_L_c0_C10) + 0.03327*ptrc_c (C1:ptrc_c0_C0 C2:ptrc_c0_C1 C3:ptrc_c0_C2 C4:ptrc_c0_C3) + 0.006388*fe2_c + 0.000223*thf_c + 0.000007*mocogdp_c + 0.000223*fad_c (C1:fad_c0_C0 C2:fad_c0_C1 C3:fad_c0_C2 C4:fad_c0_C3 C5:fad_c0_C4 C6:fad_c0_C5 C7:fad_c0_C6 C8:fad_c0_C7 C9:fad_c0_C8 C10:fad_c0_C9 C11:fad_c0_C10 C12:fad_c0_C11 C13:fad_c0_C12 C14:fad_c0_C13 C15:fad_c0_C14 C16:fad_c0_C15 C17:fad_c0_C16 C18:fad_c0_C17 C19:fad_c0_C18 C20:fad_c0_C19 C21:fad_c0_C20 C22:fad_c0_C21 C23:fad_c0_C22 C24:fad_c0_C23 C25:fad_c0_C24 C26:fad_c0_C25 C27:fad_c0_C26) + 0.004126*so4_c + 0.411184*val_DASH_L_c (C1:val_DASH_L_c0_C0 C2:val_DASH_L_c0_C1 C3:val_DASH_L_c0_C2 C4:val_DASH_L_c0_C3 C5:val_DASH_L_c0_C4) + 0.18569*k_c + 0.005381*murein4p4p_p (C1:murein4p4p_p0_C0 C2:murein4p4p_p0_C1 C3:murein4p4p_p0_C2 C4:murein4p4p_p0_C3 C5:murein4p4p_p0_C4 C6:murein4p4p_p0_C5 C7:murein4p4p_p0_C6 C8:murein4p4p_p0_C7 C9:murein4p4p_p0_C8 C10:murein4p4p_p0_C9 C11:murein4p4p_p0_C10 C12:murein4p4p_p0_C11 C13:murein4p4p_p0_C12 C14:murein4p4p_p0_C13 C15:murein4p4p_p0_C14 C16:murein4p4p_p0_C15 C17:murein4p4p_p0_C16 C18:murein4p4p_p0_C17 C19:murein4p4p_p0_C18 C20:murein4p4p_p0_C19 C21:murein4p4p_p0_C20 C22:murein4p4p_p0_C21 C23:murein4p4p_p0_C22 C24:murein4p4p_p0_C23 C25:murein4p4p_p0_C24 C26:murein4p4p_p0_C25 C27:murein4p4p_p0_C26 C28:murein4p4p_p0_C27 C29:murein4p4p_p0_C28 C30:murein4p4p_p0_C29 C31:murein4p4p_p0_C30 C32:murein4p4p_p0_C31 C33:murein4p4p_p0_C32 C34:murein4p4p_p0_C33 C35:murein4p4p_p0_C34 C36:murein4p4p_p0_C35 C37:murein4p4p_p0_C36 C38:murein4p4p_p0_C37 C39:murein4p4p_p0_C38 C40:murein4p4p_p0_C39 C41:murein4p4p_p0_C40 C42:murein4p4p_p0_C41 C43:murein4p4p_p0_C42 C44:murein4p4p_p0_C43 C45:murein4p4p_p0_C44 C46:murein4p4p_p0_C45 C47:murein4p4p_p0_C46 C48:murein4p4p_p0_C47 C49:murein4p4p_p0_C48 C50:murein4p4p_p0_C49 C51:murein4p4p_p0_C50 C52:murein4p4p_p0_C51 C53:murein4p4p_p0_C52 C54:murein4p4p_p0_C53 C55:murein4p4p_p0_C54 C56:murein4p4p_p0_C55 C57:murein4p4p_p0_C56 C58:murein4p4p_p0_C57 C59:murein4p4p_p0_C58 C60:murein4p4p_p0_C59 C61:murein4p4p_p0_C60 C62:murein4p4p_p0_C61 C63:murein4p4p_p0_C62 C64:murein4p4p_p0_C63 C65:murein4p4p_p0_C64 C66:murein4p4p_p0_C65 C67:murein4p4p_p0_C66 C68:murein4p4p_p0_C67 C69:murein4p4p_p0_C68 C70:murein4p4p_p0_C69 C71:murein4p4p_p0_C70 C72:murein4p4p_p0_C71 C73:murein4p4p_p0_C72 C74:murein4p4p_p0_C73) + 0.000223*adocbl_c + 0.005448*murein4px4p_p (C1:murein4px4p_p0_C0 C2:murein4px4p_p0_C1 C3:murein4px4p_p0_C2 C4:murein4px4p_p0_C3 C5:murein4px4p_p0_C4 C6:murein4px4p_p0_C5 C7:murein4px4p_p0_C6 C8:murein4px4p_p0_C7 C9:murein4px4p_p0_C8 C10:murein4px4p_p0_C9 C11:murein4px4p_p0_C10 C12:murein4px4p_p0_C11 C13:murein4px4p_p0_C12 C14:murein4px4p_p0_C13 C15:murein4px4p_p0_C14 C16:murein4px4p_p0_C15 C17:murein4px4p_p0_C16 C18:murein4px4p_p0_C17 C19:murein4px4p_p0_C18 C20:murein4px4p_p0_C19 C21:murein4px4p_p0_C20 C22:murein4px4p_p0_C21 C23:murein4px4p_p0_C22 C24:murein4px4p_p0_C23 C25:murein4px4p_p0_C24 C26:murein4px4p_p0_C25 C27:murein4px4p_p0_C26 C28:murein4px4p_p0_C27 C29:murein4px4p_p0_C28 C30:murein4px4p_p0_C29 C31:murein4px4p_p0_C30 C32:murein4px4p_p0_C31 C33:murein4px4p_p0_C32 C34:murein4px4p_p0_C33 C35:murein4px4p_p0_C34 C36:murein4px4p_p0_C35 C37:murein4px4p_p0_C36 C38:murein4px4p_p0_C37 C39:murein4px4p_p0_C38 C40:murein4px4p_p0_C39 C41:murein4px4p_p0_C40 C42:murein4px4p_p0_C41 C43:murein4px4p_p0_C42 C44:murein4px4p_p0_C43 C45:murein4px4p_p0_C44 C46:murein4px4p_p0_C45 C47:murein4px4p_p0_C46 C48:murein4px4p_p0_C47 C49:murein4px4p_p0_C48 C50:murein4px4p_p0_C49 C51:murein4px4p_p0_C50 C52:murein4px4p_p0_C51 C53:murein4px4p_p0_C52 C54:murein4px4p_p0_C53 C55:murein4px4p_p0_C54 C56:murein4px4p_p0_C55 C57:murein4px4p_p0_C56 C58:murein4px4p_p0_C57 C59:murein4px4p_p0_C58 C60:murein4px4p_p0_C59 C61:murein4px4p_p0_C60 C62:murein4px4p_p0_C61 C63:murein4px4p_p0_C62 C64:murein4px4p_p0_C63 C65:murein4px4p_p0_C64 C66:murein4px4p_p0_C65 C67:murein4px4p_p0_C66 C68:murein4px4p_p0_C67 C69:murein4px4p_p0_C68 C70:murein4px4p_p0_C69 C71:murein4px4p_p0_C70 C72:murein4px4p_p0_C71 C73:murein4px4p_p0_C72 C74:murein4px4p_p0_C73) + 0.004952*ca2_c + 0.000025*2fe2s_c + 0.000335*nadph_c (C1:nadph_c0_C0 C2:nadph_c0_C1 C3:nadph_c0_C2 C4:nadph_c0_C3 C5:nadph_c0_C4 C6:nadph_c0_C5 C7:nadph_c0_C6 C8:nadph_c0_C7 C9:nadph_c0_C8 C10:nadph_c0_C9 C11:nadph_c0_C10 C12:nadph_c0_C11 C13:nadph_c0_C12 C14:nadph_c0_C13 C15:nadph_c0_C14 C16:nadph_c0_C15 C17:nadph_c0_C16 C18:nadph_c0_C17 C19:nadph_c0_C18 C20:nadph_c0_C19 C21:nadph_c0_C20) + 0.000045*nadh_c (C1:nadh_c0_C0 C2:nadh_c0_C1 C3:nadh_c0_C2 C4:nadh_c0_C3 C5:nadh_c0_C4 C6:nadh_c0_C5 C7:nadh_c0_C6 C8:nadh_c0_C7 C9:nadh_c0_C8 C10:nadh_c0_C9 C11:nadh_c0_C10 C12:nadh_c0_C11 C13:nadh_c0_C12 C14:nadh_c0_C13 C15:nadh_c0_C14 C16:nadh_c0_C15 C17:nadh_c0_C16 C18:nadh_c0_C17 C19:nadh_c0_C18 C20:nadh_c0_C19 C21:nadh_c0_C20) + 0.000674*cu2_c + 0.000007*mococdp_c + 0.000223*pheme_c + 0.004439*pg161_c (C1:pg161_c0_C0 C2:pg161_c0_C1 C3:pg161_c0_C2 C4:pg161_c0_C3 C5:pg161_c0_C4 C6:pg161_c0_C5) + 0.012747*pe181_p (C1:pe181_p0_C0 C2:pe181_p0_C1 C3:pe181_p0_C2 C4:pe181_p0_C3 C5:pe181_p0_C4) + 0.282306*ile_DASH_L_c (C1:ile_DASH_L_c0_C0 C2:ile_DASH_L_c0_C1 C3:ile_DASH_L_c0_C2 C4:ile_DASH_L_c0_C3 C5:ile_DASH_L_c0_C4 C6:ile_DASH_L_c0_C5) + 0.000223*chor_c (C1:chor_c0_C0 C2:chor_c0_C1 C3:chor_c0_C2 C4:chor_c0_C3 C5:chor_c0_C4 C6:chor_c0_C5 C7:chor_c0_C6 C8:chor_c0_C7 C9:chor_c0_C8 C10:chor_c0_C9) + 0.000223*q8h2_c + 0.008151*colipa_e + 0.333448*lys_DASH_L_c (C1:lys_DASH_L_c0_C0 C2:lys_DASH_L_c0_C1 C3:lys_DASH_L_c0_C2 C4:lys_DASH_L_c0_C3 C5:lys_DASH_L_c0_C4 C6:lys_DASH_L_c0_C5) + 0.000223*enter_c + 0.000223*mlthf_c (C1:mlthf_c0_C0) + 0.000223*thmpp_c + 0.28742*arg_DASH_L_c (C1:arg_DASH_L_c0_C0 C2:arg_DASH_L_c0_C1 C3:arg_DASH_L_c0_C2 C4:arg_DASH_L_c0_C3 C5:arg_DASH_L_c0_C4 C6:arg_DASH_L_c0_C5) + 0.000002*btn_c + 0.000223*hemeO_c + 0.499149*ala_DASH_L_c (C1:ala_DASH_L_c0_C0 C2:ala_DASH_L_c0_C1 C3:ala_DASH_L_c0_C2) + 0.246506*thr_DASH_L_c (C1:thr_DASH_L_c0_C0 C2:thr_DASH_L_c0_C1 C3:thr_DASH_L_c0_C2 C4:thr_DASH_L_c0_C3) + 0.088988*cys_DASH_L_c (C1:cys_DASH_L_c0_C0 C2:cys_DASH_L_c0_C1 C3:cys_DASH_L_c0_C2) + 0.001787*nad_c (C1:nad_c0_C0 C2:nad_c0_C1 C3:nad_c0_C2 C4:nad_c0_C3 C5:nad_c0_C4 C6:nad_c0_C5 C7:nad_c0_C6 C8:nad_c0_C7 C9:nad_c0_C8 C10:nad_c0_C9 C11:nad_c0_C10 C12:nad_c0_C11 C13:nad_c0_C12 C14:nad_c0_C13 C15:nad_c0_C14 C16:nad_c0_C15 C17:nad_c0_C16 C18:nad_c0_C17 C19:nad_c0_C18 C20:nad_c0_C19 C21:nad_c0_C20) + 0.180021*phe_DASH_L_c (C1:phe_DASH_L_c0_C0 C2:phe_DASH_L_c0_C1 C3:phe_DASH_L_c0_C2 C4:phe_DASH_L_c0_C3 C5:phe_DASH_L_c0_C4 C6:phe_DASH_L_c0_C5 C7:phe_DASH_L_c0_C6 C8:phe_DASH_L_c0_C7 C9:phe_DASH_L_c0_C8) + 0.025612*dctp_c (C1:dctp_c0_C0 C2:dctp_c0_C1 C3:dctp_c0_C2 C4:dctp_c0_C3 C5:dctp_c0_C4 C6:dctp_c0_C5 C7:dctp_c0_C6 C8:dctp_c0_C7 C9:dctp_c0_C8) + 0.149336*met_DASH_L_c (C1:met_DASH_L_c0_C0 C2:met_DASH_L_c0_C1 C3:met_DASH_L_c0_C2 C4:met_DASH_L_c0_C3 C5:met_DASH_L_c0_C4) + 0.012366*pe160_c (C1:pe160_c0_C0 C2:pe160_c0_C1 C3:pe160_c0_C2 C4:pe160_c0_C3 C5:pe160_c0_C4) + 0.209121*gtp_c (C1:gtp_c0_C0 C2:gtp_c0_C1 C3:gtp_c0_C2 C4:gtp_c0_C3 C5:gtp_c0_C4 C6:gtp_c0_C5 C7:gtp_c0_C6 C8:gtp_c0_C7 C9:gtp_c0_C8 C10:gtp_c0_C9) + 0.437778*leu_DASH_L_c (C1:leu_DASH_L_c0_C0 C2:leu_DASH_L_c0_C1 C3:leu_DASH_L_c0_C2 C4:leu_DASH_L_c0_C3 C5:leu_DASH_L_c0_C4 C6:leu_DASH_L_c0_C5) + 0.007428*fe3_c + 0.092056*his_DASH_L_c (C1:his_DASH_L_c0_C0 C2:his_DASH_L_c0_C1 C3:his_DASH_L_c0_C2 C4:his_DASH_L_c0_C3 C5:his_DASH_L_c0_C4 C6:his_DASH_L_c0_C5) + 0.009618*pe161_c (C1:pe161_c0_C0 C2:pe161_c0_C1 C3:pe161_c0_C2 C4:pe161_c0_C3 C5:pe161_c0_C4) + 0.000223*10fthf_c (C1:10fthf_c0_C0) + 0.024805*datp_c (C1:datp_c0_C0 C2:datp_c0_C1 C3:datp_c0_C2 C4:datp_c0_C3 C5:datp_c0_C4 C6:datp_c0_C5 C7:datp_c0_C6 C8:datp_c0_C7 C9:datp_c0_C8 C10:datp_c0_C9) + 0.000223*5mthf_c (C1:5mthf_c0_C0) + 0.000673*murein4px4px4p_p (C1:murein4px4px4p_p0_C0 C2:murein4px4px4p_p0_C1 C3:murein4px4px4p_p0_C2 C4:murein4px4px4p_p0_C3 C5:murein4px4px4p_p0_C4 C6:murein4px4px4p_p0_C5 C7:murein4px4px4p_p0_C6 C8:murein4px4px4p_p0_C7 C9:murein4px4px4p_p0_C8 C10:murein4px4px4p_p0_C9 C11:murein4px4px4p_p0_C10 C12:murein4px4px4p_p0_C11 C13:murein4px4px4p_p0_C12 C14:murein4px4px4p_p0_C13 C15:murein4px4px4p_p0_C14 C16:murein4px4px4p_p0_C15 C17:murein4px4px4p_p0_C16 C18:murein4px4px4p_p0_C17 C19:murein4px4px4p_p0_C18 C20:murein4px4px4p_p0_C19 C21:murein4px4px4p_p0_C20 C22:murein4px4px4p_p0_C21 C23:murein4px4px4p_p0_C22 C24:murein4px4px4p_p0_C23 C25:murein4px4px4p_p0_C24 C26:murein4px4px4p_p0_C25 C27:murein4px4px4p_p0_C26 C28:murein4px4px4p_p0_C27 C29:murein4px4px4p_p0_C28 C30:murein4px4px4p_p0_C29 C31:murein4px4px4p_p0_C30 C32:murein4px4px4p_p0_C31 C33:murein4px4px4p_p0_C32 C34:murein4px4px4p_p0_C33 C35:murein4px4px4p_p0_C34 C36:murein4px4px4p_p0_C35 C37:murein4px4px4p_p0_C36 C38:murein4px4px4p_p0_C37 C39:murein4px4px4p_p0_C38 C40:murein4px4px4p_p0_C39 C41:murein4px4px4p_p0_C40 C42:murein4px4px4p_p0_C41 C43:murein4px4px4p_p0_C42 C44:murein4px4px4p_p0_C43 C45:murein4px4px4p_p0_C44 C46:murein4px4px4p_p0_C45 C47:murein4px4px4p_p0_C46 C48:murein4px4px4p_p0_C47 C49:murein4px4px4p_p0_C48 C50:murein4px4px4p_p0_C49 C51:murein4px4px4p_p0_C50 C52:murein4px4px4p_p0_C51 C53:murein4px4px4p_p0_C52 C54:murein4px4px4p_p0_C53 C55:murein4px4px4p_p0_C54 C56:murein4px4px4p_p0_C55 C57:murein4px4px4p_p0_C56 C58:murein4px4px4p_p0_C57 C59:murein4px4px4p_p0_C58 C60:murein4px4px4p_p0_C59 C61:murein4px4px4p_p0_C60 C62:murein4px4px4p_p0_C61 C63:murein4px4px4p_p0_C62 C64:murein4px4px4p_p0_C63 C65:murein4px4px4p_p0_C64 C66:murein4px4px4p_p0_C65 C67:murein4px4px4p_p0_C66 C68:murein4px4px4p_p0_C67 C69:murein4px4px4p_p0_C68 C70:murein4px4px4p_p0_C69 C71:murein4px4px4p_p0_C70 C72:murein4px4px4p_p0_C71 C73:murein4px4px4p_p0_C72 C74:murein4px4px4p_p0_C73 C75:murein4px4px4p_p0_C74 C76:murein4px4px4p_p0_C75 C77:murein4px4px4p_p0_C76 C78:murein4px4px4p_p0_C77 C79:murein4px4px4p_p0_C78 C80:murein4px4px4p_p0_C79 C81:murein4px4px4p_p0_C80 C82:murein4px4px4p_p0_C81 C83:murein4px4px4p_p0_C82 C84:murein4px4px4p_p0_C83 C85:murein4px4px4p_p0_C84 C86:murein4px4px4p_p0_C85 C87:murein4px4px4p_p0_C86 C88:murein4px4px4p_p0_C87 C89:murein4px4px4p_p0_C88 C90:murein4px4px4p_p0_C89 C91:murein4px4px4p_p0_C90 C92:murein4px4px4p_p0_C91 C93:murein4px4px4p_p0_C92 C94:murein4px4px4p_p0_C93 C95:murein4px4px4p_p0_C94 C96:murein4px4px4p_p0_C95 C97:murein4px4px4p_p0_C96 C98:murein4px4px4p_p0_C97 C99:murein4px4px4p_p0_C98 C100:murein4px4px4p_p0_C99 C101:murein4px4px4p_p0_C100 C102:murein4px4px4p_p0_C101 C103:murein4px4px4p_p0_C102 C104:murein4px4px4p_p0_C103 C105:murein4px4px4p_p0_C104 C106:murein4px4px4p_p0_C105 C107:murein4px4px4p_p0_C106 C108:murein4px4px4p_p0_C107 C109:murein4px4px4p_p0_C108 C110:murein4px4px4p_p0_C109 C111:murein4px4px4p_p0_C110) + 0.024805*dttp_c (C1:dttp_c0_C0 C2:dttp_c0_C1 C3:dttp_c0_C2 C4:dttp_c0_C3 C5:dttp_c0_C4 C6:dttp_c0_C5 C7:dttp_c0_C6 C8:dttp_c0_C7 C9:dttp_c0_C8 C10:dttp_c0_C9) + 0.000223*ribflv_c (C1:ribflv_c0_C0 C2:ribflv_c0_C1 C3:ribflv_c0_C2 C4:ribflv_c0_C3 C5:ribflv_c0_C4 C6:ribflv_c0_C5 C7:ribflv_c0_C6 C8:ribflv_c0_C7 C9:ribflv_c0_C8 C10:ribflv_c0_C9 C11:ribflv_c0_C10 C12:ribflv_c0_C11 C13:ribflv_c0_C12 C14:ribflv_c0_C13 C15:ribflv_c0_C14 C16:ribflv_c0_C15 C17:ribflv_c0_C16) + 0.000223*pydx5p_c + 0.000324*zn2_c + 0.004952*cl_c + 0.000223*sheme_c + 0.001345*murein3p3p_p (C1:murein3p3p_p0_C0 C2:murein3p3p_p0_C1 C3:murein3p3p_p0_C2 C4:murein3p3p_p0_C3 C5:murein3p3p_p0_C4 C6:murein3p3p_p0_C5 C7:murein3p3p_p0_C6 C8:murein3p3p_p0_C7 C9:murein3p3p_p0_C8 C10:murein3p3p_p0_C9 C11:murein3p3p_p0_C10 C12:murein3p3p_p0_C11 C13:murein3p3p_p0_C12 C14:murein3p3p_p0_C13 C15:murein3p3p_p0_C14 C16:murein3p3p_p0_C15 C17:murein3p3p_p0_C16 C18:murein3p3p_p0_C17 C19:murein3p3p_p0_C18 C20:murein3p3p_p0_C19 C21:murein3p3p_p0_C20 C22:murein3p3p_p0_C21 C23:murein3p3p_p0_C22 C24:murein3p3p_p0_C23 C25:murein3p3p_p0_C24 C26:murein3p3p_p0_C25 C27:murein3p3p_p0_C26 C28:murein3p3p_p0_C27 C29:murein3p3p_p0_C28 C30:murein3p3p_p0_C29 C31:murein3p3p_p0_C30 C32:murein3p3p_p0_C31 C33:murein3p3p_p0_C32 C34:murein3p3p_p0_C33 C35:murein3p3p_p0_C34 C36:murein3p3p_p0_C35 C37:murein3p3p_p0_C36 C38:murein3p3p_p0_C37 C39:murein3p3p_p0_C38 C40:murein3p3p_p0_C39 C41:murein3p3p_p0_C40 C42:murein3p3p_p0_C41 C43:murein3p3p_p0_C42 C44:murein3p3p_p0_C43 C45:murein3p3p_p0_C44 C46:murein3p3p_p0_C45 C47:murein3p3p_p0_C46 C48:murein3p3p_p0_C47 C49:murein3p3p_p0_C48 C50:murein3p3p_p0_C49 C51:murein3p3p_p0_C50 C52:murein3p3p_p0_C51 C53:murein3p3p_p0_C52 C54:murein3p3p_p0_C53 C55:murein3p3p_p0_C54 C56:murein3p3p_p0_C55 C57:murein3p3p_p0_C56 C58:murein3p3p_p0_C57 C59:murein3p3p_p0_C58 C60:murein3p3p_p0_C59 C61:murein3p3p_p0_C60 C62:murein3p3p_p0_C61 C63:murein3p3p_p0_C62 C64:murein3p3p_p0_C63 C65:murein3p3p_p0_C64 C66:murein3p3p_p0_C65 C67:murein3p3p_p0_C66 C68:murein3p3p_p0_C67) + 0.004892*pg160_p (C1:pg160_p0_C0 C2:pg160_p0_C1 C3:pg160_p0_C2 C4:pg160_p0_C3 C5:pg160_p0_C4 C6:pg160_p0_C5) + 0.129799*ctp_c (C1:ctp_c0_C0 C2:ctp_c0_C1 C3:ctp_c0_C2 C4:ctp_c0_C3 C5:ctp_c0_C4 C6:ctp_c0_C5 C7:ctp_c0_C6 C8:ctp_c0_C7 C9:ctp_c0_C8) + 0.255712*glu_DASH_L_c (C1:glu_DASH_L_c0_C0 C2:glu_DASH_L_c0_C1 C3:glu_DASH_L_c0_C2 C4:glu_DASH_L_c0_C3 C5:glu_DASH_L_c0_C4) + 0.214798*pro_DASH_L_c (C1:pro_DASH_L_c0_C0 C2:pro_DASH_L_c0_C1 C3:pro_DASH_L_c0_C2 C4:pro_DASH_L_c0_C3 C5:pro_DASH_L_c0_C4) + 0.025612*dgtp_c (C1:dgtp_c0_C0 C2:dgtp_c0_C1 C3:dgtp_c0_C2 C4:dgtp_c0_C3 C5:dgtp_c0_C4 C6:dgtp_c0_C5 C7:dgtp_c0_C6 C8:dgtp_c0_C7 C9:dgtp_c0_C8 C10:dgtp_c0_C9) + 0.000007*mobd_c + 0.255712*gln_DASH_L_c (C1:gln_DASH_L_c0_C0 C2:gln_DASH_L_c0_C1 C3:gln_DASH_L_c0_C2 C4:gln_DASH_L_c0_C3 C5:gln_DASH_L_c0_C4) + 0.001961*pg181_p (C1:pg181_p0_C0 C2:pg181_p0_C1 C3:pg181_p0_C2 C4:pg181_p0_C3 C5:pg181_p0_C4 C6:pg181_p0_C5) + 0.000658*mn2_c + 0.000223*2dmmql8_c + 0.024732*pe161_p (C1:pe161_p0_C0 C2:pe161_p0_C1 C3:pe161_p0_C2 C4:pe161_p0_C3 C5:pe161_p0_C4) + 0.000248*4fe4s_c + 0.00118*clpn181_p (C1:clpn181_p0_C0 C2:clpn181_p0_C1 C3:clpn181_p0_C2 C4:clpn181_p0_C3 C5:clpn181_p0_C4 C6:clpn181_p0_C5 C7:clpn181_p0_C6 C8:clpn181_p0_C7 C9:clpn181_p0_C8) + 0.012379*nh4_c + 0.000223*mql8_c + 0.003805*pg161_p (C1:pg161_p0_C0 C2:pg161_p0_C1 C3:pg161_p0_C2 C4:pg161_p0_C3 C5:pg161_p0_C4 C6:pg161_p0_C5) + 0.000279*accoa_c (C1:accoa_c0_C0 C2:accoa_c0_C1 C3:accoa_c0_C2 C4:accoa_c0_C3 C5:accoa_c0_C4 C6:accoa_c0_C5 C7:accoa_c0_C6 C8:accoa_c0_C7 C9:accoa_c0_C8 C10:accoa_c0_C9 C11:accoa_c0_C10 C12:accoa_c0_C11 C13:accoa_c0_C12 C14:accoa_c0_C13 C15:accoa_c0_C14 C16:accoa_c0_C15 C17:accoa_c0_C16 C18:accoa_c0_C17 C19:accoa_c0_C18 C20:accoa_c0_C19 C21:accoa_c0_C20 C22:accoa_c0_C21 C23:accoa_c0_C22) + 54.119975*atp_c (C1:atp_c0_C0 C2:atp_c0_C1 C3:atp_c0_C2 C4:atp_c0_C3 C5:atp_c0_C4 C6:atp_c0_C5 C7:atp_c0_C6 C8:atp_c0_C7 C9:atp_c0_C8 C10:atp_c0_C9) + 0.133993*tyr_DASH_L_c (C1:tyr_DASH_L_c0_C0 C2:tyr_DASH_L_c0_C1 C3:tyr_DASH_L_c0_C2 C4:tyr_DASH_L_c0_C3 C5:tyr_DASH_L_c0_C4 C6:tyr_DASH_L_c0_C5 C7:tyr_DASH_L_c0_C6 C8:tyr_DASH_L_c0_C7 C9:tyr_DASH_L_c0_C8) + 0.006744*spmd_c (C1:spmd_c0_C0 C2:spmd_c0_C1 C3:spmd_c0_C2 C4:spmd_c0_C3 C5:spmd_c0_C4 C6:spmd_c0_C5 C7:spmd_c0_C6) + 0.002944*clpn160_p (C1:clpn160_p0_C0 C2:clpn160_p0_C1 C3:clpn160_p0_C2 C4:clpn160_p0_C3 C5:clpn160_p0_C4 C6:clpn160_p0_C5 C7:clpn160_p0_C6 C8:clpn160_p0_C7 C9:clpn160_p0_C8) + 0.000116*bmocogdp_c + 0.00229*clpn161_p (C1:clpn161_p0_C0 C2:clpn161_p0_C1 C3:clpn161_p0_C2 C4:clpn161_p0_C3 C5:clpn161_p0_C4 C6:clpn161_p0_C5 C7:clpn161_p0_C6 C8:clpn161_p0_C7 C9:clpn161_p0_C8) -> 0.749831*ppi_c + 53.95*adp_c (C1:atp_c0_C0 C2:atp_c0_C1 C3:atp_c0_C2 C4:atp_c0_C3 C5:atp_c0_C4 C6:atp_c0_C5 C7:atp_c0_C6 C8:atp_c0_C7 C9:atp_c0_C8 C10:atp_c0_C9) + 0.005707*Ec_biomass_iJO1366_WT_53p95M_pg160_c_0.balance (C1:pg160_c0_C0 C2:pg160_c0_C1 C3:pg160_c0_C2 C4:pg160_c0_C3 C5:pg160_c0_C4 C6:pg160_c0_C5) + 0.000168*Ec_biomass_iJO1366_WT_53p95M_coa_c_1.balance (C1:coa_c0_C0 C2:coa_c0_C1 C3:coa_c0_C2 C4:coa_c0_C3 C5:coa_c0_C4 C6:coa_c0_C5 C7:coa_c0_C6 C8:coa_c0_C7 C9:coa_c0_C8 C10:coa_c0_C9 C11:coa_c0_C10 C12:coa_c0_C11 C13:coa_c0_C12 C14:coa_c0_C13 C15:coa_c0_C14 C16:coa_c0_C15 C17:coa_c0_C16 C18:coa_c0_C17 C19:coa_c0_C18 C20:coa_c0_C19 C21:coa_c0_C20) + 0.000055*Ec_biomass_iJO1366_WT_53p95M_udcpdp_c_2.balance (C1:udcpdp_c0_C0 C2:udcpdp_c0_C1 C3:udcpdp_c0_C2 C4:udcpdp_c0_C3 C5:udcpdp_c0_C4 C6:udcpdp_c0_C5 C7:udcpdp_c0_C6 C8:udcpdp_c0_C7 C9:udcpdp_c0_C8 C10:udcpdp_c0_C9 C11:udcpdp_c0_C10 C12:udcpdp_c0_C11 C13:udcpdp_c0_C12 C14:udcpdp_c0_C13 C15:udcpdp_c0_C14 C16:udcpdp_c0_C15 C17:udcpdp_c0_C16 C18:udcpdp_c0_C17 C19:udcpdp_c0_C18 C20:udcpdp_c0_C19 C21:udcpdp_c0_C20 C22:udcpdp_c0_C21 C23:udcpdp_c0_C22 C24:udcpdp_c0_C23 C25:udcpdp_c0_C24 C26:udcpdp_c0_C25 C27:udcpdp_c0_C26 C28:udcpdp_c0_C27 C29:udcpdp_c0_C28 C30:udcpdp_c0_C29 C31:udcpdp_c0_C30 C32:udcpdp_c0_C31 C33:udcpdp_c0_C32 C34:udcpdp_c0_C33 C35:udcpdp_c0_C34 C36:udcpdp_c0_C35 C37:udcpdp_c0_C36 C38:udcpdp_c0_C37 C39:udcpdp_c0_C38 C40:udcpdp_c0_C39 C41:udcpdp_c0_C40 C42:udcpdp_c0_C41 C43:udcpdp_c0_C42 C44:udcpdp_c0_C43 C45:udcpdp_c0_C44 C46:udcpdp_c0_C45 C47:udcpdp_c0_C46 C48:udcpdp_c0_C47 C49:udcpdp_c0_C48 C50:udcpdp_c0_C49 C51:udcpdp_c0_C50 C52:udcpdp_c0_C51 C53:udcpdp_c0_C52 C54:udcpdp_c0_C53 C55:udcpdp_c0_C54) + 0.004957*Ec_biomass_iJO1366_WT_53p95M_pe181_c_3.balance (C1:pe181_c0_C0 C2:pe181_c0_C1 C3:pe181_c0_C2 C4:pe181_c0_C3 C5:pe181_c0_C4) + 0.000112*Ec_biomass_iJO1366_WT_53p95M_nadp_c_4.balance (C1:nadp_c0_C0 C2:nadp_c0_C1 C3:nadp_c0_C2 C4:nadp_c0_C3 C5:nadp_c0_C4 C6:nadp_c0_C5 C7:nadp_c0_C6 C8:nadp_c0_C7 C9:nadp_c0_C8 C10:nadp_c0_C9 C11:nadp_c0_C10 C12:nadp_c0_C11 C13:nadp_c0_C12 C14:nadp_c0_C13 C15:nadp_c0_C14 C16:nadp_c0_C15 C17:nadp_c0_C16 C18:nadp_c0_C17 C19:nadp_c0_C18 C20:nadp_c0_C19 C21:nadp_c0_C20) + 0.140101*Ec_biomass_iJO1366_WT_53p95M_utp_c_5.balance (C1:utp_c0_C0 C2:utp_c0_C1 C3:utp_c0_C2 C4:utp_c0_C3 C5:utp_c0_C4 C6:utp_c0_C5 C7:utp_c0_C6 C8:utp_c0_C7 C9:utp_c0_C8) + 0.234232*Ec_biomass_iJO1366_WT_53p95M_asp_DASH_L_c_6.balance (C1:asp_DASH_L_c0_C0 C2:asp_DASH_L_c0_C1 C3:asp_DASH_L_c0_C2 C4:asp_DASH_L_c0_C3) + 0.002288*Ec_biomass_iJO1366_WT_53p95M_pg181_c_7.balance (C1:pg181_c0_C0 C2:pg181_c0_C1 C3:pg181_c0_C2 C4:pg181_c0_C3 C5:pg181_c0_C4 C6:pg181_c0_C5) + 0.154187*Ec_biomass_iJO1366_WT_53p95M_glycogen_c_8.balance (C1:glycogen_c0_C0 C2:glycogen_c0_C1 C3:glycogen_c0_C2 C4:glycogen_c0_C3 C5:glycogen_c0_C4 C6:glycogen_c0_C5) + 0.000098*Ec_biomass_iJO1366_WT_53p95M_succoa_c_9.balance (C1:succoa_c0_C0 C2:succoa_c0_C1 C3:succoa_c0_C2 C4:succoa_c0_C3 C5:succoa_c0_C4 C6:succoa_c0_C5 C7:succoa_c0_C6 C8:succoa_c0_C7 C9:succoa_c0_C8 C10:succoa_c0_C9 C11:succoa_c0_C10 C12:succoa_c0_C11 C13:succoa_c0_C12 C14:succoa_c0_C13 C15:succoa_c0_C14 C16:succoa_c0_C15 C17:succoa_c0_C16 C18:succoa_c0_C17 C19:succoa_c0_C18 C20:succoa_c0_C19 C21:succoa_c0_C20 C22:succoa_c0_C21 C23:succoa_c0_C22 C24:succoa_c0_C23 C25:succoa_c0_C24) + 0.031798*Ec_biomass_iJO1366_WT_53p95M_pe160_p_10.balance (C1:pe160_p0_C0 C2:pe160_p0_C1 C3:pe160_p0_C2 C4:pe160_p0_C3 C5:pe160_p0_C4) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_gthrd_c_11.balance (C1:gthrd_c0_C0 C2:gthrd_c0_C1 C3:gthrd_c0_C2 C4:gthrd_c0_C3 C5:gthrd_c0_C4 C6:gthrd_c0_C5 C7:gthrd_c0_C6 C8:gthrd_c0_C7 C9:gthrd_c0_C8 C10:gthrd_c0_C9) + 0.000031*Ec_biomass_iJO1366_WT_53p95M_malcoa_c_12.balance (C1:malcoa_c0_C0 C2:malcoa_c0_C1 C3:malcoa_c0_C2 C4:malcoa_c0_C3 C5:malcoa_c0_C4 C6:malcoa_c0_C5 C7:malcoa_c0_C6 C8:malcoa_c0_C7 C9:malcoa_c0_C8 C10:malcoa_c0_C9 C11:malcoa_c0_C10 C12:malcoa_c0_C11 C13:malcoa_c0_C12 C14:malcoa_c0_C13 C15:malcoa_c0_C14 C16:malcoa_c0_C15 C17:malcoa_c0_C16 C18:malcoa_c0_C17 C19:malcoa_c0_C18 C20:malcoa_c0_C19 C21:malcoa_c0_C20 C22:malcoa_c0_C21 C23:malcoa_c0_C22 C24:malcoa_c0_C23) + 0.209684*Ec_biomass_iJO1366_WT_53p95M_ser_DASH_L_c_13.balance (C1:ser_DASH_L_c0_C0 C2:ser_DASH_L_c0_C1 C3:ser_DASH_L_c0_C2) + 0.234232*Ec_biomass_iJO1366_WT_53p95M_asn_DASH_L_c_14.balance (C1:asn_DASH_L_c0_C0 C2:asn_DASH_L_c0_C1 C3:asn_DASH_L_c0_C2 C4:asn_DASH_L_c0_C3) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_amet_c_15.balance (C1:amet_c0_C0 C2:amet_c0_C1 C3:amet_c0_C2 C4:amet_c0_C3 C5:amet_c0_C4 C6:amet_c0_C5 C7:amet_c0_C6 C8:amet_c0_C7 C9:amet_c0_C8 C10:amet_c0_C9 C11:amet_c0_C10 C12:amet_c0_C11 C13:amet_c0_C12 C14:amet_c0_C13 C15:amet_c0_C14) + 0.595297*Ec_biomass_iJO1366_WT_53p95M_gly_c_16.balance (C1:gly_c0_C0 C2:gly_c0_C1) + 0.000605*Ec_biomass_iJO1366_WT_53p95M_murein3px4p_p_17.balance (C1:murein3px4p_p0_C0 C2:murein3px4p_p0_C1 C3:murein3px4p_p0_C2 C4:murein3px4p_p0_C3 C5:murein3px4p_p0_C4 C6:murein3px4p_p0_C5 C7:murein3px4p_p0_C6 C8:murein3px4p_p0_C7 C9:murein3px4p_p0_C8 C10:murein3px4p_p0_C9 C11:murein3px4p_p0_C10 C12:murein3px4p_p0_C11 C13:murein3px4p_p0_C12 C14:murein3px4p_p0_C13 C15:murein3px4p_p0_C14 C16:murein3px4p_p0_C15 C17:murein3px4p_p0_C16 C18:murein3px4p_p0_C17 C19:murein3px4p_p0_C18 C20:murein3px4p_p0_C19 C21:murein3px4p_p0_C20 C22:murein3px4p_p0_C21 C23:murein3px4p_p0_C22 C24:murein3px4p_p0_C23 C25:murein3px4p_p0_C24 C26:murein3px4p_p0_C25 C27:murein3px4p_p0_C26 C28:murein3px4p_p0_C27 C29:murein3px4p_p0_C28 C30:murein3px4p_p0_C29 C31:murein3px4p_p0_C30 C32:murein3px4p_p0_C31 C33:murein3px4p_p0_C32 C34:murein3px4p_p0_C33 C35:murein3px4p_p0_C34 C36:murein3px4p_p0_C35 C37:murein3px4p_p0_C36 C38:murein3px4p_p0_C37 C39:murein3px4p_p0_C38 C40:murein3px4p_p0_C39 C41:murein3px4p_p0_C40 C42:murein3px4p_p0_C41 C43:murein3px4p_p0_C42 C44:murein3px4p_p0_C43 C45:murein3px4p_p0_C44 C46:murein3px4p_p0_C45 C47:murein3px4p_p0_C46 C48:murein3px4p_p0_C47 C49:murein3px4p_p0_C48 C50:murein3px4p_p0_C49 C51:murein3px4p_p0_C50 C52:murein3px4p_p0_C51 C53:murein3px4p_p0_C52 C54:murein3px4p_p0_C53 C55:murein3px4p_p0_C54 C56:murein3px4p_p0_C55 C57:murein3px4p_p0_C56 C58:murein3px4p_p0_C57 C59:murein3px4p_p0_C58 C60:murein3px4p_p0_C59 C61:murein3px4p_p0_C60 C62:murein3px4p_p0_C61 C63:murein3px4p_p0_C62 C64:murein3px4p_p0_C63 C65:murein3px4p_p0_C64 C66:murein3px4p_p0_C65 C67:murein3px4p_p0_C66 C68:murein3px4p_p0_C67 C69:murein3px4p_p0_C68 C70:murein3px4p_p0_C69 C71:murein3px4p_p0_C70) + 0.055234*Ec_biomass_iJO1366_WT_53p95M_trp_DASH_L_c_18.balance (C1:trp_DASH_L_c0_C0 C2:trp_DASH_L_c0_C1 C3:trp_DASH_L_c0_C2 C4:trp_DASH_L_c0_C3 C5:trp_DASH_L_c0_C4 C6:trp_DASH_L_c0_C5 C7:trp_DASH_L_c0_C6 C8:trp_DASH_L_c0_C7 C9:trp_DASH_L_c0_C8 C10:trp_DASH_L_c0_C9 C11:trp_DASH_L_c0_C10) + 0.03327*Ec_biomass_iJO1366_WT_53p95M_ptrc_c_19.balance (C1:ptrc_c0_C0 C2:ptrc_c0_C1 C3:ptrc_c0_C2 C4:ptrc_c0_C3) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_fad_c_20.balance (C1:fad_c0_C0 C2:fad_c0_C1 C3:fad_c0_C2 C4:fad_c0_C3 C5:fad_c0_C4 C6:fad_c0_C5 C7:fad_c0_C6 C8:fad_c0_C7 C9:fad_c0_C8 C10:fad_c0_C9 C11:fad_c0_C10 C12:fad_c0_C11 C13:fad_c0_C12 C14:fad_c0_C13 C15:fad_c0_C14 C16:fad_c0_C15 C17:fad_c0_C16 C18:fad_c0_C17 C19:fad_c0_C18 C20:fad_c0_C19 C21:fad_c0_C20 C22:fad_c0_C21 C23:fad_c0_C22 C24:fad_c0_C23 C25:fad_c0_C24 C26:fad_c0_C25 C27:fad_c0_C26) + 0.411184*Ec_biomass_iJO1366_WT_53p95M_val_DASH_L_c_21.balance (C1:val_DASH_L_c0_C0 C2:val_DASH_L_c0_C1 C3:val_DASH_L_c0_C2 C4:val_DASH_L_c0_C3 C5:val_DASH_L_c0_C4) + 0.005381*Ec_biomass_iJO1366_WT_53p95M_murein4p4p_p_22.balance (C1:murein4p4p_p0_C0 C2:murein4p4p_p0_C1 C3:murein4p4p_p0_C2 C4:murein4p4p_p0_C3 C5:murein4p4p_p0_C4 C6:murein4p4p_p0_C5 C7:murein4p4p_p0_C6 C8:murein4p4p_p0_C7 C9:murein4p4p_p0_C8 C10:murein4p4p_p0_C9 C11:murein4p4p_p0_C10 C12:murein4p4p_p0_C11 C13:murein4p4p_p0_C12 C14:murein4p4p_p0_C13 C15:murein4p4p_p0_C14 C16:murein4p4p_p0_C15 C17:murein4p4p_p0_C16 C18:murein4p4p_p0_C17 C19:murein4p4p_p0_C18 C20:murein4p4p_p0_C19 C21:murein4p4p_p0_C20 C22:murein4p4p_p0_C21 C23:murein4p4p_p0_C22 C24:murein4p4p_p0_C23 C25:murein4p4p_p0_C24 C26:murein4p4p_p0_C25 C27:murein4p4p_p0_C26 C28:murein4p4p_p0_C27 C29:murein4p4p_p0_C28 C30:murein4p4p_p0_C29 C31:murein4p4p_p0_C30 C32:murein4p4p_p0_C31 C33:murein4p4p_p0_C32 C34:murein4p4p_p0_C33 C35:murein4p4p_p0_C34 C36:murein4p4p_p0_C35 C37:murein4p4p_p0_C36 C38:murein4p4p_p0_C37 C39:murein4p4p_p0_C38 C40:murein4p4p_p0_C39 C41:murein4p4p_p0_C40 C42:murein4p4p_p0_C41 C43:murein4p4p_p0_C42 C44:murein4p4p_p0_C43 C45:murein4p4p_p0_C44 C46:murein4p4p_p0_C45 C47:murein4p4p_p0_C46 C48:murein4p4p_p0_C47 C49:murein4p4p_p0_C48 C50:murein4p4p_p0_C49 C51:murein4p4p_p0_C50 C52:murein4p4p_p0_C51 C53:murein4p4p_p0_C52 C54:murein4p4p_p0_C53 C55:murein4p4p_p0_C54 C56:murein4p4p_p0_C55 C57:murein4p4p_p0_C56 C58:murein4p4p_p0_C57 C59:murein4p4p_p0_C58 C60:murein4p4p_p0_C59 C61:murein4p4p_p0_C60 C62:murein4p4p_p0_C61 C63:murein4p4p_p0_C62 C64:murein4p4p_p0_C63 C65:murein4p4p_p0_C64 C66:murein4p4p_p0_C65 C67:murein4p4p_p0_C66 C68:murein4p4p_p0_C67 C69:murein4p4p_p0_C68 C70:murein4p4p_p0_C69 C71:murein4p4p_p0_C70 C72:murein4p4p_p0_C71 C73:murein4p4p_p0_C72 C74:murein4p4p_p0_C73) + 0.005448*Ec_biomass_iJO1366_WT_53p95M_murein4px4p_p_23.balance (C1:murein4px4p_p0_C0 C2:murein4px4p_p0_C1 C3:murein4px4p_p0_C2 C4:murein4px4p_p0_C3 C5:murein4px4p_p0_C4 C6:murein4px4p_p0_C5 C7:murein4px4p_p0_C6 C8:murein4px4p_p0_C7 C9:murein4px4p_p0_C8 C10:murein4px4p_p0_C9 C11:murein4px4p_p0_C10 C12:murein4px4p_p0_C11 C13:murein4px4p_p0_C12 C14:murein4px4p_p0_C13 C15:murein4px4p_p0_C14 C16:murein4px4p_p0_C15 C17:murein4px4p_p0_C16 C18:murein4px4p_p0_C17 C19:murein4px4p_p0_C18 C20:murein4px4p_p0_C19 C21:murein4px4p_p0_C20 C22:murein4px4p_p0_C21 C23:murein4px4p_p0_C22 C24:murein4px4p_p0_C23 C25:murein4px4p_p0_C24 C26:murein4px4p_p0_C25 C27:murein4px4p_p0_C26 C28:murein4px4p_p0_C27 C29:murein4px4p_p0_C28 C30:murein4px4p_p0_C29 C31:murein4px4p_p0_C30 C32:murein4px4p_p0_C31 C33:murein4px4p_p0_C32 C34:murein4px4p_p0_C33 C35:murein4px4p_p0_C34 C36:murein4px4p_p0_C35 C37:murein4px4p_p0_C36 C38:murein4px4p_p0_C37 C39:murein4px4p_p0_C38 C40:murein4px4p_p0_C39 C41:murein4px4p_p0_C40 C42:murein4px4p_p0_C41 C43:murein4px4p_p0_C42 C44:murein4px4p_p0_C43 C45:murein4px4p_p0_C44 C46:murein4px4p_p0_C45 C47:murein4px4p_p0_C46 C48:murein4px4p_p0_C47 C49:murein4px4p_p0_C48 C50:murein4px4p_p0_C49 C51:murein4px4p_p0_C50 C52:murein4px4p_p0_C51 C53:murein4px4p_p0_C52 C54:murein4px4p_p0_C53 C55:murein4px4p_p0_C54 C56:murein4px4p_p0_C55 C57:murein4px4p_p0_C56 C58:murein4px4p_p0_C57 C59:murein4px4p_p0_C58 C60:murein4px4p_p0_C59 C61:murein4px4p_p0_C60 C62:murein4px4p_p0_C61 C63:murein4px4p_p0_C62 C64:murein4px4p_p0_C63 C65:murein4px4p_p0_C64 C66:murein4px4p_p0_C65 C67:murein4px4p_p0_C66 C68:murein4px4p_p0_C67 C69:murein4px4p_p0_C68 C70:murein4px4p_p0_C69 C71:murein4px4p_p0_C70 C72:murein4px4p_p0_C71 C73:murein4px4p_p0_C72 C74:murein4px4p_p0_C73) + 0.000335*Ec_biomass_iJO1366_WT_53p95M_nadph_c_24.balance (C1:nadph_c0_C0 C2:nadph_c0_C1 C3:nadph_c0_C2 C4:nadph_c0_C3 C5:nadph_c0_C4 C6:nadph_c0_C5 C7:nadph_c0_C6 C8:nadph_c0_C7 C9:nadph_c0_C8 C10:nadph_c0_C9 C11:nadph_c0_C10 C12:nadph_c0_C11 C13:nadph_c0_C12 C14:nadph_c0_C13 C15:nadph_c0_C14 C16:nadph_c0_C15 C17:nadph_c0_C16 C18:nadph_c0_C17 C19:nadph_c0_C18 C20:nadph_c0_C19 C21:nadph_c0_C20) + 0.000045*Ec_biomass_iJO1366_WT_53p95M_nadh_c_25.balance (C1:nadh_c0_C0 C2:nadh_c0_C1 C3:nadh_c0_C2 C4:nadh_c0_C3 C5:nadh_c0_C4 C6:nadh_c0_C5 C7:nadh_c0_C6 C8:nadh_c0_C7 C9:nadh_c0_C8 C10:nadh_c0_C9 C11:nadh_c0_C10 C12:nadh_c0_C11 C13:nadh_c0_C12 C14:nadh_c0_C13 C15:nadh_c0_C14 C16:nadh_c0_C15 C17:nadh_c0_C16 C18:nadh_c0_C17 C19:nadh_c0_C18 C20:nadh_c0_C19 C21:nadh_c0_C20) + 0.004439*Ec_biomass_iJO1366_WT_53p95M_pg161_c_26.balance (C1:pg161_c0_C0 C2:pg161_c0_C1 C3:pg161_c0_C2 C4:pg161_c0_C3 C5:pg161_c0_C4 C6:pg161_c0_C5) + 0.012747*Ec_biomass_iJO1366_WT_53p95M_pe181_p_27.balance (C1:pe181_p0_C0 C2:pe181_p0_C1 C3:pe181_p0_C2 C4:pe181_p0_C3 C5:pe181_p0_C4) + 0.282306*Ec_biomass_iJO1366_WT_53p95M_ile_DASH_L_c_28.balance (C1:ile_DASH_L_c0_C0 C2:ile_DASH_L_c0_C1 C3:ile_DASH_L_c0_C2 C4:ile_DASH_L_c0_C3 C5:ile_DASH_L_c0_C4 C6:ile_DASH_L_c0_C5) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_chor_c_29.balance (C1:chor_c0_C0 C2:chor_c0_C1 C3:chor_c0_C2 C4:chor_c0_C3 C5:chor_c0_C4 C6:chor_c0_C5 C7:chor_c0_C6 C8:chor_c0_C7 C9:chor_c0_C8 C10:chor_c0_C9) + 0.333448*Ec_biomass_iJO1366_WT_53p95M_lys_DASH_L_c_30.balance (C1:lys_DASH_L_c0_C0 C2:lys_DASH_L_c0_C1 C3:lys_DASH_L_c0_C2 C4:lys_DASH_L_c0_C3 C5:lys_DASH_L_c0_C4 C6:lys_DASH_L_c0_C5) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_mlthf_c_31.balance (C1:mlthf_c0_C0) + 0.28742*Ec_biomass_iJO1366_WT_53p95M_arg_DASH_L_c_32.balance (C1:arg_DASH_L_c0_C0 C2:arg_DASH_L_c0_C1 C3:arg_DASH_L_c0_C2 C4:arg_DASH_L_c0_C3 C5:arg_DASH_L_c0_C4 C6:arg_DASH_L_c0_C5) + 0.499149*Ec_biomass_iJO1366_WT_53p95M_ala_DASH_L_c_33.balance (C1:ala_DASH_L_c0_C0 C2:ala_DASH_L_c0_C1 C3:ala_DASH_L_c0_C2) + 0.246506*Ec_biomass_iJO1366_WT_53p95M_thr_DASH_L_c_34.balance (C1:thr_DASH_L_c0_C0 C2:thr_DASH_L_c0_C1 C3:thr_DASH_L_c0_C2 C4:thr_DASH_L_c0_C3) + 0.088988*Ec_biomass_iJO1366_WT_53p95M_cys_DASH_L_c_35.balance (C1:cys_DASH_L_c0_C0 C2:cys_DASH_L_c0_C1 C3:cys_DASH_L_c0_C2) + 0.001787*Ec_biomass_iJO1366_WT_53p95M_nad_c_36.balance (C1:nad_c0_C0 C2:nad_c0_C1 C3:nad_c0_C2 C4:nad_c0_C3 C5:nad_c0_C4 C6:nad_c0_C5 C7:nad_c0_C6 C8:nad_c0_C7 C9:nad_c0_C8 C10:nad_c0_C9 C11:nad_c0_C10 C12:nad_c0_C11 C13:nad_c0_C12 C14:nad_c0_C13 C15:nad_c0_C14 C16:nad_c0_C15 C17:nad_c0_C16 C18:nad_c0_C17 C19:nad_c0_C18 C20:nad_c0_C19 C21:nad_c0_C20) + 0.180021*Ec_biomass_iJO1366_WT_53p95M_phe_DASH_L_c_37.balance (C1:phe_DASH_L_c0_C0 C2:phe_DASH_L_c0_C1 C3:phe_DASH_L_c0_C2 C4:phe_DASH_L_c0_C3 C5:phe_DASH_L_c0_C4 C6:phe_DASH_L_c0_C5 C7:phe_DASH_L_c0_C6 C8:phe_DASH_L_c0_C7 C9:phe_DASH_L_c0_C8) + 0.025612*Ec_biomass_iJO1366_WT_53p95M_dctp_c_38.balance (C1:dctp_c0_C0 C2:dctp_c0_C1 C3:dctp_c0_C2 C4:dctp_c0_C3 C5:dctp_c0_C4 C6:dctp_c0_C5 C7:dctp_c0_C6 C8:dctp_c0_C7 C9:dctp_c0_C8) + 0.149336*Ec_biomass_iJO1366_WT_53p95M_met_DASH_L_c_39.balance (C1:met_DASH_L_c0_C0 C2:met_DASH_L_c0_C1 C3:met_DASH_L_c0_C2 C4:met_DASH_L_c0_C3 C5:met_DASH_L_c0_C4) + 0.012366*Ec_biomass_iJO1366_WT_53p95M_pe160_c_40.balance (C1:pe160_c0_C0 C2:pe160_c0_C1 C3:pe160_c0_C2 C4:pe160_c0_C3 C5:pe160_c0_C4) + 0.209121*Ec_biomass_iJO1366_WT_53p95M_gtp_c_41.balance (C1:gtp_c0_C0 C2:gtp_c0_C1 C3:gtp_c0_C2 C4:gtp_c0_C3 C5:gtp_c0_C4 C6:gtp_c0_C5 C7:gtp_c0_C6 C8:gtp_c0_C7 C9:gtp_c0_C8 C10:gtp_c0_C9) + 0.437778*Ec_biomass_iJO1366_WT_53p95M_leu_DASH_L_c_42.balance (C1:leu_DASH_L_c0_C0 C2:leu_DASH_L_c0_C1 C3:leu_DASH_L_c0_C2 C4:leu_DASH_L_c0_C3 C5:leu_DASH_L_c0_C4 C6:leu_DASH_L_c0_C5) + 0.092056*Ec_biomass_iJO1366_WT_53p95M_his_DASH_L_c_43.balance (C1:his_DASH_L_c0_C0 C2:his_DASH_L_c0_C1 C3:his_DASH_L_c0_C2 C4:his_DASH_L_c0_C3 C5:his_DASH_L_c0_C4 C6:his_DASH_L_c0_C5) + 0.009618*Ec_biomass_iJO1366_WT_53p95M_pe161_c_44.balance (C1:pe161_c0_C0 C2:pe161_c0_C1 C3:pe161_c0_C2 C4:pe161_c0_C3 C5:pe161_c0_C4) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_10fthf_c_45.balance (C1:10fthf_c0_C0) + 0.024805*Ec_biomass_iJO1366_WT_53p95M_datp_c_46.balance (C1:datp_c0_C0 C2:datp_c0_C1 C3:datp_c0_C2 C4:datp_c0_C3 C5:datp_c0_C4 C6:datp_c0_C5 C7:datp_c0_C6 C8:datp_c0_C7 C9:datp_c0_C8 C10:datp_c0_C9) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_5mthf_c_47.balance (C1:5mthf_c0_C0) + 0.000673*Ec_biomass_iJO1366_WT_53p95M_murein4px4px4p_p_48.balance (C1:murein4px4px4p_p0_C0 C2:murein4px4px4p_p0_C1 C3:murein4px4px4p_p0_C2 C4:murein4px4px4p_p0_C3 C5:murein4px4px4p_p0_C4 C6:murein4px4px4p_p0_C5 C7:murein4px4px4p_p0_C6 C8:murein4px4px4p_p0_C7 C9:murein4px4px4p_p0_C8 C10:murein4px4px4p_p0_C9 C11:murein4px4px4p_p0_C10 C12:murein4px4px4p_p0_C11 C13:murein4px4px4p_p0_C12 C14:murein4px4px4p_p0_C13 C15:murein4px4px4p_p0_C14 C16:murein4px4px4p_p0_C15 C17:murein4px4px4p_p0_C16 C18:murein4px4px4p_p0_C17 C19:murein4px4px4p_p0_C18 C20:murein4px4px4p_p0_C19 C21:murein4px4px4p_p0_C20 C22:murein4px4px4p_p0_C21 C23:murein4px4px4p_p0_C22 C24:murein4px4px4p_p0_C23 C25:murein4px4px4p_p0_C24 C26:murein4px4px4p_p0_C25 C27:murein4px4px4p_p0_C26 C28:murein4px4px4p_p0_C27 C29:murein4px4px4p_p0_C28 C30:murein4px4px4p_p0_C29 C31:murein4px4px4p_p0_C30 C32:murein4px4px4p_p0_C31 C33:murein4px4px4p_p0_C32 C34:murein4px4px4p_p0_C33 C35:murein4px4px4p_p0_C34 C36:murein4px4px4p_p0_C35 C37:murein4px4px4p_p0_C36 C38:murein4px4px4p_p0_C37 C39:murein4px4px4p_p0_C38 C40:murein4px4px4p_p0_C39 C41:murein4px4px4p_p0_C40 C42:murein4px4px4p_p0_C41 C43:murein4px4px4p_p0_C42 C44:murein4px4px4p_p0_C43 C45:murein4px4px4p_p0_C44 C46:murein4px4px4p_p0_C45 C47:murein4px4px4p_p0_C46 C48:murein4px4px4p_p0_C47 C49:murein4px4px4p_p0_C48 C50:murein4px4px4p_p0_C49 C51:murein4px4px4p_p0_C50 C52:murein4px4px4p_p0_C51 C53:murein4px4px4p_p0_C52 C54:murein4px4px4p_p0_C53 C55:murein4px4px4p_p0_C54 C56:murein4px4px4p_p0_C55 C57:murein4px4px4p_p0_C56 C58:murein4px4px4p_p0_C57 C59:murein4px4px4p_p0_C58 C60:murein4px4px4p_p0_C59 C61:murein4px4px4p_p0_C60 C62:murein4px4px4p_p0_C61 C63:murein4px4px4p_p0_C62 C64:murein4px4px4p_p0_C63 C65:murein4px4px4p_p0_C64 C66:murein4px4px4p_p0_C65 C67:murein4px4px4p_p0_C66 C68:murein4px4px4p_p0_C67 C69:murein4px4px4p_p0_C68 C70:murein4px4px4p_p0_C69 C71:murein4px4px4p_p0_C70 C72:murein4px4px4p_p0_C71 C73:murein4px4px4p_p0_C72 C74:murein4px4px4p_p0_C73 C75:murein4px4px4p_p0_C74 C76:murein4px4px4p_p0_C75 C77:murein4px4px4p_p0_C76 C78:murein4px4px4p_p0_C77 C79:murein4px4px4p_p0_C78 C80:murein4px4px4p_p0_C79 C81:murein4px4px4p_p0_C80 C82:murein4px4px4p_p0_C81 C83:murein4px4px4p_p0_C82 C84:murein4px4px4p_p0_C83 C85:murein4px4px4p_p0_C84 C86:murein4px4px4p_p0_C85 C87:murein4px4px4p_p0_C86 C88:murein4px4px4p_p0_C87 C89:murein4px4px4p_p0_C88 C90:murein4px4px4p_p0_C89 C91:murein4px4px4p_p0_C90 C92:murein4px4px4p_p0_C91 C93:murein4px4px4p_p0_C92 C94:murein4px4px4p_p0_C93 C95:murein4px4px4p_p0_C94 C96:murein4px4px4p_p0_C95 C97:murein4px4px4p_p0_C96 C98:murein4px4px4p_p0_C97 C99:murein4px4px4p_p0_C98 C100:murein4px4px4p_p0_C99 C101:murein4px4px4p_p0_C100 C102:murein4px4px4p_p0_C101 C103:murein4px4px4p_p0_C102 C104:murein4px4px4p_p0_C103 C105:murein4px4px4p_p0_C104 C106:murein4px4px4p_p0_C105 C107:murein4px4px4p_p0_C106 C108:murein4px4px4p_p0_C107 C109:murein4px4px4p_p0_C108 C110:murein4px4px4p_p0_C109 C111:murein4px4px4p_p0_C110) + 0.024805*Ec_biomass_iJO1366_WT_53p95M_dttp_c_49.balance (C1:dttp_c0_C0 C2:dttp_c0_C1 C3:dttp_c0_C2 C4:dttp_c0_C3 C5:dttp_c0_C4 C6:dttp_c0_C5 C7:dttp_c0_C6 C8:dttp_c0_C7 C9:dttp_c0_C8 C10:dttp_c0_C9) + 0.000223*Ec_biomass_iJO1366_WT_53p95M_ribflv_c_50.balance (C1:ribflv_c0_C0 C2:ribflv_c0_C1 C3:ribflv_c0_C2 C4:ribflv_c0_C3 C5:ribflv_c0_C4 C6:ribflv_c0_C5 C7:ribflv_c0_C6 C8:ribflv_c0_C7 C9:ribflv_c0_C8 C10:ribflv_c0_C9 C11:ribflv_c0_C10 C12:ribflv_c0_C11 C13:ribflv_c0_C12 C14:ribflv_c0_C13 C15:ribflv_c0_C14 C16:ribflv_c0_C15 C17:ribflv_c0_C16) + 0.001345*Ec_biomass_iJO1366_WT_53p95M_murein3p3p_p_51.balance (C1:murein3p3p_p0_C0 C2:murein3p3p_p0_C1 C3:murein3p3p_p0_C2 C4:murein3p3p_p0_C3 C5:murein3p3p_p0_C4 C6:murein3p3p_p0_C5 C7:murein3p3p_p0_C6 C8:murein3p3p_p0_C7 C9:murein3p3p_p0_C8 C10:murein3p3p_p0_C9 C11:murein3p3p_p0_C10 C12:murein3p3p_p0_C11 C13:murein3p3p_p0_C12 C14:murein3p3p_p0_C13 C15:murein3p3p_p0_C14 C16:murein3p3p_p0_C15 C17:murein3p3p_p0_C16 C18:murein3p3p_p0_C17 C19:murein3p3p_p0_C18 C20:murein3p3p_p0_C19 C21:murein3p3p_p0_C20 C22:murein3p3p_p0_C21 C23:murein3p3p_p0_C22 C24:murein3p3p_p0_C23 C25:murein3p3p_p0_C24 C26:murein3p3p_p0_C25 C27:murein3p3p_p0_C26 C28:murein3p3p_p0_C27 C29:murein3p3p_p0_C28 C30:murein3p3p_p0_C29 C31:murein3p3p_p0_C30 C32:murein3p3p_p0_C31 C33:murein3p3p_p0_C32 C34:murein3p3p_p0_C33 C35:murein3p3p_p0_C34 C36:murein3p3p_p0_C35 C37:murein3p3p_p0_C36 C38:murein3p3p_p0_C37 C39:murein3p3p_p0_C38 C40:murein3p3p_p0_C39 C41:murein3p3p_p0_C40 C42:murein3p3p_p0_C41 C43:murein3p3p_p0_C42 C44:murein3p3p_p0_C43 C45:murein3p3p_p0_C44 C46:murein3p3p_p0_C45 C47:murein3p3p_p0_C46 C48:murein3p3p_p0_C47 C49:murein3p3p_p0_C48 C50:murein3p3p_p0_C49 C51:murein3p3p_p0_C50 C52:murein3p3p_p0_C51 C53:murein3p3p_p0_C52 C54:murein3p3p_p0_C53 C55:murein3p3p_p0_C54 C56:murein3p3p_p0_C55 C57:murein3p3p_p0_C56 C58:murein3p3p_p0_C57 C59:murein3p3p_p0_C58 C60:murein3p3p_p0_C59 C61:murein3p3p_p0_C60 C62:murein3p3p_p0_C61 C63:murein3p3p_p0_C62 C64:murein3p3p_p0_C63 C65:murein3p3p_p0_C64 C66:murein3p3p_p0_C65 C67:murein3p3p_p0_C66 C68:murein3p3p_p0_C67) + 0.004892*Ec_biomass_iJO1366_WT_53p95M_pg160_p_52.balance (C1:pg160_p0_C0 C2:pg160_p0_C1 C3:pg160_p0_C2 C4:pg160_p0_C3 C5:pg160_p0_C4 C6:pg160_p0_C5) + 0.129799*Ec_biomass_iJO1366_WT_53p95M_ctp_c_53.balance (C1:ctp_c0_C0 C2:ctp_c0_C1 C3:ctp_c0_C2 C4:ctp_c0_C3 C5:ctp_c0_C4 C6:ctp_c0_C5 C7:ctp_c0_C6 C8:ctp_c0_C7 C9:ctp_c0_C8) + 0.255712*Ec_biomass_iJO1366_WT_53p95M_glu_DASH_L_c_54.balance (C1:glu_DASH_L_c0_C0 C2:glu_DASH_L_c0_C1 C3:glu_DASH_L_c0_C2 C4:glu_DASH_L_c0_C3 C5:glu_DASH_L_c0_C4) + 0.214798*Ec_biomass_iJO1366_WT_53p95M_pro_DASH_L_c_55.balance (C1:pro_DASH_L_c0_C0 C2:pro_DASH_L_c0_C1 C3:pro_DASH_L_c0_C2 C4:pro_DASH_L_c0_C3 C5:pro_DASH_L_c0_C4) + 0.025612*Ec_biomass_iJO1366_WT_53p95M_dgtp_c_56.balance (C1:dgtp_c0_C0 C2:dgtp_c0_C1 C3:dgtp_c0_C2 C4:dgtp_c0_C3 C5:dgtp_c0_C4 C6:dgtp_c0_C5 C7:dgtp_c0_C6 C8:dgtp_c0_C7 C9:dgtp_c0_C8 C10:dgtp_c0_C9) + 0.255712*Ec_biomass_iJO1366_WT_53p95M_gln_DASH_L_c_57.balance (C1:gln_DASH_L_c0_C0 C2:gln_DASH_L_c0_C1 C3:gln_DASH_L_c0_C2 C4:gln_DASH_L_c0_C3 C5:gln_DASH_L_c0_C4) + 0.001961*Ec_biomass_iJO1366_WT_53p95M_pg181_p_58.balance (C1:pg181_p0_C0 C2:pg181_p0_C1 C3:pg181_p0_C2 C4:pg181_p0_C3 C5:pg181_p0_C4 C6:pg181_p0_C5) + 0.024732*Ec_biomass_iJO1366_WT_53p95M_pe161_p_59.balance (C1:pe161_p0_C0 C2:pe161_p0_C1 C3:pe161_p0_C2 C4:pe161_p0_C3 C5:pe161_p0_C4) + 0.00118*Ec_biomass_iJO1366_WT_53p95M_clpn181_p_60.balance (C1:clpn181_p0_C0 C2:clpn181_p0_C1 C3:clpn181_p0_C2 C4:clpn181_p0_C3 C5:clpn181_p0_C4 C6:clpn181_p0_C5 C7:clpn181_p0_C6 C8:clpn181_p0_C7 C9:clpn181_p0_C8) + 0.003805*Ec_biomass_iJO1366_WT_53p95M_pg161_p_61.balance (C1:pg161_p0_C0 C2:pg161_p0_C1 C3:pg161_p0_C2 C4:pg161_p0_C3 C5:pg161_p0_C4 C6:pg161_p0_C5) + 0.000279*Ec_biomass_iJO1366_WT_53p95M_accoa_c_62.balance (C1:accoa_c0_C0 C2:accoa_c0_C1 C3:accoa_c0_C2 C4:accoa_c0_C3 C5:accoa_c0_C4 C6:accoa_c0_C5 C7:accoa_c0_C6 C8:accoa_c0_C7 C9:accoa_c0_C8 C10:accoa_c0_C9 C11:accoa_c0_C10 C12:accoa_c0_C11 C13:accoa_c0_C12 C14:accoa_c0_C13 C15:accoa_c0_C14 C16:accoa_c0_C15 C17:accoa_c0_C16 C18:accoa_c0_C17 C19:accoa_c0_C18 C20:accoa_c0_C19 C21:accoa_c0_C20 C22:accoa_c0_C21 C23:accoa_c0_C22) + 54.119975*Ec_biomass_iJO1366_WT_53p95M_atp_c_63.balance (C1:atp_c0_C0 C2:atp_c0_C1 C3:atp_c0_C2 C4:atp_c0_C3 C5:atp_c0_C4 C6:atp_c0_C5 C7:atp_c0_C6 C8:atp_c0_C7 C9:atp_c0_C8 C10:atp_c0_C9) + 0.133993*Ec_biomass_iJO1366_WT_53p95M_tyr_DASH_L_c_64.balance (C1:tyr_DASH_L_c0_C0 C2:tyr_DASH_L_c0_C1 C3:tyr_DASH_L_c0_C2 C4:tyr_DASH_L_c0_C3 C5:tyr_DASH_L_c0_C4 C6:tyr_DASH_L_c0_C5 C7:tyr_DASH_L_c0_C6 C8:tyr_DASH_L_c0_C7 C9:tyr_DASH_L_c0_C8) + 0.006744*Ec_biomass_iJO1366_WT_53p95M_spmd_c_65.balance (C1:spmd_c0_C0 C2:spmd_c0_C1 C3:spmd_c0_C2 C4:spmd_c0_C3 C5:spmd_c0_C4 C6:spmd_c0_C5 C7:spmd_c0_C6) + 0.002944*Ec_biomass_iJO1366_WT_53p95M_clpn160_p_66.balance (C1:clpn160_p0_C0 C2:clpn160_p0_C1 C3:clpn160_p0_C2 C4:clpn160_p0_C3 C5:clpn160_p0_C4 C6:clpn160_p0_C5 C7:clpn160_p0_C6 C8:clpn160_p0_C7 C9:clpn160_p0_C8) + 0.00229*Ec_biomass_iJO1366_WT_53p95M_clpn161_p_67.balance (C1:clpn161_p0_C0 C2:clpn161_p0_C1 C3:clpn161_p0_C2 C4:clpn161_p0_C3 C5:clpn161_p0_C4 C6:clpn161_p0_C5 C7:clpn161_p0_C6 C8:clpn161_p0_C7 C9:clpn161_p0_C8) + 53.945874*pi_c + 53.95*h_c ';
    def write_isotopomerExperiment_INCA(self, modelReaction_data_I,modelMetabolite_data_I,
                                        measuredFluxes_data_I,experimentalMS_data_I,tracer_I,
                                        parallel_I = 'experiment_id'):
        '''Write matlab script file that describes the fluxomics experiment for INCA1.1'''

        mat_script = 'clear functions\n';

        ##1. Define the model:

        ## debug reaction equations
        #tmp_script = ''
        #for rxn in modelReaction_data_I:
        #    #TODO check on how the reactions are named  
        #    tmp_script = tmp_script + 'r = reaction({...\n';
        #    tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
        #    tmp_script = tmp_script + '});\n';
        #mat_script = mat_script + tmp_script;

        # write out reaction equations
        tmp_script = ''
        tmp_script = tmp_script + 'r = reaction({...\n';
        rxn_ids_INCA = {};
        cnt = 0
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                rxn_ids_INCA[rxn['rxn_id']] = ('R'+str(cnt+1));
                cnt+=1;
                if rxn['rxn_id'] == 'Ec_biomass_iJO1366_WT_53p95M':
                    tmp_script = tmp_script + "'" + self.biomass_INCA + "';...\n"
                    #tmp_script = tmp_script + "'" + self.biomass_INCA_iJS2012 + "';...\n"
                else:
                    tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
                #tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
            #else:
            #    print 'rxn_id ' + rxn['rxn_id'] + ' will be excluded from INCA' 
        tmp_script = tmp_script + '});\n';
        mat_script = mat_script + tmp_script;

        # setup the model
        mat_script = mat_script + 'm = model(r);\n'

        # Take care of symmetrical metabolites if not done so in the reaction equations
        tmp_script = ''
        for met in modelMetabolite_data_I:
            if met['met_symmetry_atompositions']:
                tmp_script = tmp_script + "m.mets{'" + met['met_id'] + "'}.sym = list('rotate180',map('";
                for cnt,atompositions in enumerate(met['met_atompositions']):
                    tmp_script = tmp_script + met['met_elements'][cnt] + str(atompositions+1) + ':' + met['met_symmetry_elements'][cnt] + str(met['met_symmetry_atompositions'][cnt]+1) + ' ';
                tmp_script = tmp_script[:-1];
                tmp_script = tmp_script + "'));\n";
        mat_script = mat_script + tmp_script;

        # Add in the metabolite states (balance), value, and lb/ub)
        tmp_script = ''
        # specify reactions that should be forcible unbalanced
        #NOTE: hard-coded for now until a better workaround can be done
        metabolites_all = [x['met_id'] for x in modelMetabolite_data_I];
        for met in ['co2_e','h2o_e','h_e','na1_e']:
            if met in metabolites_all:
                tmp_script = tmp_script + "m.states{'" + met + ".EX" + "'}.bal = false";
                tmp_script = tmp_script + "'));\n";
        mat_script = mat_script + tmp_script;

        # Add in initial fluxes (values lb/ub) and define the reaction ids
        tmp_script = ''
        tmp_script = tmp_script + 'm.rates.flx.lb = [...\n';
        # lower bounds
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                if measuredFluxes_data_I:
                    for flux in measuredFluxes_data_I:
                        if rxn['rxn_id'] == flux['rxn_id']:
                            tmp_script = tmp_script + str(flux['flux_lb']) + ',...\n'
                            break;
                        else:
                            tmp_script = tmp_script + str(rxn['lower_bound']) + ',...\n'
                            break;
                else: tmp_script = tmp_script + str(rxn['lower_bound']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.flx.ub = [...\n';
        # upper bounds
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                if measuredFluxes_data_I:
                    for flux in measuredFluxes_data_I:
                        if rxn['rxn_id'] == flux['rxn_id']:
                            tmp_script = tmp_script + str(flux['flux_ub']) + ',...\n'
                            break;
                        else:
                            tmp_script = tmp_script + str(rxn['upper_bound']) + ',...\n'
                            break;
                else: tmp_script = tmp_script + str(rxn['upper_bound']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.flx.val = [...\n';
        # intial flux values
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                tmp_script = tmp_script + str(rxn['flux_val']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.on = [...\n';
        # include/exclude a reaction from the simulation
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            if rxn['flux_val']==0.0 and rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0:
                #tmp_script = tmp_script + 'm.rates.on(' + str(rxn_cnt) + ') = 0;\n'
                tmp_script = tmp_script + 'false' + ',...\n'
            else:
                #tmp_script = tmp_script + 'm.rates.on(' + str(rxn_cnt) + ') = 1;\n'
                tmp_script = tmp_script + 'true' + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.id = {...\n';
        # rxn_ids
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                tmp_script = tmp_script + "'" + rxn['rxn_id'] + "',...\n"
        tmp_script = tmp_script + '};\n';
        tmp_script = tmp_script + 'm.rates.id = {...\n';

        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            tmp_script = tmp_script + "'" + rxn['rxn_id'] + "',...\n"
        tmp_script = tmp_script + '};\n';
        mat_script = mat_script + tmp_script;

        ## Check that fluxes are feasible
        #mat_script = mat_script + "m.rates.flx.val = mod2stoich(m)';\n"

        ## Add in the metabolite states (value and lb/ub)
        ##TODO: decide on met_equations structure
        ##NOTE: lb, ub, val = 0 for steady-state
        #mat_script = mat_script + 'm.states.flx.lb = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['lower_bound'] + ',...\n'
        #mat_script = mat_script + '];\n';
        #mat_script = mat_script + 'm.states.flx.ub = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['upper_bound'] + ',...\n'
        #mat_script = mat_script + '];\n';
        #mat_script = mat_script + 'm.states.flx.ub = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['flux'] + ',...\n'
        #mat_script = mat_script + '];\n';

        ##2. Set simulation options

        # Specify simulation parameters (non-stationary only!)
        '''% simulate MS measurements
        nmts = 8;                               % number of total measurements
        samp = 8/60/60;                         % spacing between measurements in hours
        m.options.int_tspan = 0:samp:(samp*nmts);   % time points in hours
        m.options.sim_tunit = 'h';              % hours are unit of time
        m.options.fit_reinit = true;
        m.options.sim_ss = false;
        m.options.sim_sens = true;'''

        tmp_script = ''
        tmp_script = tmp_script + 'm.options.fit_starts = 10;\n' #10 restarts during the estimation procedure
        mat_script = mat_script + tmp_script;
        
        ##3. Define the experiment

        # write out the measured fragment information
        # (actual MS measurements will be written to the script later)

        if parallel_I == 'experiment_id':
            experiments_all = [x['experiment_id'] for x in experimentalMS_data_I];
            experiments = list(set(experiments_all));
            experiments.sort();
        
            fragments_all = [x['fragment_id'] for x in experimentalMS_data_I];
            fragments = list(set(fragments_all));
            fragments.sort();
            mets_all = [x['met_id'] for x in experimentalMS_data_I];
            mets = list(set(mets_all));
            mets.sort();
            times_all = [x['time_point'] for x in experimentalMS_data_I];
            times = list(set(times_all));
            times.sort();

            for experiment_cnt,experiment in enumerate(experiments):
                tmp_script = ''
                tmp_script = tmp_script + 'd = msdata({...\n';
                for fragment in fragments:
                    for ms_data in experimentalMS_data_I:
                        if ms_data['fragment_id'] == fragment and ms_data['experiment_id'] == experiment:
                            tmp_script = tmp_script + "'" + ms_data['fragment_id'] + ': ' + ms_data['met_id'] + ' @ ';
                            for pos_cnt,pos in enumerate(ms_data['met_atompositions']):
                                    tmp_script = tmp_script + ms_data['met_elements'][pos_cnt] + str(pos+1) + ' ';
                            tmp_script = tmp_script[:-1];
                            tmp_script = tmp_script + "';\n"
                            break;
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 'd.mdvs = mdv;\n';
                mat_script = mat_script + tmp_script;

                ## write substrate labeling (i.e. tracer) information
                tmp_script = ''
                tmp_script = tmp_script + 't = tracer({...\n';
                for tracer in tracer_I:
                    if tracer['experiment_id'] == experiment:
                        tmp_script = tmp_script + "'" + tracer['met_name'] + ': ' + tracer['met_id'] + '.EX' + ' @ '
                        for cnt,met_atompositions in enumerate(tracer['met_atompositions']):
                                tmp_script = tmp_script + tracer['met_elements'][cnt]+str(met_atompositions) + ' '
                        tmp_script = tmp_script[:-1]; #remove extra white space
                        tmp_script = tmp_script + "';...\n";
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 't.frac = [';
                for tracer in tracer_I:
                    if tracer['experiment_id'] == experiment:
                        tmp_script = tmp_script + str(tracer['ratio']) + ',';
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + '];\n'; #remove extra ,
                mat_script = mat_script + tmp_script;
                    
                ## write flux measurements
                tmp_script = ''
                tmp_script = tmp_script + "f = data('";
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment:
                        ## Temporary fix until reactions can be properly named
                        #tmp_script = tmp_script + rxn_ids_INCA[flux['rxn_id']] + " ";
                        tmp_script = tmp_script + flux['rxn_id'] + " ";
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + "');\n";
                tmp_script = tmp_script + 'f.val = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment: tmp_script = tmp_script + str(flux['flux_average']) + ',...\n';
                tmp_script = tmp_script + '];\n';
                tmp_script = tmp_script + 'f.std = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment: tmp_script = tmp_script + str(flux['flux_stdev']) + ',...\n';
                tmp_script = tmp_script + '];\n';

                tmp_script = tmp_script + 'x = experiment(t);\n'
                tmp_script = tmp_script + 'x.data_flx = f;\n'
                tmp_script = tmp_script + 'x.data_ms = d;\n'
                tmp_script = tmp_script + ('m.expts(%d) = x;\n' %(experiment_cnt+1));
                tmp_script = tmp_script + ("m.expts(%d).id = {'%s'};\n" %(experiment_cnt+1,experiment));
                mat_script = mat_script + tmp_script;

            # Add in ms data or Write ms data to separate file
            for experiment_cnt,experiment in enumerate(experiments):
                tmp_script = ''
                for i,fragment in enumerate(fragments):
                    for j,time in enumerate(times):
                        # Pad the data file:
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,1,j+1,'NaN'));
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,1,j+1,'NaN'));
                        for ms_data in experimentalMS_data_I:
                            if ms_data['fragment_id']==fragments[i] and \
                                ms_data['time_point']==times[j] and \
                                ms_data['experiment_id'] == experiment:
                                for cnt,intensity in enumerate(ms_data['intensity_normalized_average']):
                                    # each column is a seperate time point
                                    # each row is a seperate mdv
                                    # Assign names and times
                                    name = fragment + '_' + str(cnt) + '_' + str(j) + '_' + str(experiment);
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.id(%d,%d) = {'%s'};\n" %(experiment_cnt+1,i+1,1,j+1,name));
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.time(%d,%d) = %s;\n" %(experiment_cnt+1,i+1,1,j+1,time));
                                    # Assign values
                                    ave = ms_data['intensity_normalized_average'][cnt]
                                    stdev = ms_data['intensity_normalized_stdev'][cnt]
                                    # remove 0.0000 values and replace with NaN
                                    if ave < 1e-6: 
                                        ave = 'NaN';
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,ave));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %f;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,ave));
                                    if stdev < 1e-3:
                                        # check if the ave is NaN
                                        if ave=='NaN': stdev = 'NaN';
                                        elif stdev == 0.0: stdev = 0.05;
                                        else: stdev = 0.001;
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,stdev));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %f;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,stdev));
                mat_script = mat_script + tmp_script;
        elif parallel_I == 'sample_name_abbreviation':
            snas_all = [x['sample_name_abbreviation'] for x in experimentalMS_data_I];
            snas = list(set(snas_all));
            snas.sort();
        
            fragments_all = [x['fragment_id'] for x in experimentalMS_data_I];
            fragments = list(set(fragments_all));
            fragments.sort();
            mets_all = [x['met_id'] for x in experimentalMS_data_I];
            mets = list(set(mets_all));
            mets.sort();
            times_all = [x['time_point'] for x in experimentalMS_data_I];
            times = list(set(times_all));
            times.sort();
        
            for sna_cnt,sna in enumerate(snas):
                tmp_script = ''
                tmp_script = tmp_script + 'd = msdata({...\n';
                for fragment in fragments:
                    for ms_data in snaalMS_data_I:
                        if ms_data['fragment_id'] == fragment and ms_data['sample_name_abbreviation'] == sna:
                            tmp_script = tmp_script + "'" + ms_data['fragment_id'] + ': ' + ms_data['met_id'] + ' @ ';
                            for pos_cnt,pos in enumerate(ms_data['met_atompositions']):
                                    tmp_script = tmp_script + ms_data['met_elements'][pos_cnt] + str(pos+1) + ' ';
                            tmp_script = tmp_script[:-1];
                            tmp_script = tmp_script + "';\n"
                            break;
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 'd.mdvs = mdv;\n';
                mat_script = mat_script + tmp_script;

                ## write substrate labeling (i.e. tracer) information
                tmp_script = ''
                tmp_script = tmp_script + 't = tracer({...\n';
                for tracer in tracer_I:
                    if tracer['sample_name_abbreviation'] == sna:
                        tmp_script = tmp_script + "'" + tracer['met_name'] + ': ' + tracer['met_id'] + '.EX' + ' @ '
                        for cnt,met_atompositions in enumerate(tracer['met_atompositions']):
                                tmp_script = tmp_script + tracer['met_elements'][cnt]+str(met_atompositions) + ' '
                        tmp_script = tmp_script[:-1]; #remove extra white space
                        tmp_script = tmp_script + "';...\n";
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 't.frac = [';
                for tracer in tracer_I:
                    if tracer['sample_name_abbreviation'] == sna:
                        tmp_script = tmp_script + str(tracer['ratio']) + ',';
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + '];\n'; #remove extra ,
                mat_script = mat_script + tmp_script;
                    
                ## write flux measurements
                tmp_script = ''
                tmp_script = tmp_script + "f = data('";
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna:
                        ## Temporary fix until reactions can be properly named
                        #tmp_script = tmp_script + rxn_ids_INCA[flux['rxn_id']] + " ";
                        tmp_script = tmp_script + flux['rxn_id'] + " ";
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + "');\n";
                tmp_script = tmp_script + 'f.val = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna: tmp_script = tmp_script + str(flux['flux_average']) + ',...\n';
                tmp_script = tmp_script + '];\n';
                tmp_script = tmp_script + 'f.std = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna: tmp_script = tmp_script + str(flux['flux_stdev']) + ',...\n';
                tmp_script = tmp_script + '];\n';

                tmp_script = tmp_script + 'x = sna(t);\n'
                tmp_script = tmp_script + 'x.data_flx = f;\n'
                tmp_script = tmp_script + 'x.data_ms = d;\n'
                tmp_script = tmp_script + ('m.expts(%d) = x;\n' %(sna_cnt+1));
                tmp_script = tmp_script + ("m.expts(%d).id = {'%s'};\n" %(sna_cnt+1,sna));
                mat_script = mat_script + tmp_script;

            # Add in ms data or Write ms data to separate file
            for sna_cnt,sna in enumerate(snas):
                tmp_script = ''
                for i,fragment in enumerate(fragments):
                    for j,time in enumerate(times):
                        # Pad the data file:
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(sna_cnt+1,i+1,1,j+1,'NaN'));
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(sna_cnt+1,i+1,1,j+1,'NaN'));
                        for ms_data in snaalMS_data_I:
                            if ms_data['fragment_id']==fragments[i] and \
                                ms_data['time_point']==times[j] and \
                                ms_data['sample_name_abbreviation'] == sna:
                                for cnt,intensity in enumerate(ms_data['intensity_normalized_average']):
                                    # each column is a seperate time point
                                    # each row is a seperate mdv
                                    # Assign names and times
                                    name = fragment + '_' + str(cnt) + '_' + str(j) + '_' + str(sna);
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.id(%d,%d) = {'%s'};\n" %(sna_cnt+1,i+1,1,j+1,name));
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.time(%d,%d) = %s;\n" %(sna_cnt+1,i+1,1,j+1,time));
                                    # Assign values
                                    ave = ms_data['intensity_normalized_average'][cnt]
                                    stdev = ms_data['intensity_normalized_stdev'][cnt]
                                    # remove 0.0000 values and replace with NaN
                                    if ave < 1e-4: 
                                        ave = 'NaN';
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(sna_cnt+1,i+1,cnt+1,j+1,ave));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %f;\n' %(sna_cnt+1,i+1,cnt+1,j+1,ave));
                                    if stdev < 1e-3:
                                        # check if the ave is NaN
                                        if ave=='NaN': stdev = 'NaN';
                                        elif stdev == 0.0: stdev = 0.05;
                                        else: stdev = 0.001;
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(sna_cnt+1,i+1,cnt+1,j+1,stdev));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %f;\n' %(sna_cnt+1,i+1,cnt+1,j+1,stdev));
                mat_script = mat_script + tmp_script;

        return mat_script;
    #Matlab Scripts for INCA
    def writeScript_model_INCA(self, modelReaction_data_I,modelMetabolite_data_I,
                                        measuredFluxes_data_I,experimentalMS_data_I,tracer_I):
        '''Generate the model information for INCA'''

        mat_script = 'clear functions\n';

        ##1. Define the model:

        ## debug reaction equations
        #tmp_script = ''
        #for rxn in modelReaction_data_I:
        #    #TODO check on how the reactions are named  
        #    tmp_script = tmp_script + 'r = reaction({...\n';
        #    tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
        #    tmp_script = tmp_script + '});\n';
        #mat_script = mat_script + tmp_script;

        # write out reaction equations
        tmp_script = ''
        tmp_script = tmp_script + 'r = reaction({...\n';
        rxn_ids_INCA = {};
        cnt = 0
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                rxn_ids_INCA[rxn['rxn_id']] = ('R'+str(cnt+1));
                cnt+=1;
                if rxn['rxn_id'] == 'Ec_biomass_iJO1366_WT_53p95M':
                    tmp_script = tmp_script + "'" + self.biomass_INCA + "';...\n"
                #    #tmp_script = tmp_script + "'" + self.biomass_INCA_iJS2012 + "';...\n"
                else:
                    tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
                #tmp_script = tmp_script + "'" + rxn['rxn_equation'] + "';...\n"
            #else:
            #    print 'rxn_id ' + rxn['rxn_id'] + ' will be excluded from INCA' 
        tmp_script = tmp_script + '});\n';
        mat_script = mat_script + tmp_script;

        # setup the model
        mat_script = mat_script + 'm = model(r);\n'

        # Take care of symmetrical metabolites if not done so in the reaction equations
        tmp_script = ''
        for met in modelMetabolite_data_I:
            if met['met_symmetry_atompositions']:
                tmp_script = tmp_script + "m.mets{'" + met['met_id'] + "'}.sym = list('rotate180',map('";
                for cnt,atompositions in enumerate(met['met_atompositions']):
                    tmp_script = tmp_script + met['met_elements'][cnt] + str(atompositions+1) + ':' + met['met_symmetry_elements'][cnt] + str(met['met_symmetry_atompositions'][cnt]+1) + ' ';
                tmp_script = tmp_script[:-1];
                tmp_script = tmp_script + "'));\n";
        mat_script = mat_script + tmp_script;

        # Add in the metabolite states (balance), value, and lb/ub)
        tmp_script = ''
        # specify metabolites that should be forcible unbalanced
        #NOTE: hard-coded for now until a better workaround can be done
        metabolites_all = [x['met_id'] for x in modelMetabolite_data_I];
        for met in ['co2_e','h2o_e','h_e','na1_e']:
            if met in metabolites_all:
                tmp_script = tmp_script + "m.states{'" + met + ".EX" + "'}.bal = false";
                tmp_script = tmp_script + ";\n";
        for met in modelMetabolite_data_I:
            if '.balance' in met['met_id']:
                tmp_script = tmp_script + "m.states{'" + met['met_id']  + "'}.bal = false";
                tmp_script = tmp_script + ";\n";
        mat_script = mat_script + tmp_script;

        # Add in initial fluxes (values lb/ub) and define the reaction ids
        tmp_script = ''
        tmp_script = tmp_script + 'm.rates.flx.lb = [...\n';
        # lower bounds
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                if measuredFluxes_data_I:
                    for flux in measuredFluxes_data_I:
                        if rxn['rxn_id'] == flux['rxn_id']:
                            tmp_script = tmp_script + str(flux['flux_lb']) + ',...\n'
                            break;
                        else:
                            tmp_script = tmp_script + str(rxn['lower_bound']) + ',...\n'
                            break;
                else: tmp_script = tmp_script + str(rxn['lower_bound']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.flx.ub = [...\n';
        # upper bounds
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                if measuredFluxes_data_I:
                    for flux in measuredFluxes_data_I:
                        if rxn['rxn_id'] == flux['rxn_id']:
                            tmp_script = tmp_script + str(flux['flux_ub']) + ',...\n'
                            break;
                        else:
                            tmp_script = tmp_script + str(rxn['upper_bound']) + ',...\n'
                            break;
                else: tmp_script = tmp_script + str(rxn['upper_bound']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.flx.val = [...\n';
        # intial flux values
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                tmp_script = tmp_script + str(rxn['flux_val']) + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.on = [...\n';
        # include/exclude a reaction from the simulation
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            if rxn['flux_val']==0.0 and rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0:
                #tmp_script = tmp_script + 'm.rates.on(' + str(rxn_cnt) + ') = 0;\n'
                tmp_script = tmp_script + 'false' + ',...\n'
            else:
                #tmp_script = tmp_script + 'm.rates.on(' + str(rxn_cnt) + ') = 1;\n'
                tmp_script = tmp_script + 'true' + ',...\n'
        tmp_script = tmp_script + '];\n';
        tmp_script = tmp_script + 'm.rates.id = {...\n';
        # rxn_ids
        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            #if not(rxn['upper_bound']==0.0 and rxn['lower_bound']==0.0):
                tmp_script = tmp_script + "'" + rxn['rxn_id'] + "',...\n"
        tmp_script = tmp_script + '};\n';
        tmp_script = tmp_script + 'm.rates.id = {...\n';

        for rxn_cnt,rxn in enumerate(modelReaction_data_I):
            tmp_script = tmp_script + "'" + rxn['rxn_id'] + "',...\n"
        tmp_script = tmp_script + '};\n';
        mat_script = mat_script + tmp_script;

        ## Check that fluxes are feasible
        #mat_script = mat_script + "m.rates.flx.val = mod2stoich(m)';\n"

        ## Add in the metabolite states (value and lb/ub)
        ##TODO: decide on met_equations structure
        ##NOTE: lb, ub, val = 0 for steady-state
        #mat_script = mat_script + 'm.states.flx.lb = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['lower_bound'] + ',...\n'
        #mat_script = mat_script + '];\n';
        #mat_script = mat_script + 'm.states.flx.ub = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['upper_bound'] + ',...\n'
        #mat_script = mat_script + '];\n';
        #mat_script = mat_script + 'm.states.flx.ub = [...';
        #for met in modelMetabolite_data_I:
        #    #TODO check on how the metabolites are named
        #    mat_script = mat_script + met['flux'] + ',...\n'
        #mat_script = mat_script + '];\n';

        return mat_script;
    def writeScript_experiment_INCA(self, modelReaction_data_I,modelMetabolite_data_I,
                                        measuredFluxes_data_I,experimentalMS_data_I,tracer_I,
                                        parallel_I = 'experiment_id'):
        '''Generate the experimental information for INCA'''
        
        ##3. Define the experiment

        # write out the measured fragment information
        # (actual MS measurements will be written to the script later)
        mat_script = '';

        if parallel_I == 'experiment_id':
            experiments_all = [x['experiment_id'] for x in experimentalMS_data_I];
            experiments = list(set(experiments_all));
            experiments.sort();
        
            fragments_all = [x['fragment_id'] for x in experimentalMS_data_I];
            fragments = list(set(fragments_all));
            fragments.sort();
            mets_all = [x['met_id'] for x in experimentalMS_data_I];
            mets = list(set(mets_all));
            mets.sort();
            times_all = [x['time_point'] for x in experimentalMS_data_I];
            times = list(set(times_all));
            times.sort();

            for experiment_cnt,experiment in enumerate(experiments):
                tmp_script = ''
                tmp_script = tmp_script + 'd = msdata({...\n';
                for fragment in fragments:
                    for ms_data in experimentalMS_data_I:
                        if ms_data['fragment_id'] == fragment and ms_data['experiment_id'] == experiment:
                            tmp_script = tmp_script + "'" + ms_data['fragment_id'] + ': ' + ms_data['met_id'] + ' @ ';
                            for pos_cnt,pos in enumerate(ms_data['met_atompositions']):
                                    tmp_script = tmp_script + ms_data['met_elements'][pos_cnt] + str(pos+1) + ' ';
                            tmp_script = tmp_script[:-1];
                            tmp_script = tmp_script + "';\n"
                            break;
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 'd.mdvs = mdv;\n';
                mat_script = mat_script + tmp_script;

                ## write substrate labeling (i.e. tracer) information
                tmp_script = ''
                tmp_script = tmp_script + 't = tracer({...\n';
                for tracer in tracer_I:
                    if tracer['experiment_id'] == experiment:
                        tmp_script = tmp_script + "'" + tracer['met_name'] + ': ' + tracer['met_id'] + '.EX' + ' @ '
                        for cnt,met_atompositions in enumerate(tracer['met_atompositions']):
                                tmp_script = tmp_script + tracer['met_elements'][cnt]+str(met_atompositions) + ' '
                        tmp_script = tmp_script[:-1]; #remove extra white space
                        tmp_script = tmp_script + "';...\n";
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 't.frac = [';
                for tracer in tracer_I:
                    if tracer['experiment_id'] == experiment:
                        tmp_script = tmp_script + str(tracer['ratio']) + ',';
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + '];\n'; #remove extra ,
                mat_script = mat_script + tmp_script;
                    
                ## write flux measurements
                tmp_script = ''
                tmp_script = tmp_script + "f = data('";
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment:
                        ## Temporary fix until reactions can be properly named
                        #tmp_script = tmp_script + rxn_ids_INCA[flux['rxn_id']] + " ";
                        tmp_script = tmp_script + flux['rxn_id'] + " ";
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + "');\n";
                tmp_script = tmp_script + 'f.val = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment: tmp_script = tmp_script + str(flux['flux_average']) + ',...\n';
                tmp_script = tmp_script + '];\n';
                tmp_script = tmp_script + 'f.std = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['experiment_id'] == experiment: tmp_script = tmp_script + str(flux['flux_stdev']) + ',...\n';
                tmp_script = tmp_script + '];\n';

                tmp_script = tmp_script + 'x = experiment(t);\n'
                tmp_script = tmp_script + 'x.data_flx = f;\n'
                tmp_script = tmp_script + 'x.data_ms = d;\n'
                tmp_script = tmp_script + ('m.expts(%d) = x;\n' %(experiment_cnt+1));
                tmp_script = tmp_script + ("m.expts(%d).id = {'%s'};\n" %(experiment_cnt+1,experiment));
                mat_script = mat_script + tmp_script;

            # Add in ms data or Write ms data to separate file
            for experiment_cnt,experiment in enumerate(experiments):
                tmp_script = ''
                for i,fragment in enumerate(fragments):
                    for j,time in enumerate(times):
                        # Pad the data file:
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,1,j+1,'NaN'));
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,1,j+1,'NaN'));
                        for ms_data in experimentalMS_data_I:
                            if ms_data['fragment_id']==fragments[i] and \
                                ms_data['time_point']==times[j] and \
                                ms_data['experiment_id'] == experiment:
                                if not ms_data['intensity_normalized_average']: continue; #measurements will need to be added/simulated later
                                for cnt,intensity in enumerate(ms_data['intensity_normalized_average']):
                                    # each column is a seperate time point
                                    # each row is a seperate mdv
                                    # Assign names and times
                                    name = fragment + '_' + str(cnt) + '_' + str(j) + '_' + str(experiment);
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.id(%d,%d) = {'%s'};\n" %(experiment_cnt+1,i+1,1,j+1,name));
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.time(%d,%d) = %s;\n" %(experiment_cnt+1,i+1,1,j+1,time));
                                    # Assign values
                                    ave = ms_data['intensity_normalized_average'][cnt]
                                    stdev = ms_data['intensity_normalized_stdev'][cnt]
                                    # remove 0.0000 values and replace with NaN
                                    if ave < 1e-6: 
                                        ave = 'NaN';
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,ave));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %f;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,ave));
                                    if stdev < 1e-3:
                                        # check if the ave is NaN
                                        if ave=='NaN': stdev = 'NaN';
                                        elif stdev == 0.0: stdev = 0.05;
                                        else: stdev = 0.001;
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,stdev));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %f;\n' %(experiment_cnt+1,i+1,cnt+1,j+1,stdev));
                mat_script = mat_script + tmp_script;
        elif parallel_I == 'sample_name_abbreviation':
            snas_all = [x['sample_name_abbreviation'] for x in experimentalMS_data_I];
            snas = list(set(snas_all));
            snas.sort();
        
            fragments_all = [x['fragment_id'] for x in experimentalMS_data_I];
            fragments = list(set(fragments_all));
            fragments.sort();
            mets_all = [x['met_id'] for x in experimentalMS_data_I];
            mets = list(set(mets_all));
            mets.sort();
            times_all = [x['time_point'] for x in experimentalMS_data_I];
            times = list(set(times_all));
            times.sort();
        
            for sna_cnt,sna in enumerate(snas):
                tmp_script = ''
                tmp_script = tmp_script + 'd = msdata({...\n';
                for fragment in fragments:
                    for ms_data in experimentalMS_data_I:
                        if ms_data['fragment_id'] == fragment and ms_data['sample_name_abbreviation'] == sna:
                            tmp_script = tmp_script + "'" + ms_data['fragment_id'] + ': ' + ms_data['met_id'] + ' @ ';
                            for pos_cnt,pos in enumerate(ms_data['met_atompositions']):
                                    tmp_script = tmp_script + ms_data['met_elements'][pos_cnt] + str(pos+1) + ' ';
                            tmp_script = tmp_script[:-1];
                            tmp_script = tmp_script + "';\n"
                            break;
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 'd.mdvs = mdv;\n';
                mat_script = mat_script + tmp_script;

                ## write substrate labeling (i.e. tracer) information
                tmp_script = ''
                tmp_script = tmp_script + 't = tracer({...\n';
                for tracer in tracer_I:
                    if tracer['sample_name_abbreviation'] == sna:
                        tmp_script = tmp_script + "'" + tracer['met_name'] + ': ' + tracer['met_id'] + '.EX' + ' @ '
                        for cnt,met_atompositions in enumerate(tracer['met_atompositions']):
                                tmp_script = tmp_script + tracer['met_elements'][cnt]+str(met_atompositions) + ' '
                        tmp_script = tmp_script[:-1]; #remove extra white space
                        tmp_script = tmp_script + "';...\n";
                tmp_script = tmp_script + '});\n';
                tmp_script = tmp_script + 't.frac = [';
                for tracer in tracer_I:
                    if tracer['sample_name_abbreviation'] == sna:
                        tmp_script = tmp_script + str(tracer['ratio']) + ',';
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + '];\n'; #remove extra ,
                mat_script = mat_script + tmp_script;
                    
                ## write flux measurements
                tmp_script = ''
                tmp_script = tmp_script + "f = data('";
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna:
                        ## Temporary fix until reactions can be properly named
                        #tmp_script = tmp_script + rxn_ids_INCA[flux['rxn_id']] + " ";
                        tmp_script = tmp_script + flux['rxn_id'] + " ";
                tmp_script = tmp_script[:-1]; #remove extra ,
                tmp_script = tmp_script + "');\n";
                tmp_script = tmp_script + 'f.val = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna: tmp_script = tmp_script + str(flux['flux_average']) + ',...\n';
                tmp_script = tmp_script + '];\n';
                tmp_script = tmp_script + 'f.std = [...\n';
                for flux in measuredFluxes_data_I:
                    if flux['sample_name_abbreviation'] == sna: tmp_script = tmp_script + str(flux['flux_stdev']) + ',...\n';
                tmp_script = tmp_script + '];\n';
                
                tmp_script = tmp_script + 'x = experiment(t);\n'
                tmp_script = tmp_script + 'x.data_flx = f;\n'
                tmp_script = tmp_script + 'x.data_ms = d;\n'
                tmp_script = tmp_script + ('m.expts(%d) = x;\n' %(sna_cnt+1));
                tmp_script = tmp_script + ("m.expts(%d).id = {'%s'};\n" %(sna_cnt+1,sna));
                mat_script = mat_script + tmp_script;

            # Add in ms data or Write ms data to separate file
            for sna_cnt,sna in enumerate(snas):
                tmp_script = ''
                for i,fragment in enumerate(fragments):
                    for j,time in enumerate(times):
                        # Pad the data file:
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(sna_cnt+1,i+1,1,j+1,'NaN'));
                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(sna_cnt+1,i+1,1,j+1,'NaN'));
                        for ms_data in experimentalMS_data_I:
                            if ms_data['fragment_id']==fragments[i] and \
                                ms_data['time_point']==times[j] and \
                                ms_data['sample_name_abbreviation'] == sna:
                                if not ms_data['intensity_normalized_average']: continue; #measurements will need to be added/simulated later
                                for cnt,intensity in enumerate(ms_data['intensity_normalized_average']):
                                    # each column is a seperate time point
                                    # each row is a seperate mdv
                                    # Assign names and times
                                    name = fragment + '_' + str(cnt) + '_' + str(j) + '_' + str(sna);
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.id(%d,%d) = {'%s'};\n" %(sna_cnt+1,i+1,1,j+1,name));
                                    tmp_script = tmp_script + ("m.expts(%d).data_ms(%d).mdvs.time(%d,%d) = %s;\n" %(sna_cnt+1,i+1,1,j+1,time));
                                    # Assign values
                                    ave = ms_data['intensity_normalized_average'][cnt]
                                    stdev = ms_data['intensity_normalized_stdev'][cnt]
                                    # remove 0.0000 values and replace with NaN
                                    if ave < 1e-6: 
                                        ave = 'NaN';
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %s;\n' %(sna_cnt+1,i+1,cnt+1,j+1,ave));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.val(%d,%d) = %f;\n' %(sna_cnt+1,i+1,cnt+1,j+1,ave));
                                    if stdev < 1e-3:
                                        # check if the ave is NaN
                                        if ave=='NaN': stdev = 'NaN';
                                        elif stdev == 0.0: stdev = 0.05;
                                        else: stdev = 0.001;
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %s;\n' %(sna_cnt+1,i+1,cnt+1,j+1,stdev));
                                    else:
                                        tmp_script = tmp_script + ('m.expts(%d).data_ms(%d).mdvs.std(%d,%d) = %f;\n' %(sna_cnt+1,i+1,cnt+1,j+1,stdev));
                mat_script = mat_script + tmp_script;
        return mat_script
    def writeScript_experimentFromMSData_INCA(self, modelReaction_data_I,modelMetabolite_data_I,
                                        measuredFluxes_data_I,experimentalMS_data_I,tracer_I,
                                        experiment_index_I=1,experiment_name_I=None):
        '''Generate the experimental MS data for INCA'''
        return
    def writeScript_experimentFromSimulation_INCA(self, modelReaction_data_I,modelMetabolite_data_I,
                                        measuredFluxes_data_I,experimentalMS_data_I,tracer_I,
                                        experiment_index_I=1,experiment_name_I=None):
        '''Generate simulated MS data for INCA'''

        mat_script = ''

        # check that the fluxes are feasible
        mat_script += "m.rates.flx.val = mod2stoich(m)';\n";

        # simulate measurements
        mat_script += "s = simulate(m);\n"
        mat_script += "m = sim2mod(m,s);\n" # copy simulated measurements into model

        # GC/MS standard error ranges linearly
        mat_script += "x0 = 0.005; e0 = 0.003; x1 = 0.25; e1 = 0.01;\n"
        mat_script += ("for i = 1:length(m.expts(%d).data_ms)\n" %(experiment_index_I))
        mat_script += ("\tm.expts(%d).data_ms(i).mdvs.std = max(min((e1-e0)/(x1-x0)*(m.expts(%d)data_ms(i).mdvs.val-x1)+e1,e1),e0);\n" %(experiment_index_I,experiment_index_I))
        mat_script += "end\n"

        # Introduce randomly distributed error into measurements
        mat_script += ("for i = 1:length(m.expts(%d).data_ms)\n" %(experiment_index_I))
        mat_script += ("m.expts(%d).data_ms(i).mdvs.val = normrnd(m.expts(%d).data_ms(i).mdvs.val,m.expts(%d).data_ms(i).mdvs.std" %(experiment_index_I,experiment_index_I,experiment_index_I))
        mat_script += "end\n"

        return mat_script
    def writeScript_simulationOptions_Inca(self,stationary_I=True):
        '''Generate parameters for isotopomer simulation for INCA1.1'''

        # Specify simulation parameters (non-stationary only!)
        '''% simulate MS measurements
        nmts = 8;                               % number of total measurements
        samp = 8/60/60;                         % spacing between measurements in hours
        m.options.int_tspan = 0:samp:(samp*nmts);   % time points in hours
        m.options.sim_tunit = 'h';              % hours are unit of time
        m.options.fit_reinit = true;
        m.options.sim_ss = false;
        m.options.sim_sens = true;'''

        mat_script = ''
        mat_script += 'm.options.fit_starts = 10;\n' #10 restarts during the estimation procedure

        return mat_script
    def writeScript_parameterEstimation_Inca(self):
        '''Run parameter estimations INCA1.1'''

        mat_script = ''
        mat_script += "f=estimate(m,10);\n" #10 restarts

        return mat_script
    def writeScript_parameterContinuation_Inca(self):
        '''Run parameter continuation INCA1.1'''

        mat_script = ''
        mat_script += "f=continuate(f,m);\n"

        return mat_script
    def make_isotopomerRxnEquations_INCA(self,reactants_ids_I = [],
                                        products_ids_I = [],
                                        reactants_stoichiometry_I = [],
                                        products_stoichiometry_I = [],
                                        reversibility_I = True,
                                        reactants_stoichiometry_tracked_I = [],
                                        products_stoichiometry_tracked_I = [],
                                        reactants_ids_tracked_I = [],
                                        products_ids_tracked_I = [],
                                        reactants_elements_tracked_I = [[]],
                                        products_elements_tracked_I = [[]],
                                        reactants_positions_tracked_I = [[]],
                                        products_positions_tracked_I = [[]],
                                        reactants_mapping_I = [],
                                        products_mapping_I = []):
        '''Generate string represention of reactions equations for INCA1.1'''

        #e.g. A (AabB) + H (c) -> B (AbBc) + H (a)
        #e.g. A (C1:A C2:B H1R:a H1S:b) + H (H1:c) -> B (C1:A C2:B H1:a H2:c) + H (H1:a)

        rxn_equations_INCA = '';

        #add balance dummy metabolite for an exchange reaction
        if len(reactants_ids_I)==0 and len(products_ids_I)==1:
            reactants_ids_I=[products_ids_I[0] + '.EX'];
            reactants_stoichiometry_I=[products_stoichiometry_I[0]]
            if products_ids_tracked_I and products_ids_tracked_I[0]:
                reactants_stoichiometry_tracked_I=[products_stoichiometry_tracked_I[0]]
                reactants_ids_tracked_I=[products_ids_tracked_I[0] + '.EX']
                reactants_elements_tracked_I=[products_elements_tracked_I[0]]
                reactants_positions_tracked_I=[products_positions_tracked_I[0]]
                reactants_mapping_I=[products_mapping_I[0]]
        elif len(products_ids_I)==0 and len(reactants_ids_I)==1:
            products_ids_I=[reactants_ids_I[0] + '.EX'];
            products_stoichiometry_I=[reactants_stoichiometry_I[0]]
            if reactants_ids_tracked_I and reactants_ids_tracked_I[0]:
                products_stoichiometry_tracked_I=[reactants_stoichiometry_tracked_I[0]]
                products_ids_tracked_I=[reactants_ids_tracked_I[0] + '.EX']
                products_elements_tracked_I=[reactants_elements_tracked_I[0]]
                products_positions_tracked_I=[reactants_positions_tracked_I[0]]
                products_mapping_I=[reactants_mapping_I[0]]

        # pseudo metabolites
        pseudo_mets = [];

        #build the string for the reactants
        for reactant_cnt,reactant in enumerate(reactants_ids_I):
            reactants_stoichiometry = abs(reactants_stoichiometry_I[reactant_cnt]);
            if reactant in reactants_ids_tracked_I:
                for reactant_tracked_cnt, reactant_tracked in enumerate(reactants_ids_tracked_I):
                    if reactant_tracked == reactant and \
                        reactants_stoichiometry == abs(reactants_stoichiometry_tracked_I[reactant_tracked_cnt]):
                        # if the tracked reactant matches and the stoichiometry aggrees, 
                        # combine the information for the tracked_reactant and the reactant
                        #rxn_equations_INCA += ('%f' % reactants_stoichiometry) + '*' + reactant + ' ';
                        rxn_equations_INCA += str(reactants_stoichiometry) + '*' + reactant + ' ';
                        rxn_equations_INCA += '(';
                        if reactants_mapping_I[reactant_tracked_cnt]:
                            reactants_mapping = reactants_mapping_I[reactant_tracked_cnt];
                            if '[' in reactants_mapping_I[reactant_tracked_cnt]:
                                reactants_mapping = reactants_mapping_I[reactant_tracked_cnt].split('][');
                                reactants_mapping = [m.replace('[','') for m in reactants_mapping];
                                reactants_mapping = [m.replace(']','') for m in reactants_mapping];
                            for mapping_cnt, mapping in enumerate(reactants_mapping):
                                rxn_equations_INCA += reactants_elements_tracked_I[reactant_tracked_cnt][mapping_cnt] + \
                                    str(reactants_positions_tracked_I[reactant_tracked_cnt][mapping_cnt]+1) + \
                                    ':' + mapping + ' ';
                            rxn_equations_INCA = rxn_equations_INCA[:-1];
                            rxn_equations_INCA += ') ';
                            rxn_equations_INCA += '+ ';
                            reactants_stoichiometry -= abs(reactants_stoichiometry_tracked_I[reactant_tracked_cnt]); # subtract out the stoichiometry of the tracked reactant from the reactant
                    elif reactant_tracked == reactant and \
                        reactants_stoichiometry > abs(reactants_stoichiometry_tracked_I[reactant_tracked_cnt]):
                        # if the tracked reactant matches and the stoichiometry of the reactant is greater than the tracked_reactant, 
                        # combine the information for the tracked_reactant and the reactant after substracking out the stoichiometry of the tracked_reactant
                        rxn_equations_INCA += str(abs(reactants_stoichiometry_tracked_I[reactant_tracked_cnt])) + '*' + reactant + ' ';
                        rxn_equations_INCA += '(';
                        if reactants_mapping_I[reactant_tracked_cnt]:
                            reactants_mapping = reactants_mapping_I[reactant_tracked_cnt];
                            if '[' in reactants_mapping_I[reactant_tracked_cnt]:
                                reactants_mapping = reactants_mapping_I[reactant_tracked_cnt].split('][');
                                reactants_mapping = [m.replace('[','') for m in reactants_mapping];
                                reactants_mapping = [m.replace(']','') for m in reactants_mapping];
                            for mapping_cnt, mapping in enumerate(reactants_mapping):
                                rxn_equations_INCA += reactants_elements_tracked_I[reactant_tracked_cnt][mapping_cnt] + \
                                    str(reactants_positions_tracked_I[reactant_tracked_cnt][mapping_cnt]+1) + \
                                    ':' + mapping + ' ';
                            rxn_equations_INCA = rxn_equations_INCA[:-1];
                            rxn_equations_INCA += ') ';
                            rxn_equations_INCA += '+ ';
                            reactants_stoichiometry -= abs(reactants_stoichiometry_tracked_I[reactant_tracked_cnt]) # subtract out the stoichiometry of the tracked reactant from the reactant and continue
                    elif not reactant_tracked in reactants_ids_I:
                        print('unaccounted for reactant_tracked: ' + reactant_tracked);
                        if reactants_stoichiometry_tracked_I[reactant_tracked_cnt]==-1e-13 and not reactant_tracked in pseudo_mets:
                            #add in the pseudo-metabolite used to complete the atom mapping
                            rxn_equations_INCA += '0.0000000000001' + '*' + reactant_tracked + ' ';
                            rxn_equations_INCA += '(';
                            if reactants_mapping_I[reactant_tracked_cnt]:
                                reactants_mapping = reactants_mapping_I[reactant_tracked_cnt];
                                if '[' in reactants_mapping_I[reactant_tracked_cnt]:
                                    reactants_mapping = reactants_mapping_I[reactant_tracked_cnt].split('][');
                                    reactants_mapping = [m.replace('[','') for m in reactants_mapping];
                                    reactants_mapping = [m.replace(']','') for m in reactants_mapping];
                                for mapping_cnt, mapping in enumerate(reactants_mapping):
                                    rxn_equations_INCA += reactants_elements_tracked_I[reactant_tracked_cnt][mapping_cnt] + \
                                        str(reactants_positions_tracked_I[reactant_tracked_cnt][mapping_cnt]+1) + \
                                        ':' + mapping + ' ';
                                rxn_equations_INCA = rxn_equations_INCA[:-1];
                                rxn_equations_INCA += ') ';
                                rxn_equations_INCA += '+ ';
                                pseudo_mets.append(reactant_tracked); # only 1 unique pseudo_met per reaction
                if reactants_stoichiometry>0.0: # check if there is a remainder of the reactant stoichiometry that has not yet been accounted for
                                                # after iterating through all reactants
                                                # if so, the molecule is not tracked
                    rxn_equations_INCA += str(reactants_stoichiometry) + '*' + reactant + ' ';
                    rxn_equations_INCA += '+ ';
            else:
                rxn_equations_INCA += str(reactants_stoichiometry) + '*' + reactant + ' ';
                rxn_equations_INCA += '+ ';
        rxn_equations_INCA = rxn_equations_INCA[:-2];
        if reversibility_I:
            rxn_equations_INCA += '<->  ';
        else:
            rxn_equations_INCA += '->  ';
        #build the string for the products
        for product_cnt,product in enumerate(products_ids_I):
            if product_cnt == 0: rxn_equations_INCA = rxn_equations_INCA[:-1];
            #rxn_equations_INCA += str(abs(reactants_stoichiometry_I[product_cnt])) + '*' + product + ' ';
            products_stoichiometry = abs(products_stoichiometry_I[product_cnt]);
            if product in products_ids_tracked_I:
                for product_tracked_cnt, product_tracked in enumerate(products_ids_tracked_I):
                    if product_tracked == product and \
                        products_stoichiometry == abs(products_stoichiometry_tracked_I[product_tracked_cnt]):
                        # if the tracked product matches and the stoichiometry aggrees, 
                        # combine the information for the tracked_product and the product
                        rxn_equations_INCA += str(products_stoichiometry) + '*' + product + ' ';
                        rxn_equations_INCA += '(';
                        if products_mapping_I[product_tracked_cnt]:
                            products_mapping = products_mapping_I[product_tracked_cnt];
                            if '[' in products_mapping_I[product_tracked_cnt]:
                                products_mapping = products_mapping_I[product_tracked_cnt].split('][');
                                products_mapping = [m.replace('[','') for m in products_mapping];
                                products_mapping = [m.replace(']','') for m in products_mapping];
                            for mapping_cnt, mapping in enumerate(products_mapping):
                                rxn_equations_INCA += products_elements_tracked_I[product_tracked_cnt][mapping_cnt] + \
                                    str(products_positions_tracked_I[product_tracked_cnt][mapping_cnt]+1) + \
                                    ':' + mapping + ' ';
                            rxn_equations_INCA = rxn_equations_INCA[:-1];
                            rxn_equations_INCA += ') ';
                            rxn_equations_INCA += '+ ';
                            products_stoichiometry -= abs(products_stoichiometry_tracked_I[product_tracked_cnt]);
                    elif product_tracked == product and \
                        products_stoichiometry > abs(products_stoichiometry_tracked_I[product_tracked_cnt]):
                        # if the tracked product matches and the stoichiometry of the product is greater than the tracked_product, 
                        # combine the information for the tracked_product and the product after substracking out the stoichiometry of the tracked_product
                        rxn_equations_INCA += str(products_stoichiometry_tracked_I[product_tracked_cnt]) + '*' + product + ' ';
                        rxn_equations_INCA += '(';
                        if products_mapping_I[product_tracked_cnt]:
                            products_mapping = products_mapping_I[product_tracked_cnt];
                            if '[' in products_mapping_I[product_tracked_cnt]:
                                products_mapping = products_mapping_I[product_tracked_cnt].split('][');
                                products_mapping = [m.replace('[','') for m in products_mapping];
                                products_mapping = [m.replace(']','') for m in products_mapping];
                            for mapping_cnt, mapping in enumerate(products_mapping):
                                rxn_equations_INCA += products_elements_tracked_I[product_tracked_cnt][mapping_cnt] + \
                                    str(products_positions_tracked_I[product_tracked_cnt][mapping_cnt]+1) + \
                                    ':' + mapping + ' ';
                            rxn_equations_INCA = rxn_equations_INCA[:-1];
                            rxn_equations_INCA += ') ';
                            rxn_equations_INCA += '+ ';
                            products_stoichiometry -= abs(products_stoichiometry_tracked_I[product_tracked_cnt])
                    elif not product_tracked in products_ids_I:
                        print('unaccounted for product_tracked: ' + product_tracked);
                        if '.balance' in product_tracked and not product_tracked in pseudo_mets:
                            #add in the pseudo-metabolite used to complete the atom mapping
                            rxn_equations_INCA += str(products_stoichiometry_tracked_I[product_tracked_cnt]) + '*' + product_tracked + ' ';
                            rxn_equations_INCA += '(';
                            if products_mapping_I[product_tracked_cnt]:
                                products_mapping = products_mapping_I[product_tracked_cnt];
                                if '[' in products_mapping_I[product_tracked_cnt]:
                                    products_mapping = products_mapping_I[product_tracked_cnt].split('][');
                                    products_mapping = [m.replace('[','') for m in products_mapping];
                                    products_mapping = [m.replace(']','') for m in products_mapping];
                                for mapping_cnt, mapping in enumerate(products_mapping):
                                    rxn_equations_INCA += products_elements_tracked_I[product_tracked_cnt][mapping_cnt] + \
                                        str(products_positions_tracked_I[product_tracked_cnt][mapping_cnt]+1) + \
                                        ':' + mapping + ' ';
                                rxn_equations_INCA = rxn_equations_INCA[:-1];
                                rxn_equations_INCA += ') ';
                                rxn_equations_INCA += '+ ';
                                pseudo_mets.append(product_tracked); # only 1 unique pseudo_met per reaction
                if products_stoichiometry>0.0:
                    rxn_equations_INCA += str(products_stoichiometry) + '*' + product + ' ';
                    rxn_equations_INCA += '+ ';
            else:
                rxn_equations_INCA += str(products_stoichiometry) + '*' + product + ' ';
                rxn_equations_INCA += '+ ';
        rxn_equations_INCA = rxn_equations_INCA[:-2];

        #add in unbalanced metabolites to the products

        return rxn_equations_INCA;