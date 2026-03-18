function [linStabOut] = linearizeFC(FC_point,desPerturb,aeroData,massData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This function will compute stability derivatives that
% represent the linearized system about a specified point. When used within
% a trim routine, this will linearize the system about the trim point. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Input Data for Sim Run 
% Store number of sim runs, *2 because +/- for sensitivities 
numRuns = 2*numel(fieldnames(desPerturb));

% Define data structure for all input parameters needed for sim run 
% Note: Requires input data to have sub parameters labeled as below: 

ALT     = ones(numRuns,1)*FC_point.alt_ft;
MACH    = ones(numRuns,1)*FC_point.Mach;
ALPHA   = ones(numRuns,1)*FC_point.alpha_deg;
BETA    = ones(numRuns,1)*FC_point.beta_deg;
PB      = ones(numRuns,1)*FC_point.pb_rps;
QB      = ones(numRuns,1)*FC_point.qb_rps;
RB      = ones(numRuns,1)*FC_point.rb_rps;
DA      = ones(numRuns,1)*FC_point.da_deg;
DE      = ones(numRuns,1)*FC_point.de_deg;
DR      = ones(numRuns,1)*FC_point.dr_deg;
DLEF    = ones(numRuns,1)*FC_point.dlef_deg;
XCG     = ones(numRuns,1)*FC_point.xcg_ft;
YCG     = ones(numRuns,1)*FC_point.ycg_ft;
ZCG     = ones(numRuns,1)*FC_point.zcg_ft;
THR     = ones(numRuns,1)*FC_point.thr;

%% Define sensitivities by replacing specific runs with sensitivity parameter
ALPHA(1)    = ALPHA(1) + desPerturb.dAlpha_deg;
ALPHA(2)    = ALPHA(2) - desPerturb.dAlpha_deg;
BETA(3)     = BETA(3) + desPerturb.dBeta_deg;
BETA(4)     = BETA(4) - desPerturb.dBeta_deg;
MACH(5)     = MACH(5) + (desPerturb.dVair/FC_point.SoS_fps);
MACH(6)     = MACH(6) - (desPerturb.dVair/FC_point.SoS_fps);
THR(7)      = THR(7) + desPerturb.dthr;
THR(8)      = THR(8) - desPerturb.dthr;
DA(9)       = DA(9) + desPerturb.da_deg;
DA(10)      = DA(10) - desPerturb.da_deg;
DE(11)      = DE(11) + desPerturb.de_deg;
DE(12)      = DE(12) - desPerturb.de_deg;
DR(13)      = DR(13) + desPerturb.dr_deg;
DR(14)      = DR(14) - desPerturb.dr_deg;
PB(15)      = PB(15) + desPerturb.dp_rps;
PB(16)      = PB(16) - desPerturb.dp_rps;
QB(17)      = QB(17) + desPerturb.dq_rps;
QB(18)      = QB(18) - desPerturb.dq_rps;
RB(19)      = RB(19) + desPerturb.dr_rps;
RB(20)      = RB(20) - desPerturb.dr_rps;

% Create time vector to run all array points
timeVec = 0:1:numRuns-1;

%% Define input vector for Sim 
sncInVec = [timeVec',ALPHA,BETA,MACH,THR,DA,DE,DR,PB,QB,RB,ALT,DLEF,XCG,YCG,ZCG];

%% Define SnC Model for Calculation
SnC_mdl = 'SnC.slx';
load_system(SnC_mdl);
set_param(SnC_mdl(1:end-4),'StopTime',num2str(timeVec(end)));
% Run sim
simOut = sim(SnC_mdl,'','',sncInVec);
% Force close model
bdclose(SnC_mdl);
% Process SnC Outputs
sncOut = processSimData(simOut);

%% Compute Non-Dimensional Stability Derivative
d2r = pi/180;
% Alpha (1/rad)
linStabOut.CD_alpha = (sncOut.CDS_total(1)-sncOut.CDS_total(2))/(2*desPerturb.dAlpha_deg*d2r);
linStabOut.CL_alpha = (sncOut.CLS_total(1)-sncOut.CLS_total(2))/(2*desPerturb.dAlpha_deg*d2r);
linStabOut.Cm_alpha = (sncOut.CMS_total(1)-sncOut.CMS_total(2))/(2*desPerturb.dAlpha_deg*d2r);
% Beta (1/rad)
linStabOut.CY_beta  = (sncOut.CYS_total(3)-sncOut.CYS_total(4))/(2*desPerturb.dBeta_deg*d2r);
linStabOut.Cn_beta  = (sncOut.CNS_total(3)-sncOut.CNS_total(4))/(2*desPerturb.dBeta_deg*d2r);
linStabOut.Cr_beta  = (sncOut.CRS_total(3)-sncOut.CRS_total(4))/(2*desPerturb.dBeta_deg*d2r);
% Speed (unitless)
linStabOut.CD_U     = FC_point.Vair_fps*((sncOut.CDS_total(5)-sncOut.CDS_total(6))/(2*desPerturb.dVair));
linStabOut.CL_U     = FC_point.Vair_fps*((sncOut.CLS_total(5)-sncOut.CLS_total(6))/(2*desPerturb.dVair));
linStabOut.Cm_U     = FC_point.Vair_fps*((sncOut.CMS_total(5)-sncOut.CMS_total(6))/(2*desPerturb.dVair));
% Dynamic Derivatives (zero for now)
linStabOut.CD_alpdot = 0;
linStabOut.CL_alpdot = 0;
linStabOut.Cm_alpdot = 0;
linStabOut.CY_betdot = 0;
linStabOut.Cn_betdot = 0;
linStabOut.Cr_betdot = 0;
% Roll rate derivative
linStabOut.CY_p     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CYS_total(15)-sncOut.CYS_total(16))/(2*desPerturb.dp_rps));
linStabOut.Cn_p     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CNS_total(15)-sncOut.CNS_total(16))/(2*desPerturb.dp_rps));
linStabOut.Cr_p     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CRS_total(15)-sncOut.CRS_total(16))/(2*desPerturb.dp_rps));
% Pitch rate derivative 
linStabOut.CD_q     = ((2*FC_point.Vair_fps)/aeroData.CREF)*((sncOut.CDS_total(17)-sncOut.CDS_total(18))/(2*desPerturb.dq_rps));
linStabOut.CL_q     = ((2*FC_point.Vair_fps)/aeroData.CREF)*((sncOut.CLS_total(17)-sncOut.CLS_total(18))/(2*desPerturb.dq_rps));
linStabOut.Cm_q     = ((2*FC_point.Vair_fps)/aeroData.CREF)*((sncOut.CMS_total(17)-sncOut.CMS_total(18))/(2*desPerturb.dq_rps));
% Yaw rate derivative
linStabOut.CY_r     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CYS_total(19)-sncOut.CYS_total(20))/(2*desPerturb.dr_rps));
linStabOut.Cn_r     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CNS_total(19)-sncOut.CNS_total(20))/(2*desPerturb.dr_rps));
linStabOut.Cr_r     = ((2*FC_point.Vair_fps)/aeroData.BREF)*((sncOut.CRS_total(19)-sncOut.CRS_total(20))/(2*desPerturb.dr_rps));
% Throttle 
linStabOut.CD_dthr  = (sncOut.CDS_total(7)-sncOut.CDS_total(8))/(2*desPerturb.dthr);
linStabOut.CL_dthr  = (sncOut.CLS_total(7)-sncOut.CLS_total(8))/(2*desPerturb.dthr);
linStabOut.Cm_dthr  = (sncOut.CMS_total(7)-sncOut.CMS_total(8))/(2*desPerturb.dthr);
linStabOut.CY_dthr  = (sncOut.CYS_total(7)-sncOut.CYS_total(8))/(2*desPerturb.dthr);
linStabOut.Cn_dthr  = (sncOut.CNS_total(7)-sncOut.CNS_total(8))/(2*desPerturb.dthr);
linStabOut.Cr_dthr  = (sncOut.CRS_total(7)-sncOut.CRS_total(8))/(2*desPerturb.dthr);
% Aileron 
linStabOut.CD_da    = (sncOut.CDS_total(9)-sncOut.CDS_total(10))/(2*desPerturb.da_deg*d2r);
linStabOut.CL_da    = (sncOut.CLS_total(9)-sncOut.CLS_total(10))/(2*desPerturb.da_deg*d2r);
linStabOut.Cm_da    = (sncOut.CMS_total(9)-sncOut.CMS_total(10))/(2*desPerturb.da_deg*d2r);
linStabOut.CY_da    = (sncOut.CYS_total(9)-sncOut.CYS_total(10))/(2*desPerturb.da_deg*d2r);
linStabOut.Cn_da    = (sncOut.CNS_total(9)-sncOut.CNS_total(10))/(2*desPerturb.da_deg*d2r);
linStabOut.Cr_da    = (sncOut.CRS_total(9)-sncOut.CRS_total(10))/(2*desPerturb.da_deg*d2r);
% Elevator
linStabOut.CD_de    = (sncOut.CDS_total(11)-sncOut.CDS_total(12))/(2*desPerturb.de_deg*d2r);
linStabOut.CL_de    = (sncOut.CLS_total(11)-sncOut.CLS_total(12))/(2*desPerturb.de_deg*d2r);
linStabOut.Cm_de    = (sncOut.CMS_total(11)-sncOut.CMS_total(12))/(2*desPerturb.de_deg*d2r);
linStabOut.CY_de    = (sncOut.CYS_total(11)-sncOut.CYS_total(12))/(2*desPerturb.de_deg*d2r);
linStabOut.Cn_de    = (sncOut.CNS_total(11)-sncOut.CNS_total(12))/(2*desPerturb.de_deg*d2r);
linStabOut.Cr_de    = (sncOut.CRS_total(11)-sncOut.CRS_total(12))/(2*desPerturb.de_deg*d2r);
% Rudder
linStabOut.CD_dr    = (sncOut.CDS_total(13)-sncOut.CDS_total(14))/(2*desPerturb.dr_deg*d2r);
linStabOut.CL_dr    = (sncOut.CLS_total(13)-sncOut.CLS_total(14))/(2*desPerturb.dr_deg*d2r);
linStabOut.Cm_dr    = (sncOut.CMS_total(13)-sncOut.CMS_total(14))/(2*desPerturb.dr_deg*d2r);
linStabOut.CY_dr    = (sncOut.CYS_total(13)-sncOut.CYS_total(14))/(2*desPerturb.dr_deg*d2r);
linStabOut.Cn_dr    = (sncOut.CNS_total(13)-sncOut.CNS_total(14))/(2*desPerturb.dr_deg*d2r);
linStabOut.Cr_dr    = (sncOut.CRS_total(13)-sncOut.CRS_total(14))/(2*desPerturb.dr_deg*d2r);

%% Compute MOI in Stability Axis 
Ixxs_slft2 = (massData.ixx_slft2*cosd(FC_point.alpha_deg)*cosd(FC_point.alpha_deg))+...
                (massData.izz_slft2*sind(FC_point.alpha_deg)*sind(FC_point.alpha_deg))-...
                (massData.ixz_slft2*sind(2*FC_point.alpha_deg));
Iyys_slft2 = massData.iyy_slft2;
Izzs_slft2 = (massData.ixx_slft2*sind(FC_point.alpha_deg)*sind(FC_point.alpha_deg))+...
                (massData.izz_slft2*cosd(FC_point.alpha_deg)*cosd(FC_point.alpha_deg))+...
                (massData.ixz_slft2*sind(2*FC_point.alpha_deg));
Ixzs_slft2 = (0.5*(massData.ixx_slft2-massData.iyy_slft2)*sind(2*FC_point.alpha_deg))+...
                (massData.ixz_slft2*cosd(2*FC_point.alpha_deg));

%% Compute Dimensional Stability Derivatives 
% X-Force
linStabOut.dim.X_U      = -1*((FC_point.qbar_psf*aeroData.SREF)/(FC_point.mass_sl*FC_point.Vair_fps))*...
                            ((2*FC_point.CDS_total + linStabOut.CD_U));
linStabOut.dim.X_alpha  = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            (FC_point.CLS_total - linStabOut.CD_alpha);
linStabOut.dim.X_q      = -1*((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CD_q;
linStabOut.dim.X_alpdot = -1*((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CD_alpdot;
linStabOut.dim.X_dthr   = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CD_dthr;
linStabOut.dim.X_da     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CD_da;
linStabOut.dim.X_de     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CD_de;
linStabOut.dim.X_dr     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CD_dr;
% Z-Force
linStabOut.dim.Z_U      = -1*((FC_point.qbar_psf*aeroData.SREF)/(FC_point.mass_sl*FC_point.Vair_fps))*...
                            ((2*FC_point.CLS_total + linStabOut.CL_U));
linStabOut.dim.Z_alpha  = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            (FC_point.CDS_total + linStabOut.CL_alpha);
linStabOut.dim.Z_q      = -1*((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CL_q;
linStabOut.dim.Z_alpdot = -1*((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CL_alpdot;
linStabOut.dim.Z_dthr   = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CL_dthr;
linStabOut.dim.Z_da     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CL_da;
linStabOut.dim.Z_de     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CL_de;
linStabOut.dim.Z_dr     = -1*((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CL_dr;
% Pitch Moment 
linStabOut.dim.M_U      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(Iyys_slft2*FC_point.Vair_fps))*...
                            ((2*FC_point.CMS_total + linStabOut.Cm_U));
linStabOut.dim.M_alpha  = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/Iyys_slft2)*...
                            linStabOut.Cm_alpha;
linStabOut.dim.M_q      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(Iyys_slft2))*(aeroData.CREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cm_q;
linStabOut.dim.M_alpdot = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/(Iyys_slft2))*(aeroData.CREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cm_alpdot;
linStabOut.dim.M_dthr   = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/Iyys_slft2)*...
                            linStabOut.Cm_dthr;
linStabOut.dim.M_da     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/Iyys_slft2)*...
                            linStabOut.Cm_da;
linStabOut.dim.M_de     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/Iyys_slft2)*...
                            linStabOut.Cm_de;
linStabOut.dim.M_dr     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.CREF)/Iyys_slft2)*...
                            linStabOut.Cm_dr;
% Y-Force
linStabOut.dim.Y_beta   = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CY_beta;
linStabOut.dim.Y_p      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CY_p;
linStabOut.dim.Y_r      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/(2*FC_point.mass_sl*FC_point.Vair_fps))*...
                            linStabOut.CY_r;
linStabOut.dim.Y_betdot = 0;
linStabOut.dim.Y_dthr   = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CY_dthr;
linStabOut.dim.Y_da     = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CY_da;
linStabOut.dim.Y_de     = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CY_de;
linStabOut.dim.Y_dr     = ((FC_point.qbar_psf*aeroData.SREF)/FC_point.mass_sl)*...
                            linStabOut.CY_dr;
% Roll Moment 
linStabOut.dim.R_beta   = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*...
                            linStabOut.Cr_beta;
linStabOut.dim.R_p      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*(aeroData.BREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cr_p;
linStabOut.dim.R_r      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*(aeroData.BREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cr_r;
linStabOut.dim.R_betdot = 0;
linStabOut.dim.R_dthr   = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*...
                            linStabOut.Cr_dthr;
linStabOut.dim.R_da     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*...
                            linStabOut.Cr_da;
linStabOut.dim.R_de     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*...
                            linStabOut.Cr_de;
linStabOut.dim.R_dr     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Ixxs_slft2)*...
                            linStabOut.Cr_dr;
% Yaw Moment 
linStabOut.dim.N_beta   = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*...
                            linStabOut.Cn_beta;
linStabOut.dim.N_p      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*(aeroData.BREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cn_p;
linStabOut.dim.N_r      = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*(aeroData.BREF/(2*FC_point.Vair_fps))*...
                            linStabOut.Cn_r;
linStabOut.dim.N_betdot = 0;
linStabOut.dim.N_dthr   = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*...
                            linStabOut.Cn_dthr;
linStabOut.dim.N_da     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*...
                            linStabOut.Cn_da;
linStabOut.dim.N_de     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*...
                            linStabOut.Cn_de;
linStabOut.dim.N_dr     = ((FC_point.qbar_psf*aeroData.SREF*aeroData.BREF)/Izzs_slft2)*...
                            linStabOut.Cn_dr;
%% Compute Primed Derivatives 
den = 1-((Ixzs_slft2^2)/(Ixxs_slft2*Izzs_slft2));
% Roll
linStabOut.dim.Rp_beta  = (linStabOut.dim.R_beta + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_beta))/den;
linStabOut.dim.Rp_p     = (linStabOut.dim.R_p + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_p))/den;
linStabOut.dim.Rp_r     = (linStabOut.dim.R_r + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_r))/den;
linStabOut.dim.Rp_betdot= (linStabOut.dim.R_betdot + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_betdot))/den;
linStabOut.dim.Rp_dthr  = (linStabOut.dim.R_dthr + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_dthr))/den;
linStabOut.dim.Rp_da    = (linStabOut.dim.R_da + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_da))/den;
linStabOut.dim.Rp_de    = (linStabOut.dim.R_de + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_de))/den;
linStabOut.dim.Rp_dr    = (linStabOut.dim.R_dr + ((Ixzs_slft2/Ixxs_slft2)*linStabOut.dim.N_dr))/den;
% Yaw
linStabOut.dim.Np_beta  = (linStabOut.dim.N_beta + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_beta))/den;
linStabOut.dim.Np_p     = (linStabOut.dim.N_p + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_p))/den;
linStabOut.dim.Np_r     = (linStabOut.dim.N_r + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_r))/den;
linStabOut.dim.Np_betdot= (linStabOut.dim.N_betdot + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_betdot))/den;
linStabOut.dim.Np_dthr  = (linStabOut.dim.N_dthr + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_dthr))/den;
linStabOut.dim.Np_da    = (linStabOut.dim.N_da + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_da))/den;
linStabOut.dim.Np_de    = (linStabOut.dim.N_de + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_de))/den;
linStabOut.dim.Np_dr    = (linStabOut.dim.N_dr + ((Ixzs_slft2/Izzs_slft2)*linStabOut.dim.R_dr))/den;

%% Formulate Longitudinal Linear System 
% Compute A-matrix Components
lon_A11     = linStabOut.dim.Z_alpha/(FC_point.Vair_fps-linStabOut.dim.Z_alpdot);
lon_A12     = (FC_point.Vair_fps+linStabOut.dim.Z_q)/(FC_point.Vair_fps-linStabOut.dim.Z_alpdot);
lon_A13     = linStabOut.dim.Z_U/(FC_point.Vair_fps-linStabOut.dim.Z_alpdot);
lon_A14     = (-1*FC_point.gd_fps2*sind(FC_point.gamma_deg))/(FC_point.Vair_fps-linStabOut.dim.Z_alpdot);
lon_A21     = linStabOut.dim.M_alpha + ((linStabOut.dim.M_alpdot*linStabOut.dim.Z_alpha)/den);
lon_A22     = linStabOut.dim.M_q + (linStabOut.dim.M_alpdot*(FC_point.Vair_fps+linStabOut.dim.Z_q)/den);
lon_A23     = linStabOut.dim.M_U + ((linStabOut.dim.M_alpdot*linStabOut.dim.Z_U)/den);
lon_A24     = (-1*linStabOut.dim.M_alpdot*FC_point.gd_fps2*sind(FC_point.gamma_deg))/den;
lon_A31     = linStabOut.dim.X_alpha;
lon_A32     = 0;
lon_A33     = linStabOut.dim.X_U;
lon_A34     = -1*FC_point.gd_fps2*cosd(FC_point.gamma_deg);
lon_A41     = 0;
lon_A42     = 1;
lon_A43     = 0;
lon_A44     = 0;
% Formulate A-Matrix
linStabOut.sys.lon.A    = [lon_A11 lon_A12 lon_A13 lon_A14;...
                           lon_A21 lon_A22 lon_A23 lon_A24;...
                           lon_A31 lon_A32 lon_A33 lon_A34;...
                           lon_A41 lon_A42 lon_A43 lon_A44];
% Compute Stability Characteristics 
linStabOut.sys.lon.chars = modalAnalysis(linStabOut.sys.lon.A);
% Compute B-Matrix Components
lon_B11     = linStabOut.dim.Z_dthr/den;
lon_B12     = linStabOut.dim.Z_da/den;
lon_B13     = linStabOut.dim.Z_de/den;
lon_B14     = linStabOut.dim.Z_dr/den;
lon_B21     = linStabOut.dim.M_dthr+((linStabOut.dim.M_alpdot*linStabOut.dim.Z_dthr)/den);
lon_B22     = linStabOut.dim.M_da+((linStabOut.dim.M_alpdot*linStabOut.dim.Z_da)/den);
lon_B23     = linStabOut.dim.M_de+((linStabOut.dim.M_alpdot*linStabOut.dim.Z_de)/den);
lon_B24     = linStabOut.dim.M_dr+((linStabOut.dim.M_alpdot*linStabOut.dim.Z_dr)/den);
lon_B31     = linStabOut.dim.X_dthr;
lon_B32     = linStabOut.dim.X_da;
lon_B33     = linStabOut.dim.X_de;
lon_B34     = linStabOut.dim.X_dr;
lon_B41     = 0;
lon_B42     = 0;
lon_B43     = 0;
lon_B44     = 0;
% Formulate B-Matrix
linStabOut.sys.lon.B    = [lon_B11 lon_B12 lon_B13 lon_B14;...
                           lon_B21 lon_B22 lon_B23 lon_B24;...
                           lon_B31 lon_B32 lon_B33 lon_B34;...
                           lon_B41 lon_B42 lon_B43 lon_B44];
%% Formulate Lat/Dir Linear System 
% Compute A-Matrix Components 
ld_A11      = linStabOut.dim.Y_beta/FC_point.Vair_fps;
ld_A12      = FC_point.gd_fps2*cosd(FC_point.theta_deg)/FC_point.Vair_fps;
ld_A13      = linStabOut.dim.Y_p/FC_point.Vair_fps;
ld_A14      = (linStabOut.dim.Y_r - FC_point.Vair_fps)/FC_point.Vair_fps;
ld_A21      = 0;
ld_A22      = 0;
ld_A23      = cosd(FC_point.gamma_deg)/cosd(FC_point.theta_deg);
ld_A24      = sind(FC_point.gamma_deg)/cosd(FC_point.theta_deg);
ld_A31      = linStabOut.dim.Rp_beta;
ld_A32      = 0;
ld_A33      = linStabOut.dim.Rp_p;
ld_A34      = linStabOut.dim.Rp_r;
ld_A41      = linStabOut.dim.Np_beta;
ld_A42      = 0;
ld_A43      = linStabOut.dim.Np_p;
ld_A44      = linStabOut.dim.Np_r;
% Formulate A-Matrix 
linStabOut.sys.lat.A    = [ld_A11 ld_A12 ld_A13 ld_A14;...
                           ld_A21 ld_A22 ld_A23 ld_A24;...
                           ld_A31 ld_A32 ld_A33 ld_A34;...
                           ld_A41 ld_A42 ld_A43 ld_A44];
% Compute Stability Characteristics 
linStabOut.sys.lat.chars = modalAnalysis(linStabOut.sys.lat.A);
% Compute B-Matrix Components 
ld_B11      = linStabOut.dim.Y_dthr/FC_point.Vair_fps;
ld_B12      = linStabOut.dim.Y_da/FC_point.Vair_fps;
ld_B13      = linStabOut.dim.Y_de/FC_point.Vair_fps;
ld_B14      = linStabOut.dim.Y_dr/FC_point.Vair_fps;
ld_B21      = 0;
ld_B22      = 0;
ld_B23      = 0;
ld_B24      = 0;
ld_B31      = linStabOut.dim.Rp_dthr;
ld_B32      = linStabOut.dim.Rp_da;
ld_B33      = linStabOut.dim.Rp_de;
ld_B34      = linStabOut.dim.Rp_dr;
ld_B41      = linStabOut.dim.Np_dthr;
ld_B42      = linStabOut.dim.Np_da;
ld_B43      = linStabOut.dim.Np_de;
ld_B44      = linStabOut.dim.Np_dr;
% Formulate B-Matrix 
linStabOut.sys.lat.B    = [ld_B11 ld_B12 ld_B13 ld_B14;...
                           ld_B21 ld_B22 ld_B23 ld_B24;...
                           ld_B31 ld_B32 ld_B33 ld_B34;...
                           ld_B41 ld_B42 ld_B43 ld_B44];
end
