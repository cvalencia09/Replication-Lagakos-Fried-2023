%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: tables.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Creates tables for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear all;

%Load paths 
load paths

%% Table 1: Externally calibrated parameters
load cal_NGA; 
cd(tablePath);
FID = fopen('direct.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Parameter Values: Externally Calibrated Parameters} \\label{tab:direct}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{c c c c c c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '  	$\\eta$ & $\\alpha$ & $A^T$ & $\\phi$ &  $\\delta$ & $A^G$ & $\\psi$ & $\\sigma$ & $\\beta$   \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' %8.2f & %8.2f & %8.0f & %8.1f & %8.2f & %8.0f & %8.1f & %8.0f & %8.2f     \\\\  '...
    , eta, alpha, A, phi, delta, Ag, psi, sigma, beta );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Table 2: Internally calibrated parameters
load cal_NGA; 
cd(tablePath);
FID = fopen('mom.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Parameter Values: Internally Calibrated Parameters} \\label{tab:mom}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{c c c c c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '  	$\\mu$ & $\\chi$ & $A^S$ & $P^G$ & $\\Omega$  & $\\lambda$ & $\\zeta$ &$\\rho$  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\  '...
    , mu, chi, As, Pg, Omega, lambda, zeta, gamma);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Table 3: Semi-Elasticity of Generator Ownership With Respect to Firm Size
%Data For "Data" Column are stored in /output/results_micro.xlsx

load ss1_GHA; 
load ss1_UGA; 
load ss1_NGA; 
load ss1_TZA; 

cd(tablePath);
FID = fopen('elasticity.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Semi-Elasticity of Generator Ownership With Respect to Firm Size} \\label{tab:elasticity}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '  	& Data & Model  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Ghana & 0.12 & %8.2f  \\\\[-0.5ex]', elasticity_GHA1);
fprintf(FID, ' & [0.07, 0.17]  \\\\');
fprintf(FID, ' Nigeria & 0.04 & %8.2f  \\\\[-0.5ex]', elasticity_NGA1);
fprintf(FID, ' & 	[-0.04,   0.12]  \\\\');
fprintf(FID, ' Tanzania & 0.06 & %8.2f  \\\\[-0.5ex]', elasticity_TZA1);
fprintf(FID, ' &  [-0.004,   0.13]  \\\\');
fprintf(FID, ' Uganda & 0.05 & %8.2f  \\\\[-0.5ex]', elasticity_UGA1);
fprintf(FID, ' & [-0.07,  0.16] \\\\');
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);



%% Table 4: Grid electricity price, capital and supply
load ss1_GHA; load ss2_GHA;
load ss1_UGA; load ss2_UGA;
load ss1_NGA; load ss2_TZA;
load ss1_TZA; load ss2_NGA;

cd(tablePath);
FID = fopen('Kg_pg.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Grid Electricity Price, Capital, and Supply \\\\ (Value relative to initial steady state)} \\label{Kg_pg}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '  & Grid price: $P^G$ & Grid capital: $K^G$ & Grid electricity: $E^G$   \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Ghana &  %8.2f & %8.2f & %8.2f    \\\\  '...
    , Pg_GHA2/Pg_GHA1, Kg_GHA2/Kg_GHA1, Eg_GHA2/Eg_GHA1 );
fprintf(FID, 'Nigeria  &  %8.2f & %8.2f & %8.2f    \\\\  ' ...
    , Pg_NGA2/Pg_NGA1, Kg_NGA2/Kg_NGA1, Eg_NGA2/Eg_NGA1 );
fprintf(FID, 'Tanzania &  %8.2f & %8.2f & %8.2f   \\\\ '...
    , Pg_TZA2/Pg_TZA1, Kg_TZA2/Kg_TZA1, Eg_TZA2/Eg_TZA1 );
fprintf(FID, 'Uganda &   %8.2f & %8.2f & %8.2f \\\\  '...
    , Pg_UGA2/Pg_UGA1, Kg_UGA2/Kg_UGA1, Eg_UGA2/Eg_UGA1 );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Table 5: Effects of Eliminating Outages on Macro Aggregates

load ss1_GHA; load ss2_GHA;
load ss1_UGA; load ss2_UGA;
load ss1_NGA; load ss2_TZA;
load ss1_TZA; load ss2_NGA;

TFP_GHA2 = Y_GHA2/(K_GHA2^alpha); 
TFP_NGA1 = Y_NGA1/(K_NGA1^alpha); TFP_NGA2 = Y_NGA2/(K_NGA2^alpha); 
TFP_TZA1 = Y_TZA1/(K_TZA1^alpha); TFP_TZA2 = Y_TZA2/(K_TZA2^alpha); 
TFP_UGA1 = Y_UGA1/(K_UGA1^alpha); TFP_UGA2 = Y_UGA2/(K_UGA2^alpha); 

TFP2_GHA1 = Y_GHA1/((K_GHA1 + Nm_GHA1*Omega)^alpha); TFP2_GHA2 = Y_GHA2/((K_GHA2 + Nm_GHA2*Omega)^alpha); 
TFP2_NGA1 = Y_NGA1/((K_NGA1 + Nm_NGA1*Omega)^alpha); TFP2_NGA2 = Y_NGA2/((K_NGA2 + Nm_NGA2*Omega)^alpha); 
TFP2_TZA1 = Y_TZA1/((K_TZA1 + Nm_TZA1*Omega)^alpha); TFP2_TZA2 = Y_TZA2/((K_TZA2 + Nm_TZA2*Omega)^alpha); 
TFP2_UGA1 = Y_UGA1/((K_UGA1 + Nm_UGA1*Omega)^alpha); TFP2_UGA2 = Y_UGA2/((K_UGA2 + Nm_UGA2*Omega)^alpha); 

cd(tablePath);
FID = fopen('agg.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Effects of Eliminating Outages on Macro Aggregates\\\\ (percent change from the initial steady state)} \\label{agg}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, '   & Ghana & Nigeria & Tanzania & Uganda    \\\\  \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Modern entrepreneurs: $Q^M$ & %8.0f & %8.0f & %8.0f  & %8.0f    \\\\  '...
    , (Nm_GHA2/Nm_GHA1-1)*100, (Nm_NGA2/Nm_NGA1-1)*100, ...
    (Nm_TZA2/Nm_TZA1-1)*100, (Nm_UGA2/Nm_UGA1-1)*100 );
fprintf(FID, ' Modern labor: $N^M$ & %8.0f & %8.0f & %8.0f  & %8.0f    \\\\  '...
    ,(Lm_GHA2/Lm_GHA1-1)*100, (Lm_NGA2/Lm_NGA1-1)*100, ...
    (Lm_TZA2/Lm_TZA1-1)*100, (Lm_UGA2/Lm_UGA1-1)*100 );
fprintf(FID, ' Wage rate: $W$ & %8.0f & %8.0f & %8.0f  & %8.0f     \\\\  '...
    , (W_GHA2/W_GHA1-1)*100, (W_NGA2/W_NGA1-1)*100, ...
    (W_TZA2/W_TZA1-1)*100, (W_UGA2/W_UGA1-1)*100 );
fprintf(FID, 'Modern productive capital: $K^M$ & %8.0f & %8.0f & %8.0f  & %8.0f     \\\\ '...
    ,  (Km_GHA2/Km_GHA1-1)*100, (Km_NGA2/Km_NGA1-1)*100, ...
    (Km_TZA2/Km_TZA1-1)*100, (Km_UGA2/Km_UGA1-1)*100 );
fprintf(FID, 'Traditional capital: $K^T$ & %8.0f & %8.0f & %8.0f  & %8.0f     \\\\  '...
    ,  (Kt_GHA2/Kt_GHA1-1)*100, (Kt_NGA2/Kt_NGA1-1)*100, ...
    (Kt_TZA2/Kt_TZA1-1)*100, (Kt_UGA2/Kt_UGA1-1)*100 );
fprintf(FID, 'Aggregate capital stock: $K$ & %8.0f & %8.0f & %8.0f  & %8.0f     \\\\  '...
    ,  (K_GHA2/K_GHA1-1)*100, (K_NGA2/K_NGA1-1)*100, ...
    (K_TZA2/K_TZA1-1)*100, (K_UGA2/K_UGA1-1)*100 );
fprintf(FID, 'Measured TFP: $Y/(\\tilde{K}^\\alpha N^{1-\\alpha})$ & %8.0f & %8.0f & %8.0f  & %8.0f    \\\\  '...
    ,  (TFP2_GHA2/TFP2_GHA1-1)*100, (TFP2_NGA2/TFP2_NGA1-1)*100, ...
    (TFP2_TZA2/TFP2_TZA1-1)*100, (TFP2_UGA2/TFP2_UGA1-1)*100 );
fprintf(FID, 'Consumption: $C$ & %8.0f & %8.0f & %8.0f  & %8.0f    \\\\  '...
    , (C_GHA2/C_GHA1-1)*100, (C_NGA2/C_NGA1-1)*100, ...
    (C_TZA2/C_TZA1-1)*100, (C_UGA2/C_UGA1-1)*100 );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Table 6: Sensitivity to Parameter Values and Targets

load partial_eq
load ss1_NGA; load ss2_NGA;
load ss1_NGA_psi6; load ss2_NGA_psi6; 
load ss1_NGA_psi8; load ss2_NGA_psi8;
load ss1_NGA_sp; load ss2_NGA_sp; 
load ss1_NGA_bp; load ss2_NGA_bp;
load ss1_NGA_hs; load ss2_NGA_hs; 
load ss1_NGA_ls; load ss2_NGA_ls; 
load ss1_NGA_lc; load ss2_NGA_lc; 
load ss1_NGA_hc; load ss2_NGA_hc;
load ss1_NGA_q3; load ss2_NGA_q3; 
load ss1_NGA_q1; load ss2_NGA_q1;

cd(tablePath);
FID = fopen('sensitivity.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Sensitivity to Parameter Values and Targets} \\label{sensitivity}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, ' & \\multicolumn{3}{c}{Percent increase in output from initial steady state} \\\\ \n');
fprintf(FID, '\\cline{2-4} \n');
fprintf(FID, '  & Long run, G.E.& Short run, P.E.   & Ratio\\\\ [.5ex] \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\emph{Capital share in grid electricity production: $\\psi$}} \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' $\\psi = 0.6$ &  %8.1f & %8.1f   & %8.1f  \\\\  '  , (Y_NGA_psi62/Y_NGA_psi61-1)*100, SR_NGA_psi6,  ...
    (Y_NGA_psi62/Y_NGA_psi61-1)*100/SR_NGA_psi6);
fprintf(FID, '$\\psi = 0.7$ & %8.1f & %8.1f   & %8.1f   \\\\  ' , (Y_NGA2/Y_NGA1-1)*100, SR_NGA, ...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA);
fprintf(FID, '$\\psi = 0.8$ & %8.1f & %8.1f   & %8.1f \\\\  '  , (Y_NGA_psi82/Y_NGA_psi81-1)*100, SR_NGA_psi8,...
    (Y_NGA_psi82/Y_NGA_psi81-1)*100/SR_NGA_psi8);
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\emph{Modern productivity boost: $\\phi$}} \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' $\\phi = 0.3$ &  %8.1f & %8.1f  & %8.1f  \\\\  '  , (Y_NGA_sp2/Y_NGA_sp1-1)*100, SR_NGA_sp, ...
    (Y_NGA_sp2/Y_NGA_sp1-1)*100/SR_NGA_sp);
fprintf(FID, '$\\phi = 0.4$ & %8.1f & %8.1f & %8.1f \\\\  ' , (Y_NGA2/Y_NGA1-1)*100, SR_NGA, ...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA);
fprintf(FID, '$\\phi = 0.5$ & %8.1f & %8.1f & %8.1f \\\\  '  , (Y_NGA_bp2/Y_NGA_bp1-1)*100, SR_NGA_bp,...
    (Y_NGA_bp2/Y_NGA_bp1-1)*100/SR_NGA_bp);
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\emph{Electricity share of modern output}} \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Share = 0.06 & %8.1f & %8.1f & %8.1f \\\\  '  , (Y_NGA_ls2/Y_NGA_ls1-1)*100, SR_NGA_ls, ...
     (Y_NGA_ls2/Y_NGA_ls1-1)*100/SR_NGA_ls);
fprintf(FID, 'Share = 0.074 &  %8.1f & %8.1f & %8.1f  \\\\  ' , (Y_NGA2/Y_NGA1-1)*100, SR_NGA, ...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA);
fprintf(FID, 'Share = 0.09 & %8.1f & %8.1f & %8.1f  \\\\  '  , (Y_NGA_hs2/Y_NGA_hs1-1)*100, SR_NGA_hs,...
    (Y_NGA_hs2/Y_NGA_hs1-1)*100/SR_NGA_hs);
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\emph{Self-generation costs}} \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '10 percent decrease in costs &  %8.1f & %8.1f & %8.1f  \\\\  '  , (Y_NGA_lc2/Y_NGA_lc1-1)*100, SR_NGA_lc,...
    (Y_NGA_lc2/Y_NGA_lc1-1)*100/SR_NGA_lc);
fprintf(FID, 'Baseline costs & %8.1f & %8.1f & %8.1f  \\\\  ' , (Y_NGA2/Y_NGA1-1)*100, SR_NGA, ...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA);
fprintf(FID, '10 percent increase in costs &  %8.1f & %8.1f & %8.1f \\\\  '  , (Y_NGA_hc2/Y_NGA_hc1-1)*100, SR_NGA_hc,...
     (Y_NGA_hc2/Y_NGA_hc1-1)*100/SR_NGA_hc);
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multicolumn{3}{l}{\\emph{Fraction of modern firms with no outages}} \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Fraction = 0.1 &  %8.1f & %8.1f & %8.1f  \\\\  '  , (Y_NGA_q12/Y_NGA_q11-1)*100, SR_NGA_q1,...
     (Y_NGA_q12/Y_NGA_q11-1)*100/SR_NGA_q1);
fprintf(FID, 'Fraction  = 0.2 & %8.1f & %8.1f & %8.1f \\\\  ' , (Y_NGA2/Y_NGA1-1)*100, SR_NGA, ...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA);
fprintf(FID, 'Fraction = 0.3 &  %8.1f & %8.1f & %8.1f  \\\\  '  , (Y_NGA_q32/Y_NGA_q31-1)*100, SR_NGA_q3,...
    (Y_NGA_q32/Y_NGA_q31-1)*100/SR_NGA_q3);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);



%% Table C.1 Model fit
load cal_NGA
%find_elas =0;
%solveModel

%Notes: hard-coded in share_modern_labor and share_modern_labor_targ to fix
%rounding issue. To view that they are the same run:
%share_modern_labor-share_modern_labor_targ

cd(tablePath);
FID = fopen('model_fit.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Model Fit} \\label{model_fit}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c }\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' {Moment} & Model  & Target \\\\ [.5ex] \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '(variable cost of self-generation)/(grid electricity price)  & %8.2f & %8.2f    \\\\ ', ps_pg, ps_pg_targ);
fprintf(FID, '(average cost of self-generation)/(grid electricity price)  & %8.2f & %8.2f    \\\\ ', ac_self_grid, ac_self_grid_targ);
fprintf(FID, 'Fraction of self-generated electricity  & %8.2f & %8.2f    \\\\  ', frac_self, frac_self_targ);
fprintf(FID, 'Modern electricity share & %8.2f & %8.2f    \\\\  ', elec_share_modern_v928, elec_share_modern_targ);
fprintf(FID, 'Fraction of modern labor  & %8.2f & %8.2f    \\\\  ', 0.63, 0.63);
fprintf(FID, 'Fraction of modern entrepreneurs  & %8.2f & %8.2f    \\\\  ', frac_modern_firms, frac_modern_firm_targ);
fprintf(FID, 'Fraction of modern entrepreneurs with a generator & %8.2f & %8.2f    \\\\  ', Nm_gen_Nm_q0,  frac_modern_generators_target);
fprintf(FID, 'Fraction of modern entrepreneurs without outages & %8.2f & %8.2f    \\\\  ', Nm_q1_Nm,  Nm_q1_Nm_targ);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Table C.2: Elasticities of Moments to Parameters
load perturbations; 

cd(tablePath);
FID = fopen('cal_matrix.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\renewcommand{\\arraystretch}{2} \n');
fprintf(FID, '\\setlength{\\tabcolsep}{3mm}\n');
fprintf(FID, '\\vspace{3mm} \n');
fprintf(FID, '\\caption{Elasticities of Moments to Parameters} \\label{tab:cal_matrix}\n');
fprintf(FID, '\\begin{tabular}{l c c c c c c c c }\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '& $\\frac{P_S}{P_G}$ & Avg cost$\\left(\\frac{self}{grid}\\right)$  & $\\frac{Self}{Grid}$ & $\\frac{L_m}{L}$ & $\\frac{N_m}{N}$ & $\\frac{N_{m,gen}}{N_m}$ & $\\frac{N_{m,q1}}{N_m}$ & elec share \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '  $\\chi$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(1,:) );
fprintf(FID, '  $A^S$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(2,:) );
fprintf(FID, '  $P^G$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(3,:) );
fprintf(FID, '  $\\lambda$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(4,:) );
fprintf(FID, '  $\\Omega$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(5,:) );
fprintf(FID, '  $\\zeta$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(6,:) );
fprintf(FID, '  $\\rho$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(7,:) );
fprintf(FID, '  $\\mu$ & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f\\\\ ',temp_sorted2(8,:) );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Table C.3: Country-specific targets

M = xlsread(strcat(datapath, '\cross_country_moments.xlsx'));
frac_no_outages_vec = M(1:end,1); 
frac_self_vec =M(1:end,2); 
Nm_gen_Nm_q0_vec =M(1:end,3); 
ps_pg_vec = M(1:end,4); ac_pg_vec = M(1:end, 5); 
rel_y_pop_vec = M(1:end, 6); 

cd(tablePath);
FID = fopen('country_moments.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Country Specific Targets} \\label{country_moments}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' {Country} & AC$^{self}$/AC$^{grid}$  & MC$^{self}$/MC$^{grid}$ & $E^s/E$ & $Y/Y_{NGA}$ & $N_{m,q=1}/N_m$ &$N_{m,gen}/N_{m,q=0}$ \\\\  \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Ghana  & %8.2f & %8.2f & %8.2f   & %8.2f  & %8.2f& %8.2f\\\\  ',...
    ac_pg_vec(2), ps_pg_vec(2), frac_self_vec(2), rel_y_pop_vec(2), frac_no_outages_vec(2),Nm_gen_Nm_q0_vec(2) );
fprintf(FID, 'Nigeria  & %8.2f & %8.2f & %8.2f  & %8.2f   & %8.2f& %8.2f\\\\  ',...
    ac_pg_vec(1), ps_pg_vec(1), frac_self_vec(1), rel_y_pop_vec(1), frac_no_outages_vec(1), Nm_gen_Nm_q0_vec(1));
fprintf(FID, 'Tanzania & %8.2f & %8.2f & %8.2f  & %8.2f  & %8.2f & %8.2f\\\\  ', ...
    ac_pg_vec(3), ps_pg_vec(3), frac_self_vec(3), rel_y_pop_vec(3), frac_no_outages_vec(3), Nm_gen_Nm_q0_vec(3));
fprintf(FID, 'Uganda & %8.2f & %8.2f & %8.2f  & %8.2f  & %8.2f& %8.2f\\\\  ', ...
    ac_pg_vec(4), ps_pg_vec(4), frac_self_vec(4), rel_y_pop_vec(4), frac_no_outages_vec(4), Nm_gen_Nm_q0_vec(4));
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Table C.4: Country-specific parameter values

cd(tablePath);
FID = fopen('country_params.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Country Specific Parameters} \\label{country_params}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' \\multirow{2}{*}{Country} & Price of grid  & Self-gen  &Self-gen   & Trad. TFP: & Gumbel  & Fraction with  \\\\ [-.5ex] \n');
fprintf(FID, '  & electricity: $P^G$ & Leontief: $\\chi$  & TFP: $A^S$ & $A^T$ & scale: $\\zeta$ & no outages: $\\rho$ \\\\ [.5ex] \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Ghana  & %8.2f & %8.2f & %8.2f  & %8.2f   & %8.2f  & %8.2f  \\\\  ',...
    Pg_GHA1,  chi_GHA1,  As_GHA1,  A_GHA1,  zeta_GHA1, gamma_GHA1);
fprintf(FID, 'Nigeria  & %8.2f & %8.2f & %8.2f  & %8.2f  & %8.2f  & %8.2f   \\\\  ',...
    Pg_NGA1,  chi_NGA1, As_NGA1, A_NGA1, zeta_NGA1, gamma_NGA1);
fprintf(FID, 'Tanzania  & %8.2f & %8.2f & %8.2f  & %8.2f  & %8.2f  & %8.2f  \\\\  ',...
    Pg_TZA1,  chi_TZA1, As_TZA1, A_TZA1, zeta_TZA1, gamma_TZA1);
fprintf(FID, 'Uganda  & %8.2f & %8.2f & %8.2f  & %8.2f  & %8.2f   & %8.2f  \\\\  ',...
    Pg_UGA1,chi_UGA1, As_UGA1,   A_UGA1, zeta_UGA1, gamma_UGA1);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Table D.1: Equilibrium values of macro aggregates in the initial steady state

load ss1_GHA; 
load ss1_UGA; 
load ss1_NGA;
load ss1_TZA; 

cd(tablePath);
FID = fopen('agg_initial.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Equilibrium Values of Macro-Aggregates in the Initial Steady State} \\label{agg_initial}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, '  & Ghana & Nigeria & Tanzania & Uganda    \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Grid-electricity price: $P^G$ & %8.2f & %8.2f & %8.2f  & %8.2f   \\\\  ',...
 Pg_GHA1, Pg_NGA1, Pg_TZA1, Pg_UGA1);
fprintf(FID, 'Grid-electricity capital: $K^G$ & %8.2f & %8.2f & %8.2f  & %8.2f     \\\\  ',...
  Kg_GHA1, Kg_NGA1, Kg_TZA1, Kg_UGA1);
fprintf(FID, 'Grid-electricity supply: $E^G$ & %8.2f & %8.2f & %8.2f  & %8.2f     \\\\ ',...
 Eg_GHA1, Eg_NGA1, Eg_TZA1, Eg_UGA1);
fprintf(FID, 'Fraction of modern labor: $N^M$ & %8.2f & %8.2f & %8.2f  & %8.2f     \\\\  ',...
  Lm_GHA1, Lm_NGA1, Lm_TZA1, Lm_UGA1);
fprintf(FID, 'Fraction of modern entrepreneurs: $Q^M$ & %8.2f & %8.2f & %8.2f  & %8.2f     \\\\  ',...
   Nm_GHA1, Nm_NGA1, Nm_TZA1, Nm_UGA1);
fprintf(FID, 'Probability of grid power: $v$ & %8.2f & %8.2f & %8.2f  & %8.2f     \\\\ ',...
   v_GHA1, v_NGA1, v_TZA1, v_UGA1);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Table E.1. Sensitivity to Generator Fixed Cost
cd(bpath);
load partial_eq
load ss1_NGA; load ss2_NGA;
load ss1_NGA_ups05; load ss2_NGA_ups05;
load ss1_NGA_ups10; load ss2_NGA_ups10;
load ss1_NGA_ups15; load ss2_NGA_ups15;
load ss1_NGA_ups20; load ss2_NGA_ups20;

cd(tablePath);
FID = fopen('generator_fixed_cost.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Sensivity to Generator Fixed Cost: Nigeria} \\label{generator_fixed_cost}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, ' & \\multicolumn{5}{c}{Generator fixed cost} \\\\ \n');
fprintf(FID, '\\cline{2-6} \n');
fprintf(FID, '  & 0 & 0.05 & 0.1 & 0.15 & 0.2    \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Short-run partial-equilibrium effect & %8.1f & %8.1f & %8.1f  & %8.1f & %8.1f  \\\\  ',...
    SR_NGA, SR_NGA_ups05, SR_NGA_ups10, SR_NGA_ups15, SR_NGA_ups20);
fprintf(FID, 'Long-run general-equilibrium effect & %8.1f & %8.1f & %8.1f  & %8.1f & %8.1f  \\\\  ',...
    (Y_NGA2/Y_NGA1-1)*100, (Y_NGA_ups052/Y_NGA_ups051-1)*100,(Y_NGA_ups102/Y_NGA_ups101-1)*100, (Y_NGA_ups152/Y_NGA_ups151-1)*100, ...
    (Y_NGA_ups202/Y_NGA_ups201-1)*100);
fprintf(FID, 'Ratio & %8.1f & %8.1f & %8.1f  & %8.1f & %8.1f  \\\\  ',...
    (Y_NGA2/Y_NGA1-1)*100/SR_NGA, (Y_NGA_ups052/Y_NGA_ups051-1)*100/SR_NGA_ups05, ...
    (Y_NGA_ups102/Y_NGA_ups101-1)*100/SR_NGA_ups10, (Y_NGA_ups152/Y_NGA_ups151-1)*100/SR_NGA_ups15, ...
    (Y_NGA_ups202/Y_NGA_ups201-1)*100/SR_NGA_ups20);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Table E.2: Effects of higher capital and productivity
load decomp_2
cd(tablePath);
FID = fopen('decomp2.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\caption{Effects of Higher Capital and Productivity \\\\ (percent increase in output per worker from the initial steady state)} \\label{decomp2}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c}\\hline \n');
fprintf(FID, '   & Ghana & Nigeria & Tanzania & Uganda    \\\\  \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Short-run partial equilibrium & %8.1f & %8.1f & %8.1f  & %8.1f    \\\\  ', PE_country);
fprintf(FID, ' Add exogenous capital increase & %8.1f & %8.1f & %8.1f  & %8.1f    \\\\  ', Y_SR_K);
fprintf(FID, ' Add exogenous productivity increase & %8.1f & %8.1f & %8.1f  & %8.1f    \\\\  ', Y_SR_prod);
fprintf(FID, ' Long-run general equilibrium & %8.1f & %8.1f & %8.1f  & %8.1f    \\\\  ', GE_country);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);
















