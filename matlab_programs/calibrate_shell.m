%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: calibrate_shell.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Calibrates model for each country, for sensitivity exercise, and
% extension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
clc;
close all;

load paths
cd(bpath);

options = optimset('Algorithm','trust-region-dogleg','TolFun',10e-15,'TolX',10e-15, ...
    'MaxFunEvals', 10000, 'MaxIter', 10000, 'Display','none');

%% Baseline Calibration for Nigeria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exogenously determined parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.33;       % capital share
beta =  0.96;       % discount factor 
delta = 0.06;       % depreciation rate
eta = 0.85;         % Span of control
sigma = 2;          % CRAA utility coefficient
N =  1;             % Number of firms
psi =0.7;           % Capital share in grid electricity
nz = 10000;         % Number of simulated entrepreneurs
zmax = 50;          % Max entrepreneurial productivity
tau =0;             % Tax on grid electricity
ups =0.0;           % Fixed cost of operating a generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Technologies and endowments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ag = 1;            % Grid electricity productivity
A = 1;             % Traditional Productivity  
Ls = 1;            % Labor supply
phi = .4;          % Modern productivity boost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess for internal calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pg= 0.1056;          % Grid capital price
lambda = 1.7330;     % Exponent on the pareto distribution
Omega = 0.5373;      % Modern entry cost
chi = 2.1444;        % Fuel efficiency  
mu =  1.4769;        % Electricity efficiency
As =  1.02;          % Self-generation productivity 
zeta = 0.2435;       % Gumbel variance 
gamma =  0.0755;     % Probability an entrepreneur never expereinces an outage (rho in paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Targets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
share_modern_labor_targ = 0.625;      % Fraction of labor in the modern sector (/data/programs/Nigeria_firms.xlsx)
frac_modern_firm_targ = 0.297;        % Fraction of firms in the modern sector(/data/programs/Nigeria_firms.xlsx)
ps_pg_targ = 4.33;                    % Ratio of marginal costs: self-generated electricity to grid (/data/programs/ac_mc.xlsx)
frac_self_targ  =.59;                 % Fraction of electricity modern firms with generators generate themselves (/data/output/results_micro.xlsx)
elec_share_modern_targ = 0.074;       % Electricity share of modern output when v =0.928 (ACO 2016)
ac_self_grid_targ = 5.51;             % Ratio of average costs: self-generated electricity to grid (/data/programs/ac_mc.xlsx)
frac_modern_generators_target =0.82;  % Fraction of firms that experience outages that have access to a generator (/data/output/results_micro.xlsx)
Nm_q1_Nm_targ = 0.20;                  % Fraction of modern firms that never experience an outage (/data/output/results_micro.xlsx)

targets  =[ps_pg_targ, frac_self_targ, share_modern_labor_targ, ...
    frac_modern_firm_targ,  ac_self_grid_targ, frac_modern_generators_target,...
    Nm_q1_Nm_targ, elec_share_modern_targ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flags for Nigerian calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indicate = 0; % Equals 1 for non-Nigaria calibration
decomp = 0;   % Equals 1 for the decomposition exercise in Figure 6, 
                % equals 2 for the decomposition exercise in Appendix Table
                % E.2
find_elas =0; % Equals 1 to calculate the elasticity of generator ownership with respect to firm size

% Initial guess for errors_ss [W,v]; 
x0 =[1.0892    0.4100];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose parameters to so that moments in model match their empirical
% targets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vector of parameter values
params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
         ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];
     
%Internally determined parameters
p0 = [Pg, lambda, As, zeta, gamma, Omega, chi, mu];
paramsFixed = params(9:end);

%Calibrate model
[dist1, pStar, sol] = calibrate(p0, paramsFixed, targets, x0, options, indicate);

% Report calibrated moments and targets
p0  = pStar; 
params =[pStar, paramsFixed];
[sol, fval, exitflag] = fsolve(@errors_ss,sol,options,params, indicate);
[errorsTry, moments1, guess1] = errors_ss(sol,params,  indicate);
[xstar1, fval, exitflag] = fsolve(@errors_sspg, guess1,options,params, indicate, .928);
[errorsTry, moments2] = errors_sspg(xstar1,params, indicate, 0.928);
moments = [moments1(1:7), moments2];
elec_share_modern_v928 = moments2;  
fprintf('Moments and targets')
[moments; targets]

%Find elasticity of generator ownership with respect to firm size
W = sol(1); v  = sol(2);
find_elas =1; 
solveModel;

%Save intital guess
x0_initial = sol; 
x0 = x0_initial; 

%Save output for country-specific calibration
Y_NGA = Y; 

%Save Nigerian calibration
save cal_NGA; 

%% Re-calibrate for each country
%Recalibration parameters: chi, Ag, gamma, Es/E
clear
load cal_NGA
countries ={'GHA','TZA', 'UGA'};
%M = xlsread('C:\Users\sdfried\Dropbox (ASU)\Research\David\IGC_energy\replication_package\data\derived\cross_country_moments.xlsx');
M = xlsread(strcat(reppath, '/data/derived/cross_country_moments.xlsx'));
frac_no_outages_vec = M(2:end,1);   % Fraction of modern firms that never experience an outage
frac_self_vec =M(2:end,2);          % Fraction of electricity modern firms with generators generate themselves
Nm_gen_Nm_q0_vec =M(2:end,3);       % Fraction of firms that experience outages that have access to a generator
ps_pg_vec = M(2:end,4);             % Ratio of marginal costs: self-generated electricity to grid
ac_pg_vec = M(2:end, 5);            % Ratio of average costs: self-generated electricity to grid
rel_y_pop_vec = M(2:end, 6);        % Ratio of output per worker in country relative to Nigeria

%Initial guesses p0 for country specific calibration
p0guess = [0.1248    0.8698    1.5237    0.2493    0.0396;
    0.1080    0.4909    2.2775    1.6075    0.1250;
    0.0998    0.4286    3.0000    0.0773    0.1324];

%Loop over countries
for ci =1:1:length(countries)
       
    %Clear variables from previous calibration
    clearvars -except p0guess  gamma sol Nm_gen_Nm_q0_vec x0 R frac_no_outages_vec A ups Ag v_NGA rel_y_pop_vec Y_NGA ci Pg frac_self_vec options ps_pg_vec ac_pg_vec lambda  mu N alpha beta delta eta phi psi chi Omega Ls q countries nz zmax zeta
    
    % Set flags for country-specific calibration
    indicate = 1;
    decomp=0;
    tau =0; 
    find_elas=0; 

    %Calculate chi directly from moments
    chi = (ac_pg_vec(ci)/ps_pg_vec(ci)*0.8 - 0.8)/R;

    %Vector of targets; 
    targets  =[frac_self_vec(ci), rel_y_pop_vec(ci), ps_pg_vec(ci),  Nm_gen_Nm_q0_vec(ci), frac_no_outages_vec(ci)];
    
    p0 = p0guess(ci, :);
    params =[p0, Omega, chi, mu,....
        ups, eta,phi , alpha, lambda, Ag, Ls, beta, psi, N, delta, nz, zmax, Y_NGA];
    paramsFixed = params(6:end);
    
    [dist1, pStar, sol] = calibrate(p0, paramsFixed, targets, x0, options, indicate);

    %Calculate moments and targets
    p0  = pStar;
    params =[pStar, paramsFixed];
    [sol, fval, exitflag] = fsolve(@errors_ss,sol,options,params, indicate);
    [errorsTry, moments] = errors_ss(sol,params,  indicate);
    [moments; targets]
    
    %Calculate elasticity
    W = sol(1);  v = sol(2);
    find_elas=1;
    solveModel;
    
    %Save intital guess
    x0_initial = sol;
    x0 = x0_initial;
    
    %save country specific calibration
    indicate =0;
    clear indicate  kgclear decomp find_elas
    save(strcat('cal_', countries{ci}));    
end

% Check country-specific calibration
for ci = 1:length(countries)
    load(strcat('cal_', countries{ci})); 
    countries{ci}
    [moments;targets]
end

%% Sensitivity calibration for Table BLANK

%vector of names for sensitivity calibrations
sens_names = {'NGA_hs', 'NGA_ls', 'NGA_bp', 'NGA_sp',  'NGA_lc', 'NGA_hc', 'NGA_q3', 'NGA_q1' ,'NGA_psi6', 'NGA_psi8'};

%initial guess of p0 for each sensitivity calibration
p0guess =[0.1220    0.6170    0.8830    0.1829    0.0773    0.3818    2.1444    1.1630;
    0.0963    2.4966    1.1184    0.3094    0.0814    0.7165    2.1444    1.7934;
    0.1080    3.2125    0.9971    0.3163    0.0676    0.7446    2.1444    1.4814;
    0.1061    0.7174    1.0155    0.1846    0.0889    0.3807    2.1443    1.4845;
    0.1039    2.0130    1.1516    0.2749    0.0799    0.6216    2.1444    1.5042;
    0.1079    1.4448    0.9074    0.2112    0.0719    0.4514    2.1444    1.4372;
    0.1119    1.8308    0.9621    0.2155    0.1101    0.4880    2.1444    1.4661;
    0.0969    1.7147    1.1112    0.2810    0.0391    0.6099    2.1444    1.4979;
    0.0962    2.1893    1.1195    0.2863    0.0791    0.6534    2.1444    1.6177;
    0.1104    1.2684    0.9756    0.1988    0.0728    0.4184    2.1444    1.2641];


for ci =1:length(sens_names)

    clearvars -except p0guess sens_names ci

    %Load baseline calibration for Nigeria
    load cal_NGA

    %Perturbation
    if ci ==1
        elec_share_modern_targ =0.09;
    elseif ci ==2
        elec_share_modern_targ =0.06;
    elseif ci ==3
        phi = 0.5;
    elseif ci ==4
        phi =0.3;
    elseif ci ==5
        ps_pg_targ = 0.9*4.33;
        ac_self_grid_targ = 0.9*5.51;
    elseif ci ==6
        ps_pg_targ = 1.1*4.33;
        ac_self_grid_targ = 1.1*5.51;
    elseif ci ==7
        Nm_q1_Nm_targ = 0.3;
    elseif ci ==8
        Nm_q1_Nm_targ = 0.1;
    elseif ci ==9
        psi =0.6;
    elseif ci ==10
        psi =0.8;
    end

    targets  =[ps_pg_targ, frac_self_targ, share_modern_labor_targ, ...
        frac_modern_firm_targ,  ac_self_grid_targ, frac_modern_generators_target,...
        Nm_q1_Nm_targ, elec_share_modern_targ];

    %Initial guess for internally determined parameter values
    p0 = p0guess(ci, :);

    %Vector of parameter values
    params =[p0, ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];
    paramsFixed = params(9:end);

    %Calibrate model
    [dist1, pStar, sol] = calibrate(p0, paramsFixed, targets, x0, options, indicate);

    % Report calibrated moments and targets
    p0  = pStar;
    params =[pStar, paramsFixed];
    [sol, fval, exitflag] = fsolve(@errors_ss,sol,options,params, indicate);
    [errorsTry, moments1, guess1] = errors_ss(sol,params,  indicate);
    [xstar1, fval, exitflag] = fsolve(@errors_sspg, guess1,options,params, indicate, .928);
    [errorsTry, moments2] = errors_sspg(xstar1,params, indicate, 0.928);
    moments = [moments1(1:7), moments2];
    fprintf('Moments and targets')
    [moments; targets]

    %Find elasticity of generator ownership with respect to firm size
    W = sol(1); v  = sol(2);
    find_elas =1;
    solveModel;

    save(strcat('cal_', sens_names{ci}));

end

% Check sensitivity calibrations
for cj = 1:length(sens_names)
    load(strcat('cal_', sens_names{cj})); 
    sens_names{cj}
    [moments;targets]
end

%% Generator fixed cost calibrations

%Names for each calibration of upsilon
ups_names = {'NGA_ups05', 'NGA_ups10', 'NGA_ups15', 'NGA_ups20'};

%Different values of upsilon
upsvec =[0.05, 0.1, 0.15, 0.2];

%Initial guesses for calibration 
p0guess = [  0.1048    1.7871    1.0275    0.2110    0.0697    0.4498    2.1444    1.4819;
    0.1042    1.8492    1.0339    0.1775    0.0648    0.3588    2.1444    1.4852;
    0.1037    1.9239    1.0386    0.1424    0.0614    0.2627    2.1444    1.4861;
    0.1034    2.0242    1.0412    0.1053    0.0597    0.1609    2.1444    1.4849];

for upsi =1:length(ups_names)

    clearvars -except p0guess ups_names upsi upsvec

    %Load baseline calibration for Nigeria
    load cal_NGA

    ups = upsvec(upsi);
    %Initial guess for internally determined parameter values
    p0 = p0guess(upsi, :);

    %Vector of parameter values
    params =[p0, ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];
    paramsFixed = params(9:end);

     %Calibrate model
    [dist1, pStar, sol] = calibrate(p0, paramsFixed, targets, x0, options, indicate);

    % Report calibrated moments and targets
    p0  = pStar;
    params =[pStar, paramsFixed];
    [sol, fval, exitflag] = fsolve(@errors_ss,sol,options,params, indicate);
    [errorsTry, moments1, guess1] = errors_ss(sol,params,  indicate);
    [xstar1, fval, exitflag] = fsolve(@errors_sspg, guess1,options,params, indicate, .928);
    [errorsTry, moments2] = errors_sspg(xstar1,params, indicate, 0.928);
    moments = [moments1(1:7), moments2];
    fprintf('Moments and targets')
    [moments; targets]

    %Find elasticity of generator ownership with respect to firm size
    W = sol(1); v  = sol(2);
    find_elas =1;
    solveModel;

    save(strcat('cal_', ups_names{upsi}));
end

% Check upsilon calibrations
for cj = 1:length(ups_names)
    load(strcat('cal_', ups_names{cj})); 
    ups_names{cj}
    [moments;targets]
    ups/Omega
end


%% Perturbation Matrix (For Table C.2)
clear; 
load cal_NGA
perturb_mat = zeros(8, 8);

param_vec = {'chi', 'As', 'Pg', 'lambda', 'Omega', 'zeta', 'gamma', 'mu'};

for i = 1:length(param_vec)
    load cal_NGA

    %Perturb parameters 1 at a time by 1%
    eval(strcat(param_vec{i},'= 1.01*',param_vec{i},';'));
    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];

    [sol, fval, exitflag] = fsolve(@errors_ss,sol,options,params, indicate);
    [errorsTry, moments1, guess1] = errors_ss(sol,params,  indicate);
    [xstar1, fval, exitflag] = fsolve(@errors_sspg, guess1,options,params, indicate, .928);
    [errorsTry, moments2] = errors_sspg(xstar1,params, indicate, 0.928);
    moments = [moments1(1:7), moments2];
    %Percent difference between moments and targets
    perturb_mat(i, :) = (moments - targets)./targets*100;
end


% Sort
order = [1, 3, 4, 5, 2, 6, 7, 8];
temp = [order; perturb_mat];
temp_sorted = sortrows(temp',1);
temp_sorted2 =temp_sorted(:, 2:end)';

save('perturbations', 'temp_sorted2', 'param_vec');




 
 
 


