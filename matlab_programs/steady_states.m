%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: steady_states.m
% By:  Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Calculates initial and no-outages steady state, decomposes
% gains in output per worker in the SRPE effect, firm entry channel, and
% the firm expansion channel, computes the weak links and tax experiements
% and the proxmiate causes decomposition. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear all;

options = optimset('Algorithm','Levenberg-Marquardt','TolFun',10e-15,'TolX',10e-15, ...
    'MaxFunEvals', 10000, 'MaxIter', 10000, 'Display','none');

%List of countries and robustness exercieses for which to calculate the
%SRPE and LRGE effects of eliminating outages
countries1 = {'NGA', 'GHA', 'TZA', 'UGA', 'NGA_hs', 'NGA_ls', 'NGA_bp', 'NGA_sp',  'NGA_lc', 'NGA_hc', ...
    'NGA_q3', 'NGA_q1' ,'NGA_psi6', 'NGA_psi8', 'NGA_ups05', 'NGA_ups10', 'NGA_ups15', 'NGA_ups20'};

%Variable to save for each steady state
vars ={'x0', 'xstar', 'v', 'W', 'R', 'Lm', 'Kg', 'Y', 'Km', 'K', 'I','Lt' 'Ym', 'Yt', 'A', 'Kt', 'Ls', 'r','Y', 'C', ...
    'Nm', 'Nt', 'Ks','N', 'Es', 'Eg', 'E', 'Pg', 'Ag', 'chi', 'psi', 'As', 'ac_self_grid', 'Ps', 'fval', 'Leontief', ...
    'A', 'Ag', 'K', 'zeta', 'ups', 'Omega', 'gamma', 'frac_self', 'elec_share_modern', 'phi', 'mu','Ys', 'Ym_on', ...
    'prob_m_gen_q0_z', 'vbig',  'prob_m_ngen_q0_z', 'prob_m_q1_z', 'elasticity', 'TFP'};


%% Compare initital steady state with steady state with zero outages

for count = 1:length(countries1)

    %Initial steady state
    clearvars -except countries1 count options vars
    country1 = countries1{count}
    load(strcat('cal_', country1));

    decomp =0;
    indicate =0;
    x0=sol;

    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];


    [xstar, fval, exitflag] = fsolve(@errors_ss,x0,options,params, indicate);


    if max(abs(fval))< 1e-12 && isreal(fval)
        x0 = xstar;
    else
        fval
        disp('Model did not solve')
    end


    W = xstar(1); v = xstar(2);
    find_elas =1;
    solveModel;

    %Save the firm-size generator-ownership elasticity in the initial
    %steady state
    vars1 = vars;
    vars1{end+1} ='elasticity';

    filename = strcat('ss1_', country1);
    for i = 1:length(vars1)
        str = strcat(vars1{i}, '_', country1, '1', '=', vars1{i}, ';');
        eval(str);
        str2 =  strcat(vars1{i}, '_', country1, '1');

        if i ==1 save(filename, str2);
        else
            save(filename, str2, '-append');
        end
    end
    save(strcat('ss1_all_', country1));

    %%%  No-outages steady state
    v = 1;
    x0=[W, Pg];
    [xstar, fval, exitflag] = fsolve(@errors_sspg,x0,options,params, indicate, v);

    if max(abs(fval))< 1e-12 && isreal(fval)
        x0 = xstar;
    else
        fval
        disp('SS2 Model did not solve')
    end

    W = xstar(1); Pg= xstar(2);
    params(1) = Pg;
    solveModel;

    filename = strcat('ss2_', country1);
    for i = 1:length(vars)
        str = strcat(vars{i}, '_', country1, '2', '=', vars{i}, ';');
        eval(str);
        str2 =  strcat(vars{i}, '_', country1, '2');

        if i ==1 save(filename, str2);
        else
            save(filename, str2, '-append');
        end
    end

    save(strcat('ss2_all_', country1));

end

%% Short-run partial equilibrium effects of eliminating power outages
for count = 1:length(countries1)
    country1 = countries1{count};
    load(strcat('ss1_all_', country1));
    Y_SR = Ym_on + Yt;
    SR = (Y_SR/Y-1)*100; %Percent increase in output per worker from producing on the power-on line the whole period
    vars2 = {'SR'};

    filename = strcat('partial_eq');
    for i = 1:length(vars2)
        str = strcat(vars2{i}, '_', country1, '=', vars2{i}, ';');
        eval(str);
        str2 =  strcat(vars2{i}, '_', country1);
        if (i ==1 && count ==1)
            save(filename, str2);
        else
            save(filename, str2, '-append');
        end
    end
end

%% Decomposition: Contribution of each general equilibrium channel
load ss1_all_NGA;

for count = 1:length(countries1)

    clearvars -except options count countries1
    load partial_eq;
    country1 = countries1{count};
    load(strcat('cal_', country1));
    load(strcat('ss1_', country1));
    load(strcat('ss2', '_', country1));
    indicate =0;
    decomp =1;
    find_elas =0;

    %Add Firm Expansion: Only change power outages, sector probabilities
    %are fixed at initial level.
    v=1;
    prob_m_gen_q0_z = eval(strcat('prob_m_gen_q0_z_', country1, '1'));
    prob_m_ngen_q0_z = eval(strcat('prob_m_ngen_q0_z_', country1, '1'));
    prob_m_q1_z = eval(strcat('prob_m_q1_z_', country1, '1'));

    %Initial guess
    W = eval(strcat('W_', country1, '1'));
    Pg = eval(strcat('Pg_', country1, '1'));

    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];

    guess = [W, Pg];
    [xstar, fval, exitflag] = fsolve(@errors_noentry, guess ,options,params, indicate,...
        prob_m_gen_q0_z, prob_m_ngen_q0_z, prob_m_q1_z,  v);
    if max(abs(fval))> 1e-12 || isreal(fval) ==0
        fval
        disp('Model did not solve')
        count
    end
    W = xstar(1);
    Pg = xstar(2);
    params(1) = Pg;
    solveModel;

    SR = eval(strcat('SR_', country1));
    firm_expand = (Y/eval(strcat('Y_', country1, '1')) -1)*100 - SR;

    total = (eval(strcat('Y_', country1, '2'))/eval(strcat('Y_', country1, '1')) -1)*100;

    firm_entry = total - firm_expand - SR;

    vars2 = {'firm_entry', 'firm_expand', 'total'};

    filename = strcat('decomp');
    for i = 1:length(vars2)
        str = strcat(vars2{i}, '_', country1,  '=', vars2{i}, ';');
        eval(str);
        str2 =  strcat(vars2{i}, '_', country1);
        if (i==1 && count ==1)
            save(filename, str2);
        else
            save(filename, str2, '-append');
        end

    end
end

%% Size of the LR GE Effect compared to the SRPE Effect
clear;
load decomp;
load partial_eq;

PE = [SR_GHA, SR_NGA, SR_TZA, SR_UGA];
GE = [total_GHA, total_NGA, total_TZA, total_UGA] - PE;
total = [total_GHA, total_NGA, total_TZA, total_UGA];

mean(total./PE)

%% Weak Links and tax experiments (Appendix E)
clear;
load ss1_all_NGA;
load ss1_NGA; load ss2_NGA;
find_elas =0;

%Outages economy: lower Pg until we get to the initial steady state.
npg = 200; pgmax = Pg_NGA2; pgmin = Pg_NGA1; pginc = (pgmax - pgmin)/(npg-1);
pgvec = [pgmin:pginc:pgmax]';

x0 = x0_NGA1;
Y1vec = zeros(size(pgvec));
Y1vec(1) = Y_NGA1; Y1vec(end) = Y_NGA2;
E1vec = zeros(size(pgvec));
E1vec(1) = E_NGA1; E1vec(end) = E_NGA2;
Egvec = zeros(size(pgvec));
Egvec(1) = Eg_NGA1; Egvec(end) = Eg_NGA2;
Kg1vec = zeros(size(pgvec));
Kg1vec(1) = Kg_NGA1; Kg1vec(end) = Kg_NGA2;
Ks1vec = zeros(size(pgvec));
Ks1vec(1) = Ks_NGA1; Ks1vec(end) = Ks_NGA2;
Km1vec = zeros(size(pgvec));
Km1vec(1) = Km_NGA1; Km1vec(end) = Km_NGA2;
leovec = zeros(size(pgvec));
leovec(1) = Leontief_NGA1; leovec(end) = Leontief_NGA2;
vvec = zeros(size(pgvec));
vvec(1) = v_NGA1; vvec(end) = v_NGA2;
Nmvec = zeros(size(pgvec));
Nmvec(1) = Nm_NGA1; Nmvec(end) = Nm_NGA2;

cvec = zeros(size(pgvec));
cvec(1) = C_NGA1; cvec(end) = C_NGA2;

for i = 2:1:length(pgvec)-1

    Pg = pgvec(i);

    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];

    [xstar, fval, exitflag] = fsolve(@errors_ss,x0,options,params, indicate);

    if max(abs(fval))< 1e-12 && isreal(fval)
        x0 = xstar;
    else
        fval
        disp('Model did not solve')
    end

    W = xstar(1); v  =xstar(2);
    solveModel;

    Y1vec(i) = Y;
    E1vec(i) = E;
    Egvec(i) = Eg;
    Kg1vec(i) = Kg;
    Ks1vec(i) = Ks;
    Km1vec(i) = Km;
    leovec(i) = Leontief;
    vvec(i) = v;
    Nmvec(i) = Nm;
    wvec(i) = W;
    cvec(i) = C;
end

%Remove duplicate regimes
ind1 = find(leovec > 1.0e-15, 1, 'first')-1;
ind2 = find (Y1vec(ind1+1:end) >= Y1vec(ind1), 1, 'first') +ind1;


%Competitive economy: decrease Ag by E_NGA1/E_NGA2, no outages
agvec = zeros(size(pgvec));
agvec(npg) = 1;
for i = 1:npg-1
    agvec(i) = 1*E1vec(i)/E_NGA2;
end

x0 = x0_NGA2;
v = 1;
Y2vec = zeros(size(pgvec));
Y2vec(end) = Y_NGA2;
E2vec = zeros(size(pgvec));
E2vec(1) = E_NGA1; E2vec(end) = E_NGA2;
Kg2vec = zeros(size(pgvec));
Kg2vec(1) = Kg_NGA2; Kg2vec(end) = Kg_NGA2;

for i = 1:1:length(Y2vec)
    i
    Ag = agvec(i);
    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];

    [xstar, fval, exitflag] = fsolve(@errors_sspg,x0,options,params, indicate, v);

    if max(abs(fval))< 1e-12 && isreal(fval)
        x0 = xstar;
    else
        fval
        disp('SS2 Model did not solve')
    end

    W = xstar(1); Pg  =xstar(2);
    params(1) = Pg;
    solveModel;

    Y2vec(i) = Y;
    E2vec(i) = E;
    Kg2vec(i) = Kg;
end

% Tax experiment
tau_vec = pgmax - pgvec;

Y3vec = zeros(size(pgvec));
Y3vec(end) = Y_NGA2;
E3vec = zeros(size(pgvec));
E3vec(1) = E_NGA1; E3vec(end) = E_NGA2;
Kg3vec = zeros(size(pgvec));
Kg3vec(1) = Kg_NGA2; Kg3vec(end) = Kg_NGA2;
Ag = Ag_NGA1;
x0 = x0_NGA2;
v =1;

for i = length(Y3vec):-1:1
    i
    tau = tau_vec(i);
    params =[Pg, lambda , As, zeta, gamma, Omega, chi, mu,....
        ups, eta,phi , alpha, A, Ag, Ls, beta, psi, N, delta, nz, zmax];

    [xstar, fval, exitflag] = fsolve(@errors_sspg_tax,x0,options,params, indicate, v, tau);

    if max(abs(fval))< 1e-12 && isreal(fval)
        x0 = xstar;
    else
        fval
        disp('SS3 Model did not solve')
    end

    W = xstar(1); Pg  =xstar(2);
    params(1) = Pg;
    solveModel;

    Y3vec(i) = Y;
    E3vec(i) = E;
    Kg3vec(i) = Kg;
end

save weak_links_plus_tax

%% Capital versus productivity decomposition (Appendix E)

%Subset of countries for second decomposition
countries2= {'GHA', 'NGA', 'TZA', 'UGA'};

PE_country = zeros(4, 1);
Y_SR_K =zeros(4,1);
Y_SR_prod = zeros(4,1);
GE_country = zeros(4,1);

for j1 = 1:length(countries2)
    country1 = countries2{j1}
    load(strcat('ss1_all_', country1));
    load(strcat('ss2_', country1));
    Km2 = eval(strcat('Km_', country1, '2'));
    Kt2 = eval(strcat('Kt_', country1, '2'));
    Km1 = eval(strcat('Km_', country1, '1'));
    Kt1 = eval(strcat('Kt_', country1, '1'));
    Yt2 = eval(strcat('Yt_', country1, '2'));
    Yt1 = eval(strcat('Yt_', country1, '1'));
    Y1 = eval(strcat('Y_', country1, '1'));
    Y2 = eval(strcat('Y_', country1, '2'));
    deltaK_prod = (Km2+ Kt2)/(Km1+ Kt1);
    Y_SR = Ym_on + Yt;
    PE_country(j1) = (Y_SR/Y1-1)*100; %Percent increase in output per worker from producing on the power-on line the whole period

    decomp =2;
    deltaK = deltaK_prod;
    solveModel;
    Y_SR_K(j1) = ((Ym_on + Yt)/Y1 -1)*100;

    %Reallocate (raise productivity in the traditional sector)
    Y_SR_prod(j1) = ((Ym_on + 1.4*Yt*(1-Yt2/Yt1) +Yt*(Yt2/Yt1))/Y1 -1)*100 ;

    %Long run GE effect
    GE_country(j1) = (Y2/Y1-1)*100;
end

save('decomp_2', 'PE_country', 'Y_SR_K', 'Y_SR_prod', 'GE_country');

%% Check resource constraints and FOCS
clc;
for cj = 1:length(countries1)
    country1 = countries1{cj};
    load(strcat('ss1_all_', country1));
    if ( max(abs([resource, foc_km, foc_lm, foc_ks, foc_kg])) > 1.0e-14)
        fprintf('Problem in %s \n', country1);
    end
end

for cj = 1:length(countries1)
    country1 = countries1{cj};
    load(strcat('ss2_all_', country1));
    if ( max(abs([resource, foc_km, foc_lm, foc_kg])) > 1.0e-14)
        fprintf('Problem in %s \n', country1);
    end
end




