%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: solveModel.m
% By:  Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Computes equilibrium conditions aggregate variables, moments,
% and elasticity of generator ownership with respect for firm size, given
% a value of the wage and the probability of grid power. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Indicate = 1 during the country-sepcific calibration
if (indicate ==0)
    lambda = params(2);
    A = params(13);
    
elseif (indicate ==1)
    lambda = params(13);
    A = params(2);
end

Pg = params(1); As = params(3);  zeta = params(4);  gamma =params(5);
Omega = params(6);
chi = params(7);  mu = params(8); ups = params(9);
eta = params(10);  phi = params(11);  alpha = params(12); Ag = params(14);

Ls = params(15); beta = params(16); psi= params(17);
N = params(18); nz =params(20); zmax = params(21); delta = params(19);

%ups = ups_omega*Omega; 

%Grid over modern values of z
zgrid = [1:(zmax - 1)/(nz-1): zmax];
Am = (1+phi)*A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Prices and constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 1/beta - 1+delta; %rental rate
r = R - delta;
cons1 = (1-alpha)/alpha*R/W;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Grid electricity capital and supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Grid capital
Kg = (psi*(Pg-tau)*Ag/R)^(1/(1-psi));

%Grid electricity supply
Eg_supply = Ag*Kg^psi;

%Grid electricity profits
Pig = (Pg-tau)*Eg_supply - R*Kg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Traditional sector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Capital and labor demands per unit of productivity
Kt_til = (eta*alpha*A*cons1^(eta*(1-alpha))/R)^(1/(1-eta));
Lt_til = cons1*Kt_til;
if (decomp ==2)
    Kt_til = deltaK*Kt_til;
end

%Capital per unit of productivity
Yt_til = A*(Kt_til^alpha*Lt_til^(1-alpha))^eta; 


%Capital, labor, output and profits for entrepreneurs with each value of z
Kt_z = Kt_til.*zgrid; 
Lt_z = Lt_til.*zgrid; 
Yt_z = Yt_til.*zgrid; 
Pit_z = (1-eta)*Yt_z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Modern sector: q=0 (doesn't always get electricity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Firm generates her own electricity 

%Capital, labor, and electricity demands per unit of productivity if firm
%idles capital when the power is out
Ks_til_gen = (eta*Am *(mu*As)^eta/ ( R/(1-v) + 1/chi) )^(1/(1-eta));
Km_til_gen= ( v*alpha*eta*cons1^( eta*(1-alpha))*Am/( R + v*alpha*cons1^(1-alpha)*Pg*(1/mu)))^(1/(1-eta));
Lm_til_gen = cons1*Km_til_gen;
if (decomp ==2)
    Km_til_gen = deltaK*Km_til_gen; 
    Ks_til_gen = deltaK*Ks_til_gen; 
end
Xm_til_gen = Km_til_gen^alpha*Lm_til_gen^(1-alpha);
Eg_til_gen = v*Xm_til_gen/mu;
Es_til_gen = As*Ks_til_gen*(1-v);

%Check if the firm idles X when the power is out.
Leontief = Xm_til_gen - mu*As*Ks_til_gen; % Greater than zero if the producer idles capital and labor when the power is out.

if (Leontief >0)
    vbig =1; %Indicates v >vstar, so firm idles capital and labor when the power is out
else
    vbig =0; 
end

if Leontief <= 0  % Firm does not idle capital or labor when the power is out (impose that AsKs = X)
    %Capital, labor, and electricity  demands per unit of productivity
    Km_til_gen = ( eta*Am *cons1^((1-alpha)*(eta-1))/( R*( 1/(alpha*cons1^(1-alpha)) + 1/(mu *As)) ...
        + v *Pg*(1/mu) + (1-v)*(1/(chi*mu*As))))^(1/(1-eta));
    Lm_til_gen = cons1*Km_til_gen;

    if (decomp ==2)
        Km_til_gen = deltaK*Km_til_gen;
    end

    Xm_til_gen = Km_til_gen^alpha*Lm_til_gen^(1-alpha);
    Ks_til_gen = Xm_til_gen/(mu*As);

    if (decomp ==2)
        Ks_til_gen = deltaK*Ks_til_gen;
    end

    Es_til_gen = As*Ks_til_gen*(1-v);
    Eg_til_gen = v*Xm_til_gen/mu;
    Leontief = Xm_til_gen - mu*As*Ks_til_gen; %Should equal zero
end

%Self-generated electricity and modern output per unit of productivity
Ys_til_gen = (1-v)*(1/chi)*Ks_til_gen;
Ym_til_gen =Am*(v*Xm_til_gen^eta +(1-v)*(mu*As*Ks_til_gen)^eta);
Ym_til_gen_on = Am*Xm_til_gen^eta; 

%Profits, capital, labor, output, electricity, for entrepreneurs with each
%value of z 
Pim_gen_z = (1-eta)*Ym_til_gen.*zgrid- A*ups - A*Omega; 
Km_gen_z = Km_til_gen.*zgrid; 
Lm_gen_z = Lm_til_gen.*zgrid; 
Ym_gen_z = Ym_til_gen.*zgrid; 
Eg_gen_z = Eg_til_gen.*zgrid; 
Es_gen_z = Es_til_gen.*zgrid; 
Ks_gen_z = Ks_til_gen.*zgrid;

%Modern output for entrepreneurs with each value of z when the power is on
%ex-post
Ym_gen_on_z = Ym_til_gen_on.*zgrid; 

%%%% Firm does not generate its own electricity

%Capital, labor, and electricity demands per unit of productivity
Km_til_ngen =( v*alpha*eta*cons1^( eta*(1-alpha))*Am/( R + v*alpha*cons1^(1-alpha)*Pg*(1/mu)))^(1/(1-eta));

%%
Lm_til_ngen = cons1*Km_til_ngen;
if (decomp ==2)
    Km_til_ngen = deltaK*Km_til_ngen;
end
Xm_til_ngen = Km_til_ngen^alpha*Lm_til_ngen^(1-alpha);
Eg_til_ngen = v*Xm_til_ngen/mu;

%Modern output and profits per unit of productivity
Ym_til_ngen = Am*v*Xm_til_ngen^eta;
Pim_til_ngen = (1-eta)*Ym_til_ngen - A*Omega; 

%Modern output when the power is on ex-post
Ym_til_ngen_on = Am*Xm_til_ngen^eta; 

%Profits, capital, labor, output, electricity, for entrepreneurs with each
%value of z 
Pim_ngen_z = (1-eta)*Ym_til_ngen.*zgrid -A*Omega; 
Km_ngen_z =  Km_til_ngen.*zgrid; 
Lm_ngen_z = Lm_til_ngen.*zgrid; 
Ym_ngen_z = Ym_til_ngen.*zgrid; 
Ym_ngen_on_z = Ym_til_ngen_on.*zgrid; 
Eg_ngen_z = Eg_til_ngen.*zgrid; 
Es_ngen_z = 0.*zgrid; 
Ks_ngen_z = 0.*zgrid; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Modern sector: q=1 (always gets grid electricity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Capital, labor, and electricity demands per unit of productivity 
Km_til_q1 = (1*alpha*eta*cons1^(eta*(1-alpha))*Am/( R + 1*alpha*cons1^(1-alpha)*Pg*(1/mu)))^(1/(1-eta));
Lm_til_q1 = cons1*Km_til_q1;
if (decomp ==2)
    Km_til_q1 = deltaK*Km_til_q1;
end
Xm_til_q1 = Km_til_q1^alpha*Lm_til_q1^(1-alpha);
Ks_til_q1 = 0;
Es_til_q1 = 0;
Eg_til_q1 = 1*Xm_til_q1/mu;

%Self-generated electricity and modern output per unit of productivity
Ys_til_q1 = 0;
Ym_til_q1 =Am*(1*Xm_til_q1^eta +(1-1)*(mu*As*Ks_til_q1)^eta);

%Modern output when the power is on ex-post
Ym_til_q1_on = Am*Xm_til_q1^eta; 

%Profits, capital, labor, output, electricity, for entrepreneurs with each
%value of z 
Pim_til_q1 = (1-eta)*Ym_til_q1 - A*Omega; %Profits exclusive of the gumbel taste shock
Pim_q1_z = (1-eta)*Ym_til_q1.*zgrid- A*Omega; 
Km_q1_z = Km_til_q1.*zgrid; 
Lm_q1_z = Lm_til_q1.*zgrid; 
Ym_q1_z = Ym_til_q1.*zgrid; 
Eg_q1_z = Eg_til_q1.*zgrid; 
Es_q1_z = 0.*zgrid; 
Ks_q1_z = 0.*zgrid;
Ym_on_q1_z = Ym_til_q1_on.*zgrid; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Aggregate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (decomp ==0)  %if decomp is greater than zero, then the distribution of firms across generators is exogenous
   
    %Probability firm with q=0 is in modern sector and has a generator
    prob_m_gen_q0_z = exp(zeta^(-1)*Pim_gen_z)./...
        (exp(zeta^(-1)*Pim_gen_z)+ exp(zeta^(-1)*Pim_ngen_z)+exp(zeta^(-1)*Pit_z));
    
    %Probability firm with q=0 is in modern sector and does not have a generator
    prob_m_ngen_q0_z = exp(zeta^(-1)*Pim_ngen_z)./(exp(zeta^(-1)*Pim_gen_z)+exp(zeta^(-1)*Pim_ngen_z)+exp(zeta^(-1)*Pit_z));
   
    %Probability firm with q=1 is in modern sector
    prob_m_q1_z = (exp(zeta^(-1)*Pim_q1_z) + exp(zeta^(-1)*Pim_q1_z))./...
        (exp(zeta^(-1)*Pim_q1_z)+exp(zeta^(-1)*Pim_q1_z)+exp(zeta^(-1)*Pit_z));
end

%Probability of each value of z (comes from the pdf of the pareto distribution)
sn =0.0001; %Find the probability of each point on the zgrid as the difference in the cdf for an interval of length 2*sn
prob_z = (1-(zgrid  + sn).^(-lambda)) - (1-(zgrid - sn).^(-lambda)); 

%Traditional capital, labor, and output
Kt = (1-gamma)*prob_z*((1-prob_m_gen_q0_z - prob_m_ngen_q0_z).*Kt_z)'/(sum(prob_z)) + ...
    gamma*prob_z*((1-prob_m_q1_z).*Kt_z)'/(sum(prob_z));
Lt = (1-gamma)*prob_z*((1-prob_m_gen_q0_z - prob_m_ngen_q0_z).*Lt_z)'/(sum(prob_z)) + ...
    gamma*prob_z*((1-prob_m_q1_z).*Lt_z)'/(sum(prob_z));
Yt = (1-gamma)* prob_z*((1-prob_m_gen_q0_z - prob_m_ngen_q0_z).*Yt_z)'/(sum(prob_z))+ ...
    gamma* prob_z*((1-prob_m_q1_z).*Yt_z)'/(sum(prob_z));

%Modern capital
Km_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Km_gen_z)'/(sum(prob_z)); 
Km_ngen_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Km_ngen_z)'/(sum(prob_z));
Km_q1=  gamma*prob_z*(prob_m_q1_z.*Km_q1_z)'/(sum(prob_z)); 
Km = Km_gen_q0 + Km_ngen_q0 + Km_q1;
Km_q0 = Km_gen_q0 + Km_ngen_q0; 

%Modern labor
Lm_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Lm_gen_z)'/(sum(prob_z)); 
Lm_ngen_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Lm_ngen_z)'/(sum(prob_z));
Lm_q1= gamma*prob_z*(prob_m_q1_z.*Lm_q1_z)'/(sum(prob_z));  
Lm = Lm_gen_q0 + Lm_ngen_q0 + Lm_q1;

%Modern output
Ym_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Ym_gen_z)'/(sum(prob_z)); 
Ym_ngen_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Ym_ngen_z)'/(sum(prob_z)); 
Ym_q1= gamma*prob_z*(prob_m_q1_z.*Ym_q1_z)'/(sum(prob_z));  
Ym = Ym_gen_q0 + Ym_ngen_q0 + Ym_q1;

%Modern output when the power is on ex-post
Ym_gen_on_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Ym_gen_on_z)'/(sum(prob_z)); 
Ym_ngen_on_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Ym_ngen_on_z)'/(sum(prob_z)); 
Ym_on_q1= gamma*prob_z*(prob_m_q1_z.*Ym_on_q1_z)'/(sum(prob_z));  
Ym_on = Ym_gen_on_q0 + Ym_ngen_on_q0 + Ym_on_q1;

%Demand for grid electricity
Eg_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Eg_gen_z)'/(sum(prob_z)); 
Eg_ngen_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Eg_ngen_z)'/(sum(prob_z)); 
Eg_q1= gamma*prob_z*(prob_m_q1_z.*Eg_q1_z)'/(sum(prob_z));  
Eg = Eg_gen_q0 + Eg_ngen_q0 + Eg_q1;
Eg_demand_rationed = Eg; 

%Self generated electricity
Es_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Es_gen_z)'/(sum(prob_z)); 
Es_ngen_q0 = (1-gamma)*prob_z*(prob_m_ngen_q0_z.*Es_ngen_z)'/(sum(prob_z));
Es_q1 =0; 
Es = Es_gen_q0 + Es_ngen_q0 + Es_q1;

%Generator capital
Ks_gen_q0 = (1-gamma)*prob_z*(prob_m_gen_q0_z.*Ks_gen_z)'/(sum(prob_z)); 
Ks = Ks_gen_q0;

%Numbers of  firms
Nm= (1-gamma)*prob_z*(prob_m_gen_q0_z + prob_m_ngen_q0_z)'/(sum(prob_z)) + ...
    gamma*prob_z*(prob_m_q1_z)'/(sum(prob_z));
Nm_q1 =  gamma*prob_z*(prob_m_q1_z)'/(sum(prob_z));
Nm_gen= (1-gamma)*prob_z*(prob_m_gen_q0_z)'/(sum(prob_z));
Nm_ngen = Nm - Nm_gen;
Nt = 1 -Nm;
Nm_q0 = Nm - Nm_q1; 
Nm_ngen_q0 = Nm_ngen - Nm_q1;

%Aggregate output, electricity and capital
E = Es + Eg;
K = Km + Kg + Kt + Ks;
Ys = (1-v)*(1/chi)*Ks;
Y = Yt + Ym - Ys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Market clearing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labor_market = Lt + Lm -Ls; 
elec_market = Eg_supply - Eg_demand_rationed;
Eg = Eg_supply;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Moments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable cost ratio
Ps = 1/(chi*As); %Variable cost of self generated energy
ps_pg = Ps/Pg;

%Average cost ratio
Ys = (1-v)*(1/chi)*Ks;
ac_grid = R*Kg/Eg;
ac_self_cf_8 = (R*Ks + 0.8*(1/chi)*Ks)/(0.8*As*Ks); % average cost of self with capacity factor = 0.8, consistent with data
ac_self_grid =ac_self_cf_8/Pg;

%Electricity share of modern output
elec_share_modern = (Pg*Eg + Ys)/Ym; 

%Fraction of modern labor and firms
share_modern_labor = Lm/Ls;
frac_modern_firms = Nm/N;

%Fraction of electricity frims with generaters produce themselves
frac_self = Es_gen_q0/(Es_gen_q0 + Eg_gen_q0); 

%Fraction of modern firms that experience outages (Nm- Nm_q1) that have
%generators
Nm_gen_Nm_q0 = Nm_gen/(Nm - Nm_q1); 

%Fraction of modern firms that never exeperience outages
Nm_q1_Nm = Nm_q1/Nm; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Measured TFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ktilde = K + Nm*A*Omega + Nm_gen*A*ups; 
TFP = Y/Ktilde^alpha; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Aggregate Resource Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = delta*(Kt + Km + Ks + Kg);
Pit = (1-eta)*Yt; 
Pim = (1-eta)*Ym;
C = W*(Lt + Lm) + R*(Km + Kt + Ks+ Kg) - I + Pit +  Pim + Pig- A*Nm*Omega - Nm_gen*A*ups;
resource = Yt + Ym  - Ys- I - A*Nm*Omega - Nm_gen*A*ups - C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Algebra checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profit_check_t = Yt - W*Lt - R*Kt - Pit;
profit_check_m = Ym - W*Lm - R*Km - R*Ks -Ys - Pg*Eg - Pim;
profit_check_g = (Pg-tau)*Eg_supply - R*Kg  - Pig;

zi =1.05;
Kmi = Km_til_gen*zi; Lmi = Lm_til_gen*zi;

if abs(Leontief)<1.0e-15 %Firm does not idle capital and labor when the power is out (v < vstar)
    foc_km = alpha*(Kmi/Lmi)^(alpha-1)*(eta *Am*zi^(1-eta)*(Kmi^alpha*Lmi^(1-alpha))^(eta-1)...
        - R*(1/(mu*As)) - v*Pg*(1/mu) - (1-v)*(1/(chi*mu*As))) - R;
    foc_lm = (1-alpha)*(Kmi/Lmi)^(alpha)*(eta *Am*zi^(1-eta)*(Kmi^alpha*Lmi^(1-alpha))^(eta-1)...
        - R*(1/(mu*As)) - v*Pg*(1/mu) - (1-v)*(1/(chi*mu*As))) - W;
    foc_ks = 0;
    foc_kg = psi*(Pg-tau)*Ag*Kg^(psi-1) - R;
else %producer idles capital and labor whenever the power is out (v > vstar)
    Ksi = Ks_til_gen*zi;
    foc_km = v*alpha*(Kmi/Lmi)^(alpha-1)*(eta*Am*zi^(1-eta)*(Kmi^alpha*Lmi^(1-alpha))^(eta-1) - Pg*(1/mu)) - R;
    foc_lm = v*(1-alpha)*(Kmi/Lmi)^(alpha)*(eta*Am*zi^(1-eta)*(Kmi^alpha*Lmi^(1-alpha))^(eta-1) - Pg*(1/mu)) - W;
    foc_ks = (1-v)*(eta*Am*zi^(1-eta)*(mu*As *Ksi)^(eta-1)*mu*As- 1/chi)-R;
    foc_kg = psi*(Pg-tau)*Ag*Kg^(psi-1) - R;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Semi-elasticity of generator ownership with respect to firm size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (find_elas ==1)
    avg_firm_size = (Lm + Lt)/(Nm + Nt);
    
    firm_size_vec = [1.05*Lt_til:0.01:Lm_til_gen*5];
    frac_gen =zeros(size(firm_size_vec));
    
    for si = 1:length(firm_size_vec)
        Li = firm_size_vec(si);

        %Find the number of firms with Li workers
        z_gen_Li = Li/Lm_til_gen;
        z_ngen_Li = Li/Lm_til_ngen;
        z_t_Li = Li/Lt_til;
        z_m_q1_Li = Li/Lm_til_q1;
        z_t_q1_Li = Li/Lt_til;
        
        if (z_gen_Li < 1)
            frac_gen(si) =0;
        elseif (z_gen_Li > zmax)
            frac_gen(si) =1;
        else
            zgrid_aux = abs(zgrid - z_gen_Li);
            ind = find((zgrid_aux) == min((zgrid_aux)));
            num_gen_q0 =(1-gamma)*((1-(z_gen_Li  + sn)^(-lambda)) - (1-(z_gen_Li - sn)^(-lambda)))*prob_m_gen_q0_z(ind);
            
            zgrid_aux = abs(zgrid - z_ngen_Li);
            ind = find((zgrid_aux) == min((zgrid_aux)));
            num_ngen_q0 =(1-gamma)*((1-(z_ngen_Li  + sn)^(-lambda)) - (1-(z_ngen_Li - sn)^(-lambda)))*prob_m_ngen_q0_z(ind);
            
            zgrid_aux = abs(zgrid - z_t_Li);
            ind = find(zgrid_aux == min((zgrid_aux)));
            num_t_q0 =(1-gamma)*((1-(z_t_Li  + sn)^(-lambda)) - (1-(z_t_Li - sn)^(-lambda)))*(1-prob_m_ngen_q0_z(ind) - prob_m_gen_q0_z(ind));
            
            zgrid_aux = abs(zgrid - z_m_q1_Li);
            ind = find(zgrid_aux == min((zgrid_aux)));
            num_m_q1 = gamma*((1-(z_m_q1_Li  + sn)^(-lambda)) - (1-(z_m_q1_Li - sn)^(-lambda)))*prob_m_q1_z(ind);
            
            zgrid_aux = abs(zgrid - z_t_q1_Li);
            ind = find(zgrid_aux == min((zgrid_aux)));
            num_t_q1 = gamma*((1-(z_t_q1_Li  + sn)^(-lambda)) - (1-(z_t_q1_Li - sn)^(-lambda)))*(1-prob_m_q1_z(ind));
                        
            frac_gen(si) = num_gen_q0/(num_gen_q0 + num_ngen_q0 + num_t_q0 + num_m_q1 + num_t_q1);
        end
    end
    
      
    Xvar = [ones(size(firm_size_vec))', log(firm_size_vec)'];
    Yvar = frac_gen';
    
    reg_coeffs = (Xvar'*Xvar)^(-1)*Xvar'*Yvar;
    elasticity = reg_coeffs(2); 
    clear si
    
end












