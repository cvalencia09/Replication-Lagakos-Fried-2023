%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: errors_ss.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Reports the value of the equilibrium and market clearing
% conditions and the model values of the relevant moments. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errors_vec, moments_vec, guess2] = errors_ss(guess,params, indicate) 

decomp=0; 
tau =0; 
find_elas =0; 

W= guess(1);
v = guess(2);

solveModel;

M = 1e5;
w_penalty = (W< 0);
v_penalty = ((v< 0) || v>1);
penalty = M*(w_penalty + v_penalty);

errors_vec = zeros(1, 2);

%Minimization vector
errors_vec(1)= labor_market + penalty;
errors_vec(2) = elec_market + penalty;

moments_vec =[ps_pg, frac_self, share_modern_labor, frac_modern_firms, ac_self_grid, Nm_gen_Nm_q0, ...
    Nm_q1_Nm, elec_share_modern];

if indicate ==0
    moments_vec = moments_vec;
elseif indicate ==1
    Y_NGA = params(22); 
   moments_vec = [frac_self, Y/Y_NGA, ps_pg,  Nm_gen_Nm_q0, Nm_q1_Nm];
end

guess2 = [W, Pg];

