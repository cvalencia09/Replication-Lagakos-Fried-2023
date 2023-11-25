%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: errors_noentry.m
% By: Stephie Fried and David Lagakos
% Date: Winter 2022
% Purpose: calculates the labor market clearing condition and the grid electricity
% market clearing condition for a given value of the wage and price of grid 
% electricity, when the distribution of entrepreneurs is fixed at its value
% in the initial steady state.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errors_vec] = errors_noentry(guess,params, indicate,   prob_m_gen_q0_z,prob_m_ngen_q0_z, prob_m_q1_z, v)
%%
decomp =1; 
tau =0; 
find_elas =0; 

W= guess(1);
Pg = guess(2);

params(1) = Pg; 
solveModel; 

M = 1e5;
w_penalty = (W< 0);
penalty = M*(w_penalty); 

errors_vec = zeros(2,1); 
%Minimization vector
errors_vec(1)= labor_market + penalty;
errors_vec(2) = elec_market + penalty;
