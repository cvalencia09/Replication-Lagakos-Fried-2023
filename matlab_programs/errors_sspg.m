%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: errors_sspg.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: calculates the labor market clearing condition 
% and the grid electricity market clearing condition for a given value of 
% the wage and price of grid electricity.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errors_vec, moments_vec] = errors_sspg(guess,params, indicate, v) 
%%
decomp =0; 
tau =0; 
find_elas =0; 

W= guess(1);
Pg = guess(2);

params(1) = Pg; 

solveModel;

M = 1e5;
w_penalty = (W< 0);
penalty = M*(w_penalty);

%Minimization vector
errors_vec = zeros(1, 2);
errors_vec(1)= labor_market + penalty;
errors_vec(2) = elec_market + penalty;

moments_vec = elec_share_modern;

end

