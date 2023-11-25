%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: nelder.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Forms a new point for a reflection, contraction, expansion in
% nelder-meade loop and calculates the squared distance between the moments and targets 
% at that new point. Adapted to Matlab from Numerical Recipes routine amotry (Press et al. 1988) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ytry, x0]= nelder2(ptry, x0, paramsFixed, targets, options, indicate)



%evaluate distance between moments and targets at trial point
params = [ptry, paramsFixed];

[x,fval,EXITFLAG] = fsolve(@errors_ss,x0,options,params, indicate);

error = 0; 
if EXITFLAG<1 || isreal(x) ==0 || max(abs(fval))>1e-13
    display(['Model did not solve.'])
    beep;
    error = 1; 
    
else 
 %display(['Model did solve. x0 updated'])
  x0 =x;  
  
end

[errorsTry, moments1, guess2] = errors_ss(x0,params, indicate);

temp = errors_sspg(guess2, params, indicate, 0.928);

if max(isnan(temp) ==1)
    guess2(2) =0.5;
end

if indicate ==0
    [xstar, fval, exitflag] = fsolve(@errors_sspg, guess2,options,params, indicate, .928);
    [errorsTry, moments2] = errors_sspg(xstar,params, indicate, 0.928);
    moments = [moments1(1:7), moments2];
else
    moments = moments1;
end

ytry = (targets - moments)*(targets -moments)'; %distance between targets and moments

%if ytry is better than the highest point, replace the hightest point
%with ytry

%in the case of the expansion, replace reflected point with the
%expanded point if the expanded point is better.

end 