%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: nelder.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Forms a new point for a reflection, contraction, expansion in
% nelder-meade loop and calculates the squared distance between the moments and targets 
% at that new point. Adapted to Matlab from Numerical Recipes routine amotry (Press et al. 1988) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ytry, x0, yNew, pNew, psumNew,problem]= nelder(p, y, psum, ndim, ihi, fac, x0, paramsFixed,...
     targets, options, indicate)

%form the new point (based on whether it is a reflection, contraction,
%expansion, etc)
fac1 = (1-fac)/ndim;
fac2 = fac1 - fac;
ptry = zeros(1, ndim);


for j = 1:ndim
    ptry(j) = psum(j)*fac1 - p(ihi, j)*fac2; %for each minimized parameter
end

%evaluate distance between moments and targets at trial point

parameters = [ptry, paramsFixed];

[x,fval,EXITFLAG] = fsolve(@errors_ss,x0,options,parameters, indicate);

problem = 0; 
if EXITFLAG<1 || isreal(x) ==0 || max(abs(fval))>1e-10
    display(['Model did not solve.'])
    beep;
    problem = 1; 
    
else 
 %display(['Model did solve. x0 updated'])
  x0 =x;  
end

[errorsTry, moments1, guess2] = errors_ss(x0,parameters, indicate);

temp = errors_sspg(guess2, parameters, indicate, 0.928);
if max(isnan(temp) ==1)
    guess2(2) =0.5;
    temp = errors_sspg(guess2, parameters, indicate, 0.928);
    if max(isnan(temp) ==1)
        guess2(2) =0.75; 
        temp = errors_sspg(guess2, parameters, indicate, 0.928);
        if max(isnan(temp) ==1)
            guess2(2) =0.25;
        end
    end 
end

if indicate ==0
    [xstar, fval, exitflag] = fsolve(@errors_sspg, guess2,options,parameters, indicate, .928);
    [errorsTry, moments2] = errors_sspg(xstar,parameters, indicate, 0.928);
    moments = [moments1(1:7), moments2];
    %moments =moments1;
else
    moments = moments1;
end

ytry = (targets - moments)*(targets -moments)'; %distance between targets and moments

%if ytry is better than the highest point, replace the hightest point
%with ytry

%in the case of the expansion, replace reflected point with the
%expanded point if the expanded point is better.

if ytry < y(ihi) && problem ==0
    y(ihi) = ytry;
    p(ihi, :) = ptry;
end
psumNew = sum(p,1);
yNew = y;
pNew =p;

end 