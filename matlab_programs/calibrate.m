%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: calibrate.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Uses Nelder-Mead minnimization routine to calibrates model.
% Adopted from Numerical Recipes amoeba (Press et al. 1988, pg 411)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, pStar, x0] = calibrate(p0, paramsFixed, targets, x0,options, indicate)
%%
%COMPUTATIONAL PARAMETERS
tiny = 1e-10;
ftol  = 1e-10;

%DEFINE INITIAL SIMPLEX
pp=p0; 
ndim = length(pp); %number of parameters
step = -.01*pp; %step size is 5%
p = [pp; repmat(pp, ndim, 1)+diag(step)]; %Each row is a vertex of the starting simplex
y = zeros(length(p0) +1, 1);

%INIATILIZE NELDER MEAD LOOP
for i = 1:length(p0) +1;
    i
   [y(i), x0]= nelder2(p(i, :), x0, paramsFixed, targets, options, indicate);
end

done =0;
psum = sum(p, 1);
iter =0;
%%%%%%%% NELDER MEAD LOOP %%%%%%%%%
%%
while(done ==0)
    %%
    
    iter = iter +1;
    
    %FIND THE HIGHEST, SECOND HIGHEST AND LOWEST POINTS ON THE SIMPLEX
    ihi = find(y == max(y), 1, 'first'); %highest 
    y1 = y(y ~= max(y));
    inhi = find(y == max(y1), 1, 'first'); %second highest
    ilo = find(y == min(y), 1, 'first');  %lowest
    
    rtol = 2*(y(ihi) - y(ilo))/(y(ihi) + y(ilo) + tiny); %stopping criteria
    if rtol < ftol || y(ilo) < 1e-10 %squared distance between targets and moments is less than 1e-10
        done=1; 
    end 
    
    %BEGIN A NEW ITERATION
    %reflect the simplex from ihi
    fac =-1;
    [ytry, x0, y, p, psum] =nelder(p, y, psum, ndim, ihi,...
           fac, x0, paramsFixed, targets, options, indicate); 
            
    %if ytry is better than the best point, try an expansion by a factor of
    %two about ihi (which is the new point).  (expansion)
    
    if ytry <= y(ilo)
        fac =2;
        [ytry, x0, y, p, psum, problem] =nelder(p, y, psum, ndim, ihi,...
           fac, x0, paramsFixed, targets, options, indicate);
           
        %if ytry is worse than the second highest point,do a contraction

    elseif ytry >= y(inhi) 
        fac = .5;
        ysave = y(ihi);
        [ytry, x0, y, p, psum] =nelder(p, y, psum, ndim, ihi,...
           fac, x0, paramsFixed, targets, options, indicate);
        
        %If the new point is worse than ysave, contract around the best
        %point
        
        if ytry >= ysave
            disp('Many D contraction')
            for i =1:ndim+1
                for j =1:ndim
                    p(i,j) = .5*(p(i,j) + p(ilo, j));
                end
                [y(i), x0]= nelder2(p(i,:), x0, paramsFixed, targets, options, indicate);
                           
            end
            psum = sum(p,1);
        end
    end

  [rtol, y(ilo)]
  iter
  if (iter > 1000) done =1; 
  
end

pStar = p(ilo, :); 
dist = y(ilo); 

end 