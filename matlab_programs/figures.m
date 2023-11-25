%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: figures.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Creates  figres for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear all;

%Load paths 
load paths

options = optimset('Algorithm','trust-region-dogleg','TolFun',10e-15,'TolX',10e-15, ...
    'MaxFunEvals', 10000, 'MaxIter', 10000, 'Display','none');

%Set pic = 1 to save figures
pic =1; 

%% Section 3: Figures 4 and 5

%Set parameter values for plotting
eta = 0.6;
R = 1;
Omega_m =10; 
Am = 9.98;
v = 1;
gamma = 0.5;
lambda = 1.7325; 

%Find vstar
fun = @(v) (v/R - (1-v)/(R + 1 - v));
vstar = fzero(fun,0.5);

[zstar fval, exitflag] = fsolve(@findzstar,2,options,lambda, eta, gamma, v, Am, vstar, Omega_m, R);


vvec = [.01:0.01:1];
ind1 = find(vvec < vstar, 1,  'last');

Km0 = zeros(size(vvec)); 
Km1_l = Km0; 
Km1_h = Km0; 
Ks1_l = Km0;
Ks1_h = Km0; 

for i = 1:length(vvec)
   x0 = 40;
    [zstar fval, exitflag] = fsolve(@findzstar,x0,options,lambda, eta, gamma, vvec(i), Am, vstar, Omega_m, R);
    zstar_vec(i) = zstar;
    N =1;
    Zbarm_vec(i) = (lambda/(lambda-1))*zstar;
    Nm_N = zstar^(-lambda); %Fraction of modern firms
    Nm_vec(i) = Nm_N*N;
    Km0(i) = (1-gamma)*Nm_vec(i)*Zbarm_vec(i)*(vvec(i)*eta*Am/R).^(1/(1-eta));
    Km1_l(i) = gamma*Nm_vec(i)*Zbarm_vec(i)*(eta*Am./(R + R  + 1-vvec(i))).^(1/(1-eta));
    Km1_h(i) = gamma*Nm_vec(i)*Zbarm_vec(i)*(vvec(i)*eta*Am/R).^(1/(1-eta));
    Ks1_l(i) = Km1_l(i);
    Ks1_h(i) = gamma*Nm_vec(i)*Zbarm_vec(i)*((1-vvec(i))*eta*Am./(R + 1-vvec(i))).^(1/(1-eta));
end


Km1 = [Km1_l(1:ind1), Km1_h(ind1+1:end)]; 
Ks1 = [Ks1_l(1:ind1), Ks1_h(ind1+1:end)];

%Aggregate output across all modern firms
Ym_on_agg = (gamma*Zbarm_vec.*Nm_vec*Am).^(1-eta).*Km1.^eta + ((1-gamma)*Zbarm_vec.*Nm_vec*Am).^(1-eta).*Km0.^eta; 
Ym_off_agg =(gamma*Zbarm_vec.*Nm_vec*Am).^(1-eta).*Ks1.^eta;
Ym_agg  = vvec.*Ym_on_agg + (1-vvec).*Ym_off_agg; 

%For firm with a particular z. 
Nm =1;
Zbarm = 1; 
Km0 = (1-gamma)*Nm*Zbarm*(vvec*eta*Am/R).^(1/(1-eta));
Km1_l = gamma*Nm*Zbarm*(eta*Am./(R + R  + 1-vvec)).^(1/(1-eta));
Km1_h = gamma*Nm*Zbarm*(vvec*eta*Am/R).^(1/(1-eta));
Ks1_l = Km1_l;
Ks1_h = gamma*Nm*Zbarm*((1-vvec)*eta*Am./(R + 1-vvec)).^(1/(1-eta));

Km1 = [Km1_l(1:ind1), Km1_h(ind1+1:end)]; 
Ks1 = [Ks1_l(1:ind1), Ks1_h(ind1+1:end)];

Ym_on_1 = (Am)^(1-eta)*(Km1/(gamma*Nm*Zbarm)).^eta; 
Ym_off_1 =(Am)^(1-eta)*(Ks1/(gamma*Nm*Zbarm)).^eta;
Ym_on_0 = (Am)^(1-eta)*(Km0/((1-gamma)*Nm*Zbarm)).^eta; 
Ym_off_0 = zeros(size(Ym_on_0)); 
s2 = 1; 
s3 =61;

% Figure 4: Solution to the modern firm's problem
figure(1)
hax = axes; 
set(gca, 'FontSize', 12)
hold on 
plot(vvec(s2:s3), Ym_on_1(s2:s3) ,'Color', [0, 0.4, 0], 'LineWidth', 1.5);
plot(vvec(s2:s3), [Ym_off_1(s2:ind1), [Ym_off_1(ind1):-(Ym_off_1(ind1) - 0)/(s3 - ind1-1):0]], 'o','Color', [0, 0.4, 0], 'LineWidth', 1.5, 'MarkerIndices',1:5:length(vvec(s2:s3))); 
plot(vvec(s2:s3), Ym_on_0(s2:s3) , '--', 'Color', [0, 0, 0.8], 'LineWidth', 2.5);
plot(vvec(s2:s3), Ym_off_0(s2:s3), 'x','Color', [0, 0, 0.8], 'LineWidth',  1.5, 'MarkerIndices',1:5:length(vvec(s2:s3))); 
plot(vvec(19),Ym_on_1(19),'.', 'MarkerSize', 20, 'Color', [0,0,0])
plot(vvec(19),Ym_on_0(19),'.', 'MarkerSize', 20, 'Color', [0,0,0])
plot(vvec(19),Ym_off_0(19),'.', 'MarkerSize', 20, 'Color', [0,0,0])
plot(vvec(s3),Ym_on_0(s3),'.', 'MarkerSize', 20, 'Color', [0,0,0])
text(.95*vvec(19),Ym_on_1(19)*1.1,'A', 'FontSize', 12)
text(.95*vvec(19),Ym_on_0(19)*1.3,'B', 'FontSize', 12)
text(.95*vvec(19),Ym_off_0(19)+1,'C', 'FontSize', 12)
text(.95*vvec(s3),Ym_on_0(s3)+.5,'D', 'FontSize', 12)
legend boxoff
xlim([vvec(s2), vvec(s3)])
xlabel('Probability of grid power: v')
set(gca,'YTick',[]);
set(gca,'YTick',[]);
xticks([vvec(s2), vvec(19),  vstar, vvec(s3)])
xticklabels({'0', 'v_1' 'v*', '1'})
line([vstar, vstar], get(hax, 'YLim'), 'Color', [0,0,0], 'LineStyle', '--')
legend('Power on: generator', 'Power off: generator', 'Power on: no generator', 'Power off: no generator', 'Location', 'NW')
legend boxoff
ylabel('Output')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'model_graph2b'],'epsc');
end


% Figure 5: Short-Run Partial-Equilibrium and Long-Run General-Equilibrium
% Effects of Outages
xvec = [0:0.01:1]; 
yvec2 =xvec.^4; 
yvec1 = xvec.^2;

figure(2)
hax = axes; 
set(gca, 'FontSize', 12)
hold on 
plot(xvec, yvec2, '--','Color', [0, 0, 0.8], 'LineWidth', 1.5); 
plot(xvec, yvec1 ,'Color', [0, 0.4, 0], 'LineWidth', 1.5);
legend boxoff
xlim([0, 1])
xlabel('Probability of grid power: v')
set(gca,'YTick',[]);
set(gca,'YTick',[]);
xticks([0, xvec(60),  1])
xticklabels({'0','v_1', '1'})
x1 = [0.59 0.59];
y1 = [0.25 0.39];
annotation('line',x1,y1)
x2 = [0.4 0.58];
y2 = [0.3 0.3];
annotation('textarrow',x2,y2, 'String', 'S.R. P.E.')
x3 = [0.6 .91];
y3 = [0.23 0.23];
annotation('line',x3,y3)
x4 = [.91 .91];
y4 = [0.235 .91];
annotation('line',x4,y4)
x5 = [0.5 0.9];
y5 = [0.6 0.6];
annotation('textarrow',x5,y5, 'String', 'L.R. G.E.')
plot(xvec(60),yvec1(60),'.', 'MarkerSize', 25, 'Color', [0,0,0])
plot(xvec(60),yvec2(60),'.', 'MarkerSize', 25, 'Color', [0,0,0])
plot(xvec(end),yvec1(end),'.', 'MarkerSize', 25, 'Color', [0,0,0])
text(.92*xvec(60),yvec1(60)*1.1,'C', 'FontSize', 13)
text(.92*xvec(60),yvec2(60)*1.1,'A', 'FontSize', 13)
text(.95*xvec(end),yvec1(end),'B', 'FontSize', 13)
legend('Aggregate modern output', 'Aggregate modern output, power on ex-post', 'Location', 'NW')
legend boxoff
ylabel('Output')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'model_graph3'],'epsc');
end



%% Figure 6: Long Run General Equilibrium Effects of Eliminating Outages
load decomp
load partial_eq;
countries2 = {'GHA', 'NGA', 'TZA', 'UGA'};
xnames ={'GHA', 'NGA', 'TZA', 'UGA'};
X = 1:1:length(countries2);

vars1 ={'firm_expand', 'firm_entry',  'total', 'SR'};

for k = 1:1
    for i = 1:length(countries2)
        country1 = countries2{i};
        for j = 1:length(vars1)
            eval(strcat(vars1{j} ,'(i) =',  vars1{j}, '_', country1 ,';'));
        end
    end
end

figure(9)
bar([SR; firm_expand; firm_entry]','stacked', 'FaceColor', 'flat');
legend('SR, PE', 'LR, GE, firm-expansion channel', 'LR, GE, firm-expansion and firm-entry channels', 'Location', 'NE')
labels = arrayfun(@(value) num2str(value,'%2.1f'),total, 'UniformOutput',false);
text(X,total,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
xticklabels(xnames)
ylim([0, 30])
xlabel('Country')
box off
legend boxoff
ylabel('Percent Increase in Output per Worker')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'decomp'],'epsc');
end

%% Weak Links and tax comparison comparisons (Figures E.1 and E.2)
load  weak_links_plus_tax

%Weak links (Figure E.1)
delta_ag = (agvec./agvec(npg) -1)*100;
figure(3)
set(gca, 'Fontsize', 14)
hold on
plot(([E1vec(1:ind1); E1vec(ind2:end)]./E1vec(end)-1)*100, [Y1vec(1:ind1); Y1vec(ind2:end)]./Y1vec(end), 'Color', [0, 0.4, 0], 'LineWidth', 1.5)
plot(delta_ag, Y2vec./Y2vec(end), '--', 'Color', [0, 0, 0.8], 'LineWidth', 1.5)
xlabel('Electricity Supply Shock')
%ylabel('Output per worker')
legend({'Outages economy', 'Competitive economy'}, 'Location', 'SE', 'FontSize', 15)
box off
legend boxoff
title('Aggregate Output per Worker')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'weak_y'],'epsc');
end

figure(4)
set(gca, 'Fontsize', 14)
hold on
plot(([E1vec(1:ind1); E1vec(ind2:end)]./E1vec(end)-1)*100, ([Kg1vec(1:ind1); Kg1vec(ind2:end)])./(Kg1vec(end)), 'Color', [0, 0.4, 0], 'LineWidth', 1.5)
%plot((E1vec./E1vec(end)-1)*100, Ks1vec./Ks1vec(end-1), '-.', 'Color', [0.3, 0, 0.6], 'LineWidth', 1.5)
plot(delta_ag, Kg2vec./Kg2vec(end), '--', 'Color', [0, 0, 0.8], 'LineWidth', 1.5)
xlabel('Electricity Supply Shock')
%legend('Outages economy', 'Competitive economy', 'Location', 'SE')
box off
title('Grid-Electricity Capital per Worker')
%title('Electricity Capital per Worker')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'weak_k'],'epsc');
end

%Tax (Figure E.2)
figure(5)
set(gca, 'Fontsize', 14)
hold on
plot(([pgvec(1:ind1); pgvec(ind2:end)]./pgvec(end)-1)*100, [Y1vec(1:ind1); Y1vec(ind2:end)]./Y1vec(end), 'Color', [0, 0.4, 0], 'LineWidth', 1.5)
plot(([pgvec(1:end)]./pgvec(end)-1)*100, Y3vec./Y3vec(end), '--', 'Color', [0, 0, 0.8], 'LineWidth', 1.5)
xlabel('Electricity Price Shock')
%ylabel('Output per worker')
legend({'Outages economy', 'Competitive economy with tax'}, 'Location', 'SE', 'FontSize', 15)
box off
legend boxoff
xlim([-35, 0])
xticks([-45:5:0])
title('Aggregate Output per Worker')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'tax_y'],'epsc');
end

figure(6)
set(gca, 'Fontsize', 14)
hold on
plot(([pgvec(1:ind1); pgvec(ind2:end)]./pgvec(end)-1)*100, ([Kg1vec(1:ind1); Kg1vec(ind2:end)])./(Kg1vec(end)), 'Color', [0, 0.4, 0], 'LineWidth', 1.5)
plot(([pgvec(1:end)]./pgvec(end)-1)*100, Kg3vec./Kg3vec(end), '--', 'Color', [0, 0, 0.8], 'LineWidth', 1.5)
xlabel('Electricity Price Shock')
%legend('Outages economy', 'Competitive economy', 'Location', 'SE')
box off
xlim([-35, 0])
xticks([-45:5:0])
title('Grid-Electricity Capital per Worker')
%title('Electricity Capital per Worker')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'tax_k'],'epsc');
end







