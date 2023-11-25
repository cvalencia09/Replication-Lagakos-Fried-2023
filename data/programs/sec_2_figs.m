%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: sec_2_figs.m
% By: Stephie Fried and David Lagakos
% Date: Winter 2022
% Purpose: Creates Figures 1-3 for Section 2 of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

% Run the programs from the outermost folder in the replication package
% (click "add to path" after hitting run)
reppath = pwd; %Outermost folder in replication pacakage
datapath1 = strcat(reppath, '/data/output'); %path to folder with data for figure1
datapath2 = strcat(reppath, '/data/derived'); %path to folder with data for figures 2 and 3
figpath = strcat(reppath, '/figures'); % path to folder where figures are stored
bpath = strcat(reppath, '/matlab_programs'); %path to folder with matlab programs
cd(bpath);

%Flag set equal to 1 to save the figures
pic =1;
%% Figure 1 Percent of firms experiencing electricity outages
cd(datapath1)
[~, names] = xlsread('fig1_data.xlsx', 'A2:A81');
gdp_pop = xlsread('fig1_data.xlsx', 'B2:B81');
ougages = xlsread('fig1_data.xlsx', 'C2:C81');

x = log2(gdp_pop); 
y =ougages; 

X = [ones(size(x)), x];
bmat = (X'*X)^(-1)*X'*y;

yfit = bmat(2)*x + bmat(1);
fsize=8; 
color=[0 .4 0];
axis([-1 6 0 100]);
for i=1:size(y,2);
  text(x,y(:,i),names,'FontSize',fsize,'Color',color,'LineWidth',1);
end;
xticklabels({' ','1','2','4','8','16','32', '64'})
xlabel('GDP per Capita x $1000')
ylabel('Percent')
hold on
plot(x, yfit, 'LineWidth', 1, 'Color', [0, 0, .8])
box off
if pic ==1
    saveas(gcf,[figpath ,filesep, 'fig1'],'epsc');
end

%% Figure 2 and 3
countries2 = {'CIV','ETH', 'GHA', 'KEN', 'MDG', 'MOZ', 'NER', 'NGA', 'TZA', 'UGA'};

%Import data
cd(datapath2);
M1 = xlsread('ps_pg_generators_matlab.xlsx');
pg_data = M1(:, 1); ac_data = M1(:,2); gen_data = M1(:,3);

%Figure 2: Percent of firms that own or share a generator
figure(2)
bar(gen_data, 0.8,'EdgeColor', 'flat')
xt = get(gca, 'XTick');
labels = arrayfun(@(value) num2str(value,'%2.0f'),gen_data,'UniformOutput',false);
text(xt, gen_data, labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
box off
set(gca,'YLim', [0, 100], 'XLim', [0.5, 10.5])
ylabel('Percent')
xticklabels(countries2)
xlabel('Country')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'fig2'],'epsc');
end

%Figure 3: Average Cost of Electricity: Grid vs. Generator
X = 1:1:length(countries2);
figure(10)
b = bar([pg_data'; ac_data']'*100, 1,'FaceColor', 'flat');
b(1).CData =	[0, 0.4470, 0.7410];
b(2).CData =[0.4660, 0.7740, 0.1880];
xt = get(gca, 'XTick');
labels = arrayfun(@(value) num2str(value,'%2.0f'),pg_data*100,'UniformOutput',false);
text(xt-.18, pg_data*100, labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
labels = arrayfun(@(value) num2str(value,'%2.0f'),ac_data*100,'UniformOutput',false);
text(xt+.15, ac_data*100, labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
box off
set(gca,'YLim', [0, 55], 'XLim',[0 10.5])
ylabel('Cents/kwh (2014 $)')
xticklabels(countries2)
legend boxoff
xlabel('Country')
legend('Grid', 'Generator', 'Location', 'NE', 'Orientation', 'horizontal')
if pic ==1
    saveas(gcf,[figpath ,filesep, 'fig3'],'epsc');
end

