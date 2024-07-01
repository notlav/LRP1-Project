% Mechanistic subnetwork identification
% Written by Taylor Eggertsen
% Last updated 10/25/2022

%% Inputs for simulations
clear all; close all; clc

model = 'exampleNet.xlsx';

formattedReactions = table;
% Species/Reaction information from network, 'species'/'reactions' tab 
networkReactions = readtable(model, 'Sheet', 'reactions');
if ismember('=',cell2mat(networkReactions{1,3}))==0
    networkReactions(1,:)=[];  
end
% Formats network reactions to only show the product/output. 
for i = 1:height(networkReactions)
    reaction = string(networkReactions{i,3});
    nodeOfReaction = extractAfter(reaction, '=>'); 
    formattedReactions{i,1} = strtrim(nodeOfReaction);
end

% Import parameter and ode files
[params,y0] = NetfluxODE_Toy_loadParams;
[rpar,tau,ymax,speciesNames]=params{:};
w = rpar(1,:);
n = rpar(2,:);
EC50 = rpar(3,:);

% Baseline condition (set if different than excel file values)
w(1) = 0.5;
w(2) = 0.5;

% Perturbation 
input = 'ymax(2)=0;';

% Phenotype measured
phenotype = 'E';
phen = find(contains(speciesNames,phenotype)); 

%% Control Knockouts

rpar = [w;n;EC50];
tspan = [0 50]; 
options = [];
params = {rpar,tau,ymax,speciesNames};
[t,y] = ode15s(@NetfluxODE_Toy,tspan,y0,options,params);
y0 = real(y(end,:)');

% simulate knockdowns
for r = 0:length(speciesNames)
    disp(num2str(r))
    ymax_Knockout = ymax;
    if r > 0
        ymax_Knockout(r) = 0.01;
    end

    params = {rpar,tau,ymax_Knockout,speciesNames};
    tspan = [0 50]; 
    options = []; 
    [t2,y2] = ode15s(@NetfluxODE_Toy,tspan,y0,options,params); 

    % Calculate final values 
    if r == 0
        ySim = real(y2(:,phen))';
        [c2] = ySim;
        c2End_control = c2(end);
        c2End_knockout = c2(end);
        label = 'Control';
    else
        ySim = real(y2(:,phen))';
        [c2] = ySim;
        c2End_knockout = c2(end);
        label = speciesNames{r};
    end

    % Calculates the change in cell area from the knockout (in the
    % presence of stimulus alone) and control (also in presence
    % of stimulus alone)
    dataVector_control = {label, c2End_knockout-c2End_control};
    knockoutData_control(r+1, :) = dataVector_control;
end

%% Perturbation knockdowns

eval(input);
params = {rpar,tau,ymax_Knockout,speciesNames};

% simulate knockdowns
for r = 0:length(speciesNames)
    disp(num2str(r));
    ymax_Knockout = ymax;
    if r > 0
        ymax_Knockout(r) = 0.01;
    end
    params = {rpar,tau,ymax_Knockout,speciesNames};
    tspan = [0 50]; 
    options = []; 
    [t2,y2] = ode15s(@NetfluxODE_Toy,tspan,y0,options,params); 

    % Calculate final values
    if r == 0
        T_Fin = real(y2(end,:)'); T0 = real(y2(1,:)');
        T_Diff = T_Fin-T0;
        ySim = real(y2(:,phen))';
        [c1] = ySim;
        c1End_control = c1(end);
        c1End_knockout = c1(end);
        label = 'Control';
    else
        ySim = real(y2(:,phen))';
        [c1] = ySim;
        c1End_knockout = c1(end);
        label = speciesNames{r};
    end
    
    % Calculates the change in phenotype from the knockout (in the
    % presence of stimulus and drug) and control (also in presence
    % of stimulus and drug)
    dataVector = {label, c1End_knockout-c1End_control};
    knockoutData_perturb(r+1, :) = dataVector;
end

%% Identify subnetwork 

% Calculate the 'difference of the difference' - This finds how the
% changes in phenotype in each analysis differ from each other,
% resulting in a measurement of how each node contributes to the
% mechanism of the perturbation
knockoutData = cell2mat(knockoutData_perturb(:,2))-cell2mat(knockoutData_control(:,2));

%Filter by nodes that were affected by drug action
for ii=1:length(T_Diff)
    if abs(T_Diff(ii)) < 0.01
        knockoutData(ii+1) = 0;
    end
end

%Filter by effect size (identifies the most important nodes)
thresh = 0.1*max(abs(knockoutData));
for ii = 1:length(knockoutData)
    if abs(knockoutData(ii)) < thresh
        knockoutData(ii) = 0;
    end
end

%% Figures and data export

X = categorical(knockoutData_control(:,1)); 
X = reordercats(X,string(X));

figure
bar(X,cell2mat(knockoutData_control(:,2)))
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', X,'FontSize', 14);
xtickangle(90)
ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' (KD - Baseline)'))
xlabel('Knockdowns')

figure
bar(X,cell2mat(knockoutData_perturb(:,2)))
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', X,'FontSize', 14);
xtickangle(90)
ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' ((Drug+KD) - (Drug))'))
xlabel('Knockdowns')

figure
bar(X,knockoutData)
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', X,'FontSize', 14);
xtickangle(90)
ylabel(strcat('Difference of KD Effect on',{' '},char(speciesNames(phen))))
xlabel('Knockdowns')

%write data to text file
T = table(X,knockoutData,'VariableNames',{'Species_name','KDdata'}); 
T(1,:) = []; filename = strcat('B_to_',char(speciesNames(phen)),'_KDdata.txt');
writetable(T, filename)
%write data to text file
X(1) = [];
S = table(X,T0,T_Fin,T_Fin-T0,'VariableNames',{'Species_name','T0','T_Fin','T_Diff'}); 
filename = strcat('B_to_',char(speciesNames(phen)),'.txt');
writetable(S, filename)

%% Cytoscape instructions

% import network structure into cytoscape
% import data (txt) files into cytoscape

% use KDdata to filter critical nodes involved in perturbed activity:
    % in Filter sidebar add
        % Node: KDdata is not between (range of critical values for nodes)
        % Node: Type doesn't contain 'connector'
% Hide selected nodes

% use T_Diff to color remaining nodes
