% Automated Validation
function [percentAgreeLit, percentChangeAct, resultChart] = autovalidation(...
    filename, controlCond, threshold, y0, params)

% Outputs:
% percentAgreeLit = percent of agreement with literature
% percentChangeAct = percent change in activity of the output with the
%    stimulus in the relationship being validated
% results = a summary of all the validation information
[~, txt, raw] = xlsread(filename);

ids = txt(2:end, 1);
inputs = txt(2:end, 2); % Input species
inputCode = txt(2:end, 3); % Code for input
outputSpec = txt(2:end, 4); % Output species to be compared
measurement = txt(2:end, 6); % Increase/decrease/no change

% Output variables
percentChangeAct = zeros(1, length(outputSpec));
percentChangeActStr = cell(1, length(outputSpec));
prediction = cell(1, length(outputSpec));
numMatching = 0;
match = cell(1, length(outputSpec));

outputSpecIndex = zeros(1, length(measurement));
[rpar, tau, ymax, speciesNames] = params{:};   %this line unpackages the structure of params
w = rpar(1,:);  % this is the weight parameter now as an array
n = rpar(2,:);  % this is the hill coefficient parameter now as an array
EC50 = rpar(3,:);

for i = 1:length(outputSpecIndex)
    [~, outputSpecIndex(i)] = ismember(outputSpec{i}, speciesNames);
end

% Run each simulation
for i=1:length(inputCode) 
    disp(i);
    % Unpack parameters
    [params,y0] = NetfluxODE_loadParams();
    [rpar, tau, ymax, speciesNames] = params{:}; 
    w = rpar(1,:);
    n = rpar(2,:);
    EC50 = rpar(3,:);
    
    % Control input
    eval(controlCond);
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,speciesNames};
    tspan = [0 300];
    options = [];
    [~,y] = ode15s(@NetfluxODE,tspan,y0,options,params);
    control_y = y(end,:);
    
    % Experimental input
    eval(inputCode{i});
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,speciesNames};
    tspan = [0 300];
    options = [];
    [~,y] = ode15s(@NetfluxODE,tspan,y0,options,params);
    exp_y = y(end,:);
    
    % Calculate change in output
    percentChangeAct(i) = real(exp_y(outputSpecIndex(i)))/...
        real(control_y(outputSpecIndex(i)));
    
    if percentChangeAct(i) > (threshold + 1)
        prediction{i} = 'Increase';
    elseif percentChangeAct(i) < (1 - threshold)
        prediction{i} = 'Decrease';
    else
        prediction{i} = 'No Change';
    end
    
    percentChangeActStr{i} = num2str(percentChangeAct(i));
    if isequal(prediction{i}, measurement{i})
        numMatching = numMatching + 1;
        match{i} = 1;
    else
        match{i} = 0;
    end
end

disp('end of loop');
percentAgreeLit = numMatching/length(measurement)*100;
resultChart = {inputs, outputSpec, measurement, prediction', percentChangeActStr', match'};
resultChart = horzcat(resultChart{:});
header = { 'input' , 'output', 'measurement', 'prediction','predictedChange', 'match'};
resultChart = vertcat(header, resultChart);
end