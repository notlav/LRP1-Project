% script to run validation code across model

% for Mac
Modelfile_path = '/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/model/'; % Change the path to the file current path
Validationfile_path = '/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/validation/';% Change the path to the file current path
Int_time = 40;
Steady_time = 40;
Threshold = 0.1; % percent of changes
Model_version = 1; % 1= Original 2= Modified
[percentMatch, resultChart, BMatch, byClass] = Automated_Validation_V1(Modelfile_path, Validationfile_path, Int_time, Steady_time,Threshold, Model_version);

% accuracy
percentMatch
