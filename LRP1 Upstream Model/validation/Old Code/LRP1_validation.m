% NetfluxODE_run.m

 model = 'paper-apoptosis-model.xlsx'
 a = generateODEs(model)

 filename = 'paper-validation-model.xlsx'; 
 
 % Load parameters
 [params,y0] = NetfluxODE_loadParams(); %will need to update when change model
 %disp(params{1})
 %make sure validation indices are right while updating model

 %Finding steady state
 temp = y0;
 options = [];

 %for i=0:1:4 % incubation is 4 hours NOT 1 day but have also seen 10min?
 %let run for steady state for 
 %for i = 1:4 %similar time point to ML
     %no steady state
     %[t, y] = ode23(@NetfluxODE,[0:0.001:0.10],temp,options,params); %steady state, 4 hrs, ML does 10 minutes
     [t, y] = ode15s(@NetfluxODE,[0 4],temp,options,params); %this is what using for ML data not 10 min
     y = real(y);
     %possible may want to change for murine, human
     %remember to change length stimulation for human and neutrophils
     %  y_diff = y(end,:) - temp;
   %  y_prime = y(end,:) .* 0.001;
   %  if isequal((y_diff < y_prime), ones(size(y_prime)))
    %    break
   %  end
 %    temp = y(end,:);
% end

 y0 = y(end, :); %no y0 if no steady state
 
 %filename = 'neutrophil_validation_Spring_Clean.xlsx';
 %filename = 'Human neutrophil validation.xlsx'; %check incubation times and stimulation times
 %and stim times for these, humans stim for about 3 hours
 %filename = 'Murine neutrophil validation.xlsx'; %check incubation times, %stimulate mice for 3h
 %and stimulation times
 %and stim times for these
 %add folder with validation excel files
 %threshold 0.1
 %threshold may want to lower

 %reactionVars = params(1,1); %get weight, n, EC50
 reactionWeightsOriginal = params{1,1}(1,:); %index all weights of reaction
 %now can see what happens to validation when take out one connection at a
 %time

 %for loop not quite working
 % tic
 % for i = 1:length(reactionWeightsOriginal)
 %     disp(i)
 %     reactionWeights = reactionWeightsOriginal;  %set params back to normal start of each loop
 %     reactionWeights(i) = 0;
 %     params{1,1}(1,:) = reactionWeights; %index and change just weight

 %  Auto validate model
     [percentAgreeLit, percentChangeAct, resultChart] = autovalidation(filename, '', 0.1, y0, params);

%      accuracy(i) = percentAgreeLit;
%      percentChangeActivity{i} = percentChangeAct; %storing a matrix every iteration
%      rawresults{i} = resultChart; %storing a matrix, not single double every iter
% 
%  end
% toc

 %all inputs 0, find values of nodes
%params{4} should work but calls fibroblast model?
     
%      tspan = [0 4];
%      [params, y0] = NetfluxODE_loadParams();
%      [rpar, tau, ymax, speciesNames] = params{:}; %unpack parameters
%      w = rpar(1,:); n = rpar(2,:); EC50 = rpar(3,:);
%      %w(1) = 0; 
%      %w(2) = 0; 
%     % w(3) = 0; 
%      w(4) = 1;
%      rpar = [w;n;EC50]; %have to repack it
%      paramsNew = {rpar, tau, ymax, speciesNames};
%      
%      [tZero,yZero] = ode23(@NetfluxODE,tspan,y0,options,paramsNew);

 
 function output = drugApp(cmd1, cmd2, timeInt1, timeInt2, params, y0)
 [rpar,tau,ymax,speciesNames]=params{:}; 
 w = rpar(1,:);
 n = rpar(2,:);
 EC50 = rpar(3,:);
 eval(cmd1);
 rpar = [w;n;EC50];
 params = {rpar,tau,ymax,speciesNames};
 options = [];
 [t,y] = ode23(@NetfluxODE,timeInt1,y0,options,params);
 new_y = real(interp1(t, y, (0:1:timeInt1(2))));
 
 y0 = new_y(end,:);
 eval(cmd2);
 rpar = [w;n;EC50];
 params = {rpar,tau,ymax,speciesNames};
 options = [];
 [t1,y1] = ode23(@NetfluxODE,timeInt2,y0,options,params);
 new_y1 = real(interp1(t1, y1, (0:1:timeInt2(2))));

 output = [new_y.' new_y1.'];
 end
 
 function im = createImage(output, speciesNames)
 % Angela's method
 figure
 colormap(flip(bone));
 i1 = imagesc(output*100);
 set(gca, 'XTick', 0:24:144);
 set(gca, 'XTickLabel', 0:24:144);
 set(gca, 'YTick', 1:length(speciesNames));
 set(gca, 'YTickLabel', speciesNames);
 c1 = colorbar('southoutside');
 c1.Label.String = 'Activity %';
 c1.Label.FontSize = 8;
 c1.Position = [0.75 0.05 0.15 0.03];
 end
 
function im1 = createSensAnalysis(output, speciesNames)
% Angela's method
figure
cmaprange = [0.5:0.01:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
maxVal = 0.25;
caxis([-maxVal, maxVal]);
imagesc(output,[-maxVal,maxVal]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',8);
xticklabel_rotate;
xlabel('Perturbation');
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',speciesNames,'fontsize',7);
ylabel('Sensitivity of Outputs');
title(['Sensitivity Analysis']);
colorbar;
end
