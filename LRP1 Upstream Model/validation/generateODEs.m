function a = generateODEs(model)

%throw +util package in path, may have to update "Z","Y", etc
addpath('Z:\Malcolm\Neutrophil Model\Netflux-master\Netflux')
addpath('Y:\Malcolm\Neutrophil Model\Netflux-master\Netflux'); %depends on if remote or not

% Generate the ODE and parameter files from NETFLUX
if exist([pwd '\NetfluxODE.m'],'file') == 2
    delete('NetfluxODE.m');
end
if exist([pwd '\NetfluxODE_loadParams.m'],'file') == 2
    delete('NetfluxODE_loadParams.m');
end
namepos = findstr('.xls', model); namestr = cellstr(model(1:namepos-1));

[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,model);
commandODE = util.exportODE(specID,paramList,ODElist, 'NetfluxODE');
[a,commandPARAM,b] = util.exportODE(specID,paramList,ODElist,'NetfluxODE');
util.textwrite('NetfluxODE.m',commandODE);
util.textwrite('NetfluxODE_loadParams.m',commandPARAM);

a = 'complete'
end

