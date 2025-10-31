function [] = z06_MakeRPs(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir.folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];

rpsFolder = [dataDir,filesep,'RPs'];
if ~exist(rpsFolder,'dir')
    mkdir(rpsFolder);
end

censorThresh = 0.5; % mm FWD

rpList = dir([epiDir,filesep,'*.txt']);
nRuns = numel(rpList);

for iR = 1:nRuns
    % Get data
    FileName_RpText = [rpList(iR).folder,filesep,rpList(iR).name];
    RPs = readmatrix(FileName_RpText);
    nEpis = size(RPs,1);

    % Compute derivatives
    dRPs = diff(RPs,1,1);
    qt = (0:1:(nEpis-1))';
    t = mean([qt(1:end-1),qt(2:end)],2);
    dRPs = cell2mat(...
        cellfun(@(y)sincInterp(qt,t,y),mat2cell(dRPs,nEpis-1,ones(1,6)),...
        'UniformOutput',false));

    % Compute frame-wise displacement (FWD)
    translations = sum(abs(dRPs(:,1:3)),2);
    rotations = sum(abs(dRPs(:,4:6)),2);
    fwd = translations + 50*rotations;

    % Compute censors
    volIdxs = find(fwd > censorThresh);
    censors = zeros(nEpis,numel(volIdxs));
    for ii = 1:numel(volIdxs)
        censors(volIdxs(ii),ii) = 1;
    end

    % Construct R
    R = [RPs,dRPs,fwd];
    R = zscore(R,[],1);
    R = [R,censors];

    % Save R
    save(sprintf('%s%sRP%i.mat',epiDir,filesep,iR),'R');
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return