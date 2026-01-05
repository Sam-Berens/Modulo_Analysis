function [DataTabel01] = getWAIC(model)

G ='G1';
dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.vnM = fullfile(dirs.Data,'_Group','A00',model);
subjectIds = getSubjectIds(G);

strct = load([dir.vnM,filesep,'Output.mat'],'PostSamps');
allPSs = strct.PostSamps;

%For each sub, make a matrix which is 8000 (samples) by nTrials

%in subject loop you select out just that subject:
for iSub=1:numel(subjectIds)
    %we need the samples for the parameters
    cPSs = allPSs(allPSs.subjectId==cSubId,:);
    %and the trial data 
    trialData = fullfile(dirs.Data,cSubId,'Analysis','A00','InputData.json');
    trialData = jsondecode(fileread(trialData));
    LLs = makeLls(trialData,cPSs);
    %TO DO ..rest of WAIC steps
    
end




return

function [LLs] = makeLls(dt,PSs)
nSamples = size(PSs,1);
fh = @(x,c) mod((x-c)*(pi/3),2*pi);
LLs = nan(nSamples,dt.nTrials);
for iTrial = 1:dt.nTrials
    cPairId = dt.pairId(iTrial);
    column = ['a_',int2str(cPairId)];
    cASamples = PSs.(column);
    column = ['b_',int2str(cPairId)];
    cBSamples = PSs.(column);
    attempts = dt.Y(iTrial);
    c = dt.c(iTrial);
    theta = fh(attempts,c);
    xSup = dt.x;
    nlls = arrayfun(@(a,b) spth_nll([a,b],xSup,theta),cASamples,cBSamples);
    LLs(iTrial,:) = (-nlls)';   %TO DO check if this works !!
end 
return







