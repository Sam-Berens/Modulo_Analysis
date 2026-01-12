function [modelTb] = getWAIC(model)
%add nicer checking of below
model = string(model);
%construct model-specific function names:
if strcmp(model,'vonMises')
    f = @vM_ll;
elseif strcmp(model,'Binomial')
    f = @BN_ll;
end
model = char(model);
G ='G0';
dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.model = fullfile(dirs.Data,'_Group','A00',model);
subjectIds = getSubjectIds(G);

strct = load([dirs.model,filesep,'Output.mat'],'PostSamps');
allPSs = strct.PostSamps;
nSubs = numel(subjectIds);
nSamps = 8000;
WAIC = nan(nSubs,1);
seElpd = nan(nSubs,1);
seWAIC = nan(nSubs,1);

for iSub=1:nSubs
    cSubId = subjectIds(iSub);
    %WAIC step 1: calculate sample x trial mat of loglikelihoods 
    %we need the samples for the parameters
    cPSs = allPSs(allPSs.subjectId==cSubId,:);
    %and the trial data 
    tDataFn = fullfile(dirs.Data,char(cSubId),'Analysis','A00','InputData.json');
    trialData = jsondecode(fileread(tDataFn));
    Lls = makeLls(f,trialData,cPSs);
    
    %save likelihoods
    llFn= extractBefore(string(tDataFn),"Input");
    llFn = [char(llFn),model,filesep,'LLs.mat'];
    save(llFn,"Lls");

    %step 2: pointwise lppd (log of the mean of likelihood)
    %either we unlog everything again or we dont log in the last step?
    Ls = exp(Lls);
    lppds = log(mean(Ls,1)); %samples are dim 1, trials are dim 2
    %step 3: get the point-wise penalty (variance of the log likelihood for
    %every sample, on a given trial)
    waicPs = var(Lls,[],1);
    %step 4: grand lppd
    lppd = sum(lppds);
    %step 5: grand penalty 
    waicP = sum(waicPs);
    %step 6: Expected Log Pointwise Predictive Density 
    elpd = lppd - waicP;
    WAIC(iSub) = -2 * elpd;
    %step 7: Stadard error 
    elpds = lppds - waicPs;
    seElpd(iSub) = sqrt((nSamps*var(elpds)));
    seWAIC(iSub) = 2 * seElpd(iSub);
end

%Format /save table , add subIDs etc
modelTb = table(subjectIds,WAIC,seElpd,seWAIC);
modelTb.Properties.VariableNames = {'subjectId',...
[model,'_WAIC'],[model,'_seElpd'],[model,'_seWAIC']};

return

function [LLs] = makeLls(f,dt,PSs)
nSamples = size(PSs,1);
getTheta = @(x,c) mod((x-c).*(pi/3),2*pi);
LLs = nan(nSamples,dt.nTrials);
for iTrial = 1:dt.nTrials
    cPairId = dt.pairId(iTrial);
    column = ['a_',int2str(cPairId)];
    cASamples = PSs.(column);
    column = ['b_',int2str(cPairId)];
    cBSamples = PSs.(column);
    attempts = dt.Y(iTrial,:);
    c = dt.c(iTrial);
    theta = getTheta(attempts,c);
    xSup = dt.x(iTrial);
    params = [cASamples,cBSamples];
    %get the models log-likelihoods via model-specific function
    clls = f(params,xSup,theta); 
    %instead of doing this we need to do the operation to the vectors at
    %the same time
    LLs(:,iTrial) = clls';
end 
return




