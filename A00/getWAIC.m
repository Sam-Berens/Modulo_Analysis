function [wAIC] = getWAIC(G,model)

% Get subject list
subjectIds = getSubjectIds(G);

% Set the model-specific log-likelihood
if strcmp(model,'vonMises')
    llFunc = @spvm_ll;
elseif strcmp(model,'Binomial')
    llFunc = @spbn_ll;
end

% Set some dirs
dirs.Data = fullfile('..','..','Data');
dirs.Model = fullfile(dirs.Data,'_Group','A00',model);

stanOut = load([dirs.Model,filesep,'Output.mat']);
PostSamps = stanOut.PostSamps;
nSubjects = numel(subjectIds);
wAIC = nan(nSubjects,1);
wAIC_se = nan(nSubjects,1);
elppd_se = nan(nSubjects,1);

for iSubject = 1:nSubjects

    cSubjectId = subjectIds(iSubject);
    inputDataFn = fullfile(dirs.Data,char(cSubjectId),...
        'Analysis','A00','InputData.json');
    inputData = jsondecode(fileread(inputDataFn));

    % Step 1: Calculate sample x trial matrix of loglikelihoods
    cSamples = PostSamps(PostSamps.subjectId==cSubjectId,:);
    Ll = makeLl(llFunc,inputData,cSamples);

    % Step 2: Compute pointwise lppd (log of the mean of likelihood)
    % dim1 = Samples
    % dim2 = Trials
    m = max(Ll,[],1);
    lppd = m + log(mean(exp(Ll - m),1));

    % Step 3: Compute the point-wise penalty ...
    % ... (variance of the log likelihood for every sample on each trial).
    ppp = var(Ll,[],1);

    % Step 4: Expected Log Pointwise Predictive Density
    elppd = sum(lppd) - sum(ppp);

    % Step 5: Compute the WAIC
    wAIC(iSubject) = -2 * elppd;

    % Step 6: Compute the stadard errors
    elpds = lppd - ppp;
    nTrials = numel(elpds);
    elppd_se(iSubject) = sqrt(nTrials * var(elpds));
    wAIC_se(iSubject) = 2 * elppd_se(iSubject);
end

% Format the output table
wAIC = table(subjectIds,wAIC,wAIC_se,elppd_se);
wAIC.Properties.VariableNames = {'subjectId',...
    [model,'_wAIC'],[model,'_wAIC_se'],[model,'_elppd_se']};

return

function [L] = makeLl(llFun,inputData,PosSamps)
nSamples = size(PosSamps,1);
x = inputData.x;
nTrials = numel(x);
L = nan(nSamples,nTrials);
getTheta = @(yy,t) mod((yy-t).*(pi/3),2*pi);

colName = arrayfun(...
    @(ii)sprintf('a_%i',ii),...
    inputData.pairId,...
    'UniformOutput',false);
colIdx = cellfun(...
    @(x)find(strcmp(x,PosSamps.Properties.VariableNames)),...
    colName);
A = PosSamps{:,colIdx};

colName = arrayfun(...
    @(ii)sprintf('b_%i',ii),...
    inputData.pairId,...
    'UniformOutput',false);
colIdx = cellfun(...
    @(x)find(strcmp(x,PosSamps.Properties.VariableNames)),...
    colName);
B = PosSamps{:,colIdx};


for iTrial = 1:numel(x)
    x_ii = x(iTrial,:);
    target = inputData.c(iTrial);
    y = inputData.Y(iTrial,:);
    theta = getTheta(y,target);
    P = [A(:,iTrial),B(:,iTrial)];
    L(:,iTrial) = llFun(P,x_ii,theta);
end
return