function [PatternSim] = getPatternSim(G,roiId)
% [PatternSim] = getPatternSim(G,roiId)
% PatternSim output is a table containing:
%  - subjectId: [nSubj,1]
%  - zTemplate: [nSubj,1]
%  - pCover: [nSubj,1]
%  - patternSim: {nSubj,1}

%% Construct H (hypothesised similarity matrix)
simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
H = nan(6);
for ii = 1:36
    [x,y] = ind2sub([6,6],ii);
    H(ii) = simFun(x,y);
end

%% Get Subjets
subjectId = getSubjectIds(G);
nSubjects = numel(subjectId);

%% Do stuff
lowerS = tril(true(6),-1);
zTemplate = nan(nSubjects,1);
pCover = nan(nSubjects,1);
patternSim = cell(nSubjects,1);
fh = waitbar(0,sprintf('Computing similarity stats for %s',roiId));
for iSubject = 1:nSubjects
    cSubjectId = subjectId(iSubject);
    [Data,pCover(iSubject)] = getTpatterns_EpiRes(G,cSubjectId,roiId);
    R = corr(Data);
    zTemplate(iSubject) = atanh(corr(H(lowerS),R(lowerS)));
    patternSim{iSubject} = R;
    waitbar(iSubject/nSubjects,fh);
end
close(fh);

%% Make table
PatternSim = table(subjectId,zTemplate,pCover,patternSim);
return
