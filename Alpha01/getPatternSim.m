function [PatternSim] = getPatternSim(G,roiId)
% Outputs is a table containing:
%  - subjectId: [nSubj*2,1]
%  - zTemplate: [nSubj*2,2]
%  - colocation: [nSubj*2,1]
%  - pCover: [nSubj*2,1]
%  - patternSim: {nSubj*2,1}


% Construct H (hypothesised similarity matrix)
simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
H = nan(6);
for ii = 1:36
    [x,y] = ind2sub([6,6],ii);
    H(ii) = simFun(x,y);
end


% Get Subjects
subjectId = getSubjectIds(G);
nSubjects = numel(subjectId);

lowerS = tril(true(6),-1);
rmvDiag = boolean(tril(ones(6),-1) + triu(ones(6),1));
zTemplate1 = nan(nSubjects,1);
zTemplate2 = nan(nSubjects,1);
coLocation = [ones(nSubjects,1);(-1 * ones(nSubjects,1))];
pCover = nan(nSubjects,1);
patternSim = cell(nSubjects,1);
for iSubject=1:numel(subjectId)
    cSubjectId = subjectId(iSubject);
    %reminder that data is [12,nvox], with the 1st 6 rows being as
    [Data,pCover(iSubject)] = getTpatterns_EpiRes(G,cSubjectId,roiId);

    %reminder that this M is (nStims,nVoxels)
    % Compute pairwise correlation of stimuli across voxels
    %D = pdist(Data,'euclidean');
    %D = squareform(D); 
    R = corr(Data);
    patternSim{iSubject} = R;

    %first chop off top triangle (just make zeros)
    lowerBig = tril(true(12),-1);
    lR = lowerBig .* R;
    aa = lR(1:6,1:6);
    ab_ba = R(7:12,1:6);
    bb = lR(7:12,7:12);

    r1a = corr(H(lowerS),aa(lowerS));
    r1b = corr(H(lowerS),bb(lowerS));
    z1a = atanh(r1a);
    z1b = atanh(r1b);
    z1 = mean([z1a,z1b]);
    r2 = corr(H(rmvDiag),ab_ba(rmvDiag)); 
    z2 = atanh(r2);
    patternSim{iSubject,1} = R;
    zTemplate1(iSubject,1) = z1;
    zTemplate2(iSubject,1) = z2;
end

zTemplate = [zTemplate1;zTemplate2];
subjectId = repmat(subjectId,[2,1]);
pCover = repmat(pCover,[2,1]);
patternSim =  repmat(patternSim,[2,1]);
PatternSim = table(coLocation,subjectId,pCover,zTemplate,patternSim);
return