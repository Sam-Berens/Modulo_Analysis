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
patternSim1 = cell((nSubjects),1);
patternSim2 = cell((nSubjects),1);
for iSubject=1:numel(subjectId)
    cSubjectId = subjectId(iSubject);
    %reminder that data is [12,nvox], with the 1st 6 rows being As
    [Data,pCover(iSubject)] = getTpatterns_EpiRes(G,cSubjectId,roiId);
    % Compute pairwise correlation of stimuli across voxels
    R = corr(Data);
    patternSim{iSubject} = R;

    %first chop off top triangle (just make zeros)
    lowerBig = tril(true(12),-1);
    lR = lowerBig .* R;
    aa = lR(1:6,1:6);
    ab_ba = R(7:12,1:6);
    bb = lR(7:12,7:12);
    %flip these lower triangle values to be in the top
    bb = permute(bb,[2,1]);
    % add together colocal triangles
    aa_bb = aa + bb;

    r1a = corr(H(lowerS),aa(lowerS));
    r1b = corr(H(lowerS),bb(lowerS));
    z1a = atanh(r1a);
    z1b = atanh(r1b);
    z1 = mean([z1a,z1b]);
    r2 = corr(H(rmvDiag),ab_ba(rmvDiag)); 
    z2 = atanh(r2);
    zTemplate1(iSubject,1) = z1;
    zTemplate2(iSubject,1) = z2;
    patternSim1 = aa_bb;
    patternSim2 = ab_ba;
end

zTemplate = [zTemplate1;zTemplate2];
subjectId = repmat(subjectId,[2,1]);
pCover = repmat(pCover,[2,1]);
patternSim = cat(1,patternSim1,patternSim2); 
PatternSim = table(coLocation,subjectId,pCover,zTemplate,patternSim);
return