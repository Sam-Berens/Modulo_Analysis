function [PatternSim] = getPatternSim(G,roiId)
% Outputs is a table containing:
%  - subjectId: [nSubj*2,1]
%  - colocation: [nSubj*2,1]
%  - zTemplate: [nSubj*2,2]
%  - pCover: [nSubj*2,1]
%  - patternSim: {nSubj*2,1}

%% Construct H1 and H2 (hypothesised similarity matrices)
simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
H = nan(6);
for ii = 1:36
    [x,y] = ind2sub([6,6],ii);
    H(ii) = simFun(x,y);
end

% Hp1 is for the colocation = +1
Hp1 = kron(eye(2),H);
Hp1(triu(true(12))) = NaN;
Hp1(logical(kron([0,0;1,0],true(6)))) = NaN;

% Hn1 is for the colocation = -1
Hn1 = kron([0,0;1,0],H);
Hn1(triu(true(12))) = NaN;
Hn1(~isnan(Hp1)) = NaN; % Remove colocation = +1
Hn1(logical(kron([0,0;1,0],eye(6)))) = NaN; % Remove the visual effect

% Finally, remove upper triangle of H, including the diagonal
H(triu(true(6))) = NaN;

%% Get the subject list
subjectId = getSubjectIds(G);
nSubjects = numel(subjectId);

%% Preallocate
colocation = nan(nSubjects*2,1);
zTemplate = nan(nSubjects*2,1);
pCover = nan(nSubjects*2,1);
patternSim = cell(nSubjects*2,1);

%% Loopy loo
fh = waitbar(0,['Getting pattern similarity for ',roiId]);
iIn = 0;
for iSubject = 1:nSubjects
    cSubjectId = subjectId(iSubject);

    % Data is [nVox,12], with the 1st 6 cols being the 'A' position
    [Data,cpCover] = getTpatterns_EpiRes(G,cSubjectId,roiId);

    % Compute pairwise correlation of stimuli across voxels
    R = corr(Data);

    % Compute zTemplate and extract patternSim for each colocatoin
    for cl = -1:2:1
        iIn = iIn + 1;
        colocation(iIn) = cl;
        pCover(iIn) = cpCover;
        if cl == -1
            zTemplate(iIn) = nanzcorr(Hn1,R);
            Ps = nan(6);
            Ps(xor(true(6),logical(eye(6)))) = R(~isnan(Hn1));
            patternSim{iIn} = Ps;
        else
            idx = find(~isnan(Hp1));
            Ps_aa = nan(6);
            Ps_aa(tril(true(6),-1)) = R(idx(1:15));
            Ps_bb = nan(6);
            Ps_bb(tril(true(6),-1)) = R(idx(16:30));
            zTemplate(iIn) = (nanzcorr(H,Ps_aa) + nanzcorr(H_,Ps_bb))/2;
            Ps = sum(cat(3,Ps_aa,Ps_bb'),3,'omitmissing');
            Ps(logical(eye(6))) = NaN;
            patternSim{iIn} = Ps;
        end
    end

    waitbar(iSubject/nSubjects,fh);
end
close(fh);

%% Double up subjectId
subjectId = [subjectId,subjectId];
subjectId = subjectId';
subjectId = subjectId(:);

%% Make the data table
PatternSim = table(subjectId,colocation,zTemplate,zTemplate_,pCover,patternSim);
return

function [z] = nanzcorr(H,R)
S = ~isnan(H);
z = atanh(corr(H(S),R(S)));
return