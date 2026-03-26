function [DataTables05] = getDataTable05(G)

if exist('DataTables05.mat','file')
    temp = load('DataTables05.mat');
    DataTables05 = temp.DataTables05;
    return
end

%% Get the visual similarity imgPerm=[1,2,3,4,5,6]
temp = load('VisualSim.mat');
VisSim = temp.VisualSim;

%% Cd out
wd = pwd;
cd ..;

%% Get pNonc
pNonc = get_pNonc(G);
pNonc.cpNonc = pNonc.pNonc - mean(pNonc.pNonc);
subjectIds = pNonc.subjectId;

%% Get ROI list
roiNames = arrayfun(...
    @(ii)sprintf('N17P200_R%03d',ii),...
    (1:200)',...
    'UniformOutput',false);

%% ROI loop
DataTables05 = cell(size(roiNames));
parfor iRoi = 1:numel(roiNames)
    DataTables05{iRoi} = getRoiData(VisSim,G,subjectIds,roiNames{iRoi});
end

%% Do the outer join with pNonc
for iRoi = 1:numel(roiNames)
    DT = DataTables05{iRoi};
    DT = outerjoin(pNonc,DT);

    % Rename subjectId after outerjoin
    DT.subjectId_DT = [];
    DT.Properties.VariableNames{1} = 'subjectId';
    DataTables05{iRoi} = DT;
end

%% Cd back and save
cd(wd);
save('DataTables05.mat','DataTables05');
return

function [PatternSim] = getRoiData(VisSim,G,subjectList,roiId)
% Outputs is a table containing:
%  - subjectId: [nSubj*2,1]
%  - colocation: [nSubj*2,1]
%  - zTemplate: [nSubj*2,2]
%  - zVisual: [nSubj*2,2]
%  - pCover: [nSubj*2,1]

nSubjects = numel(subjectList);

%% Preallocate
subjectId = cell(nSubjects*2,1);
colocation = nan(nSubjects*2,1);
zTemplate = nan(nSubjects*2,1);
zVisual = nan(nSubjects*2,1);
pCover = nan(nSubjects*2,1);

%% Loopy loop
fh = waitbar(0,['Getting pattern similarity for ',roiId]);
iIn = 0;
for iSubject = 1:nSubjects
    cSubjectId = subjectList(iSubject);

    % Data is [nVox,12], with the 1st 6 cols being the 'A' position
    [Data,cpCover] = getTpatterns_EpiRes(G,cSubjectId,roiId);

    % Get subject's image perm, undo zero-ordering
    imgPerm = getImgPerm(char(cSubjectId));
    imgPerm = imgPerm + 1;
    z = mdl5Func(VisSim,Data,imgPerm);
    % z is a 3x1 vector of {-ve, +ve, vis} Fisher-transformed stats

    % Populate colocation, zTemplate and pCover
    for cl = -1:2:1
        iIn = iIn + 1;
        subjectId{iIn} = char(cSubjectId);
        colocation(iIn) = cl;
        if cl == -1
            zTemplate(iIn) = z(1);
        else
            zTemplate(iIn) = z(2);
        end
        zVisual(iIn) = z(3);
        pCover(iIn) = cpCover;
    end

    waitbar(iSubject/nSubjects,fh);
end
close(fh);
subjectId = categorical(subjectId); 

%% Make the data table
PatternSim = table(subjectId,colocation,zTemplate,zVisual,pCover);
return


function [z] = mdl5Func(VisSim,M,imgPerm)
nanzscore = @(v) (v-mean(v,'omitmissing')) ./ std(v,'omitmissing');

% Construct the basic 6x6 similarity hypothesis
simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
H_ = nan(6);
for ii = 1:36
    [x,y] = ind2sub([6,6],ii);
    H_(ii) = simFun(x,y);
end

% Hypothesis for colocation=-1
Hn = kron(ones(2),H_);
Hn(logical(kron([1,1;0,1],ones(6)))) = NaN;
Hn(logical(kron([0,0;1,0],eye(6)))) = NaN;

Haa = kron(ones(2),H_);
Haa(logical(kron([0,1;1,1],ones(6)))) = NaN;
Haa(triu(true(12))) = NaN;

Hbb = kron(ones(2),H_);
Hbb(logical(kron([1,1;1,0],ones(6)))) = NaN;
Hbb(triu(true(12))) = NaN;

Hp = sum(cat(3,Haa,Hbb),3,'omitmissing');
Hp(isnan(Haa) & isnan(Hbb)) = NaN;

% Make the selectors for the real data RSM
S.n = ~isnan(Hn);
S.aa = ~isnan(Haa);
S.bb = ~isnan(Hbb);
S.p = ~isnan(Hp);
S.all = S.n | S.p;

% Load in predicted visual similarity from denseNet-169
V = VisSim(imgPerm,imgPerm);
V = kron(ones(2),V);

% Construct the design matrix
X = [...
    nanzscore(Hn(S.all)),...
    nanzscore(Hp(S.all)),...
    nanzscore(V(S.all))];

% Replace N/A values with mean
X(isnan(X)) = 0;

% Invert for regression
IX = pinv(X);

% Get the 12x12 RSM for the searchlight image
R = corr(M);

% Zscore within aa/bb/ba pairs
R(S.n) = nanzscore(R(S.n));
R(S.aa) = nanzscore(R(S.aa));
R(S.bb) = nanzscore(R(S.bb));
y = R(S.all);

% Do regression and Fischer-transform the results
z = atanh(IX*y);
return