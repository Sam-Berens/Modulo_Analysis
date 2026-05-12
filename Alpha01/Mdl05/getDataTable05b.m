function [DataTables05b] = getDataTable05b(G,roiId,UseParfor)
% Builds datatable from which Mdl5 can be estimated on native roi data,
% either all schaefer parcells (default) or a specific roi
DataTables05b = struct;
if nargin < 1
    error('At least one argument (GroupId) required.');

elseif ~startsWith(G,'G')
    error('First argument must be GroupId.');

elseif nargin == 1
    fn = [G,'_DataTables05b.mat'];
    if exist(fn,'file')
        DataTables05b = temp.DataTables05;
        return
    end
    roiId = arrayfun(...
        @(ii)sprintf('N17P200_R%03d',ii),...
        (1:200)',...
        'UniformOutput',false);

elseif (nargin > 1) && ~iscell(roiId)
    roiId = {roiId};

elseif nargin > 3
    error('Maximum of 3 arguments: GroupId and (optional) roiPattern');
end

if ~exist('UseParfor','var')
    UseParfor = false;
end

%% Add .. to path
addpath('../');

%% Get pNonc
Pnonc = get_pNonc(G);
Pnonc.cpNonc = Pnonc.pNonc - mean(Pnonc.pNonc);
subjectIds = Pnonc.subjectId;

%% Loop through and only get the data datatables we need
s = ~ismember(roiId,fieldnames(DataTables05b));
toDo = roiId(s);
DT = cell(size(roiId));
if UseParfor
    parfor iRoi = 1:numel(toDo)
        DT{iRoi} = loopFunction(G,Pnonc,subjectIds,toDo{iRoi});
    end
else
    for iRoi = 1:numel(toDo)
        DT{iRoi} = loopFunction(G,Pnonc,subjectIds,toDo{iRoi});
    end
end

%% Unpack
for iRoi = 1:numel(toDo)
    DataTables05b.(toDo{iRoi}) = DT{iRoi};
end

%% Save
fn = [G,'_DataTables05b.mat'];
save(fn,'DataTables05b');

%% Remove .. from path
rmpath('../');

%% Select out the data tables we need
[~,idx] = ismember(roiId,fieldnames(DataTables05b));
DT = struct2cell(DataTables05b);
if isscalar(nonzeros(idx))
    DataTables05b = DT{nonzeros(idx)};
else
    DataTables05b = struct;
    for iRoi = 1:numel(roiId)
        DataTables05b.(roiId{iRoi}) = DT{iRoi};
    end

end

return

function [PatternSim] = getPatternSim(G,subjectList,roiId)
% Outputs is a table containing:
%  - subjectId: [nSubj*2,1]
%  - colocation: [nSubj*2,1]
%  - zTemplate: [nSubj*2,2]
%  - zVisual: [nSubj*2,2]
%  - pCover: [nSubj*2,1]
%  - voxSim: {nSubj*2,1}

nSubjects = numel(subjectList);

%% Preallocate
subjectId = cell(nSubjects*2,1);
colocation = nan(nSubjects*2,1);
zTemplate = nan(nSubjects*2,1);
zVisual = nan(nSubjects*2,1);
pCover = nan(nSubjects*2,1);
voxSim = cell(nSubjects*2,1);

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
    z = mdl5Func(Data,imgPerm);
    % z is a 3x1 vector of {-ve, +ve, vis} Fisher-transformed stats
    
    R = corr(Data);
    R(logical(kron(ones(2),eye(6)))) = NaN;
    Ps_n1 = R(7:end,1:6);
    A = R(1:6,1:6);
    B = R(7:end,7:end);
    A(triu(true(6),1)) = 0;
    B(triu(true(6),1)) = 0;
    Ps_p1 = A + B';

    % Populate colocation, zTemplate and pCover
    for cl = -1:2:1
        iIn = iIn + 1;
        subjectId{iIn} = char(cSubjectId);
        colocation(iIn) = cl;
        if cl == -1
            zTemplate(iIn) = z(1);
            voxSim{iIn} = Ps_n1;
        else
            zTemplate(iIn) = z(2);
            voxSim{iIn} = Ps_p1;
        end
        zVisual(iIn) = z(3);
        pCover(iIn) = cpCover;
    end

    waitbar(iSubject/nSubjects,fh);
end
close(fh);
subjectId = categorical(subjectId);

%% Make the data table
PatternSim = table(subjectId,colocation,zTemplate,zVisual,pCover,voxSim);
return

function [DT] = loopFunction(G,Pnonc,subjectIds,roiId)
DT = getPatternSim(G,subjectIds,roiId);
DT = outerjoin(Pnonc,DT,'MergeKeys',true);
return