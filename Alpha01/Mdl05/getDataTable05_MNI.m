function [DataTables05] = getDataTable05_MNI(mniCoords)
%% Takes MNI coordintes as a 3 x n vector (x,y,z)
%% Returns a cell array of tables which is n(Sample locations) x 1

%TO DO fix names as inconsistent and is probs also a bit inefficient
if size(mniCoords,1)~=3 
    mniCoords = mniCoords';
elseif numel(size(mniCoords)) ~= 2
    error('input must be 2d')
end 

G = 'G1';

dirs.Data = '../../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl5 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl05a');

%get 'scan' names in correct order for R rows
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
dSubjectIds = repelem(subjectIds,2,1);
%colocation=-1 goes first
srchLghtFn = categorical({'wzTemplate_colocation=-1.nii';...
    'wzTemplate_colocation=+1.nii'});
srchLghtFn = repmat(srchLghtFn,[nSubs,1]);
srchLghtFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),dSubjectIds,srchLghtFn,...
    'UniformOutput',false);
v = spm_vol(char(srchLghtFn));
zTemplate = spm_read_vols(v);

% get voxel coords from mni coords 
%(affine should be the same for all normalised images)
affine = v(1).mat;
dims = v(1).dim';
% mni = [mniCoords; 1];

%round mni coordinates to nearest real voxel
nSamps = size(mniCoords,2);
mniCoords = mat2cell(mniCoords,3,ones(1,nSamps));
mniCoords = cellfun(@(x) spm_XYZreg('RoundCoords',x,affine,dims),...
    mniCoords,'UniformOutput',false);
%unpack them
mniCoords =  [cell2mat(mniCoords);ones(1,nSamps)];
vxIdxs = affine \ mniCoords;
vxIdxs = vxIdxs(1:3,:)';
zTs = getSamps(zTemplate,vxIdxs);


visFn = categorical({'wzTemplate_visSim.nii'});
visFn = repmat(visFn,[nSubs,1]);
visFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),subjectIds,visFn,...
    'UniformOutput',false);
v = spm_vol(char(visFn));
visSim = spm_read_vols(v);

vSs = getSamps(visSim,vxIdxs);

hoodFn = categorical({'whoodSize.nii'});
hoodFn = repmat(hoodFn,[nSubs,1]);
hoodFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),subjectIds,hoodFn,...
    'UniformOutput',false);
v = spm_vol(char(hoodFn));
hoodSize = spm_read_vols(v);

hSs = getSamps(hoodSize,vxIdxs);

DataTables05 = compileTbls(zTs,vSs,hSs,G); 

return


function [Samps] = getSamps(M,vxIdxs)
nSamps = size(vxIdxs,1);
nVox = numel(M(:,:,:,1));
nPages = size(M,4);
Samps = nan(nSamps,nPages);
for ii=1:nSamps
    vxIdx = num2cell(vxIdxs(ii,:));
    lVxIdx = sub2ind(size(M,[1,2,3]),vxIdx{:});
    idxRead = nVox*(0:(nPages-1)) + lVxIdx;
    Samps(ii,:) = M(idxRead);
end
return


function [DataTables05] = compileTbls(zTemplate,viSim,hoodSize,G)
nSamps = size(zTemplate,1);
%% Get pNonc
pNonc = get_pNonc(G);
pNonc.cpNonc = pNonc.pNonc - mean(pNonc.pNonc);
subjectIds = pNonc.subjectId;

%% ROI loop
DataTables05 = cell(size(nSamps));
for iSamp = 1:nSamps   %change back to parfor?
    czT = zTemplate(iSamp,:);
    cvS = viSim(iSamp,:);
    chS = hoodSize(iSamp,:);
    DataTables05{iSamp} = buildTbl(subjectIds,czT,cvS,chS);
end

%% Do the outer join with pNonc
for iSamp = 1:nSamps
    DT = DataTables05{iSamp};
    DT = outerjoin(pNonc,DT);
    % Rename subjectId after outerjoin
    DT.subjectId_DT = [];
    DT.Properties.VariableNames{1} = 'subjectId';
    DataTables05{iSamp} = DT;
end

return

function [PatternSim] = buildTbl(subjectList,zT,vS,hS)
nSubjects = numel(subjectList);

%% Preallocate
subjectId = cell(nSubjects*2,1);
colocation = nan(nSubjects*2,1);
zTemplate = nan(nSubjects*2,1);
zVisual = nan(nSubjects*2,1);
hoodSize = nan(nSubjects*2,1); %change name to hood size or keep same?

%% Loopy loop
iIn = 0;
for iSubject = 1:nSubjects
    cSubjectId = subjectList(iSubject);
    iIn = iIn + 1;
    % Populate colocation, zTemplate, visSim and hoodSize
    colocation(iIn) = -1;
    zTemplate(iIn) = zT(iIn); %both coloc conditions are already in here
    zVisual(iIn) = vS(iSubject);
    hoodSize(iIn) = hS(iSubject);
    subjectId{iIn} = char(cSubjectId);
    % Populate colocation, zTemplate, visSim and hoodSize
    iIn = iIn + 1;
    colocation(iIn) = +1;
    zTemplate(iIn) = zT(iIn);
    zVisual(iIn) = vS(iSubject);
    hoodSize(iIn) = hS(iSubject);
    subjectId{iIn} = char(cSubjectId);
end
subjectId = categorical(subjectId);

%% Make the data table
PatternSim = table(subjectId,colocation,zTemplate,zVisual,hoodSize);
return
