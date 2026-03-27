function [DataTables05,mniCoords] = getDataTable05_MNI(mniCoords)
%% Takes 3 x n vector (x,y,z) of MNI coords
% Returns a cell array of tables which is n(Sample locations) x 1,
% and the rounded MNI coords sampled from

%TO DO fix names as inconsistent and is probs also a bit inefficient
if size(mniCoords,1)~=3 
    mniCoords = mniCoords';
elseif numel(size(mniCoords)) ~= 2
    error('input must be 2d')
elseif all(size(mniCoords) == [3,3])
    disp(['Warning: uncheckable dimensions,'...
        'make sure the rows of your input are xyz'])
end 

G = 'G1';
dirs.Data = '../../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl5 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl05a');

%make a list of the filenames for the volumes to load (in the same order as
%the datatable will be in)
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
dSubjectIds = repelem(subjectIds,2,1);
%colocation=-1 goes first
zTempFn = categorical({'wzTemplate_colocation=-1.nii';...
    'wzTemplate_colocation=+1.nii'});
zTempFn = repmat(zTempFn,[nSubs,1]);
zTempFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),dSubjectIds,zTempFn,...
    'UniformOutput',false);
V = spm_vol(char(zTempFn));
ztM = spm_read_vols(V);

% get voxel coords from mni coords 
%(affine should be the same for all normalised images)
affine = V(1).mat;
dims = V(1).dim';

%round mni coordinates to nearest real voxel
nSamps = size(mniCoords,2);
mniCoords = mat2cell(mniCoords,3,ones(1,nSamps));
mniCoords = cellfun(@(x) spm_XYZreg('RoundCoords',x,affine,dims),...
    mniCoords,'UniformOutput',false);
%unpack them
mniCoords =  [cell2mat(mniCoords);ones(1,nSamps)];
%Find the corresponding voxel locations
vxIdxs = affine \ mniCoords;
vxIdxs = vxIdxs(1:3,:)';
%reformat rounded mniCoords to be returned (so you know exactly where
%you've sampled from)
mniCoords = mniCoords(1:3,:);
%sample the 4d image across subjects and coloc conditions
ztS = getSamps(ztM,vxIdxs);

visFn = categorical({'wzTemplate_visSim.nii'});
visFn = repmat(visFn,[nSubs,1]);
visFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),subjectIds,visFn,...
    'UniformOutput',false);
V = spm_vol(char(visFn));
visM = spm_read_vols(V);
%sample the 4d image across subjects 
visS = getSamps(visM,vxIdxs);

hoodFn = categorical({'whoodSize.nii'});
hoodFn = repmat(hoodFn,[nSubs,1]);
hoodFn = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),subjectIds,hoodFn,...
    'UniformOutput',false);
V = spm_vol(char(hoodFn));
hsM = spm_read_vols(V);
%sample the 4d image across subjects 
hsS = getSamps(hsM,vxIdxs);
%Build a cell array of tables (nSample_locations x 1)
DataTables05 = compileTbls(ztS,visS,hsS,G); 

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


function [DataTables05] = compileTbls(ztS,visS,hsS,G)
nSamps = size(ztS,1);
%% Get pNonc
pNonc = get_pNonc(G);
pNonc.cpNonc = pNonc.pNonc - mean(pNonc.pNonc);
subjectIds = pNonc.subjectId;

%% ROI loop
DataTables05 = cell(size(nSamps));
for iSamp = 1:nSamps
    %select out the data for the current sample location
    czT = ztS(iSamp,:);
    cvS = visS(iSamp,:);
    chS = hsS(iSamp,:);
    %construct table 
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

function [PatternSim] = buildTbl(subjectList,ztS,visS,hsS)
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
    zTemplate(iIn) = ztS(iIn); %both coloc conditions are already in here
    zVisual(iIn) = visS(iSubject);
    hoodSize(iIn) = hsS(iSubject);
    subjectId{iIn} = char(cSubjectId);
    % Populate colocation, zTemplate, visSim and hoodSize
    iIn = iIn + 1;
    colocation(iIn) = +1;
    zTemplate(iIn) = ztS(iIn);
    zVisual(iIn) = visS(iSubject);
    hoodSize(iIn) = hsS(iSubject);
    subjectId{iIn} = char(cSubjectId);
end
subjectId = categorical(subjectId);

%% Make the data table
PatternSim = table(subjectId,colocation,zTemplate,zVisual,hoodSize);
return
