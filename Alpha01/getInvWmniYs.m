function [yIms] = getInvWmniYs()
matName = 'invW_mmY.mat';
if exist(matName,"file")
    strct = load(matName);
    yIms = strct.yIms;
    return
end

G = 'G1';
dirs.data = '../../Data';
dirs.template = fullfile(dirs.data,'_Group',G,'Structural','DARTEL_Templates');

g.V = spm_vol(fullfile(dirs.template,'G1_Template_6.nii'));
dims = g.V.dim;
strct = load(fullfile(dirs.template,'G1_Template_6_2mni.mat'));
affine = strct.mni.affine;
y.M = nan(dims); %I think this means this voxels matrix is technically in group space because of the dimensions
[X,Y,Z] = ndgrid(1:dims(1),1:dims(2),1:dims(3));
voxCoords = [X(:) Y(:) Z(:) ones(numel(X),1)]';
mmXYZ = affine * voxCoords;
mmY = mmXYZ(2,:)';
y.M(isnan(y.M)) = mmY; %but contains the y coordinate of each voxel in mni space
y.V = g.V(1);
y.V.fname = fullfile(dirs.template,'_mmY.nii');
y.V.descrip = 'MNI-mm y coordinates for group-space voxel matrix';
y.V.dt = [16,0];
spm_write_vol(y.V,y.M);
inverseWarpY(G,y.V.fname);

%load new images in
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
yIms = struct('V',[],'M',[]);
for iSub = 1:nSubs
    cid = char(subjectIds(iSub));
    ydir = fullfile(dirs.data,cid,'Structural',G,...
        'NativeRois','w');
    ydir = dir([ydir,filesep,'w_mmY*.nii']);
    yFname = fullfile(ydir.folder,ydir.name);
    yIms(iSub).V = spm_vol(yFname);
    yIms(iSub).M = spm_read_vols(y.V);
end 

save(matName,'yIms');

return




function [] = inverseWarpY(G,fname)
% Inverse warp the ycoord image from group space to native space (for each participant).
% This procedure uses DARTEL tools and so requires a DARTEL group template,
% along with flowfield images per participant. This routine involves three
% distinct steps:
%   1) Calling the "tools.dartel.crt_iwarped" utility in SPM to perform the
%      non-linear warp.
%   2) Moving the inverse-warped files into a dedicated Native directory
%      (see code for details).
% Note that the resulting inverse normalised images have the same
% dimensions as the flowfields.
%

%% Input checks
if nargin < 1
    error('You must provide a group ID (string).');
end
if ischar(fname)
    fname = {fname};
elseif ~iscell(fname)
    error(...
    'fname input must be either a char or a cell array of chars.');
end

%% Get subjectsIds
subjectIds = getSubjectIds(G);

%% Set the paths structure
paths.data = ['..',filesep,'..',filesep,'Data'];
paths.templates = [paths.data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];

% Loop through the participants to get the locations of the flowfields ...
% ... and the target location of the warped native ROIs.
paths.flowfields = cell(size(subjectIds));
paths.targets = cell(size(subjectIds));
for iSubject = 1:numel(subjectIds)
    paths.flowfields{iSubject} = [paths.data,...
        filesep,char(subjectIds(iSubject)),...
        filesep,'Structural',...
        filesep,G];
    paths.targets{iSubject} = [paths.flowfields{iSubject},...
        filesep,'NativeRois',...
        filesep,'w'];
end

%% Set the fns structure
fns.template = dir([paths.templates,filesep,'*_Template_6.nii']);
fns.template = [fns.template.folder,filesep,fns.template.name];

% Loop through participants to get the names of the flowfields
fns.flowfields = cell(size(subjectIds));
for iSubject = 1:numel(subjectIds)
    dirList = dir(...
        [paths.flowfields{iSubject},filesep,'u_*.nii']);
    fns.flowfields{iSubject} = [dirList.folder,filesep,dirList.name];
end


%% Create and execute SPM batch
spmBatch = {};
spmBatch{1}.spm.tools.dartel.crt_iwarped.flowfields = fns.flowfields;
spmBatch{1}.spm.tools.dartel.crt_iwarped.images = fname;
spmBatch{1}.spm.tools.dartel.crt_iwarped.K = 9;
spmBatch{1}.spm.tools.dartel.crt_iwarped.interp = 7;
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

%% Loop through participants to move the inverse warped ROIs
for iSubject = 1:numel(subjectIds)
    if ~exist(paths.targets{iSubject},'dir')
        mkdir(paths.targets{iSubject});
    end
    movefile(...
        [paths.flowfields{iSubject},filesep,'w_mmY*.nii'],...
        paths.targets{iSubject});
end

return