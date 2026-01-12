function [] = z13_InverseWarpRois(G,roisToWarp)
% Inverse warp ROIs from MNI space to native space (for each participant).
% This procedure uses DARTEL tools and so requires a DARTEL group template,
% along with flowfield images per participant. This routine involves three
% distinct steps:
%   1) Performing an affine transformation of the ROI images so they are
%      in registration with the DARTEL group template (by manipulating the
%      image headers alone).
%   2) Calling the "tools.dartel.crt_iwarped" utility in SPM to perform the
%      non-linear warp.
%   3) Moving the inverse-warped files into a dedicated NativeRoi directory
%      (see code for details).
% Note that the resulting inverse normalised images have the same
% dimensions as the flowfields.
%

%% Input checks
if nargin < 2
    roisToWarp = {};
    % This will be populated will all the native ROIs that match the 
    % following filename schema:
    %    ../../Data/_Group/MniRois/MniRois/_*.nii
elseif nargin < 1
    error('You must provide a group ID (string).');
end
if ischar(roisToWarp)
    roisToWarp = {roisToWarp};
elseif ~iscell(roisToWarp)
    error(...
    'roisToWarp input must be either a char or a cell array of chars.');
end

%% Get subjectsIds
subjectIds = getSubjectIds(G);

%% Set the paths structure
paths.data = ['..',filesep,'..',filesep,'Data'];
paths.source = [paths.data,...
    filesep,'_Group',...
    filesep,'MniRois'];
paths.templates = [paths.data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];
paths.groupRois = [paths.data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'GroupRois'];

% Loop through the participants to get the locations of the flowfields ...
% ... and the target location of the warped native ROIs.
paths.flowfields = cell(size(subjectIds));
paths.targets = cell(size(subjectIds));
for iSubject = 1:numel(subjectIds)
    paths.flowfields{iSubject} = [paths.data,...
        filesep,subjectIds{iSubject},...
        filesep,'Structural',...
        filesep,G];
    paths.targets{iSubject} = [paths.flowfields{iSubject},...
        filesep,'NativeRois',...
        filesep,'w'];
end

%% Set the fns structure
fns.template = dir([paths.templates,filesep,'*_Template_6.nii']);
fns.template = [fns.template.folder,filesep,fns.template.name];
fns.group2Mni = dir([paths.templates,filesep,'*_Template_6_2mni.mat']);
fns.group2Mni = [fns.group2Mni.folder,filesep,fns.group2Mni.name];
if isempty(roisToWarp)
    dirList = dir([paths.source,filesep,'_*.nii']);
    roisToWarp = {dirList.name}';
end
fns.roisToWarp = cellfun(@(s)[paths.source,filesep,s],...
    roisToWarp,'UniformOutput',false);
fns.groupRois = cellfun(@(s)[paths.groupRois,filesep,s],...
    roisToWarp,'UniformOutput',false);

% Loop through participants to get the names of the flowfields
fns.flowfields = cell(size(subjectIds));
for iSubject = 1:numel(subjectIds)
    dirList = dir(...
        [paths.flowfields{iSubject},filesep,'u_*.nii']);
    fns.flowfields{iSubject} = [dirList.folder,filesep,dirList.name];
end

%% 1) Affine transform the MNI ROIs
% M1 is the direction cosine matrix of the ROI in MNI space:
%   VxROI -> MNI
% M2 is the *rotation matrix* that moves the DARTEL template to MNI space:
%   GroupSpace -> MNI (we want the inverse of this);
% M3 is the direction cosine matrix of the DARTEL template in group space:
%   VxTemplate -> GroupSpace
% M4 is the direction cosine matrix of the ROI in group space:
%   VxROI -> GroupSpace
%   M4 = M3*inv(M2)*M1;

% First we load M2 and M3
M2 = load(fns.group2Mni);
M2 = M2.mni.affine;
iM2 = inv(M2);
M3 = spm_vol(fns.template);
M3 = M3.mat;

% Now loop through ROIs to transform .mat and write
if ~exist(paths.groupRois,'dir')
    mkdir(paths.groupRois);
end
for iRoi = 1:numel(roisToWarp)
    V1 = spm_vol(fns.roisToWarp{iRoi});
    Y1 = spm_read_vols(V1);
    V1.fname = fns.groupRois{iRoi};
    V1.descrip = ['Affined2Group ',V1.descrip];
    M1 = V1.mat;
    M4 = M3*iM2*M1; %#ok<MINV>
    V1.mat = M4;
    spm_write_vol(V1,Y1);
end

%% Create and execute SPM batch
spmBatch = {};
spmBatch{1}.spm.tools.dartel.crt_iwarped.flowfields = fns.flowfields;
spmBatch{1}.spm.tools.dartel.crt_iwarped.images = fns.groupRois;
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
        [paths.flowfields{iSubject},filesep,'w_*.nii'],...
        paths.targets{iSubject});
end

return