function [] = norm2MNI(G,imgPath)
%NORM2MNI  Normalise subject images to MNI space using DARTEL
%
%   NORM2MNI(G, IMGPATH) normalises subject-level NIfTI images into
%   MNI space using precomputed DARTEL flowfields and a group-specific
%   DARTEL template. The function constructs and runs an SPM DARTEL
%   normalisation batch, then moves the resulting normalised images
%   into their corresponding target directories.
%
%   INPUTS
%   ------
%   G : char | string
%       Group identifier (e.g. 'G0', 'G1'). This is used to:
%         - obtain subject IDs via GETSUBJECTIDS
%         - locate group-specific DARTEL templates
%         - define group-level directory paths
%
%   IMGPATH : char | string
%       Relative path (from each subject directory) to the images
%       that should be normalised to MNI space. These images are
%       assumed to already be in DARTEL (template) space.
%
%   DIRECTORY ASSUMPTIONS
%   ---------------------
%   The function assumes the following directory structure:
%
%     Source images to be normalised:
%      ../../Data/<subjectID>/<IMGPATH>/*.nii
%
%     Final DARTEL template:
%      ../../Data/_Group/<G>/Structural/DARTEL_Templates/*_Template_6.nii
%
%     Subject-specific flowfields:
%      ../../Data/<subjectID>/Structural/<G>/u_*.nii
%
%   BEHAVIOUR
%   ---------
%   - Uses SPM's DARTEL MNI normalisation (spm.tools.dartel.mni_norm)
%   - Voxel size and bounding box are left as NaN (SPM defaults)
%   - Jacobian modulation is disabled (preserve = 0)
%   - No smoothing is applied (fwhm = [0 0 0])
%   - Output images are prefixed with 'w'
%   - Normalised images are moved into subject-specific target folders
%
%   OUTPUTS
%   -------
%   None. This function is called for its side effects (file generation
%   and file movement).
%
%   DEPENDENCIES
%   ------------
%   - SPM (DARTEL toolbox)
%   - GETSUBJECTIDS
%
%   NOTES
%   -----
%   - Existing files matching 'w*.nii' in the source directories may be
%     overwritten or moved.
%   - This function assumes that DARTEL templates and flowfields have
%     already been created.
%
%   See also: GETSUBJECTIDS, SPM_JOBMAN

%% Get the list of subjectIds
subjectIds = getSubjectIds(G);
nSubjects = numel(subjectIds);

%% Set fns.template
paths.data = '../../Data';
paths.templates = [...
    paths.data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];
fns.template = dir([paths.templates,filesep,'*_Template_6.nii']);
if numel(fns.template) ~= 1
    error('Unexpected number of templates!');
end
fns.template = [fns.template.folder,filesep,fns.template.name];

%% Get the subject dirs
paths.subject = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.subject{iSubject} = [...
        paths.data,...
        filesep,char(subjectIds(iSubject))];
end

%% Get the names of the flowfields
paths.flowfields = cell(size(subjectIds));
fns.flowfields = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.flowfields{iSubject} = [...
        paths.subject{iSubject},...
        filesep,'Structural',...
        filesep,G];
    dirList = dir(...
        [paths.flowfields{iSubject},filesep,'u_*.nii']);
    if numel(dirList) ~= 1
        error('Unexpected number of flowfields for %s!',...
            char(subjectIds{iSubject}));
    end
    fns.flowfields{iSubject} = [dirList.folder,filesep,dirList.name];
end

%% Get the input/output paths
paths.sources = cell(size(subjectIds));
paths.targets = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.sources{iSubject} = [...
        paths.subject{iSubject},...
        filesep,imgPath];
    paths.targets{iSubject} = [...
        paths.sources{iSubject},...
        filesep,G];
    if ~exist(paths.targets{iSubject},'dir')
        mkdir(paths.targets{iSubject});
    end
end

%% Loop through participants to get fns.sources
for iSubject = 1:nSubjects
    dirList = dir([paths.sources{iSubject},filesep,'*.nii']);
    [~,ord] = sort({dirList.name});
    dirList = dirList(ord);
    fnCount = numel(dirList);
    if iSubject == 1
        fnCount0 = fnCount;
        fns.sources = cell(1,fnCount);
        for ii = 1:fnCount
            fns.sources{ii} = cell(size(subjectIds));
        end
    elseif fnCount ~= fnCount0
        error('Inconsistent number of source files!');
    end
    for ii = 1:fnCount
        fns.sources{ii}{iSubject} = [dirList(ii).folder,...
            filesep,dirList(ii).name];
    end
end

%% Run DARTEL
SpmBatch = {};
SpmBatch{1}.spm.tools.dartel.mni_norm.template = {fns.template};
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = fns.flowfields;
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = fns.sources;
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Move the output files
for iSubject = 1:nSubjects
    movefile([paths.sources{iSubject},filesep,'w*.nii'],...
        paths.targets{iSubject});
end
return