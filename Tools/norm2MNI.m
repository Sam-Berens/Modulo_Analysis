function [] = norm2MNI(G,imgPath,fwhm, voxSize)
%NORM2MNI  Normalise subject images to MNI space using DARTEL
%
%   NORM2MNI(G, IMGPATH) normalises subject-level NIfTI images into
%   MNI space using precomputed DARTEL flowfields and a group-specific
%   DARTEL template.
%
%   NORM2MNI(G, IMGPATH, FWHM) additionally applies Gaussian smoothing
%   during normalisation, where FWHM is a 1×3 vector specifying the
%   kernel size in mm (e.g. [8 8 8]).
%
%   NORM2MNI(G, IMGPATH, FWHM, VOXSIZE) additionally controls the voxel
%   size of the output images. VOXSIZE is a 1×3 vector in mm (e.g.
%   [2 2 2]). Use [NaN NaN NaN] to inherit the resolution of the DARTEL
%   template.
%
%   If FWHM is omitted, it defaults to [0 0 0] (no smoothing) and a
%   warning is issued.
%
%   If VOXSIZE is omitted, it defaults to [NaN NaN NaN] (template
%   resolution) and a warning is issued.
%
%   INPUTS
%   ------
%   G : char | string
%       Group identifier (e.g. 'G0', 'G1').
%
%   IMGPATH : char | string
%       Relative path to the images to be normalised.
%
%   FWHM : 1×3 numeric vector, optional
%       Smoothing kernel size in mm. Default = [0 0 0].
%
%   VOXSIZE : 1×3 numeric vector, optional
%       Output voxel size in mm. Default = [NaN NaN NaN].
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
%   - Bounding box is left as NaN (SPM defaults)
%   - Jacobian modulation is disabled (preserve = 0)
%   - Output images are prefixed with 'w' or 'sw'
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
%   - Existing files matching 'w*.nii' or 'sw*.nii' in the source 
%     directories may be overwritten or moved.
%   - This function assumes that DARTEL templates and flowfields have
%     already been created.
%
%   See also: GETSUBJECTIDS, SPM_JOBMAN

%% Handle smoothing kernel input
if nargin < 3 || isempty(fwhm)
    fwhm = [0 0 0];
    warning([ ...
        'No smoothing kernel specified.%c', ...
        'Defaulting to FWHM = [0 0 0] mm (no smoothing).'],10);
end

% Validate kernel
if ~isnumeric(fwhm) || numel(fwhm) ~= 3 || any(fwhm < 0)
    error('FWHM must be a 1×3 numeric vector of non-negative values (in mm).');
end
fwhm = double(fwhm(:)');  % enforce row vector

%% Handle voxel size input
if nargin < 4 || isempty(voxSize)
    voxSize = [NaN NaN NaN];
    warning([ ...
        'No voxel size specified.%c', ...
        'Defaulting to VOXSIZE = [NaN NaN NaN] (template resolution).'],10);
end

% Validate voxel size
if ~isnumeric(voxSize) || numel(voxSize) ~= 3
    error('VOXSIZE must be a 1×3 numeric vector (in mm or NaN).');
end
if any(voxSize <= 0 & ~isnan(voxSize))
    error('VOXSIZE values must be positive or NaN.');
end
voxSize = double(voxSize(:)');  % enforce row vector

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
            char(subjectIds(iSubject)));
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
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = voxSize;
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = fwhm;
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Move the output files
for iSubject = 1:nSubjects
    if isequal(fwhm,[0,0,0])
        movefile([paths.sources{iSubject},filesep,'w*.nii'],...
            paths.targets{iSubject});
    else
        movefile([paths.sources{iSubject},filesep,'sw*.nii'],...
            paths.targets{iSubject});
    end
end
return