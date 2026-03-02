function [] = z16_Norm2MNI_Epi(G ,fwhm, voxSize)
%Z16_NORM2MNI_EPI  Normalise EPI images to MNI space using DARTEL
%
%   Z16_NORM2MNI_EPI(G) normalises all EPI NIfTI images for the subject
%   group G into MNI space using precomputed DARTEL flowfields and a
%   group-specific DARTEL template. Images are written at EPI resolution
%   with default smoothing applied.
%
%   Z16_NORM2MNI_EPI(G, FWHM) additionally specifies the Gaussian smoothing
%   kernel applied during normalisation. FWHM must be a 1×3 vector in mm
%   (e.g. [5 5 5]). Use [0 0 0] for no smoothing.
%
%   Z16_NORM2MNI_EPI(G, FWHM, VOXSIZE) additionally specifies the voxel
%   size of the output images in mm (e.g. [2 2 2]). NaN values may be used
%   to inherit the corresponding resolution from the DARTEL template.
%
%   INPUTS
%   ------
%   G : char | string
%       Group identifier (e.g. 'G0', 'G1'), used to locate group-specific
%       DARTEL templates and subject flowfields.
%
%   FWHM : 1×3 numeric vector, optional
%       Full-width at half-maximum (mm) of the Gaussian smoothing kernel
%       applied during normalisation. Default = [5 5 5].
%
%   VOXSIZE : 1×3 numeric vector, optional
%       Output voxel size in mm. Values must be positive or NaN.
%       Default = [2 2 2] (EPI resolution).
%
%   DIRECTORY ASSUMPTIONS
%   ---------------------
%   The function assumes the following directory structure:
%
%     Subject data:
%       ../../Data/<subjectID>/EPI/2_Temporal/R*/ *.nii
%
%     Subject-specific DARTEL flowfields:
%       ../../Data/<subjectID>/Structural/<G>/u_*.nii
%
%     Group DARTEL template:
%       ../../Data/_Group/<G>/Structural/DARTEL_Templates/*_Template_6.nii
%
%   Normalised images are moved to:
%       ../../Data/<subjectID>/EPI/<G>/kXX/R*/
%   where XX corresponds to the mean FWHM smoothing value.
%
%   BEHAVIOUR
%   ---------
%   - Uses SPM DARTEL MNI normalisation (spm.tools.dartel.mni_norm)
%   - Processes all EPI volumes across all runs for each subject
%   - Bounding box is left unspecified (SPM defaults)
%   - Jacobian modulation is disabled (preserve = 0)
%   - Output files are prefixed with:
%       'w'  if FWHM = [0 0 0]
%       'sw' otherwise
%   - Output files are moved into run-specific target directories
%
%   OUTPUTS
%   -------
%   None. This function is called for its side effects (SPM batch
%   execution, file creation, and file movement).
%
%   DEPENDENCIES
%   ------------
%   - SPM12 (DARTEL toolbox)
%   - getSubjectIds
%
%   NOTES
%   -----
%   - Existing files matching 'w*.nii' or 'sw*.nii' may be overwritten
%     or moved.
%   - This function assumes that DARTEL templates and flowfields have
%     already been generated.
%
%   See also: GETSUBJECTIDS, SPM_JOBMAN

%% Handle smoothing kernel input
if nargin < 3 || isempty(fwhm)
    fwhm = [5 5 5];
    warning([ ...
        'No smoothing kernel specified.%c', ...
        'Defaulting to FWHM = [5 5 5] mm.'],10);
end

% Validate kernel
if ~isnumeric(fwhm) || numel(fwhm) ~= 3 || any(fwhm < 0)
    error('FWHM must be a 1×3 numeric vector of non-negative values (in mm).');
end
fwhm = double(fwhm(:)');  % enforce row vector

%% Handle voxel size input
if nargin < 4 || isempty(voxSize)
    voxSize = [2, 2, 2];
    warning([ ...
        'No voxel size specified.%c', ...
        'Defaulting to VOXSIZE = [2 2 2] (EPI resolution).'],10);
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

%% Set the EPI paths
paths.epi = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.epi{iSubject} = [...
        paths.subject{iSubject},...
        filesep,'EPI'];
end

%% Set the runLists and Temporal + sources paths
runLists = cell(size(subjectIds));
paths.temporal = cell(size(subjectIds));
paths.sources = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.temporal{iSubject} = [...
        paths.epi{iSubject},...
        filesep,'2_Temporal'];
    dirList = dir(paths.temporal{iSubject});
    s = [dirList.isdir]' & ...
        cellfun(@(s)strcmp(s(1),'R'),{dirList.name}');
    runLists{iSubject} = {dirList(s).name}';
    paths.sources{iSubject} = cell(size(runLists{iSubject}));
    for iRun = 1:numel(runLists{iSubject})
        paths.sources{iSubject}{iRun} = [...
            paths.temporal{iSubject},...
            filesep,runLists{iSubject}{iRun}];
    end
end

%% Set the G, K, and target paths
paths.G = cell(size(subjectIds));
paths.k = cell(size(subjectIds));
paths.targets = cell(size(subjectIds));
for iSubject = 1:nSubjects
    paths.G{iSubject} = [...
        paths.epi{iSubject},...
        filesep,G];
    paths.k{iSubject} = [...
        paths.G{iSubject},...
        filesep,sprintf('k%02d',round(mean(fwhm)))];
    paths.targets{iSubject} = cell(size(runLists{iSubject}));
    for iRun = 1:numel(runLists{iSubject})
        paths.targets{iSubject}{iRun} = [...
            paths.k{iSubject},...
            filesep,runLists{iSubject}{iRun}];
        if ~exist(paths.targets{iSubject}{iRun},'dir')
            mkdir(paths.targets{iSubject}{iRun});
        end
    end
end

%% Loop through participants to get fns.sources
fns.sources = cell(size(subjectIds));
for iSubject = 1:nSubjects
    for iRun = 1:numel(paths.sources{iSubject})
        dirList = dir([paths.sources{iSubject}{iRun},filesep,'*.nii']);
        [~,ord] = sort({dirList.name});
        dirList = dirList(ord);
        fnCount = numel(dirList);
        for ii = 1:fnCount
            fns.sources{iSubject} = [fns.sources{iSubject};...
                {[dirList(ii).folder,...
                filesep,dirList(ii).name]}];
        end
    end
end

%% Run DARTEL
SpmBatch = {};
SpmBatch{1}.spm.tools.dartel.mni_norm.template = {fns.template};
for iSubj = 1:nSubjects
    SpmBatch{1}.spm.tools.dartel.mni_norm.data.subj(iSubj).flowfield = ...
        fns.flowfields(iSubj);
    SpmBatch{1}.spm.tools.dartel.mni_norm.data.subj(iSubj).images = ...
        fns.sources{iSubj};
end
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = voxSize;
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = fwhm;
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Move the output files
for iSubject = 1:nSubjects
    for iRun = 1:numel(paths.sources{iSubject})
        if isequal(fwhm,[0,0,0])
            movefile([paths.sources{iSubject}{iRun},filesep,'w*.nii'],...
                paths.targets{iSubject}{iRun});
        else
            movefile([paths.sources{iSubject}{iRun},filesep,'sw*.nii'],...
                paths.targets{iSubject}{iRun});
        end
    end
end
return