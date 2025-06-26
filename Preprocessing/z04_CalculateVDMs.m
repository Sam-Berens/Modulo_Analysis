function [] = z04_CalculateVDMs(subjectId)

spmRoot = getSpmHome();
preprocDir = pwd;
% Cd into the data folder (from the pre-processingc folder)
cd(['..',filesep,'..',filesep,'Data']);
% Cd into subject dir
cd(subjectId);

%% Get FileNames from /Fieldmap
cd Fieldmap;

% Make FileName_ORI variable:
FileList = dir('*ORI.nii');
if sum(double(size(FileList) == [1,1]),2) ~= 2
    error('There appears to be multiple (%i) ORIs for %s',...
        max(size(FileList)),subjectId);
end
FilePath_ORI= sprintf('%s%s%s',pwd,filesep,FileList.name);

% Make FileName_FmMag variable:
FileList = dir('*FmMag.nii');
if sum(double(size(FileList) == [1,1]),2) ~= 2
    error('There appears to be multiple (%i) FmMags for %s',...
        max(size(FileList)),subjectId);
end
FilePath_FmMag= sprintf('%s%s%s',pwd,filesep,FileList.name);


cd ..;

%% Cd into the EPI/0_Raw folder
cd EPI/0_Raw;

% Get filepaths for first volume of each functional run
FilePath_EPIs = cell(5,1);
for iRun = 1:5
    cd(['R',int2str(iRun)]);
    dirList = dir('*001.nii');
    cd ..;
    if isempty(dirList)
        error('The first EPI volume of for %s, run %i is missing!',...
            subjectId,iRun);
    end
    FilePath_EPIs{iRun} = [dirList.folder, filesep, dirList.name];
end

%% Create and execute batch

% Select the off-resonance image (in Hz)
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap...
    .precalcfieldmap = {FilePath_ORI};

% Select the magnitude image (for co-registration to the EPI data)
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap...
    .magfieldmap = {FilePath_FmMag};

% Set a stand-in short and long echo time to appease the spm_job manager
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .et = [0,0];

% Set the mask brain option
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .maskbrain = 1;

% Set the blip direction of the EPI images, SPM defines P->A as positive
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .blipdir = -1;

% Set the total readout time in milliseconds (SPM definition)
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .tert = 69.68;

% Specify whether the fieldmaps are EPI based
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .epifm = 0;

% Specify whether to use Jacobian modulation
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .ajm = 0;

% Specify the unwrapping options
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .uflags.method = 'Mark3D';
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .uflags.fwhm = 10;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .uflags.pad = 0;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .uflags.ws = 1;

% Specify the segmentation options
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.template = {[spmRoot,'/toolbox/FieldMap/T1.nii']};
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.fwhm = 5;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.nerode = 2;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.ndilate = 4;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.thresh = 0.5;
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval...
    .mflags.reg = 0.02;

% Loop through to add EPI volumes for each
for iRun = 1:5
    SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(iRun)...
        .epi = FilePath_EPIs(iRun);
end

% Specify whether the to co-registerc the VDM to the first EPI in each run
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;

% Specify the run prefix
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'R';

% Specify whether to write a sample-unwarpped
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;

% QA options
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
SpmBatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Move the sample unwarpped
for iRun = 1:5
    movefile(sprintf('.%sR%i%su_*.nii',filesep,iRun,filesep),...
        './SampleUnwarp/');
end

%% Cd out of subject and then out of data
cd(preprocDir);

return

function [spmRoot] = getSpmHome()
matlabPath = path;
% Get all paths in the MATLAB path
allPaths = strsplit(matlabPath, pathsep)';
spmRoot = allPaths{cellfun(@(x) strcmp(x(end-4:end), 'SPM25'), allPaths)};
return