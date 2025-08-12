function [] = z06_Segment(subjectId)

preprocDir = pwd;
spmRoot = getSpmHome();

% Cd into the data folder (from the preprocessing folder)
cd(['..',filesep,'..',filesep,'Data']);
% Cd into subject dir
cd(subjectId);

%% Make FileName_Structural var:
cd Structural;
FileList = dir('*.nii');
if sum(double(size(FileList) == [1,1]),2) ~= 2
    error('There appears to be multiple (%i) structural scans for %s',...
        max(size(FileList)),subjectId);
end
FileName_Structural = sprintf('%s%s%s',pwd,filesep,FileList.name);
cd ..;

clearvars -except 'subjectId' 'FileName_Structural' 'spmRoot' 'preprocDir';

%% Create and execute SPM batch:
SpmBatch = {};
SpmBatch{1}.spm.spatial.preproc.channel.vols = {FileName_Structural};
SpmBatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; 
SpmBatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
SpmBatch{1}.spm.spatial.preproc.channel.write = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(1).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,1']};
SpmBatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
SpmBatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(2).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,2']};
SpmBatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
SpmBatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(3).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,3']};
SpmBatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
SpmBatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(4).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,4']};
SpmBatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
SpmBatch{1}.spm.spatial.preproc.tissue(4).native = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(5).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,5']};
SpmBatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
SpmBatch{1}.spm.spatial.preproc.tissue(5).native = [1 1];
SpmBatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(6).tpm = ...
    {[spmRoot,'\tpm\TPM.nii,6']};
SpmBatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
SpmBatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
SpmBatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
SpmBatch{1}.spm.spatial.preproc.warp.mrf = 1;
SpmBatch{1}.spm.spatial.preproc.warp.cleanup = 1;
SpmBatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
SpmBatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
SpmBatch{1}.spm.spatial.preproc.warp.fwhm = 0;
SpmBatch{1}.spm.spatial.preproc.warp.samp = 1;
SpmBatch{1}.spm.spatial.preproc.warp.write = [1 1];

% The voxel size (isotropic, in mm) of the written images.
SpmBatch{1}.spm.spatial.preproc.warp.vox = 1;

% The bounding box (in mm) of the volume which is to be written ...
% ... (relative to the anterior commissure).
SpmBatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN; NaN NaN NaN];

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Cd out of subject and then out of data
cd(preprocDir);

return

function [spmRoot] = getSpmHome()
matlabPath = path;
% Get all paths in the MATLAB path
allPaths = strsplit(matlabPath, pathsep)'; 
spmRoot = allPaths{cellfun(@(x) strcmp(x(end-4:end), 'SPM25'), allPaths)};
return