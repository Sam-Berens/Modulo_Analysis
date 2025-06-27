function [] = z07_Coregister(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
subjDir = [dataDir,filesep,subjectId];

%% Get the MP-RAGE
t1Dir = dir([subjDir,'/Structural/m_*.nii']);
FileName_T1 = [t1Dir.folder,filesep,t1Dir.name];

%% Get the mean EPI image
muDir = dir([subjDir,'/EPI/meanu_*.nii']);
FileName_Mu = [muDir.folder,filesep,muDir.name];

%% Get the file paths for the EPI data
epiDir = [subjDir,filesep,'EPI/1_Realigned'];

FilePath_EPIs = cell(0,1);
for iRun = 1:5
    sourceDir = [epiDir,filesep,'R',int2str(iRun)];
    dirList = dir(sprintf('%s%su_*.nii',sourceDir,filesep));
    FilePath_EPIs = [FilePath_EPIs;
        cellfun(...
        @(x,y)[x,filesep,y],...
        {dirList.folder}',...
        {dirList.name}',...
        'UniformOutput',false)]; %#ok<AGROW>
end

%% Create and execute the SPM batch
SpmBatch = {};
SpmBatch{1}.spm.spatial.coreg.estimate.ref = {FileName_T1};
SpmBatch{1}.spm.spatial.coreg.estimate.source = {FileName_Mu};
SpmBatch{1}.spm.spatial.coreg.estimate.other = FilePath_EPIs;
SpmBatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
SpmBatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2 1];
SpmBatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
SpmBatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

return