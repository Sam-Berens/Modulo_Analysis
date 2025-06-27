function [] = z05_RealignUnwarp(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
subjDir = [dataDir,filesep,subjectId];

%% Get FilePath_VDMs
fmDir = [subjDir,filesep,'Fieldmap'];

% Get filepaths for first volume of each functional run
FilePath_VDMs = cell(5,1);
for iRun = 1:5
    dirList = dir(sprintf('%s%svdm5__*_FmORI_R%i.nii',fmDir,filesep,iRun));
    FilePath_VDMs{iRun} = [dirList.folder, filesep, dirList.name];
end

%% Get the file paths for the EPI data
epiRawDir = [subjDir,filesep,'EPI/0_Raw'];

FilePath_EPIs = cell(5,1);
for iRun = 1:5
    sourceDir = [epiRawDir,filesep,'R',int2str(iRun)];
    dirList = dir(sprintf('%s%s_*.nii',sourceDir,filesep));
    FilePath_EPIs{iRun} = cellfun(...
        @(x,y)[x,filesep,y],...
        {dirList.folder}',...
        {dirList.name}',...
        'UniformOutput',false);
end

%% Create and execute the SPM batch
for iRun = 1:5
    SpmBatch{1}.spm.spatial.realignunwarp.data(iRun).scans = ...
        FilePath_EPIs{iRun};
    SpmBatch{1}.spm.spatial.realignunwarp.data(iRun).pmscan = ...
        FilePath_VDMs(iRun);
end

SpmBatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.sep = 1;
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 1;
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 1;
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 1 0];
SpmBatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 2;
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
SpmBatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.jm = 1;
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 1 0];
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
SpmBatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Move new files
for iRun = 1:5
    sourceDir = [epiRawDir,filesep,'R',int2str(iRun)];
    targetDir = [subjDir,'/EPI/1_Realigned/R',int2str(iRun)];
    mkdir(targetDir);
    movefile(...
        [sourceDir,'/u_*.nii'],...
        targetDir);
    movefile(...
        [sourceDir,'/_*_uw.mat'],...
        targetDir);
    movefile([sourceDir,'/rp__*.txt'],[subjDir,'/EPI']);
    if iRun == 1
        movefile([sourceDir,'/meanu_*.nii'],[subjDir,'/EPI']);
    end
end

return