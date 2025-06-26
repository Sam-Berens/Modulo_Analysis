function [] = z05_RealignUnwarp(subjectId)

preprocDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);
cd(subjectId);

%% Get FilePath_VDMs
cd Fieldmap;

% Get filepaths for first volume of each functional run
FilePath_VDMs = cell(5,1);
for iRun = 1:5
    dirList = dir(sprintf('vdm5__*_FmROI_R%i.nii',iRun));
    FilePath_VDMs{iRun} = [dirList.folder, filesep, dirList.name];
end

cd ..;

%% Get the file paths for the EPI data
cd EPI/0_Raw;

FilePath_EPIs = cell(5,1);
for iRun = 1:5
    cd(['R',int2str(iRun)]);
    dirList = dir('_*.nii');
    cd ..;
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
    cd(['R',int2str(iRun)]);
    mkdir(['../../1_Realigned/R',int2str(iRun)]);
    movefile(...
        'u_*.nii',...
        ['../../1_Realigned/R',int2str(iRun),filesep]);
    cd ..;
end

%% Cd out of subject and then out of data
cd(preprocDir);

return