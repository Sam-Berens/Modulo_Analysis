function [] = z08_SliceTimeCor(subjectId)
% useful details:
% TR = 2200 ms
% TA = 2112.5 ms
% deadtime = 87.5 ms
% tHalfofMax_TA = 2112.5/2 = 1,056.25ms
% closest t for real slice = 1057.5ms (slice 33 and 66) 
% slice 33 time using interslice-gap of 66ms: 1056ms

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
subjDir = [dataDir,filesep,subjectId];

%% Get the file paths for the EPI data
epiFolder = [subjDir,filesep,'EPI',filesep,'1_Realigned'];
epiDir = dir([epiFolder,filesep,'R*']);
runNames = {epiDir.name}';
nRuns = numel(runNames);
FilePath_EPIs = cell(nRuns,1);
for iRun = 1:nRuns
    sourceDir = [epiFolder,filesep,runNames{iRun}]; 
    dirList = dir(sprintf('%s%su_*.nii',sourceDir,filesep));
    FilePath_EPIs{iRun} = cellfun(...
        @(x,y)[x,filesep,y],...
        {dirList.folder}',...
        {dirList.name}',...
        'UniformOutput',false);
end

%% Create and execute the SPM batch
for iRun = 1:nRuns
    spmBatch{1}.spm.temporal.st.scans{:,iRun} = ...
        FilePath_EPIs{iRun};
end

spmBatch{1}.spm.temporal.st.nslices = 0; % to appease the spm jobman, not used
spmBatch{1}.spm.temporal.st.tr = 2.2; %needs to be in seconds
spmBatch{1}.spm.temporal.st.ta = 0; % to appease the spm jobman, not used
spmBatch{1}.spm.temporal.st.so = getSliceT(); %SPM will know you are providing slice time for sliceOrder arg as long as ~isequal(1:nslices,sort(sliceOrder))
spmBatch{1}.spm.temporal.st.refslice = 1056; %this is ms, SPM knows its not an index because of above condition 
spmBatch{1}.spm.temporal.st.prefix = 'a';
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);


%% Move new files
for iRun = 1:nRuns
    sourceDir = [epiFolder,filesep,runNames{iRun}];
    targetDir = [subjDir,'/EPI/2_Temporal/',runNames{iRun}];
    mkdir(targetDir);
    movefile(...
        [sourceDir,filesep,'a_*.nii'],...
        targetDir);
    % movefile(...
    %     [sourceDir,'/_*_uw.mat'],...
    %     targetDir);
end

return
