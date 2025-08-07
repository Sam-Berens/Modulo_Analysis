function [] = z08_SliceTime(subjectId)
% Useful details:
% TR = 2200 ms
% TA = 2112.5 ms
% Deadtime = 87.5 ms
% 0.5*TA = 1056.25 ms
% Slice numbers corresponding to 0.5*TA = (33 and 66)
% Time of slices 33 and 66 = 1057.5 ms (variable?)
% 0.5*TA assuming an average 66 ms inter-slice gap: 1056 ms

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
subjDir = [dataDir,filesep,subjectId];

%% Get the file paths for the EPI data
epiFolder = [subjDir,filesep,'EPI'];
realinedFolder = [epiFolder,filesep,'1_Realigned'];
realinedDir = dir([realinedFolder,filesep,'R*']);
runNames = {realinedDir.name}';
nRuns = numel(runNames);
FilePath_EPIs = cell(nRuns,1);
for iRun = 1:nRuns
    sourceDir = [realinedDir(iRun).folder,filesep,runNames{iRun}];
    dirList = dir(sprintf('%s%su_*.nii',sourceDir,filesep));
    FilePath_EPIs{iRun} = cellfun(...
        @(x,y)[x,filesep,y],...
        {dirList.folder}',...
        {dirList.name}',...
        'UniformOutput',false);
end

%% Create and execute the SPM batch

% scans is a [1,nRuns] cell array of cell arrays containing file paths.
for iRun = 1:nRuns
    SpmBatch{1}.spm.temporal.st.scans{iRun} = FilePath_EPIs{iRun};
end
% nslices set to zero to appease the spm jobman (Unused).
SpmBatch{1}.spm.temporal.st.nslices = 0;

SpmBatch{1}.spm.temporal.st.tr = 2.2; % In seconds.

% TA set to zero to appease the spm jobman (Unused).
SpmBatch{1}.spm.temporal.st.ta = 0;

% Specify the in-volume timing of every slice in ms.
% SPM will know you are providing slice time for sliceOrder as long as:
% ~isequal(1:nSlices,sort(sliceOrder))
SpmBatch{1}.spm.temporal.st.so = getSliceT();

% Specify the timing of the "reference slice" in ms.
% SPM knows this parameter is not an index because of the above condition.
SpmBatch{1}.spm.temporal.st.refslice = 1056;

SpmBatch{1}.spm.temporal.st.prefix = 'a';

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);


%% Move new files
for iRun = 1:nRuns
    sourceDir = [realinedFolder,filesep,runNames{iRun}];
    targetDir = [epiFolder,filesep,'2_Temporal',filesep,runNames{iRun}];
    mkdir(targetDir);
    movefile(...
        [sourceDir,filesep,'au_*.nii'],...
        targetDir);
end

return

function [tSlice] = getSliceT()
% gap: The mean inter-slice gap in ms
gap = 66;
tSignal = (0:32)'*gap;
% Slice order is interleaved, starting on odd numbered slices (1-ordered)
sliceOrder = reshape([(1:17);[(18:33),NaN]],34,1);
sliceOrder = sliceOrder(~isnan(sliceOrder));
sliceOrder = [sliceOrder;sliceOrder];
tSlice = tSignal(sliceOrder);
return