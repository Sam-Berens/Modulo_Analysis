function []= estimL1_Alpha00(groupFunc)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;

subjectIds = groupFunc(); %this will be e.g. getG1()
nSubs = numel(subjectIds);


for ii=1:nSubs
    subDir = [dataDir,filesep,subjectIds{ii}];
    anlDir = [subDir,filesep,'Analysis'];
    alphDir = [anlDir,'Alpha00'];
    epiDir = [subDir,filesep,'EPI'];
    epi2Dir = [epiDir,filesep,'2_Temporal'];
    runList  = dir([epi2Dir,filesep,'R*']);
    runList = unique([runList.name]');
    runList = runList(1:end-1);
    runIds = (arrayfun(@(x) str2double(x),runList))';
    nRuns = numel(runList);

    %% EPIs
    FilePath_EPIs = cell(nRuns,1);
    for iRun = runIds
        epiFileList = dir([epi2Dir,filesep,'R',int2str(runIds(iRun)),'au_*.nii']);
        FilePath_EPIs{iRun} = cellfun(...
            @(x,y)[x,filesep,y],...
            {epiFileList.folder}',...
            {epiFileList.name}',...
            'UniformOutput',false);
    end

    %% RPS
    rpFileList = dir([epiDir,filesep,'RP*.mat']);
    FileName_RPs = cellfun(...
        @(x,y)[x,filesep,y],...
        {rpFileList.folder}',...
        {rpFileList.name}',...
        'UniformOutput',false);

    %set TR
    for iRun = runIds
        %% Stim loop
        for stimId = 0:5
            specAndEstim(alphDir,FilePath_EPIs{iRun},runIds(iRun),stimId);
        end
    end

end

return

function [] = specAndEstim(alphDir,FilePath_EPIs,runId,stimId)
TR = 2.2;

outDir = [alphDir,filesep,sprintf('i%i_R%i',stimId,runId)];
spmMatFile = [outDir,filesep,'SPM.mat'];
eventSpec = dir([outDir,filesep,'EventSpec.mat']);
eventSpec = [eventSpec.folder,filesep,eventSpec.name];

SpmBatch = [{},{}];
SpmBatch{1}.spm.stats.fmri_spec.dir = {outDir};
SpmBatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
SpmBatch{1}.spm.stats.fmri_spec.timing.RT = TR;
SpmBatch{1}.spm.stats.fmri_spec.timing.fmri_t = 66; % aquisition times were start of each aquistion not halfway through so we need double the resolution to specfiy reference slice's point
SpmBatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 33 ; %....and 12
SpmBatch{1}.spm.stats.fmri_spec.sess(1).scans = FilePath_EPIs;
SpmBatch{1}.spm.stats.fmri_spec.sess(1).cond = struct(...
    'name', {}, ... %THIS IS JUST TO MAKE THE GUI HAPPY?
    'onset', {}, ...
    'duration', {}, ...
    'tmod', {}, ...
    'pmod', {}, ...
    'orth', {});
SpmBatch{1}.spm.stats.fmri_spec.sess(1).multi = eventSpec; %THEN PUT IN THE ACTUAL VALUES HERE
SpmBatch{1}.spm.stats.fmri_spec.sess(1).regress = struct( ...
    'name', {}, 'val', {}); %THIS IS JUST TO MAKE THE GUI HAPPY
SpmBatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = FileName_RPs;  %THEN PUT IN THE ACTUAL VALUES HERE
SpmBatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; %defaults

SpmBatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
SpmBatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % Do we not use time/dispersion derivatives in the first level model because we're not estimating the patterns at a group level?
SpmBatch{1}.spm.stats.fmri_spec.volt = 1; %(Do we not use voltera/2nd order interaction of HRFs because our smallest stim seperation time is 1s which is > 800ms at which nonlinearity has been observed?
SpmBatch{1}.spm.stats.fmri_spec.global = 'None';
SpmBatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
SpmBatch{1}.spm.stats.fmri_spec.mask = {''}; %MAKE DECISION ABOUT YOUR MASK VS SPM MASK
SpmBatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; %default

%% TO DO find out if this way of specifying batches works?
%Estimate the model once the spmMatFile is made
%the batch numbers are such that for someone with 5 runs, 1-30 is the model
%specification job for producing the spm.mat files and then
%31-60 is the estimation jobs for each of the models.
SpmBatch{2}.spm.stats.fmri_est.spmmat = {spmMatFile};
SpmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
SpmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

%%
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);
return