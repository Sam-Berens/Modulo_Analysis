function []= specifyModel_l1(groupFunc)

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
    runList = unique([runs.name]');
    runList = runs(1:end-1);
    runs = (arrayfun(@(x) str2double(x),runList))';
    nRuns = numel(runList);

%% EPIs
FilePath_EPIs = cell(nRuns,1);
for iRun = runs
    epiFileList = dir([epi2Dir,filesep,'R',int2str(runs(iRun)),'au_*.nii']);
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
TR = 2.2;
nSessions = nRuns*6;
    for iRun = runs
            %% Stim loop
        for iStim = 1:6
            outDir = [alphaDir,filesep,sprintf('i%i_R%i', (iStim-1), runs(iRun))];
            spmMatFile = [outDir,filesep,'SPM.mat'];
            eventSpec = dir([outDir, filesep,'EventSpec.mat']);
            eventSpec = [eventSpec.folder, filesep, eventSpec.name];

            iSession = (iRun-1)*6 + iStim;
            SpmBatch = [{} {}];
            SpmBatch{iSession}.spm.stats.fmri_spec.dir = {outDir};
            SpmBatch{iSession}.spm.stats.fmri_spec.timing.units = 'scans';
            SpmBatch{iSession}.spm.stats.fmri_spec.timing.RT = TR;
            SpmBatch{iSession}.spm.stats.fmri_spec.timing.fmri_t = 16; %defaults 
            SpmBatch{iSession}.spm.stats.fmri_spec.timing.fmri_t0 = 8; %defaults 
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).scans = FilePath_EPIs{iRun};
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).cond = struct(...
                'name', {}, ... %THIS IS JUST TO MAKE THE GUI HAPPY?
                'onset', {}, ...
                'duration', {}, ...
                'tmod', {}, ...
                'pmod', {}, ...
                'orth', {});
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).multi = eventSpec; %THEN PUT IN THE ACTUAL VALUES HERE
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).regress = struct( ...
                'name', {}, 'val', {}); %THIS IS JUST TO MAKE THE GUI HAPPY
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).multi_reg = FileName_RPs(iRun);  %THEN PUT IN THE ACTUAL VALUES HERE
            SpmBatch{iSession}.spm.stats.fmri_spec.sess(1).hpf = 128; %defaults 

            SpmBatch{iSession}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            SpmBatch{iSession}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            SpmBatch{iSession}.spm.stats.fmri_spec.volt = 1;
            SpmBatch{iSession}.spm.stats.fmri_spec.global = 'None';
            SpmBatch{iSession}.spm.stats.fmri_spec.mthresh = 0.8;
            SpmBatch{iSession}.spm.stats.fmri_spec.mask = {''};
            SpmBatch{iSession}.spm.stats.fmri_spec.cvi = 'AR(1)';

            %% TO DO find out if this way of specifying batches works?
            %Estimate the model once the spmMatFile is made
            SpmBatch{nSessions+iSession}.spm.stats.fmri_est.spmmat = {spmMatFile};
            SpmBatch{nSessions+iSession}.spm.stats.fmri_est.write_residuals = 0;
            SpmBatch{nSessions+iSession}.spm.stats.fmri_est.method.Classical = 1;



        end
    end 
   
end 

% TO DO:
%% CHECK WHETHER YOUR SYSTEM OF CHECKING WHETHER THE JOBMAN HAS FINISHED WITH THE FIRST LOT OF BATCHES BEFORE THE SECOND LOT GETS DONE ACTUALLY WORKS 
spm_jobman('initcfg');
result = spm_jobman('run',SpmBatch(1:nSessions));
specDone = 0;
while ~specDone
 specDone = exist('result','var');
 if specDone
     spm_jobman('run',SpmBatch(nSessions:(nSessions*2)));
 end 
end 
    

