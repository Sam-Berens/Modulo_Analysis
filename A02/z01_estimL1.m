function [] = z01_estimL1(G)
% Estimate first-level GLMs for A02.
%
%   z01_estimL1(G) estimates first-level SPM models for A02 across all
%   subjects in the group defined by G.
%
%   Inputs:
%     G  - Group identifier consumed by getSubjectIds(G).
%
%   Requirements / Assumptions:
%     • Event specification files exist at:
%         */[SubjectId]/Analysis/A02/EventSpec_R*.mat
%     • Realignment parameter files exist at:
%         */[SubjectId]/EPI/RP*.mat
%     • EPI images exist in per-run folders:
%        */[SubjectId]/EPI/[Group]/k05/R*/<*.nii>
%     • Run indices encoded in EventSpec_R*.mat match those in the EPI
%       run folders; the function checks and errors if mismatched.
%
%   Model details:
%     • Units: scans (TR = 2.2 s)
%     • Basis functions: canonical HRF, no time/dispersion derivatives
%     • Microtime resolution: 66
%     • Global normalisation: none
%     • Masking: implicit, threshold = 0; no explicit mask
%     • High-pass filter: 128 s
%     • Serial correlations: AR(1)
%     • Estimation method: Classical (ReML)
%
%   Outputs (per subject):
%     • SPM.mat and parameter estimate images (beta_*.nii) saved to:
%         */[Subject]/Analysis/A02/
%
%   Notes:
%     • Residual images are not written (write_residuals = 0).
%     • No contrasts are specified here.

% Set some constants
tr = 2.2;
dirs_Data = ['..',filesep,'..',filesep,'Data'];
% epiMask = [dirs_Data,filesep,'_Group',filesep,'G1',...
%     filesep,'Structural',filesep,'GrpEpiMask00',...
%     filesep,'G1_GrpEpiMask00.nii'];
subjectIds = getSubjectIds(G);
G = 'G1';
for iSubject = 1:numel(subjectIds)
    dirs = struct;

    % Set required directory paths
    subjectId = char(subjectIds(iSubject));
    dirs.Subject = [dirs_Data,filesep,subjectId];
    dirs.EPI = [dirs.Subject,filesep,'EPI'];
    dirs.G = [dirs.EPI,filesep,G];
    dirs.Y = [dirs.G,filesep,'k05'];
    dirs.A02 = [dirs.Subject,filesep,'Analysis',filesep,'A02'];
    %Set the normed EPI mask filename
    epiMask = [dirs.EPI,filesep,'G1',filesep,sprintf('w_%s_epiMask00.nii',subjectId)];

    % Set the filenames of the realignment parameters
    rpsFns = getRpsFns(dirs.EPI);

    % Set the EPI filenames
    epiFns = getEpiFns(dirs.Y);

    % Loop through stimIds to estimate
    if ~exist([dirs.A02,filesep,'ResMS.nii'],'file')
        estimL1(tr,dirs.A02,epiFns,rpsFns,epiMask);
    end

end
return

function [fileList] = getRpsFns(path2data)
fullpath = @(s)[s.folder,filesep,s.name];
fileList  = dir([path2data,filesep,'RP*.mat']);
[~,ord] = sort({fileList.name});
fileList = fileList(ord);
fileList = arrayfun(fullpath,fileList,'UniformOutput',false);
return

function [epiFns] = getEpiFns(path2data)
fullpath = @(s)[s.folder,filesep,s.name];
runList  = dir([path2data,filesep,'R*']);
[~,ord]  = sort({runList.name});
runList  = runList(ord);
runList = arrayfun(fullpath,runList,'UniformOutput',false);
epiFns = cell(size(runList));
for iRun = 1:numel(runList)
    temp = dir([runList{iRun},filesep,'*.nii']);
    epiFns{iRun} = arrayfun(fullpath,temp,'UniformOutput',false);
end
return

function [] = estimL1(tr,outDir,epiFns,rpsFns,epiMaskFn)

% Set some constants
fullpath = @(s)[s.folder,filesep,s.name];
empty0 = struct(...
    'name',{},'onset',{},'duration',{},'tmod',{},'pmod',{},'orth',{});
empty1 = struct('name',{},'val',{});

% Set the event spec file names
eventSpecs = dir([outDir,filesep,'EventSpec_R*.mat']);
[~,ord]  = sort({eventSpecs.name});
eventSpecs  = eventSpecs(ord);
eventSpecs = arrayfun(fullpath,eventSpecs,'UniformOutput',false);

%% Check there is no mismatch in the RunIds ...
% ... between epiFns and eventSpecs (rpsFns may be misaligned by design)
runIds.eventSpecs = cellfun(...
    @(s,ii)str2double(s(ii)),...
    eventSpecs,num2cell(cellfun(@(ii)ii+2,strfind(eventSpecs,'_R'))));
firstNames = cellfun(@(c)c{1},epiFns,'UniformOutput',false);
runIds.epiFns = cellfun(...
    @(s,ii)str2double(s(ii)),...
    firstNames,num2cell(cellfun(@(ii)ii+1,strfind(firstNames,'R'))));
if  var([numel(epiFns);numel(rpsFns);numel(eventSpecs)]) ~= 0
    error('Mismatching runs between epiFns, rpsFns, and/or eventSpecs.');
end
if ~all(runIds.eventSpecs==runIds.epiFns)
    error('Mismatch between Run IDs!');
end

%% Job definition: Specify
SpmJob{1}.spm.stats.fmri_spec.dir = {outDir};
SpmJob{1}.spm.stats.fmri_spec.timing.units = 'scans';
SpmJob{1}.spm.stats.fmri_spec.timing.RT = tr;
SpmJob{1}.spm.stats.fmri_spec.timing.fmri_t = 66;
SpmJob{1}.spm.stats.fmri_spec.timing.fmri_t0 = 33;

% Run loop
for iRun = 1:numel(eventSpecs)
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).scans = epiFns{iRun};
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).cond = empty0;
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).multi = eventSpecs(iRun);
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).regress = empty1;
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = rpsFns(iRun);
    SpmJob{1}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
end

% No factorial structure
SpmJob{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});

% No time or dispersion derivatives
SpmJob{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

% No Voltera interaction or 2nd order effect for the HRFs
SpmJob{1}.spm.stats.fmri_spec.volt = 1;

% No per-voxel scaling
SpmJob{1}.spm.stats.fmri_spec.global = 'None';

% No implicit masking
SpmJob{1}.spm.stats.fmri_spec.mthresh = 0;

% Explicit masking with custom group epi mask
SpmJob{1}.spm.stats.fmri_spec.mask = {epiMaskFn};

% AR(1)
SpmJob{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Job definition: Estimate
SpmJob{2}.spm.stats.fmri_est.spmmat = {[outDir,filesep,'SPM.mat']};
SpmJob{2}.spm.stats.fmri_est.write_residuals = 0;
SpmJob{2}.spm.stats.fmri_est.method.Classical = 1;

%% Job execution
spm_jobman('initcfg');
spm_jobman('run',SpmJob);
return