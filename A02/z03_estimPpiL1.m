function [] = z03_estimL1PPI(G,roiId)
% Estimate first-level GLMs for A02, with PP1 predictors for a given roi.
%
%   z01_estimL1(G) estimates first-level SPM models for A02 across all
%   subjects in the group defined by G.
%
%   Inputs:
%     G  - Group identifier consumed by getSubjectIds(G).
%
%   Requirements / Assumptions:
%     • PPI structure files exist at:
%         */[SubjectId]/Analysis/A02/PPI_[roiId]/PPI_[iRun]*[stimX].mat
%     • Realignment parameter files exist at:
%         */[SubjectId]/EPI/RP*.mat
%     • EPI images exist in per-run folders:
%         */[SubjectId]/EPI/[Group]/k05/R*/<*.nii>
%
%   Model details:
%     • Units: scans (TR = 2.2 s)
%     • Basis functions: canonical HRF, no time/dispersion derivatives
%     • Microtime resolution: 66
%     • Global normalisation: none
%     • Masking: implicit, threshold = 0; explicit mask
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
% EG:
% roiId = 'rVisual';

%% 
dirs_Data = ['..',filesep,'..',filesep,'Data'];
dirs_Data = dir(dirs_Data);
dirs_Data = dirs_Data(1).folder;
subjectIds = getSubjectIds(G);
parfor iSubject = 1:numel(subjectIds)

    dirs = struct;
    subjectId = char(subjectIds(iSubject));
    dirs.Subject = fullfile(dirs_Data,subjectId);

    % Set dirs.PPI (the destination)
    dirs.A02 = fullfile(dirs.Subject,'Analysis','A02');
    dirs.PPI = fullfile(dirs.A02,sprintf('PPI_%s',roiId));
    toDel = [
        dir(fullfile(dirs.PPI,'*.nii'));
        dir(fullfile(dirs.PPI,'SPM.mat'))];
    for iToDel = 1:numel(toDel)
        delete(fullfile(toDel(iToDel).folder,toDel(iToDel).name));
    end

    % Set dirs.Y (the source)
    dirs.EPI = fullfile(dirs.Subject,'EPI');
    dirs.G = fullfile(dirs.EPI,G);
    dirs.Y = fullfile(dirs.G,'k05');

    epiMask = fullfile(...
        dirs.EPI,...
        'G1',...
        sprintf('w_%s_epiMask00.nii',subjectId));

    % Get the number of runs
    A02SPMFn = fullfile(dirs.A02,'SPM.mat');
    A02SPM = load(A02SPMFn);
    nRuns = size(A02SPM.SPM.Sess,2);

    % Set the EPI filenames
    epiFns = getEpiFns(dirs.Y);

    % Set the movement regressor mat filenames
    rpFns = getRpsFns(dirs.EPI);

    % Load the PPI structures into two cell arrays, for a and b seperately
    a = cell(nRuns,1);
    b = cell(nRuns,1);
    for iRun = 1:nRuns
        regressor_A = load(...
            sprintf('%s%sPPI_R%i-%s*a.mat',dirs.PPI,filesep,iRun,roiId));
        a{iRun} = regressor_A.PPI;
        regressor_B = load(...
            sprintf('%s%sPPI_R%i-%s*b.mat',dirs.PPI,filesep,iRun,roiId));
        b{iRun} = regressor_B.PPI;
    end
    
    % Estimate 1st level GLM on whole brain data, using PPI predictors
    estimPPImodel(dirs,epiMask,roiId,epiFns,rpFns,a,b);

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

function [] = estimPPImodel(dirs,epiMask,roiId,epiFns,rpFns,a,b)

tr = 2.2;

% Job definition: Specify
spmBatch{1}.spm.stats.fmri_spec.dir = {dirs.PPI};
spmBatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
spmBatch{1}.spm.stats.fmri_spec.timing.RT = tr;
spmBatch{1}.spm.stats.fmri_spec.timing.fmri_t = 66;
spmBatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 33;
    
for iRun = 1:numel(rpFns)
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).scans = epiFns{iRun};
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).cond = struct('name',...
        {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod',...
        {}, 'orth', {});
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).multi = {''};

    % Psychological
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(1).name = ...
        'a';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(1).val = ...
        a{iRun}.P;
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(2).name = ...
        'b';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(2).val = ...
        b{iRun}.P;

    % Physiological
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(3).name = ...
        roiId;
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(3).val = ...
        a{iRun}.Y;
    % We only take the Physiological regressor for a but it is indentical
    % for both a and b.

    % Interaction
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(4).name = ...
        sprintf('a:%s',roiId);
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(4).val = ...
        a{iRun}.ppi;
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(5).name = ...
        sprintf('b:%s',roiId);
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(5).val = ...
        b{iRun}.ppi;

    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = rpFns(iRun);
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
end

% No factorial structure
spmBatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% No time or dispersion derivatives
spmBatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
% No Voltera interaction or 2nd order effect for the HRFs
spmBatch{1}.spm.stats.fmri_spec.volt = 1;
% No per-voxel scaling
spmBatch{1}.spm.stats.fmri_spec.global = 'None';
% No implicit masking
spmBatch{1}.spm.stats.fmri_spec.mthresh = 0;
% Explicit masking with custom group mask
spmBatch{1}.spm.stats.fmri_spec.mask = {epiMask};
% AR(1)
spmBatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Estimate
spmMatFn = [dirs.PPI,filesep,'SPM.mat'];
spmBatch{2}.spm.stats.fmri_est.spmmat = {spmMatFn};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Make contrasts
spmBatch{3}.spm.stats.con.spmmat = {spmMatFn};
spmBatch{3}.spm.stats.con.delete = 1;

spmBatch{3}.spm.stats.con.consess{1}.tcon.name = 'a';
spmBatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
spmBatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';

spmBatch{3}.spm.stats.con.consess{2}.tcon.name = 'b';
spmBatch{3}.spm.stats.con.consess{2}.tcon.weights = [0,1];
spmBatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';

spmBatch{3}.spm.stats.con.consess{3}.tcon.name = roiId;
spmBatch{3}.spm.stats.con.consess{3}.tcon.weights = [0,0,1];
spmBatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';

spmBatch{3}.spm.stats.con.consess{4}.tcon.name = sprintf('a:%s',roiId);
spmBatch{3}.spm.stats.con.consess{4}.tcon.weights = [0,0,0,1,0];
spmBatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';

spmBatch{3}.spm.stats.con.consess{5}.tcon.name = sprintf('b:%s',roiId);
spmBatch{3}.spm.stats.con.consess{5}.tcon.weights = [0,0,0,0,1];
spmBatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';

% The following contrast compute a mean
spmBatch{3}.spm.stats.con.consess{6}.tcon.name = ...
    sprintf('a:%s + b:%s',roiId,roiId);
spmBatch{3}.spm.stats.con.consess{6}.tcon.weights = [0,0,0,0.5,0.5];
spmBatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';

spmBatch{3}.spm.stats.con.consess{7}.tcon.name = ...
    sprintf('a:%s - b:%s',roiId,roiId);
spmBatch{3}.spm.stats.con.consess{7}.tcon.weights = [0,0,0,1,-1];
spmBatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';

%% Run Job
spm_jobman('initcfg');
spm_jobman('run',spmBatch);
return