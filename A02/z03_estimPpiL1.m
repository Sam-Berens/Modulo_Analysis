% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_rVisual.nii';
% roiInfo.id = 'rVisual';

function [] = z03_estimPpiL1(G,roiInfo)
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

% Set some constants
tr = 2.2;
dirs_Data = ['..',filesep,'..',filesep,'Data'];
% epiMask = [dirs_Data,filesep,'_Group',filesep,'G1',...
%     filesep,'Structural',filesep,'GrpEpiMask00',...
%     filesep,'G1_GrpEpiMask00.nii'];

subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject = 1:nSubs
    subjectId = char(subjectIds(iSubject));
    %get subject-specific filenames and info
    dirs.Subject = fullfile(dirs_Data,subjectId);
    dirs.A02 = fullfile(dirs.Subject,'Analysis','A02');
    dirs.PPI = fullfile(dirs.A02,sprintf('PPI_%s',roiInfo.id));
    dirs.EPI = [dirs.Subject,filesep,'EPI'];
    dirs.G = [dirs.EPI,filesep,G];
    dirs.Y = [dirs.G,filesep,'k05'];
    epiMask = [dirs.EPI,filesep,'G1',filesep,sprintf('w_%s_epiMask00.nii',subjectId)];
    A02spmFn = fullfile(dirs.A02,'SPM.mat');
    tmp = load(A02spmFn);
    nRuns = size(tmp.SPM.Sess,2);
    % Set the EPI filenames
    epiFns = getEpiFns(dirs.Y);
    % Set the movement regressor mat filenames
    [rpFns] = getRpsFns(dirs.EPI);

    % need to load all runs for each PPI column and concatinate them (in
    % the same order as the scans you're loading (Y in this model is just BOLD of the whole brain so I think this means we litterally just want the smoothed EPI images?)
    a = cell(nRuns,1);
    b = cell(nRuns,1);
    for iRun=1:nRuns
        tmp = load(sprintf('%s%sPPI_R%i-rVisual*a.mat',dirs.PPI,filesep,iRun));
        a{iRun} = tmp.PPI;
        tmp = load(sprintf('%s%sPPI_R%i-rVisual*b.mat',dirs.PPI,filesep,iRun));
        b{iRun} = tmp.PPI;
    end
    clear tmp;
    
    %Estimate 1st level GLM on whole brain data, using PPI predictors
    estimPPImodel(dirs,tr,epiMask,roiInfo,epiFns,rpFns,a,b)

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

function [] = estimPPImodel(dirs,tr,epiMask,roiInfo,epiFns,rpFns,a,b)
%% Job definition: Specify
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

    %Psychological
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(1).name = 'Psych-a';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(1).val = a{iRun}.P;
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(2).name = 'Psych-b';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(2).val = b{iRun}.P;

    %Physiological
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(3).name = 'rVisual-BOLD';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(3).val = a{iRun}.Y; %this is indentical for both

    %Interaction
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(4).name = 'PPI-interaction-a';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(4).val = a{iRun}.ppi;
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(5).name = 'PPI-interaction-b';
    spmBatch{1}.spm.stats.fmri_spec.sess(iRun).regress(5).val = b{iRun}.ppi;

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

spmMatFn = [dirs.PPI,filesep,'SPM.mat'];
spmBatch{2}.spm.stats.fmri_est.spmmat = {spmMatFn};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Make contrasts
conNames = cell(7,1);
%main effect of a
conNames{1} = 'a';
H.a = 1;
%main effect of b
conNames{2} = 'b';
H.b = [0,1];
%main effect of roi
conNames{3} = roiInfo.id;
H.Y = [0,0,1];
%PPI interaction effect (a)
conNames{4} = sprintf('a:%s',roiInfo.id);
H.axY =  [0,0,0,1];
%PPI interaction effect (b)
conNames{5} = sprintf('b:%s',roiInfo.id);
H.bxY =  [0,0,0,0,1];
%'Main effect' of x:roi (where x is either a or b)
conNames{6} = sprintf('a:%s + b:%s',roiInfo.id);
H.axY_plus_bxY =  [0,0,0,1,1];
%'Interaction' of interaction effects (diff between a:roi vs b:roi)
conNames{7} = sprintf('a:%s - b:%s',roiInfo.id);
H.axY_minus_bxY =  [0,0,0,1,-1];
fields = fieldnames(H);

spmBatch{3}.spm.stats.con.spmmat = {spmMatFn};
for iH = 1:numel(conNames)
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.weights = H.(fields{iH});
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'replsc';
end
spmBatch{3}.spm.stats.con.delete = 1;

spm_jobman('initcfg');
spm_jobman('run',spmBatch);
return