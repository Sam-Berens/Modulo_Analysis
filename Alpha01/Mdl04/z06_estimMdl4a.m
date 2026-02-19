function [] = z06_estimMdl4a(G)
% Estimate second-level GLM for Mdl04.
%
%   z06_FullInteraction2ndLvlMdl(G) estimates second-level SPM model for
%   Mdl04 across all subjects in the group defined by G.
%
%   Inputs:
%     G  - Group identifier consumed by getSubjectIds(G).
%
%   Requirements / Assumptions:
%     • Design matrix file exist at:
%         */_Group/[G]/Alpha01/Mdl4/X.mat
%
%   Model details:
%     • Formula: zTemplate ~ (1|SubjectId) + colocation + colocation:zPnonc
%
%   Outputs (per group):
%     • SPM.mat and parameter estimate images (beta_*.nii) saved to:
%         */_Group/[G]/Alpha01/Mdl4/
%   Notes:
%     • Residual images are not written (write_residuals = 0).
%     • No contrasts are specified here.

%% Cd out
wd = pwd;
cd ..;

%% Load the design mat
dirs.Data = '../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl4 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl04_r3');
Rfn = fullfile(dirs.Mdl4,'X.mat');
loadStrct = load(Rfn);
names = loadStrct.names;

%% Get group mask
maskFn = fullfile(dirs.Group, 'Structural','GrpEpiMask00',...
    'G1_GrpEpiMask00.nii');
%get 'scan' names in correct order for R rows
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
names = categorical(names(1:nSubs))';
order = arrayfun(@(x) find(x == subjectIds),names);
subjectIds = subjectIds(order);
subjectIds = repelem(subjectIds,2,1);
%very important that the colocation=-1 goes first
srchLghtFn = categorical({'wzTemplate_colocation=-1.nii';...
    'wzTemplate_colocation=+1.nii'});
srchLghtFn = repmat(srchLghtFn,[nSubs,1]);
scans = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Searchlight',G,char(y)),subjectIds,srchLghtFn,...
    'UniformOutput',false);

%% Job definition: Specify
SpmBatch{1}.spm.stats.factorial_design.dir = {dirs.Mdl4};
SpmBatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;
SpmBatch{1}.spm.stats.factorial_design.des.mreg.mcov = ...
    struct('c', {}, 'cname', {}, 'iCC', {});
SpmBatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
SpmBatch{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
SpmBatch{1}.spm.stats.factorial_design.multi_cov.files = {Rfn};
SpmBatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
SpmBatch{1}.spm.stats.factorial_design.multi_cov.iCC = 5;
SpmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
% No implicit masking
SpmBatch{1}.spm.stats.factorial_design.masking.im = 1;
% Explicit masking with group EPI mask
SpmBatch{1}.spm.stats.factorial_design.masking.em = {maskFn};
SpmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
SpmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
SpmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Job definition: Estimate
SpmBatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(dirs.Mdl4,'SPM.mat')};
SpmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
SpmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Job execution
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

% Cd back in
cd(wd);
return