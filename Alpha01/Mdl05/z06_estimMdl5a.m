function [] = z06_estimMdl5a(G)
% Estimate second-level GLM for Mdl05.
%
%   z06_FullInteraction2ndLvlMdl(G) estimates second-level SPM model for
%   Mdl05 across all subjects in the group defined by G.
%
%   Inputs:
%     G  - Group identifier consumed by getSubjectIds(G).
%
%   Requirements / Assumptions:
%     • Design matrix file exist at:
%         */_Group/[G]/Alpha01/Mdl5/X.mat
%
%   Model details:
%     • Formula: zTemplate ~ (1|SubjectId) + colocation + colocation:zPnonc
%
%   Outputs (per group):
%     • SPM.mat and parameter estimate images (beta_*.nii) saved to:
%         */_Group/[G]/Alpha01/Mdl5/
%   Notes:
%     • Residual images are not written (write_residuals = 0).
%     • T contrasts are written.

%% Cd out
wd = pwd;
cd ..;

%% Load the design mat
dirs.Data = '../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl5 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl05a');
Rfn = fullfile(dirs.Mdl5,'X.mat');
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
    'Alpha01','Mdl05',G,char(y)),subjectIds,srchLghtFn,...
    'UniformOutput',false);

%% Job definition: Specify
spmBatch{1}.spm.stats.factorial_design.dir = {dirs.Mdl5};
spmBatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;
spmBatch{1}.spm.stats.factorial_design.des.mreg.mcov = ...
    struct('c', {}, 'cname', {}, 'iCC', {});
spmBatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
spmBatch{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
spmBatch{1}.spm.stats.factorial_design.multi_cov.files = {Rfn};
spmBatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
spmBatch{1}.spm.stats.factorial_design.multi_cov.iCC = 5;
spmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
% No implicit masking
spmBatch{1}.spm.stats.factorial_design.masking.im = 1;
% Explicit masking with group EPI mask
spmBatch{1}.spm.stats.factorial_design.masking.em = {maskFn};
spmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Job definition: Estimate
spmBatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(dirs.Mdl5,'SPM.mat')};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrasts
spmMatfn = fullfile(dirs.Mdl5,'SPM.mat');
spmBatch{3}.spm.stats.con.spmmat(1) = {spmMatfn};

pNonc = get_pNonc('G1');
zPnonc = zscore(pNonc.pNonc);
zPnonc = zPnonc(order);
zPnonc_d = zPnonc; 

conNames = cell(5,1);
%intercept
conNames{1} = 'intercept';
H.intercept = [ones(1,nSubs).*(1/nSubs),0,0];
%main effects
conNames{2} = 'zPnonc';
H.zPnonc = [zPnonc_d',0,0];
conNames{3} = 'coloc';
H.colocation = [zeros(1,nSubs),1,0];
%interaction effect
conNames{4} = 'zPnonc:coloc';
H.zPnoncXcoloc =  [zeros(1,nSubs),0,1];
%simple effect
conNames{5} = 'zPnonc_coloc=+1';
H.zPnoncInColocP1 = [zPnonc_d',0,1];
fields = fieldnames(H);

for iH = 1:numel(conNames)
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.weights = H.(fields{iH});
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'none';
end
spmBatch{3}.spm.stats.con.delete = 1;

%% Job execution
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

% Cd back in
cd(wd);
return