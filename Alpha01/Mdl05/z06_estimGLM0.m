function [] = z06_estimGLM0(G)
% Estimate Mdl05/GLM0.
%
%   z06_estimGLM0(G) estimates a second-level SPM model (Mdl05/GLM0) for
%   all subjects in the group defined by G.
%
%   Inputs:
%     G  - Group identifier consumed by getSubjectIds(G).
%
%   Requirements / Assumptions:
%     • Design matrix file exist at:
%         */_Group/[G]/Alpha01/Mdl05/GlM0/X.mat
%
%   Model details:
%     • Formula: zTemplate ~ (1|SubjectId) + colocation + colocation:zPnonc
%
%   Outputs (per group):
%     • SPM.mat and parameter estimate images (beta_*.nii) saved to:
%         */_Group/[G]/Alpha01/Mdl05/GlM0
%   Notes:
%     • Residual images are not written (write_residuals = 0).
%     • No contrasts are specified here.


%% Load the design mat
dirs.Data = fullfile('..','..','..','Data');
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl05 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl05');
dirs.GLM0 = fullfile(dirs.Mdl05,'GLM0');
X = fullfile(dirs.GLM0,'X.mat');
temp = load(X);
names = temp.names;

%% Get group mask
maskFn = fullfile(dirs.Group, 'Structural','GrpEpiMask00',...
    'G1_GrpEpiMask00.nii');

%% Get Y filenames in correct order for rows in R
pNonc = get_pNonc(G);
subjectIds = pNonc.subjectId;
nSubjects = numel(subjectIds);
if ~isequal(categorical(names(1:nSubjects)'),pNonc.subjectId)
    error('Something is wrong.');
end

subjectIds = cellstr(repelem(subjectIds,2,1));
yFilenames = repmat({....
    'wzTemplate_colocation=-1.nii';... % Important this goes first
    'wzTemplate_colocation=+1.nii'},...
    nSubjects,1);
yFilenames = cellfun(@(x,y) fullfile(...
    dirs.Data,...
    x,...
    'Analysis',...
    'Alpha01',...
    'Mdl05',...
    G,...
    y),...
    subjectIds,yFilenames,...
    'UniformOutput',false);

%% Job definition: Specify
spmBatch{1}.spm.stats.factorial_design.dir = {dirs.GLM0};
spmBatch{1}.spm.stats.factorial_design.des.mreg.scans = yFilenames;
spmBatch{1}.spm.stats.factorial_design.des.mreg.mcov = ...
    struct('c', {}, 'cname', {}, 'iCC', {});
spmBatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
spmBatch{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
spmBatch{1}.spm.stats.factorial_design.multi_cov.files = {X};
spmBatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
spmBatch{1}.spm.stats.factorial_design.multi_cov.iCC = 5;
spmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;

% No implicit masking
spmBatch{1}.spm.stats.factorial_design.masking.im = 0;

% Explicit masking with group EPI mask
spmBatch{1}.spm.stats.factorial_design.masking.em = {maskFn};
spmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Job definition: Estimate
spmBatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(dirs.GLM0,'SPM.mat')};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Job definition: Contrast spec
zPnonc = zscore(pNonc.pNonc);
conNames = cell(4,1);

conNames{1} = 'intercept';
H.intercept = [ones(1,nSubjects).*(1/nSubjects),0,0];

conNames{2} = 'zPnonc';
H.zPnonc = [zPnonc',0,0];

conNames{3} = 'coloc';
H.coloc = [zeros(1,nSubjects),1,0];

conNames{4} = 'zPnonc:coloc';
H.zPnoncXcoloc =  [zeros(1,nSubjects),0,1];

fields = fieldnames(H);
spmBatch{3}.spm.stats.con.spmmat(1) = {fullfile(dirs.GLM0,'SPM.mat')};
for iH = 1:numel(conNames)
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.weights = H.(fields{iH});
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'none';
end
spmBatch{3}.spm.stats.con.delete = 1;

%% Job execution
spm_jobman('initcfg');
spm_jobman('run',spmBatch);
return