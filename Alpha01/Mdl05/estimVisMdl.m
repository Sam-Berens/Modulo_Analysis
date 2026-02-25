function [] = estimVisMdl(G)

%% Cd out
wd = pwd;
cd ..;

%get initial folder locations
dirs.Data = '../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group',G);
dirs.Mdl5c = fullfile(dirs.Group,'Analysis','Alpha01','Mdl05c');
if ~exist(dirs.Mdl5c,'dir')
    mkdir(dirs.Mdl5c);
end 

% Get group mask
maskFn = fullfile(dirs.Group, 'Structural','GrpEpiMask00',...
    'G1_GrpEpiMask00.nii');
%get 'scan' names and prep covariate vector
pNonc = get_pNonc('G1');
zPnonc = zscore(pNonc.pNonc);
subjectIds = pNonc.subjectId;
nSubs = numel(subjectIds);
srchLghtFn = categorical({'wzTemplate_visSim.nii'});
srchLghtFn = repmat(srchLghtFn,[nSubs,1]);
scans = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl05',G,char(y)),subjectIds,srchLghtFn,...
    'UniformOutput',false);

%% Job definition: Specify
spmBatch{1}.spm.stats.factorial_design.dir = {dirs.Mdl5c};
spmBatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
spmBatch{1}.spm.stats.factorial_design.cov.c = zPnonc;
spmBatch{1}.spm.stats.factorial_design.cov.cname = 'zPnonc';
spmBatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
spmBatch{1}.spm.stats.factorial_design.cov.iCC = 5;
spmBatch{1}.spm.stats.factorial_design.multi_cov = ...
    struct('files', {}, 'iCFI', {}, 'iCC', {});
% No implicit masking
spmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
spmBatch{1}.spm.stats.factorial_design.masking.im = 1;
% Explicit masking with group EPI mask
spmBatch{1}.spm.stats.factorial_design.masking.em = {maskFn};
spmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Job definition: Estimate
spmMatfn = fullfile(dirs.Mdl5c,'SPM.mat');
spmBatch{2}.spm.stats.fmri_est.spmmat(1) = {spmMatfn};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrasts
spmBatch{3}.spm.stats.con.spmmat(1) = {spmMatfn};

conNames = cell(1,1);
%intercept
conNames{1} = 'intercept';
H.intercept = 1;
%main effect
conNames{2} = 'zPnonc';
H.zPnonc = [0,1];

for iH = 1:numel(conNames)
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.weights = H.(conNames{iH});
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'none';
end
spmBatch{3}.spm.stats.con.delete = 1;

%% Job execution
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

% Cd back in
cd(wd);
return