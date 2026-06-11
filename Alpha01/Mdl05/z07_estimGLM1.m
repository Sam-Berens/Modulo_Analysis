function [] = z07_estimGLM1(G)

%% Get a table that includes all subjectIds and performance stats
pNonc = get_pNonc(G);
pNonc.zpNonc = zscore(pNonc.pNonc);
nSubjects = size(pNonc,1);
pNonc.fnY = cell(nSubjects,1);

%% Set the path to the Data directory
dirs.Data = fullfile('..','..','..','Data');

%% Set the output path
dirs.Output = [dirs.Data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Analysis',...
    filesep,'Alpha01',...
    filesep,'Mdl05',...
    filesep,'GLM1'];
if ~exist(dirs.Output,'dir')
    mkdir(dirs.Output);
end

%% Get group mask
maskFn = fullfile(dirs.Data,'_Group',G,'Structural','GrpEpiMask00',...
    'G1_GrpEpiMask00.nii');

%% Get the input filenames
for iSubject = 1:nSubjects
    pNonc.fnY{iSubject} = [dirs.Data,...
        filesep,char(pNonc.subjectId(iSubject)),...
        filesep,'Analysis',...
        filesep,'Alpha01',...
        filesep,'Mdl05',...
        filesep,G,...
        filesep,'wzTemplate_visSim.nii'];
end

%% Create and run the batch job
SpmJob = [{},{}];

% Output dir
SpmJob{1}.spm.stats.factorial_design.dir = {dirs.Output};

% Target fns
SpmJob{1}.spm.stats.factorial_design.des.mreg.scans = ...
    pNonc.fnY;

% Covariate vector
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.c = ...
    pNonc.zpNonc;

% Covariate name
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'zpNonc';

% Centering
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;

% Other stuff
SpmJob{1}.spm.stats.factorial_design.des.mreg.incint = 1;
SpmJob{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
SpmJob{1}.spm.stats.factorial_design.multi_cov = ...
    struct('files', {}, 'iCFI', {}, 'iCC', {});
SpmJob{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
SpmJob{1}.spm.stats.factorial_design.masking.im = 0; % No implicit mask
SpmJob{1}.spm.stats.factorial_design.masking.em = {maskFn};
SpmJob{1}.spm.stats.factorial_design.globalc.g_omit = 1;
SpmJob{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
SpmJob{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Estimation
spmMatfn = [dirs.Output,filesep,'SPM.mat'];
SpmJob{2}.spm.stats.fmri_est.spmmat = {spmMatfn};
SpmJob{2}.spm.stats.fmri_est.write_residuals = 0;
SpmJob{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrasts
SpmJob{3}.spm.stats.con.spmmat(1) = {spmMatfn};

conNames = cell(1,1);
%intercept
conNames{1} = 'intercept';
H.intercept = 1;
%main effect
conNames{2} = 'zPnonc';
H.zPnonc = [0,1];

for iH = 1:numel(conNames)
    SpmJob{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    SpmJob{3}.spm.stats.con.consess{iH}.tcon.weights = H.(conNames{iH});
    SpmJob{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'none';
end
SpmJob{3}.spm.stats.con.delete = 1;

% Run th job
spm_jobman('initcfg');
spm_jobman('run',SpmJob);

return