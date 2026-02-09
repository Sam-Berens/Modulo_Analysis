function [] = z07_estimMdl4b()

% Cd out
wd = pwd;
cd ..;

%% Get a table that includes all subjectIds and performance stats
InputTable = get_pNonc(G);
InputTable.cpNonc = InputTable.pNonc - mean(InputTable.pNonc);
nSubjects = size(InputTable,1);
InputTable.fn_colocPos = cell(nSubjects,1);

%% Set the path to the /Data directory
dirs.Data = ['..',filesep,'..',filesep,'Data'];

%% Set the output path
dirs.Output = [dirs.Data,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Analysis',...
    filesep,'Alpha01',...
    filesep,'colocation=+1'];
if ~exist(dirs.Output,'dir')
    mkdir(dirs.Output);
end

%% Get the input filenames
for iSubject = 1:nSubjects
    InputTable.fn_colocPos{iSubject} = [dirs.Data,...
        filesep,char(InputTable.subjectId(iSubject)),...
        filesep,'Analysis',...
        filesep,'Alpha01',...
        filesep,'Searchlight',...
        filesep,G,...
        filesep,'wzTemplate_colocation=+1.nii'];
end

%% Create and run the batch job
SpmJob = [{},{}];

% Output dir
SpmJob{1}.spm.stats.factorial_design.dir = {dirs.Output};

% Input fns
SpmJob{1}.spm.stats.factorial_design.des.mreg.scans = ...
    InputTable.fn_colocPos;

% Covariate vector
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.c = ...
    InputTable.cpNonc;

% Covariate name
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'cpNonc';

% Centering
SpmJob{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;

% Other stuff
SpmJob{1}.spm.stats.factorial_design.des.mreg.incint = 1;
SpmJob{1}.spm.stats.factorial_design.cov = ...
    struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
SpmJob{1}.spm.stats.factorial_design.multi_cov = ...
    struct('files', {}, 'iCFI', {}, 'iCC', {});
SpmJob{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
SpmJob{1}.spm.stats.factorial_design.masking.im = 1;
SpmJob{1}.spm.stats.factorial_design.masking.em = {''};
SpmJob{1}.spm.stats.factorial_design.globalc.g_omit = 1;
SpmJob{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
SpmJob{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Estimation
SpmJob{2}.spm.stats.fmri_est.spmmat = {[dirs.Output,filesep,'SPM.mat']};
SpmJob{2}.spm.stats.fmri_est.write_residuals = 0;
SpmJob{2}.spm.stats.fmri_est.method.Classical = 1;

% Run th job
spm_jobman('initcfg');
spm_jobman('run',SpmJob);

% Cd back in
cd(wd);

return


