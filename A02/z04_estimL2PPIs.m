function [] = z04_estimL2PPIs(G,roiId)
%% Estimates a seperate 2nd level model each for 4 contrast, for a given roi
% cons are the two PPI effects, the main effect of PPI (irrespective of
% a or b condition) and the difference between the two PPI effects

% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_rVisual.nii';
% roiInfo.id = 'rVisual';
% Pick which contrast numbers you want to estimate models for and build a
% vector of them
conNames = {
    sprintf('a:%s',roiId);
    sprintf('b:%s',roiId);
    sprintf('a:%s + b:%s',roiId,roiId);
    sprintf('a:%s - b:%s',roiId,roiId);
    };
for iCon = 1:numel(conNames)
    estim2ndLvl(G,roiId,conNames{iCon});
end

return

function [] = estim2ndLvl(G,roiId,conName)

subjectIds = sort(getSubjectIds(G));
dirs.Data = fullfile('..','..','Data');
dirs.G = fullfile(dirs.Data,'_Group',G);

%% Get the contrast number
conSource = fullfile(...
    dirs.Data,...
    char(subjectIds(1)),...
    'Analysis',...
    'A02',...
    sprintf('PPI_%s',roiId));

dirList = dir(fullfile(conSource,'con*.nii'));
conList = fullfile({dirList.folder}',{dirList.name}');
conList = char(conList);
V = spm_vol(conList);
iCon = find(contains({V.descrip}',[conName,' - All Sessions']));
%TO DO FIX
if ~isscalar(iCon)
    [~,headEndIdx] = regexp(V(1).descrip, 'Contrast 1: ');
    iCon = find(arrayfun(@(x) matches({x.descrip(headEndIdx+1:end)}',[conName,' - All Sessions']),V));
    if ~isscalar(iCon)
        error('Ambiguous contrast name');
    end
end

%%
epiMask = fullfile(...
    dirs.G,...
    'Structural',...
    'GrpEpiMask00',...
    sprintf('%s_GrpEpiMask00.nii',G));

dirs.output = fullfile(...
    dirs.Data,...
    '_Group',...
    G,...
    'Analysis',...
    'A02',...
    sprintf('PPI_%s',roiId),...
    conName);

if ~exist(dirs.output,"dir")
    mkdir(dirs.output);
end

%% Get a full list of contrast file names
conList = cellstr(conList);
conFns = strrep(conList{iCon},char(subjectIds(1)),cellstr(subjectIds));

%% Get pNonc
pNonc = get_pNonc(G);
pNonc.zpNonc = zscore(pNonc.pNonc);
pNonc = sortrows(pNonc,'subjectId');
zpNonc = pNonc.zpNonc;

%% Job specification
spmBatch{1}.spm.stats.factorial_design.dir = {dirs.output};
spmBatch{1}.spm.stats.factorial_design.des.t1.scans = conFns;
spmBatch{1}.spm.stats.factorial_design.cov.c = zpNonc;
spmBatch{1}.spm.stats.factorial_design.cov.cname = 'zpNonc';
spmBatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
spmBatch{1}.spm.stats.factorial_design.cov.iCC = 5;
spmBatch{1}.spm.stats.factorial_design.multi_cov = ...
    struct('files', {}, 'iCFI', {}, 'iCC', {});

spmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
spmBatch{1}.spm.stats.factorial_design.masking.im = 1;
spmBatch{1}.spm.stats.factorial_design.masking.em = {epiMask};
spmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Estimation
spmBatch{2}.spm.stats.fmri_est.spmmat = {fullfile(dirs.output,'SPM.mat')};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Specify the contrasts for term in the model
spmBatch{3}.spm.stats.con.spmmat = {fullfile(dirs.output,'SPM.mat')};
spmBatch{3}.spm.stats.con.delete = 1;

spmBatch{3}.spm.stats.con.consess{1}.tcon.name = '+Intercept';
spmBatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
spmBatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

spmBatch{3}.spm.stats.con.consess{2}.tcon.name = '+zPnonc';
spmBatch{3}.spm.stats.con.consess{2}.tcon.weights = [0,1];
spmBatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

spmBatch{3}.spm.stats.con.consess{3}.tcon.name = '-Intercept';
spmBatch{3}.spm.stats.con.consess{3}.tcon.weights = -1;
spmBatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

spmBatch{3}.spm.stats.con.consess{4}.tcon.name = '-zPnonc';
spmBatch{3}.spm.stats.con.consess{4}.tcon.weights = [0,-1];
spmBatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

spmBatch{3}.spm.stats.con.consess{5}.fcon.name = 'Intercept';
spmBatch{3}.spm.stats.con.consess{5}.fcon.weights = 1;
spmBatch{3}.spm.stats.con.consess{5}.fcon.sessrep = 'none';

spmBatch{3}.spm.stats.con.consess{6}.fcon.name = 'zPnonc';
spmBatch{3}.spm.stats.con.consess{6}.fcon.weights = [0 1];
spmBatch{3}.spm.stats.con.consess{6}.fcon.sessrep = 'none';

%% Run th job
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

return