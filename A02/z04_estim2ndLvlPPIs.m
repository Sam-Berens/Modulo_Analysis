function [] = z04_estim2ndLvlPPIs(G,roiInfo)
%% Estimates a seperate 2nd level model each for 4 contrast, for a given roi 
% cons are the two PPI effects, the main effect of PPI (irrespective of
% a or b condition) and the difference between the two PPI effects

% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_rVisual.nii';
% roiInfo.id = 'rVisual';
% Pick which contrast numbers you want to estimate models for and build a
% vector of them 
for ii=4:7
    estim2ndLvl(G,roiInfo,ii);
end
return

function [] = estim2ndLvl(G,roiInfo,con)

subjectIds = sort(getSubjectIds(G));
dir_Data = ['..',filesep,'..',filesep,'Data'];
dir_G = fullfile(dir_Data,'_Group','G1');
epiMask = fullfile(dir_G,'Structural','GrpEpiMask00',...
    sprintf('%s_GrpEpiMask00.nii',G));
dir_outPut = fullfile(dir_Data,'_Group',G,'Analysis','A02',sprintf('con_%03d',con));
%check if A02 group dir exists, and then if the subdir for that contrast
%model exists before making them
if ~exist(dir_outPut(1:end-8),"dir")
    mkdir(dir_outPut(1:end-8))
end
if ~exist(dir_outPut,"dir")
    mkdir(dir_outPut)
end

%construct scans cell array 
dirs = arrayfun(@(x) dirFnc(x,dir_Data,roiInfo), subjectIds,'UniformOutput',false);
scans = cellfun(@(x) sprintf('%s%scon_%04d.nii',x.PPI,...
    filesep,con),dirs,'UniformOutput',false);
%get covariate
pNonc = get_pNonc(G);
pNonc.zpNonc = zscore(pNonc.pNonc);
%make sure its ordered correctly
pNonc.Properties.RowNames = cellstr([pNonc.subjectId]');
pNonc = sortrows(pNonc);
zpNonc = pNonc{:,'zpNonc'};

spmBatch{1}.spm.stats.factorial_design.dir = {dir_outPut};
spmBatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
spmBatch{1}.spm.stats.factorial_design.cov.c = zpNonc;
spmBatch{1}.spm.stats.factorial_design.cov.cname = 'zpNonc';
spmBatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
spmBatch{1}.spm.stats.factorial_design.cov.iCC = 5;
spmBatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

spmBatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
spmBatch{1}.spm.stats.factorial_design.masking.im = 1; %from the sounds of it this is just saying if its already nan in the image, take it to be part of the mask (what happens if nans arent shared across peopel?)
spmBatch{1}.spm.stats.factorial_design.masking.em = {epiMask}; %explicit masking?
spmBatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
spmBatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Estimation
spmBatch{2}.spm.stats.fmri_est.spmmat = {fullfile(dir_outPut,'SPM.mat')};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Make contrasts
conNames = cell(2,1);
%effect of intercept
conNames{1} = 'intercept';
H.intercept = 1;
%main effect of pNonc
conNames{2} = 'zPnonc';
H.zPnonc = [0,1];
fields = fieldnames(H);

spmBatch{3}.spm.stats.con.spmmat = {fullfile(dir_outPut,'SPM.mat')};
for iH = 1:numel(conNames)
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.name = conNames{iH};
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.weights = H.(fields{iH});
    spmBatch{3}.spm.stats.con.consess{iH}.tcon.sessrep = 'replsc';
end
spmBatch{3}.spm.stats.con.delete = 1;

% Run th job
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

return


function [dirs] = dirFnc(subjectId,dir_Data,roiInfo)
subjectId = char(subjectId);
dir_Data = char(dir_Data);
%get subject-specific filenames and info
dirs.Subject = fullfile(dir_Data,subjectId);
dirs.A02 = fullfile(dirs.Subject,'Analysis','A02');
dirs.PPI = fullfile(dirs.A02,sprintf('PPI_%s',roiInfo.id));
return