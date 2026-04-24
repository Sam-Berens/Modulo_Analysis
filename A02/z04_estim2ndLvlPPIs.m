function [] = z04_estim2ndLvlPPIs(G)
% Pick which contrast numbers you want to estimate models for and build a
% vector of them
% cons 4-7 are the two PPI effects, the main effect of PPI (irrespective of
% a or b condition) and the difference between the two PPI effects 
for ii=4:7
    estim2ndLvl(ii);
end

return


%TO DO - finish below!!

function [] = estim2ndLvl(con,G)
%construct scans cell aray doing like 
subjectIds = sort(getSubjectIds(G));
nSubs = numel(subjectIds);

dirs = arrayfun(@(x) dirFnc(x), subjectIds,'UniformOutput',false);
scans = cellfun(@(x) sprintf('%s%s%0d3.nii',x,...
    filesep,con),dirs.PPI,'UniformOutput',false);
%get covariate
pNonc = get_pNonc(G);
pNonc.zpNonc = zscore(pNonc.pNonc);
%make sure its ordered correctly
pNonc.Properties.RowNames = cellstr([pNonc.subjectId]');
pNonc = sortrows(pNonc);
zpNonc = pNonc{:,'zpNonc'};

spmBatch{1}.spm.stats.factorial_design.dir = {'/mnt/Erebus/Modulo/Data/_Group/G1/Analysis'};
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
spmBatch{2}.spm.stats.fmri_est.spmmat = {[dirs.Output,filesep,'SPM.mat']};
spmBatch{2}.spm.stats.fmri_est.write_residuals = 0;
spmBatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Run th job
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

return


function [dirs] = dirFnc(subjectId)
subjectId = char(subjectId);
%get subject-specific filenames and info
dirs.Subject = fullfile(dirs_Data,subjectId);
dirs.A02 = fullfile(dirs.Subject,'Analysis','A02');
dirs.PPI = fullfile(dirs.A02,sprintf('PPI_%s',roiInfo.id));
return