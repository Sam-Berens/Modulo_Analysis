function [] = z03_MakeFmMag(subjectId)

preprocDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);
cd(subjectId);
cd('Fieldmap');

%% Load the raw fieldmap images
fileList = [dir('_*_FmAP-*.nii');dir('_*_FmPA-*.nii')];
filenames = {fileList.name}';
V = spm_vol(filenames);
M = cellfun(@(s)spm_read_vols(s),V,'UniformOutput',false);

%% Load the movement parameters from topup
moveParFn = dir('*_FmTopupCoefs-Movpar.txt');
T = readtable(moveParFn.name);
T.Properties.VariableNames = {'tX','tY','tZ','rX','rY','rZ'};
R = T.Variables;

%% Move all the images (by adjusting the direction cosines matrix)
newV = V;
for ii = 1:numel(newV)
    newV{ii}.fname = ['r',newV{ii}.fname];
    newV{ii}.mat = spm_matrix(R(ii,:))*V{ii}.mat;
    spm_write_vol(newV{ii},M{ii});
end

%% Set-up and run the SPM batch
SpmBatch{1}.spm.util.imcalc.input = cellfun(@(s)s.fname,newV,...
    'UniformOutput',false);
SpmBatch{1}.spm.util.imcalc.output = sprintf('_%s_FmMag.nii',subjectId);
SpmBatch{1}.spm.util.imcalc.outdir = {''};
SpmBatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8)/8';
SpmBatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
SpmBatch{1}.spm.util.imcalc.options.dmtx = 0;
SpmBatch{1}.spm.util.imcalc.options.mask = 0;
SpmBatch{1}.spm.util.imcalc.options.interp = -4;
SpmBatch{1}.spm.util.imcalc.options.dtype = 512;

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%% Cd out of subject and then out of data
cd(preprocDir);

return