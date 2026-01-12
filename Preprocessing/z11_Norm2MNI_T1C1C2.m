function [] = z11_Norm2MNI_T1C1C2(G)

% Get the subjectIds
subjectIds = getSubjectIds(G);
nSubjects = numel(subjectIds);

% Set pathData
pathData = '../../Data';

% Set the template output path
pathTemplates = [...
    pathData,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];

% Set filenameTemplate
filenameTemplate = dir([pathTemplates,filesep,'*_Template_6.nii']);
filenameTemplate = [filenameTemplate.folder,filesep,filenameTemplate.name];

% Preallocate filename arrays for flow-fields and images to norm
fnsFFs = cell(nSubjects,1);
fnsT1s = cell(nSubjects,1);
fnsC1s = cell(nSubjects,1);
fnsC2s = cell(nSubjects,1);

% Populate filename arrays for flow-fields and images to norm
for iSubject = 1:nSubjects
    pathSubject = [pathData,filesep,char(subjectIds(iSubject))];
    pathStructural = [pathSubject,filesep,'Structural'];
    pathFF = [pathStructural,filesep,G];

    dirFF = dir([pathFF,filesep,'u_rc1_*.nii']);
    fnsFFs{iSubject} = [dirFF.folder,filesep,dirFF.name];

    dirT1 = dir([pathStructural,filesep,'m_*_T1*.nii']);
    fnsT1s{iSubject} = [dirT1.folder,filesep,dirT1.name];

    dirC1 = dir([pathStructural,filesep,'c1_*.nii']);
    fnsC1s{iSubject} = [dirC1.folder,filesep,dirC1.name];

    dirC2 = dir([pathStructural,filesep,'c2_*.nii']);
    fnsC2s{iSubject} = [dirC2.folder,filesep,dirC2.name];
end

% Create and execute the SPM batch
SpmBatch = {};
SpmBatch{1}.spm.tools.dartel.mni_norm.template = {filenameTemplate};
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = fnsFFs;
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {fnsT1s,fnsC1s,fnsC2s};
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

% Loop to move new images
for iSubject = 1:nSubjects
    pathSubject = [pathData,filesep,char(subjectIds(iSubject))];
    pathStructural = [pathSubject,filesep,'Structural'];
    pathFF = [pathStructural,filesep,G];
    movefile([pathStructural,filesep,'w*.nii'],pathFF);
end

return