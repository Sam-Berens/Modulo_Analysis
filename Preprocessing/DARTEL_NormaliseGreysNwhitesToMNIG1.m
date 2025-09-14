function [] = DARTEL_NormaliseGreysNwhitesToMNIG1()

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
g1Dir = [dataDir,filesep,'_Group',filesep,'G1'];
g1StrDir = [g1Dir,filesep,'Structural'];

%find group template
Dartel_Template = dir([g1StrDir,filesep,'*G1_Template_6.nii']);
Dartel_Template = [Dartel_Template.folder,filesep,Dartel_Template.name];

%get subject list
subjectIds = getG1();
nSubs = numel(subjectIds);

%collate all T1s and flowfields
allT1s = cell(nSubs,1);
allC1s = cell(nSubs,1);
allC2s = cell(nSubs,1);
allFlFds = cell(nSubs,1);

for iSubject=1:nSubs
    subjDir = [dataDir,filesep,subjectIds{iSubject}];
    struDir = [subjDir,filesep,'Structural'];

    %populate T1 and tissue filenames
    srcT1 = dir([struDir,filesep,'m_*_T1*.nii']); %want the bias corrected T1s
    allT1s{iSubject,1} = [srcT1.folder,filesep,srcT1.name];
    src1 = dir([struDir,filesep,'rc1_*.nii']);
    allC1s{iSubject,1} = [src1.folder,filesep,src1.name];
    src2 = dir([struDir,filesep,'rc2_*.nii']);
    allC2s{iSubject,1} = [src2.folder,filesep,src2.name];
    %populate flowfields
    flowfield = 'u_rc1_*.nii';
    srcFf = dir([struDir,filesep,flowfield]);
    allFlFds{iSubject,1} = [srcFf.folder,filesep,srcFf.name];
end

%% Create SpmBatch and hand over to JobMan.
SpmBatch = [{} {} {}]; %what does this bit do??
SpmBatch{1}.spm.tools.dartel.mni_norm.template = {Dartel_Template};
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = allFlFds;  %these get combined with the group->MNI affine transform
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {allT1s,allC1s,allC2s};
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN]; %will default to the resolution of the template
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN]; % No bounding box extension needed
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0; %preserve concentration as not doing VBM
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0]; %aliasing is fine

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);


%% Loop to move new images:
%file location:
% warped T1s and chained deformation .mats should both be saved in the
% strucutural folder of each subject
for iSubject = 1:nSubs
    subjectDir = [dataDir,filesep,subjectIds{iSubject}];
    struDir = [subjectDir,filesep,'Structural'];
    g1Dir = [struDir,filesep,'G1'];
    % movefile('w*',g1Dir);
    movefile([struDir,filesep,'w*'],g1Dir)
    movefile([struDir,filesep,'u_rc1*.nii'],g1Dir)

end

return