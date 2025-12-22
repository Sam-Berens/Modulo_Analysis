function [targetPath] = norm2MNI(G,subjectId,imagePaths)
%% Takes a cell of full filepaths for a set of images for a subject,...
% norms them to using the flowfield to group space G, 
% and then replaces their affine matrix with one which
% aligns the voxels to mni space 

%% NOTE: it will save all normed images inside a G folder within the source
%folder of the first image specified in the input!

% Set pathData
pathData = '../../Data';

% Set the template path
pathTemplates = [...
    pathData,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];

% Set filenameTemplate
filenameTemplate = dir([pathTemplates,filesep,'*_Template_6.nii']);
filenameTemplate = [filenameTemplate.folder,filesep,filenameTemplate.name];

% Populate filename for flow-field
    pathSubject = [pathData,filesep,subjectId];
    pathStructural = [pathSubject,filesep,'Structural'];
    pathFF = [pathStructural,filesep,G];
    dirFF = dir([pathFF,filesep,'u_rc1_*.nii']);
    fnFF = [dirFF.folder,filesep,dirFF.name];
% make a new G directory to move the normed results to inside src folder
    [imFolders,imNames] = cellfun(@(x) fileparts(x), imagePaths, 'UniformOutput',false); 
    targetPath = [imFolders{1},filesep,G];
    if ~exist("targetPath","dir")
        mkdir(targetPath);
    end 

% Create and execute the SPM batch
% SpmBatch = [{} {} {}];
SpmBatch{1}.spm.tools.dartel.mni_norm.template = {filenameTemplate};
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {fnFF};
SpmBatch{1}.spm.tools.dartel.mni_norm.data.subj.images = imagePaths; %If it goes wrong it might be that imageFns needs to be transposed?
SpmBatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
SpmBatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
SpmBatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);

%Move new images into the subject's group folder created inside the first
%src image's folder

srcPaths = cellfun(@(x,y) [x,filesep,'w',y,'.nii'],imFolders,imNames,'UniformOutput',false);
for ii=1:numel(imagePaths)
    movefile(srcPaths{ii},targetPath);
end

return
