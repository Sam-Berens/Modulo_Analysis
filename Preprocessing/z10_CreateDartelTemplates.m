function [] = z10_CreateDartelTemplates(G)

% Get the subjectIds
subjectIds = getSubjectIds(G);
nSubjects = numel(subjectIds);

% Preallocate the DARTEL-imported TPMs
FileNames_rc1 = cell(nSubjects,1);
FileNames_rc2 = cell(nSubjects,1);
FileNames_rc3 = cell(nSubjects,1);

% Set pathData
pathData = '../../Data';

% Set the template output path
pathTemplates = [...
    pathData,...
    filesep,'_Group',...
    filesep,G,...
    filesep,'Structural',...
    filesep,'DARTEL_Templates'];
if ~exist(pathTemplates,'dir')
    mkdir(pathTemplates);
else
    error('It appears that DARTEL templates for the requested group ID already exist!');
end

% Subject loop to populate inputs and create folder for the flow-fields
for iSubject = 1:nSubjects

    pathSubject = [pathData,filesep,subjectIds{iSubject}];
    pathStructural = [pathSubject,filesep,'Structural'];

    % Populate FileNames_rc*:
    Src1 = dir([pathStructural,filesep,'rc1_*.nii']);
    FileNames_rc1{iSubject,1} = [Src1.folder,filesep,Src1.name];
    Src2 = dir([pathStructural,filesep,'rc2_*.nii']);
    FileNames_rc2{iSubject,1} = [Src2.folder,filesep,Src2.name];
    Src3 = dir([pathStructural,filesep,'rc3_*.nii']);
    FileNames_rc3{iSubject,1} = [Src3.folder,filesep,Src3.name];

    % Make the named group folder with in the subject's Structural dir
    pathFlowFields = [pathStructural,filesep,G];
    if ~exist(pathFlowFields,'dir')
        mkdir(pathFlowFields);
    else
        error('It appears that this group ID has already been used!');
    end

end

% Create and execute the SPM batch
SpmBatch = {};
SpmBatch{1}.spm.tools.dartel.warp.images = {FileNames_rc1, FileNames_rc2, FileNames_rc3};
SpmBatch{1}.spm.tools.dartel.warp.settings.template = [G,'_Template'];
SpmBatch{1}.spm.tools.dartel.warp.settings.rform = 0;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
SpmBatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
SpmBatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
SpmBatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 8;
SpmBatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);


% Move the flow-fields images and the Templates
movedTemplate = false;
for iSubject = 1:nSubjects

    pathSubject = [pathData,filesep,subjectIds{iSubject}];
    pathStructural = [pathSubject,filesep,'Structural'];
    pathFlowFields = [pathStructural,filesep,G];

    % Flow-fields
    movefile([pathStructural,filesep,'u_rc1_*.nii'],...
        pathFlowFields);

    % Templates
    dirTemplates = dir([pathStructural,filesep,'*_Template_*']);
    if numel(dirTemplates) > 0 && ~movedTemplate
        move([pathStructural,filesep,'*_Template_*'],pathTemplates);
        movedTemplate = true;
    elseif numel(dirTemplates) > 0
        error('Something bad has happened!');
    end

end

return