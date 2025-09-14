function [] = DARTEL_CreateTemplatesG1()

dirLits = dir('../../Data');
g1Dir = [dirLits(1).folder,filesep,'_Group',filesep,'G1'];
g1StrDir = [g1Dir,filesep,'Structural'];
if ~exist(g1StrDir,'dir')
    mkdir(g1StrDir)
    drtlDir = [g1StrDir,filesep,'DARTEL_Templates'];
    mkdir(drtlDir)
end

%get subject list
subjectIds = getG1();
nSubs = numel(subjectIds);

%collate all rigid alignment tissue maps for white,grey and csf?
FileNames_rc1 = cell(nSubs,1);
FileNames_rc2 = cell(nSubs,1);
FileNames_rc3 = cell(nSubs,1);

for iSubject=1:nSubs

    dataDir = dir(['..',filesep,'..',filesep,'Data']);
    dataDir = dataDir(1).folder;
    subjDir = [dataDir,filesep,subjectIds{iSubject}];
    struDir = [subjDir,filesep,'Structural'];

    %% Populate FileNames_rc*:
    Src1 = dir([struDir,filesep,'rc1_*.nii']);
    FileNames_rc1{iSubject,1} = [Src1.folder,filesep,Src1.name];
    Src2 = dir([struDir,filesep,'rc2_*.nii']);
    FileNames_rc2{iSubject,1} = [Src2.folder,filesep,Src2.name];
    Src3 = dir([struDir,filesep,'rc3_*.nii']);
    FileNames_rc3{iSubject,1} = [Src3.folder,filesep,Src3.name];

    %make the g1 folder to put the flowfields in within each sub's dir
    g1Dir = [struDir,filesep,'G1'];
    if ~exist(g1Dir,'dir')
        mkdir(g1Dir)
    end

end

SpmBatch = {};
% SpmBatch{1}.spm.tools.dartel.warp.output.option = 'allin'; %the spm_dartel_warp job appears to have output as field so we should send the images here to the group folder to be safe
% SpmBatch{1}.spm.tools.dartel.warp.output.outDir = g1StrDir;
SpmBatch{1}.spm.tools.dartel.warp.images = {FileNames_rc1, FileNames_rc2, FileNames_rc3};
SpmBatch{1}.spm.tools.dartel.warp.settings.template = 'G1_Template'; %name of template to be saved - expect it to be appended to flowimages too
SpmBatch{1}.spm.tools.dartel.warp.settings.rform = 0;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
SpmBatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06]; % everything is same as  default
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
SpmBatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 8; % default is 3 but this is just a speed accuracy trade off and higher is more acc
SpmBatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

spm_jobman('initcfg');
spm_jobman('run',SpmBatch);


%find where the group templates /affine mat have gone and move them out
%that subjects folder
for iSubject=1:nSubs
    dataDir = dir(['..',filesep,'..',filesep,'Data']);
    dataDir = dataDir(1).folder;
    subjDir = [dataDir,filesep,subjectId];
    struDir = [subjDir,filesep,'Structural'];
    g1Dir = [struDir,filesep,'G1'];
    src = dir([g1Dir,filesep,'G1_Template_*']);
    if numel(src)>0
        movefile('G1_Template_*',drtlDir)
    end

end

return