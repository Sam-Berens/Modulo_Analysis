function [] = z10_InverseWarpRois()
% Create inverse normalised versions of some image(s). The image that is inverse-normalised should
% be in alignment with the template (generated during the warping procedure). Note that the results
% have the same dimensions as the “flow fields”, but are mapped to the original images via the affine
% transformations in their headers.

%ergo in order to warp our MNI ROIs to native space we need to get them in
%alignment with the group template before using the dartel.inverseWarp
%procedure (using each subjects flowfields)


%to do this we could just do coregistration on the ROIs with the template
%as the reference image? - or we could use the formula below

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
mRoiDir = [dataDir,filesep,'_Group',filesep,'MniRois'];
drtlDir = [dataDir,filesep,'_Group',filesep,'G1',filesep,'Structural', filesep,'DARTEL_templates'];
grpRoiDir = [drtlDir, filesep, 'GroupSpaceRois'];

if ~exist(grpRoiDir, "dir")
    mkdir(grpRoiDir)
end

%find group template
srcDartelTemplate = dir([drtlDir,filesep,'*G1_Template_6.nii']);
dartelTemplate = [srcDartelTemplate.folder,filesep,srcDartelTemplate.name];


subjectIds = getG1();
nSubs = numel(subjectIds);
nativeDirs = cell(nSubs);
for ii=1:nSubs
    subjDir = [dataDir,filesep,subjectIds{ii}];
    struDir = [subjDir,filesep,'Structural'];
    nativeDirs{ii} = [struDir,filesep,'NativeRois'];
    if ~exist(nativeDirs{ii},'dir')
        mkdir([nativeDirs{ii},filesep,'w'])
        mkdir([nativeDirs{ii},filesep,'rw']) %for later sliced rois to go in
        mkdir([nativeDirs{ii},filesep,'mrw']) %for later masked rois to go in
    end
end

%M is the mapping from the voxel grid of the ROI to group space
% it is calculated as the path from
% roiVoxelGrid->MNIspace(M1) to...
% MNIspace->groupTemplateVoxelGrid(M3^-1) to...
% groupTemplateVoxelGrid->groupSpace(M2)

%it is neccesary to * by M2 in order to make M map the ROI voxel grid map
%to group space, otherwise M1*inv(M3) would just tell spm to align the ROIs
%according to the voxel grid (e.g. the intensity of voxel [1,1,1](which
%should be black, hailing from the bottom left corner below your head)being
%presented near the origin of the image.

%the reason that these .mat files are important despite us not displaying
%the images is that the volumes need to be aligned to group space for the
%flowfields to work.

% for formula discussion see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1203&L=spm&P=R89160&1=spm&9=A&J=on&d=No+Match%3BMatch%3BMatches&z=4

%load the native2group2MNI .mat file (M3)
affineFn = [drtlDir,filesep,'G1_Template_6_2mni.mat'];

%load in the rois (for M1)
mRoisrc = dir([mRoiDir,filesep,'_*.nii']);
mRois = cellfun(@(x) x,{mRoisrc.name}','UniformOutput',false);

%M2 is the mapping of the group template voxel grid to group space
M2 = spm_get_space(dartelTemplate);
%M3 is the mapping of the group template voxel grid to MNI space
load(affineFn);
M3 = mni.affine;

%make initial group-space ROIS
for ii=1:numel(mRois)
    roiFilename = [mRoiDir,filesep,mRois{ii}];
    gspaceRoiFn = [grpRoiDir,filesep,'w',mRois{ii}]; 
    copyfile(roiFilename,gspaceRoiFn);
%M1 is the mapping of the roi voxel grid to MNI space
    M1 = spm_get_space(roiFilename);
%M is the mapping from the voxel grid of the ROI to group space   
    M = M2*inv(M3)*M1;
%write the affine .mat to the group-space copy of the given roi volume
    spm_get_space(gspaceRoiFn,M)
end 

%% Individually warp the rois to native space for each person
%using their flowfield

%collect the grpRois for batch warping
srcGrpRoiDir = dir([grpRoiDir,filesep,'w_*.nii']);% should all be named this after creation above
grpRois = cellfun(@(x,y) [x,filesep,y],{srcGrpRoiDir.folder}',{srcGrpRoiDir.name}','UniformOutput',false);

for ii=1:nSubs
    subjDir = [dataDir,filesep,subjectIds{ii}];
    sStruDir = [subjDir,filesep,'Structural'];
    sG1Dir = [sStruDir,filesep,'G1'];
    srcFf = dir([sG1Dir,filesep,'u_rc1_*.nii']);
    flowField = [srcFf.folder,filesep,srcFf.name];

    spmBatch{1}= {};
    spmBatch{1}.spm.tools.dartel.crt_iwarped.flowfields = {flowField};
    spmBatch{1}.spm.tools.dartel.crt_iwarped.images = grpRois;
    spmBatch{1}.spm.tools.dartel.crt_iwarped.K = 6; %default
    spmBatch{1}.spm.tools.dartel.crt_iwarped.interp = 7; %putting it on highest for accuracy

    spm_jobman('initcfg');
    spm_jobman('run',spmBatch);

    targetfld =  [nativeDirs{ii}, filesep, 'w'];
    %move all the native-warped ROIs for that subject from where they get sent(their G1 folder because thats where the flowfield is?) to their native ROI folder
    movefile([sG1Dir,filesep,'ww_*.nii'],targetfld)
end


%instead of reslicing ROIs into EPI resolution now, we will leave these unbinarised
% and will toggle a threshold for inclusion in the final analaysis rois
% masks in a later script (will need to use interpolation but can do this
% with spm function e.g. reslice?


