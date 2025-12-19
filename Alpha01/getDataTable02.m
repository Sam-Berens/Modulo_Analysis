function [l,r] = getDataTable02()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
Omega.l = nan(nSubs,1);
Omega.r = nan(nSubs,1);
% loop through subject
for iSubject=1:nSubs
        qFname = getQfileList(subjectId);
        %warp to MNI space
        targetFolder = norm2MNI(G,subjectId,qFname);
        %load q image for left and right
        qDir = dir(targetFolder,filesep,'wq*.nii');
        qPaths = cellfun(@(x) fullfile(x.folder,x.name),qDir,'UniformOutput',false);
hemi = {'l'; 'r'};
    for iHem=1:2
        cHemi = hemi{iHem};
        roiId = [cHemi,'HippC']; %this is the combined head body and tail
        %this is both right and left
        %TO DO read in q image and mask
        
        Q.(cHemi) = M(mniRoi.idx);
        Y.(cHemi) = mniRoi.idx * V.mat%% TO DO this is not how to do it just place holder !! 
        % TO DO - 
        % zX0 - make q
        % zX1_... = make norm images loop into seperar function which is like
        % getDataTable01... =  %TO DO load and apply roi masks! (using roi
        % filtering helper function)(this is where you do left and right
        % seperately ), also collecting up qs and ys and testing their corr
        % put that in the datatable and then ttest on datatable across
        % subjects
    end

        Omega.l(iSubject) = corr(ql,yl);
        Omega.r(iSubject) = corr(qr,yr);
   
end
%TO DO BUILD OMEGA INTO TABLE INSTEAD OF STRUCT
return

%% TO DO: 

%ADD to the markdown for the analysis folder explaining that mdl00 is
%prereg, mdl01 is searchlight version of mdl00 and mdl02 is this here. 
%inside mdl00.mlx file use .1 etc to denote added terms


function [qFnames] = getQfileList(subjectId)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl02 = [dirs.Alpha01,filesep,'Mdl02'];
qs = dir([dirs.Mdl02,filesep,'q*.nii']);
qFnames = cellfun(@(x)[x.folder,filesep,x.name],qs,'UniformOutput',false);
return
