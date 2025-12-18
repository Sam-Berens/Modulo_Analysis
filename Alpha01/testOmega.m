function [l,r] = testOmega()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
Omega.l = nan(nSubs,1);
Omega.r = nan(nSubs,1);
% loop through subject
for iSubject=1:nSubs
    hemi = {'l'; 'r'};
    for iHem=1:2
        cHemi = hemi{iHem};
        roiId = [cHemi,'HippC']; %this is the combined head body and tail
        makeQ(G,subjectId,roiId);
        %get filenames of files to warp
    end 
        %this is both right and left
        qFname = getQfileList(subjectId);
        %warp to MNI space
        targetFolder = norm2MNI(G,subjectId,qFname);
        %load q image for left and right
        qDir = dir(targetFolder,filesep,'wq*.nii');
        qPaths = cellfun(@(x) fullfile(x.folder,x.name),qDir,'UniformOutput',false);
        l = qPaths(contains(qPaths,'r'));
        r = qPaths(contains(qPaths,'l')); %this is probs dangerous if the roi name changes
        lV =  spm_vol(l);
        lM = spm_read_vols(rl);
        rV = spm_vol(r);
        rM = spm_read_vols(rV);

        %TO DO load and apply roi masks!


        %then i guess you can collect up an array of q vals and Y coords for
        %each hemisphere
        ql = lM(mniRoi.idx);
        yl = mniRoi.idx * lV.mat %% TO DO this is not how to do it just place holder !! 
        qr = rM(mniRoi.idx);
        yr = mniRoi.idx * rV.mat
        Omega.l(iSubject) = corr(ql,yl);
        Omega.r(iSubject) = corr(qr,yr);
end

[l.h,l.p,l.ci] = ttest(Omega.l);
[r.h,r.p,r.ci] = ttest(Omega.l);
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
