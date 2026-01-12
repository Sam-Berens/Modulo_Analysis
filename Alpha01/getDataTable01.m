function [dt1] = getDataTable01()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%get MNI masks
dirs.Data = '../../Data';
rois = dir([fullfile(dirs.Data,'_Group','MniRois'),filesep,'*HippC*.nii']);
hemi = {'l'; 'r'}; %NOTE l=1 r=2
for iHem=1:2
    cHemi = hemi{iHem};
    cDir = dir([fullfile(dirs.Data,'_Group','MniRois'),...
        filesep,'*',[cHemi,'HippC*.nii']]);
    fName = [cDir.folder,filesep,cDir.name];
    roi.V{iHem} = spm_vol(fName);
    roi.M{iHem} = spm_read_vols(roi.V{iHem});
    roi.idx{iHem} = find(roi.M{iHem}>0.5);
end

%preallocate columns for table
nVx.l = numel(roi.idx{1});
nVx.r = numel(roi.idx{2});
y.l = nan(nSubs,nVx.l);
y.r = nan(nSubs,nVx.r);
q.l = nan(nSubs,nVx.l);
q.r = nan(nSubs,nVx.r);

%to do : need to filter out the idx that there is no Q for?


% loop through subject
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject);
    for iHem=1:2
        cHemi = hemi{iHem};
        roiId = [cHemi,'HippC']; %this is the combined head body and tail
        %load in q and mask out with roi
        cq = getNormdQ(char(subjectId)).M;
        q.(cHemi)(iSubject,:)  = cq(roi.idx{iHem});
        %TO DO - obtain list of mm y coords for each q values
        [~,y,~] = ind2sub(size(roi.M{iHem}),roi.idx{iHem});
        y.(cHemi)(iSubject,:) = y .* roi.V{iHem}.mat;
        %do corr on each sub
    end
end

% dt1 = table(subjectIds,corVals);
%save('Datatable01.mat',"dt1");

return

function [q] = getNormdQ(subjectId)
dirs.Data = '../../Data';
qFname = fullfile(dirs.Data,subjectId,'Analysis','Alpha01','Mdl01','G1','wq.nii');
q.V = spm_vol(qFname);
q.M = spm_read_vols(q.V);
return
