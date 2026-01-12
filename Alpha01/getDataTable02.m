function [dt2] = getDataTable02()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%preallocate columns for table
y.l = nan(nSubs);
y.r = nan(nSubs);
q.l = nan(nSubs);
q.r = nan(nSubs);

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
end

% loop through subject
for iSubject=1:nSubs
    for iHem=1:2
        cHemi = hemi{iHem};
        roiId = [cHemi,'HippC']; %this is the combined head body and tail
        %load in q and mask out with roi
        q.(cHemi)(iSubject)  = (getNormdQ(char(subjectId))).*roi.M{iHem};
        %flatten q values

        %TO DO - obtain list of mm y coords for each q values
        roiIdx = find(roi.M{iHem}>0.5);
        [~,y,~] = ind2sub(size(roi.M{iHem}),roiIdx);
        y.(cHemi)(iSubject) = y .* roi.V{iHem}.mat;
    end
end

dt2 = table(subjectIds,y.l,q.l,y.r,q.r);
save('Datatable02.mat',"dt2");

return

function [q] = getNormdQ(subjectId)
dirs.Data = '../../Data';
qFname = fullfile(dirs.Data,subjectId,'Analysis','Alpha01','Mdl02','G1','wq.nii');
q.V = spm_vol(qFname);
q.M = spm_read_vols(q.V);
return
