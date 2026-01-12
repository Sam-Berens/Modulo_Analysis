function [dt01] = getDataTable01()
dtName = 'Datatable01.mat';
if exist("dtName","dir")
    strct = load(dtName);
    dt01 = strct.dt01;
    return
end

G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%get MNI masks
dirs.Data = '../../Data';
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
mmY.l = cell(nSubs,1);
mmY.r = cell(nSubs,1);
lq.l = cell(nSubs,1);
lq.r = cell(nSubs,1);
rho.l = nan(nSubs,1);
rho.r = nan(nSubs,1);

% loop through subject
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject);
    for iHem=1:2
        cHemi = hemi{iHem};
        %load in log(q) and mask out with roi
        clq = getIm(char(subjectId),'wlQ.nii').M;
        %mask out the values of q which did not meet the inequality...
        %(mask needs thresholding because its been normed to mni)
        qMask =  getIm(char(subjectId),'wqMask.nii').M > 0.5;
        clq(qMask<1)= nan;
        %store log(q) values from inside the mni roi for that subject
        roilq = clq(roi.idx{iHem});
        lq.(cHemi)(iSubject)  = {roilq};
        %TO DO - obtain list of mm y coords for each q values
        [x,y,z] = ind2sub(size(roi.M{iHem}),roi.idx{iHem});
        voxCoords = [x y z ones(numel(x),1)]';   % 4 × N
        mmCoords = roi.V{iHem}.mat * voxCoords;
        mmY.(cHemi)(iSubject) = {mmCoords(2,:)'};
        %do corr on each sub
        %ignoring the q values which have been masked with nans
        rho.(cHemi)(iSubject) = corr(mmCoords(2,:)',roilq, 'Rows', 'complete');
    end
end

z_rho.l = atanh(rho.l);
z_rho.r = atanh(rho.r);

dt01 = table(subjectIds,mmY.l,lq.l,rho.l,z_rho.l,mmY.r,lq.r,rho.r,z_rho.r);
dt01.Properties.VariableNames = {'subjectId','left_mmY','left_logQ',...
    'left_Corr','left_zCorr','right_mmY','right_logQ',...
    'right_Corr','right_zCorr'};
save(dtName,"dt01");

return

function [im] = getIm(subjectId,fname)
dirs.Data = '../../Data';
qFname = fullfile(dirs.Data,subjectId,'Analysis','Alpha01','Mdl01','G1',fname);
im.V = spm_vol(qFname);
im.M = spm_read_vols(im.V);
return
