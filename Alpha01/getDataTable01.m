function [dt01,mmY] = getDataTable01()
dtName = 'DataTable01.mat';
if exist(dtName,"file")
    strct = load(dtName);
    dt01 = strct.dt01;
    mmY = strct.mmY;
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
    [x,y,z] = ind2sub(size(roi.M{iHem}),roi.idx{iHem});
    voxCoords = [x y z ones(numel(x),1)]';   % 4 × N
    mmCoords = roi.V{iHem}.mat * voxCoords;
    mmY.(cHemi) = mmCoords(2,:)';
end
mmY.comb = [mmY.l; mmY.r];

%preallocate columns for table
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
        cLq = getIm(char(subjectId),'wlQ.nii').M;
        %mask out the values of q which did not meet the inequality...
        %(mask needs thresholding because its been normed to mni)
        qMask =  getIm(char(subjectId),'wqMask.nii').M > 0.5;
        cLq(qMask<1)= nan;
        %store log(q) values from inside the mni roi for that subject
        roiLq = cLq(roi.idx{iHem});
        lq.(cHemi)(iSubject)  = {roiLq};
        %do corr on each sub
        %ignoring the q values which have been masked with nans
        rho.(cHemi)(iSubject) = corr(mmY.(cHemi),roiLq, 'Rows','complete'); 
    end
end

%make combined cols
lq.comb = arrayfun(@(l,r) {[l{:};r{:}]}, lq.l,lq.r);
rho.comb = arrayfun(@(a) corr(a{:}, mmY.comb,'Rows', 'complete'), lq.comb);
%Fischer Z-transform
z_rho.l = atanh(rho.l);
z_rho.r = atanh(rho.r);
z_rho.comb = atanh(rho.comb);

dt01 = table(subjectIds,z_rho.comb, z_rho.l,z_rho.r,...
    rho.comb,rho.l,rho.r,lq.comb,lq.l,lq.r);
dt01.Properties.VariableNames = {'subjectId','zCorr_comb','zCorr_l','zCorr_r',...
    'Corr_comb','Corr_l','Corr_r','logQ_comb','logQ_l','logQ_r'};
%get performance on nonCom
dt02 = get_pNonc(G);
dt02.mcPnonc = dt02.pNonc - mean( dt02.pNonc);
dt01 = join(dt01,dt02);

save(dtName,"dt01","mmY");
return

function [im] = getIm(subjectId,fname)
dirs.Data = '../../Data';
qFname = fullfile(dirs.Data,subjectId,'Analysis','Alpha01','Mdl01','G1',fname);
im.V = spm_vol(qFname);
im.M = spm_read_vols(im.V);
return
