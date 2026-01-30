function [dt03,mmY] = getDataTable03()
dtName = 'DataTable03.mat';
if exist(dtName,"file")
    strct = load(dtName);
    dt03 = strct.dt03;
    mmY = strct.mmY;
    return
end

G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%get MNI masks
dirs.Data = '../../Data';
hName = {'l'; 'r'}; 

for iHem=1:2
    cHemi = hName{iHem};
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
nHems = 2;
nColocs = 2;


%what we're going to do is make a seperate table for coloc =- 1 and coloc
%=+1 and then join them at the end so we're nto doing crazy indexing

% loop through subject
for iColoc = -1:2:1
    imIdx = (iColoc >0) + 1;
    for iHem=-1:2:1
        cHemi = hName{((cHemi >0) + 1)};
        epsilon = nan((nSubs*nHems*nColocs),1);
        coLocation = ones(nSubs,1)*iColoc;
        hemisphere = ones(nSubs,1)*iHem;
        rho = nan((nSubs*nHems*nColocs),1);
        for iSubject=1:nSubs
            subjectId = subjectIds(iSubject);
            %load in log(q) and mask out with roi
            cE = getIm(char(subjectId),'wepsilon.nii').M;
            cE = cE(:,:,:,imIdx);%select out the slice for the current coloc condition
            %mask out the values of q which did not meet the inequality...
            %(mask needs thresholding because its been normed to mni)
            eMask =  getIm(char(subjectId),'wepsilonMask.nii').M > 0.5;
            eMask = eMask(:,:,:,imIdx);
            cE(eMask<1)= nan;
            %store log(q) values from inside the mni roi for that subject
            roiE = cE(roi.idx{iHem});
            epsilon(iSubject)  = {roiE};
            %do corr on each sub
            %ignoring the epsilon values which have been masked with nans
            rho(iSubject) = corr(mmY.(cHemi),roiE, 'Rows','complete');
        end
        t = table(subjectIds,hemisphere,coLocation,epsilon,rho);
        if iHem==-1 && iColoc==-1
            dt03 = t;
        else 
            dt03 = [dt03;t];%#ok<AGROW>
        end 
    end
end


%Fischer Z-transform
dt03.zCorr = atanh(dt03.rho);

%get performance on nonCom
dtB = get_pNonc(G);
dtB.mcPnonc = dtB.pNonc - mean( dtB.pNonc);
dt03 = join(dt03,dtB);

save(dtName,"dt03","mmY");
return

function [im] = getIm(subjectId,fname)
dirs.Data = '../../Data';
fname = fullfile(dirs.Data,subjectId,...
    'Analysis','Alpha01','Mdl03','G1',fname);
im.V = spm_vol(fname);
im.M = spm_read_vols(im.V);
return
