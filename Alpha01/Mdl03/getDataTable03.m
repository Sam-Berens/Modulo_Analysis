function [dt03,mmY] = getDataTable03()
dtName = 'DataTable03.mat';
if exist(dtName,"file")
    strct = load(dtName);
    dt03 = strct.dt03;
    mmY = strct.mmY;
    return
end

% Cd out
wd = pwd;
cd ..;

G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%get MNI masks
dirs.Data = '../../Data';
hName = {'l'; 'r'}; 

roi = struct('V',[],'M',[],'idx',[]);
for iHem=1:2
    cHemi = hName{iHem};
    cDir = dir([fullfile(dirs.Data,'_Group','MniRois'),...
        filesep,'*',[cHemi,'HippC*.nii']]);
    fName = fullfile(cDir.folder,cDir.name);
    roi(iHem) = getIm(fName);
    [x,y,z] = ind2sub(size(roi(iHem).M),roi(iHem).idx);
    voxCoords = [x y z ones(numel(x),1)]';   % 4 × N
    mmCoords = roi(iHem).V.mat * voxCoords;
    mmY.(cHemi) = mmCoords(2,:)';
end

invW_mmY = getInvWmniYs();


%make a seperate table for coloc =- 1 and coloc =+1 
% same for each hemisphere and then join them at the end 

% loop through subject
fPart = fullfile('Analysis','Alpha01','Mdl03');
for iColoc = -1:2:1
    imIdx = (iColoc >0) + 1;
    for iHem=-1:2:1
        hemIdx = ((iHem >0) + 1);
        cHemi = hName{hemIdx};
        coLocation = ones(nSubs,1)*iColoc;
        hemisphere = ones(nSubs,1)*iHem;
        corr_method1 = nan(nSubs,1);
        corr_method2 = nan(nSubs,1);
        for iSubject=1:nSubs
            subjectId = subjectIds(iSubject);
            %% Method 1: corr on each sub, using mni coords and mni-normed epsilons
            %load in epsilon and mask out with roi
            cFldr = fullfile(dirs.Data,char(subjectId),fPart);
            fName = fullfile(cFldr,G,'wepsilon_new.nii');
            strct = getIm(fName);
            cE = strct(imIdx).M;
            %mask out the values of epsilon which did not meet the inequality...
            %nans from native image were turned to zeros in normed image
            fName = fullfile(cFldr,G,'wepMask_new.nii');
            strct = getIm(fName);
            eMask = strct(imIdx);
            eMask = eMask.M > 0.5;
            cE(eMask<1)= nan;
            %filter epsilon values for inside the mni roi 
            roiE = cE(roi(hemIdx).idx);

            if  sum(~isnan(roiE)) < 3 % correlation is invalid in this context
                corr_method1(iSubject) = nan;
            else
                corr_method1(iSubject) = corr(mmY.(cHemi),...
                    roiE, 'Rows','complete');
            end

            %% Method 2: corr using native epsilon and mni coords inverse
            fName = fullfile(cFldr,'epsilon_new.nii');
            strct = getIm(fName);
            %select out the M and V for the current coloc condition
            cE = strct(imIdx);
            roiId = [cHemi,'HippC'];
            %return epsilon and y image values inside the epi-res native roi
            cInvY = invW_mmY(iSubject);
            [roiInvY,roiE] = getMskdIms_EpiRes(G,subjectId,roiId,cInvY,cE);
            if  sum(~isnan(roiE)) < 3
                corr_method2(iSubject) = nan;
            else
                corr_method2(iSubject) = corr(roiInvY,roiE,...
                    'Rows','complete');
            end
        end
        t = table(subjectIds,hemisphere,coLocation,corr_method1,corr_method2);
        if iHem==-1 && iColoc==-1
            dt03 = t;
        else 
            dt03 = [dt03;t];%#ok<AGROW>
        end 
    end
end


%Fischer Z-transform
dt03.zCorr_m1 = atanh(dt03.corr_method1);
dt03.zCorr_m2 = atanh(dt03.corr_method2);

%get performance on nonCom
dtB = get_pNonc(G);
dtB.mcPnonc = dtB.pNonc - mean( dtB.pNonc);
dt03 = renamevars(dt03,"subjectIds","subjectId");
dt03 = join(dt03,dtB);

% Cd back and save
cd(wd);
save(dtName,"dt03","mmY");
return

function [im] = getIm(fname)
V = spm_vol(fname);
for ii=1:numel(V)
    im(ii).V = V(ii);
    im(ii).M = spm_read_vols(im(ii).V);
    im(ii).idx = find(im(ii).M>0.5);
end
return



