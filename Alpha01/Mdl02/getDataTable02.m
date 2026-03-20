function [DataTable02] = getDataTable02()

filename = 'DataTable02.mat';
if false%exist(filename,"file")
    X = load(filename);
    DataTable02 = X.DataTable02;
    return
end


%% Get MNI Hippocampal masks
dirs.Data = '../../../Data';

fileList = dir([fullfile(dirs.Data,'_Group','MniRois'), ...
    filesep,'*','lHippC*.nii']);
filename = fullfile(fileList.folder,fileList.name);
leftHippo = loadMask(filename);

fileList = dir([fullfile(dirs.Data,'_Group','MniRois'), ...
    filesep,'*','rHippC*.nii']);
filename = fullfile(fileList.folder,fileList.name);
rightHippo = loadMask(filename);


%%
hemiName = {'l'; 'r'};

hippoMasks = struct('V',[],'M',[],'idx',[]);
for iHemi = 1:2
    cHemi = hemiName{iHemi};
    cDir = dir([fullfile(dirs.Data,'_Group','MniRois'),...
        filesep,'*',[cHemi,'HippC*.nii']]);
    fName = fullfile(cDir.folder,cDir.name);
    volume = loadMask(fName);
    [x,y,z] = ind2sub(volume.V.dim,volume.idx');
    VxIdx = [x; y; z; ones(size(x))];   % 4 × N
    Mni = volume.V.mat * VxIdx;
    mmY.(cHemi) = Mni(2,:)';
    hippoMasks(iHemi) = volume;
end

Y_left = nan(size(hippoMasks(1).M));
Y_left(hippoMasks(1).idx) = Mni(2,:);

Y_right = nan(size(hippoMasks(2).M));
Y_right(hippoMasks(2).idx) = Mni(2,:);

%%
G = 'G1';

% CD out
wd = pwd;
cd ..;

subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);

%make a seperate table for coloc =- 1 and coloc =+1
% same for each hemisphere and then join them at the end

% loop through subject
fPart = fullfile('Analysis','Alpha01','Mdl02');
for iColoc = -1:2:1
    imIdx = (iColoc >0) + 1;
    for iHemi=-1:2:1
        hemIdx = ((iHemi >0) + 1);
        cHemi = hemiName{hemIdx};
        coLocation = ones(nSubs,1)*iColoc;
        hemisphere = ones(nSubs,1)*iHemi;
        corr_method1 = nan(nSubs,1);
        for iSubject=1:nSubs
            subjectId = subjectIds(iSubject);
            %% Method 1: corr on each sub, using mni coords and mni-normed log(q)
            %load in log(q) and mask out with roi
            cFldr = fullfile(dirs.Data,char(subjectId),fPart);
            fName = fullfile(cFldr,G,'wlQ.nii');
            X = load4dVol(fName);
            clQ = X(imIdx).M;
            %mask out the values of q which did not meet the inequality...
            %(mask needs thresholding because its been normed to mni)
            fName = fullfile(cFldr,G,'wqMask.nii');
            X = load4dVol(fName);
            qMask = X(imIdx);
            qMask = qMask.M > 0.5;
            clQ(qMask<1)= nan;
            %store log(q) values from inside the mni roi for that subject
            roiLq = clQ(hippoMasks(hemIdx).idx);
            if sum(~isnan(roiLq)) < 3
                corr_method1(iSubject) = nan;
            else
                corr_method1(iSubject) = corr(mmY.(cHemi),roiLq, 'Rows','complete');
            end
        end
        t = table(subjectIds,hemisphere,coLocation,corr_method1);
        if iHemi==-1 && iColoc==-1
            DataTable02 = t;
        else
            DataTable02 = [DataTable02;t];%#ok<AGROW>
        end
    end
end


%Fischer Z-transform
DataTable02.zCorr_m1 = atanh(DataTable02.corr_method1);

%get performance on nonCom
dtB = get_pNonc(G);
dtB.mcPnonc = dtB.pNonc - mean( dtB.pNonc);
DataTable02 = renamevars(DataTable02,"subjectIds","subjectId");
DataTable02 = join(DataTable02,dtB);

% Cd back and save
cd(wd);
save(filename,"DataTable02","mmY");
return

function [mask] = loadMask(filename)
V = spm_vol(filename);
mask.V = V;
mask.M = spm_read_vols(mask.V);
mask.idx = find(mask.M > 0.5);
[x,y,z] = ind2sub(V.dim,mask.idx');
VxIdx = [x; y; z; ones(size(x))]; % 4 × N
Mni = V.mat * VxIdx;
MniY = nan(V.dim);
MniY(mask.idx) = Mni(2,:);
mask.MniY = MniY;
return

function [Vols] = load4dVol(filename)
V = spm_vol(filename);
nVol = numel(V);
Vols = struct;
for ii = 1:nVol
    Vols(ii).V = V;
    Vols(ii).M = spm_read_vols(Vols(ii).V);
end
return