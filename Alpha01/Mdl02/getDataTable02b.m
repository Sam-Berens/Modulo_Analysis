function [dt02b,mmY] = getDataTable02b()
dtName = 'DataTable02b.mat';
if exist(dtName,"file")
    strct = load(dtName);
    dt02b = strct.DataTable02b;
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
hemiName = {'l'; 'r'};
% hippoMasks = struct('V',[],'M',[],'idx',[]);
for iHemi = 1:2
    cHemi = hemiName{iHemi};
    cDir = dir([fullfile(dirs.Data,'_Group','MniRois'),...
        filesep,'*',[cHemi,'HippC*.nii']]);
    fName = fullfile(cDir.folder,cDir.name);
    mask = loadMask(fName);
    %these are an output because needed for scatter graphs
    mmY.(cHemi) = mask.MniY(mask.idx); 
    %these are used in selecting out the lQ vals
   hippoMasks(iHemi) = mask;%#ok<AGROW
end

% loop through subject
fPart = fullfile('Analysis','Alpha01','Mdl02');
for iColoc = -1:2:1
    imIdx = (iColoc >0) + 1;
    for iHemi=-1:2:1
        hemIdx = ((iHemi >0) + 1);
        colocation = ones(nSubs,1)*iColoc;
        hemisphere = ones(nSubs,1)*iHemi;
        lQ = cell(nSubs,1);
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
            lQ{iSubject} = clQ(hippoMasks(hemIdx).idx);
        end
        t = table(subjectIds,hemisphere,colocation,lQ);
        if iHemi==-1 && iColoc==-1
            DataTable02b = t;
        else
            DataTable02b = [DataTable02b;t];%#ok<AGROW>
        end
    end
end
DataTable02b.Properties.VariableNames{1} = 'subjectId';
%get performance on nonCom
dtP = get_pNonc(G);
dtP.mcPnonc = dtP.pNonc - mean( dtP.pNonc);
DataTable02b = join(DataTable02b,dtP);

% Cd back and save
cd(wd);
save(dtName,"DataTable02b","mmY");
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