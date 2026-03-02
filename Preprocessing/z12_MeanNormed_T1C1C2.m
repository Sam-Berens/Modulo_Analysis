function [] = z12_MeanNormed_T1C1C2(G)

% Get the subjectIds
subjectIds = getSubjectIds(G);
nSubjects = numel(subjectIds);

% Set paths.Data
paths.Data = fullfile('..','..','Data');

% Set paths.GrgStruct
paths.GrpStruct = fullfile(paths.Data,'_Group',G,'Structural');

% Preallocate and populate filename arrays for images to norm
fns.T1s = cell(nSubjects,1);
fns.C1s = cell(nSubjects,1);
fns.C2s = cell(nSubjects,1);
for iSubject = 1:nSubjects
    paths.Subject = [paths.Data,filesep,char(subjectIds(iSubject))];
    paths.Structural = [paths.Subject,filesep,'Structural'];
    paths.SG1 = [paths.Structural,filesep,'G1'];

    dirList = dir([paths.SG1,filesep,'wm_*_T1*.nii']);
    fns.T1s{iSubject} = [dirList.folder,filesep,dirList.name];

    dirList = dir([paths.SG1,filesep,'wc1_*.nii']);
    fns.C1s{iSubject} = [dirList.folder,filesep,dirList.name];

    dirList = dir([paths.SG1,filesep,'wc2_*.nii']);
    fns.C2s{iSubject} = [dirList.folder,filesep,dirList.name];
end

% Mean the imags
meanT1s(G,paths.GrpStruct,fns.T1s);
meanCxs(G,paths.GrpStruct,fns.C1s,'C1');
meanCxs(G,paths.GrpStruct,fns.C2s,'C2');

return

function [] = meanT1s(G,grpStruct,fns)
V = spm_vol(char(string(fns)));
M = spm_read_vols(V);
Mu1 = mean(M,[1,2,3]);
Sd1 = std(M,1,[1,2,3]);
Z = (M - Mu1) ./ Sd1;
MuZ = mean(Z,4);
Vnew = V(1);
Vnew.fname = fullfile(grpStruct,'MeanT1',[G,'_MeanT1.nii']);
Vnew.dt(1) = 64;
Vnew.descrip = [G,' mean T1 image;'];
if ~exist(fileparts(Vnew.fname),'dir')
    mkdir(fileparts(Vnew.fname));
end
spm_write_vol(Vnew,MuZ);
return

function [] = meanCxs(G,grpStruct,fns,C)
V = spm_vol(char(string(fns)));
M = spm_read_vols(V);
MuM = mean(M,4);
Vnew = V(1);
Vnew.fname = fullfile(grpStruct,['Mean',C],[G,'_MeanT1.nii']);
Vnew.dt(1) = 64;
Vnew.descrip = [G,' mean ',C,' image;'];
if ~exist(fileparts(Vnew.fname),'dir')
    mkdir(fileparts(Vnew.fname));
end
spm_write_vol(Vnew,MuM);
return