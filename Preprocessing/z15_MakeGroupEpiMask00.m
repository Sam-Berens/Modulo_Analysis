function [] = z15_MakeGroupEpiMask00(G)

% Get the canonical mask
fns.spmMask = fullfile(getSpmHome(),'tpm','mask_ICV.nii');
spmMask.V = spm_vol(fns.spmMask);

% Get the subjectIds
subjectIds = getSubjectIds(G);
nSubjects = numel(subjectIds);

% Set paths.Data
paths.Data = fullfile('..','..','Data');

% Set paths.GrgStruct
paths.GrpStruct = fullfile(paths.Data,'_Group',G,'Structural');

% Preallocate and populate filename arrays for images to norm
fns.epiMask = cell(nSubjects,1);
for iSubject = 1:nSubjects
    paths.Subject = [paths.Data,filesep,char(subjectIds(iSubject))];
    paths.EPI = [paths.Subject,filesep,'EPI'];
    paths.SG1 = [paths.EPI,filesep,'G1'];

    dirList = dir([paths.SG1,filesep,'w_*_epiMask00.nii']);
    fns.epiMask{iSubject} = [dirList.folder,filesep,dirList.name];

end

% Read in the subject specific normed EPI masks
grpMask.V = spm_vol(char(string(fns.epiMask)));
grpMask.M = spm_read_vols(grpMask.V);

% Produce a group mask
grpMask.M = sum(grpMask.M,4,'omitmissing');
grpMask.M = grpMask.M >= 16;
grpMask.idx = find(grpMask.M);

% Sample from the canonical image
[x,y,z] = ind2sub(grpMask.V(1).dim,grpMask.idx);
Xyz = [x,y,z,ones(size(x))]';
mm = grpMask.V(1).mat * Xyz;
Xyz = spmMask.V.mat \ mm;
Xyz = Xyz(1:3,:)';
Xyz = mat2cell(Xyz,size(Xyz,1),ones(1,3));
s = spm_sample_vol(spmMask.V,Xyz{:},-7);

% Clip the group mask that is outside the canonical image
grpMask.idx = grpMask.idx(s > 0.5);
grpMask.M = zeros(grpMask.V(1).dim);
grpMask.M(grpMask.idx) = 1;

Vnew = grpMask.V(1);
Vnew.fname = fullfile(paths.GrpStruct,'GrpEpiMask00',...
    [G,'_GrpEpiMask00.nii']);
Vnew.dt(1) = 2;
Vnew.descrip = [G,' group EPI mask 00;'];
if ~exist(fileparts(Vnew.fname),'dir')
    mkdir(fileparts(Vnew.fname));
end
spm_write_vol(Vnew,grpMask.M);
return

function [spmRoot] = getSpmHome()
matlabPath = path;
% Get all paths in the MATLAB path
allPaths = strsplit(matlabPath, pathsep)'; 
spmRoot = allPaths{cellfun(@(x) strcmp(x(end-4:end), 'SPM25'), allPaths)};
return