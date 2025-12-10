function [Data,pCover] = getTpatterns_EpiRes(G,subjectId,roiId)
%Data is a [nVox,(nStim*nPos)] matrix of tStatistics;

if iscategorical(subjectId)
    subjectId = char(subjectId);
end

dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.EPI = [dirs.Subject,filesep,'EPI'];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

% Get ROI mask
roiPath = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiPath);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

%Spec the ROI coverage in ROI space
roiCoverage = zeros(size(roiMask.M));

% Get EPI mask
epiMask.name = sprintf('%s%s_%s_epiMask00.nii',dirs.EPI,filesep,subjectId);
epiMask.V = spm_vol(epiMask.name);

epiMask.M = spm_read_vols(epiMask.V);
epiMask.idx = find(epiMask.M > 0.5);
epiMask.size = size(epiMask.M);

% Spec the ROI mask in EPI space!
epiRoi.M = zeros(size(epiMask.M));

% Loop through the epiMask to:
%   1) populate epiRoi by samppling from roiMask
%   2) populate roiCoverage by keeping track of where we have sampled from
for iVoxel = 1:numel(epiMask.idx)

    % Step 1:
    [vx,vy,vz] = ind2sub(epiMask.size,epiMask.idx(iVoxel));
    vxEpi = [vx,vy,vz,1]';
    mm = epiMask.V.mat * vxEpi;
    vxRoi = roiMask.V.mat \ mm;
    sample = spm_sample_vol(roiMask.V, ...
        vxRoi(1),vxRoi(2),vxRoi(3),-5);
    epiRoi.M(vx,vy,vz) = sample;
    % Step 2:
    if sample > 0.5
        [w,idxToWeight] = getLerpW(vxRoi(1:3)',roiMask.size);
        roiCoverage(idxToWeight) = roiCoverage(idxToWeight) + w;
    end
end

% Get the masked EPI voxel indecies
epiRoi.idx = find(epiRoi.M>0.5);

% Get a matrix of t-stats
nStim = 6;
nPos = 2;
positions = 'ab';
Data = nan(numel(epiRoi.idx),nStim,nPos);
for iPos = 1:nPos
    cPos = positions(iPos);
    for iStim = 1:nStim
        stimId = [cPos,num2str(iStim-1)];
        tFn = [dirs.Alpha01,filesep,stimId,filesep,'spmT_0001.nii'];
        tV = spm_vol(tFn);
        tM = spm_read_vols(tV);
        Data(:,iStim,iPos) = tM(epiRoi.idx);
    end
end

%stack a stim ontop of b stim
% we want stim to be the columns as thats what corr expects
Data = [Data(:,:,1),Data(:,:,2)];

% Compute coverage stats
pCover = sum(roiCoverage>(1e-6),'all')/sum(roiMask.M,'all');
return








