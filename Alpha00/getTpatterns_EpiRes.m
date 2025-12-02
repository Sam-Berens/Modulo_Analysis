function [Data,pCover] = getTpatterns_EpiRes(G,subjectId,roiId)
if iscategorical(subjectId)
    subjectId = char(subjectId);
end

dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.EPI = [dirs.Subject,filesep,'EPI'];
dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];


% Get ROI mask
roiPath = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiPath);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

% Spec the ROI coverage in ROI space
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
Data = nan(numel(epiRoi.idx),6);
for iStim = 1:6
    stimId = ['i',num2str(iStim-1)];
    tFn = [dirs.Alpha00,filesep,stimId,filesep,'spmT_0001.nii'];
    tV = spm_vol(tFn);
    tM = spm_read_vols(tV);
    Data(:,iStim) = tM(epiRoi.idx);
end

% Compute coverage stats
pCover = sum(roiCoverage>(1e-6),'all')/sum(roiMask.M,'all');
return