function [mskd_Yim, mskd_eIm] = getMskdIms_EpiRes(G,subjectId,roiId,yIm,eIm)

%% Make sure subjectId is a char
if iscategorical(subjectId)
    subjectId = char(subjectId);
end

%% Get the ROI mask
roiMask.name = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiMask.name);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

%% Get EPI mask
epiMask = getEpiMask(subjectId);

%% Preallocate the ROI mask in EPI space!
epiRoi.M = zeros(size(epiMask.M));
epiYim.M = zeros(size(epiMask.M));
epiRoi.size = epiMask.size;

%% Loop through the epiMask to:
%  1) populate epiRoi by sampling from roiMask
for iVoxel = 1:numel(epiMask.idx)
    % Step 1:
    [vx,vy,vz] = ind2sub(epiMask.size,epiMask.idx(iVoxel));
    vxEpi = [vx,vy,vz,1]';
    mm = epiMask.V.mat * vxEpi;
    vxRoi = roiMask.V.mat \ mm;
    sample1 = spm_sample_vol(roiMask.V, ...
        vxRoi(1),vxRoi(2),vxRoi(3),-5);
    sample2 = spm_sample_vol(yIm.V, ...
        vxRoi(1),vxRoi(2),vxRoi(3),-5);
    epiRoi.M(vx,vy,vz) = sample1;
    epiYim.M(vx,vy,vz) = sample2;
end

%% Get the masked EPI voxel indecies
epiRoi.idx = find(epiRoi.M>0.5);

eIm.size = size(eIm.M);
%% Get the stat-image
if ~isequal(epiRoi.size,eIm.size)
    error('The EPI mask and the stat-image are not the same size!');
end
if isfield(epiMask,'V') && isfield(epiMask.V,'mat')
    if ~isequal(epiMask.V.mat, eIm.V.mat)
        error('The EPI mask and the stat-image are not aligned!');
    end
end

%% Format the data
%mask both images with nans outside the roi
mskd_eIm = eIm.M(epiRoi.idx);
mskd_Yim = epiYim.M(epiRoi.idx);

return