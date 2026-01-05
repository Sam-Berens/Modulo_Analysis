function [Data,pCover] = getTpatterns_EpiRes(G,subjectId,roiId)

% Data is a [nVox,nPattern] matrix of tStatistics

%% Make sure subjectId is a char
if iscategorical(subjectId)
    subjectId = char(subjectId);
end

%% Get the ROI mask
roiMask.name = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiMask.name);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

%% Preallocate the ROI coverage in ROI space
roiCoverage = zeros(roiMask.size);

%% Get EPI mask
epiMask = getEpiMask(subjectId);

%% Preallocate the ROI mask in EPI space!
epiRoi.M = zeros(size(epiMask.M));
epiRoi.size = epiMask.size;

%% Loop through the epiMask to:
%  1) populate epiRoi by sampling from roiMask
%  2) populate roiCoverage by keeping track of where we have sampled from
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

%% Compute coverage stats
pCover = sum(roiCoverage>(1e-6),'all')/sum(roiMask.M,'all');

%% Get the masked EPI voxel indecies
epiRoi.idx = find(epiRoi.M>0.5);

%% Get the t-images
Timgs = getTimgs(subjectId); % Note we don't mask here

%% Format the data
n0 = prod(epiRoi.size);
nT = numel(Timgs.P);
Idx = n0*(0:(nT-1)) + epiRoi.idx;
Data = Timgs.M(Idx);
return