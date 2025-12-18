function [roiTData,pCover] = getRoiTs_EpiRes(G,subjectId,roiId)
%Data is a [nVox,(nStim*nPos)] matrix of tStatistics;

if iscategorical(subjectId)
    subjectId = char(subjectId);
end

% Get ROI mask
roiPath = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiPath);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

%Spec the ROI coverage in ROI space
roiCoverage = zeros(size(roiMask.M));

% Get EPI mask
[epiMask] = getEpiMask(subjectId);

% Spec the ROI mask in EPI space!
epiRoi.M = zeros(size(epiMask.M));

% Loop through the epiMask to:
%   1) populate epiRoi by sampling from roiMask
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
gnStim = nPos*nStim;
roiTData = nan(numel(epiRoi.idx),gnStim);
[tData] = getAlpha01Ts(subjectId,epiMask);
for iStim = 1:gnStim
    cTData = tData(:,:,:,iStim);%this is 3d
    roiTData(:,iStim) = cTData(epiRoi.idx);
end

% Compute coverage stats
pCover = sum(roiCoverage>(1e-6),'all')/sum(roiMask.M,'all');
return








