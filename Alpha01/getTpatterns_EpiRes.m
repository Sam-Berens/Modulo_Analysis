function [Data,pCover] = getTpatterns_EpiRes(G,subjectId,roiId)
% GETTPATTERNS_EPIRES  Extract Alpha01 T-image patterns at EPI resolution.
%
%   [Data, pCover] = GETTPATTERNS_EPIRES(G, subjectId, roiId)
%
%   Loads a fixed set of first-level SPM T-statistic images for a subject
%   (typically 12 images: two places × six integers) and extracts
%   voxelwise T-values within an ROI represented in EPI space ("EPI
%   resolution").
%
%   This function is intended for pipelines where an ROI defined in native
%   space is mapped/resampled into the subject's EPI grid (or where an EPI
%   ROI already exists). The output is a matrix with one row per in-ROI EPI
%   voxel and one column per T-image/condition.
%
%   INPUTS
%   ------
%   G : char | string
%       Group identifier used by the data pipeline (e.g., 'G0', 'G1').
%
%   subjectId : char | string
%       Subject identifier corresponding to a subject folder in ../../Data.
%
%   roiId : char | string
%       ID of the ROI. The ROI mask is expected to be located via the
%       pipeline helper (e.g., getNativeRoiPath) and mapped into EPI space.
%
%   OUTPUTS
%   -------
%   Data : double [nVox × nT]
%       Matrix of T-values extracted within the EPI-space ROI. Rows
%       correspond to in-ROI voxels, columns correspond to T-images/
%       conditions.
%
%   pCover : double
%       ROI coverage metric: the proportion of ROI voxels for which
%       coverage is non-zero (or above a small threshold) when intersected
%       with the EPI mask.
%
%   NOTES
%   -----
%   - This function assumes that the ROI and T-images are defined on
%     compatible voxel grids. TGhe first three dimensions of the T-image
%     volumes must match the EPI grid size.
%   - Requires SPM on the MATLAB path.
%
%   EXAMPLE
%   -------
%   [Data, pCover] = getTpatterns_EpiRes('G1','eade18a5','lHippAnt');
%
%   See also GETTIMGS, GETEPIMASK, SPM_VOL, SPM_READ_VOLS.

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
pCover = sum(roiCoverage>(1e-6),'all')/sum(roiMask.M>0.5,'all');

%% Get the masked EPI voxel indecies
epiRoi.idx = find(epiRoi.M>0.5);

%% Get the t-images
Timgs = getTimgs(subjectId); % Note we don't mask here
if ~isequal(epiRoi.size,Timgs.size)
    error('The EPI mask and the t-images are not the same size!');
end
if isfield(epiMask,'V') && isfield(epiMask.V,'mat')
    if ~isequal(epiMask.V.mat, Timgs.V(1).mat)
        error('The EPI mask and the t-images are not aligned!');
    end
end

%% Format the data
n0 = prod(epiRoi.size);
nT = numel(Timgs.P);
Idx = n0*(0:(nT-1)) + epiRoi.idx;
Data = Timgs.M(Idx);
return