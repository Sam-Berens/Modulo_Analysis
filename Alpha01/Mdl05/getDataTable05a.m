function [DataTables05a] = getDataTable05a(G,roiId)
%GETDATATABLE05A Build ROI-wise data tables for Mdl05 template statistics.
%
%   DataTables05a = GETDATATABLE05A(G, roiId) loads subject-level
%   non-commutable scores and ROI-averaged statistical images for the
%   group/model identified by G, returning one analysis table per ROI.
%
%   INPUTS
%   ------
%   G
%       Group identifier, specified as a character vector or string scalar.
%       This must begin with 'G' and is used to locate subject-level
%       statistical images in:
%
%           ../../../Data/<subjectId>/Analysis/Alpha01/Mdl05/<G>/
%
%   roiId
%       ROI identifier, or a cell array of ROI identifiers. Each ROI
%       identifier is passed to getMNIRoiPath, loaded as an ROI mask, and
%       used as a field name in the returned structure.
%
%   OUTPUT
%   ------
%   DataTables05a
%       Structure containing one table per ROI. The field names are taken
%       directly from roiId. For example:
%
%           DataTables05a.<roiId>
%
%       Each table contains two rows per subject: one row for the
%       colocation = -1 contrast and one row for the colocation = +1
%       contrast. The table is formed by outer-joining the pNonc table
%       returned by get_pNonc(G) with the ROI-averaged image statistics
%       described below.
%
%   VARIABLES IN EACH ROI TABLE
%   ---------------------------
%   subjectId
%       Subject identifier. This is used as the key when joining the
%       pNonc data with the ROI-averaged image statistics.
%
%   pNonc
%       Subject-level non-commutable score, as returned by get_pNonc(G).
%
%   cpNonc
%       Mean-centred non-commutable score:
%
%           cpNonc = pNonc - mean(pNonc)
%
%   colocation
%       Condition/contrast code for the colocation effect.
%       This takes values:
%
%           -1  corresponding to wzTemplate_colocation=-1.nii
%           +1  corresponding to wzTemplate_colocation=+1.nii
%
%   zTemplate
%       Mean statistic value within the ROI for the relevant colocation
%       template image. There is one value per subject and colocation
%       condition.
%
%   zVisSim
%       Mean statistic value within the ROI for:
%
%           wzTemplate_visSim.nii
%
%       This value is computed once per subject and repeated across that
%       subject's two colocation rows.
%
%   NOTES
%   -----
%   The function first checks whether all ROI masks and statistical images
%   are defined on the same voxel grid. If they are, ROI means are computed
%   directly using linear indexing into the loaded image volumes. If not,
%   each statistical image is sampled at the world-coordinate locations of
%   the ROI voxels using spm_sample_vol, and finite sampled values are
%   averaged.
%
%   This function requires SPM and the helper functions get_pNonc,
%   getMNIRoiPath, and the nested ROI-loading/sampling utilities defined
%   below.
%
%   EXAMPLE
%   -------
%   DataTables05a = getDataTable05a('G1', {'HPC', 'EC'});
%
%   T_HPC = DataTables05a.HPC;
%   T_EC  = DataTables05a.EC;

DataTables05a = struct;
if nargin < 1
    error('At least one argument (GroupId) required.');

elseif ~startsWith(G,'G')
    error('First argument must be GroupId.');

elseif nargin == 1

elseif (nargin > 1) && ~iscell(roiId)
    roiId = {roiId};

elseif nargin > 3
    error('Maximum of 3 arguments: GroupId and (optional) roiPattern');
end


%% Add .. to path
addpath('../');

%% Load the ROI masks and determine if they are on the same grid
roi = cell(size(roiId));
for iRoi = 1:numel(roiId)
    roi{iRoi} = loadMask(getMNIRoiPath(roiId{iRoi}));
end
affine = cellfun(@(s)s.V.mat,roi,'UniformOutput',false);
maskDim = cellfun(@(s)s.V.dim(1:3),roi,'UniformOutput',false);
allMasksTheSame = isscalar(affine) || ...
    (isequal(affine{:}) && isequal(maskDim{:}));

%% Get pNonc
Pnonc = get_pNonc(G);
Pnonc.cpNonc = Pnonc.pNonc - mean(Pnonc.pNonc);
subjectIds = Pnonc.subjectId;

%% Load the statistical imges
a = cellstr(repelem(Pnonc.subjectId,3,1));
b = repmat({ ...
    'wzTemplate_colocation=-1.nii'; ...
    'wzTemplate_colocation=+1.nii'; ...
    'wzTemplate_visSim.nii'}, ...
    size(subjectIds,1),1);

zTemplateFn = cellfun(@(x,y) fullfile(...
    '..','..','..','Data',x,'Analysis',...
    'Alpha01','Mdl05',G,y),...
    a,b,'UniformOutput',false);

V = spm_vol(char(zTemplateFn));

%% Test if all images and masks are on the same grid
imgAffine = arrayfun(@(s)s.mat,V,'UniformOutput',false);
imgDim = arrayfun(@(s)s.dim(1:3),V,'UniformOutput',false);
allStatImgsTheSame = isscalar(V) || ...
    (isequal(imgAffine{:}) && isequal(imgDim{:}));
allImgsTheSame = allMasksTheSame && allStatImgsTheSame && ...
    isequal(affine{1},V(1).mat) && isequal(maskDim{1},V(1).dim(1:3));

%% Get/Sample the stats
zTemplate_Cell = cell(size(roiId));
if allImgsTheSame
    M = spm_read_vols(V);
    nVox = prod(V(1).dim);
end
for iRoi = 1:numel(roiId)
    if allImgsTheSame
        idx = (nVox * (0:(size(M,4)-1))') + roi{iRoi}.idx';
        zTemplate_Cell{iRoi} = mean(M(idx),2);
    else
        zTemplate_Cell{iRoi} = sampleV(roi{iRoi},V);
    end
end

%% Make the data tables!
for iRoi = 1:numel(roiId)
    zTemplate = zTemplate_Cell{iRoi};

    colocation = repmat([-1;1],numel(subjectIds),1);
    s = false(size(zTemplate));
    s(3:3:end) = true;
    zVisSim = repelem(zTemplate(s),2,1);
    zTemplate = zTemplate(~s);

    subjectId = repelem(subjectIds,2,1);
    T = table(subjectId,colocation,zTemplate,zVisSim);
    DataTables05a.(roiId{iRoi}) = outerjoin(Pnonc,T,"MergeKeys",true);
end

%% Remove .. from path
rmpath('../');
return

function [mask] = loadMask(filename)
V = spm_vol(filename);
mask.V = V;
mask.M = spm_read_vols(mask.V);
mask.idx = find(mask.M > 0.5);
return

function [zTemplate] = sampleV(roi,V)
% The ROI mask and the statistical images are not in the same voxel
% grid.  Sample each statistical image at the centre of every ROI
% voxel after transforming those voxel centres into the statistical
% image's voxel coordinate system.
zTemplate = nan(numel(V),1);

[x,y,z] = ind2sub(roi.V.dim(1:3),roi.idx(:)');
roiVoxXYZ = [x;y;z;ones(1,numel(x))];
roiWorldXYZ = roi.V.mat * roiVoxXYZ;

for iImg = 1:numel(V)
    imgVoxXYZ = V(iImg).mat \ roiWorldXYZ;
    sampledVals = spm_sample_vol(V(iImg), ...
        imgVoxXYZ(1,:),imgVoxXYZ(2,:),imgVoxXYZ(3,:),1);

    % spm_sample_vol can return NaNs where the ROI extends outside
    % the statistical image's field of view.  Exclude only those
    % invalid samples from the ROI average.
    sampledVals = sampledVals(isfinite(sampledVals));
    if ~isempty(sampledVals)
        zTemplate(iImg) = mean(sampledVals);
    end
end
return