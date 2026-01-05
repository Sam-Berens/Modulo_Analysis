function [Timgs] = getTimgs(subjectId, Mask)
% GETTIMGS Load and optionally mask first-level SPM T-images for a subject.
%
%   Timgs = GETTIMGS(subjectId)
%   Timgs = GETTIMGS(subjectId, Mask)
%
%   Loads a fixed set of SPM T-statistic images for a given subject from
%   the Alpha01 first-level analysis directory and returns them stacked
%   into a single 4D array.
%
%   The function assumes a factorial structure with two places (a, b) and
%   six intervals (0–5), yielding 12 T-images in total. Images are loaded
%   from (e.g.):
%    ../../Data/subjectId/Analysis/Alpha01/a0/spmT_0001.nii
%
%   INPUTS
%   ------
%   subjectId : char | string
%       Subject identifier corresponding to a directory name in the Data
%       folder (e.g., 'eade18a5').
%
%   Mask : struct | numeric array (optional)
%       Brain mask used to exclude voxels outside the region of interest.
%       Can be provided in one of two forms:
%
%       1) Numeric 3D array:
%          - Voxels > 0.5 are treated as in-mask.
%
%       2) Struct with fields:
%          - M    : 3D logical or numeric mask volume
%          - idx  : Linear indices of in-mask voxels
%          - size : Size of the mask volume
%          - V    : (optional) SPM volume struct for alignment checking
%
%       If provided, voxels outside the mask are set to NaN in the output.
%
%   OUTPUT
%   ------
%   Timgs : struct
%       Structure containing the loaded T-images and metadata:
%
%       Timgs.M
%           4D array of T-statistic images
%           [X × Y × Z × N], where N = 12.
%
%       Timgs.V
%           1 × N struct array of SPM volume headers corresponding to the
%           T-images.
%
%       Timgs.P
%           N × 1 cell array of full file paths to each T-image.
%
%       Timgs.idx
%           Linear indices of voxels retained by the mask (only present if
%           a mask is supplied).
%
%       Timgs.size
%           Size of the 3D image volume [X Y Z].
%
%   NOTES
%   -----
%   - The function assumes that all T-images are in the same space and
%     share identical dimensions and affine transforms. 
%   - If a mask struct includes an SPM volume field (Mask.V.mat), alignment
%     between the mask and T-images is explicitly checked.
%   - Voxels outside the mask are set to NaN rather than zero to prevent
%     contamination of downstream statistical analyses.
%
%   DEPENDENCIES
%   ------------
%   Requires SPM on the MATLAB path:
%       - spm_vol
%       - spm_read_vols
%
%   EXAMPLE
%   -------
%   % Load unmasked T-images
%   Timgs = getTimgs('eade18a5');
%
%   % Load masked T-images
%   Mask = getEpiMask('eade18a5');
%   Timgs = getTimgs('eade18a5', Mask);
%
%   See also SPM_VOL, SPM_READ_VOLS.

%% Set some vars
nInts = 6;
nPlcs = 2;
nTimg = nInts * nPlcs;

%% Check inputs
maskGiven = (nargin == 2);
if nargin > 2
    error('Too many inputs!');
elseif nargin < 1
    error('Too few inputs!');
end

%% Normalise Mask if it is provided as a matrix (not a struct)
if maskGiven && ~isstruct(Mask)
    Mraw = Mask;
    Mask = struct();
    Mask.M = Mraw > 0.5;
    Mask.idx = find(Mask.M);
    Mask.size = size(Mraw);
elseif maskGiven && isstruct(Mask) && ~isfield(Mask,'idx')
    % "Mask-like" struct but missing fields
    Mraw = Mask;
    Mask = struct();
    Mask.M = Mraw > 0.5;
    Mask.idx = find(Mask.M);
    Mask.size = size(Mraw);
end

%% Preallocate the output
Timgs = struct();
if maskGiven
    Timgs.idx  = Mask.idx;
    Timgs.size = Mask.size(1:3);
end

%% Set directories
dirs.Data    = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

%% Construct filenames
Timgs.P = cell(nTimg,1);
for iT = 1:nTimg
    plc = char(floor((iT-1)/nInts) + double('a'));
    int = num2str(mod((iT-1),nInts));
    Timgs.P{iT} = [dirs.Alpha01,filesep,plc,int,filesep,'spmT_0001.nii'];
end

%% Load headers and data
Timgs.V = spm_vol(char(Timgs.P));

% Read volumes into 4D array
vol1 = spm_read_vols(Timgs.V(1));
M = nan([size(vol1), nTimg], 'like', vol1);
M(:,:,:,1) = vol1;
for iT = 2:nTimg
    M(:,:,:,iT) = spm_read_vols(Timgs.V(iT));
end
Timgs.M = M;

%% Check mask size/alignment (only if info exists)
if maskGiven
    if ~isequal(Timgs.V(1).dim, Mask.size)
        error('The mask and the t-images are not the same size!');
    end
    if isfield(Mask,'V') && isfield(Mask.V,'mat')
        if ~isequal(Mask.V.mat, Timgs.V(1).mat)
            error('The mask and the t-images are not aligned!');
        end
    end
end

%% Apply mask
if maskGiven
    Mmask = repmat(logical(Mask.M), [1,1,1,nTimg]);
    Timgs.M(~Mmask) = NaN;
else
    sz = size(Timgs.M);
    Timgs.size = sz(1:3);
end
return