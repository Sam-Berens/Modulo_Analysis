function [Y,N] = searchlight3D(r,f,Mask,X,yDepth,searchID)

%% Check inputs
if nargin < 5
    error('Not enough inputs.');
elseif nargin == 5
    waitmssg = 'Running searchlight ...';
elseif nargin == 6
    waitmssg = ['Running searchlight: ',searchID];
else
    error('Too many inputs');
end

% Check whether the Mask has been loaded by getEpiMask.m ...
% ... and adapt if not.
if ~isfield(Mask,'idx') || ~isfield(Mask,'size')
    Mraw = Mask;
    Mask = struct();
    Mask.M = Mraw > 0.5;
    Mask.idx = find(Mask.M);
    Mask.size = size(Mraw,[1,2,3]);
end
Mask.nIdx = numel(Mask.idx);

% Check whether X has been loaded by (e.g.) getTimgs.m ...
% ... and adapt it if not.
if ~isfield(X,'size') || ~isfield(X,'M')
    Xraw = X;
    X = struct();
    X.M = Xraw;
    X.size = size(Xraw,[1,2,3]);
end
X.nVox = prod(X.size);
X.nPttrns = size(X.M,4);

% Check that the size of epiMask is consistent with X
if any(Mask.size ~= X.size)
    error('Size mismatch between mask and X');
end

% Check mask size/alignment (only if info exists)
if isfield(Mask,'V') && isfield(X,'V')
    if ~isequal(Mask.V.mat, X.V(1).mat)
        error('The mask and the data in X are not aligned!');
    end
end

%% Make the searchlight (Ball)
Ball = makeBall(r);

%% Loop through the Mask to apply the searchlight
Y = nan([Mask.size,yDepth]);
N = nan(Mask.size);
fh = waitbar(0,waitmssg);
for ii = 1:Mask.nIdx

    % Get the central voxel
    [x,y,z] = ind2sub(Mask.size,Mask.idx(ii));
    xyz = [x,y,z];

    % Construct the neighbourhood
    Hood = Ball + xyz;

    % Contain the neighbourhood within the FOV
    Hood = constrain(Hood,Mask.size);

    % Convert to linear indices
    Hood = mat2cell(Hood,size(Hood,1),[1,1,1]);
    Hood = sub2ind(Mask.size,Hood{:});

    % Ensure all voxels are in the mask
    selector = ismember(Hood,Mask.idx);
    Hood = Hood(selector);

    % Ensure the neighbourhood has at least 3 voxels
    n = numel(Hood);
    if n < 3
        Y(Mask.idx(ii)) = NaN;
        N(Mask.idx(ii)) = 0;
        continue
    end

    % Index all patterns (in M)
    IdxRead = X.nVox*(0:(X.nPttrns-1)) + Hood;
    M = X.M(IdxRead);

    % Apply the function
    IdxWrit = X.nVox*(0:(yDepth-1)) + Mask.idx(ii);
    Y(IdxWrit) = f(M);
    N(Mask.idx(ii)) = n;

    % Update the waitbar
    if mod(ii,97) == 96
        waitbar(ii/Mask.nIdx,fh);
    end
end
close(fh);

return

function [Ball] = makeBall(r)
d = ceil(r)*2 + 1;
offset = ceil(r) + 1;
n = d^3;
Ball = nan(n,3);
for ii = 1:n
    [x,y,z] = ind2sub([d,d,d],ii);
    xyz = [x,y,z] - offset;
    if norm(xyz) <= r
        Ball(ii,:) = xyz;
    end
end
Ball = Ball(~any(isnan(Ball),2),:);
return

function [Xyz] = constrain(Xyz,sz)
c1 = Xyz > 0;
c2 = Xyz <= sz;
selector = all(c1&c2,2);
Xyz = Xyz(selector,:);
return