function [srchIm] = searchLight(image,mask,fh,radius)
%image is a 4D volume with stimulus condition in the 4th dim
%mask is the mask for which you want to sample the image
%fh is the figure handle for the function you want to apply to the
%searchlight
%radius is the size of the searchlight you want to look with
nCentres = numel(mask.idx);
%this gives us the bounds of the actual volume (should be 104x104x66)
[maxX,maxY,maxZ] = size(mask.M);

%% Make ball:
[dx,dy,dz] = meshgrid(-radius:radius,-radius:radius,-radius:radius);
dxdydz = cat(4,dx,dy,dz);
dxdydz = mat2cell(dxdydz,ones(1,7),ones(1,7),ones(1,7),radius);
ball = cellfun(@(c)norm(squeeze(c))<=radius,dxdydz);
[sx,sy,sz] = ind2sub(size(ball),find(ball));
S = [sx,sy,sz];
%make cords relative to centre
S = S - ((size(ball)-1)./2) - 1;

nSlices = nargout(fh);
srchIm = nan([size(image,1,2,3),nSlices]);

for iSlice=1:numel(nSlices)
result.(['x',(int2str(iSlice))]) = nan(size(image,1,2,3));
end

for iVox=1:nCentres %TO DO refactor to searchlight function which takes arg of im and function to apply, and a radius
    %TO DO change testOmega() so it just loads the q images once for both
    %hemisphere
    %get sub coords for your central voxel
    [x,y,z] = ind2sub(mask.size,mask.idx(iVox));
    xyz = [x,y,z];
    %get the ball subscripts when centrered over current vox
    v = xyz + S;
    %remove voxel which are outside the field of view
    inFOV = ...
        v(:,1) >= 1 & v(:,1) <= maxX & ...
        v(:,2) >= 1 & v(:,2) <= maxY & ...
        v(:,3) >= 1 & v(:,3) <= maxZ;
    v = v(inFOV,:);

    %% Extract tStatistics for relevant voxels (remember image should be [x,y,z,nConditions]
    cBall = image(v(:,1),v(:,2),v(:,3),:);
    %apply your function to the searchlight
    [varargout{1:nSlices}] = fh(cBall);
    if isnan(varargout{:})
        continue
    end
    %fill up all condition images at that location with result of function
    for ir=1:numel(nSlices)
        result.(iSlice)(x,y,z) = varargout{ir};
    end
end

for ir=1:numel(nSlices)
    srchIm(:,:,:,iSlice) = result.(['x',(int2str(iSlice))]);
end

return