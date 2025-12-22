function [srchIm] = searchLight(imageM,mask,fh,radius)
%image is a 4D volume [x,y,z,nConditions]
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
srchIm = nan([size(imageM,1,2,3),nSlices]);

%flatten imageM so that we can easily index the search light area for all
%conditions at once
imageM = reshape(imageM, [], size(imageM,4));
for iVox=1:nCentres
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

    %% Extract searchlight-voxels from the input image,across all conditions
    % cVals = arrayfun(@(x,y,z) imageM(x,y,z,:),v(:,1), v(:,2), v(:,3));
    linIdx = sub2ind([maxX maxY maxZ], v(:,1), v(:,2), v(:,3));
    cVals = imageM(linIdx,:);
    %apply your function to the searchlight
    fOfVals = fh(cVals);
    
    if iVox==1
        nSlices = size(fOfVals,1);
        %preallocate output based on how many result images the fh produces
        srchIm = nan([size(imageM,1,2,3),nSlices]);
    end 
   
    %this skips searchlights whos stats don't meet the p inequality
    if isnan(fOfVals(1,1))
        continue 
    end

    %fill up all condition images at that location with result of function
    for iSlice=1:numel(nSlices)
        srchIm(x,y,z,iSlice) = fOfVals(iSlice,1);
    end
end

return
