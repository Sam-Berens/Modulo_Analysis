function [w,idxToWeight] = getLerpW(vxRoi,roiSize)
x = vxRoi(1);
y = vxRoi(2);
z = vxRoi(3);

x0 = floor(x);
x1 = x0 + 1;
y0 = floor(y);
y1 = y0 + 1;
z0 = floor(z);
z1 = z0 + 1;

dx = x - x0;
dy = y - y0;
dz = z - z0;

w000 = (1-dx)*(1-dy)*(1-dz);
w100 = dx    *(1-dy)*(1-dz);
w010 = (1-dx)*dy    *(1-dz);
w110 = dx    *dy    *(1-dz);
w001 = (1-dx)*(1-dy)*dz;
w101 = dx    *(1-dy)*dz;
w011 = (1-dx)*dy    *dz;
w111 = dx    *dy    *dz;

idxVox = [
    x0 y0 z0;
    x1 y0 z0;
    x0 y1 z0;
    x1 y1 z0;
    x0 y0 z1;
    x1 y0 z1;
    x0 y1 z1;
    x1 y1 z1];

w = [w000; w100; w010; w110; w001; w101; w011; w111];

% Bounds check (edge cases)
inside = idxVox(:,1) >= 1 & idxVox(:,1) <= roiSize(1) & ...
    idxVox(:,2) >= 1 & idxVox(:,2) <= roiSize(2) & ...
    idxVox(:,3) >= 1 & idxVox(:,3) <= roiSize(3);

idxVox = idxVox(inside,:);
w = w(inside);

idxToWeight = sub2ind(roiSize, idxVox(:,1), idxVox(:,2), idxVox(:,3));
return