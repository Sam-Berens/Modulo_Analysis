function [] = makeQ(subjectId,roiId)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl02 = [dirs.Alpha01,filesep,'Mdl02'];
dirs.G1 = [dirs.Mdl02,filesep,'G1'];
if ~exist(dirs.Mdl02 ,'var')
    mkdir(dirs.Mdl02);
end

epiMask = getEpiMask(subjectId);
nCentres = numel(epiMask.idx);
%this gives us the bounds of the actual volume (should be 104x104x66)
[maxX,maxY,maxZ] = size(epiMask.M);

%% Get a matrix of t-stats for all conditions
%remember that a stims are stacked ontop of b stim
[tM] = getAlpha01Ts(subjectId,epiMask);

%% Make ball:
[dx,dy,dz] = meshgrid(-3:3,-3:3,-3:3);
dxdydz = cat(4,dx,dy,dz);
dxdydz = mat2cell(dxdydz,ones(1,7),ones(1,7),ones(1,7),3);
ball = cellfun(@(c)norm(squeeze(c))<=3,dxdydz);
[sx,sy,sz] = ind2sub(size(ball),find(ball));
S = [sx,sy,sz];
%make cords relative to centre
S = S - ((size(ball)-1)./2) - 1;

%% loop through searchlight centres to test produce q term from precursor to mdl02
%reminder: these vectors are not the full size of the image matrix, they 
%just correspond to the idxs which are 1 in the epiMask
Q = nan(nCentres,1);
B0 = nan(nCentres,1);
B1 = nan(nCentres,1);
rH1 = nan(nCentres,1);
rH2 = nan(nCentres,1);
rH3 = nan(nCentres,1);
err = nan(nCentres,1);

for iVox=1:nCentres
    %get sub coords for your central voxel
    [x,y,z] = ind2sub(epiMask.size,epiMask.idx(iVox));
    xyz = [x,y,z];
    %get the ball subscripts when centrered over current vox
    v = xyz + S;
    %remove voxel which are outside the field of view
    inFOV = ...
        v(:,1) >= 1 & v(:,1) <= maxX & ...
        v(:,2) >= 1 & v(:,2) <= maxY & ...
        v(:,3) >= 1 & v(:,3) <= maxZ;
    v = v(inFOV,:);

    %% Extract tStatistics for relevant voxels (remember t stat is[x,y,z,(iStim*pos)]
    tStats = tM(v(:,1),v(:,2),v(:,3),:);
    %squeeze into 2d so the column is the stim condition
    tStats = reshape(tStats, [], size(tStats,4));
  
    %estimate model which predicts the similarity score of a given dist bin
    [pHat,cRhoHat,cErr] = estimQ(tStats);
    %dont add Q if inequality did not hold 
    %(output will be set to false in this case)
    if isnan(pHat)
        continue
    end  
 
    B0(iVox,1) = pHat(1,1);
    B1(iVox,1) = pHat(2,1);
    Q(iVox,1) = pHat(3,1);
    %these are the similarity scores for each distance bin
    rH1(iVox,1) = cRhoHat(1,1); 
    rH2(iVox,1) = cRhoHat(2,1);
    rH3(iVox,1) = cRhoHat(3,1);
    err(iVox,1) = cErr;
end
%we're going to save all stats maps as nifti so we can norm q images
%to mni space and so that we can do QA checks on the others  
images = {Q,B0,B1,rH1,rH2,rH3,err};
names = {'q','b0','b1','rH1','rH2','rH3','error'};
for iIm=1:numel(images)
    im.M = nan(size(epiMask.M,[1,2,3]));
    im.M(epiMask.idx) = images{iIm};
    im.V = epiMask.V;
    label = names{iIm};
    %save to their mdl02 folder
    im.V.fname = [dirs.Mdl02,filesep,label,'_',roiId,'.nii'];
    im.V.descrip = 'Statistic map of nonlinear pattern–RDM fit model';
    spm_write_vol(im.V,im.M);
end
return

function [varargout] = estimQ(tPatterns)
D = corr(tPatterns);
lower = tril(true(6),-1);
D = D(lower); %have checked that these are ordered the same as nchoose 
%we will compute the arithmetic mean of all similarity scores r_(x,y)
% that have the same predicted value according to our structural hypothesis
nrows =6;
f = @(x,y)min(mod((x-y),6),mod((y-x),6));
distPairs = nchoosek(1:nrows,2);
x = distPairs(:,1);
y = distPairs(:,2);
distBin = f(x,y); 
%pi=b0+b1(1-i3)q
%pi is the mean simiarity scores(rho) for
%p1, p2, and p3(where i is the modular distance)
P = nan(3,1);
I = nan(3,1);
for ii=1:3
    idxs = find(distBin==ii);
    %PUT p into the design Matrix
    P(ii,1) = mean(D(idxs));
    %I is the modular distance between pairs of sparks
    %transform I to everything inside the brackets of the model
    I(ii,1) = ii;
end
%check whether p1>p2>p3 inequality holds
validP = P(1) > P(2) & P(2) > P(3);
if validP
    [pHat,rhoHat,err] = fitNonlinMdl(I,P);
    [varargout{1,1}, varargout{2,1}, varargout{3,1}]...
        = deal(pHat,rhoHat,err); 
else 
     varargout(1:3,1) = {nan};
end
return

function [pHat,rhoHat,err] = fitNonlinMdl(dist,rho)
similarity = 1-(dist./3);
rhoHatFnc = @(p) p(1) + p(2).*(similarity.^p(3));
costFnc = @(p) sum((rhoHatFnc(p) - rho).^2);
A = [
    0,-1,0;
    0,0,-1;
    0,0,1];
c = [0;0;7];
p0 = [0.2;1;1];
[pHat,err] = fmincon(costFnc,p0,A,c);
rhoHat = rhoHatFnc(pHat);
return