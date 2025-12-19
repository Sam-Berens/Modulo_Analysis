function [] = zX0_makeQ()
dirs.Data = '../../Data';
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject);
    makeQ(subjectId)
end
return


function [] = makeQ(subjectId,G)
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl02 = [dirs.Alpha01,filesep,'Mdl02'];
dirs.G1 = [dirs.Mdl02,filesep,G];
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

srchIm = searchLight(tM,epiMask,@estimQ, 3);

B0 = srchIm(:,:,:,1);
B1 = srchIm(:,:,:,2);
Q = srchIm(:,:,:,2);
%these are the similarity scores for each distance bin
rH1 = srchIm(:,:,:,3);
rH2 = srchIm(:,:,:,4);
rH3 = srchIm(:,:,:,5);
err = srchIm(:,:,:,6);

%we're going to save all stats maps as nifti so we can norm q images
%to mni space and so that we can do QA checks on the others
images = {Q,B0,B1,rH1,rH2,rH3,err};
names = {'q','b0','b1','rH1','rH2','rH3','error'};
for iIm=1:numel(images)
    im.M = nan(size(epiMask.M,[1,2,3]));
    im.M(epiMask.idx) = images{iIm};
    im.V = epiMask.V; %TO DO - change data type to be approrpiate not binary
    im.V.dt(1) = 64;
    label = names{iIm};
    %save to their mdl02 folder
    im.V.fname = [dirs.Mdl02,filesep,label,'.nii'];
    im.V.descrip = 'Statistic map of nonlinear pattern–RDM fit model';
    spm_write_vol(im.V,im.M);
end
return

function [varargout] = estimQ(tPatterns)
tPatterns = reshape(tPatterns, [], size(tPatterns,4));
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
    [pHat, rhoHat, err] = fitNonlinMdl(I,P);
    allVals = [pHat(:); rhoHat(:); err(:)];
    varargout = num2cell(allVals);
else 
     varargout(1:7,1) = {nan};
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