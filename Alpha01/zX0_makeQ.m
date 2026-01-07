function [] = zX0_makeQ()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
subjectIds = char(subjectIds);
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject,:);
    makeQ(subjectId,G);
end
return


function [] = makeQ(subjectId,G)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl02 = [dirs.Alpha01,filesep,'Mdl02'];
dirs.G1 = [dirs.Mdl02,filesep,G];
if ~exist(dirs.Mdl02 ,'dir')
    mkdir(dirs.Mdl02);
end

epiMask = getEpiMask(subjectId);

%% Get a matrix of t-stats for all conditions
%remember that a stims are stacked ontop of b stim
[tImgs] = getTimgs(subjectId,epiMask);

%% loop through searchlight centres to test produce q term from precursor to mdl02
[Y,N] = searchlight3D(3,@estimQ,epiMask,tImgs); %reminder N is the volume of each searchlight
clear tImgs; %this is to save memory since we dont need it anymore
%we're going to save all stats maps as nifti so we can norm q images
%to mni space and so that we can do QA checks on the others
B0 = srchIm(:,:,:,1);
B1 = srchIm(:,:,:,2);
Q = srchIm(:,:,:,3);
%these are the similarity scores for each distance bin
rH1 = srchIm(:,:,:,4);
rH2 = srchIm(:,:,:,5);
rH3 = srchIm(:,:,:,6);
%error of the model
err = srchIm(:,:,:,7);

images = {Q,B0,B1,rH1,rH2,rH3,err};
names = {'q','b0','b1','rH1','rH2','rH3','error'};
for iIm=1:numel(images)
    im.M = nan(size(epiMask.M,[1,2,3]));
    im.M(epiMask.idx) = images{iIm};
    im.V = epiMask.V;
    im.V.dt(1) = 64;
    label = names{iIm};
    %save to their mdl02 folder
    im.V.fname = [dirs.Mdl02,filesep,label,'.nii'];
    im.V.descrip = 'Statistic map of nonlinear pattern–RDM fit model';
    spm_write_vol(im.V,im.M);
end
return

function [allVals] = estimQ(tPatterns)
%% expects flattened tPatterns which are [linIdx,condition]
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
    
else 
     allVals = nan(7,1);
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

opts = optimoptions('fmincon','Display','off','Algorithm','sqp');
[pHat,err] = fmincon(costFnc,p0,A,c,[],[],[],[],[],opts);
rhoHat = rhoHatFnc(pHat);
return
