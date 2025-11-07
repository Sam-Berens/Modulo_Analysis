function [mdl] = fitV1Mdl(G,subjectId,roiPattern,doSVD,preW,tauOffset)
% FITV1MDL Fit a simple visual-stimulus GLM to native EPI data across runs.
%
%   mdl = FITV1MDL(G, subjectId, roiPattern, doSVD, preW, tauOffset)
%
% DESCRIPTION
%   Loads run-wise EPI time series from a native-space ROI, builds a design
%   matrix for each run/offset, optionally extracts the 1st PC (SVD) from
%   the ROI, applies optional pre-whitening (SPM REML with AR(1) or FAST
%   covariance components), and fits an OLS model. Returns per-run,
%   per-offset estimates and diagnostics.
%
% INPUTS (all optional; defaults in brackets)
%   G           : Experiment/Group identifier used by getNativeRoiPath
%   subjectId   : Subject identifier subfolder under Data/
%   roiPattern  : ROI name/pattern understood by getNativeRoiPath
%   doSVD       : If true, replace ROI voxel matrix with 1st PC timecourse
%   preW        : Pre-whitening mode: 0=none, 1=AR(1), 2=FAST
%   tauOffset   : Scalar/vector of TR offsets (samples) to shift X by
%
% FILES / FOLDERS (relative to ../../Data/<subjectId>/)
%   Behavioural/ScanTaskIO.mat        -> TaskIO struct
%   EPI/RP*.mat                       -> motion (cell R with one per run)
%   EPI/*_epiMask00.nii               -> brain/ROI mask in native space
%   EPI/2_Temporal/R*/<*.nii>         -> per-run 4D EPI split to 3D volumes
%
% OUTPUT
%   mdl(iRun,iOff) struct with fields:
%     X, Y                       : design and data (time x channels)
%     W, WX, WY (if preW>0)      : whitening operator and transformed mats
%     B, C                       : OLS betas and covariance (k x 1 x ch)
%     P, M                       : hat and residual-forming matrices
%     df1, df2, mse, ll          : dof, per-channel MSE, Gaussian loglik
%     Err or WErr                : residuals (pre- or post-whitened)
%
% DEPENDENCIES (external)
%   getNativeRoiPath, getNativeData, getX_V1
%   SPM: spm_reml, spm_sqrtm, spm_inv
%
% NOTES
%   • TR is assumed to be 2.2 s (change 'tr' below to override).
%   • With doSVD=true, the ROI is summarised to a single timecourse via the
%     first principal component (u*s*mean(v) from svds).
%   • preW==1 uses an AR(1)-like component dictionary; preW==2 uses FAST.
%
% EXAMPLE
%   mdl = fitV1Mdl('G1','329e92e2','R001',true,1,0:2);
%
% SEE ALSO: svds, pinv, spm_reml

%% Set default inputs
if (nargin < 1) || isempty(G)
    G = 'G1';
end
if (nargin < 2) || isempty(subjectId)
    subjectId = '329e92e2';
end
if (nargin < 3) || isempty(roiPattern)
    roiPattern = 'R001';
end
if (nargin < 4) || isempty(doSVD)
    doSVD = true;
end
if (nargin < 5) || isempty(preW)
    preW = 1;
end
if (nargin < 6) || isempty(tauOffset)
    tauOffset = 0;
end

%% Constants
tr = 2.2;
% For AR1: maxQ dictates the lookahead time for temporal correlations
% For FAST: fastP dictates the number of temporal correlation components
maxQ = 32;
fastP = 6;

%%

dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.Bhv = [dirs.Subject,filesep,'Behavioural'];
dirs.EPI = [dirs.Subject,filesep,'EPI'];
dirs.Y = [dirs.EPI,filesep,'2_Temporal'];

TaskIO = getTaskIO(dirs.Bhv);
R = getR(dirs.EPI);
maskFn = dir([dirs.EPI,filesep,'_*_epiMask00.nii']);
maskFn = [maskFn.folder,filesep,maskFn.name];
epiFns = getEpiFns(dirs.Y);
roiPath = getNativeRoiPath(G,subjectId,roiPattern);

maxRun = numel(R);
mdl = struct();
for iRun = 1:maxRun

    Y = getNativeData(maskFn,epiFns{iRun},roiPath);
    n = size(Y,1);
    tau = (0:n-1)';
    switch preW
        case 0
            % No whitening
        case 1
            Q = GetCovComps_AR1(tr, tau, maxQ);
        case 2
            Q = GetCovComps_FAST(tr, tau, fastP);
        otherwise
            error(['preW must be 0 (none), 1 (AR1), or 2 (FAST).%c',...
                'Got %g.'],10,preW);
    end

    %% Extract the first principle component if requested
    if doSVD
        [u,s,v] = svds(Y,1);
        y = u*s*mean(v,1)';
        Y = y;
    end

    %%
    for iOff = 1:numel(tauOffset)

        X = getX_V1(TaskIO,R{iRun},iRun,tauOffset(iOff));

        %% Pre-whitening
        if preW == 0
            W = eye(size(Y,1));
        else
            YY = Y*Y';
            V = spm_reml(YY,X,Q);
            W = spm_sqrtm(spm_inv(V));
            %wScale = full(sum(W,'all'))/size(W,1);
            wScale = mean(diag(W));
            W = W./wScale;
        end
        WY = W*Y;
        WX = W*X;

        %% OLS
        B = pinv(WX) * WY;
        IXX = pinv(WX'*WX);
        P = WX * pinv(WX);
        M = eye(size(P,1)) - P;
        Err = M * WY;
        df1 = trace(P);
        df2 = trace(M);
        mse = diag(Err'*Err)./df2;
        C = IXX.*reshape(mse,1,1,[]);
        % ll = sum(log(normpdf(...
        %     Err,...
        %     zeros(n,numel(mse)),....
        %     repmat(sqrt(mse'),n,1)...
        %     )),1)';
        sse = sum(Err.^2, 1)';
        ll  = -0.5*(n*log(2*pi*mse) + sse./mse);

        %% Save
        mdl(iRun,iOff).X = X;
        mdl(iRun,iOff).Y = Y;
        if preW > 0
            mdl(iRun,iOff).W = W;
            mdl(iRun,iOff).WX = WX;
            mdl(iRun,iOff).WY = WY;
        end
        mdl(iRun,iOff).df1 = df1;
        mdl(iRun,iOff).df2 = df2;
        mdl(iRun,iOff).B = B;
        mdl(iRun,iOff).C = C;
        mdl(iRun,iOff).P = P;
        mdl(iRun,iOff).M = M;
        if preW == 0
            mdl(iRun,iOff).Err = Err;
        else
            mdl(iRun,iOff).WErr = Err;
        end
        mdl(iRun,iOff).mse = mse;
        mdl(iRun,iOff).ll = ll;
    end

end

return

function [TaskIO] = getTaskIO(path2data)
fullpath = @(s)[s.folder,filesep,s.name];
fileList  = dir([path2data,filesep,'ScanTaskIO.mat']);
fileList = arrayfun(fullpath,fileList,'UniformOutput',false);
syncData = load(fileList{1});
TaskIO = syncData.TaskIO;
return

function [R] = getR(path2data)
fullpath = @(s)[s.folder,filesep,s.name];
fileList  = dir([path2data,filesep,'RP*.mat']);
[~,ord] = sort({fileList.name});
fileList = fileList(ord);
fileList = arrayfun(fullpath,fileList,'UniformOutput',false);
R = cell(size(fileList));
for iRun = 1:numel(fileList)
    RPs = load(fileList{iRun});
    R{iRun} = RPs.R;
end
return

function [epiFns] = getEpiFns(path2data)
fullpath = @(s)[s.folder,filesep,s.name];
runList  = dir([path2data,filesep,'R*']);
[~,ord]  = sort({runList.name});
runList  = runList(ord);
runList = arrayfun(fullpath,runList,'UniformOutput',false);
epiFns = cell(size(runList));
for iRun = 1:numel(runList)
    vols = dir([runList{iRun},filesep,'*.nii']);
    [~,vord] = sort({vols.name});
    vols = vols(vord);
    epiFns{iRun} = arrayfun(fullpath,vols,'UniformOutput',false);
end
return

function [Q] = GetCovComps_AR1(tr,tau,maxQ)
% Q = GetCovComps_AR1()
% Get a dictionary of covariance components for the AR(1) model.

ts = tau*tr; % Time in seconds

% e: Time constants (seconds)
e = 2.^((floor(log2(tr/4)):log2(maxQ))');

Q = cell(1,numel(e)*2);
for ii = 1:numel(e)
    for jj = 0:1
        Q{(ii-1)*2 + jj + 1} = toeplitz((ts.^jj).*exp(-ts/e(ii)));
    end
end
return

function [Q] = GetCovComps_FAST(tr,tau,p)
% Q = GetCovComps_FAST()
% Get a dictionary of covariance components for the FAST model.
% See https://doi.org/10.1002/hbm.24218
% By Sam Berens (s.berens@sussex.ac.uk)
% I know it's not that elegant but it helped me understand.

nScans = numel(tau);
q = (1:p)';
alpha = 8./(2.^q);

ts = tau*tr; % Time in seconds
Q = cell(1,p*3);
for ialpha = 1:p
    for n = 0:2
        C = nan(nScans);
        for ii = 1:nScans
            for jj = 1:nScans
                if (ii == jj) && (n==0)
                    C(ii,jj) = 1;
                else
                    C(ii,jj) = abs(ts(jj)-ts(ii))^n * ...
                        exp(-alpha(ialpha)*abs(ts(jj)-ts(ii)));
                end
            end
        end
        Q{((ialpha-1)*3)+n+1} = C;
    end
end
return