function [X] = getX_V1(TaskIO,RPs,iRun,tauOffset,filterPeriod)
%
%

%% Default values
if nargin < 2 || isempty(TaskIO) || isempty(RPs)
    error('Feed me!');
end

if nargin < 3 || isempty(iRun)
    iRun = 1;
end

if nargin < 4 || isempty(tauOffset)
    tauOffset = 0;
end

if nargin < 5 || isempty(filterPeriod)
    filterPeriod = 128; % Seconds
end

%% Filter for a specific run
TaskIO = TaskIO(TaskIO.iRun == iRun,:);

%% Set some constants
n = size(RPs,1);
tr = 2.2;
microRes = 12^2;
m = n*microRes; % Number of microtime states
microTimeOnset = 1056/2200; % microTimeOnset (in factions of tau)

%% Preallocate X
X = zeros(m,3);

%% Sparks
onset1 = sort([TaskIO.tauShowA;TaskIO.tauShowB]) + tauOffset;
onset1 = onset1(~isnan(onset1));
idx1 = round(onset1*microRes + (0:(((3/tr)*microRes)-1)));
idx1 = idx1(:);
idx1 = idx1((idx1>0)&(idx1<=m));
X(idx1,1) = tr/microRes;

%% Symbol
onset2 = sort(TaskIO.tauShowB)+(3/tr) + tauOffset;
onset2 = onset2(~isnan(onset2));
idx2 = round(onset2*microRes + (0:(((1/tr)*microRes)-1)));
idx2 = idx2(:);
idx2 = idx2((idx2>0)&(idx2<=m));
X(idx2,2) = tr/microRes;

%% Response period
onset3 = sort(TaskIO.tauArray)+ tauOffset;
onset3 = onset3(~isnan(onset3));
idx3 = round(onset3*microRes + (0:(((6/tr)*microRes)-1)));
idx3 = idx3(:);
idx3 = idx3((idx3>0)&(idx3<=m));
X(idx3,3) = tr/microRes;

%% Convolution
hrf = spm_hrf(tr/microRes);

% Normalise the amplitude of the hrf to have a unit peak
hrf = hrf ./ max(hrf);

XX = zeros(m+numel(hrf)-1,3);
for iX = 1:3
    XX(:,iX) = conv(X(:,iX),hrf);
end

% Trim
X = XX(1:m,:);
clear XX;

%% Downsample
tau = linspace(...
    0 - microTimeOnset,... Lower limit
    n - microTimeOnset,... Upper limit
    m)';
[~,idxDown] = min((tau-(0:n-1)).^2,[],1);
idxDown = idxDown';
X = X(idxDown,:);

%% Add RPs
X = [X,RPs];

%% Add Cosine basis set 
 k = fix((2*n*tr/filterPeriod)+1); % Order of high-pass filter
 CosBasis = spm_dctmtx(n,k);
 CosBasis = flip(CosBasis(:,2:end),2);
 CosBasis = CosBasis./max(CosBasis,[],1);
 X = [X,CosBasis,ones(n,1)];

return