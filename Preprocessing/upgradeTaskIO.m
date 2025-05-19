function [] = upgradeTaskIO(subjectId)

scriptsDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);
% Cd into subject dir
cd([subjectId,filesep,'Behavioural']);

%% Read in data, extract:
% taskSpikes: A timeseries object that records IOport pulses.
% scanSpikes: A collection of scan times (in seconds, not samples).
% pulseSpikes: A timeseries object that records detected heartbeats.
matFileNames = dir('*_202*.mat');
spikeData = getSpikeData(subjectId,matFileNames);
taskSpikes = getTaskSpikes(spikeData);
scanSpikes = getScanSpikes(spikeData);
pulseSpikes = getPulseSpikes(spikeData);

%% Produce a table of all heart beats
[heartEvents] = getHeartEvents(pulseSpikes);

%% Loop over runs
scanDiagnostics = struct;
for iRun = 1:5

    % Load the TaskIO etc.
    runSel = cellfun(@(s)~isempty(s),strfind({matFileNames.name}',...
        sprintf('_R%i_',iRun)));
    runData = load(matFileNames(runSel).name);

    % Create a hypothesised timeserise for expected taskSpikes (h)
    h = getExpetedTaskSpikes(taskSpikes,runData);

    % Compute the dot d the timings of each scoll button pressproduct between h and taskSpikes at each sample
    p = serialDot(taskSpikes.Data,h);

    % Find sample indicies that approximate the start and end of this run 
    [~,iRunStart] = max(p);
    iRunStart = iRunStart - 144;
    iRunEnd = iRunStart + numel(h) + 2*144;

    % Produce a table of all task events (stimuli and button presses)
    taskEvents = getTaskEvents(taskSpikes,runData,iRunStart,iRunEnd);

    % Get and filter the scan times from this run
    scanTimes = getScanTimes(scanSpikes,taskEvents,subjectId,iRun);
    [scanTimes,scanIdx,scanDiagnostics(iRun,1)] = trimScanTimes(scanTimes);

    % Use cubic interpolation to convert seconds into scanNum ...
    % ... for each task event.
    taskEvents.scanNum = cubicInterp(taskEvents.t,scanTimes,scanIdx);

    % Use cubic interpolation to estimate the heart phase ...
    % ... for each task event.
    taskEvents.heartPhase = mod(...
        cubicInterp(taskEvents.t,heartEvents.t,heartEvents.n),...
        1).*2*pi;

    % Extract and re-write TaskIO:
    % 1) Replace all the timings with the scanNum metric computed above...
    %    ... (including tArray);
    % 2) Add the timings of each scoll button press;
    % 3) Add response times;
    % 4) Add the heart phases for each event;
    % 5) Recode the responses;

end

%% Save each new TaskIO structure into a new table...
% ... with runs indexed by a specific variable

% Also save scanDiagnostics;

%% CD back
cd(scriptsDir);

return

function [spikeData] = getSpikeData(subjectId,matFileNames)
strIdx = strfind(matFileNames(1).name,'_202');
dateTime = matFileNames(1).name((strIdx+1):(strIdx+8));
spikeData = load([subjectId,'_',dateTime,'.mat']);
spikeData.DateTimeStart = datetime(spikeData.file.start');
return

function [taskSpikes] = getTaskSpikes(spikeData)
spikeDataFns = fieldnames(spikeData);
taskSpikes = spikeData.(...
    spikeDataFns{structfun(@(s)strcmp(s.title,'taskeven'),...
    rmfield(spikeData,{'file','DateTimeStart'}))});
si = taskSpikes.interval; % Sample interval
t = (0:(taskSpikes.length-1))' * si;
taskSpikes = timeseries(taskSpikes.values,t,'Name','TaskSpikes');
taskSpikes.DataInfo.Units = 'Volts';
taskSpikes = setuniformtime(taskSpikes,'Interval',si);
return

function [scanSpikes] = getScanSpikes(spikeData)
spikeDataFns = fieldnames(spikeData);
scanSpikes = spikeData.(...
    spikeDataFns{structfun(@(s)strcmp(s.title,'Scan Vol'),...
    rmfield(spikeData,{'file','DateTimeStart'}))}).times;
return

function [pulseSpikes] = getPulseSpikes(spikeData)
spikeDataFns = fieldnames(spikeData);
pulseSpikes = spikeData.(...
    spikeDataFns{structfun(@(s)strcmp(s.title,'pulseEvt'),...
    rmfield(spikeData,{'file','DateTimeStart'}))});
si = pulseSpikes.interval; % Sample interval
t = (0:(pulseSpikes.length-1))' * si;
pulseSpikes = timeseries(pulseSpikes.values,t,'Name','PulseSpikes');
pulseSpikes.DataInfo.Units = 'Volts';
pulseSpikes = setuniformtime(pulseSpikes,'Interval',si);
return

function [heartEvents] = getHeartEvents(pulseSpikes)
state = double(pulseSpikes.Data(1) > 1);
if state == 1
    count = 1;
end
c = pulseSpikes.TimeInfo.Increment;
jj = 0; % Index of each event
M = zeros(0,3);
for ii = 2:numel(pulseSpikes.Data)
    if (state==0) && (pulseSpikes.Data(ii) > 1)
        state = 1;
        count = 1;
        iStartEvn = ii;
    elseif (state==1) && pulseSpikes.Data(ii) > 1
        count = count + 1;
    elseif (state==1) && (pulseSpikes.Data(ii) < 1)
        state = 0;
        dur = count * c;
        jj = jj + 1;
        M(jj,1) = dur;
        M(jj,2) = iStartEvn;
        M(jj,3) = (iStartEvn-1)*c;
    end
end
dBeat = diff(M(:,3));
dBeat = [dBeat;NaN];
filter = dBeat < 0.4; % Filter out beat faster than 150BPM
M = M(~filter,:);
M = M(:,2:3);
heartEvents = array2table(M,'VariableNames',{'SampleIdx','t'});
N = size(heartEvents,1);
heartEvents.n = (0:(N-1))';
return

function [h] = getExpetedTaskSpikes(taskSpikes,runData)
tS = runData.TaskIO(1).tShowA;
tE = runData.TaskIO(end).tRespo;
dur = tE - tS;
nSamps = ceil(dur / taskSpikes.TimeInfo.Increment);
h = zeros(nSamps,1);
for iTrial = 1:numel(runData.TaskIO)

    a = runData.TaskIO(iTrial).a;
    tAs = runData.TaskIO(iTrial).tShowA - tS; % Start time (secs from t=0)
    if ~isnan(tAs)
        tAe = tAs + (a+9)*runData.globals.portUnitLength;
        iAs = round(tAs / taskSpikes.TimeInfo.Increment) + 1;
        iAe = round(tAe / taskSpikes.TimeInfo.Increment) + 1;
        h(iAs:iAe) = 1;
    end

    b = runData.TaskIO(iTrial).b;
    tBs = runData.TaskIO(iTrial).tShowB - tS; % Start time (secs from t=0)
    if ~isnan(tBs)
        tBe = tBs + (b+17)*runData.globals.portUnitLength;
        iBs = round(tBs / taskSpikes.TimeInfo.Increment) + 1;
        iBe = round(tBe / taskSpikes.TimeInfo.Increment) + 1;
        h(iBs:iBe) = 1;
    end

    tRs = runData.TaskIO(iTrial).tRespo - tS;
    if ~isnan(tRs)
        tRe = tRs + 2*runData.globals.portUnitLength;
        iRs = round(tRs / taskSpikes.TimeInfo.Increment) + 1;
        iRe = round(tRe / taskSpikes.TimeInfo.Increment) + 1;
        h(iRs:iRe) = 1;
    end

end
return

function [p] = serialDot(v,k)
v = v./norm(v);
k = k./norm(k);
m = numel(k);
n = numel(v)-m+1;
p = zeros(n,1);
parfor ii = 1:n
    s = v(ii:(ii+m-1)); %#ok<PFBNS>
    p(ii) = s'*k;
end
return

function [taskEvents] = getTaskEvents(taskSpikes,runData,iRunStart,iRunEnd)
state = 0;
c = taskSpikes.TimeInfo.Increment / runData.globals.portUnitLength;
jj = 0; % Index of each event
M = zeros(0,3);
for ii = iRunStart:iRunEnd
    if (state==0) && (taskSpikes.Data(ii) > 1)
        state = 1;
        count = 1;
        iStartEvn = ii;
    elseif (state==1) && taskSpikes.Data(ii) > 1
        count = count + 1;
    elseif (state==1) && (taskSpikes.Data(ii) < 1)
        state = 0;
        id = round(count * c);
        jj = jj + 1;
        M(jj,1) = id;
        M(jj,2) = iStartEvn;
        M(jj,3) = (iStartEvn-1)*taskSpikes.TimeInfo.Increment;
    end
end
taskEvents = array2table(M,'VariableNames',{'EventId','SampleIdx','t'});
return

function [scanTimes] = getScanTimes(scanSpikes,taskEvents,subjectId,iRun)
D = scanSpikes'-scanSpikes;
D(tril(true(numel(scanSpikes),numel(scanSpikes)))) = NaN;
L = abs(D-(292*2.2))<1e-2;
[iStart,iEnd] = ind2sub(size(L),find(L));
tStart = scanSpikes(iStart);
tEnd = scanSpikes(iEnd);
iScanSeq = NaN;
for ii = 1:numel(tStart)
    if (taskEvents.t(1)>tStart(ii)) && (taskEvents.t(end)<tEnd(ii))
        iScanSeq = ii;
    end
end
if isnan(iScanSeq)
    error('Unable to find a scan sequence for run %i and subject %s.',...
        iRun,subjectId);
end
scanTimes = scanSpikes(iStart(iScanSeq):iEnd(iScanSeq));
return

function [scanTimes,scanIdx,diagnostics] = trimScanTimes(scanTimes)

% Fit the scan observation model
startParams = [0;2.2;1e-3];
opts = optimoptions(@fmincon,'Display','iter');
problem = createOptimProblem('fmincon',...
    'objective',@(p)func2min(p,scanTimes),...
    'x0',startParams,...
    'lb',[-1.1;2.1;1e-8],...
    'ub',[+1.1;2.3;1e-1],...
    'options',opts);
ms = MultiStart;
params = run(ms,problem,50);
[scanIdx,resids] = getScanIdx(params,scanTimes);

% Remove scans with a large residual
scanSelect = abs(resids) < 1e-2;
scanTimes = scanTimes(scanSelect);
scanIdx = scanIdx(scanSelect);

diagnostics.params = params;
diagnostics.scanSelect = scanSelect;
diagnostics.resids = resids;
return

function [nll] = func2min(params,scanTimes)
b0 = params(1);
b1 = params(2);
sigma = params(3);
cauchy = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',1);
d = (scanTimes-min(scanTimes)+b0)./b1;
resids = computeResids(d);
nll = -sum(log(pdf(cauchy,resids)));
return

function [idx,resids] = getScanIdx(params,scanTimes)
b0 = params(1);
b1 = params(2);
d = (scanTimes-min(scanTimes)+b0)./b1;
idx = round(d);
resids = computeResids(d);
return

function [resids] = computeResids(d)
resids = mod(d,1);
signedResids = [resids-1,resids];
unsinedResids = abs(signedResids);
[~,sel] = min(unsinedResids,[],2);
resids = nan(size(resids));
for ii = 1:numel(resids)
    resids(ii) = signedResids(ii,sel(ii));
end
return

function [vy] = cubicInterp(vx,x,y)
vy = nan(size(vx));
for ii = 1:numel(vx)
    [~,index] = sort(abs(x-vx(ii)));
    sindex = sort(index(1:4));
    sx = x(sindex);
    [sz,mu,sigma] = zscore(sx);
    vz = (vx(ii)-mu)/sigma;
    sy = y(sindex);
    p = polyfit(sz,sy,3);
    vy(ii) = polyval(p,vz);
end
return