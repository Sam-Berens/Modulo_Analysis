function [] = z00_makeEventSpec(G)

% Set some constants
tr = 2.2;
stimIds = (0:6)';
stimDur = 3/tr; % Duration of Spark stimuli
decisiDur = 6/tr; % Duration of Decision period
dirs.Data = dir(['..',filesep,'..',filesep,'Data']);
subjectIds = getSubjectIds(G);

% Subject loop
for iSubject = 1:numel(subjectIds)

    % Set required directory paths
    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Behav = [dirs.Subject,filesep,'Behavioural'];
    dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];

    % Get ScanTaskIO
    try
        STIO = load([dirs.Behav,filesep,'ScanTaskIO.mat']);
    catch
        error('Cannot load ScanTaskIO.mat for %s!',cId);
    end

    % Loop through stimuli and runs
    TaskIO = STIO.TaskIO;
    runIds = unique(TaskIO.iRun)';
    for iStim = stimIds
        signalName = sprintf('i%i',iStim);

        % Set and create the output directory (if needed)
        dirs.Target = [dirs.Alpha00,filesep,signalName];
        if ~exist(dirs.Target, "dir")
            mkdir(dirs.Target);
        end

        for iRun = runIds
            % Specify the names (e.g., [i0, ###, Decision])
            names = {signalName,'###','Decision'};

            % Specify the onsets (tauSignal, tauResidu, tauDecision)
            T = TaskIO(TaskIO.iRun==iRun,:);
            tauAB = [T.tauShowA,T.tauShowB];
            sAB = [T.a==iStim,T.b==iStim];
            sNan = isnan(tauAB);
            tauSignal = sort(tauAB(sAB));
            tauResidu = sort(tauAB((~sAB)&(~sNan)));
            tauDecisi = T.tauArray(~isnan(T.tauArray));
            onsets = {tauSignal,tauResidu,tauDecisi};

            % Specify the durations [stimDur, stimDur, decisiDur]
            durations = {...
                ones(size(tauSignal)).*stimDur,...
                ones(size(tauResidu)).*stimDur,...
                ones(size(tauDecisi)).*decisiDur};

            % Save the EventSpec
            outFn = sprintf('%s%sEventSpec_R%i.mat',...
                dirs.Target,filesep,iRun);
            save(outFn,'names','onsets','durations');
        end
    end
end
return