function [] = z00_makeEventSpec(G)
% Generate and save first-level model specification files for A02.
%
%   z00_makeEventSpec(G) creates first-level model specification files for
%   A02 across all subjects in the group defined by the input groupId (G).
%   These files contain stimulus names, onsets, and durations.
%
%   The function performs the following steps:
%     1. Loads behavioural task data (ScanTaskIO.mat) for each subject.
%     2. Extracts event onsets for:
%        - All stimuli in the "a" position
%        - All stimuli in the "b" position
%        - Decision periods ("Decision")
%     3. Assigns durations in TR units:
%        - Spark stimuli: 3/tr
%        - Decision periods: 6/tr
%     4. Saves an EventSpec_R*.mat file per run containing:
%        • names      - cell array of event names
%        • onsets     - cell array of onset times (in TR units)
%        • durations  - cell array of durations (in TR units)
%
%   Outputs:
%     EventSpec_R*.mat files are saved in:
%       */[SubjectId]/Analysis/A02/EventSpec_R*.mat


% Set some constants
tr = 2.2;
stimDur = 3/tr; % Duration of Spark stimuli
decisiDur = 6/tr; % Duration of Decision period
dirs.Data = ['..',filesep,'..',filesep,'Data'];
subjectIds = getSubjectIds(G);

% Subject loop
for iSubject = 1:numel(subjectIds)

    % Set required directory paths
    cId = char(subjectIds(iSubject));
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Behav = [dirs.Subject,filesep,'Behavioural'];
    dirs.A02 = [dirs.Subject,filesep,'Analysis',filesep,'A02'];

    % Get ScanTaskIO
    try
        STIO = load([dirs.Behav,filesep,'ScanTaskIO.mat']);
    catch
        error('Cannot load ScanTaskIO.mat for %s!',cId);
    end

    % Create the output directory (if needed)
    if ~exist(dirs.A02, "dir")
        mkdir(dirs.A02);
    end

    % Loop through stimuli and runs
    TaskIO = STIO.TaskIO;
    runIds = unique(TaskIO.iRun)';

    for iRun = runIds
        % Specify the names
        names = {'a','b','Decision'};

        % Specify the onsets (tauA, tauB, tauDecision)
        T = TaskIO(TaskIO.iRun==iRun,:);
        tauA = T.tauShowA(~isnan(T.tauShowA));
        tauB = T.tauShowB(~isnan(T.tauShowB));
        tauDecisi = T.tauArray(~isnan(T.tauArray));
        onsets = {tauA,tauB,tauDecisi};

        % Specify the durations [stimDur, stimDur, decisiDur]
        durations = {...
            ones(size(tauA)).*stimDur,...
            ones(size(tauB)).*stimDur,...
            ones(size(tauDecisi)).*decisiDur};

        % Save the EventSpec
        outFn = sprintf('%s%sEventSpec_R%i.mat',...
            dirs.A02,filesep,iRun);
        save(outFn,'names','onsets','durations');
    end
end
return