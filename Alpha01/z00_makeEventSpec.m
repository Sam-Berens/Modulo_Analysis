function [] = z00_makeEventSpec(G)
% Generate and save first-level model specification files for Alpha01.
%
%   z00_makeEventSpec(G) creates first-level model specification files for
%   Alpha01 across all subjects in the group defined by the input groupId
%   (G). These files contain stimulus names, onsets, and durations,
%   organised to allow for Least-Squares-Separate decoding in SPM.
%
%   The function performs the following steps:
%     1. Loads behavioural task data (ScanTaskIO.mat) for each subject.
%     2. Iterates through all stimuli (0–5), both positions (A and B),
%        and task runs.
%     3. Extracts event onsets for:
%        - The stimulus of interest (now split by position; e.g. "a0"/"b0" ... "a5"/"b5")
%        - Residual stimuli ("###")
%        - Decision periods ("Decision")
%     4. Assigns durations in TR units:
%        - Spark stimuli: 3/tr
%        - Decision periods: 6/tr
%     5. Saves an EventSpec_R*.mat file per run containing:
%        • names      - cell array of event names
%        • onsets     - cell array of onset times (in TR units)
%        • durations  - cell array of durations (in TR units)
%
%   Important difference from Alpha00:
%     - Alpha00 combined A- and B-position occurrences into a single regressor
%       for each stimulus. Alpha01 instead creates regressors of interest that
%       include only the onsets/durations for that stimulus when shown in a
%       specific position (either A or B).
%
%   Outputs:
%     EventSpec_R*.mat files are saved in:
%       */[SubjectId]/Analysis/Alpha01/*(from a,b)*(from 1:5)/EventSpec_R*.mat


% Set some constants
tr = 2.2;
stimIds = (0:5)';
stimDur = 3/tr; % Duration of Spark stimuli
decisiDur = 6/tr; % Duration of Decision period
dirs.Data = ['..',filesep,'..',filesep,'Data'];
subjectIds = getSubjectIds(G);

% Subject loop
for iSubject = 1:numel(subjectIds)

    % Set required directory paths
    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Behav = [dirs.Subject,filesep,'Behavioural'];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

    % Get ScanTaskIO
    try
        STIO = load([dirs.Behav,filesep,'ScanTaskIO.mat']);
    catch
        error('Cannot load ScanTaskIO.mat for %s!',cId);
    end

    % Loop through stimuli and runs
    TaskIO = STIO.TaskIO;
    runIds = unique(TaskIO.iRun)';
    for iStim = stimIds'
        %positions to loop through
        position = 'AB';
        for iPos=1:2
            X = position(iPos);
            x = lower(position(iPos));
            %rgi stands for regressor of intrest
            rgi = sprintf('%s%i',x,iStim);
            % Set and create the output directory (if needed)
            dirs.Target = [dirs.Alpha01,filesep,rgi];
            if ~exist(dirs.Target, "dir")
                mkdir(dirs.Target);
            end
            for iRun = runIds
                % Specify the names (e.g., [i0, ###, Decision])
                names = {rgi,'###','Decision'};

                % Specify the onsets (tauSignal, tauResidu, tauDecision)
                T = TaskIO(TaskIO.iRun==iRun,:);
                %this is either tauShowA or tauShow B
                tauRgi = [T.(['tauShow',X])];
                tauAB = [T.tauShowA,T.tauShowB];
                sRgi = T.(x)==iStim;
                tauSignal =  sort(tauRgi(sRgi));
                %A is the 1st column and B is the 2nd column
                otherPos = mod(iPos,2)+1;
                tauResidu = sort([tauAB(~sRgi, iPos); tauAB(:, otherPos)]);
                tauResidu = tauResidu(~isnan(tauResidu));
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
end
return