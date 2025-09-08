function []= makeEventSpec_Alpha00(groupFunc)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;

subjectIds = groupFunc(); %this will be e.g. getG1()
nSubs = numel(subjectIds);

for ii=1:nSubs
    subDir = [dataDir,filesep,subjectIds{ii}];
    anlDir = [subDir,filesep,'Analysis'];
    alphDir = [anlDir,filesep,'Alpha00'];
    if ~exist(alphDir, "dir")
        mkdir(alphDir);
    end
    %get TaskIO
    behaDir = [subDir,filesep,'Behavioural'];
    % Load the ScanTaskIO
    try
        X = load([behaDir,filesep,'ScanTaskIO.mat']);
    catch
        fprintf(['Either you specified the filename of TaskIO wrong or'...
            'the file is corrupt/too big \n Filename:%s\n'],[behaDir,filesep,'ScanTaskIO.mat'])
        return
    end
    taskIO = X.TaskIO;
    stimDur = 3/2.2;
    stims = 0:5;
    runs = unique(taskIO.iRun)';
    for iRun = runs
        %% Stim loop
        for iStim = 0:5
            %probs unecess but just clear variables to be sure
            onsets =  [];%maximum num rows should be 
            durations = [];
            names = {[],[]};

            names = {sprintf('i%i',iStim),'other_stims'}; %change to better name?
            tau = [...
                taskIO.tauShowA((taskIO.a==iStim)&(taskIO.iRun==iRun));
                taskIO.tauShowB((taskIO.b==iStim)&(taskIO.iRun==iRun))];
            tau = sort(tau); %so this is all the onsents for the rgi in this run
         
            nuis = stims(~ismember(stims,iStim))';
            decisTau = taskIO.tauArray((taskIO.iRun==iRun) & ~isnan(taskIO.tauArray));
            nuisTau = [...
                taskIO.tauShowA((ismember(taskIO.a,nuis))&(taskIO.iRun==iRun)); %all other stims at A
                taskIO.tauShowB((ismember(taskIO.b,nuis))&(taskIO.iRun==iRun)); %all other stims at B
                decisTau; %all decision periods
                ];

            %preacllocate mats
            onsets = cell(1,2);
            durations = cell(1,2);
            %load in values for RGI
            onsets{1,1} = tau; %load in the scan times for regressor of interest
            durations{1,1} = ones(size(tau))*stimDur;
            %load in values for RGI and nuis
            onsets{1,2} = sort(nuisTau);
            nuisDur = ones(size(nuisTau))*stimDur;
            nuisDur(ismember(sort(nuisTau), decisTau)) = 6/2.2; %change the durs for decision periods
            durations{1,2} =nuisDur;

            %% save event spec for 1 run and 1 stim at a time
            folderName = [alphaDir,filesep,sprintf('i%i_R%i', iStim, runs(iRun))];
            if ~exist(folderName, "dir")
                mkdir(folderName);
            end
            eventSpecFn = [folderName, filesep, 'EventSpec.mat'];
            save(eventSpecFn, "names","durations","onsets", "-mat"); %does order of saving matter?


        end

    end
end

return

