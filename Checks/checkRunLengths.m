function [runLengths] = checkRunLengths()

cd(['..',filesep,'..',filesep,'Data']);

dirList = dir();
subjectIds = {dirList(cellfun(@(s)numel(s)==8,{dirList.name}')).name}';
runLengths = cell(size(subjectIds));

for iSubject = 1:numel(subjectIds)
    cd(subjectIds{iSubject});
    cd Behavioural/;

    matFileNames = dir('*_202*.mat');
    spikeData = getSpikeData(subjectIds{iSubject},matFileNames);
    tScans = getScanTimes(spikeData);

    dtScans = diff(tScans);
    isSameRun = double(abs(dtScans-2.2) < 1e-2);
    iRun = 1;
    runLengths{iSubject} = 1;
    for ii = 1:numel(isSameRun)
        if isSameRun(ii)
            runLengths{iSubject}(iRun,1) = runLengths{iSubject}(iRun,1) + 1;
        else
            iRun = iRun + 1;
            runLengths{iSubject}(iRun,1) = 1;
        end
    end

    cd ../..;

end

runLengths = cell2struct(runLengths,cellfun(@(s)['s_',s],subjectIds,'UniformOutput',false));

return


function [spikeData] = getSpikeData(subjectId,matFileNames)
strIdx = strfind(matFileNames(1).name,'_202');
dateTime = matFileNames(1).name((strIdx+1):(strIdx+8));
spikeData = load([subjectId,'_',dateTime,'.mat']);
spikeData.DateTimeStart = datetime(spikeData.file.start');
return

function [scanTimes] = getScanTimes(spikeData)
spikeDataFns = fieldnames(spikeData);
scanVols = spikeData.(...
    spikeDataFns{structfun(@(s)strcmp(s.title,'Scan Vol'),...
    rmfield(spikeData,{'file','DateTimeStart'}))});
scanTimes = scanVols.times;
return