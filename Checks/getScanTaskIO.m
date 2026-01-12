function [TaskIO,ScanDiag,sIdsToDo] = getScanTaskIO(G)

subjectIds = getSubjectIds(G);
TaskIO = [];
ScanDiagTable = [];
for iSubject = 1:numel(subjectIds)
    [T,S] = loadData(char(subjectIds(iSubject)));
    if ~isempty(T)
        TaskIO = [TaskIO;T]; %#ok<*AGROW>
        ScanDiagTable = [ScanDiagTable;S];
    end
end
scanTauResids = cell2mat(cellfun(@(r,s)r(s),...
    ScanDiagTable.resids,ScanDiagTable.scanSelect,...
    'UniformOutput',false));
ScanDiag.Table = ScanDiagTable;
ScanDiag.tauResids = scanTauResids;
sIdsToDo = subjectIds(...
    ~ismember(subjectIds,unique(TaskIO.SubjectId)));

return

function [TaskIO,ScanDiagnostics] = loadData(subjectId)
pathToData = sprintf('..%s..%sData%s%s%sBehavioural',...
    filesep,filesep,filesep,subjectId,filesep);
dirList = dir(pathToData);
fileList = {dirList.name}';
if any(strcmp(fileList,'ScanTaskIO.mat'))
    Data = load([pathToData,filesep,'ScanTaskIO.mat']);
    TaskIO = Data.TaskIO;
    ScanDiagnostics = ...
        [table(...
        repmat(categorical({subjectId}),size(Data.scanDiagnostics,1),1),...
        (1:size(Data.scanDiagnostics,1))',...
        'VariableNames',{'subjectId','run'}),...
        struct2table(Data.scanDiagnostics)];

else
    TaskIO = [];
    ScanDiagnostics = [];
end
return