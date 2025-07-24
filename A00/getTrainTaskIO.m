function [TaskIO] = getTrainTaskIO(SubjectId)
% getTrainTaskIO.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/07/2025
%
% Syntax:  Data = getTrainTaskIO(SubjectId)
%
% Description:
%    Retrieves and processes subject-specific task data from a web server
%    (or local file if avaible). The function downloads TaskIO data,
%    converts the structure to a table, assigns proper data types to each
%    field, and organizes trials based on predefined pair type indices. It
%    also creates unique trial identifiers.
%
% Inputs:
%    SubjectId - Character vector (or string) specifying the subject
%                identifier.
%
% Outputs:
%    Data - A table containing the processed subject data with the
%           following variables:
%               SubjectId, SessionId, TrialId, tSup, PairId, PairType, 
%               TrialType, FieldIdx_A, FieldIdx_B, FieldIdx_C, FieldIdx_R,
%               and RT.
%
% Example:
%    Data = getTrainTaskIO('subject123');
%

%% Check to see whether data has been saved locally
pathToTaskIO = sprintf('..%s..%sData%s%s%sBehavioural',...
    filesep,filesep,filesep,SubjectId,filesep);
dirList = dir(pathToTaskIO);
fileList = {dirList.name}';
if any(contains(fileList,'TrainTaskIO.mat'))
    TaskIO = load([pathToTaskIO,filesep,'TrainTaskIO.mat']);
    TaskIO = TaskIO.TaskIO;
    return
end

%% Get TaskIO from the web-server
TaskIO = webwrite(...
    'https://b01.learningandinference.org/GetTaskIO.php',...
    'SubjectId',SubjectId);
TaskIO = struct2table(TaskIO);
TaskIO.SubjectId = categorical(TaskIO.SubjectId);
TaskIO.FieldSize = [];
TaskIO.SessionId = str2double(TaskIO.SessionId);
TaskIO.TrialId = str2double(TaskIO.TrialId);
TaskIO.PairId = str2double(TaskIO.PairId);
TaskIO.TrialType = categorical(TaskIO.TrialType);
TaskIO.OppId = [];
TaskIO.FieldIdx_A = str2double(TaskIO.FieldIdx_A);
TaskIO.FieldIdx_B = str2double(TaskIO.FieldIdx_B);
TaskIO.FieldIdx_C = str2double(TaskIO.FieldIdx_C);
TaskIO.AttemptNum = str2double(TaskIO.AttemptNum);
TaskIO.FieldIdx_R = str2double(TaskIO.FieldIdx_R);
TaskIO.Correct = str2double(TaskIO.Correct);
TaskIO.RT = str2double(TaskIO.RT);
TaskIO.DateTime_Write = datetime(TaskIO.DateTime_Write);

%% Specify all PairTypes and their indices
AllIdx = (0:35)';
ZeropuIdx = [(0:5)';(6:6:30)'];
ComunsIdx = [3;6;10;23;26;31;34];
ComunsIdx = ComunsIdx(~ismember(ComunsIdx,ZeropuIdx));
NonunsIdx = [8,13,17,21,28,32]';
NonunsIdx = NonunsIdx(~ismember(NonunsIdx,ZeropuIdx));
SupervIdx = AllIdx(~ismember(AllIdx,[ZeropuIdx;ComunsIdx;NonunsIdx]));
TypeIdxs = struct();
TypeIdxs.Superv = SupervIdx;
TypeIdxs.Zeropu = ZeropuIdx;
TypeIdxs.Comuns = ComunsIdx;
TypeIdxs.Nonuns = NonunsIdx;
PairTypes = categorical(fieldnames(TypeIdxs));

%% Sort the data
TaskIO = sortrows(TaskIO,...
    {'DateTime_Write','AttemptNum'},{'ascend','ascend'});

%% Check for post-scan entries
[bool] = checkForPostScanEvents(pathToTaskIO,TaskIO.DateTime_Write);
if any(bool)
    if contains(SubjectId,...
            {'3d00293e';'b777a363';'b8282457';'d0d844d6';'eec99b44'})
        % This is the list of people who did extra training
        TaskIO = TaskIO(~bool,:);
    else
        error('You should probably be aware of this: %s',SubjectId);
    end
end

%% Construct trial-unique IDs for each trial
M = [TaskIO.SessionId,TaskIO.TrialId];
uM = unique(M,'rows');
M = mat2cell(M',2,ones(1,size(M,1)))';
uM = mat2cell(uM',2,ones(1,size(uM,1)))';
TaskIO.TrialId = cellfun(@(v)find(cellfun(@(u)all(v==u),uM)),M);

%% Loop through each trial to construct the output
Data = struct;
PairCounts = zeros(6);
for iTrial = 1:max(TaskIO.TrialId)
    T = TaskIO(TaskIO.TrialId==iTrial,:);
    Data.SubjectId(iTrial,1) = T.SubjectId(1);
    Data.SessionId(iTrial,1) = T.SessionId(1);
    Data.TrialId(iTrial,1) = iTrial - 1;
    Data.tSup(iTrial,1) = mean(PairCounts(SupervIdx+1));
    Data.PairId(iTrial,1) = T.PairId(1);
    Data.PairType(iTrial,1) = PairTypes(...
        structfun(@(v)ismember(T.PairId(1),v),TypeIdxs));
    Data.TrialType(iTrial,1) = T.TrialType(1);
    Data.FieldIdx_A(iTrial,1) = T.FieldIdx_A(1);
    Data.FieldIdx_B(iTrial,1) = T.FieldIdx_B(1);
    Data.FieldIdx_C(iTrial,1) = T.FieldIdx_C(1);
    Data.FieldIdx_R(iTrial,1) = {T.FieldIdx_R};
    Data.RT(iTrial,1) = {T.RT};
    % Update pair counts
    PairId = Data.PairId(iTrial);
    PairCounts(PairId+1) = PairCounts(PairId+1) + 1;
end

%% Redefine TaskIO and save the data locally
TaskIO = struct2table(Data);
save([pathToTaskIO,filesep,'TrainTaskIO.mat'],'TaskIO');

return

function [bool] = checkForPostScanEvents(pathToTaskIO,trainTimes)
fileList = dir([pathToTaskIO,filesep,'*_R*.mat']);
fn = fileList(end).name;
iS = strfind(fn,'_202')+1;
iE = iS + 14;
scanTime = datetime(fn(iS:iE));
bool = trainTimes > scanTime;
return