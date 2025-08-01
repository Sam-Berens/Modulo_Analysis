function [] = z01_makeInputs()
% z01_makeInputs
% Sam Berens (s.berens@sussex.ac.uk)
% 24/07/2025
%
% Syntax: z01_makeInputs()
%
% Description:
%    Loads TrainTaskIO for each subject and extract key variables needed
%    to estimate a model in Stan. These variables are recorded in a local 
%    JSON file.
%
% Inputs:
%    NONE
%
% Outputs:
%    NONE
%
% Example:
%    z01_makeInputs();
%

%% Make paths structure
dirStruct = struct;
dirStruct.Data = ['..',filesep,'..',filesep,'..',filesep,'Data',filesep];

%% Get a list of subject IDs
dirLits = dir(dirStruct.Data);
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);

%%  Loop through each subject
for iSubject = 1:numel(subjectIds)
    sId = char(subjectIds(iSubject));
    dirStruct.Target = [dirStruct.Data,sId,filesep,...
        'Analysis',filesep,'A00',filesep];
    if ~exist(dirStruct.Target,'dir')
        mkdir(dirStruct.Target);
    end
    TaskIO = getTrainTaskIO(sId);
    % Get the Json string
    DataStruct = getDataStruct(TaskIO);
    Json = jsonencode(DataStruct,...
        'PrettyPrint',true,'ConvertInfAndNaN',false);
    %% Save the Json string
    fid = fopen([dirStruct.Target,'InputData.json'],'w');
    fprintf(fid,'%s',Json);
    fclose(fid);
end
return

function [DataStruct] = getDataStruct(TaskIO)

% nTrials: Trial count, Int[1,Inf]
nTrials = size(TaskIO,1);

% pairId: Identity of the pair tested on each trial, Int[1,36]
pairId = TaskIO.PairId + 1;
nPairs = numel(unique(pairId));

% pairType: Type of pair tested on each trial, Int[1,4]
pairType = dummyvar(TaskIO.PairType)*[3;4;2;1];
nTypes = numel(unique(pairType));
id2type = nan(nPairs,1);
for ii = 1:numel(pairId)
    id2type(pairId(ii)) = pairType(ii);
end

% x: Learning state predictor, Real[0,Inf]
x = TaskIO.tSup;
maxX = max(x);

% c: Target response on each trial, Int[0,5]
c = nan(numel(x),1);

% Y: Response per trial and attempt, Array of Ints[0,5]
Y = nan(numel(x),6); % The actual responses
for ii = 1:size(TaskIO,1)
    c(ii) = TaskIO.FieldIdx_C(ii);
    r = TaskIO.FieldIdx_R{ii};
    r = unique(r,'stable');
    if numel(r) < 6
        r = [r;nan(6-numel(r),1)]; %#ok<AGROW>
    end
    Y(ii,:) = r';
end

%% Package
DataStruct = struct();
DataStruct.nTrials = nTrials;
DataStruct.nPairs = nPairs;
DataStruct.pairId = pairId;
DataStruct.nTypes = nTypes;
DataStruct.id2type = id2type;
DataStruct.maxX = maxX;
DataStruct.x = x;
DataStruct.c = c;
DataStruct.Y = Y;
return