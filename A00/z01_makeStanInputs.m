function [] = z01_makeStanInputs()
% z01_makeStanInputs.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/07/2025
%
% Syntax:  z01_makeStanInputs()
%
% Description:
%    Loads TrainTaskIO for each subject and extract key variables needed
%    to estimate a model in Stan. These variables are recorded in a set of
%    local JSON files.
%
% Inputs:
%    NONE
%
% Outputs:
%    NONE
%
% Example:
%    z01_makeStanInputs();
%
%% Make paths structure
dirPaths = struct;
dirPaths.Data = ['..',filesep,'..',filesep,'Data',filesep];

%% Get a list of subject IDs
dirLits = dir(dirPaths.Data);
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);

%%  Loop through each subject
for iSubject = 1:numel(subjectIds)
    sId = char(subjectIds(iSubject));
    dirPaths.Target = [dirPaths.Data,sId,filesep,...
        'Behavioural',filesep,'A00',filesep,'InputData',filesep];
    if ~exist(dirPaths.Target,'dir')
        mkdir(dirPaths.Target);
    end
    TaskIO = getTrainTaskIO(sId);
    for pairId = 0:(6^2-1)
        % Get the Json string
        DataStruct = getPairStruct(TaskIO,pairId);
        Json = jsonencode(DataStruct,...
            'PrettyPrint',true,'ConvertInfAndNaN',false);
        %% Save the Json string
        fid = fopen(sprintf('%sP%02d.json',dirPaths.Target,pairId),'w');
        fprintf(fid,'%s',Json);
        fclose(fid);
    end
end
return

function [DataStruct] = getPairStruct(TaskIO,pairId)

% Select data from TaskIO related to the queried pair
PairIO = TaskIO(TaskIO.PairId==pairId,:);

%% Construct variables
n = size(PairIO,1);
x = PairIO.tSup;
maxX = max(x);
c = nan(numel(x),1); % The targets
Y = nan(numel(x),6); % The actual responses
for ii = 1:size(PairIO,1)
    c(ii) = PairIO.FieldIdx_C(ii);
    r = PairIO.FieldIdx_R{ii};
    r = unique(r,'stable');
    if numel(r) < 6
        r = [r;nan(6-numel(r),1)]; %#ok<AGROW>
    end
    Y(ii,:) = r';
end

%% Package
DataStruct = struct();
DataStruct.n = n;
DataStruct.x = x;
DataStruct.maxX = maxX;
DataStruct.c = c;
DataStruct.Y = Y;
return