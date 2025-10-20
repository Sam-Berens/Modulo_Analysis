function [] = z00_makeTrainTaskIO()
% z00_makeTrainTaskIO.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/07/2025
%
% Syntax:  z00_makeTrainTaskIO()
%
% Description:
%    Ensures that subject-specific training data (in TrainTaskIO) have been
%    retrieved from a web server, pre-processed and saved to a local file.
%
% Inputs:
%    NONE
%
% Outputs:
%    NONE
%
% Example:
%    z00_makeTrainTaskIO();
%
dirLits = dir(['..',filesep,'..',filesep,'Data']);
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);
fh = waitbar(0,'Saving data...');
for iSubject = 1:numel(subjectIds)
    sId = char(subjectIds(iSubject));
    [~] = getTrainTaskIO(sId);
    waitbar(iSubject/numel(subjectIds),fh);
end
close(fh);
return