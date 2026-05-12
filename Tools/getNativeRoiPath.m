function [roiPath] = getNativeRoiPath(G,subjectId,roiStrPattern)
% Helper function to return the full file path for an ROI
dirs.Subject = [...
    fileparts(mfilename('fullpath')),filesep,...
    '..',filesep,'..',filesep,...
    'Data',filesep,subjectId];

dirs.Roi = [...
    dirs.Subject,filesep,...
    'Structural',filesep,...
    G,filesep,...
    'NativeRois',filesep,...
    'w'];

roiPath = dir([dirs.Roi,filesep,'*',roiStrPattern,'*']);
if ~numel(roiPath)
    roiPath = dir([dirs.Roi,filesep,roiStrPattern,'*']);
end

if numel(roiPath) == 0
    error('Bad ROI pattern.');
end

if numel(roiPath) > 1
    error('Ambiguous ROI pattern.');
end

roiPath = [roiPath.folder, filesep, roiPath.name];
return