function [roiPath] = getNativeRoiPath(G,subjectId,roiPattern)
% Helper function to return the full file path for an ROI

dirs.Subject = [...
    '..',filesep,'..',filesep,...
    'Data',filesep,subjectId];

dirs.Roi = [...
    dirs.Subject,filesep,...
    'Structural',filesep,...
    G,filesep,...
    'NativeRois',filesep,...
    'w'];

roiPath = dir([dirs.Roi,filesep,'*',roiPattern,'*']);
if ~numel(roiPath)
    roiPath = dir([dirs.Roi,filesep,roiPattern,'*']);
end

if numel(roiPath) == 0
    error('Bad ROI pattern.');
end

if numel(roiPath) > 1
    error('Ambiguous ROI pattern.');
end

roiPath = [roiPath.folder, filesep, roiPath.name];

return