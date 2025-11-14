function [DataTable] = z03_Alpha00_0(G)

subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];

%% Loop through rois
roiNames=['lHip','rHip','lEnt','rEnt','lMPFC','lVis','rVis'];
roiIds = ['lHip','rHip','rEnt','082','187','001','101'];

%check this is the way im supposed to stack them
patternSim = nan(6,6,numel(subjectIds),numel(roiIds));
for ii=1:numel(roiIds)

    [zTemplate,patternSim(:,:,:,ii)] = getPatternSim(G,roiIds(ii));

    %meant to arrange stuff such that you have a model object which is a struct
    %with rows for roi and fields incld name of roi, estimated result of
    %running model, and so on

    %% TO DO estim glme/lme
    %input data is zTemplate
    %need to compute the profiency stuff before being fully ready

end

return

