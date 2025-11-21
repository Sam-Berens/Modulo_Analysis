function [mdlObj,patternSim] = z03_Alpha00_0(G)

subjectIds = getSubjectIds(G);
subjectIds = categorical(subjectIds);
%naming convention is that captilised variables are table columns!!!

nSubs = numel(subjectIds);
dirs.Data = ['..',filesep,'..',filesep,'Data'];

%% Loop through rois
roiNames={'lHip','rHip','lEnt','rEnt','lMPFC','rMPFC','lVis','rVis'};
%fyi odd indexes are left hemisphere, even indexes are right hemisphere
roiIds = ["lHi";"rHi";"lEn";"rEn";"082";"187";"001";"101"];

%check this is the way im supposed to stack them
patternSim = nan(6,6,nSubs,numel(roiIds));
nHems = 2;
SubjectIds = repmat(subjectIds,[nHems,1]);
nRows = nSubs * nHems;
%collect predictors 
p_noncom = getP_noncom(G);

%this is going to create one model per roi (L&R comb)
for ii=1:numel(roiIds)
    cRoiId = roiIds(ii,1);
    cRoiNm = roiNames{1,ii};
    [zTemplate,patternSim(:,:,:,ii)] = getPatternSim(G,cRoiId);

    %% TO DO estim glme/lme
    %input data is zTemplate
    %isOdd = mod(ii,2);
    isOdd = contains(cRoiNm,'l');
    if isOdd
        Ztemplate = nan(nRows,1);
        Hemisphere = nan(nRows,1);
        P_noncom = nan(nRows,1);
        hemisphereCode = ones(1,nSubs);
        cRows = (1:nSubs)';
    else
        hemisphereCode = -1 *ones(1,nSubs);
        cRows = ((nSubs+1):(nSubs*nHems))';
    end
    
    Ztemplate(cRows,1) = zTemplate;
    Hemisphere(cRows,1) = hemisphereCode;
    P_noncom(cRows,1) = p_noncom;
    
    %we collect the columns for the left and then right hemisphere 
    %and only construct the table on the second hemisphere so theres
    %one table per comb roi
    if ~isOdd
      T = table(Ztemplate,Hemisphere,P_noncom,SubjectIds);
      mdlObj.(cRoiNm(2:end)) = fitlme(T, 'Ztemplate ~ Hemisphere + P_noncom + (1|SubjectIds)');
    end 
    %no Colocation in this model because its not Alpha01!

end

return

