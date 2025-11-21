function [mdlObj] = z03_Alpha01_0(G)

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
patternSim = nan(12,12,nSubs,numel(roiIds));
nHems = 2;
nPos = 2;
SubjectIds = repmat(subjectIds,[(nHems*nPos),1]);
%nRows = nSubs * nHems * nPos;
%collect predictors
p_noncom = getP_noncom(G);

%this is going to create one model per roi (L&R comb)
for ii=1:numel(roiIds)
    cRoiId = roiIds(ii,1);
    cRoiNm = roiNames{1,ii};
    %zTemplate has nSubjects going down, and colocation variable going
    %across
    [zTemplate,patternSim(:,:,:,ii)] = getPatternSim(G,cRoiId);
    isLeft = contains(cRoiNm,'l');

    for jj=1:2
        isSamePos = mod(jj,2);
        if isSamePos
            coLocationCode = ones(nSubs,1);
            cZtemplate = zTemplate(:,1);
        else
            coLocationCode = -1 * ones(nSubs,1);
            cZtemplate = zTemplate(:,2);
        end

        if isLeft
            hemisphereCode = ones(nSubs,1);
            if  isSamePos
                %make empty arrays on the first go only
                Ztemplate = cZtemplate;
                Hemisphere = hemisphereCode;
                P_noncom = p_noncom;
                CoLocation = coLocationCode;
            else
                Ztemplate = [Ztemplate;cZtemplate]; %#ok<AGROW>
                Hemisphere = [Hemisphere;hemisphereCode]; %#ok<AGROW>
                P_noncom = [P_noncom;p_noncom]; %#ok<AGROW>
                CoLocation = [CoLocation;coLocationCode];%#ok<AGROW>
            end
        else
            hemisphereCode = -1 *ones(nSubs,1);
            Ztemplate = [Ztemplate;cZtemplate]; %#ok<AGROW>
            Hemisphere = [Hemisphere;hemisphereCode]; %#ok<AGROW>
            P_noncom = [P_noncom;p_noncom]; %#ok<AGROW>
            CoLocation = [CoLocation;coLocationCode];%#ok<AGROW>
            if ~isSamePos
                T = table(Ztemplate,CoLocation,Hemisphere,P_noncom,...
                    SubjectIds);
                mdlObj.(cRoiNm(2:end)) = fitlme(T,...
                    'Ztemplate ~ CoLocation + Hemisphere + P_noncom + (1|SubjectIds)');
            end
        end

    end
end

return
