function [DataTable01] = getDataTable01(G)

if exist('DataTable01.mat','file')
    DataTable01 = load('DataTable01.mat');
    DataTable01 = DataTable01.DataTable01;
    return
end

pNonc = get_pNonc(G);
pNonc.cpNonc = pNonc.pNonc - mean(pNonc.pNonc);
[roiIds,roiNames] = getRoiList();
DataTable01 = [];
for iRoi = 1:numel(roiIds)
    cRoiId = roiIds{iRoi};
    cRoiName = roiNames{iRoi};
    PatternSim = getPatternSim(G,cRoiId);
    hemi = double(strcmp(cRoiName(1),'r'))*2 - 1;
    roiInfo = table(...
        repmat(hemi,size(PatternSim,1),1),...
        repmat(categorical({cRoiName}),size(PatternSim,1),1),...
        'VariableNames',{'hemisphere','roiName'});
    PatternSim = [roiInfo,PatternSim]; %#ok<AGROW>
    PatternSim = outerjoin(pNonc,PatternSim);
    % Rename subjectId after outerjoin
    PatternSim.subjectId_PatternSim = [];
    PatternSim.Properties.VariableNames{1} = 'subjectId';
    DataTable01 = [DataTable01;PatternSim]; %#ok<AGROW>
end

% Sort by subjectId
DataTable01 = sortrows(DataTable01,'subjectId');

save('DataTable01.mat','DataTable01');
return