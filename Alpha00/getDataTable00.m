function [DataTable00] = getDataTable00(G)

if exist('DataTable00.mat','file')
    DataTable00 = load('DataTable00.mat');
    DataTable00 = DataTable00.DataTable00;
    return
end

pNonc = get_pNonc(G);
pNonc.cpNonc = pNonc.pNonc - mean(pNonc.pNonc);
[roiIds,roiNames] = getRoiList();
DataTable00 = [];
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
    PatternSim = join(pNonc,PatternSim);
    DataTable00 = [DataTable00;PatternSim]; %#ok<AGROW>
end
DataTable00 = sortrows(DataTable00,'subjectId');

save('DataTable00.mat','DataTable00');
return