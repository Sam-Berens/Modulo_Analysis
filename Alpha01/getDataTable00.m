function [DataTable00] = getDataTable00(G)
if exist('DataTable00.mat','file')
    DataTable00 = load('DataTable00.mat');
    DataTable00 = DataTable00.DataTable00;
    return
end

pNonc = get_pNonc(G);
pNonc.mCpNonc = pNonc.pNonc - mean(pNonc.pNonc);
[roiIds,roiNames] = getRoiList();
DataTable00 = [];
for iRoi = 1:numel(roiIds)
    cRoiId = roiIds{iRoi};
    cRoiName = roiNames{iRoi};
    PatternSim = getPatternSim(G,cRoiId);
    roiInfo = table(...
        repmat(categorical({cRoiName}),size(PatternSim,1),1),...
        'VariableNames',{'roiName'});
    PatternSim = [roiInfo,PatternSim]; %#ok<AGROW>
    [~,idxPS] = ismember(PatternSim.subjectId,pNonc.subjectId);
    PatternSim.pNonc = nan(height(PatternSim),1);
    PatternSim.pNonc = pNonc.pNonc(idxPS);
    PatternSim.mCpNonc = nan(height(PatternSim),1);
    PatternSim.mCpNonc = pNonc.mCpNonc(idxPS);
    DataTable00 = [DataTable00;PatternSim]; %#ok<AGROW>
end

% Make Hemi var
hemi = double(contains(roiNames,'l'));
hemi(~hemi) = -1;
hemiInfo = table(categorical(roiNames),hemi,'VariableNames',{'roiName','hemisphere'});
DataTable00 = join(DataTable00,hemiInfo);
return
