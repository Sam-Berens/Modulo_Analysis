function [DataTable] = z03_CorrEuc(G)

subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
H = (3-min(cat(3,mod((0:5)'-(0:5),6),mod((0:5)-(0:5)',6)),[],3))./3;
lower = tril(true(6),-1);

for iSubject=1:numel(subjectIds)
    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];
    %% Get roi paths 
    [roiPaths] = getROIpaths(dirs, G);
    for iRoi = 1:numel(roiPaths)
    %% Sample t-stat data using mask ROI coordinates
    data.M = sampleTs(dirs, roiPaths{iRoi});
    %reminder that this M is (nStims,nVoxels)
    % Compute pairwise Euclidian distances across conditions
    D = pdist(data.M,'euclidean');
    D = squareform(D);
    r = corr(H(lower),D(lower),'Type','Kendall');
    z = atanh(r);
    DataTable(iSubject,1).(roiName) = z;
    end

    DataTable = struct2table(DataTable,'RowNames',subjectIds);
    varNames = DataTable.Properties.VariableNames;

    for ii = 1:numel(varNames)
        v = DataTable.(varNames{ii});
        s = cellfun(@isempty,v);
        v(s) = {NaN};
        v = cell2mat(v);
        DataTable.(varNames{ii}) = v;
    end

end

[~,lE.p,lE.ci,lE.stats] = ttest(DataTable.rw_Neuromorph_lEntorhinal);
[~,rE.p,rE.ci,rE.stats] = ttest(DataTable.rw_Neuromorph_rEntorhinal);
[~,lH.p,lH.ci,lH.stats] = ttest(DataTable.rw_RitcheyProb_lHippComb);
[~,rH.p,rH.ci,rH.stats] = ttest(DataTable.rw_RitcheyProb_rHippComb);
[~,lV.p,lV.ci,lV.stats] = ttest(DataTable.rw_Schaefer2018_N17P200_R001);
[~,rV.p,rV.ci,rV.stats] = ttest(DataTable.rw_Schaefer2018_N17P200_R101);
[~,lM.p,lM.ci,lM.stats] = ttest(DataTable.rw_Schaefer2018_N17P200_R082);
[~,rM.p,rM.ci,rM.stats] = ttest(DataTable.rw_Schaefer2018_N17P200_R187);

return

function [roiList] = getROIpaths(dirs,G)
dirs.NtvRois = [dirs.Subject,filesep,'Structural',filesep,G,...
    filesep,'NativeRois',filesep,'w'];
roiList = dir([dirs.NtvRois,filesep,'*.nii']);
roiList = arrayfun(@(x) [x.folder,filesep,x.name],roiList,...
    'UniformOutput',false);
return