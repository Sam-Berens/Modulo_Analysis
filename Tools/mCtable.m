function [mC] = mCtable(mdl,roiIdx)

fields = fieldnames(mdl);

% correcting across the first 4 rois, not V1!
if any(contains(fields,'rois'))
    nCols = size(mdl.FFX.pVal,2);
    newT = mdl.FFX.pVal(roiIdx,:);
    rjctNll = zeros(numel(roiIdx),nCols);
    for ii=1:nCols
        pVals = newT{roiIdx,ii};
        [~,cRjctNll] = holmBonf(0.05,pVals);
        rjctNll(:,ii) = cRjctNll;
        toKeep = any(rjctNll,2);
    end
    mC = mdl.FFX.pVal(roiIdx(toKeep),:);
else
    pVals = dataset2table(mdl.Coefficients(:,'pValue'));
    coefficient = dataset2table(mdl.Coefficients(:,'Name'));
    [~,rejectNll] = holmBonf(0.05,pVals{:,1});
    mC = table(coefficient,pVals);
    mC.Properties.VariableNames = {'coefficient','pVal'};
    mC = mC(rejectNll,:);
end


return