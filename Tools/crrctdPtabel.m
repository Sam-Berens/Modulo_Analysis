function [mC] = crrctdPtabel(mdl)

fields = fieldnames(mdl);

% correcting across the first 4 rois, not V1!
if any(contains(fields,'rois'))
    nCols = size(mdl.FFX.pVal,2);
    newT = mdl.FFX.pVal(1:4,:);
    for ii=1:nCols
        pVals = newT{1:4,ii};
        alphas = holmBonf(0.05,pVals);
        newT{1:4,ii} = alphas;
        varName = newT.Properties.VariableNames{ii};
        newT = renamevars(newT,varName,[varName,'_alpha']);
    end
    mC = join(mdl.FFX.pVal(1:4,:),newT,'Keys', 'RowNames');
    order = nan(1,16);
    for ii=1:8
        col = ii+(ii-1);
        order(col) = ii;
        order(col+1) = 8+ii;
    end
    mC = mC(:,order);
else
    pVals = dataset2table(mdl.Coefficients(:,'pValue'));
    coefficient = dataset2table(mdl.Coefficients(:,'Name'));
    alphas = holmBonf(0.05,pVals{:,1});
    alphas = table(alphas);
    mC = [coefficient,pVals,alphas];
    mC.Properties.VariableNames = {'coefficient','pVal','corrected_Alpha_Level'};
end


return