function [fgh] = plotMdl5AtMNI(dt5,coordList)

formula = 'zTemplate ~ colocation * cpNonc + (1|subjectId)';
mdl05 = fitRoiLmes(coordList,formula,dt5);
names = unique(dt5.roiName);
fgh = nan(numel(names),1);
for ii=1:numel(names)
    cRoi = names{ii};
    cMdl = mdl05.mdl.mdl{ii};
    cDT = dt5(dt5.roiName == cRoi,:);
   fgh(ii) = plotRoiEffect(cDT,cMdl,cRoi);
end
return