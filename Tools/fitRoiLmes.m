function [Mdl] = fitRoiLmes(rois,formula,T)

%% Check the input table
if ~any(ismember(T.Properties.VariableNames,'roiName'))
    error('The input table T must contain a variable named "roiName".');
end

%% Prep the output
Mdl = struct();
Mdl.formula = formula;
Mdl.rois = rois;

%% Loop through to fit the models
fits = cell(size(rois));
for iRoi = 1:numel(rois)

    s = contains(cellstr(T.roiName),rois{iRoi});
    cT = T(s,:);
    m = fitlme(cT,formula);

    if iRoi == 1
        nFfx = m.NumCoefficients;
        sFfx = m.CoefficientNames;

        estm = nan(numel(rois),nFfx);
        tVal = nan(numel(rois),nFfx);
        pVal = nan(numel(rois),nFfx);
        dfe  = nan(numel(rois),nFfx);
        ciLo = nan(numel(rois),nFfx);
        ciHi = nan(numel(rois),nFfx);

        ffxRawsz = strings(numel(rois), nFfx);
        ffxTests = strings(numel(rois), nFfx);
    end

    % Core stats
    estm(iRoi,:) = m.Coefficients.Estimate';
    tVal(iRoi,:) = m.Coefficients.tStat';
    pVal(iRoi,:) = m.Coefficients.pValue';
    dfe(iRoi,:) = repmat(m.DFE, 1, nFfx);
    ciLo(iRoi,:) = m.Coefficients.Lower';
    ciHi(iRoi,:) = m.Coefficients.Upper';

    % Build printable strings for each coefficient
    for jj = 1:nFfx
        ffxRawsz(iRoi,jj) = sprintf(...
            '%s [%s, %s];', ...
            fmtNum(estm(iRoi,jj)), ...
            fmtNum(ciLo(iRoi,jj)), ...
            fmtNum(ciHi(iRoi,jj)));
        ffxTests(iRoi,jj) = sprintf(...
            't(%s) = %s, p = %s;', ...
            fmtInt(dfe(iRoi,jj)), ...
            fmtNum(tVal(iRoi,jj)), ...
            fmtP(pVal(iRoi,jj)) );
    end

    fits{iRoi} = m;
end

%% Set the display tables
Mdl.mdl = table(fits,'RowNames',rois,'VariableNames',{'mdl'});

Mdl.FFX.estm = array2table(estm,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.tVal = array2table(tVal,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.pVal = array2table(pVal,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.ciLo = array2table(ciLo,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.ciHi = array2table(ciHi,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.dfe = array2table(dfe, 'RowNames',rois,'VariableNames',sFfx);

Mdl.FFX.rawsz = array2table(ffxRawsz,'RowNames',rois,'VariableNames',sFfx);
Mdl.FFX.tests = array2table(ffxTests,'RowNames',rois,'VariableNames',sFfx);

return

function s = fmtNum(x)
% Compact numeric formatting for estimates, t, CI, etc.
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.3g', x);
end
return

function s = fmtInt(x)
% DFE is typically integer-like; handle NaN gracefully
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%d', round(x));
end
return

function s = fmtP(p)
% p-value formatting with a common convention
if isnan(p)
    s = 'NaN';
elseif p < 1e-4
    s = '<1e-4';
else
    s = sprintf('%.3g', p);
end
return