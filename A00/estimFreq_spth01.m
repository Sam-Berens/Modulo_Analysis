function [spth01,Data,fh] = estimFreq_spth01(SubjectId,dispPlot)
% estimFreq_spth01.m
% Sam Berens (s.berens@sussex.ac.uk)
% 24/03/2025
%
% Syntax:  [spth01, Data, fh] = estimFreq_spth01(SubjectId, dispPlot)
%
% Description:
%    Estimates model parameters for a Softplus-tanh model using the von
%    Mises distribution for a given subject. The function retrieves subject
%    data using getSubjectData, processes the data, estimates the model
%    parameters for different pair types via negative log-likelihood
%    minimization (using fmincon with MultiStart), predicts probability
%    mass functions (PMFs) for discrete response angles, and optionally
%    plots the prediction curves.
%
% Inputs:
%    SubjectId - (optional) Character vector specifying the subject
%                           identifier.
%                Default is 'eec99b44'.
%    dispPlot  - (optional) Boolean flag indicating whether to display
%                           prediction plots.
%                Default is true.
%
% Outputs:
%    spth01 - A table containing the estimated model parameters (b0, b1) 
%             and negative log-likelihood (nll) for each pair type.
%    Data   - A table containing the subject data along with predicted
%             outcome fields.
%    fh     - If dispPlot is true, fh is the figure handle of the plotted
%             prediction curves. If dispPlot is false, fh returns as NaN.
%
% Example:
%    [spth01, Data, fh] = estimFreq_spth01('subject123', true);
%
% See also: getSubjectData, spth_nll, spth_pred
%

%% Get the data
if nargin < 1
    SubjectId = 'eec99b44';
    dispPlot = true;
elseif nargin < 2
    dispPlot = true;
elseif ~dispPlot
    fh = nan;
end
Data = getSubjectData(SubjectId);
Data.AttemptIds = [];

%% Add empty variables to Data that will store outcomes and predicted PMFs
Data.Y = nan(size(Data,1),6);
for iTheta = 0:3
    Data.(['pmf',num2str(iTheta)]) = nan(size(Data,1),1);
end

%% Specify PairTypes
PairType = categorical({'Superv';'Zeropu';'Comuns';'Nonuns'});

%% Loop through each trial type to estimate the model
if dispPlot
    fh = figure;
end
bHat = cell(numel(PairType),1);
nll = nan(numel(PairType),1);
for iPairType = 1:numel(PairType)
    
    % Extract rows of Data that correspond to the current pair type
    rowSel = Data.PairType==PairType(iPairType);
    sData = Data(rowSel,:);
    
    %% Extract the predictor and response variables
    x = sData.tSup;
    maxX = max(x);
    Y = nan(numel(x),6);
    for ii = 1:size(sData,1)
        r = sData.FieldIdx_R{ii};
        r = unique(r,'stable');
        if numel(r) < 6
            r = [r;nan(6-numel(r),1)]; %#ok<AGROW>
        end
        Y(ii,:) = wrapTo2Pi((r'-sData.FieldIdx_C(ii)).*(pi/3));
    end
    Data.Y(rowSel,:) = Y;
    
    %% Estimate the model
    A = [-1,0;1,0;0,1];
    b = [0;maxX+1;3];
    problem = createOptimProblem(...
        'fmincon',...
        'x0',[maxX/2;0.01],...
        'objective',@(p) spth_nll(p,x,Y),...
        'Aineq',A,'bineq',b);
    [bHat{iPairType},nll(iPairType)] = run(MultiStart,problem,50);
    
    %% Predict
    [pmf,angles] = spth_pred(x,bHat{iPairType});
    for iTheta = 0:3
        sTheta = abs(abs(wrapToPi(angles))-(pi*iTheta/3))<(1e-6);
        Data.(['pmf',num2str(iTheta)])(rowSel,:) = sum(pmf(:,sTheta),2);
    end
    
    %% Plot
    if dispPlot
        plotPreds(x,pmf,angles);
    end
end

%% Save the outputs
b0 = cellfun(@(v)v(1),bHat);
b1 = cellfun(@(v)v(2),bHat);
spth01 = table(PairType,b0,b1,nll);
return

function [] = plotPreds(x,pmf,angles)
% plotPreds plots the predicted PMFs for each discrete angle bin.
%
% Inputs:
%    x      - Vector of predictor values.
%    pmf    - Matrix of predicted probability mass functions.
%    angles - Vector of angle values (in radians).
%
% This function creates a 2x2 subplot, one for each angle bin (0 to
% 3*pi/3), and plots the predicted probability curves.

persistent lineCount;
if isempty(lineCount)
    lineCount = 1;
else
    lineCount = lineCount + 1;
end

for iTheta = 0:3
    subplot(2,2,iTheta+1);
    hold on;
    sTheta = abs(abs(wrapToPi(angles))-(pi*iTheta/3))<(1e-6);
    p = sum(pmf(:,sTheta),2);
    plot(x,p);
    title(sprintf('%s=%f',char(952),mean(abs(wrapToPi(angles(sTheta))))));
    if (iTheta==0) && (lineCount==4)
        legend({'Superv','Zeropu','Comuns','Nonuns'},'Location','best');
    end
end
return