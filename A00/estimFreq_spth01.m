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
%    parameters for different pairs via negative log-likelihood
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
%             and negative log-likelihood (nll) for each pair.
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

%% Add empty variables to Data that will store outcomes and predicted PMFs
Data.Y = nan(size(Data,1),6);
for iTheta = 0:3
    Data.(['pmf',num2str(iTheta)]) = nan(size(Data,1),1);
end

%% Specify PairId and preallocate PairType
PairId = unique(Data.PairId);
PairType = cell(size(PairId));

%% Loop through each trial type to estimate the model
bHat = cell(numel(PairId),1);
nll = nan(numel(PairId),1);
for iPairId = 1:numel(PairId)

    % Extract rows of Data that correspond to the current pair ID
    rowSel = Data.PairId==PairId(iPairId);
    sData = Data(rowSel,:);
    PairType{iPairId} = char(sData.PairType(1));

    %% Extract the predictor and response variables
    x = sData.tSup;
    maxX = max(x);
    c = nan(numel(x),1); % The targets
    Y = nan(numel(x),6); % The actual responses
    for ii = 1:size(sData,1)
        c(ii) = sData.FieldIdx_C(ii);
        r = sData.FieldIdx_R{ii};
        r = unique(r,'stable');
        if numel(r) < 6
            r = [r;nan(6-numel(r),1)]; %#ok<AGROW>
        end
        Y(ii,:) = r';
    end
    Data.Y(rowSel,:) = Y;

    %% Estimate the model
    A = [-1,0;1,0;0,1];
    b = [0;maxX+1;13];
    problem = createOptimProblem(...
        'fmincon',...
        'x0',[maxX/2;0.01],...
        'objective',@(p) spth_nll(p,x,c,Y),...
        'Aineq',A,'bineq',b);
    [bHat{iPairId},nll(iPairId)] = run(MultiStart,problem,50);

    %% Predict
    [pmf,angles] = spth_pred(x,bHat{iPairId});
    for iTheta = 0:3
        sTheta = abs(abs(wrapToPi(angles))-(pi*iTheta/3))<(1e-6);
        Data.(['pmf',num2str(iTheta)])(rowSel,:) = sum(pmf(:,sTheta),2);
    end
end

%% Save the outputs
b0 = cellfun(@(v)v(1),bHat);
b1 = cellfun(@(v)v(2),bHat);
PairType = categorical(PairType);
spth01 = table(PairId,PairType,b0,b1,nll);

%% Plot
if dispPlot
    fh = figure;
    hold on;
    for iPairId = 1:numel(PairId)
        rowSel = Data.PairId==PairId(iPairId);
        sData = Data(rowSel,:);
        x = sData.tSup;
        y = sData.pmf0;
        if sData.PairType(1)==categorical({'Zeropu'})
            p1 = plot(x,y,'Color','black','LineStyle','-', ...
                'DisplayName','Zeropu');
        elseif sData.PairType(1)==categorical({'Superv'})
            p2 = plot(x,y,'Color','red','LineStyle','-', ...
                'DisplayName','Superv');
        elseif sData.PairType(1)==categorical({'Comuns'})
            p3 = plot(x,y,'Color','green','LineStyle','-', ...
                'DisplayName','Comuns');
        else
            p4 = plot(x,y,'Color','blue','LineStyle','-', ...
                'DisplayName','Nonuns');
        end
    end
    legend([p1,p2,p3,p4]);
end
return