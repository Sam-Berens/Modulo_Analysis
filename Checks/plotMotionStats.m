function [] = plotMotionStats()

censorThresh = 0.9;

% Set dataDir and subjList
dataList = dir(['..',filesep,'..',filesep,'Data']);
dataList = dataList(cellfun(@(s)numel(s)==8,{dataList.name}));
dataDir = dataList.folder;
subjList = {dataList.name}';

rpsFolder = [dataDir,filesep,'RPs'];
if ~exist(rpsFolder,'dir')
    mkdir(rpsFolder);
end

iIn = 0;
for iSubject = 1:numel(subjList)

    % Set subjectId, rpList, nRuns
    subjectId = subjList{iSubject};
    subjDir = [dataDir,filesep,subjectId];
    epiDir = [subjDir,filesep,'EPI'];
    rpList = dir([epiDir,filesep,'*.txt']);
    nRuns = numel(rpList);

    % Run loop
    for iR = 1:nRuns
        iIn = iIn + 1;

        % Set RPs and nEpis
        FileName_RpText = [rpList(iR).folder,filesep,rpList(iR).name];
        RPs = readmatrix(FileName_RpText);
        nEpis = size(RPs,1);

        % Compute derivatives
        dRPs = diff(RPs,1,1);
        qt = (0:1:(nEpis-1))';
        t = mean([qt(1:end-1),qt(2:end)],2);
        dRPs = cell2mat(...
            cellfun(@(y)sincInterp(qt,t,y),...
            mat2cell(dRPs,nEpis-1,ones(1,6)),...
            'UniformOutput',false));

        % Compute frame-wise displacement (FWD)
        translations = sum(abs(dRPs(:,1:3)),2);
        rotations = sum(abs(dRPs(:,4:6)),2);
        fwd = translations + 50*rotations;

        % Compute censors
        volIdxs = find(fwd > censorThresh);
        censors = zeros(nEpis,numel(volIdxs));
        for ii = 1:numel(volIdxs)
            censors(volIdxs(ii),ii) = 1;
        end

        RPs(:,4:end) = rad2deg(RPs(:,4:end));
        dRPs(:,4:end) = rad2deg(dRPs(:,4:end));

        %% Plots
        % Translations
        fh1 = figure('Units','normalized','OuterPosition',[0,0,1,1]);
        subplot(2,1,1);
        hold on;
        plot(RPs(:,1),'-r');
        plot(RPs(:,2),'-g');
        plot(RPs(:,3),'-b');
        plot(dRPs(:,1),'--','Color',[1.0,0.5,0.5]);
        plot(dRPs(:,2),'--','Color',[0.5,1.0,0.5]);
        plot(dRPs(:,3),'--','Color',[0.5,0.5,1.0]);
        Legend = legend('X','Y','Z','\DeltaX','\DeltaY','\DeltaZ');
        set(Legend,'Location','eastoutside');
        title(sprintf('Translations: %s iRun=%i',subjectId,iR));
        xlabel('Scans \rightarrow');
        ylabel('Displacement /mm');
        ax = gca;
        ax.YAxis.Exponent = 0;   % turn off the ×10^n scaling
        ax.YRuler.TickLabelFormat = '%.2f';  % control decimal places
        ax.YTickLabelRotation = 45;
        ylim([-1.5,1.5]);
        hold off;

        % Rotations
        subplot(2,1,2);
        hold on;
        plot(RPs(:,4),'-m');
        plot(RPs(:,5),'-y');
        plot(RPs(:,6),'-c');
        plot(dRPs(:,4),'--','Color',[1.0,0.5,1.0]);
        plot(dRPs(:,5),'--','Color',[1.0,1.0,0.5]);
        plot(dRPs(:,6),'--','Color',[0.5,1.0,1.0]);
        Legend = legend('Pitch','Roll','Yaw',...
            '\DeltaPitch','\DeltaRoll','\DeltaYaw');
        set(Legend,'Location','eastoutside');
        title(sprintf('Rotations: %s iRun=%i',subjectId,iR));
        xlabel('Scans \rightarrow');
        ylabel('Rotation /rad');
        ax = gca;
        ax.YAxis.Exponent = 0;   % turn off the ×10^n scaling
        ax.YRuler.TickLabelFormat = '%.2f';  % control decimal places
        ax.YTickLabelRotation = 45;
        ylim([-2,2]);
        hold off;

        print(fh1,'-dpng',sprintf('%s%s%s-R%i_RPs.png',...
            rpsFolder,filesep,subjectId,iR));
        close(fh1);

        % create FWD graph
        fh2 = figure('Units','normalized','OuterPosition',[0,0,1,1]);
        hold on;
        plot(fwd,'--','Color',[0,0.5,1]);
        scatter(volIdxs,zeros(size(volIdxs)),'black','filled');
        ylim([0,1]);
        title(sprintf('Framewise displacement: %s iRun=%i',subjectId,iR));
        xlabel('Scans \rightarrow');
        ylabel('Displacement /mm');
        print(fh2,'-dpng',sprintf('%s%s%s-R%i_FWD.png',...
            rpsFolder,filesep,subjectId,iR));
        close(fh2);

    end
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return