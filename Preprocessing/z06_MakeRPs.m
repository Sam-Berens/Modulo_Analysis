function [] = z06_MakeRPs(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir.folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];

rpsFolder = [dataDir,filesep,'RPs'];
if ~exist(rpsFolder,'dir')
    mkdir(rpsFolder);
end

censorThresh = 0.9; % mm FWD

rpList = dir([epiDir,filesep,'*.txt']);
nRuns = numel(rpList);

for iR = 1:nRuns
    % Get data
    FileName_RpText = [rpList(iR).folder,filesep,rpList(iR).name];
    RPs = readmatrix(FileName_RpText);
    nEpis = size(RPs,1);

    % Compute derivatives
    dRPs = diff(RPs,1,1);
    qt = (0:1:(nEpis-1))';
    t = mean([qt(1:end-1),qt(2:end)],2);
    dRPs = cell2mat(...
        cellfun(@(y)sincInterp(qt,t,y),mat2cell(dRPs,nEpis-1,ones(1,6)),...
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

    % Construct R
    R = [RPs,dRPs,RPs.^2,dRPs.^2,fwd]; %%%%%%%%%%%%%%%%%%%%%%%% TO CONSIDER
    R = zscore(R,[],1);
    R = [R,censors];

    % Save R
    save(sprintf('%s%sRP%i.mat',epiDir,filesep,iR),'R'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Print plots. THIS NEEDS SOME WORK TOO
    % translations
    subplot(2,1,1)
    plot(RPs(:,1),'-r');
    hold on;
    plot(RPs(:,2),'-g');
    plot(RPs(:,3),'-b');
    plot(dRPs(:,1),'--','Color',[1.0,0.5,0.5]);
    plot(dRPs(:,2),'--','Color',[0.5,1.0,0.5]);
    plot(dRPs(:,3),'--','Color',[0.5,0.5,1.0]);
    Legend = legend('X','Y','Z','\DeltaX','\DeltaY','\DeltaZ');
    set(Legend,'Location','eastoutside');
    title(sprintf('Translations for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Displacement /mm');
    ax = gca;
    ax.YAxis.Exponent = 0;   % turn off the ×10^n scaling
    ax.YRuler.TickLabelFormat = '%.2f';  % control decimal places
    ax.YTickLabelRotation = 45; 
    ymin = -1;
    ymax = 1.5 ;%hard coded based on inspecting the extremes
    ylim([ymin,ymax]);
    pause(.1);
    hold off;

    %rotations
    subplot(2,1,2)
    plot(RPs(:,4),'-m');
    hold on;
    plot(RPs(:,5),'-y');
    plot(RPs(:,6),'-c');
    plot(dRPs(:,4),'--','Color',[1.0,0.5,1.0]);
    plot(dRPs(:,5),'--','Color',[1.0,1.0,0.5]);
    plot(dRPs(:,6),'--','Color',[0.5,1.0,1.0]);
    Legend = legend('Pitch','Roll','Yaw','\DeltaPitch','\DeltaRoll','\DeltaYaw');
    set(Legend,'Location','eastoutside');
    title(sprintf('Rotations for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Rotation /rad');
    ax = gca;
    ax.YAxis.Exponent = 0;   % turn off the ×10^n scaling
    ax.YRuler.TickLabelFormat = '%.2f';  % control decimal places
    ax.YTickLabelRotation = 45;
    ymin = -0.025;
    % graphs are easily comparible between people
    ymax = 0.05; %this is hard coded based on inspecting the extremes
    ylim([ymin,ymax]);
    print('-dpng',sprintf('%s%s%s_R%i_TsAndRs.png',rpsFolder,filesep,subjectId,iR));
    pause(.1);
    hold off;
    close(gcf);

    % create FWD graph
    plot(fwd,'--','Color',[0,0.5,1]);
    title(sprintf('Framewise displacement for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Displacement /mm');
    %save FWD figure
    print('-dpng',sprintf('%s%s%s_R%i_FramewiseDisplacement.png',...
        rpsFolder,filesep,subjectId,iR));
    pause(.1);
    close(gcf)

    fnExpr = sprintf('%s%sR%i_%s%s%s',rpsFolder,filesep,...
        iR,subjectId,'ExtremeParams_R','.mat');
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return