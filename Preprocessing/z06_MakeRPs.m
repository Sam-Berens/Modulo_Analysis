function [MotionStats] = z06_MakeRPs(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir.folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];

rpsFolder = [dataDir,filesep,'RPs'];
if ~exist(rpsFolder,"dir")
    mkdir(rpsFolder);
end

censorThresh = 0.9; % mm fwd

rpFiles = dir([epiDir,filesep,'*.txt']);
runN = numel(rpFiles);
MotionStats = struct;

for iR = 1:runN
    % Get data
    FileName_RpText = [rpFiles(iR).folder,filesep,rpFiles(iR).name];
    RPs = readmatrix(FileName_RpText);
    nEpis = size(RPs,1);

    % Compute derivatives
    dRPs = diff(RPs,1,1);
    qt = (0:1:(nEpis-1))';
    t = mean([qt(1:end-1),qt(2:end)],2);
    dRPs = cell2mat(...
        cellfun(@(y)sincInterp(qt,t,y),mat2cell(dRPs,nEpis-1,ones(1,6)),...
        'UniformOutput',false));

    % Compute frame-wise displacement
    fwd = nan(nEpis-1,1);
    for ii = 1:(nEpis-1)
        translations = sum(abs(dRPs(ii,1:3)));
        rotations = sum(abs(dRPs(ii,4:6)));
        fwd(ii,1) = translations + 50*rotations;
    end
    fwd = sincInterp(qt,t,fwd);

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
    save(sprintf('%s%sRP%i.mat',epiDir,filesep,iR),'R');

    %% Populate MotionStats
    MotionStats(iR,1).subjectId = categorical({subjectId});
    MotionStats(iR,1).runIdx = iR;
    MotionStats(iR,1).disRang = ...
        max(RPs(:,1:3),[],1) - min(RPs(:,1:3),[],1);
    MotionStats(iR,1).totTrav = norm(MotionStats(iR,1).maxDisp);
    MotionStats(iR,1).rotRang = ...
        rad2deg(max(RPs(:,4:6),[],1) - min(RPs(:,4:6),[],1));
    % https://chatgpt.com/share/68a0a8dc-5028-8003-b235-dbab491e5fd6
    MotionStats(iR,1).rotTrav

    %% THIS NEEDS SOME WORK
    % MaxAbs:
    Extremes_T = zeros(2,3);
    Extremes_T(1,:) = min(RPs(:,1:3));
    Extremes_T(2,:) = max(RPs(:,1:3));
    Extremes_R = zeros(2,3);
    Extremes_R(1,:) = min(RPs(:,4:6));
    Extremes_R(2,:) = max(RPs(:,4:6));
    Extremes_T = abs(Extremes_T);
    Extremes_R = abs(Extremes_R);
    MaxAbs_T = max(max(Extremes_T(:,1:3)));
    MaxAbs_R = max(max(Extremes_R(:,1:3)));
    ExtremeParams(1,1).MaxAbsT = MaxAbs_T;
    ExtremeParams(1,1).MaxAbsR = MaxAbs_R;

    % MaxDer:
    Extremes_dT = zeros(2,3);
    Extremes_dT(1,:) = min(dRPs(:,1:3));
    Extremes_dT(2,:) = max(dRPs(:,1:3));
    Extremes_dR = zeros(2,3);
    Extremes_dR(1,:) = min(dRPs(:,4:6));
    Extremes_dR(2,:) = max(dRPs(:,4:6));
    Extremes_dT = abs(Extremes_dT);
    Extremes_dR = abs(Extremes_dR);
    MaxDer_dT = max(max(Extremes_dT(:,1:3)));
    MaxDer_dR = max(max(Extremes_dR(:,1:3)));
    ExtremeParams(1,1).MaxDerT = MaxDer_dT;
    ExtremeParams(1,1).MaxDerR = MaxDer_dR;

    % Max framewise displacement:
    ExtremeParams.MaxFwd = max(fwd);

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
    save(fnExpr,'ExtremeParams'); % save to group location
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return
