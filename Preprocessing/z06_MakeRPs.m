function [] = z06_MakeRPs(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir.folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];
rpsFolder = [dataDir,filesep,'RPs'];
if ~exist(rpsFolder,"dir")
    mkdir(rpsFolder);
end
anlysFolder = [subjDir,filesep,'Analysis',filesep,'Alpha00'];
if ~exist(anlysFolder,"dir")
    mkdir(anlysFolder);
end

ExtremeParams(1,1).MaxAbsT = NaN;
ExtremeParams(1,1).MaxAbsR = NaN;
ExtremeParams(1,1).MaxDerT = NaN;
ExtremeParams(1,1).MaxDerR = NaN;
ExtremeParams(1,1).MaxFwd = NaN;

rpFiles = dir([epiDir,filesep,'*.txt']);
runN = numel(rpFiles);




for iR = 1:runN
    % Get data:
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

    % Combine, make censor regressor and then z-score
    R = [RPs,dRPs];
    R = [R,R.^2]; % Add squares of og params and the derivs

    % Compute frame-wise displacement
    fwd = nan(nEpis-1,1);
    for ii= 1:(nEpis-1)
        ts = sum(abs(dRPs(ii,1:3)));
        rs = sum(abs(dRPs(ii,4:6)));
        fwd(ii,1) = ts + 50*rs;
    end


    %%
    % MaxAbs:
    Extremes_T = zeros(2,3);
    Extremes_T(1,:) = min(R(:,1:3));
    Extremes_T(2,:) = max(R(:,1:3));
    Extremes_R = zeros(2,3);
    Extremes_R(1,:) = min(R(:,4:6));
    Extremes_R(2,:) = max(R(:,4:6));
    Extremes_T = abs(Extremes_T);
    Extremes_R = abs(Extremes_R);
    MaxAbs_T = max(max(Extremes_T(:,1:3)));
    MaxAbs_R = max(max(Extremes_R(:,1:3)));
    ExtremeParams(1,1).MaxAbsT = MaxAbs_T;
    ExtremeParams(1,1).MaxAbsR = MaxAbs_R;

    % MaxDer:
    Extremes_dT = zeros(2,3);
    Extremes_dT(1,:) = min(R(:,7:9));
    Extremes_dT(2,:) = max(R(:,7:9));
    Extremes_dR = zeros(2,3);
    Extremes_dR(1,:) = min(R(:,10:12));
    Extremes_dR(2,:) = max(R(:,10:12));
    Extremes_dT = abs(Extremes_dT);
    Extremes_dR = abs(Extremes_dR);
    MaxDer_dT = max(max(Extremes_dT(:,1:3)));
    MaxDer_dR = max(max(Extremes_dR(:,1:3)));
    ExtremeParams(1,1).MaxDerT = MaxDer_dT;
    ExtremeParams(1,1).MaxDerR = MaxDer_dR;
    %Max framewise displacement:
    ExtremeParams.MaxFwd = max(fwd);

    threshold = 0.1 ;% pick meaningful threshold!!
    linIdx = find(abs(dRPs)>=threshold); %volume<->volume changes so only want the derivatives(columns 7-12)
    [volId,~] = ind2sub(size(dRPs),linIdx);
    volId = unique(volId); %this is needed as technically you could get multiple columns being above threshold for the same volume
    mask = zeros(nEpis,numel(volId));
    for mm=1:numel(volId)
        mask(volId(mm),mm)= 1;
    end

    R = zscore(R,[],1); %after using the raw values to find the vols to put in the censor regresor - now we can normalise
    R = [R,mask]; % append the censor regressor onto the other nuisance regressors
    %% Save nuissance regressors - maybe these should be concatenated?:
    save(sprintf('%s%sNuisance_R%i.mat',anlysFolder,filesep,iR),'R');

    %% Print plots:
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
    set(Legend,'Location','NorthWest');
    title(sprintf('Translations for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Displacement /mm');
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
    set(Legend,'Location','NorthWest');
    title(sprintf('Rotations for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Rotation /rad');
    print('-dpng',sprintf('%s%s%s_%i_TsAndRs.png',rpsFolder,filesep,subjectId,iR));
    pause(.1);
    hold off;
    close(gcf);

    % create FWD graph
    plot(fwd,'--','Color',[0,0.5,1]);
    title(sprintf('Framewise displacement for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Displacement /mm');
    %save FWD figure
    print('-dpng',sprintf('%s%s%s_%i_FramewiseDisplacement.png',...
        rpsFolder,filesep,subjectId,iR));
    pause(.1);
    close(gcf)

    fnExpr = sprintf('%s%s%s%i%s%s',rpsFolder,filesep,...
        'ExtremeParams_R',iR,subjectId,'.mat');
    save(fnExpr,'ExtremeParams'); % save to group location
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return
