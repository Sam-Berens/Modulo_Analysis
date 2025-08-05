function [] = z06_MakeRPs(subjectId)

dataDir = ['..',filesep,'..',filesep,'Data'];
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];

ExtremeParams(1,1).MaxAbsT = NaN;
ExtremeParams(1,1).MaxAbsR = NaN;
ExtremeParams(1,1).MaxDerT = NaN;
ExtremeParams(1,1).MaxDerR = NaN;

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

    % Combine and z-score
    R = [RPs,dRPs]; % TODO: Add squares????
    %               % TODO: Add censor regressors!!!
    R = zscore(R,[],1);

    % Compute mean frame-wise displacement
    % TODO!

    %%
    % TODO: Pass Extremes out for collation?
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

    %% Save RPS:
    save(sprintf('%s%sRPs_%i.mat',epiDir,filesep,iR),'R');

    %% Print plots:
    % TODO: Extremes and plots should be saved in one place for all ppants

    plot(TX,'-r');
    hold on;
    plot(TY,'-g');
    plot(TZ,'-b');
    plot(ndTX,'--','Color',[1.0,0.5,0.5]);
    plot(ndTY,'--','Color',[0.5,1.0,0.5]);
    plot(ndTZ,'--','Color',[0.5,0.5,1.0]);

    Legend = legend('X','Y','Z','\DeltaX','\DeltaY','\DeltaZ');
    set(Legend,'Location','NorthWest');
    title(sprintf('Translations for %s',subjectId));
    xlabel('Scans \rightarrow');
    ylabel('Displacement /mm');
    %print('-dpng',sprintf('%s%s%s_Translations.png',epiDir,filesep,subjectId));
    pause(.1);
    hold off;
    %rotations
    plot(RX,'-m');
    hold on;
    plot(RY,'-y');
    plot(RZ,'-c');
    plot(ndRX,'--','Color',[1.0,0.5,1.0]);
    plot(ndRY,'--','Color',[1.0,1.0,0.5]);
    plot(ndRZ,'--','Color',[0.5,1.0,1.0]);

    Legend = legend('Pitch','Roll','Yaw','\DeltaPitch','\DeltaRoll','\DeltaYaw');
    set(Legend,'Location','NorthWest');
    title(sprintf('Rotations for %s',subjectId));
    xlabel('Scans \rightarrow');
    %ylabel('Rotation /\circ');
    %print('-dpng',sprintf('%s%s%s_Rotations.png',epiDir,filesep,subjectId));
    pause(.1);
    hold off;

    %save([epiDir,filesep,'ExtremeParams.mat'],'ExtremeParams');
end

return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return
