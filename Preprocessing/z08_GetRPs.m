function [] = MakeRPs(subjectId)

dataDir = dir(['..',filesep,'..',filesep,'Data']);
dataDir = dataDir(1).folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];

ExtremeParams(1,1).MaxAbsT = NaN;
ExtremeParams(1,1).MaxAbsR = NaN;
ExtremeParams(1,1).MaxDerT = NaN;
ExtremeParams(1,1).MaxDerR = NaN;
AFun_CalcDerivs = @(V) V - [V(1);V(1:end-1)];

rpFiles = dir([epiDir,filesep,'*.txt']);
runN = numel(rpFiles);

for iR = 1:runN

    %% Get data:
    FileName_RpText = [rpFiles(iR).folder,filesep,rpFiles(iR).name];

    RPs = readmatrix(FileName_RpText);
    nEpis = size(RPs,1);

    TX = RPs(:,1);
    TY = RPs(:,2);
    TZ = RPs(:,3);
    RX = RPs(:,4);
    RY = RPs(:,5);
    RZ = RPs(:,6);

    %% Convert from radians to degrees:
    RX = RX.*(360/(2*pi));
    RY = RY.*(360/(2*pi));
    RZ = RZ.*(360/(2*pi));

    %% Calculate derivatives:
    dTX = AFun_CalcDerivs(TX);  %these are your inputs for sinc interp
    dTX = dTX(2:end); %chop off the start value which is just change with respect to itself
    dTY = AFun_CalcDerivs(TY);
    dTY = dTY(2:end);
    dTZ = AFun_CalcDerivs(TZ);
    dTZ = dTZ(2:end);
    dRX = AFun_CalcDerivs(RX);
    dRX = dRX(2:end);
    dRY = AFun_CalcDerivs(RY);
    dRY = dRY(2:end);
    dRZ = AFun_CalcDerivs(RZ);
    dRZ = dRZ(2:end);

    %define your old and new timepoints to sample from
    ogSampSpace = (0.5:1:292); %the sample space which calculating the derivatives will produce (n-1)
    newSampSpace = (0:292); %this is the new sample space we want so that there is a derivative for each scan (size n)

    %sinc Interpolate derivatives to their correct locations (half way between scans)
    ndTX = sincInterp(ogSampSpace, newSampSpace, dTX); %overwrite derivs with those corresponding to the scan vol sample space
    ndTY = sincInterp(ogSampSpace, newSampSpace, dTY);
    ndTZ = sincInterp(ogSampSpace, newSampSpace, dTZ);
    ndRX = sincInterp(ogSampSpace, newSampSpace, dRX);
    ndRY = sincInterp(ogSampSpace, newSampSpace, dRY);
    ndRZ = sincInterp(ogSampSpace, newSampSpace, dRZ);


    %% Fill the matrix R:
    R = zeros(nEpis,12);
    R(:,01) = TX;
    R(:,02) = TY;
    R(:,03) = TZ;
    R(:,04) = RX;
    R(:,05) = RY;
    R(:,06) = RZ;
    R(:,07) = ndTX;
    R(:,08) = ndTY;
    R(:,09) = ndTZ;
    R(:,10) = ndRX;
    R(:,11) = ndRY;
    R(:,12) = ndRZ;

    %% MaxAbs:
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

    %% MaxDer:
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

    %% Zscore and save:
    R = zscore(R,[],1);
    R = [R,zscore(R.^2,[],1)]; %#ok<AGROW>
    save(sprintf('%s%sRPs_%i.mat',epiDir,filesep,iR),'R'); %remeber to change this back to linux filesep
    %% Print plots:
    %translations
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
    print('-dpng',sprintf('%s%s%s_Translations.png',epiDir,filesep,subjectId));
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
    ylabel('Rotation /\circ');
    print('-dpng',sprintf('%s%s%s_Rotations.png',epiDir,filesep,subjectId));
    pause(.1);
    hold off;

    save([epiDir,filesep,'ExtremeParams.mat'],'ExtremeParams');

end


return


function [yNew] = sincInterp(tOld, tNew, yOld)
[TNew, TOld] = ndgrid(tNew,tOld); %sinc expects a grid so it can compare each old t with each new t
S = sinc(TNew-TOld); %tIntval is 1 so not needed
yNew = S * yOld; 
return
