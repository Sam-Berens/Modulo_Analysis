function [MotionStats] = getMotionStats()

censorThresh = 0.9; % mm FWD
MotionStats = struct;

% Set dataDir and subjList
dataList = dir(['..',filesep,'..',filesep,'Data']);
dataList = dataList(cellfun(@(s)numel(s)==8,{dataList.name}));
dataDir = dataList.folder;
subjList = {dataList.name}';

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

        %% Populate MotionStats
        MotionStats(iIn,1).subjectId = categorical({subjectId});
        MotionStats(iIn,1).runIdx = iR;

        % Ovarial displacement
        MotionStats(iIn,1).tra = ...
            max(RPs(:,1:3),[],1) - min(RPs(:,1:3),[],1);
        MotionStats(iIn,1).traSum = norm(MotionStats(iIn,1).tra);

        % Overall rotation
        rot = max(RPs(:,4:6),[],1) - min(RPs(:,4:6),[],1);
        MotionStats(iIn,1).rot = rad2deg(rot);
        %%% Total rotation is computed using the axis-angle representation.
        %%% theta = acos((tr(R)-1)/2);
        %%% Note that the -2 below accounts of the 4th dimension in R.
        MotionStats(iIn,1).rotSum = ...
            rad2deg(acos((trace(spm_matrix([zeros(1,3),rot]))-2)/2));

        % Summary stats for the norm of the translation derivatives
        dTra = sqrt(sum(dRPs(:,1:3).^2,2));
        MotionStats(iIn,1).dTraMea = mean(dTra);
        MotionStats(iIn,1).dTraMed = median(dTra);
        MotionStats(iIn,1).dTraMax = max(dTra);

        % Summary stats for the total rotation derivative
        dRot = nan(nEpis,1);
        for ii = 1:nEpis
            M = spm_matrix([zeros(1,3),dRPs(ii,4:6)]);
            dRot(ii) = abs(rad2deg(acos((trace(M)-2)/2)));
        end
        MotionStats(iIn,1).dRotMea = mean(dRot);
        MotionStats(iIn,1).dRotMed = median(dTra);
        MotionStats(iIn,1).dRotMax = max(dRot);

        % Summary stats for FWD
        MotionStats(iIn,1).fwdMea = mean(fwd);
        MotionStats(iIn,1).fwdMed = median(fwd);
        MotionStats(iIn,1).fwdMax = max(fwd);

        % Number of censored volumes
        MotionStats(iIn,1).nCen = numel(volIdxs);
    end
end

MotionStats = struct2table(MotionStats);
return

function [qy] = sincInterp(qt,t,y)
[Tq, T] = ndgrid(qt,t);
qy = sinc(Tq-T) * y;
return