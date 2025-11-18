function [zTemplate,patternSim] = getPatternSim(G,roiId)

% zTemplate = [nSubjects,1 OR 2 depending on alpha00 or alpha01];
% patternSim = [6,6,nSubjects];

subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
H = (3-min(cat(3,mod((0:5)'-(0:5),6),mod((0:5)-(0:5)',6)),[],3))./3;
lower = tril(true(6),-1);
patternSim = nan([6,6,numel(subjectIds)]);
zTemplate = nan(numel(subjectIds),1);
for iSubject=1:numel(subjectIds)
    cId = subjectIds{iSubject};
    %% prep for getTpatterns
    %check about whether this is meant to come out or not because it could
    %happen higher up but then we'd be passing dirs in instead of subjectId??
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];
    dirs.EPI = [dirs.Subject,filesep,'EPI'];
    dirs.maskPath = sprintf('%s%s_%s_epiMask00.nii',dirs.EPI,filesep,cId);

    %% Sample t-stat data using mask ROI coordinates
    Data = getTpatterns_EpiRes(G,cId,roiId,dirs);
    %reminder that this M is (nStims,nVoxels)
    % Compute pairwise Euclidian distances across conditions
    D = pdist(Data,'euclidean');
    D = squareform(D); % TO DO - check this is in the right form - i think it has zeros in the diagonal but mb we want to remove them?
    r = corr(H(lower),D(lower),'Type','Kendall'); %this is the type of corr, between this and pearsons
    patternSim(:,:,iSubject) = D;
    z = atanh(r);
    zTemplate(iSubject,1) = z;
end
return