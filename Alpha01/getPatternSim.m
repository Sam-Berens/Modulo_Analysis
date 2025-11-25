function [zTemplate,patternSim] = getPatternSim(G,roiId)

% zTemplate = [nSubjects,2];
% patternSim = [12,12,nSubjects];

subjectIds = getSubjectIds(G);
[~,order] = sort(subjectIds);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
H = (3-min(cat(3,mod((0:5)'-(0:5),6),mod((0:5)-(0:5)',6)),[],3))./3;
rmvDiag = boolean(tril(ones(6),-1) + triu(ones(6),1));
patternSim = nan([12,12,numel(subjectIds)]);
zTemplate = nan(numel(subjectIds),2);
for iSubject=1:numel(subjectIds)
    cId = subjectIds{iSubject};
    % prep for getTpatterns
    %check about whether this is meant to come out or not because it could
    %happen higher up but then we'd be passing dirs in instead of subjectId??
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
    dirs.EPI = [dirs.Subject,filesep,'EPI'];
    dirs.maskPath = sprintf('%s%s_%s_epiMask00.nii',dirs.EPI,filesep,cId);

    % Sample t-stat data using mask ROI coordinates
    %remember that Data is organised such as [nStim,nVox,nPos]
    Data = getTpatterns_EpiRes(G,cId,roiId,dirs);
    Data = [Data(:,:,1);Data(:,:,2)];
    %reminder that this M is (nStims,nVoxels)
    % Compute pairwise Euclidian distances across conditions
    D = pdist(Data,'euclidean');
    D = squareform(D); % TO DO - check this is in the right form - i think it has zeros in the diagonal but mb we want to remove them?
    %first chop off top triangle (just make zeros)
    lowerBig = tril(true(12),-1);
    lD = lowerBig .* D;
    %then you need to take the top 6 rows and reflect these into top right
    %of an empty square which is 6 by 6 ,
    aa = lD(1:6,1:6);
    ab_ba = D(7:12,1:6);
    bb = lD(7:12,7:12);
    %flip these lower triangle values to be in the top
    bb = permute(bb,[2,1]);
   % add together colocal triangles
    aa_bb = aa + bb;
    %not sure how it makes sense to correlate a square with a square? we should be collapsing both right?  
    % TO DO - FIND OUT WHETHER I'M MEANT TO BE CORRELATING LIKE THE UPPER
    % PART OF H WITH BB AND LOWER PART OF H WITH AA OR IF THAT DOESNT MAKE
    % SENSE
    r1 = corr(H(rmvDiag),aa_bb(rmvDiag),'Type','Kendall'); %this is the type of corr, between this and pearsons
    r2 = corr(H(rmvDiag),ab_ba(rmvDiag),'Type','Kendall'); 
    patternSim(:,:,iSubject) = D;
    z1 = atanh(r1);
    z2 = atanh(r2);
    zTemplate(iSubject,1:2) = [z1,z2];
end
%make sure they are ordered the same as other subjectId ordered arrays
patternSim = patternSim(:,:,order);
zTemplate = zTemplate(order,:);
return