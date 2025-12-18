function [Data] = getAlpha01Ts(subjectId,epiMask)
dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
nStim = 6;
nPos = 2;
positions = 'ab';
Data = nan(epiMask.size(1),epiMask.size(2),epiMask.size(3),nStim, nPos);
for iPos = 1:nPos
    cPos = positions(iPos);
    for iStim = 1:nStim
        stimId = [cPos,num2str(iStim-1)];
        tFn = [dirs.Alpha01,filesep,stimId,filesep,'spmT_0001.nii'];
        tV = spm_vol(tFn);
        tM = spm_read_vols(tV);
        Data(:,:,:,iStim,iPos) = tM;
    end
end
Data = cat(4,Data(:,:,:,:,1),Data(:,:,:,:,2));
return
