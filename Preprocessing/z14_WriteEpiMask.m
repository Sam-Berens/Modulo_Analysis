function [] = z14_WriteEpiMask(subjectId)
[V,Y] = computeMask(subjectId);
spm_write_vol(V,Y);
return

function [V,Y] = computeMask(subjectId)
% Collect all the EPI data into a nVox*time matrix (M)
dataDir = dir(['..',filesep,'..',filesep, 'Data']);
dataDir = dataDir.folder;
subjDir = [dataDir,filesep,subjectId];
epiDir = [subjDir,filesep,'EPI'];
tempDir = [epiDir,filesep,'2_Temporal'];
runList = dir([tempDir,filesep,'R*']);

for iRun = 1:numel(runList)
    runN = runList(iRun).name;
    xepiList = dir([tempDir,filesep,runN,filesep,'au_*.nii']);
    epiList = arrayfun(@(x)[x.folder,filesep,x.name],xepiList,'UniformOutput',false);
    V = spm_vol(epiList);
    nM = cellfun(@(v)spm_read_vols(v),V,'UniformOutput',false);
    nM = reshape(nM,1,1,1,numel(nM));
    nM = cell2mat(nM);
    nM = reshape(nM,[],size(nM,4));
    if iRun == 1
        M = nM;
    else
        M = [M,nM]; %#ok<AGROW>
    end
end

%% Compute threshold for log-mean
lMean = log(mean(M,2));
[f,x] = ecdf(lMean);
[~,iMin] = min((f-0.6).^2);
tlMean = x(iMin);

%% Compute threshold for log-var
lVar = log(var(M,[],2));
[f,x] = ecdf(lVar);
[~,iMin] = min((f-0.6).^2);
tlVar = x(iMin);

%% Threshold
Y = (lMean>tlMean) | (lVar>tlVar);

%% Prep the volume structure
V = V{1};
V.fname = sprintf('%s%s_%s_epiMask00.nii',epiDir,filesep,subjectId);
V.dt = [2,0];
V.descrip = 'Custom EPI mask (v00): cdf(log-Mean)>0.6 | cdf(log-Var)>0.6';

%% Reshape
Y = reshape(Y,V.dim);
return